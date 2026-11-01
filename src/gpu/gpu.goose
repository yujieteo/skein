// The GPU path of the page of `skein build --target wasm` (docs/roadmap.md,
// phase 4.5). Only the program of the page imports this package: its
// externs are the WebGPU calls of the page (GPU_C in
// src/codegen/page.goose, and src/gpu/gpu.js), and the native profile run
// of the page links stubs that report no GPU.
//
// The program of the page calls these functions in place of the solves of
// the packages (src/codegen/build.goose, page_program):
//
//     central_difference   elements::central_difference: an explicit step,
//                          on the GPU whenever the page has WebGPU
//     grid_transient       grid::solve_transient (ssprk3 of fv and dg): the
//                          same
//     solve                elements::solve: an implicit step, on the GPU
//                          only with the control of the page on
//                          (SKEIN_GPU_IMPLICIT=1), else CG with AMG
//     grid_static, modal,  grid::solve_static, elements::modal,
//     buckling,            elements::buckling and meshless::solve_static:
//     meshless_static      on the CPU
//
// Each falls back to the CPU solve when the page has no WebGPU. This
// package sets up the GPU: it converts the arrays of the CPU solve to f32
// buffers, makes a bind group for each kernel of kernels.wgsl, records the
// dispatches, and reads the results back. The GPU keeps the state of the
// steps between read-backs. Each solve appends its device to the record
// `device` of the run.
//
// The implicit path is CG with the Jacobi preconditioner on the GPU, in
// f32, inside a defect correction in f64 on the CPU: each pass solves
// K d = r to INNER_RTOL, adds d to u and computes r = f - K u again, until
// the residual of the CPU solve (elements::CG_TOLERANCE) holds. The result
// thus has the accuracy of the CPU solve.
import std;
import math;
import os;
namespace gpu;
import elements from "../elements";
import grid from "../grid";
import meshless from "../meshless";
import solvers from "../solvers";

extern fn skein_gpu_ready() -> i32;
extern fn skein_gpu_reset();
extern fn skein_gpu_buffer(bytes: i64) -> i32;
extern fn skein_gpu_write(buf: i32, at: i64, data: const f32[:]);
extern fn skein_gpu_write_ints(buf: i32, at: i64, data: const i32[:]);
extern fn skein_gpu_bind(kernel: const u8[:], pairs: const i32[:]) -> i32;
extern fn skein_gpu_dispatch(group: i32, threads: i64);
extern fn skein_gpu_submit();
extern fn skein_gpu_read(buf: i32, at: i64, out: f32[:]) -> i32;

let WG = 64;            // the workgroup size of kernels.wgsl
let ROW_GROUPS = 65535; // the most workgroups in a row of a dispatch (gpu.js)
let BATCH = 256;        // the steps of central_difference in one submit
let CG_BATCH = 32;      // the cg iterations between two read-backs
let INNER_RTOL = 1.0e-3;    // the residual of a pass of cg on the GPU, relative to its start
let PASSES = 30;            // the most passes of the defect correction
let CG_MAX = 20000;         // the most cg iterations of a pass

// The bindings of kernels.wgsl.
let B_P: i32 = 0;
let B_ROWS: i32 = 1;
let B_COLS: i32 = 2;
let B_VALS: i32 = 3;
let B_U: i32 = 4;
let B_V: i32 = 5;
let B_ACC: i32 = 6;
let B_W: i32 = 7;
let B_F: i32 = 8;
let B_X0: i32 = 10;
let B_XV: i32 = 11;
let B_XO: i32 = 12;
let B_GFACE: i32 = 13;
let B_FBC: i32 = 14;
let B_GRAD: i32 = 15;
let B_PHI: i32 = 16;
let B_STATS: i32 = 17;
let B_CX: i32 = 20;
let B_CR: i32 = 21;
let B_CP: i32 = 22;
let B_CQ: i32 = 23;
let B_PART: i32 = 25;
let B_SC: i32 = 26;

var ready = -1;                 // skein_gpu_ready, or -1 before the first solve
var step: u8[..64] = [];        // the step that the program runs
var devices: u8[>..] = [];
var device_count = 0;
// The f32 copies of an upload or a read-back go through these chunks of
// CHUNK values. Reason: the native profile run of the page has no GPU, so
// it does not measure a growable array of this package, and a copy of the
// matrix of sample 01 overran its data stack in the page.
let CHUNK = 65536;
var chunk: f32[65536] = [0.0; 65536];
var ichunk: i32[65536] = [0; 65536];
var fill = 0;           // the values in chunk or ichunk
var fill_at = 0;        // the index in the buffer of chunk[0]
var fill_buf: i32 = -1;

// The program of the page calls this at the start of each step.
fn step_name(name: u8[:]) {
    step.clear();
    step.append(name[..min(name.len, 64)]);
}

fn note(solver: u8[:], on_gpu: bool) {
    format(devices, "  ", step, ",", solver, ",", if on_gpu { "gpu" } else { "cpu" }, "\n");
    device_count++;
}

// The record of the device of each solve, in TOON, after the records of
// the run.
fn records() -> u8[>..] { return str("device[", device_count, "]{step,solver,device}:\n", devices); }

fn present() -> bool {
    if ready < 0 { ready = skein_gpu_ready(); }
    return ready == 1;
}

fn implicit_on() -> bool {
    var v: u8[>..] = [];
    return env("SKEIN_GPU_IMPLICIT", v) && v == "1";
}

// ---- buffers and dispatches ----

// The writes of buffer b from its start: put or put_int each value, then
// flush or flush_ints.
fn put_start(b: i32) {
    fill_buf = b;
    fill = 0;
    fill_at = 0;
}

fn put(x: f64) {
    chunk[fill] = x as f32;
    fill++;
    if fill == CHUNK { flush(); }
}

fn flush() {
    if fill > 0 { skein_gpu_write(fill_buf, fill_at, chunk[..fill]); }
    fill_at += fill;
    fill = 0;
}

fn put_int(x: i64) {
    ichunk[fill] = x as! i32;
    fill++;
    if fill == CHUNK { flush_ints(); }
}

fn flush_ints() {
    if fill > 0 { skein_gpu_write_ints(fill_buf, fill_at, ichunk[..fill]); }
    fill_at += fill;
    fill = 0;
}

fn upload(data: f64[:]) -> i32 {
    let b = skein_gpu_buffer(4 * data.len);
    put_start(b);
    for x in data { put(x); }
    flush();
    return b;
}

fn zeros(n: i64) -> i32 { return skein_gpu_buffer(4 * n); }

fn upload_ints(data: i64[:]) -> i32 {
    let b = skein_gpu_buffer(4 * data.len);
    put_start(b);
    for x in data { put_int(x); }
    flush_ints();
    return b;
}

fn upload_cols(data: const i32[:]) -> i32 {
    let b = skein_gpu_buffer(4 * data.len);
    put_start(b);
    for x in data { put_int(x as i64); }
    flush_ints();
    return b;
}

// A bind group of `kernel`: binding pairs[2 k] is buffer pairs[2 k + 1].
fn bind(kernel: u8[:], pairs: const i32[:]) -> i32 { return skein_gpu_bind(kernel, pairs); }

fn dispatch(group: i32, threads: i64) { skein_gpu_dispatch(group, threads); }

// The workgroups of a dispatch of n threads, with the empty ones of its
// last row: the sums of the cg kernels have one pair for each.
fn workgroups(n: i64) -> i64 {
    let g = max(1, (n + WG - 1) / WG);
    let x = min(g, ROW_GROUPS);
    return x * ((g + x - 1) / x);
}

// Reads buffer b back into out. False if the GPU failed.
fn read(b: i32, out: f64[:]) -> bool {
    var at = 0;
    while at < out.len {
        let m = min(CHUNK, out.len - at);
        guard skein_gpu_read(b, at, chunk[..m]) == 1 else { return false; };
        for i in m { out[at + i] = chunk[i]; }
        at += m;
    }
    return true;
}

fn finite(v: f64[:]) -> bool {
    for x in v { if !(abs(x) < 1.0e300) { return false; } }
    return true;
}

// ---- central_difference ----

// The transient step central_difference of elements, on the GPU when the
// page has WebGPU.
fn central_difference(end: f64, dt: f64) -> bool {
    guard present() else {
        note("central_difference", false);
        return elements::central_difference(end, dt);
    };
    note("central_difference", true);
    let start = clock();
    guard elements::cd_prepare(end, dt) else { return false; };
    let n = elements::sol.len;
    // the acceleration at the start, in f64 on the CPU
    solvers::rows_start(elements::row_start, elements::col, elements::val);
    elements::multiply(elements::sol, elements::accel);
    solvers::rows_stop();
    elements::cd_accel();
    skein_gpu_reset();
    let p = upload([n * 1.0, elements::cd_dt]);
    let rows = upload_ints(elements::row_start);
    let cols = upload_cols(elements::col);
    // M⁻¹ K and M⁻¹ f, so that cd_kick fits in 8 buffers
    let vals = skein_gpu_buffer(4 * elements::val.len);
    put_start(vals);
    for d in n {
        for e in elements::row_start[d]..elements::row_start[d + 1] { put(elements::minv[d] * elements::val[e]); }
    }
    flush();
    let f = zeros(n);
    put_start(f);
    for d in n { put(elements::minv[d] * elements::force[d]); }
    flush();
    let u = upload(elements::sol);
    let v = upload(elements::velocity);
    let a = upload(elements::accel);
    let drift = bind("cd_drift", [B_P, p, B_U, u, B_V, v, B_ACC, a]);
    let kick = bind("cd_kick", [B_P, p, B_ROWS, rows, B_COLS, cols, B_VALS, vals, B_U, u, B_V, v, B_ACC, a, B_F, f]);
    var done = 0;
    while done < elements::cd_steps {
        let k = min(BATCH, elements::cd_steps - done);
        for i in k {
            dispatch(drift, n);
            dispatch(kick, n);
        }
        skein_gpu_submit();
        done += k;
    }
    let ok = read(u, elements::sol) && read(v, elements::velocity);
    skein_gpu_reset();
    elements::last_iterations = elements::cd_steps;
    elements::last_seconds = clock() - start;
    return ok && finite(elements::sol);
}

// ---- ssprk3 of fv and dg ----

// The buffers of a transient grid solve: the stage kernels of each parity
// of the step (the solution alternates between two buffers) and each
// stage, and the kernels of fv before a stage, for each stage input.
struct GridGroups {
    stages: i32[6],
    grads: i32[3],
    venkats: i32[3],
    minmax: i32[3],
}

// The parameters of a stage (c0, c1) of ssprk3 with the step k.
fn grid_params(c0: f64, c1: f64, k: f64) -> i32 {
    let n = grid::dof_count() * 1.0;
    if grid::method == grid::DG { return upload([n, c0, c1, k]); }
    let l = sqrt(grid::nx * grid::hx * grid::ny * grid::hy);
    let kh = grid::VENKAT_K * sqrt(grid::hx * grid::hy) / l;
    return upload([n, grid::nx * 1.0, grid::ny * 1.0, grid::hx, grid::hy, grid::ax, grid::ay, grid::diffusivity,
                   grid::flux * 1.0, grid::gradient * 1.0, grid::limiter * 1.0, kh * kh * kh, c0, c1, k, 1.0 / (grid::ct * grid::hx * grid::hy)]);
}

// The transient step of fv or dg (ssprk3), on the GPU when the page has
// WebGPU. The source and the boundary values of each stage come from the
// CPU (grid::stage_terms); the GPU does the rest of each stage.
fn grid_transient<SRC, BC>(end: f64, dt: f64) -> bool {
    guard present() else {
        note("ssprk3", false);
        return grid::solve_transient(end, dt, SRC, BC);
    };
    note("ssprk3", true);
    let start = clock();
    let k, steps = grid::transient_start(end, dt);
    let n = grid::dof_count();
    let fv = grid::method == grid::FV;
    skein_gpu_reset();
    let c0 = [0.0, 0.75, 1.0 / 3.0];
    let c1 = [1.0, 0.25, 2.0 / 3.0];
    let at = [0.0, 1.0, 0.5];   // the time of each stage, in steps after the time of the step
    var ps: i32[3] = [0; 3];
    for s in 3 { ps[s] = grid_params(c0[s], c1[s], k); }
    // the solution in xs[0] or xs[1], and the second stage in xs[2]
    let xs = [upload(grid::sol), zeros(n), zeros(n)];
    let f = zeros(n);
    // fv: the face values and then the source of each cell
    let gf = zeros(if fv { 5 * n } else { 1 });
    let fbc = upload_ints(grid::face_bc);
    let grad = zeros(2 * n);
    let phi = zeros(n);
    let stats = zeros(2);
    var gr = GridGroups { stages: [0; 6], grads: [0; 3], venkats: [0; 3], minmax: [0; 3] };
    if fv {
        for x in 3 {
            gr.grads[x] = bind("fv_grad", [B_P, ps[0], B_XV, xs[x], B_GFACE, gf, B_FBC, fbc, B_GRAD, grad]);
            gr.minmax[x] = bind("fv_minmax", [B_P, ps[0], B_XV, xs[x], B_STATS, stats]);
            gr.venkats[x] = bind("fv_venkat", [B_P, ps[0], B_XV, xs[x], B_GFACE, gf, B_FBC, fbc, B_GRAD, grad, B_PHI, phi, B_STATS, stats]);
        }
    }
    let rows = upload_ints(grid::row_start);
    let cols = upload_cols(grid::col);
    // dg: the rows of A scaled by the stage weight, so that rk_stage fits
    // in 8 buffers
    let vals = skein_gpu_buffer(4 * grid::val.len);
    put_start(vals);
    for d in if fv { 0 } else { n } {
        for e in grid::row_start[d]..grid::row_start[d + 1] { put(grid::stage_weight(d) * grid::val[e]); }
    }
    flush();
    for parity in 2 {
        let sol = parity;
        let other = 1 - parity;
        // stage s reads ins[s] and writes outs[s]; the new solution is in `other`
        let ins = [sol, other, 2];
        let outs = [other, 2, other];
        for s in 3 {
            let gi = 3 * parity + s;
            if fv {
                gr.stages[gi] = bind("fv_stage", [B_P, ps[s], B_X0, xs[sol], B_XV, xs[ins[s]], B_XO, xs[outs[s]], B_GFACE, gf,
                                                  B_FBC, fbc, B_GRAD, grad, B_PHI, phi]);
            } else {
                gr.stages[gi] = bind("rk_stage", [B_P, ps[s], B_ROWS, rows, B_COLS, cols, B_VALS, vals, B_X0, xs[sol],
                                                  B_XV, xs[ins[s]], B_XO, xs[outs[s]], B_F, f]);
            }
        }
    }
    let t0 = grid::time;
    for i in steps {
        let parity = i % 2;
        let ins = [parity, 1 - parity, 2];
        for s in 3 {
            put_start(if fv { gf } else { f });
            grid::stage_terms(t0 + (i + at[s]) * k, SRC, BC, put);
            flush();
            if fv && grid::flux == grid::LINEAR_UPWIND {
                dispatch(gr.grads[ins[s]], n);
                if grid::limiter == grid::VENKATAKRISHNAN {
                    dispatch(gr.minmax[ins[s]], WG);
                    dispatch(gr.venkats[ins[s]], n);
                }
            }
            dispatch(gr.stages[3 * parity + s], n);
        }
    }
    let ok = read(xs[steps % 2], grid::sol);
    skein_gpu_reset();
    let done = grid::transient_end(end, steps);
    grid::last_seconds = clock() - start;
    return ok && done;
}

// The static step of fv or dg, the modal and buckling steps of elements,
// and the static step of rkpm: on the CPU, with a record of the device.
fn grid_static<SRC, BC>(rtol: f64) -> bool {
    note("gmres ilu0", false);
    return grid::solve_static(rtol, SRC, BC);
}

fn modal(count: i64) -> bool {
    note("lobpcg", false);
    return elements::modal(count);
}

fn buckling(count: i64) -> bool {
    note("lobpcg", false);
    return elements::buckling(count);
}

fn meshless_static<SRC, BC>() -> bool {
    note("skyline lu", false);
    return meshless::solve_static(SRC, BC);
}

// ---- the implicit solve ----

// The buffers and the kernels of cg on the GPU.
struct Cg {
    x: i32,
    r: i32,
    sc: i32,
    start: i32,
    matvec: i32,
    update: i32,
    direction: i32,
    shift: i32,
    sums: i32[3],    // the reduce of cg_start, cg_matvec and cg_update
}

fn cg_setup(n: i64, dinv: f64[:]) -> Cg {
    let p = upload([n * 1.0]);
    let rows = upload_ints(elements::row_start);
    let cols = upload_cols(elements::col);
    let vals = upload(elements::val);
    let w = upload(dinv);
    let x = zeros(n);
    let r = zeros(n);
    let pp = zeros(n);
    let q = zeros(n);
    let groups = workgroups(n);
    let part = zeros(2 * groups);
    let sc = zeros(8);
    var c = Cg { x: x, r: r, sc: sc, start: 0, matvec: 0, update: 0, direction: 0, shift: 0, sums: [0; 3] };
    c.start = bind("cg_start", [B_P, p, B_W, w, B_CR, r, B_CP, pp, B_PART, part]);
    c.matvec = bind("cg_matvec", [B_P, p, B_ROWS, rows, B_COLS, cols, B_VALS, vals, B_W, w, B_CP, pp, B_CQ, q, B_PART, part]);
    c.update = bind("cg_update", [B_P, p, B_W, w, B_CX, x, B_CR, r, B_CP, pp, B_CQ, q, B_PART, part, B_SC, sc]);
    c.direction = bind("cg_direction", [B_P, p, B_W, w, B_CR, r, B_CP, pp, B_SC, sc]);
    c.shift = bind("cg_shift", [B_SC, sc]);
    let slots = [0.0, 1.0, 2.0, 4.0, 3.0, 1.0];
    for k in 3 { c.sums[k] = bind("reduce", [B_P, upload([groups * 1.0, slots[2 * k], slots[2 * k + 1]]), B_PART, part, B_SC, sc]); }
    return c;
}

// Solves K d = r on the free dofs on the GPU to INNER_RTOL, with d in out.
// It gives the iterations, and false if the GPU failed or cg broke down.
fn cg_pass(c: Cg, n: i64, r: f64[:], out: f64[:]) -> i64, bool {
    put_start(c.r);
    for x in r { put(x); }
    flush();
    put_start(c.x);
    for x in r { put(0.0); }
    flush();
    dispatch(c.start, n);
    dispatch(c.sums[0], WG);
    var sc: f64[8] = [0.0; 8];
    guard read(c.sc, sc[..]) && sc[1] > 0.0 else { return 0, false; };
    let target = INNER_RTOL * INNER_RTOL * sc[1];
    var its = 0;
    while sc[1] > target && its < CG_MAX {
        for i in CG_BATCH {
            dispatch(c.matvec, n);
            dispatch(c.sums[1], WG);
            dispatch(c.update, n);
            dispatch(c.sums[2], WG);
            dispatch(c.direction, n);
            dispatch(c.shift, 1);
        }
        its += CG_BATCH;
        guard read(c.sc, sc[..]) && sc[1] >= 0.0 && sc[1] < 1.0e300 else { return its, false; };
    }
    return its, read(c.x, out);
}

// The static solve of elements: on the GPU with the control of the page
// on, else elements::solve.
fn solve() -> bool {
    let wave = elements::physics == elements::WAVE;
    guard present() && implicit_on() && !wave && elements::rad_count() == 0 else {
        note(if wave { "gmres ilu0" } else { "cg amg" }, false);
        return elements::solve();
    };
    note("cg jacobi", true);
    let start = clock();
    elements::assemble("simplex cg amg", 5);
    elements::hold_dead_nodes();
    let n = elements::sol.len;
    let sizes = elements::state_sizes();
    var r: f64[>..] = [];
    var d: f64[>..] = [];
    var q: f64[>..] = [];
    var dinv: f64[>..] = [];
    solvers::zeros(r, n);
    solvers::zeros(d, n);
    solvers::zeros(q, n);
    solvers::zeros(dinv, n);
    var f2 = 0.0;
    for e in n {
        elements::sol[e] = if elements::fixed[e] != 0 { elements::prescribed[e] } else { 0.0 };
        if elements::fixed[e] != 0 { continue; }
        let fe = elements::load[e] + elements::assembled_source[e];
        f2 += fe * fe;
        let k = elements::val[elements::entry(e, e)];
        dinv[e] = if k > 0.0 { 1.0 / k } else { 0.0 };
    }
    solvers::rows_start(elements::row_start, elements::col, elements::val);
    skein_gpu_reset();
    let c = cg_setup(n, dinv);
    var its = 0;
    var ok = true;
    var passes = 0;
    var rr = residual(r, q);
    let target = elements::CG_TOLERANCE * sqrt(max(f2, rr));
    while ok && sqrt(rr) > target && passes < PASSES {
        let k, good = cg_pass(c, n, r, d);
        its += k;
        ok = good;
        for e in n { elements::sol[e] += d[e]; }
        rr = residual(r, q);
        passes++;
    }
    solvers::rows_stop();
    skein_gpu_reset();
    guard elements::state_sizes() == sizes else { abort("gpu solve: a solve grows no array of the state (rule 3)"); }
    elements::last_iterations = its;
    elements::last_seconds = clock() - start;
    return ok && sqrt(rr) <= target;
}

// r = f - K u on the free dofs and 0 on the fixed dofs, in f64, with q for
// K u. It gives r·r.
fn residual(r: f64[:], q: f64[:]) -> f64 {
    elements::multiply(elements::sol, q);
    var rr = 0.0;
    for e in r.len {
        r[e] = if elements::fixed[e] != 0 { 0.0 } else { elements::load[e] + elements::assembled_source[e] - q[e] };
        rr += r[e] * r[e];
    }
    return rr;
}
