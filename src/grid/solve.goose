// The solves: a static step by defect correction with GMRES and ILU(0)
// (src/solvers/), and a transient step by the explicit SSP Runge-Kutta
// scheme of order 3 (Shu and Osher 1988).
//
// A static step solves R(u) = 0, with R the residual of fv (fv.goose) or
// A u - b of dg (dg.goose). Each iteration solves J δ = -R(u) and adds δ
// to u. J is the matrix of fv_matrix or A. It is the exact Jacobian of
// the linear forms, and then one iteration solves the step; with
// `linear_upwind`, the iterations converge to the residual of the second
// order flux.
//
// A transient step integrates ct M du/dt = -R(u, t), with M the volume of
// each cell (fv) or the diagonal mass matrix (dg).
import std;
import math;
import os;
namespace grid;
import solvers from "../solvers";
import .mesh;
import .fv;
import .dg;

let DEFECT_MAX = 200;           // the most iterations of the defect correction
// The most iterations with no new lowest residual. Reason: the residual of
// minmod does not fall at each iteration, and the residual of dg with
// p = 3 stops at the round-off of the matrix above an rtol of 1e-12.
let DEFECT_STALL = 20;
let LINEAR_RTOL = 1.0e-8;       // GMRES, relative to the norm of the residual
let LINEAR_RESTART = 100;
let LINEAR_MAX = 5000;
let CFL = 0.5;                  // the fraction of the stable step that `dt: stable` takes

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];
var lu: f64[>..] = [];
var diag: i64[>..] = [];
var grad_u: f64[>..] = [];      // fv: the gradient of each cell
var phi: f64[>..] = [];         // fv: the limiter of each cell
var rhs: f64[>..] = [];         // dg: b(t)
var r: f64[>..] = [];
var du: f64[>..] = [];
var u1: f64[>..] = [];
var u2: f64[>..] = [];

var iterations = 0;             // of the last solve: defect corrections, or time steps
var linear_iterations = 0;
var last_seconds = 0.0;
var first_residual = 0.0;
var last_residual = 0.0;

// The pattern of the cells beside each other: the rows of the dofs of a
// cell have the columns of the dofs of the cell and of its neighbours.
fn pattern() {
    row_start.clear();
    col.clear();
    row_start.push(0);
    let n = per_cell();
    for c in cell_count() {
        let cells = [neighbour(c, 2), neighbour(c, 0), c, neighbour(c, 1), neighbour(c, 3)];
        for a in n {
            for m in cells {
                if m < 0 { continue; }
                for b in n { col.push((m * n + b) as i32); }
            }
            row_start.push(col.len);
        }
    }
}

// The arrays of the solves, sized once for the method.
fn size_solver() {
    pattern();
    let n = dof_count();
    solvers::zeros(val, col.len);
    solvers::zeros(lu, col.len);
    diag.clear();
    for d in n { diag.push(0); }
    solvers::zeros(grad_u, if method == FV { 2 * n } else { 0 });
    solvers::zeros(phi, if method == FV { n } else { 0 });
    solvers::zeros(rhs, n);
    solvers::zeros(r, n);
    solvers::zeros(du, n);
    solvers::zeros(u1, n);
    solvers::zeros(u2, n);
}

// The residual R(u) at time t, in out.
fn residual<SRC, BC>(u: f64[:], t: f64, out: f64[:]) {
    if method == FV {
        fv_residual(u, t, out, SRC, BC);
        return;
    }
    solvers::multiply(row_start, col, val, u, out);
    for d in out.len { out[d] -= rhs[d]; }
}

// The matrix of the method in val: J of fv, or A of dg.
fn matrix() { if method == FV { fv_matrix(); } else { dg_matrix(); } }

// Solves the static step at the time of the state, to a residual of rtol
// times the first. False if GMRES or the iterations do not converge.
fn solve_static<SRC, BC>(rtol: f64) -> bool {
    guard method >= 0 && rtol > 0.0 && sol.len == dof_count() && val.len == col.len else { abort("solve_static: a method, clear_bcs, and rtol > 0"); }
    let start = clock();
    let sizes = state_sizes();
    solvers::memory_check("grid gmres ilu0", col.len, dof_count(), LINEAR_RESTART + 8, 0, "");
    matrix();
    if method == DG { dg_rhs(time, rhs, SRC, BC); }
    var converged = solvers::ilu0(row_start, col, val, lu, diag);
    residual(sol, time, r, SRC, BC);
    first_residual = solvers::norm(r);
    last_residual = first_residual;
    iterations = 0;
    linear_iterations = 0;
    var lowest = first_residual;
    var lowest_at = 0;
    while converged && last_residual > rtol * first_residual {
        if iterations >= DEFECT_MAX {
            converged = false;
            break;
        }
        for d in r.len {
            r[d] = -r[d];
            du[d] = 0.0;
        }
        let its, solved = solvers::gmres(row_start, col, val, r, du, LINEAR_RTOL, LINEAR_RESTART, LINEAR_MAX) { b, z =>
            solvers::ilu0_apply(row_start, col, lu, diag, b, z);
        };
        linear_iterations += its;
        converged = solved;
        for d in sol.len { sol[d] += du[d]; }
        residual(sol, time, r, SRC, BC);
        last_residual = solvers::norm(r);
        iterations++;
        if last_residual < lowest {
            lowest = last_residual;
            lowest_at = iterations;
        }
        if iterations - lowest_at > DEFECT_STALL { converged = false; }
    }
    guard state_sizes() == sizes else { abort("solve_static: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - start;
    return converged;
}

// The largest rate of the explicit operator, |a|/h plus d/h² on each
// axis, for the degree of the method, divided by ct.
fn stable_dt() -> f64 {
    guard ct > 0.0 && method >= 0 else { abort("stable_dt: ct > 0 and a method"); }
    let p = order + 1.0;
    let rate = (2.0 * order + 1.0) * (abs(ax) / hx + abs(ay) / hy) + 2.0 * p * p * p * p * diffusivity * (1.0 / (hx * hx) + 1.0 / (hy * hy));
    return if rate > 0.0 { CFL * ct / rate } else { 1.0e300 };
}

// The weight 1/(ct m) of dof d in a stage: m is the volume of the cell
// (fv) or the diagonal of the mass matrix (dg).
fn stage_weight(d: i64) -> f64 { return 1.0 / (ct * (if method == FV { hx * hy } else { mass(d % per_cell()) })); }

// out = c0 u0 + c1 (v - dt R(v, t)/(ct M)), the stage of SSP-RK3.
fn stage<SRC, BC>(u0: f64[:], v: f64[:], t: f64, dt: f64, c0: f64, c1: f64, out: f64[:]) {
    if method == DG { dg_rhs(t, rhs, SRC, BC); }
    residual(v, t, r, SRC, BC);
    for d in out.len { out[d] = c0 * u0[d] + c1 * (v[d] - dt * r[d] * stage_weight(d)); }
}

// The terms of a stage at time t that the GPU path of the page
// (src/gpu/) takes from the CPU, one value at a time through OUT: b(t)
// times the stage weight of each dof (dg), or the value of each boundary
// face with a condition, 4 for each cell, and then the source times the
// volume of each cell (fv).
fn stage_terms<SRC, BC, OUT>(t: f64) {
    if method == DG {
        dg_rhs(t, rhs, SRC, BC);
        for d in dof_count() { OUT(stage_weight(d) * rhs[d]); }
        return;
    }
    for c in cell_count() {
        for s in 4 { OUT(if face_bc[4 * c + s] >= 0 { face_value(c, s, t, BC) } else { 0.0 }); }
    }
    for c in cell_count() {
        let p = centre(c);
        OUT(SRC(p[0], p[1], t) * hx * hy);
    }
}

// The checks before a transient solve to `end` with the step dt (dt <= 0
// for the stable step), and the matrix of dg. It gives the step and the
// count of steps.
fn transient_start(end: f64, dt: f64) -> f64, i64 {
    guard method >= 0 && ct > 0.0 && end > time && sol.len == dof_count() && u1.len == sol.len else { abort("solve_transient: a method, ct > 0, end after the time of the state, and clear_bcs"); }
    solvers::memory_check("grid ssprk3", col.len, dof_count(), 6, 0, "");
    if method == DG { matrix(); }
    let h = if dt > 0.0 { dt } else { stable_dt() };
    let steps = max(1, ceil((end - time) / h - 1.0e-9) as! i64);
    return (end - time) / steps, steps;
}

// The end of a transient solve of `steps` steps to `end`. False if the
// state is not finite.
fn transient_end(end: f64, steps: i64) -> bool {
    time = end;
    iterations = steps;
    linear_iterations = 0;
    var finite = true;
    for v in sol { if !(abs(v) < 1.0e300) { finite = false; } }
    return finite;
}

// Integrates from the time of the state to `end` with the step dt, or
// with the stable step for dt <= 0. False if the state is not finite.
fn solve_transient<SRC, BC>(end: f64, dt: f64) -> bool {
    let start = clock();
    let sizes = state_sizes();
    let k, steps = transient_start(end, dt);
    let t0 = time;
    for i in steps {
        let t = t0 + i * k;
        stage(sol, sol, t, k, 0.0, 1.0, u1, SRC, BC);
        stage(sol, u1, t + k, k, 0.75, 0.25, u2, SRC, BC);
        stage(sol, u2, t + 0.5 * k, k, 1.0 / 3.0, 2.0 / 3.0, sol, SRC, BC);
    }
    let finite = transient_end(end, steps);
    guard state_sizes() == sizes else { abort("solve_transient: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - start;
    return finite;
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[15] {
    return [face_bc.len, selected.len, sol.len, row_start.len, col.len, val.len, lu.len, diag.len, grad_u.len, phi.len,
            rhs.len, r.len, du.len, u1.len, u2.len];
}
