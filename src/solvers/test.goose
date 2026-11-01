// Tests of the solvers. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
//
// 1. Newton on a manufactured nonlinear problem: the residual must fall
//    quadratically, and the discrete solution must converge to the exact
//    one at second order in h.
// 2. GMRES with ILU(0) against CG with Jacobi on the grid cantilever of
//    src/runtime/test.goose: the same tip deflection in fewer iterations.
// 3. The near null space of each physics (modes.goose), and the AMG levels
//    that a page keeps from one run to the next (warm.goose).
import std;
import math;
import runtime from "../runtime";
import solvers from "../solvers";

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

// --- Newton: -div((1 + u²) grad u) = f on the unit square ---------------------
//
// The exact solution is u = sin(πx) sin(πy), with u = 0 on the boundary,
// so f = -(1 + u²) Δu - 2u |grad u|². The discrete problem is the 5-point
// finite volume scheme on n × n cells, with the conductivity 1 + m² at the
// midpoint value m of the two nodes of each link. The Jacobian is not
// symmetric, so each Newton step needs GMRES.

var grid = 0;           // n: cells along each side
var h = 0.0;            // the size of a cell
var source: f64[>..] = [];
var exact: f64[>..] = [];
var pattern_start: i64[>..] = [];
var pattern_col: i32[>..] = [];
var jacobian: f64[>..] = [];

// The unknown of interior node (i, j), 1 ≤ i, j ≤ n - 1.
fn unknown(i: i64, j: i64) -> i64 { return (i - 1) + (grid - 1) * (j - 1); }

fn interior(i: i64, j: i64) -> bool { return i > 0 && i < grid && j > 0 && j < grid; }

// The 4 neighbours of a node, south, west, east and north: the order of
// their unknowns.
let DI = [0, -1, 1, 0];
let DJ = [-1, 0, 0, 1];

fn setup_grid(n: i64) {
    grid = n;
    h = 1.0 / n;
    source.clear();
    exact.clear();
    pattern_start.clear();
    pattern_col.clear();
    jacobian.clear();
    pattern_start.push(0);
    for j in 1..n {
        for i in 1..n {
            let x = i * h;
            let y = j * h;
            let u = sin(PI * x) * sin(PI * y);
            let ux = PI * cos(PI * x) * sin(PI * y);
            let uy = PI * sin(PI * x) * cos(PI * y);
            exact.push(u);
            source.push((1.0 + u * u) * 2.0 * PI * PI * u - 2.0 * u * (ux * ux + uy * uy));
            for q in 4 {
                if q == 2 {
                    pattern_col.push(unknown(i, j) as i32);
                    jacobian.push(0.0);
                }
                if interior(i + DI[q], j + DJ[q]) {
                    pattern_col.push(unknown(i + DI[q], j + DJ[q]) as i32);
                    jacobian.push(0.0);
                }
            }
            pattern_start.push(pattern_col.len);
        }
    }
}

fn value_at(u: f64[:], i: i64, j: i64) -> f64 {
    return if interior(i, j) { u[unknown(i, j)] } else { 0.0 };
}

fn residual(u: f64[:], f: f64[:]) {
    for j in 1..grid {
        for i in 1..grid {
            let p = unknown(i, j);
            var s = 0.0;
            for q in 4 {
                let v = value_at(u, i + DI[q], j + DJ[q]);
                let m = 0.5 * (u[p] + v);
                s += (1.0 + m * m) * (u[p] - v);
            }
            f[p] = s / (h * h) - source[p];
        }
    }
}

fn jacobian_values(u: f64[:], val: f64[:]) {
    for j in 1..grid {
        for i in 1..grid {
            let p = unknown(i, j);
            var e = pattern_start[p];
            let self_entry = solvers::entry(pattern_start, pattern_col, p, p);
            val[self_entry] = 0.0;
            for q in 4 {
                let ni = i + DI[q];
                let nj = j + DJ[q];
                let v = value_at(u, ni, nj);
                let m = 0.5 * (u[p] + v);
                val[self_entry] += ((1.0 + m * m) + m * (u[p] - v)) / (h * h);
                if interior(ni, nj) {
                    let c = solvers::entry(pattern_start, pattern_col, p, unknown(ni, nj));
                    val[c] = (-(1.0 + m * m) + m * (u[p] - v)) / (h * h);
                }
            }
        }
    }
}

// Solves on n × n cells from u = 0, prints the residual of each step and
// the rate log(r(k+1) / r(k)) / log(r(k) / r(k-1)), and gives the largest
// error at the nodes.
fn newton_on_grid(n: i64) -> f64 {
    setup_grid(n);
    var u: f64[>..] = [];
    solvers::zeros(u, exact.len);
    var history: f64[>..] = [];
    let tol = 1.0e-12 * solvers::norm(source);
    let steps, converged = solvers::newton(pattern_start, pattern_col, jacobian, u, tol, 20,
                                           1.0e-13, history, residual, jacobian_values);
    if !converged { fail(str("Newton on ", n, " x ", n, " does not converge in ", steps, " steps")); }
    print("newton: ", n, " x ", n, " cells, ", u.len, " unknowns, ", steps, " steps");
    print("  step  |F|                     rate");
    var last_rate = 0.0;
    for r, k in history {
        if k >= 2 && history[k - 1] < history[k - 2] {
            let rate = log(r / history[k - 1]) / log(history[k - 1] / history[k - 2]);
            print("  ", k, "     ", r, "  ", rate);
            // A rate counts only before the residual reaches round-off.
            if r > 1.0e3 * tol { last_rate = rate; }
        } else {
            print("  ", k, "     ", r);
        }
    }
    if abs(last_rate - 2.0) > 0.2 { fail(str("Newton on ", n, " x ", n, ": the last rate is ", last_rate, ", not 2")); }
    var error = 0.0;
    for v, p in u { error = max(error, abs(v - exact[p])); }
    print("  largest error against sin(pi x) sin(pi y): ", error);
    return error;
}

fn test_newton() {
    let coarse = newton_on_grid(32);
    let fine = newton_on_grid(64);
    let order = log(coarse / fine) / log(2.0);
    print("newton: order of the error in h ", order, " (want 2)");
    if abs(order - 2.0) > 0.1 { fail(str("Newton: the order of the error is ", order, ", not 2")); }
}

// --- GMRES with ILU(0) against CG with Jacobi on the cantilever --------------

let L = 100.0;
let B = 10.0;
let H = 10.0;
let E = 210000.0;
let NU = 0.3;
let P = 100.0;
let RTOL = 1.0e-10;         // the tolerance of runtime::solve_static

// The tip deflection of runtime::solve_static on this mesh (CG with
// Jacobi, src/runtime/test.goose).
let RUNTIME_TIP = -0.190502;

// Solves the cantilever with the linear solver and the preconditioner of
// the names, and gives the iterations and the tip deflection.
fn solve_cantilever(solver: u8[:], preconditioner: u8[:]) -> i64, f64 {
    runtime::mesh_box(0.0, 0.0, 0.0, L, B, H, 40, 4, 4);
    runtime::material(E, NU);
    assert(runtime::fix("x0"));
    assert(runtime::traction("x1", 0.0, 0.0, -P / (B * H)));
    runtime::linear_solver(solver, preconditioner);
    if !runtime::solve_static() { fail(str("cantilever: ", solver, "+", preconditioner, " does not converge")); }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    print("  ", solver, "+", preconditioner, " ", runtime::cg_iterations(), " iterations, tip uz ", uz, " mm, ",
          runtime::solve_seconds(), " s");
    // RUNTIME_TIP has 6 digits: 5e-7 mm is half of its last digit.
    if abs(uz - RUNTIME_TIP) > 5.0e-7 { fail(str("cantilever: ", solver, " tip uz ", uz, ", want ", RUNTIME_TIP)); }
    return runtime::cg_iterations(), uz;
}

fn test_cantilever() {
    print("cantilever: 40 x 4 x 4 hex27, relative tolerance ", RTOL);
    let cg_its, cg_tip = solve_cantilever("cg", "jacobi");
    let gm_its, gm_tip = solve_cantilever("gmres", "ilu0");
    if gm_its >= cg_its { fail(str("cantilever: GMRES+ILU(0) takes ", gm_its, " iterations, CG+Jacobi ", cg_its)); }
    let amg_its, amg_tip = solve_cantilever("cg", "amg");
    if amg_its * 3 >= cg_its { fail(str("cantilever: CG+AMG takes ", amg_its, " iterations, CG+Jacobi ", cg_its)); }
}

// The memory check (memory.goose): the /proc fields, and an estimate of the
// RSS now and the entries and the vectors of a path.
fn test_memory() {
    let meminfo = "MemTotal:       32762052 kB\nMemAvailable:   1234 kB\n";
    if solvers::proc_kb(meminfo, "MemAvailable:") != 1234.0 || solvers::proc_kb(meminfo, "Swap:") != -1.0 {
        fail("memory: proc_kb of /proc/meminfo");
    }
    let p = solvers::MemoryPath { name: "test", entry_bytes: 10.0, dof_bytes: 4.0, thread_bytes: 16.0 };
    if solvers::path_bytes(p, 1000, 100, 2, 50, 1) != 12000.0 || solvers::path_bytes(p, 1000, 100, 2, 50, 4) != 28400.0 + 4.0 * solvers::THREAD_BYTES {
        fail("memory: path_bytes");
    }
}

// The modes of a scalar field, of a 2D and a 3D solid, and of a thin plate.
fn test_modes() {
    // The 8 corners of a box of sides 1, 1 and t.
    var cube: f64[>..] = [];
    var plate: f64[>..] = [];
    for k in 8 {
        for t in 3 { cube.push(((k >> t) & 1) as f64); }
        plate.push((k & 1) as f64);
        plate.push(((k >> 1) & 1) as f64);
        plate.push(0.05 * ((k >> 2) & 1));
    }
    var fixed: u8[>..] = [];
    for d in 24 { fixed.push(if d < 3 { 1 } else { 0 }); }
    let m3, nb3 = solvers::near_null_space(solvers::ELASTIC_MODES, cube, 3, 3, fixed);
    let mp, nbp = solvers::near_null_space(solvers::ELASTIC_MODES, plate, 3, 3, fixed);
    let m2, nb2 = solvers::near_null_space(solvers::ELASTIC_MODES, cube[..12], 2, 2, fixed[..8]);
    let ms, nbs = solvers::near_null_space(solvers::SCALAR_MODES, cube, 3, 1, fixed[..8]);
    if nb3 != 6 || nbp != 9 || nb2 != 3 || nbs != 1 { fail(str("modes: ", nb3, ", ", nbp, ", ", nb2, " and ", nbs, " modes, want 6, 9, 3 and 1")); }
    if solvers::mode_count(solvers::ELASTIC_MODES, plate, 3, 3) != nbp || solvers::mode_count(solvers::ELASTIC_MODES, cube, 3, 3) != nb3 {
        fail("modes: mode_count differs from near_null_space");
    }
    for v in mp[..3 * nbp] { if v != 0.0 { fail("modes: a mode is not 0 on a fixed dof"); } }
    // u = (-y, x, 0) about the centre (0.5, 0.5, ...): the rotation about z at node 7 (1, 1, 1).
    if mp[nbp * 21 + 3] != -0.5 || mp[nbp * 22 + 3] != 0.5 || mp[nbp * 23 + 3] != 0.0 { fail("modes: the rotation about z"); }
    // The curvature along x of the plate: u_z = x²/2 about the centre.
    if mp[nbp * 23 + 6] != 0.125 { fail("modes: the bending mode of the plate"); }
    if ms[2] != 0.0 || ms[3] != 1.0 { fail("modes: the constant of a scalar field"); }
}

// The 5-point Laplacian on 80 x 80 nodes (3 AMG levels), for test_warm.
var lap_rs: i64[>..] = [0];
var lap_col: i32[>..] = [];
var lap_val: f64[>..] = [];

fn lap_multiply(x: f64[:], y: f64[:]) { solvers::multiply(lap_rs, lap_col, lap_val, x, y); }

fn lap_amg(r: f64[:], z: f64[:]) { solvers::amg_apply(r, z); }

// CG with the AMG levels that amg_reuse reads from the files of amg_keep
// gives the solution and the iterations of CG with the levels of
// amg_setup; the levels do not fit a matrix of another pattern. It needs
// SKEIN_WARM, which ./build.sh test sets to an empty directory.
fn test_warm() {
    if solvers::warm_dir().len == 0 { fail("warm: SKEIN_WARM is not set"); return; }
    let n = 80;
    for j in n {
        for i in n {
            let r = i + n * j;
            if j > 0 { lap_col.push((r - n) as i32); lap_val.push(-1.0); }
            if i > 0 { lap_col.push((r - 1) as i32); lap_val.push(-1.0); }
            lap_col.push(r as i32);
            lap_val.push(4.0);
            if i < n - 1 { lap_col.push((r + 1) as i32); lap_val.push(-1.0); }
            if j < n - 1 { lap_col.push((r + n) as i32); lap_val.push(-1.0); }
            lap_rs.push(lap_col.len);
        }
    }
    var b: f64[>..] = [];
    var x1: f64[>..] = [];
    var x2: f64[>..] = [];
    var free: u8[>..] = [];
    for r in n * n {
        b.push(sin(0.1 * r));
        x1.push(0.0);
        x2.push(0.0);
        free.push(0);
    }
    let key = solvers::warm_key();
    if !solvers::amg_setup(lap_rs, lap_col, lap_val, 1) { fail("warm: amg_setup"); return; }
    if solvers::amg_levels < 3 { fail(str("warm: ", solvers::amg_levels, " AMG levels, want 3 or more")); }
    let ok1, its1 = runtime::cg(b, x1, 1.0e-10, 500, lap_multiply, lap_amg);
    solvers::amg_keep(key, lap_rs, lap_col, free);
    var other: f64[>..] = [];
    for v in lap_val { other.push(if v > 0.0 { 5.0 } else { v }); }
    if !solvers::amg_setup(lap_rs, lap_col, other, 1) { fail("warm: amg_setup of another matrix"); return; }
    if !solvers::amg_reuse(key, lap_rs, lap_col, lap_val, free) { fail("warm: amg_reuse of the levels of amg_keep"); return; }
    let ok2, its2 = runtime::cg(b, x2, 1.0e-10, 500, lap_multiply, lap_amg);
    if !ok1 || !ok2 || its1 != its2 { fail(str("warm: CG takes ", its1, " iterations, and ", its2, " with the reused levels")); }
    for v, i in x1 { if v != x2[i] { fail("warm: the reused levels give another solution"); break; } }
    var moved: i32[>..] = [];
    moved.append(lap_col);
    moved[1] = 2;
    if solvers::amg_reuse(key, lap_rs, moved, lap_val, free) { fail("warm: amg_reuse of another pattern"); }
}

fn main() {
    test_newton();
    test_cantilever();
    test_memory();
    test_modes();
    test_warm();
    if failures > 0 { exit(1); }
    print("solvers: all tests pass");
}
