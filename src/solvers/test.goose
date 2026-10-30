// Tests of the solvers. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
//
// 1. Newton on a manufactured nonlinear problem: the residual must fall
//    quadratically, and the discrete solution must converge to the exact
//    one at second order in h.
// 2. GMRES with ILU(0) against CG with Jacobi on the grid cantilever of
//    src/runtime/test.goose: the same tip deflection in fewer iterations.
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
// GMRES(m) with ILU(0) needs about 220 iterations here. A restart at 50 or
// 200 makes it stall above RTOL, so m is larger than that count.
let RESTART = 300;

// The tip deflection of runtime::solve_static on this mesh (CG with
// Jacobi, src/runtime/test.goose).
let RUNTIME_TIP = -0.190502;

fn test_cantilever() {
    runtime::mesh_box(0.0, 0.0, 0.0, L, B, H, 40, 4, 4);
    runtime::material(E, NU);
    assert(runtime::fix_face(runtime::X, 0.0));
    assert(runtime::traction_face(runtime::X, L, 0.0, 0.0, -P / (B * H)));
    if !runtime::solve_static() { fail("cantilever: runtime CG does not converge"); }
    let runtime_iterations = runtime::cg_iterations();
    let rx, ry, rz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);

    // The assembled matrix, with the supports as identity rows.
    var val: f64[>..] = [];
    for v in runtime::val { val.push(v); }
    var b: f64[>..] = [];
    for v in runtime::load { b.push(v); }
    solvers::constrain(runtime::row_start, runtime::col, val, runtime::fixed, b);
    let n = b.len;

    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    if !solvers::ilu0(runtime::row_start, runtime::col, val, lu, diag) { fail("cantilever: ILU(0) has a zero pivot"); }
    var x: f64[>..] = [];
    solvers::zeros(x, n);
    let gm_its, gm_ok = solvers::gmres(runtime::row_start, runtime::col, val, b, x, RTOL, RESTART, 1000) { r, z =>
        solvers::ilu0_apply(runtime::row_start, runtime::col, lu, diag, r, z);
    };
    if !gm_ok { fail("cantilever: GMRES with ILU(0) does not converge"); }

    for v, d in x { runtime::disp[d] = v; }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    print("cantilever: ", n, " dofs, relative tolerance ", RTOL);
    print("  runtime CG+Jacobi ", runtime_iterations, " iterations, tip uz ", rz, " mm");
    print("  GMRES(", RESTART, ")+ILU(0) ", gm_its, " iterations, tip uz ", uz, " mm");
    if gm_its >= runtime_iterations { fail(str("cantilever: GMRES+ILU(0) takes ", gm_its, " iterations, CG+Jacobi ", runtime_iterations)); }
    // RUNTIME_TIP has 6 digits: 5e-7 mm is half of its last digit.
    if abs(uz - RUNTIME_TIP) > 5.0e-7 { fail(str("cantilever: GMRES tip uz ", uz, ", want ", RUNTIME_TIP)); }
    if abs(rz - RUNTIME_TIP) > 5.0e-7 { fail(str("cantilever: runtime tip uz ", rz, ", want ", RUNTIME_TIP)); }
}

fn main() {
    test_newton();
    test_cantilever();
    if failures > 0 { exit(1); }
    print("solvers: all tests pass");
}
