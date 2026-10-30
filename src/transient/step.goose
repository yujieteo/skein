// Implicit time stepping of M du/dt + K u = 0 with fixed values on some
// dofs. M and K are sparse (CSR) with the same pattern, and they hold all
// dofs, with the fixed dofs too.
//
// Each step solves A u_new = M w + K v by CG (runtime::cg) with a Jacobi
// preconditioner, where A = cm M + ck K:
//
//     scheme           cm          ck     w                      v
//     backward_euler   1/dt        1      u/dt                   0
//     crank_nicolson   1/dt        1/2    u/dt                   -u/2
//     bdf2             3/(2 dt)    1      (4 u - u_old)/(2 dt)   0
//
// BDF2 needs two old states. Its first step in a run is a backward Euler
// step. Its local error is O(dt²), so the global error stays O(dt²).
import std;
import math;
namespace transient;
import runtime from "../runtime";
import solvers from "../solvers";

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var mass: f64[>..] = [];
var stiff: f64[>..] = [];

var temp: f64[>..] = [];            // u, the temperature of each dof
var fixed: u8[>..] = [];            // 1 for each dof that is held
var fixed_value: f64[>..] = [];     // the value of each held dof
var now = 0.0;

var last_steps = 0;
var last_iterations = 0;

let CG_TOLERANCE = 1.0e-12;   // the residual relative to the right-hand side

fn clear_matrices() {
    row_start.clear();
    col.clear();
    mass.clear();
    stiff.clear();
}

// The scheme of the current step: A = cm M + ck K in a_val, and the
// inverse of its diagonal on the free dofs (0 on the fixed dofs).
var a_val: f64[>..] = [];
var inv_diag: f64[>..] = [];

// y = A x on the free dofs, and 0 on the fixed dofs.
fn free_operator(x: f64[:], y: f64[:]) {
    solvers::multiply(row_start, col, a_val, x, y);
    for d in y.len { if fixed[d] != 0 { y[d] = 0.0; } }
}

fn jacobi(r: f64[:], z: f64[:]) {
    for d in r.len { z[d] = inv_diag[d] * r[d]; }
}

// Solves (cm M + ck K) temp = b on the free dofs, with temp = fixed_value
// on the fixed dofs, by runtime::cg. Starts from the current temp. Returns
// false if CG does not converge.
fn solve_step(cm: f64, ck: f64, b: f64[>..]&) -> bool {
    let n = temp.len;
    a_val.clear();
    for e in mass.len { a_val.push(cm * mass[e] + ck * stiff[e]); }
    solvers::zeros(inv_diag, n);
    for d in n {
        if fixed[d] != 0 {
            temp[d] = fixed_value[d];
            b[d] = 0.0;
        } else {
            inv_diag[d] = 1.0 / a_val[solvers::entry(row_start, col, d, d)];
        }
    }
    let converged, it = runtime::cg(b, temp, CG_TOLERANCE, 10 * n, free_operator, jacobi);
    last_iterations += it;
    return converged;
}

// Steps from `now` to `end` with the scheme. See the table at the top.
fn advance(end: f64, dt_max: f64, scheme: i64) -> bool {
    assert(dt_max > 0.0 && scheme >= 0 && scheme <= 2);
    let steps = max(ceil((end - now) / dt_max - 1.0e-9) as i64, 0);
    let n = temp.len;
    var old: f64[>..] = [];
    var w: f64[>..] = [];
    var v: f64[>..] = [];
    var b: f64[>..] = [];
    var mw: f64[>..] = [];
    solvers::zeros(old, n);
    solvers::zeros(w, n);
    solvers::zeros(v, n);
    solvers::zeros(b, n);
    solvers::zeros(mw, n);
    last_steps = 0;
    last_iterations = 0;
    guard steps > 0 else { return true; }
    let start = now;
    let dt = (end - start) / steps;
    for s in steps {
        var cm = 1.0 / dt;
        var ck = 1.0;
        for d in n { w[d] = temp[d] / dt; v[d] = 0.0; }
        if scheme == CRANK_NICOLSON {
            ck = 0.5;
            for d in n { v[d] = -0.5 * temp[d]; }
        }
        if scheme == BDF2 && s > 0 {
            cm = 1.5 / dt;
            for d in n { w[d] = (4.0 * temp[d] - old[d]) / (2.0 * dt); }
        }
        solvers::multiply(row_start, col, mass, w, mw);
        solvers::multiply(row_start, col, stiff, v, b);
        for d in n {
            b[d] += mw[d];
            old[d] = temp[d];
        }
        guard solve_step(cm, ck, b) else { return false; }
        now = start + (s + 1) * dt;
        last_steps++;
    }
    return true;
}
