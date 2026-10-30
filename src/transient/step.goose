// Implicit time stepping of M du/dt + K u = 0 with fixed values on some
// dofs. M and K are sparse (CSR) with the same pattern, and they hold all
// dofs, with the fixed dofs too.
//
// Each step solves A u_new = M w + K v by CG with a Jacobi preconditioner,
// where A = cm M + ck K:
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

fn zeros(n: i64) -> f64[>..] {
    var v: f64[>..] = [];
    for i in n { v.push(0.0); }
    return v;
}

fn dot(a: f64[:], b: f64[:]) -> f64 {
    var s = 0.0;
    for v, i in a { s += v * b[i]; }
    return s;
}

// y = (cm M + ck K) x
fn multiply(cm: f64, ck: f64, x: f64[:], y: f64[:]) {
    for r in y.len {
        var s = 0.0;
        for e in row_start[r]..row_start[r + 1] {
            s += (cm * mass[e] + ck * stiff[e]) * x[col[e]];
        }
        y[r] = s;
    }
}

// Solves (cm M + ck K) temp = b on the free dofs, with temp = fixed_value
// on the fixed dofs. Starts from the current temp. Returns false if CG
// does not converge.
fn solve_step(cm: f64, ck: f64, b: f64[:]) -> bool {
    let n = temp.len;
    var r = zeros(n);
    var z = zeros(n);
    var p = zeros(n);
    var q = zeros(n);
    var inv_diag = zeros(n);
    for d in n {
        if fixed[d] != 0 { temp[d] = fixed_value[d]; }
    }
    multiply(cm, ck, temp, q);
    var bnorm = 0.0;
    for d in n {
        if fixed[d] != 0 { continue; }
        r[d] = b[d] - q[d];
        bnorm += b[d] * b[d];
        for e in row_start[d]..row_start[d + 1] {
            if col[e] == d { inv_diag[d] = 1.0 / (cm * mass[e] + ck * stiff[e]); }
        }
    }
    let target = CG_TOLERANCE * sqrt(bnorm);
    for d in n { z[d] = inv_diag[d] * r[d]; p[d] = z[d]; }
    var rz = dot(r, z);
    var converged = sqrt(dot(r, r)) <= target;
    var it = 0;
    while !converged && it < 10 * n {
        multiply(cm, ck, p, q);
        for d in n { if fixed[d] != 0 { q[d] = 0.0; } }
        let alpha = rz / dot(p, q);
        for d in n {
            temp[d] += alpha * p[d];
            r[d] -= alpha * q[d];
        }
        it++;
        if sqrt(dot(r, r)) <= target { converged = true; break; }
        for d in n { z[d] = inv_diag[d] * r[d]; }
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    last_iterations += it;
    return converged;
}

// Steps from `now` to `end` with the scheme. See the table at the top.
fn advance(end: f64, dt_max: f64, scheme: i64) -> bool {
    assert(dt_max > 0.0 && scheme >= 0 && scheme <= 2);
    let steps = max(ceil((end - now) / dt_max - 1.0e-9) as i64, 0);
    let n = temp.len;
    var old = zeros(n);
    var w = zeros(n);
    var v = zeros(n);
    var b = zeros(n);
    var mw = zeros(n);
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
        multiply(1.0, 0.0, w, mw);
        multiply(0.0, 1.0, v, b);
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
