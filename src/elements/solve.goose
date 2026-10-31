// CG with a Jacobi preconditioner on the free dofs of the mesh of this
// package. The solver interface of src/runtime/solve.goose and the
// solvers of src/solvers/ work on the mesh of src/runtime/fem.goose.
//
// The fixed dofs take their prescribed values first. CG then solves
// K du = f - K u0 on the free dofs: it holds the residual and the search
// direction at 0 on the fixed dofs.
import std;
import math;
import os;
namespace elements;
import .mesh;
import .assemble;

var last_iterations = 0;
var last_seconds = 0.0;

let CG_TOLERANCE = 1.0e-12;   // the residual relative to the larger of the load and the first residual

// y = K x
fn multiply(x: f64[:], y: f64[:]) {
    guard x.len == y.len && row_start.len == y.len + 1 && col.len == val.len else { abort("multiply: x and y have one entry for each row of the matrix"); }
    for r in y.len {
        var s = 0.0;
        for e in row_start[r]..row_start[r + 1] { s += val[e] * x[col[e]]; }
        y[r] = s;
    }
}

fn dot(a: f64[:], b: f64[:]) -> f64 {
    guard a.len == b.len else { abort("dot: a.len == b.len"); }
    var s = 0.0;
    for v, i in a { s += v * b[i]; }
    return s;
}

fn zeros(n: i64) -> f64[>..] {
    guard n >= 0 else { abort("zeros: n >= 0"); }
    var v: f64[>..] = [];
    for i in n { v.push(0.0); }
    return v;
}

fn solve_system() -> bool {
    guard physics >= 0 && dofs_sized() else { abort("solve_system: a physics and the arrays of the dofs"); }
    let start = clock();
    assemble();
    let n = sol.len;
    var r = zeros(n);
    var z = zeros(n);
    var p = zeros(n);
    var q = zeros(n);
    var inv_diag = zeros(n);
    let sizes = state_sizes();
    for d in n { sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 }; }
    multiply(sol, q);
    var f2 = 0.0;
    for d in n {
        if fixed[d] != 0 { continue; }
        let f = load[d] + assembled_source[d];
        f2 += f * f;
        r[d] = f - q[d];
        inv_diag[d] = 1.0 / val[entry(d, d)];
    }
    let target = CG_TOLERANCE * sqrt(max(f2, dot(r, r)));
    for d in n {
        z[d] = inv_diag[d] * r[d];
        p[d] = z[d];
    }
    var rz = dot(r, z);
    var converged = sqrt(dot(r, r)) <= target;
    var it = 0;
    while !converged && it < 10 * n {
        multiply(p, q);
        for d in n { if fixed[d] != 0 { q[d] = 0.0; } }
        let alpha = rz / dot(p, q);
        for d in n {
            sol[d] += alpha * p[d];
            r[d] -= alpha * q[d];
        }
        it++;
        if sqrt(dot(r, r)) <= target {
            converged = true;
            break;
        }
        for d in n { z[d] = inv_diag[d] * r[d]; }
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    guard state_sizes() == sizes else { abort("solve_system: a solve grows no array of the state (rule 3)"); }
    last_iterations = it;
    last_seconds = clock() - start;
    return converged;
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[11] {
    return [coords.len, conn.len, radial.len, fixed.len, prescribed.len, load.len, sol.len, row_start.len,
            col.len, val.len, assembled_source.len];
}
