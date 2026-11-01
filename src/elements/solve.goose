// CG with a smoothed aggregation AMG preconditioner (src/solvers/amg.goose)
// on the free dofs of the mesh of this package. The solver interface of src/runtime/solve.goose and the
// solvers of src/solvers/ work on the mesh of src/runtime/fem.goose.
//
// The fixed dofs take their prescribed values first. CG then solves
// K du = f - K u0 on the free dofs: it holds the residual and the search
// direction at 0 on the fixed dofs. The AMG levels are built from K with
// the rows and the columns of the fixed dofs cut off, and the nodes of the
// mesh as their blocks.
import std;
import math;
import os;
namespace elements;
import .mesh;
import .assemble;
import .stages;
import .radiation;
import solvers from "../solvers";

var last_iterations = 0;
var last_seconds = 0.0;

let CG_TOLERANCE = 1.0e-12;   // the residual relative to the larger of the load and the first residual
var cg_tolerance = CG_TOLERANCE; // of the next cg_solve: a Newton iteration needs less

// y = K x, on the workers of solvers::rows_start during a solve.
fn multiply(x: f64[:], y: f64[:]) {
    guard x.len == y.len && row_start.len == y.len + 1 && col.len == val.len else { abort("multiply: x and y have one entry for each row of the matrix"); }
    solvers::rows_multiply(row_start, col, val, x, y);
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

// z = M⁻¹ r, held at 0 on the fixed dofs.
fn precondition(r: f64[:], z: f64[:]) {
    solvers::amg_apply(r, z);
    for d in z.len { if fixed[d] != 0 { z[d] = 0.0; } }
}

// The AMG levels of K with the rows and the columns of the fixed dofs cut
// off. With `rigid`, the near null space of an elastic K is the rigid-body
// modes, 0 on the fixed dofs: LOBPCG (eigen.goose) needs it, 199 steps
// fall to 36 on a cantilever. A static solve takes the constant of each
// component: on sample 04 the rigid-body modes halve the CG iterations but
// double the time of the setup.
fn setup_preconditioner(rigid: bool) -> bool {
    var free_val: f64[>..] = [];
    free_val.append(val);
    for d in sol.len {
        for e in row_start[d]..row_start[d + 1] {
            if col[e] != d && (fixed[d] != 0 || fixed[col[e]] != 0) { free_val[e] = 0.0; }
        }
    }
    guard rigid && physics == ELASTIC else { return solvers::amg_setup(row_start, col, free_val, ncomp); }
    let nb = if dim == 3 { 6 } else { 3 };
    var modes: f64[>..] = [];
    for d in sol.len {
        let n = d / ncomp;
        let c = d % ncomp;
        let x = node_x(n, 0);
        let y = node_x(n, 1);
        let z = node_x(n, 2);
        // the translations, then the rotations about z, x and y
        var m: f64[6] = [0.0; 6];
        m[c] = 1.0;
        if c == 0 { m[dim] = -y; }
        if c == 1 { m[dim] = x; }
        if dim == 3 && c == 1 { m[4] = -z; }
        if dim == 3 && c == 2 { m[4] = y; }
        if dim == 3 && c == 0 { m[5] = z; }
        if dim == 3 && c == 2 { m[5] = -x; }
        for k in nb { modes.push(if fixed[d] != 0 { 0.0 } else { m[k] }); }
    }
    return solvers::amg_setup_modes(row_start, col, free_val, ncomp, modes, nb);
}

fn solve_system() -> bool {
    guard physics >= 0 && dofs_sized() else { abort("solve_system: a physics and the arrays of the dofs"); }
    let start = clock();
    // r, z, p, q and the source vector
    assemble("simplex cg amg", 5);
    hold_dead_nodes();
    for d in sol.len { sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 }; }
    if rad_count() > 0 { return radiation_solve(start); }
    return cg_solve(start);
}

// CG on K du = f - K u from the u in sol, for the matrix and the loads of
// the last assemble.
fn cg_solve(start: f64) -> bool {
    let n = sol.len;
    var r = zeros(n);
    var z = zeros(n);
    var p = zeros(n);
    var q = zeros(n);
    let sizes = state_sizes();
    solvers::rows_start(row_start, col, val);
    multiply(sol, q);
    var f2 = 0.0;
    for d in n {
        if fixed[d] != 0 { continue; }
        let f = load[d] + assembled_source[d];
        f2 += f * f;
        r[d] = f - q[d];
    }
    guard setup_preconditioner(false) else {
        solvers::rows_stop();
        last_iterations = 0;
        last_seconds = clock() - start;
        return false;
    }
    let target = cg_tolerance * sqrt(max(f2, dot(r, r)));
    precondition(r, z);
    for d in n { p[d] = z[d]; }
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
        precondition(r, z);
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    solvers::rows_stop();
    guard state_sizes() == sizes else { abort("solve_system: a solve grows no array of the state (rule 3)"); }
    last_iterations = it;
    last_seconds = clock() - start;
    return converged;
}

// The Picard iteration of radiation.goose: each step solves
// (K + D) T' = f - R(T) + D T with CG, from T. At the end `load` holds
// f - R(T), so the reactions and the check of the equilibrium hold the
// radiated heat.
fn radiation_solve(start: f64) -> bool {
    exchange_matrix();
    let n = sol.len;
    var base: f64[>..] = [];
    var k: f64[>..] = [];
    var r = zeros(n);
    var dg = zeros(n);
    var last = zeros(n);
    base.append(load);
    k.append(val);
    var converged = false;
    var its = 0;
    rad_steps = 0;
    while !converged && rad_steps < RAD_STEPS {
        radiated(r, dg);
        for d in n {
            val[entry(d, d)] = k[entry(d, d)] + dg[d];
            load[d] = base[d] - r[d] + dg[d] * sol[d];
            last[d] = sol[d];
        }
        guard cg_solve(start) else { return false; }
        its += last_iterations;
        rad_steps++;
        var change = 0.0;
        var top = 0.0;
        for d in n {
            change = max(change, abs(sol[d] - last[d]));
            top = max(top, abs(sol[d]));
        }
        converged = change <= RAD_TOLERANCE * top;
    }
    radiated(r, dg);
    for d in n { load[d] = base[d] - r[d]; }
    for e in val.len { val[e] = k[e]; }
    last_iterations = its;
    return converged;
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[11] {
    return [coords.len, conn.len, radial.len, fixed.len, prescribed.len, load.len, sol.len, row_start.len,
            col.len, val.len, assembled_source.len];
}
