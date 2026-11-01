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
// off, with the near null space of the physics (solvers/modes.goose): the
// rigid-body modes of an elastic K (and the bending modes of a thin body),
// the constant of a scalar field. A solve
// with a key > 0 takes the levels of the last run of a page if they fit K
// (solvers/warm.goose), and keeps its levels for the next run.
// Reason: against the constant of each component, the modes take sample 01
// (tet10, 319320 dofs, 4 threads) from 36.5 s to 30.6 s and 181 to 53 CG
// iterations, and a thin plate (tet, 91809 dofs) from 32.0 s to 5.9 s and
// 866 to 55 iterations. The memory estimate stays at or above the peak:
// on sample 01, 3123 MB against 2776 MB with 4 threads, 2646 MB against
// 2429 MB with 1 thread.
fn setup_preconditioner(key: i64) -> bool {
    var free_val: f64[>..] = [];
    free_val.append(val);
    for d in sol.len {
        for e in row_start[d]..row_start[d + 1] {
            if col[e] != d && (fixed[d] != 0 || fixed[col[e]] != 0) { free_val[e] = 0.0; }
        }
    }
    if key > 0 && solvers::amg_reuse(key, row_start, col, free_val, fixed) { return true; }
    let modes, nb = solvers::near_null_space(mode_kind(), coords, dim, ncomp, fixed);
    guard solvers::amg_setup_modes(row_start, col, free_val, ncomp, modes, nb) else { return false; }
    if key > 0 { solvers::amg_keep(key, row_start, col, fixed); }
    return true;
}

fn mode_kind() -> i64 { return if physics == ELASTIC { solvers::ELASTIC_MODES } else { solvers::SCALAR_MODES }; }

// The vectors of the near null space of the physics.
fn mode_count() -> i64 { return solvers::mode_count(mode_kind(), coords, dim, ncomp); }

fn solve_system() -> bool {
    guard physics >= 0 && dofs_sized() else { abort("solve_system: a physics and the arrays of the dofs"); }
    let start = clock();
    // r, z, p, q, the source vector, w of a warm start, and the modes of the
    // AMG with their copy
    assemble(if physics == ELASTIC { "simplex elastic cg amg" } else { "simplex cg amg" }, 6 + 2 * mode_count());
    hold_dead_nodes();
    for d in sol.len { sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 }; }
    let key = solvers::warm_key();
    let warm = !direct && solvers::warm_start(key, sol, fixed);
    let converged = if rad_count() > 0 { radiation_solve(start, key, warm) } else { linear_solve(start, key, warm) };
    if converged { solvers::warm_keep(key, sol); }
    return converged;
}

// Scales the free dofs of sol, the last solution, by the factor with the
// least energy, and r, the residual of sol, with it. z is scratch. Returns
// the square of the residual of the start without the last solution, so
// the target of CG is the target of a cold run.
fn warm_scale(r: f64[:], z: f64[:]) -> f64 {
    let n = sol.len;
    var w = zeros(n);
    for d in n { if fixed[d] == 0 { w[d] = sol[d]; } }
    multiply(w, z);
    var wr = 0.0;
    var wkw = 0.0;
    var cold = 0.0;
    for d in n {
        if fixed[d] != 0 { continue; }
        wr += w[d] * r[d];
        wkw += w[d] * z[d];
        cold += (r[d] + z[d]) * (r[d] + z[d]);
    }
    let alpha = if wkw > 0.0 { (wr + wkw) / wkw } else { 0.0 };
    for d in n {
        if fixed[d] != 0 { continue; }
        sol[d] = alpha * w[d];
        r[d] -= (alpha - 1.0) * z[d];
    }
    return cold;
}

// CG on K du = f - K u from the u in sol, for the matrix and the loads of
// the last assemble. With a key >= 0, it builds the AMG levels of K
// (setup_preconditioner); with -1, it keeps the levels of the last call.
// With `warm`, the free dofs of sol hold the solution of the last run of a
// page, and CG starts from the multiple of it with the least energy. On
// sample 01, a change of the load from 2 kN to 3 kN then needs 9
// iterations, not 53.
fn cg_solve(start: f64, key: i64, warm: bool) -> bool {
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
    let r2 = if warm { warm_scale(r, z) } else { dot(r, r) };
    guard key < 0 || setup_preconditioner(key) else {
        solvers::rows_stop();
        last_iterations = 0;
        last_seconds = clock() - start;
        return false;
    }
    let target = cg_tolerance * sqrt(max(f2, r2));
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

// The solve of K du = f - K u from the u in sol: the direct solve of a
// lattice (plan_direct), or cg_solve.
fn linear_solve(start: f64, key: i64, warm: bool) -> bool {
    return if direct { direct_static(start, key) } else { cg_solve(start, key, warm) };
}

// The direct solve of K du = f - K u (solvers/nested.goose) on the plan
// of plan_direct, as cg_solve. If the factor finds a pivot with no
// stiffness, CG with AMG takes the solve, as before the direct solve.
// Reason: a nearly incompressible solid has a pivot that the factor rejects,
// and CG with AMG solved it before the direct solve.
// With a key > 0, a page run takes the factor of the last run if K is the
// same (solvers/warm.goose), and keeps its own.
fn direct_static(start: f64, key: i64) -> bool {
    let n = sol.len;
    var r = zeros(n);
    var du = zeros(n);
    let sizes = state_sizes();
    solvers::multiply(row_start, col, val, sol, du);
    for d in n { r[d] = if fixed[d] != 0 { 0.0 } else { load[d] + assembled_source[d] - du[d] }; }
    var done = key > 0 && solvers::nested_reuse(key, row_start, col, val, fixed, r, du);
    if !done {
        done = solvers::nested_factor(row_start, col, val, fixed);
        if done {
            if key > 0 { solvers::nested_keep(key, row_start, col, val, fixed); }
            solvers::nested_solve(r, du, fixed);
        } else {
            write_stderr(str("elements: the direct solve has a pivot with no stiffness at dof ", solvers::nested_singular(), ", CG with AMG takes the solve\n"));
            solvers::nested_free();
            solvers::memory_check(direct_path[..], col.len, dof_total(), direct_vectors, worker_mesh_bytes(), "");
            direct = false;
            return cg_solve(start, key, false);
        }
    }
    if done { for d in n { if fixed[d] == 0 { sol[d] += du[d]; } } }
    guard state_sizes() == sizes else { abort("direct_static: a solve grows no array of the state (rule 3)"); }
    last_iterations = 1;
    last_seconds = clock() - start;
    return done;
}

// The Picard iteration of radiation.goose: each step solves
// (K + D) T' = f - R(T) + D T with CG, from T. The AMG levels of the first
// step stay the preconditioner of the next steps: only D changes. At the end `load` holds
// f - R(T), so the reactions and the check of the equilibrium hold the
// radiated heat.
fn radiation_solve(start: f64, key: i64, warm: bool) -> bool {
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
        // Reason: a step after the first keeps the AMG levels of step 0. On
        // radiation_annulus_tri, the CG iterations stay at 210, the time
        // falls from 1.46 s to 1.18 s, and the result agrees within 1e-10.
        guard linear_solve(start, if rad_steps == 0 { key } else { -1 }, warm && rad_steps == 0) else { return false; }
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
