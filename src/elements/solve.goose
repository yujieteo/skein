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

let CG_TOLERANCE = 1.0e-12;
let CG_ITERATIONS_PER_DOF = 10;
var cg_tolerance = CG_TOLERANCE;

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

fn precondition(r: f64[:], z: f64[:]) {
    solvers::amg_apply(r, z);
    for d in z.len { if fixed[d] != 0 { z[d] = 0.0; } }
}

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

fn mode_count() -> i64 { return solvers::mode_count(mode_kind(), coords, dim, ncomp); }

fn solve_system() -> bool {
    guard physics >= 0 && dofs_sized() else { abort("solve_system: a physics and the arrays of the dofs"); }
    let start = clock();
    assemble(if physics == ELASTIC { "simplex elastic cg amg" } else { "simplex cg amg" }, 6 + 2 * mode_count());
    hold_dead_nodes();
    let nonlinear = physics == HEAT && (curve_k.len > 0 || amb_cell.len > 0);
    for d in sol.len { sol[d] = if fixed[d] != 0 { prescribed[d] } else if nonlinear { initial_value } else { 0.0 }; }
    if nonlinear && curve_k.len > 0 { reassemble(); }
    let key = solvers::warm_key();
    let warm = !direct && !nonlinear && solvers::warm_start(key, sol, fixed);
    var converged = false;
    if nonlinear {
        converged = nonlinear_heat_solve(start, key);
    } else {
        converged = if rad_count() > 0 { radiation_solve(start, key, warm) } else { linear_solve(start, key, warm) };
    }
    if converged {
        solvers::warm_keep(key, sol);
        state_valid = true;
    }
    return converged;
}

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
    while !converged && it < CG_ITERATIONS_PER_DOF * n {
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

fn linear_solve(start: f64, key: i64, warm: bool) -> bool {
    return if direct { direct_static(start, key) } else { cg_solve(start, key, warm) };
}

fn direct_static(start: f64, key: i64) -> bool {
    let n = sol.len;
    var r = zeros(n);
    var du = zeros(n);
    let sizes = state_sizes();
    solvers::multiply(row_start, col, val, sol, du);
    for d in n { r[d] = if fixed[d] != 0 { 0.0 } else { load[d] + assembled_source[d] - du[d] }; }
    var done = key > 0 && solvers::nested_reuse(key, row_start, col, val, fixed, r, du);
    if !done && keep_factor && factor_ready {
        solvers::nested_solve(r, du, fixed);
        done = true;
    }
    if !done {
        done = solvers::nested_factor(row_start, col, val, fixed);
        factor_ready = done;
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

let NL_STEPS = 200;
var nl_tolerance = 1.0e-10;
var nl_steps = 0;

fn nonlinear_heat_solve(start: f64, key: i64) -> bool {
    let n = sol.len;
    var base: f64[>..] = [];
    base.append(load);
    var r = zeros(n);
    var last = zeros(n);
    var tangent = zeros(val.len);
    var tangent_sol = zeros(n);
    var converged = false;
    var its = 0;
    nl_steps = 0;
    while !converged && nl_steps < NL_STEPS {
        if nl_steps > 0 && curve_k.len > 0 { reassemble(); }
        for d in n { r[d] = 0.0; }
        for e in val.len { tangent[e] = 0.0; }
        ambient_residual(r);
        ambient_tangent(tangent);
        solvers::multiply(row_start, col, tangent, sol, tangent_sol);
        for e in val.len { val[e] += tangent[e]; }
        for d in n {
            load[d] = base[d] - r[d] + tangent_sol[d];
            last[d] = sol[d];
        }
        let done = linear_solve(start, if nl_steps == 0 { key } else { -1 }, false);
        for e in val.len { val[e] -= tangent[e]; }
        guard done else { return false; }
        its += last_iterations;
        nl_steps++;
        var change = 0.0;
        var top = 0.0;
        for d in n {
            change = max(change, abs(sol[d] - last[d]));
            top = max(top, abs(sol[d]));
        }
        converged = change <= nl_tolerance * top;
    }
    if curve_k.len > 0 { reassemble(); }
    for d in n { r[d] = 0.0; }
    ambient_residual(r);
    for d in n { load[d] = base[d] - r[d]; }
    last_iterations = its;
    return converged;
}

fn state_sizes() -> i64[11] {
    return [coords.len, conn.len, radial.len, fixed.len, prescribed.len, load.len, sol.len, row_start.len,
            col.len, val.len, assembled_source.len];
}
