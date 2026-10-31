// The assembly of the residual and the Jacobian, the boundary conditions,
// and Newton's method with GMRES and ILU(0) (src/solvers/).
//
// A fixed dof d has the residual u[d] - prescribed[d] and a row of the
// identity in the Jacobian. `raw` keeps the assembled residual of each dof
// before that change: on a fixed velocity dof it is the force of the
// boundary on the flow (the reaction), ∮ N_a (σ n)_i.
import std;
import math;
import os;
namespace flow;
import solvers from "../solvers";
import .mesh;
import .cell;

let NEWTON_MAX = 30;            // the most Newton steps
let NEWTON_SUFFICIENT = 1.0e-4; // the Armijo constant of the line search
let NEWTON_HALVINGS = 20;
let LINEAR_RTOL = 1.0e-6;       // GMRES, relative to the norm of the residual
// Reason: on the P2 cavity of 24174 dofs, GMRES(150) took 1122 iterations
// and GMRES(50) took 12627.
let LINEAR_RESTART = 150;
let LINEAR_MAX = 5000;

var fixed: u8[>..] = [];
var prescribed: f64[>..] = [];
var sol: f64[>..] = [];
var raw: f64[>..] = [];
var val: f64[>..] = [];
var lu: f64[>..] = [];
var pre: f64[>..] = [];     // the matrix of the preconditioner
var diag: i64[>..] = [];

var newton_steps = 0;
var linear_iterations = 0;
var last_seconds = 0.0;
var first_residual = 0.0;
var last_residual = 0.0;

// Sets the material and the equations. The values are scaled.
fn fluid(density: f64, viscosity: f64, navier_stokes: bool, stab: i64) {
    guard density > 0.0 && viscosity > 0.0 && stab >= NONE && stab <= SUPG_PSPG_LSIC else { abort("fluid: density > 0, viscosity > 0 and a stabilize value"); }
    guard order == 2 || stab >= SUPG_PSPG else { abort("fluid: order 2 for the velocity, or pressure stabilization"); }
    rho = density;
    mu = viscosity;
    convective = navier_stokes;
    stabilize = stab;
}

fn body_force(fx: f64, fy: f64) {
    force[0] = fx;
    force[1] = fy;
}

// Frees each dof and sets the unknowns to 0. The arrays of the state get
// their sizes here, before the solve.
fn clear_bcs() {
    guard row_start.len == dof_count() + 1 else { abort("clear_bcs: a mesh and its pattern"); }
    let n = dof_count();
    fixed.clear();
    prescribed.clear();
    sol.clear();
    raw.clear();
    for d in n {
        fixed.push(0);
        prescribed.push(0.0);
        sol.push(0.0);
        raw.push(0.0);
    }
    // The matrices and the solve of solve_steady: the basis of GMRES, f, du,
    // trial, r, z and the `where` of ILU(0). Each worker of the assembly
    // copies the mesh and u.
    solvers::memory_check("flow newton gmres ilu0", col.len, n, LINEAR_RESTART + 7, 8 * (coords.len + conn.len + first_dof.len + n), "");
    val.clear();
    lu.clear();
    pre.clear();
    for e in col.len {
        val.push(0.0);
        lu.push(0.0);
        pre.push(0.0);
    }
    diag.clear();
    for d in n { diag.push(0); }
}

// Fixes velocity component comp (0 or 1) on the nodes of the region.
fn fix_velocity(comp: i64, value: f64) {
    guard (comp == 0 || comp == 1) && cell_count() > 0 && fixed.len == dof_count() && region.len == node_count() else { abort("fix_velocity: comp 0 or 1, and clear_bcs on a mesh"); }
    for n in node_count() {
        if region[n] == 0 { continue; }
        let d = vdof(n, comp);
        fixed[d] = 1;
        prescribed[d] = value;
    }
}

// Fixes the pressure on the vertices of the region.
fn fix_pressure(value: f64) {
    guard cell_count() > 0 && fixed.len == dof_count() && region.len == node_count() else { abort("fix_pressure: clear_bcs on a mesh"); }
    for n in vertex_count {
        if region[n] == 0 { continue; }
        fixed[pdof(n)] = 1;
        prescribed[pdof(n)] = value;
    }
}

// The assembled residual of u in r, before the boundary conditions.
fn assemble_residual(u: f64[:], r: f64[:]) {
    guard u.len == dof_count() && r.len == u.len else { abort("assemble_residual: u and r have one entry for each dof"); }
    for d in r.len { r[d] = 0.0; }
    assemble_cells(u, false, r);
}

// F(u): the residual with the boundary conditions.
fn residual(u: f64[:], r: f64[:]) {
    assemble_residual(u, r);
    for d in r.len {
        if fixed[d] != 0 { r[d] = u[d] - prescribed[d]; }
    }
}

// The Jacobian of F at u, on the pattern of mesh.goose.
fn jacobian(u: f64[:], v: f64[:]) {
    guard u.len == dof_count() && v.len == col.len else { abort("jacobian: u has one entry for each dof, v one for each entry of the pattern"); }
    for e in v.len { v[e] = 0.0; }
    assemble_cells(u, true, v);
    for d in u.len {
        if fixed[d] != 0 { v[solvers::entry(row_start, col, d, d)] = 1.0; }
    }
}

// The u and the kind of an assembly, as the argument of the spawn of a
// worker: the residual, or with jac the Jacobian.
struct Assembly { jac: bool, u: f64[] }

// Adds the residual (out has one entry for each dof) or with jac the
// Jacobian (out has one entry for each entry of the pattern) of each cell
// at u into out, in cell order. It runs on the pool of the threads of the
// case (src/solvers/pool.goose) if it has more than one thread.
fn assemble_cells(u: f64[:], jac: bool, out: f64[:]) {
    let k = cell_dof_count();
    let size = if jac { k * k } else { k };
    let workers = min(solvers::case_threads(), (cell_count() + solvers::CELL_BLOCK - 1) / solvers::CELL_BLOCK);
    if workers < 2 {
        var d: f64[>..] = [];
        for c in cell_count() {
            d.clear();
            cell_values(c, u, jac, d);
            add_cell(c, d[..], jac, out);
        }
        return;
    }
    var ids: i64[>..] = [];
    for w in workers { ids.push(thread_spawn(cell_worker, Assembly { jac: jac, u: u })); }
    solvers::gather_cells(cell_count(), workers, size) { lo, hi, d =>
        for c in lo..hi { add_cell(c, d[(c - lo) * size..(c - lo + 1) * size], jac, out); }
    };
    for id in ids { thread_wait(id); }
}

// A worker of the assembly: the values of each cell of each task, one
// cell after the other. Its copy of the mesh is the copy of the globals
// at its spawn.
thread_fn cell_worker(a: Assembly) {
    var data: f64[>..] = [];
    loop {
        let job = qget<solvers::CellJob>();
        guard job.lo >= 0 else { break; };
        data.clear();
        for c in job.lo..job.hi { cell_values(c, a.u, a.jac, data); }
        qput(solvers::CellData { lo: job.lo, bad: -1, data: data[..] });
    }
}

// Appends to d the residual of cell c at u, or with jac its Jacobian, row
// by row: cell_dof_count() or cell_dof_count() squared values.
fn cell_values(c: i64, u: f64[:], jac: bool, d: f64[>..]&) {
    var dofs: i64[15] = [0; 15];
    var ul: f64[15] = [0.0; 15];
    var rl: f64[15] = [0.0; 15];
    var kl: f64[225] = [0.0; 225];
    let k = cell_dof_count();
    cell_dofs(c, dofs);
    for a in k { ul[a] = u[dofs[a]]; }
    if !jac {
        cell_residual(c, ul, rl);
        d.append(rl[..k]);
        return;
    }
    cell_jacobian(c, ul, rl, kl);
    for a in k { d.append(kl[15 * a..15 * a + k]); }
}

// Adds the values of cell c (cell_values) into out. The Jacobian skips
// the rows of the fixed dofs.
fn add_cell(c: i64, d: f64[:], jac: bool, out: f64[:]) {
    var dofs: i64[15] = [0; 15];
    let k = cell_dof_count();
    cell_dofs(c, dofs);
    for a in k {
        let row = dofs[a];
        if !jac {
            out[row] += d[a];
            continue;
        }
        if fixed[row] != 0 { continue; }
        for b in k { out[solvers::entry(row_start, col, row, dofs[b])] += d[k * a + b]; }
    }
}

// Solves F(u) = 0 from u = 0 with the fixed values, to the norm of F at
// most rtol times its norm at the start. The first step from u = 0 is a
// Stokes solve, since the Jacobian of the convection term is 0 there.
// It grows no array of the state (state_sizes).
fn solve_steady(rtol: f64) -> bool {
    guard cell_count() > 0 && fixed.len == dof_count() && rtol > 0.0 else { abort("solve_steady: clear_bcs, and rtol > 0"); }
    let start = clock();
    let sizes = state_sizes();
    for d in sol.len { sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 }; }
    var f: f64[>..] = [];
    var du: f64[>..] = [];
    var trial: f64[>..] = [];
    solvers::zeros(f, sol.len);
    solvers::zeros(du, sol.len);
    solvers::zeros(trial, sol.len);
    residual(sol, f);
    first_residual = solvers::norm(f);
    last_residual = first_residual;
    newton_steps = 0;
    linear_iterations = 0;
    var converged = true;
    while last_residual > rtol * first_residual {
        if newton_steps >= NEWTON_MAX || !newton_step(f, du, trial) {
            converged = false;
            break;
        }
    }
    assemble_residual(sol, raw);
    guard state_sizes() == sizes else { abort("solve_steady: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - start;
    return converged;
}

// One Newton step with a backtracking line search. False if GMRES or the
// line search fails.
fn newton_step(f: f64[:], du: f64[:], trial: f64[:]) -> bool {
    jacobian(sol, val);
    preconditioner(sol);
    if !solvers::ilu0(row_start, col, pre, lu, diag) { return false; }
    for d in f.len {
        f[d] = -f[d];
        du[d] = 0.0;
    }
    solvers::rows_start(row_start, col, val);
    let its, solved = solvers::gmres(row_start, col, val, f, du, LINEAR_RTOL, LINEAR_RESTART, LINEAR_MAX) { r, z =>
        solvers::ilu0_apply(row_start, col, lu, diag, r, z);
    };
    solvers::rows_stop();
    linear_iterations += its;
    if !solved { return false; }
    var lambda = 1.0;
    for h in NEWTON_HALVINGS + 1 {
        for d in sol.len { trial[d] = sol[d] + lambda * du[d]; }
        residual(trial, f);
        let tn = solvers::norm(f);
        if tn <= (1.0 - NEWTON_SUFFICIENT * lambda) * last_residual {
            for d in sol.len { sol[d] = trial[d]; }
            last_residual = tn;
            newton_steps++;
            return true;
        }
        lambda *= 0.5;
    }
    return false;
}

// The matrix of the preconditioner: the Jacobian, plus the PSPG term of
// the pressure rows when the equations have no PSPG. Reason: with the
// dofs of each node together (mesh.goose), the pressure rows of
// Taylor-Hood cells have a zero diagonal, and ILU(0) stops on it.
fn preconditioner(u: f64[:]) {
    guard u.len == dof_count() && pre.len == val.len else { abort("preconditioner: clear_bcs, and u has one entry for each dof"); }
    for e in val.len { pre[e] = val[e]; }
    guard stabilize < SUPG_PSPG else { return; }
    var dofs: i64[15] = [0; 15];
    var ul: f64[15] = [0.0; 15];
    var kp: f64[9] = [0.0; 9];
    let k = cell_dof_count();
    for c in cell_count() {
        cell_dofs(c, dofs);
        for a in k { ul[a] = u[dofs[a]]; }
        pressure_term(c, ul, kp);
        for a in 3 {
            let row = dofs[2 * per_cell + a];
            if fixed[row] != 0 { continue; }
            for b in 3 { pre[solvers::entry(row_start, col, row, dofs[2 * per_cell + b])] += kp[3 * a + b]; }
        }
    }
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[14] {
    return [coords.len, conn.len, on_boundary.len, row_start.len, col.len, region.len, fixed.len, prescribed.len,
            sol.len, raw.len, val.len, lu.len, pre.len, diag.len];
}

// The sum of the reactions of velocity component comp on the region: the
// force of the region on the flow, ∮ (σ n)_comp over it, scaled.
fn region_force(comp: i64) -> f64 {
    guard (comp == 0 || comp == 1) && raw.len == dof_count() && region.len == node_count() else { abort("region_force: comp 0 or 1, and a solve"); }
    var s = 0.0;
    for n in node_count() {
        if region[n] != 0 { s += raw[vdof(n, comp)]; }
    }
    return s;
}
