// The registration point of the solvers lane (src/solvers/): the linear
// solvers and preconditioners other than cg with jacobi (solve.goose). It
// solves K x = b on the free dofs, as solve_free does, and gives whether it
// knows the pair of names, whether it converged and the iterations.
//
// "cg" with "amg" is CG with the smoothed aggregation AMG of
// solvers/amg.goose, built on K with the couplings of the fixed dofs
// removed, and with the near null space of the physics (rigid_setup).
//
// "gmres" with "ilu0" is GMRES(m) with ILU(0) on K with the fixed dofs as
// identity rows and columns (solvers::constrain). It is for small models
// only: the restart is fixed at GMRES_RESTART = 300, and each solve holds
// about 301 basis vectors of n doubles. A finer mesh needs more iterations,
// and the restart can then make GMRES stall.
import std;
namespace runtime;
import .solve;
import .fem;
import solvers from "../solvers";

// m of GMRES(m). On the grid cantilever (src/solvers/test.goose) with
// ILU(0), the residual estimate of GMRES is at the tolerance after 217
// iterations, and the true residual is then 1.9 times larger. Restart
// cycles of one iteration each then take it under the tolerance at 385
// iterations. Each m from 250 to 390 gives these 385 iterations. A
// restart at 100 or 200 stalls.
let GMRES_RESTART = 300;

fn solvers_solve(solver: u8[:], preconditioner: u8[:], b: f64[:], x: f64[>..]&, tol: f64,
                 max_it: i64) -> bool, bool, i64 {
    if solver == "cg" && preconditioner == "amg" { return amg_solve(b, x, tol, max_it); }
    guard solver == "gmres" && preconditioner == "ilu0" else { return false, false, 0; }
    var a: f64[>..] = [];
    for v in val { a.push(v); }
    var rhs: f64[>..] = [];
    for v in b { rhs.push(v); }
    solvers::constrain(row_start, col, a, fixed, rhs);
    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    guard solvers::ilu0(row_start, col, a, lu, diag) else {
        write_stderr("runtime: ILU(0) has a zero pivot\n");
        return true, false, 0;
    }
    let its, converged = solvers::gmres(row_start, col, a, rhs, x, tol, GMRES_RESTART, max_it) { r, z =>
        solvers::ilu0_apply(row_start, col, lu, diag, r, z);
    };
    if !converged {
        write_stderr(str("runtime: gmres with ilu0 did not converge in ", its,
                         " iterations with a restart of ", GMRES_RESTART, "\n"));
    }
    return true, converged, its;
}

// z = M⁻¹ r of the AMG levels, and 0 on the fixed dofs.
fn amg_precondition(r: f64[:], z: f64[:]) {
    solvers::amg_apply(r, z);
    for d in z.len { if fixed[d] != 0 { z[d] = 0.0; } }
}

// A page run takes the AMG levels and the solution of the last run if they
// fit (solvers/warm.goose), and keeps its own for the next run. CG starts
// from the multiple of the last solution with the least energy: a change
// of the size of the loads alone then needs no iterations.
fn amg_solve(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, bool, i64 {
    var a: f64[>..] = [];
    a.append(val);
    for d in fixed.len {
        for e in row_start[d]..row_start[d + 1] {
            if col[e] != d && (fixed[d] != 0 || fixed[col[e]] != 0) { a[e] = 0.0; }
        }
    }
    let key = solvers::warm_key();
    if !solvers::amg_reuse(key, row_start, col, a, fixed) {
        guard rigid_setup(a) else {
            write_stderr("runtime: the coarsest level of AMG is singular\n");
            return true, false, 0;
        }
        solvers::amg_keep(key, row_start, col, fixed);
    }
    if solvers::warm_start(key, x, fixed) {
        var ax = zeros(x.len);
        free_operator(x, ax);
        let xax = dot(x, ax);
        let alpha = if xax > 0.0 { dot(x, b) / xax } else { 0.0 };
        for d in x.len { if fixed[d] == 0 { x[d] *= alpha; } }
    }
    let converged, its = cg(b, x, tol, max_it, free_operator, amg_precondition);
    if converged { solvers::warm_keep(key, x); }
    return true, converged, its;
}

// The AMG levels of the matrix a, with the near null space of the physics
// (solvers/modes.goose): the rigid-body modes of a solid (3 dofs for each
// node) and the bending modes of a thin solid, the constant of each
// component otherwise.
// Reason: on warpage_thick_copper (44541 dofs, 1 cell through 20 um of
// copper) the constant of each component needs 571 CG iterations and 21 s,
// and the rigid-body modes need 88 iterations and 6.2 s.
fn rigid_setup(a: f64[:]) -> bool {
    let modes, nb = solvers::near_null_space(mode_kind(), coords, 3, dofs_per_node, fixed);
    return solvers::amg_setup_modes(row_start, col, a, dofs_per_node, modes, nb);
}

// The near null space of the field: a solid has 3 dofs for each node.
fn mode_kind() -> i64 { return if dofs_per_node == 3 { solvers::ELASTIC_MODES } else { solvers::SCALAR_MODES }; }

// The memory check of a solve with solver_name and preconditioner_name,
// after the pattern of `entries` entries (solvers::memory_check). The
// vectors: r, z, p, q, b and inv_diag of cg, with Ax of a warm start and
// the modes of the AMG and their copy, and the basis of gmres with rhs, r, z, b, inv_diag and the
// i64 diag and `where` of ILU(0).
fn check_memory(entries: i64) {
    let modes = if preconditioner_name == "amg" { solvers::mode_count(mode_kind(), coords, 3, dofs_per_node) } else { 0 };
    let vectors = if solver_name == "gmres" { GMRES_RESTART + 8 } else { 7 + 2 * modes };
    // The coarse levels of 9 modes hold more for each entry than those of 6.
    let bending = if modes > 6 { " bending" } else { "" };
    let fix = if preconditioner_name == "jacobi" { "" } else { ", or write solve { linear: cg, precond: jacobi } in the step: it needs less memory, but it can need many more iterations" };
    // The mesh that each cell_worker copies at its spawn.
    let mesh = 8 * (coords.len + cell_start.len + cell_nodes.len + cell_young.len + cell_poisson.len) + cell_type.len;
    solvers::memory_check(str("box ", solver_name, " ", preconditioner_name, bending), entries, dof_total(), vectors, mesh, fix);
}
