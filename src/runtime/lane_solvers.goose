// The registration point of the solvers lane (src/solvers/): the linear
// solvers and preconditioners other than cg with jacobi (solve.goose). It
// solves K x = b on the free dofs, as solve_free does, and gives whether it
// knows the pair of names, whether it converged and the iterations.
//
// "cg" with "amg" is CG with the smoothed aggregation AMG of
// solvers/amg.goose, built on K with the couplings of the fixed dofs
// removed, and with the rigid-body modes of a solid (rigid_setup).
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

fn amg_solve(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, bool, i64 {
    var a: f64[>..] = [];
    a.append(val);
    for d in fixed.len {
        for e in row_start[d]..row_start[d + 1] {
            if col[e] != d && (fixed[d] != 0 || fixed[col[e]] != 0) { a[e] = 0.0; }
        }
    }
    guard rigid_setup(a) else {
        write_stderr("runtime: the coarsest level of AMG is singular\n");
        return true, false, 0;
    }
    let converged, its = cg(b, x, tol, max_it, free_operator, amg_precondition);
    return true, converged, its;
}

// The AMG levels of the matrix a, with the rigid-body modes as the near
// null space for a solid (3 dofs for each node), 0 on the fixed dofs.
// Reason: on warpage_thick_copper (44541 dofs, 1 cell through 20 um of
// copper) the constant of each component needs 571 CG iterations and 21 s,
// and the rigid-body modes need 88 iterations and 6.2 s.
fn rigid_setup(a: f64[:]) -> bool {
    guard dofs_per_node == 3 else { return solvers::amg_setup(row_start, col, a, dofs_per_node); }
    var modes: f64[>..] = [];
    for d in fixed.len {
        let n = d / 3;
        let c = d % 3;
        let x = coords[3 * n];
        let y = coords[3 * n + 1];
        let z = coords[3 * n + 2];
        // the translations, then the rotations about z, x and y
        var m: f64[6] = [0.0; 6];
        m[c] = 1.0;
        if c == 0 { m[3] = -y; m[5] = z; }
        if c == 1 { m[3] = x; m[4] = -z; }
        if c == 2 { m[4] = y; m[5] = -x; }
        for k in 6 { modes.push(if fixed[d] != 0 { 0.0 } else { m[k] }); }
    }
    return solvers::amg_setup_modes(row_start, col, a, 3, modes, 6);
}
