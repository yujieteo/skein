// The registration point of the solvers lane (src/solvers/): the linear
// solvers and preconditioners other than cg with jacobi (solve.goose). It
// solves K x = b on the free dofs, as solve_free does, and gives whether it
// knows the pair of names, whether it converged and the iterations.
//
// "gmres" with "ilu0" is GMRES(m) with ILU(0) on K with the fixed dofs as
// identity rows and columns (solvers::constrain). It is for small models
// only: the restart is fixed at GMRES_RESTART = 300, and each solve holds
// about 301 basis vectors of n doubles. A finer mesh needs more iterations,
// and the restart can then make GMRES stall.
import std;
namespace runtime;
import .solve;
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
