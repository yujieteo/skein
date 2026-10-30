// The registration point of the solvers lane (src/solvers/): the linear
// solvers and preconditioners other than cg with jacobi (solve.goose). It
// solves K x = b on the free dofs, as solve_free does, and gives whether it
// knows the pair of names, whether it converged and the iterations.
import std;
namespace runtime;
import .solve;

fn solvers_solve(solver: u8[:], preconditioner: u8[:], b: f64[:], x: f64[>..]&, tol: f64,
                 max_it: i64) -> bool, bool, i64 {
    return false, false, 0;
}
