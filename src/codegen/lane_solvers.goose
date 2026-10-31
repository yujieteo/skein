// The registration point of the solvers lane (src/solvers/) in the code
// generator: the pairs of `solve { linear: <solver>, precond: <precond> }`
// other than cg with jacobi. A known pair goes to runtime::linear_solver,
// and src/runtime/lane_solvers.goose solves with it.
import std;
namespace codegen;

fn solvers_known(solver: u8[:], preconditioner: u8[:]) -> bool {
    return (solver == "gmres" && preconditioner == "ilu0") || (solver == "cg" && preconditioner == "amg");
}
