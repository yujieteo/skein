import std;
namespace codegen;

fn solvers_known(solver: u8[:], preconditioner: u8[:]) -> bool {
    return (solver == "gmres" && preconditioner == "ilu0") || (solver == "cg" && preconditioner == "amg");
}
