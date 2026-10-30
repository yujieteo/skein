// The entry of the package solvers: the sparse linear and nonlinear
// solvers. A program imports it with
// `import solvers from "<path to src>/solvers";`.
//
// A sparse matrix is CSR in three arrays, as in the runtime
// (docs/implementation.md §6): `row_start` (n + 1 entries), `col` (the
// columns of each row in increasing order) and `val`. The solvers keep no
// state between calls.
namespace solvers;
import .sparse;
import .ilu;
import .krylov;
import .newton;
