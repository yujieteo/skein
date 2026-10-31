// Rule 7 of the Power of Ten (docs/tradeoffs.md §6): the plan checks the
// status that each solver call returns. A generated statement that is only
// a solver call drops its status, and skein run stops with an internal
// fault.
import std;
namespace codegen;

// The solver calls of the runtimes, each of which returns a status.
let SOLVER_CALLS = ["runtime::solve_static(", "elements::solve(", "shells::solve("];

// The number, from 1, of the first line of `program` that is a statement
// with only a solver call, or 0 if there is none.
fn unchecked_status(program: u8[:]) -> i64 {
    for line, i in split(program, '\n') {
        let t = trim(line);
        for call in SOLVER_CALLS {
            if starts_with(t, call) { return i + 1; }
        }
    }
    return 0;
}
