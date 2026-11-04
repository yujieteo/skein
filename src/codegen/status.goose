import std;
namespace codegen;

let SOLVER_CALLS = ["runtime::solve_static(", "elements::solve(", "shells::solve("];

fn unchecked_status(program: u8[:]) -> i64 {
    for line, i in split(program, '\n') {
        let t = trim(line);
        for call in SOLVER_CALLS {
            if starts_with(t, call) { return i + 1; }
        }
    }
    return 0;
}
