// The dispatch tables of the code generator: the physics templates of
// `physics ...: <template>` and the procedures of `step ...: <procedure>`.
// The generator of this package handles linear_elastic and static. A lane
// of phase 2 gives the generator of each other entry in its own file
// (lane_*.goose), so that the lanes add code there and do not edit this
// table.
import std;
namespace codegen;
import .lane_elements;
import .lane_mesher;
import .lane_transient;

let TEMPLATES = ["linear_elastic", "heat"];
let LINEAR_ELASTIC = 0;
let HEAT = 1;

let PROCEDURES = ["static", "transient"];
let STATIC = 0;
let TRANSIENT = 1;

fn table_index(table: const u8[:][:], word: u8[:]) -> i64 {
    for s, i in table {
        if s == word { return i; }
    }
    return -1;
}

// Reads the physics item n of template t (not linear_elastic). False if
// the template has no generator.
fn read_template(t: i64, n: i64) -> bool {
    if t == HEAT { return heat_read(n); }
    return false;
}

// Writes the calls of template t (not linear_elastic) that come after the
// mesh in main.
fn emit_template(out: u8[>..]&, t: i64) {
    if t == HEAT { heat_emit(out); }
}

// Reads the step item n of procedure p (not static). False if the
// procedure has no generator.
fn read_procedure(p: i64, n: i64) -> bool {
    if p == TRANSIENT { return transient_read(n); }
    return false;
}

// Writes the function of step si of procedure p (not static).
fn emit_procedure(out: u8[>..]&, p: i64, si: i64, results_dir: u8[:]) {
    if p == TRANSIENT { transient_emit(out, si, results_dir); }
}
