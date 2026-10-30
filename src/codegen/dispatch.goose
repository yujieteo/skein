// The dispatch tables of the code generator: the physics templates of
// `physics ...: <template>` and the procedures of `step ...: <procedure>`.
// The generator of this package handles linear_elastic and static (CORE).
// Each lane of phase 2 has one registration file (lane_*.goose). The file
// says which template and procedure names the lane owns, and it holds
// their generators, so that a lane adds entries in its own file and does
// not edit this table.
import std;
namespace codegen;
import .lane_elements;
import .lane_mesher;
import .lane_shells;
import .lane_solvers;
import .lane_transient;

let CORE = 0;
let ELEMENTS = 1;
let SHELLS = 2;
let TRANSIENT = 3;

// The owner of the template `word`, or -1 if no lane owns it.
fn template_lane(word: u8[:]) -> i64 {
    if word == "linear_elastic" { return CORE; }
    if elements_owns_template(word) { return ELEMENTS; }
    if shells_owns_template(word) { return SHELLS; }
    return -1;
}

// Reads the physics item n of the template `word` of lane `lane` (not
// CORE). False if the lane has no generator for it.
fn read_template(lane: i64, word: u8[:], n: i64) -> bool {
    if lane == ELEMENTS { return elements_read_template(word, n); }
    if lane == SHELLS { return shells_read_template(word, n); }
    return false;
}

// Writes the calls of the template of lane `lane` (not CORE) that come
// after the mesh in main.
fn emit_template(out: u8[>..]&, lane: i64) {
    if lane == ELEMENTS { elements_emit_template(out); }
    if lane == SHELLS { shells_emit_template(out); }
}

// The owner of the procedure `word`, or -1 if no lane owns it.
fn procedure_lane(word: u8[:]) -> i64 {
    if word == "static" { return CORE; }
    if transient_owns_procedure(word) { return TRANSIENT; }
    return -1;
}

// Reads the step item n of the procedure `word` of lane `lane` (not CORE).
// False if the lane has no generator for it.
fn read_procedure(lane: i64, word: u8[:], n: i64) -> bool {
    if lane == TRANSIENT { return transient_read_procedure(word, n); }
    return false;
}

// Writes the function of step si of lane `lane` (not CORE).
fn emit_procedure(out: u8[>..]&, lane: i64, si: i64, results_dir: u8[:]) {
    if lane == TRANSIENT { transient_emit_procedure(out, si, results_dir); }
}
