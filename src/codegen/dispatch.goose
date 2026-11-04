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
let MESHER = 4;

fn template_lane(word: u8[:]) -> i64 {
    if word == "linear_elastic" { return CORE; }
    if elements_owns_template(word) { return ELEMENTS; }
    if shells_owns_template(word) { return SHELLS; }
    return -1;
}

fn read_template(lane: i64, word: u8[:], n: i64) -> bool {
    if lane == ELEMENTS { return elements_read_template(word, n); }
    if lane == SHELLS { return shells_read_template(word, n); }
    return false;
}

fn emit_template(out: u8[>..]&, lane: i64) {
    if lane == ELEMENTS { elements_emit_template(out); }
    if lane == SHELLS { shells_emit_template(out); }
}

fn procedure_lane(word: u8[:]) -> i64 {
    if word == "static" { return CORE; }
    if transient_owns_procedure(word) { return TRANSIENT; }
    return -1;
}

fn read_procedure(lane: i64, word: u8[:], n: i64) -> bool {
    if lane == TRANSIENT { return transient_read_procedure(word, n); }
    return false;
}

fn emit_procedure(out: u8[>..]&, lane: i64, si: i64, results_dir: u8[:]) {
    if lane == TRANSIENT { transient_emit_procedure(out, si, results_dir); }
}
