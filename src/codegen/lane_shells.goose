// The registration point of the shells, beams and laminates lane
// (src/shells/) in the code generator: its templates (dispatch.goose).
import std;
namespace codegen;

fn shells_owns_template(word: u8[:]) -> bool { return false; }

fn shells_read_template(word: u8[:], n: i64) -> bool { return false; }

fn shells_emit_template(out: u8[>..]&) {}
