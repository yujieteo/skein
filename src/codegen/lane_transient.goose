// The registration point of the transient lane (src/transient/) in the code
// generator: the procedure transient (dispatch.goose).
import std;
namespace codegen;

fn transient_owns_procedure(word: u8[:]) -> bool { return word == "transient"; }

fn transient_read_procedure(word: u8[:], n: i64) -> bool { return false; }

fn transient_emit_procedure(out: u8[>..]&, si: i64, results_dir: u8[:]) {}
