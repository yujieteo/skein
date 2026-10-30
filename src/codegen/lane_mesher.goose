// The registration point of the mesher lane (src/mesher/) in the code
// generator: a solid that is not one box, and its mesh. mesher_read gets
// the `let` item n of the solid and its value; false if the mesher cannot
// mesh it. mesher_emit writes the calls in main that fill the mesh of
// src/runtime/fem.goose, in place of runtime::mesh_box.
import std;
namespace codegen;

fn mesher_read(n: i64, value: i64) -> bool { return false; }

fn mesher_emit(out: u8[>..]&) {}
