// The registration point of the mesher lane (src/mesher/) in the code
// generator: a solid that is not one box, and its mesh. mesher_read gets
// the `let` item n of the solid and its value; false if the mesher cannot
// mesh it. mesher_surface gets the `let` item n of a surface on that solid
// and its value, and sets `set` to the name of the node set and the face
// set of the mesh that the surface is; false if it is not one. mesher_emit
// writes the calls in main that fill the mesh of src/runtime/fem.goose and
// its named sets, in place of runtime::mesh_box.
import std;
namespace codegen;

fn mesher_read(n: i64, value: i64) -> bool { return false; }

fn mesher_surface(n: i64, value: i64, set: u8[..64]&) -> bool { return false; }

fn mesher_emit(out: u8[>..]&) {}
