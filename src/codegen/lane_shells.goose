// The registration point of the shells, beams and laminates lane
// (src/shells/) in the code generator: its templates (dispatch.goose), and
// a shell or beam geometry in a model with no solid. shells_read_geometry
// gets the `let` item n of a Surface or a Curve and its value; true if the
// lane meshes it as the geometry of the model. shells_surface gets a later
// Surface or Curve item on that geometry and sets `set` to the name of the
// node set and the face set of the mesh that it is; false if it is not
// one. shells_emit_mesh writes the calls in main that fill the mesh of
// src/runtime/fem.goose and its named sets, in place of runtime::mesh_box.
import std;
namespace codegen;

fn shells_owns_template(word: u8[:]) -> bool { return false; }

fn shells_read_template(word: u8[:], n: i64) -> bool { return false; }

fn shells_emit_template(out: u8[>..]&) {}

fn shells_read_geometry(n: i64, value: i64) -> bool { return false; }

fn shells_surface(n: i64, value: i64, set: u8[..64]&) -> bool { return false; }

fn shells_emit_mesh(out: u8[>..]&) {}
