// The registration point of the elements and templates lane
// (src/elements/) in the code generator: the template heat (dispatch.goose)
// and the cells of the mesh other than hex with order 2.
import std;
namespace codegen;

fn elements_owns_template(word: u8[:]) -> bool { return word == "heat"; }

fn elements_read_template(word: u8[:], n: i64) -> bool { return false; }

fn elements_emit_template(out: u8[>..]&) {}

// Whether the lane meshes the box with the mesh setting `cell: <cell>` and
// the method setting `order: <order>`.
fn elements_cell(cell: u8[:], order: i64) -> bool { return false; }

// Writes the calls in main that mesh the box with the cells of
// elements_cell, in place of runtime::mesh_box. False if the model uses
// hex cells of order 2.
fn elements_emit_mesh(out: u8[>..]&) -> bool { return false; }
