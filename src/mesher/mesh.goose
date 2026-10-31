// The output of the mesher: the shared unstructured mesh of the runtime
// (src/runtime/fem.goose), with tri3 cells in 2D and tet4 cells in 3D, or
// the quad4 and hex8 cells of a grid (grid.goose).
import std;
namespace mesher;
import runtime from "../runtime";

var mesh_dim = 0;   // 2 or 3

fn clear_mesh(dim: i64) {
    mesh_dim = dim;
    runtime::mesh_clear();
}

fn point_count() -> i64 { return runtime::node_count(); }
fn cell_count() -> i64 { return runtime::cell_count(); }

fn add_point(x: f64, y: f64, z: f64) -> i64 { return runtime::add_node(x, y, z); }

// Adds a cell of dim + 1 nodes in positive orientation.
fn add_simplex(nodes: i64[:]) {
    runtime::add_cell(if mesh_dim == 2 { runtime::TRI3 } else { runtime::TET4 }, nodes);
}

// Adds a cell of a grid: a quad4 in 2D, a hex8 in 3D, in VTK node order.
fn add_grid_cell(nodes: i64[:]) {
    runtime::add_cell(if mesh_dim == 2 { runtime::QUAD4 } else { runtime::HEX8 }, nodes);
}

// Node i of cell c.
fn cell_node(c: i64, i: i64) -> i64 { return runtime::cell_node(c, i); }

fn point(n: i64) -> f64[3] {
    return [runtime::coords[3 * n], runtime::coords[3 * n + 1], runtime::coords[3 * n + 2]];
}
