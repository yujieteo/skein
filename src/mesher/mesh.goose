// The unstructured mesh that the mesher makes.
//
// This is a local copy of the shared unstructured mesh type of the
// runtime interface (task skein-p2-interface). Replace it with that type
// when it lands.
import std;
namespace mesher;

let TRIANGLE = 3;      // nodes in a cell of a 2D mesh
let TETRAHEDRON = 4;   // nodes in a cell of a 3D mesh

var mesh_dim = 0;               // 2 or 3
var points: f64[>..] = [];      // x, y and z of each node (z = 0 in 2D)
var cells: i64[>..] = [];       // dim + 1 nodes of each cell, positive orientation

fn clear_mesh(dim: i64) {
    mesh_dim = dim;
    points.clear();
    cells.clear();
}

fn point_count() -> i64 { return points.len / 3; }
fn cell_count() -> i64 { return cells.len / (mesh_dim + 1); }

fn add_point(x: f64, y: f64, z: f64) -> i64 {
    points.push(x);
    points.push(y);
    points.push(z);
    return points.len / 3 - 1;
}
