// The registration point of the shells, beams and laminates lane
// (src/shells/): the registry entries of quad4_shell and beam2 (fem.goose).
// Each function gives false for a cell type that has no entry yet. The node
// order of each cell type is its VTK node order (vtu.goose).
import std;
namespace runtime;
import .fem;

fn shells_shape(t: i64, xi: f64[3], n: f64[>..]&, dn: f64[>..]&) -> bool { return false; }

fn shells_quadrature(t: i64, xi: f64[>..]&, w: f64[>..]&) -> bool { return false; }

fn shells_matrix(c: i64, ke: f64[>..]&) -> bool { return false; }

fn shells_face_weights(c: i64, f: i64, w: f64[>..]&) -> bool { return false; }
