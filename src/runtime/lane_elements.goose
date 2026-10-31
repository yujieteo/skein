// The registration point of the elements lane (src/elements/): the
// registry entries of tri3, tri6, tet4, tet10 and hex8 (fem.goose). Each
// function gives false for a cell type that has no entry yet. The node
// order of each cell type is its VTK node order, because vtu.goose writes
// the nodes in that order.
import std;
namespace runtime;
import .fem;

fn elements_shape(t: i64, xi: f64[3], n: f64[>..]&, dn: f64[>..]&) -> bool { return false; }

fn elements_quadrature(t: i64, xi: f64[>..]&, w: f64[>..]&) -> bool { return false; }

fn elements_matrix(c: i64, ke: f64[>..]&) -> bool { return false; }

fn elements_face_weights(c: i64, f: i64, w: f64[>..]&) -> bool { return false; }

fn elements_inside(t: i64, xi: f64[3]) -> bool { return false; }
