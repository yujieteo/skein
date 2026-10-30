// The registration point of the elements lane (src/elements/): the
// registry entries of tri3, tri6, tet4, tet10 and hex8 (fem.goose). Each
// function gives false for a cell type that has no entry yet. A 3D cell
// with shape functions and a quadrature here can use solid_stiffness
// (solid.goose) as its element routine.
import std;
namespace runtime;
import .fem;

fn elements_shape(t: i64, xi: f64[3], n: f64[>..]&, dn: f64[>..]&) -> bool { return false; }

fn elements_quadrature(t: i64, xi: f64[>..]&, w: f64[>..]&) -> bool { return false; }

fn elements_matrix(c: i64, ke: f64[>..]&) -> bool { return false; }
