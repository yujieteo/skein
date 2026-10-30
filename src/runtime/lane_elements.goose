// The registration point of the elements lane (src/elements/): the
// registry entries of tri3, tri6, tet4, tet10 and hex8 (fem.goose). Each
// function gives false for a cell type that has no entry yet. A 3D cell
// with shape functions and a quadrature here can use solid_stiffness
// (solid.goose) as its element routine. The node order of each cell type
// is its VTK node order, because vtu.goose writes the nodes in that order.
//
// The P1 and P2 triangles and tetrahedra are the Lagrange simplices of
// src/elements/shape.goose, with the same reference cell, node order and
// rules. They are written again here, because a package that the runtime
// imports must be declared in each root program. The test of src/elements/
// solves a linear field with these entries to round-off. A tet4 or tet10
// cell is a solid of solid_stiffness. A triangle has no element routine
// here, because the runtime has no plane solid.
import std;
namespace runtime;
import .fem;
import .solid;

// The midpoint node d + 1 + e of a P2 cell is on edge e, from vertex
// EDGES[2 e] to vertex EDGES[2 e + 1].
let TRI_EDGES: i64[6] = [0, 1, 1, 2, 2, 0];
let TET_EDGES: i64[12] = [0, 1, 1, 2, 0, 2, 0, 3, 1, 3, 2, 3];

// The nodes of local face f of a tet10 cell: TET_FACES[6 f..6 f + 6], in
// the node order of a tri6 cell. A tet4 cell uses the first 3.
let TET_FACES: i64[24] = [
    0, 1, 2, 4, 5, 6,
    0, 1, 3, 4, 8, 7,
    1, 2, 3, 5, 9, 8,
    0, 2, 3, 6, 9, 7,
];

let TET_A = 0.1381966011250105;
let TET_B = 0.5854101966249685;

// The dimension and the order of a simplex cell type, or 0, 0.
fn simplex_of(t: i64) -> i64, i64 {
    if t == TRI3 { return 2, 1; }
    if t == TRI6 { return 2, 2; }
    if t == TET4 { return 3, 1; }
    if t == TET10 { return 3, 2; }
    return 0, 0;
}

// The derivative of barycentric coordinate k by xi_r. The reference
// simplex has its vertices at the origin and at the d unit points, so
// L0 = 1 - xi_1 - ... - xi_d and Lk = xi_k.
fn bary_grad(k: i64, r: i64) -> f64 {
    if k == 0 { return -1.0; }
    return if k - 1 == r { 1.0 } else { 0.0 };
}

// The shape functions of the simplex of dimension d and order `order`.
fn simplex_shape(d: i64, order: i64, xi: f64[3], n: f64[>..]&, dn: f64[>..]&) {
    var bc: f64[4] = [1.0, 0.0, 0.0, 0.0];
    for r in d {
        bc[r + 1] = xi[r];
        bc[0] -= xi[r];
    }
    for l in d + 1 {
        let f = if order == 1 { 1.0 } else { 4.0 * bc[l] - 1.0 };
        n.push(if order == 1 { bc[l] } else { bc[l] * (2.0 * bc[l] - 1.0) });
        for r in 3 { dn.push(f * bary_grad(l, r)); }
    }
    guard order == 2 else { return; }
    for e in if d == 2 { 3 } else { 6 } {
        let a = if d == 2 { TRI_EDGES[2 * e] } else { TET_EDGES[2 * e] };
        let b = if d == 2 { TRI_EDGES[2 * e + 1] } else { TET_EDGES[2 * e + 1] };
        n.push(4.0 * bc[a] * bc[b]);
        for r in 3 { dn.push(4.0 * (bc[b] * bary_grad(a, r) + bc[a] * bary_grad(b, r))); }
    }
}

fn elements_shape(t: i64, xi: f64[3], n: f64[>..]&, dn: f64[>..]&) -> bool {
    let d, order = simplex_of(t);
    guard d > 0 else { return false; }
    simplex_shape(d, order, xi, n, dn);
    return true;
}

// Rules of degree 2: 3 points on a triangle, 4 on a tetrahedron. The
// weights sum to the size of the reference simplex.
fn simplex_quadrature(d: i64, xi: f64[>..]&, w: f64[>..]&) {
    if d == 2 {
        xi.append([1.0 / 6.0, 1.0 / 6.0, 0.0, 2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 6.0, 2.0 / 3.0, 0.0]);
        for q in 3 { w.push(1.0 / 6.0); }
        return;
    }
    for q in 4 {
        var p: f64[3] = [TET_A, TET_A, TET_A];
        if q > 0 { p[q - 1] = TET_B; }
        xi.append(p);
        w.push(1.0 / 24.0);
    }
}

fn elements_quadrature(t: i64, xi: f64[>..]&, w: f64[>..]&) -> bool {
    let d, order = simplex_of(t);
    guard d > 0 else { return false; }
    simplex_quadrature(d, xi, w);
    return true;
}

fn elements_matrix(c: i64, ke: f64[>..]&) -> bool {
    let t = cell_type[c] as i64;
    guard t == TET4 || t == TET10 else { return false; }
    return solid_stiffness(c, ke);
}

// Local face f of a tet4 or tet10 cell c is face f of TET_FACES. The rule
// of degree 2 on the face is exact for the weights of a face with straight
// edges.
var face_xi: f64[>..] = [];
var face_qw: f64[>..] = [];
var face_n: f64[>..] = [];
var face_dn: f64[>..] = [];

fn elements_face_weights(c: i64, f: i64, w: f64[>..]&) -> bool {
    let t = cell_type[c] as i64;
    let d, order = simplex_of(t);
    guard d == 3 && f >= 0 && f < 4 else { return false; }
    let nn = CELL_NODES[t];
    let fnn = if order == 1 { 3 } else { 6 };
    for l in 4 * nn { w.push(0.0); }
    // The outward side of the face: away from the centroid of the cell.
    var out: f64[3] = [0.0, 0.0, 0.0];
    for a in 4 {
        for s in 3 { out[s] -= 0.25 * coords[3 * cell_node(c, a) + s]; }
    }
    for a in 3 {
        for s in 3 { out[s] += coords[3 * cell_node(c, TET_FACES[6 * f + a]) + s] / 3.0; }
    }
    face_xi.clear();
    face_qw.clear();
    simplex_quadrature(2, face_xi, face_qw);
    for pw, q in face_qw {
        face_n.clear();
        face_dn.clear();
        simplex_shape(2, order, [face_xi[3 * q], face_xi[3 * q + 1], 0.0], face_n, face_dn);
        var ta: f64[3] = [0.0, 0.0, 0.0];
        var tb: f64[3] = [0.0, 0.0, 0.0];
        for l in fnn {
            let n = cell_node(c, TET_FACES[6 * f + l]);
            for s in 3 {
                ta[s] += face_dn[3 * l] * coords[3 * n + s];
                tb[s] += face_dn[3 * l + 1] * coords[3 * n + s];
            }
        }
        var nv: f64[3] = [ta[1] * tb[2] - ta[2] * tb[1], ta[2] * tb[0] - ta[0] * tb[2], ta[0] * tb[1] - ta[1] * tb[0]];
        if nv[0] * out[0] + nv[1] * out[1] + nv[2] * out[2] < 0.0 {
            for s in 3 { nv[s] = -nv[s]; }
        }
        let da = pw * sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
        for l in fnn {
            let a = TET_FACES[6 * f + l];
            w[4 * a] += face_n[l] * da;
            for s in 3 { w[4 * a + 1 + s] += face_n[l] * pw * nv[s]; }
        }
    }
    return true;
}

fn elements_inside(t: i64, xi: f64[3]) -> bool {
    let d, order = simplex_of(t);
    guard d == 3 else { return false; }
    let eps = 1.0e-9;
    return xi[0] >= -eps && xi[1] >= -eps && xi[2] >= -eps && xi[0] + xi[1] + xi[2] <= 1.0 + eps;
}
