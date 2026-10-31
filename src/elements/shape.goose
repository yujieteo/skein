// The shape functions and the quadrature of the Lagrange simplices: lines,
// triangles and tetrahedra of order 1 (P1) and order 2 (P2).
//
// The local nodes follow VTK: the vertices, then the midpoint of each edge
// in the sequence of `edge_node`. The reference simplex of dimension d has
// its vertices at the origin and at the d unit points, so the barycentric
// coordinates are L0 = 1 - xi_1 - ... - xi_d and Lk = xi_k.
import std;
import math;
namespace elements;

let MAX_NODES = 10;     // a 10-node tetrahedron

// The edges of the line, the triangle and the tetrahedron, as pairs of
// vertices. The midpoint node of edge e is node d + 1 + e.
let LINE_EDGES: i64[2] = [0, 1];
let TRI_EDGES: i64[6] = [0, 1, 1, 2, 2, 0];
let TET_EDGES: i64[12] = [0, 1, 1, 2, 0, 2, 0, 3, 1, 3, 2, 3];

fn edge_count(d: i64) -> i64 {
    guard d >= 1 && d <= 3 else { abort("edge_count: 1 <= d <= 3"); }
    if d == 1 { return 1; }
    if d == 2 { return 3; }
    return 6;
}

// Vertex `side` (0 or 1) of edge e of the simplex of dimension d.
fn edge_node(d: i64, e: i64, side: i64) -> i64 {
    guard e >= 0 && e < edge_count(d) && (side == 0 || side == 1) else { abort("edge_node: e < edge_count(d) and side 0 or 1"); }
    if d == 1 { return LINE_EDGES[2 * e + side]; }
    if d == 2 { return TRI_EDGES[2 * e + side]; }
    return TET_EDGES[2 * e + side];
}

fn node_count_of(d: i64, order: i64) -> i64 {
    guard d >= 1 && d <= 3 && (order == 1 || order == 2) else { abort("node_count_of: 1 <= d <= 3 and order 1 or 2"); }
    if order == 1 { return d + 1; }
    return d + 1 + edge_count(d);
}

// The derivative of barycentric coordinate k by xi_r.
fn bary_grad(k: i64, r: i64) -> f64 {
    guard k >= 0 && k <= 3 && r >= 0 && r < 3 else { abort("bary_grad: 0 <= k <= 3 and 0 <= r < 3"); }
    if k == 0 { return -1.0; }
    if k - 1 == r { return 1.0; }
    return 0.0;
}

// The values n[l] and the reference gradients dn[3 l + r] = dN_l / dxi_r of
// the shape functions at xi.
fn shape(d: i64, order: i64, xi: f64[3], n: f64[10]&, dn: f64[30]&) {
    guard node_count_of(d, order) <= n.len && 3 * n.len <= dn.len else { abort("shape: n and dn hold the nodes of the simplex"); }
    var bc: f64[4] = [1.0, 0.0, 0.0, 0.0];
    for r in d {
        bc[r + 1] = xi[r];
        bc[0] -= xi[r];
    }
    for l in d + 1 {
        if order == 1 {
            n[l] = bc[l];
            for r in d { dn[3 * l + r] = bary_grad(l, r); }
        } else {
            n[l] = bc[l] * (2.0 * bc[l] - 1.0);
            for r in d { dn[3 * l + r] = (4.0 * bc[l] - 1.0) * bary_grad(l, r); }
        }
    }
    guard order == 2 else { return; }
    for e in edge_count(d) {
        let a = edge_node(d, e, 0);
        let b = edge_node(d, e, 1);
        let l = d + 1 + e;
        n[l] = 4.0 * bc[a] * bc[b];
        for r in d { dn[3 * l + r] = 4.0 * (bc[b] * bary_grad(a, r) + bc[a] * bary_grad(b, r)); }
    }
}

// --- quadrature --------------------------------------------------------------

// Rules of degree 2 on the reference simplex, so the stiffness of a P2 cell
// with straight edges is exact: 2 Gauss points on a line, 3 points on a
// triangle and 4 points on a tetrahedron. The weights sum to the size of
// the reference simplex.
let TET_A = 0.1381966011250105;
let TET_B = 0.5854101966249685;
let LINE_G = 0.21132486540518713;   // (1 - 1/sqrt(3)) / 2

fn quad_count(d: i64) -> i64 {
    guard d >= 1 && d <= 3 else { abort("quad_count: 1 <= d <= 3"); }
    return d + 1;
}

fn quad_point(d: i64, q: i64) -> f64[3], f64 {
    guard q >= 0 && q < quad_count(d) else { abort("quad_point: 0 <= q < quad_count(d)"); }
    if d == 1 {
        let s = if q == 0 { LINE_G } else { 1.0 - LINE_G };
        return [s, 0.0, 0.0], 0.5;
    }
    if d == 2 {
        if q == 0 { return [1.0 / 6.0, 1.0 / 6.0, 0.0], 1.0 / 6.0; }
        if q == 1 { return [2.0 / 3.0, 1.0 / 6.0, 0.0], 1.0 / 6.0; }
        return [1.0 / 6.0, 2.0 / 3.0, 0.0], 1.0 / 6.0;
    }
    var p: f64[3] = [TET_A, TET_A, TET_A];
    if q > 0 { p[q - 1] = TET_B; }
    return p, 1.0 / 24.0;
}

// --- facets ------------------------------------------------------------------

// Local node l of facet f of a cell of dimension d, in the node sequence of
// the facet as a simplex of dimension d - 1: an edge of a triangle, or a
// face of a tetrahedron.
let TET_FACES: i64[24] = [
    0, 1, 2, 4, 5, 6,
    0, 1, 3, 4, 8, 7,
    1, 2, 3, 5, 9, 8,
    0, 2, 3, 6, 9, 7,
];

fn facet_count(d: i64) -> i64 {
    guard d == 2 || d == 3 else { abort("facet_count: d is 2 or 3"); }
    return d + 1;
}

fn facet_node(d: i64, f: i64, l: i64) -> i64 {
    guard f >= 0 && f < facet_count(d) && l >= 0 && l < node_count_of(d - 1, 2) else { abort("facet_node: f < facet_count(d) and l < the nodes of a facet"); }
    if d == 3 { return TET_FACES[6 * f + l]; }
    if l < 2 { return TRI_EDGES[2 * f + l]; }
    return 3 + f;
}
