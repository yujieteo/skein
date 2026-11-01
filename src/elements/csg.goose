// The simplex mesh of a shape of src/mesher/: the triangles or the
// tetrahedra of its mesher, with its refine targets, and the parts of the
// boundary on the surface of a shape.
import std;
import math;
import dictionary;
namespace elements;
import mesher from "../mesher";
import runtime from "../runtime";
import .shape;
import .mesh;

var body_shape = -1;    // the shape of the mesh of mesh_shape, or -1

// |φ| of shape g at node n, the distance of the node to its surface.
fn shape_gap(g: i64, n: i64) -> f64 { return abs(mesher::phi(g, node_x(n, 0), node_x(n, 1), node_x(n, 2))); }

// Moves node m, the midpoint of the edge a b, onto the surface of shape g
// by Newton steps along the gradient of φ. It stays if the move is more
// than a tenth of the edge, for example across a cut corner of the shape.
fn snap(g: i64, m: i64, a: i64, b: i64) {
    var l = 0.0;
    for t in 3 { l += (node_x(a, t) - node_x(b, t)) * (node_x(a, t) - node_x(b, t)); }
    l = sqrt(l);
    var q: f64[3] = [node_x(m, 0), node_x(m, 1), node_x(m, 2)];
    let e = 1.0e-7 * l;
    for it in 8 {
        let f = mesher::phi(g, q[0], q[1], q[2]);
        if abs(f) <= 1.0e-12 * l { break; }
        var grad: f64[3] = [0.0, 0.0, 0.0];
        for c in dim {
            var lo = q;
            var hi = q;
            lo[c] -= e;
            hi[c] += e;
            grad[c] = (mesher::phi(g, hi[0], hi[1], hi[2]) - mesher::phi(g, lo[0], lo[1], lo[2])) / (2.0 * e);
        }
        let g2 = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
        guard g2 > 0.0 else { return; }
        for c in 3 { q[c] -= f * grad[c] / g2; }
    }
    var moved = 0.0;
    for t in 3 { moved += (q[t] - node_x(m, t)) * (q[t] - node_x(m, t)); }
    guard sqrt(moved) <= 0.1 * l && abs(mesher::phi(g, q[0], q[1], q[2])) <= 1.0e-9 * l else { return; }
    for t in 3 { coords[3 * m + t] = q[t]; }
}

// Meshes shape g with the mesher, P1 or P2 (ord) cells of size `size` in
// dimension d, and smaller at its refine targets. For order 2, the
// midpoint node of each edge on the surface of g moves onto the surface.
// False if the mesher does not finish.
fn build_shape(d: i64, ord: i64, g: i64, size: f64) -> bool {
    guard d == 2 || d == 3 else { abort("build_shape: d = 2 or 3"); }
    let done = if d == 3 { mesher::mesh_3d(g, size) } else { mesher::mesh_2d(g, size) };
    guard done && runtime::cell_count() > 0 else { return false; }
    var vertices: i64[>..] = [];
    vertices.append(runtime::cell_nodes);
    build_simplices(d, ord, runtime::coords, vertices);
    body_shape = g;
    if ord == 2 { snap_edges(d); }
    return true;
}

// Meshes each shape of `shapes` alone, as build_shape, so that two bodies
// share no node: the bodies of a contact (contact.goose). Cell c is in
// body cell_body[c]. g is the union of the shapes, the body of the
// boundary planes. False if the mesher does not finish.
fn build_bodies(d: i64, ord: i64, shapes: i64[:], g: i64, size: f64) -> bool {
    guard (d == 2 || d == 3) && shapes.len >= 2 else { abort("build_bodies: d = 2 or 3 and two shapes or more"); }
    var xyz: f64[>..] = [];
    var vertices: i64[>..] = [];
    var owner: i64[>..] = [];
    for s, k in shapes {
        let done = if d == 3 { mesher::mesh_3d(s, size) } else { mesher::mesh_2d(s, size) };
        guard done && runtime::cell_count() > 0 else { return false; }
        let base = xyz.len / 3;
        xyz.append(runtime::coords);
        for v in runtime::cell_nodes { vertices.push(base + v); }
        for c in runtime::cell_count() { owner.push(k); }
    }
    // The mesh of the mesher holds the last body: it gets all of them, for
    // mesh_quality.
    let t = runtime::cell_type[0];
    runtime::coords.clear();
    runtime::coords.append(xyz);
    runtime::cell_nodes.clear();
    runtime::cell_nodes.append(vertices);
    runtime::cell_type.clear();
    runtime::cell_start.clear();
    runtime::cell_start.push(0);
    for c in owner.len {
        runtime::cell_type.push(t);
        runtime::cell_start.push((c + 1) * (d + 1));
    }
    build_simplices(d, ord, xyz, vertices);
    body_shape = g;
    body_shapes.append(shapes);
    cell_body.append(owner);
    for n in node_count() { node_body.push(0); }
    for c in cell_count() {
        for l in per_cell { node_body[cell_node(c, l)] = cell_body[c]; }
    }
    if ord == 2 { snap_edges(d); }
    return true;
}

// Moves the midpoint node of each edge on the surface of the body of its
// cell onto that surface.
fn snap_edges(d: i64) {
    let tol = geometry_tolerance();
    let edges = boundary_edges(d);
    let nv = node_count();
    for c in cell_count() {
        let g = if cell_body.len > 0 { body_shapes[cell_body[c]] } else { body_shape };
        for e in edge_count(d) {
            let a = cell_node(c, edge_node(d, e, 0));
            let b = cell_node(c, edge_node(d, e, 1));
            // An edge through the body, such as the third edge of a cell
            // with two edges on a round side, stays straight.
            guard edges.get(min(a, b) * nv + max(a, b)) != null else { continue; }
            if shape_gap(g, a) <= tol && shape_gap(g, b) <= tol { snap(g, cell_node(c, d + 1 + e), a, b); }
        }
    }
    // A snap across a corner of the shape can fold a cell: such a cell gets
    // its straight edges back. A straight cell is never folded, so this ends.
    for pass in 8 {
        var changed = false;
        for c in cell_count() {
            guard folded(d, c) else { continue; }
            for e in edge_count(d) { straighten(cell_node(c, d + 1 + e), cell_node(c, edge_node(d, e, 0)), cell_node(c, edge_node(d, e, 1))); }
            changed = true;
        }
        guard changed else { return; }
    }
}

// Whether the Jacobian of cell c at a vertex or at a quadrature point is
// less than 0.2 times that of its straight cell.
fn folded(d: i64, c: i64) -> bool {
    var lo = 1.0e300;
    var xi: f64[3] = [0.0, 0.0, 0.0];
    for q in quad_count(d) + d + 1 {
        xi = [0.0, 0.0, 0.0];
        if q < quad_count(d) {
            let p, w = quad_point(d, q);
            xi = p;
        } else if q > quad_count(d) {
            xi[q - quad_count(d) - 1] = 1.0;
        }
        shape(dim, order, xi, shp, dref);
        var jac: f64[9] = [0.0; 9];
        lo = min(lo, cell_jacobian(c, jac));
    }
    // the straight cell: the Jacobian at the centroid of its vertices
    shape(dim, 1, [0.0, 0.0, 0.0], shp, dref);
    var p1: f64[9] = [0.0; 9];
    for l in d + 1 {
        for r in d {
            for t in d { p1[3 * r + t] += dref[3 * l + r] * node_x(cell_node(c, l), t); }
        }
    }
    let det0 = if d == 2 { p1[0] * p1[4] - p1[1] * p1[3] } else {
        p1[0] * (p1[4] * p1[8] - p1[5] * p1[7]) - p1[1] * (p1[3] * p1[8] - p1[5] * p1[6]) + p1[2] * (p1[3] * p1[7] - p1[4] * p1[6]) };
    return lo < 0.2 * det0;
}

// Puts midpoint node m back on the middle of a and b.
fn straighten(m: i64, a: i64, b: i64) {
    for t in 3 { coords[3 * m + t] = 0.5 * (node_x(a, t) + node_x(b, t)); }
}

// The edges of the facets on the boundary of the mesh, each a key
// min(a, b) * node_count() + max(a, b) of its vertices a and b. A facet on
// the boundary is a facet of one cell.
fn boundary_edges(d: i64) -> dictionary<i64, i64> {
    let nv = node_count();
    guard nv < 2000000 else { abort("boundary_edges: fewer than 2e6 nodes"); }
    var count = dictionary<i64, i64> {};
    for c in cell_count() {
        for f in facet_count(d) { count.update(facet_key(d, c, f), 0) { it += 1; }; }
    }
    var edges = dictionary<i64, i64> {};
    for c in cell_count() {
        for f in facet_count(d) {
            let shared = count.get(facet_key(d, c, f));
            guard shared != null && shared == 1 else { continue; }
            for i in d {
                for j in i + 1..d {
                    let a = cell_node(c, facet_node(d, f, i));
                    let b = cell_node(c, facet_node(d, f, j));
                    edges.insert(min(a, b) * nv + max(a, b), 0);
                }
            }
        }
    }
    return edges;
}

// The key of facet f of cell c: its sorted vertices in base node_count().
fn facet_key(d: i64, c: i64, f: i64) -> i64 {
    let nv = node_count();
    var v: i64[3] = [0, 0, 0];
    for l in d { v[l] = cell_node(c, facet_node(d, f, l)); }
    if d == 2 { return min(v[0], v[1]) * nv + max(v[0], v[1]); }
    sort(v);
    return (v[0] * nv + v[1]) * nv + v[2];
}
