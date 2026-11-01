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
    guard ord == 2 else { return true; }
    let tol = geometry_tolerance();
    let edges = boundary_edges(d);
    let nv = node_count();
    for c in cell_count() {
        for e in edge_count(d) {
            let a = cell_node(c, edge_node(d, e, 0));
            let b = cell_node(c, edge_node(d, e, 1));
            // An edge through the body, such as the third edge of a cell
            // with two edges on a round side, stays straight.
            guard edges.get(min(a, b) * nv + max(a, b)) != null else { continue; }
            if shape_gap(g, a) <= tol && shape_gap(g, b) <= tol { snap(g, cell_node(c, d + 1 + e), a, b); }
        }
    }
    return true;
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
