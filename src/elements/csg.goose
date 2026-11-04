import std;
import math;
import dictionary;
namespace elements;
import mesher from "../mesher";
import runtime from "../runtime";
import .shape;
import .mesh;

let SNAP_MAX_MOVE = 0.1;
let FOLD_RATIO = 0.2;
let SNAP_PASSES = 8;
var body_shape = -1;

fn shape_gap(g: i64, n: i64) -> f64 { return abs(mesher::phi(g, node_x(n, 0), node_x(n, 1), node_x(n, 2))); }

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
    guard sqrt(moved) <= SNAP_MAX_MOVE * l && abs(mesher::phi(g, q[0], q[1], q[2])) <= 1.0e-9 * l else { return; }
    for t in 3 { coords[3 * m + t] = q[t]; }
}

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

fn build_hex_shape(g: i64, size: f64, planes: i64[:]) -> bool {
    guard mesher::mesh_cells(g, 3, size, planes) && runtime::cell_count() > 0 else { return false; };
    guard runtime::cell_nodes.len == 8 * runtime::cell_count() else { return false; };
    var vertices: i64[>..] = [];
    vertices.append(runtime::cell_nodes);
    start_mesh(3, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
    coords.append(runtime::coords);
    conn.append(vertices);
    hex_cells = true;
    per_cell = 8;
    for t in 3 {
        bbox_lo[t] = coords[t];
        bbox_hi[t] = coords[t];
        for n in node_count() {
            bbox_lo[t] = min(bbox_lo[t], coords[3 * n + t]);
            bbox_hi[t] = max(bbox_hi[t], coords[3 * n + t]);
        }
    }
    map_lo = bbox_lo;
    map_hi = bbox_hi;
    body_shape = g;
    return true;
}

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

fn snap_edges(d: i64) {
    let tol = geometry_tolerance();
    let edges = boundary_edges(d);
    let nv = node_count();
    for c in cell_count() {
        let g = if cell_body.len > 0 { body_shapes[cell_body[c]] } else { body_shape };
        for e in edge_count(d) {
            let a = cell_node(c, edge_node(d, e, 0));
            let b = cell_node(c, edge_node(d, e, 1));
            guard edges.get(min(a, b) * nv + max(a, b)) != null else { continue; }
            if shape_gap(g, a) <= tol && shape_gap(g, b) <= tol { snap(g, cell_node(c, d + 1 + e), a, b); }
        }
    }
    for pass in SNAP_PASSES {
        var changed = false;
        for c in cell_count() {
            guard folded(d, c) else { continue; }
            for e in edge_count(d) { straighten(cell_node(c, d + 1 + e), cell_node(c, edge_node(d, e, 0)), cell_node(c, edge_node(d, e, 1))); }
            changed = true;
        }
        guard changed else { return; }
    }
}

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
    shape(dim, 1, [0.0, 0.0, 0.0], shp, dref);
    var p1: f64[9] = [0.0; 9];
    for l in d + 1 {
        for r in d {
            for t in d { p1[3 * r + t] += dref[3 * l + r] * node_x(cell_node(c, l), t); }
        }
    }
    let det0 = if d == 2 { p1[0] * p1[4] - p1[1] * p1[3] } else {
        p1[0] * (p1[4] * p1[8] - p1[5] * p1[7]) - p1[1] * (p1[3] * p1[8] - p1[5] * p1[6]) + p1[2] * (p1[3] * p1[7] - p1[4] * p1[6]) };
    return lo < FOLD_RATIO * det0;
}

fn straighten(m: i64, a: i64, b: i64) {
    for t in 3 { coords[3 * m + t] = 0.5 * (node_x(a, t) + node_x(b, t)); }
}

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

fn facet_key(d: i64, c: i64, f: i64) -> i64 {
    let nv = node_count();
    var v: i64[3] = [0, 0, 0];
    for l in d { v[l] = cell_node(c, facet_node(d, f, l)); }
    if d == 2 { return min(v[0], v[1]) * nv + max(v[0], v[1]); }
    sort(v);
    return (v[0] * nv + v[1]) * nv + v[2];
}

fn quarter_point(x: f64, y: f64) -> i64 {
    guard dim == 2 && order == 2 else { abort("quarter_point: P2 triangles"); }
    var count = 0;
    var shortest = 1.0e300;
    for n in node_count() { shortest = min(shortest, abs(node_x(n, 0) - x) + abs(node_x(n, 1) - y)); }
    let tol = 1.0e-9 * (abs(x) + abs(y) + 1.0) + shortest;
    var moved = dictionary<i64, i64> {};
    for c in cell_count() {
        for e in edge_count(2) {
            let a = cell_node(c, edge_node(2, e, 0));
            let b = cell_node(c, edge_node(2, e, 1));
            let m = cell_node(c, 3 + e);
            let a_tip = abs(node_x(a, 0) - x) + abs(node_x(a, 1) - y) <= tol;
            let b_tip = abs(node_x(b, 0) - x) + abs(node_x(b, 1) - y) <= tol;
            guard a_tip != b_tip && moved.get(m) == null else { continue; }
            moved.insert(m, 1);
            let at = if a_tip { a } else { b };
            let other = if a_tip { b } else { a };
            for t in 2 { coords[3 * m + t] = node_x(at, t) + 0.25 * (node_x(other, t) - node_x(at, t)); }
            count++;
        }
    }
    return count;
}
