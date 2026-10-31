// `quad` cells on a 2D shape that is not a union of rects (spec §9.2):
// Blossom-Quad after Remacle et al. (2012). The triangles of mesh_2d are
// the vertices of a graph, with an edge between two neighbour triangles
// when their quad has an acceptable quality. A matching of the largest
// size in this graph (Edmonds 1965), started from the best pairs first,
// gives the quads. A triangle with no partner stays.
//
// For the swept hex grid (sweep.goose), all_quads cuts each quad into 4
// quads and each triangle into 3, so that no triangle stays.
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .tri;
import runtime from "../runtime";

// The lowest quality of a quad from two triangles: 1 - (2/π) times the
// largest difference of a corner angle from π/2 (Remacle et al. 2012).
// 0.25 accepts the corner angles from 22.5 to 157.5 degrees.
let QUAD_QUALITY = 0.25;

// The quality of the quad of points a, b, c, d in counter-clockwise order,
// or -1 if it is not convex.
fn quad_quality(q: i64[4]) -> f64 {
    var worst = 0.0;
    for i in 4 {
        let o = point(q[i]);
        let a = point(q[(i + 1) % 4]);
        let b = point(q[(i + 3) % 4]);
        let ax = a[0] - o[0];
        let ay = a[1] - o[1];
        let bx = b[0] - o[0];
        let by = b[1] - o[1];
        let cross = ax * by - ay * bx;
        guard cross > 0.0 else { return -1.0; };
        let angle = atan2(cross, ax * bx + ay * by);
        worst = max(worst, abs(0.5 * PI - angle));
    }
    return 1.0 - 2.0 / PI * worst;
}

// The quad of triangle t and its neighbour u across their shared edge, in
// counter-clockwise order. t goes r, p, q and u goes s, q, p, so the quad
// is p, s, q, r.
fn pair_quad(t: i64, u: i64) -> i64[4] {
    for k in 3 {
        let p = cell_node(t, k);
        let q = cell_node(t, (k + 1) % 3);
        for j in 3 {
            if cell_node(u, j) == q && cell_node(u, (j + 1) % 3) == p {
                return [p, cell_node(u, (j + 2) % 3), q, cell_node(t, (k + 2) % 3)];
            }
        }
    }
    abort("pair_quad: the triangles are neighbours");
    return [0, 0, 0, 0];
}

// The edges of the graph: for each triangle, its neighbours with an
// acceptable quad. adj[start[t]..start[t + 1]] are the neighbours of t.
fn quad_graph(start: i64[>..]&, adj: i64[>..]&, weight: f64[>..]&) {
    let n = point_count();
    var edge = dictionary<i64, i64> {};
    for t in cell_count() {
        for k in 3 {
            let a = cell_node(t, k);
            let b = cell_node(t, (k + 1) % 3);
            edge.insert(a * n + b, t);
        }
    }
    start.clear();
    start.push(0);
    for t in cell_count() {
        for k in 3 {
            let u = edge.get(cell_node(t, (k + 1) % 3) * n + cell_node(t, k));
            guard u else { continue; };
            let w = quad_quality(pair_quad(t, u));
            guard w >= QUAD_QUALITY else { continue; };
            adj.push(u);
            weight.push(w);
        }
        start.push(adj.len);
    }
}

// The base of the blossom that holds v and w: the first vertex that both
// alternating paths to the root reach.
fn common_base(v: i64, w: i64, mate: i64[:], parent: i64[:], base: i64[:], seen: bool[>..]&) -> i64 {
    for i in seen.len { seen[i] = false; }
    var a = v;
    loop {
        a = base[a];
        seen[a] = true;
        guard mate[a] >= 0 else { break; };
        a = parent[mate[a]];
    }
    var b = w;
    loop {
        b = base[b];
        if seen[b] { return b; }
        b = parent[mate[b]];
    }
    return -1;
}

// Marks the blossom from v down to the base b, and points the parents of
// the path to `child`.
fn mark_blossom(v0: i64, b: i64, child0: i64, mate: i64[:], parent: i64[>..]&, base: i64[:], blossom: bool[>..]&) {
    var v = v0;
    var child = child0;
    while base[v] != b {
        blossom[base[v]] = true;
        blossom[base[mate[v]]] = true;
        parent[v] = child;
        child = mate[v];
        v = parent[mate[v]];
    }
}

// Searches an augmenting path from the free vertex root (Edmonds 1965), and
// flips it. True if it found one.
fn augment(root: i64, start: i64[:], adj: i64[:], mate: i64[>..]&, parent: i64[>..]&, base: i64[>..]&,
           used: bool[>..]&, blossom: bool[>..]&, seen: bool[>..]&) -> bool {
    let n = mate.len;
    for i in n {
        parent[i] = -1;
        base[i] = i;
        used[i] = false;
    }
    used[root] = true;
    var queue: i64[>..] = [root];
    var head = 0;
    while head < queue.len {
        let v = queue[head];
        head++;
        for e in start[v]..start[v + 1] {
            let w = adj[e];
            guard base[v] != base[w] && mate[v] != w else { continue; };
            if w == root || (mate[w] >= 0 && parent[mate[w]] >= 0) {
                let b = common_base(v, w, mate[..], parent[..], base[..], seen);
                for i in n { blossom[i] = false; }
                mark_blossom(v, b, w, mate[..], parent, base[..], blossom);
                mark_blossom(w, b, v, mate[..], parent, base[..], blossom);
                for i in n {
                    guard blossom[base[i]] else { continue; };
                    base[i] = b;
                    if !used[i] {
                        used[i] = true;
                        queue.push(i);
                    }
                }
            } else if parent[w] < 0 {
                parent[w] = v;
                if mate[w] < 0 {
                    var x = w;
                    while x >= 0 {
                        let px = parent[x];
                        let next = mate[px];
                        mate[x] = px;
                        mate[px] = x;
                        x = next;
                    }
                    return true;
                }
                used[mate[w]] = true;
                queue.push(mate[w]);
            }
        }
    }
    return false;
}

// The partner of each triangle, or -1: the best pairs first, then each
// augmenting path.
fn quad_matching(start: i64[:], adj: i64[:], weight: f64[:]) -> i64[>..] {
    let n = start.len - 1;
    var order: i64[>..] = [];
    for e in adj.len { order.push(e); }
    sort(order) { a, b => weight[a] > weight[b] };
    var from: i64[>..] = [];
    for t in n {
        for e in start[t]..start[t + 1] { from.push(t); }
    }
    var mate: i64[>..] = [];
    for t in n { mate.push(-1); }
    for e in order {
        let t = from[e];
        let u = adj[e];
        if mate[t] < 0 && mate[u] < 0 {
            mate[t] = u;
            mate[u] = t;
        }
    }
    var parent: i64[>..] = [];
    var base: i64[>..] = [];
    var used: bool[>..] = [];
    var blossom: bool[>..] = [];
    var seen: bool[>..] = [];
    for t in n {
        parent.push(-1);
        base.push(t);
        used.push(false);
        blossom.push(false);
        seen.push(false);
    }
    for t in n {
        guard mate[t] < 0 && start[t + 1] > start[t] else { continue; };
        augment(t, start, adj, mate, parent, base, used, blossom, seen);
    }
    return mate;
}

// Removes the cells of the mesh and keeps its points.
fn clear_cells() {
    runtime::cell_type.clear();
    runtime::cell_start.clear();
    runtime::cell_start.push(0);
    runtime::cell_nodes.clear();
}

// Merges the triangles of the mesh into quads, where a pair has an
// acceptable quality.
fn merge_quads() {
    var start: i64[>..] = [];
    var adj: i64[>..] = [];
    var weight: f64[>..] = [];
    quad_graph(start, adj, weight);
    let mate = quad_matching(start[..], adj[..], weight[..]);
    var tris: i64[>..] = [];
    for t in cell_count() {
        for k in 3 { tris.push(cell_node(t, k)); }
    }
    var quads: i64[>..] = [];
    for t in mate.len {
        guard mate[t] > t else { continue; };
        for v in pair_quad(t, mate[t]) { quads.push(v); }
    }
    clear_cells();
    for q in quads.len / 4 { add_grid_cell(quads[4 * q..4 * q + 4]); }
    for t in mate.len {
        if mate[t] < 0 { add_simplex(tris[3 * t..3 * t + 3]); }
    }
}

// Meshes the 2D shape g with quads and the triangles that stay. False as
// mesh_2d.
fn mesh_blossom(g: i64, size: f64) -> bool {
    let done = mesh_2d(g, size);
    merge_quads();
    return done;
}

// Cuts each quad of the mesh into 4 quads, and each triangle into 3, at
// the midpoints of its edges and its centre. A midpoint of an edge on the
// boundary of g moves onto the boundary.
fn all_quads(g: i64) {
    let n = point_count();
    var cells: i64[>..] = [];
    var sides: i64[>..] = [];
    var count = dictionary<i64, i64> {};
    for c in cell_count() {
        let k = if runtime::cell_type[c] as i64 == runtime::QUAD4 { 4 } else { 3 };
        sides.push(k);
        for i in k {
            let a = cell_node(c, i);
            let b = cell_node(c, (i + 1) % k);
            cells.push(a);
            count.update(min(a, b) * n + max(a, b), 0) { it += 1; };
        }
    }
    var mid = dictionary<i64, i64> {};
    fn midpoint(a: i64, b: i64) -> i64 {
        let key = min(a, b) * n + max(a, b);
        let have = mid.get(key);
        if have { return have; }
        let pa = point(a);
        let pb = point(b);
        var x: f64[3] = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.0];
        let shared = count.get(key);
        if shared && shared == 1 { x = project(g, x); }
        let m = add_point(x[0], x[1], 0.0);
        mid.insert(key, m);
        return m;
    }
    clear_cells();
    var at = 0;
    for k in sides {
        let v = cells[at..at + k];
        var cx = 0.0;
        var cy = 0.0;
        for i in k {
            let p = point(v[i]);
            cx += p[0] / k;
            cy += p[1] / k;
        }
        let o = add_point(cx, cy, 0.0);
        for i in k {
            let before = midpoint(v[(i + k - 1) % k], v[i]);
            let after = midpoint(v[i], v[(i + 1) % k]);
            add_grid_cell([v[i], after, o, before]);
        }
        at += k;
    }
}
