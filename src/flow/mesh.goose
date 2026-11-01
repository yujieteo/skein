// The mesh of the flow solver: P1 or P2 triangles from the triangles that
// skein run meshes with src/mesher/, the dofs of the velocity and the
// pressure, the sparsity of the Jacobian, and the boundary regions.
//
// Nodes 0 to vertex_count - 1 are the vertices of the triangles. For order 2,
// each edge then gets a node at its midpoint. Local node 3 + e of a cell is
// on edge e, between local vertices EDGES[2 e] and EDGES[2 e + 1].
//
// The dofs follow the reverse Cuthill-McKee sequence of the nodes, so that
// ILU(0) is near the full LU. The dofs of node n start at first_dof[n]:
// vx, vy, and p for a vertex. Reason: with all the velocity dofs first,
// GMRES(50) with ILU(0) did not converge in 5000 iterations on a P2
// cavity of 24174 dofs. With the dofs of each node together and the
// preconditioner of solve.goose, it converged.
//
// A turbulent flow (rans.goose) has P1 cells, and each vertex then has 5
// dofs: vx, vy, p, k and w.
import std;
import math;
import dictionary;
import os;
namespace flow;

let EDGES: i64[6] = [0, 1, 1, 2, 2, 0];

var order = 1;          // the order of the velocity; the pressure has order 1
var turbulent = false;  // the k-ω SST unknowns k and w at each vertex
var per_cell = 3;
var vertex_count = 0;
var coords: f64[>..] = [];      // x and y of each node
var conn: i64[>..] = [];        // the per_cell nodes of each cell
var on_boundary: u8[>..] = [];  // 1 for a node on the boundary of the mesh
var first_dof: i64[>..] = [];   // the first dof of each node
var twin: i64[>..] = [];        // the node whose dofs a node shares (periodic), or -1
var dof_total = 0;
var diameter = 0.0;             // of the bounding box

// The sparsity of the Jacobian: CSR, with the columns of each row in
// increasing order.
var row_start: i64[>..] = [];
var col: i32[>..] = [];

// The region of the last select call: a flag for each node.
var region: u8[>..] = [];

// A turbulent flow has the unknowns k and w. Call it before the mesh.
fn set_turbulent(on: bool) { turbulent = on; }

fn node_count() -> i64 { return coords.len / 2; }
fn cell_count() -> i64 { return conn.len / per_cell; }
fn vertex_dofs() -> i64 { return if turbulent { 5 } else { 3 }; }
fn dof_count() -> i64 { return dof_total; }

fn vdof(n: i64, comp: i64) -> i64 {
    guard n >= 0 && n < first_dof.len && (comp == 0 || comp == 1) else { abort("vdof: n < node_count() and comp 0 or 1"); }
    return first_dof[n] + comp;
}

fn pdof(n: i64) -> i64 {
    guard n >= 0 && n < vertex_count && n < first_dof.len else { abort("pdof: n < vertex_count"); }
    return first_dof[n] + 2;
}

// The dof of k (which 0) or w (which 1) at vertex n of a turbulent flow.
fn turbulence_dof(n: i64, which: i64) -> i64 {
    guard turbulent && n >= 0 && n < vertex_count && n < first_dof.len && (which == 0 || which == 1) else { abort("turbulence_dof: a turbulent flow, n < vertex_count and which 0 or 1"); }
    return first_dof[n] + 3 + which;
}

fn cell_node(c: i64, l: i64) -> i64 {
    guard c >= 0 && c < cell_count() && l >= 0 && l < per_cell else { abort("cell_node: c < cell_count() and l < per_cell"); }
    return conn[per_cell * c + l];
}

// The dofs of cell c in the local sequence of cell.goose: vx and vy of
// each node, then p of each vertex, then for a turbulent flow k and then w
// of each vertex.
fn cell_dofs(c: i64, dofs: i64[15]&) {
    guard c >= 0 && c < cell_count() else { abort("cell_dofs: c < cell_count()"); }
    for l in per_cell {
        dofs[2 * l] = vdof(cell_node(c, l), 0);
        dofs[2 * l + 1] = vdof(cell_node(c, l), 1);
    }
    for k in 3 { dofs[2 * per_cell + k] = pdof(cell_node(c, k)); }
    guard turbulent else { return; }
    for k in 3 {
        dofs[2 * per_cell + 3 + k] = turbulence_dof(cell_node(c, k), 0);
        dofs[2 * per_cell + 6 + k] = turbulence_dof(cell_node(c, k), 1);
    }
}

fn cell_dof_count() -> i64 { return 2 * per_cell + (if turbulent { 9 } else { 3 }); }

// Takes the vertices `xy` (x and y of each) and the triangles `tri` (3
// vertices each), turns each triangle to positive orientation, and adds
// the edge nodes for order 2.
fn mesh_from_arrays(xy: f64[:], tri: i64[:], ord: i64) {
    guard (ord == 1 || ord == 2) && (ord == 1 || !turbulent) && xy.len % 2 == 0 && tri.len % 3 == 0 && tri.len > 0 else {
        abort("mesh_from_arrays: order 1 or 2 (1 for a turbulent flow), and triangles");
    }
    order = ord;
    per_cell = if ord == 1 { 3 } else { 6 };
    vertex_count = xy.len / 2;
    coords.clear();
    coords.append(xy);
    conn.clear();
    var mid = dictionary<i64, i64> {};
    var count = dictionary<i64, i64> {};
    for c in tri.len / 3 {
        var v: i64[3] = [tri[3 * c], tri[3 * c + 1], tri[3 * c + 2]];
        guard v[0] >= 0 && v[1] >= 0 && v[2] >= 0 && v[0] < vertex_count && v[1] < vertex_count && v[2] < vertex_count else {
            abort("mesh_from_arrays: each vertex of a triangle < xy.len / 2");
        }
        if signed_area(v) < 0.0 {
            let t = v[1];
            v[1] = v[2];
            v[2] = t;
        }
        conn.append(v[..]);
        for e in 3 {
            let a = v[EDGES[2 * e]];
            let b = v[EDGES[2 * e + 1]];
            let key = min(a, b) * vertex_count + max(a, b);
            count.update(key, 0) { it += 1; };
            guard ord == 2 else { continue; }
            let m = mid.get_or_insert(key, node_count());
            if m == node_count() {
                coords.push(0.5 * (coords[2 * a] + coords[2 * b]));
                coords.push(0.5 * (coords[2 * a + 1] + coords[2 * b + 1]));
            }
            conn.push(m);
        }
    }
    mark_boundary(count, mid);
    twin.clear();
    for n in node_count() { twin.push(-1); }
    order_nodes();
    build_pattern();
}

fn signed_area(v: i64[3]) -> f64 {
    let ax = coords[2 * v[0]];
    let ay = coords[2 * v[0] + 1];
    return 0.5 * ((coords[2 * v[1]] - ax) * (coords[2 * v[2] + 1] - ay) - (coords[2 * v[1] + 1] - ay) * (coords[2 * v[2]] - ax));
}

// Reads the mesh that skein run wrote in the directory `dir`: the files
// flow.xy and flow.tri, the bytes of the arrays of mesh_from_arrays. False
// if it cannot read them.
fn read_mesh(dir: const u8[:], ord: i64) -> bool {
    var a: u8[>..] = [];
    var b: u8[>..] = [];
    guard read_file(str(dir, "/flow.xy"), a) && read_file(str(dir, "/flow.tri"), b) else { return false; }
    let xy, xy_ok = from_bytes<f64[>..]>(a);
    let tri, tri_ok = from_bytes<i64[>..]>(b);
    guard xy_ok && tri_ok && tri.len > 0 else { return false; }
    mesh_from_arrays(xy, tri, ord);
    return true;
}

// Marks the nodes of the edges of one cell, and sets the diameter.
fn mark_boundary(count: dictionary<i64, i64>&, mid: dictionary<i64, i64>&) {
    on_boundary.clear();
    region.clear();
    for n in node_count() {
        on_boundary.push(0);
        region.push(0);
    }
    count.each() { key, k =>
        if k == 1 {
            on_boundary[key / vertex_count] = 1;
            on_boundary[key % vertex_count] = 1;
            if mid.contains(key) { on_boundary[mid.get_or_insert(key, 0)] = 1; }
        }
    };
    var lo: f64[2] = [coords[0], coords[1]];
    var hi: f64[2] = [coords[0], coords[1]];
    for n in node_count() {
        for t in 2 {
            lo[t] = min(lo[t], coords[2 * n + t]);
            hi[t] = max(hi[t], coords[2 * n + t]);
        }
    }
    diameter = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]));
}

// The node graph: each two nodes of a cell are neighbours. CSR in
// `start` and `next`.
fn node_graph(start: i64[>..]&, next: i64[>..]&) {
    let n = node_count();
    var pairs: i64[>..] = [];
    for c in cell_count() {
        for a in per_cell {
            for b in per_cell {
                if a != b { pairs.push(cell_node(c, a) * n + cell_node(c, b)); }
            }
        }
    }
    sort(pairs);
    start.clear();
    next.clear();
    start.push(0);
    var last = -1;
    for e in pairs {
        if e == last { continue; }
        last = e;
        while start.len <= e / n { start.push(next.len); }
        next.push(e % n);
    }
    while start.len <= n { start.push(next.len); }
}

// The breadth-first sequence of the nodes from node `from`, with the
// neighbours of each node in the order of their degrees (Cuthill and
// McKee 1969). When the sequence ends before the last node, it starts
// again at the first node that it did not reach.
fn cuthill_mckee(from: i64, start: i64[:], next: i64[:], seq: i64[>..]&) {
    let n = node_count();
    var seen: u8[>..] = [];
    for i in n { seen.push(0); }
    seq.clear();
    seq.push(from);
    seen[from] = 1;
    var head = 0;
    var around: i64[>..] = [];
    while seq.len < n {
        if head == seq.len {
            // another region: start at its first node
            var k = 0;
            while seen[k] != 0 { k++; }
            seq.push(k);
            seen[k] = 1;
        }
        let v = seq[head];
        head++;
        around.clear();
        for e in start[v]..start[v + 1] {
            if seen[next[e]] == 0 { around.push((start[next[e] + 1] - start[next[e]]) * n + next[e]); }
        }
        sort(around);
        for a in around {
            seen[a % n] = 1;
            seq.push(a % n);
        }
    }
}

// Sets first_dof: the reverse Cuthill-McKee sequence from the last
// node of a sequence from a node of the smallest degree, which is near
// the edge of the mesh (George and Liu 1981).
fn order_nodes() {
    let n = node_count();
    var start: i64[>..] = [];
    var next: i64[>..] = [];
    node_graph(start, next);
    var low = 0;
    for v in n {
        if start[v + 1] - start[v] < start[low + 1] - start[low] { low = v; }
    }
    var seq: i64[>..] = [];
    cuthill_mckee(low, start, next, seq);
    cuthill_mckee(seq[n - 1], start, next, seq);
    first_dof.clear();
    for v in n { first_dof.push(0); }
    var d = 0;
    for i in n {
        let v = seq[n - 1 - i];
        if twin[v] >= 0 { continue; }
        first_dof[v] = d;
        d += if v < vertex_count { vertex_dofs() } else { 2 };
    }
    for v in n {
        if twin[v] >= 0 { first_dof[v] = first_dof[twin[v]]; }
    }
    dof_total = d;
}

// Ties each node of the region to the boundary node at its position minus
// (px, py): the two share their dofs (a `bloch` symmetry with k = 0,
// spec §12). Call it before clear_bcs. False if a node of the region has
// no such node.
fn periodic(px: f64, py: f64) -> bool {
    guard region.len == node_count() && twin.len == node_count() && cell_count() > 0 else { abort("periodic: a mesh and a region"); }
    let tol = tolerance();
    for n in node_count() {
        if region[n] == 0 { continue; }
        let x = coords[2 * n] - px;
        let y = coords[2 * n + 1] - py;
        var found = -1;
        for m in node_count() {
            if m != n && on_boundary[m] != 0 && region[m] == 0 && twin[m] < 0 && abs(coords[2 * m] - x) <= tol && abs(coords[2 * m + 1] - y) <= tol { found = m; }
        }
        guard found >= 0 && (n < vertex_count) == (found < vertex_count) else { return false; }
        twin[n] = found;
    }
    order_nodes();
    build_pattern();
    return true;
}

// The CSR pattern: each dof of a cell couples with each dof of the cell.
// The pairs of dofs are sorted as row n + column.
fn build_pattern() {
    let n = dof_count();
    var pairs: i64[>..] = [];
    var dofs: i64[15] = [0; 15];
    let k = cell_dof_count();
    for c in cell_count() {
        cell_dofs(c, dofs);
        for a in k {
            for b in k { pairs.push(dofs[a] * n + dofs[b]); }
        }
    }
    sort(pairs);
    row_start.clear();
    col.clear();
    row_start.push(0);
    var last = -1;
    for e in pairs {
        if e == last { continue; }
        last = e;
        while row_start.len <= e / n { row_start.push(col.len); }
        col.push((e % n) as! i32);
    }
    while row_start.len <= n { row_start.push(col.len); }
}

// The tolerance of a test of a coordinate: 1e-6 of the diameter (spec §8.2).
fn tolerance() -> f64 { return 1.0e-6 * diameter; }

// Selects the nodes on the boundary where P(x, y) holds. An edge node also
// counts when P holds at both ends of its edge, so that a region on a
// circle keeps the edge nodes on its chords. False if no node holds.
fn select_boundary<P>() -> bool {
    guard region.len == node_count() && cell_count() > 0 else { abort("select_boundary: a mesh"); }
    var found = false;
    for n in node_count() {
        region[n] = if on_boundary[n] != 0 && P(coords[2 * n], coords[2 * n + 1]) { 1 } else { 0 };
    }
    for c in cell_count() {
        for e in per_cell - 3 {
            let m = cell_node(c, 3 + e);
            let a = cell_node(c, EDGES[2 * e]);
            let b = cell_node(c, EDGES[2 * e + 1]);
            if on_boundary[m] != 0 && region[a] != 0 && region[b] != 0 { region[m] = 1; }
        }
    }
    for n in node_count() { if region[n] != 0 { found = true; } }
    return found;
}

// Selects the vertex at (x, y). False if no vertex is within tolerance().
fn select_point(x: f64, y: f64) -> bool {
    guard region.len == node_count() && cell_count() > 0 else { abort("select_point: a mesh"); }
    var best = 0;
    var dist = 1.0e300;
    for n in node_count() {
        region[n] = 0;
        let d = sqrt((coords[2 * n] - x) * (coords[2 * n] - x) + (coords[2 * n + 1] - y) * (coords[2 * n + 1] - y));
        if n < vertex_count && d < dist {
            dist = d;
            best = n;
        }
    }
    guard dist <= tolerance() else { return false; }
    region[best] = 1;
    return true;
}
