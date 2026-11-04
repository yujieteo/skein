import std;
import math;
import dictionary;
import os;
namespace flow;

let EDGES: i64[6] = [0, 1, 1, 2, 2, 0];

var order = 1;
var turbulent = false;
var spalart_on = false;
var thermal = false;
var per_cell = 3;
var vertex_count = 0;
var coords: f64[>..] = [];
var conn: i64[>..] = [];
var on_boundary: u8[>..] = [];
var first_dof: i64[>..] = [];
var twin: i64[>..] = [];
var dof_total = 0;
let MAX_CELL_DOFS = 21;
var part: u8[>..] = [];
var edge_free: u8[>..] = [];
var on_part_boundary: u8[>..] = [];
var in_fluid: u8[>..] = [];
var region_part = 0;
var diameter = 0.0;

var row_start: i64[>..] = [];
var col: i32[>..] = [];

var region: u8[>..] = [];

fn set_turbulent(on: bool) {
    turbulent = on;
    spalart_on = false;
}

fn set_spalart(on: bool) {
    turbulent = on;
    spalart_on = on;
}

fn set_thermal(on: bool) {
    guard !(on && turbulent) else { abort("set_thermal: a flow is thermal or turbulent"); }
    thermal = on;
}

fn node_count() -> i64 { return coords.len / 2; }
fn cell_count() -> i64 { return conn.len / per_cell; }
fn vertex_dofs() -> i64 { return if turbulent { if spalart_on { 4 } else { 5 } } else if thermal { 4 } else { 3 }; }
fn edge_dofs() -> i64 { return if thermal { 3 } else { 2 }; }

fn is_solid_cell(c: i64) -> bool { return part.len > 0 && part[c] != 0; }

fn tdof(n: i64) -> i64 {
    guard thermal && n >= 0 && n < first_dof.len else { abort("tdof: a thermal flow and n < node_count()"); }
    return first_dof[n] + (if n < vertex_count { 3 } else { 2 });
}
fn dof_count() -> i64 { return dof_total; }

fn vdof(n: i64, comp: i64) -> i64 {
    guard n >= 0 && n < first_dof.len && (comp == 0 || comp == 1) else { abort("vdof: n < node_count() and comp 0 or 1"); }
    return first_dof[n] + comp;
}

fn pdof(n: i64) -> i64 {
    guard n >= 0 && n < vertex_count && n < first_dof.len else { abort("pdof: n < vertex_count"); }
    return first_dof[n] + 2;
}

fn turbulence_dof(n: i64, which: i64) -> i64 {
    guard turbulent && n >= 0 && n < vertex_count && n < first_dof.len && (which == 0 || (which == 1 && !spalart_on)) else { abort("turbulence_dof: a turbulent flow, n < vertex_count and which 0 or 1"); }
    return first_dof[n] + 3 + which;
}

fn cell_node(c: i64, l: i64) -> i64 {
    guard c >= 0 && c < cell_count() && l >= 0 && l < per_cell else { abort("cell_node: c < cell_count() and l < per_cell"); }
    return conn[per_cell * c + l];
}

fn cell_dofs(c: i64, dofs: i64[21]&) {
    guard c >= 0 && c < cell_count() else { abort("cell_dofs: c < cell_count()"); }
    for l in per_cell {
        dofs[2 * l] = vdof(cell_node(c, l), 0);
        dofs[2 * l + 1] = vdof(cell_node(c, l), 1);
    }
    for k in 3 { dofs[2 * per_cell + k] = pdof(cell_node(c, k)); }
    if thermal {
        for l in per_cell { dofs[2 * per_cell + 3 + l] = tdof(cell_node(c, l)); }
        return;
    }
    guard turbulent else { return; }
    for k in 3 {
        dofs[2 * per_cell + 3 + k] = turbulence_dof(cell_node(c, k), 0);
        if !spalart_on { dofs[2 * per_cell + 6 + k] = turbulence_dof(cell_node(c, k), 1); }
    }
}

fn cell_dof_count() -> i64 { return 2 * per_cell + (if turbulent { if spalart_on { 6 } else { 9 } } else if thermal { 3 + per_cell } else { 3 }); }

fn mesh_from_arrays(xy: f64[:], tri: i64[:], ord: i64) {
    var none: u8[>..] = [];
    mesh_from_parts(xy, tri, none[..], ord);
}

fn mesh_from_parts(xy: f64[:], tri: i64[:], parts: u8[:], ord: i64) {
    guard (ord == 1 || ord == 2) && (ord == 1 || !turbulent) && xy.len % 2 == 0 && tri.len % 3 == 0 && tri.len > 0 else {
        abort("mesh_from_arrays: order 1 or 2 (1 for a turbulent flow), and triangles");
    }
    guard parts.len == 0 || (thermal && parts.len == tri.len / 3) else { abort("mesh_from_parts: parts of a thermal flow, one for each triangle"); }
    order = ord;
    per_cell = if ord == 1 { 3 } else { 6 };
    vertex_count = xy.len / 2;
    coords.clear();
    coords.append(xy);
    conn.clear();
    part.clear();
    part.append(parts);
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
            count.update(2 * key + (if parts.len > 0 { parts[c] as! i64 } else { 0 }), 0) { it += 1; };
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

fn read_mesh(dir: const u8[:], ord: i64) -> bool {
    var a: u8[>..] = [];
    var b: u8[>..] = [];
    guard read_file(str(dir, "/flow.xy"), a) && read_file(str(dir, "/flow.tri"), b) else { return false; }
    let xy, xy_ok = from_bytes<f64[>..]>(a);
    let tri, tri_ok = from_bytes<i64[>..]>(b);
    guard xy_ok && tri_ok && tri.len > 0 else { return false; }
    var c: u8[>..] = [];
    var parts: u8[>..] = [];
    if thermal && read_file(str(dir, "/flow.part"), c) {
        let p, p_ok = from_bytes<u8[>..]>(c);
        guard p_ok && p.len == tri.len / 3 else { return false; }
        parts.append(p);
    }
    mesh_from_parts(xy, tri, parts[..], ord);
    return true;
}

fn mark_boundary(count: dictionary<i64, i64>&, mid: dictionary<i64, i64>&) {
    on_boundary.clear();
    on_part_boundary.clear();
    in_fluid.clear();
    region.clear();
    for n in node_count() {
        on_boundary.push(0);
        on_part_boundary.push(0);
        in_fluid.push(0);
        region.push(0);
    }
    edge_free.clear();
    for c in cell_count() {
        for k in per_cell { if !is_solid_cell(c) { in_fluid[cell_node(c, k)] = 1; } }
        for e in 3 {
            let a = cell_node(c, EDGES[2 * e]);
            let b = cell_node(c, EDGES[2 * e + 1]);
            let key = 2 * (min(a, b) * vertex_count + max(a, b)) + (if is_solid_cell(c) { 1 } else { 0 });
            edge_free.push(if count.get_or_insert(key, 0) == 1 { 1 } else { 0 });
        }
    }
    count.each() { key, k =>
        if k == 1 {
            let e = key / 2;
            if key % 2 == 0 {
                on_boundary[e / vertex_count] = 1;
                on_boundary[e % vertex_count] = 1;
                if mid.contains(e) { on_boundary[mid.get_or_insert(e, 0)] = 1; }
            } else {
                on_part_boundary[e / vertex_count] = 1;
                on_part_boundary[e % vertex_count] = 1;
                if mid.contains(e) { on_part_boundary[mid.get_or_insert(e, 0)] = 1; }
            }
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
        d += if v < vertex_count { vertex_dofs() } else { edge_dofs() };
    }
    for v in n {
        if twin[v] >= 0 { first_dof[v] = first_dof[twin[v]]; }
    }
    dof_total = d;
}

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

fn build_pattern() {
    let n = dof_count();
    var pairs: i64[>..] = [];
    var dofs: i64[21] = [0; 21];
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

fn tolerance() -> f64 { return 1.0e-6 * diameter; }

fn select_boundary<P>() -> bool {
    guard region.len == node_count() && cell_count() > 0 else { abort("select_boundary: a mesh"); }
    var found = false;
    region_part = 0;
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

fn select_masked<P>(mask: u8[:]) -> bool {
    guard region.len == node_count() && cell_count() > 0 && mask.len == vertex_count else { abort("select_masked: a mesh and a mask of its vertices"); }
    var found = false;
    region_part = 0;
    for n in node_count() {
        region[n] = if n < vertex_count && on_boundary[n] != 0 && mask[n] != 0 && P(coords[2 * n], coords[2 * n + 1]) { 1 } else { 0 };
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

fn select_solid_part<P>() -> bool {
    guard region.len == node_count() && cell_count() > 0 && part.len > 0 else { abort("select_solid_part: a mesh with solid cells"); }
    var found = false;
    region_part = 1;
    for n in node_count() {
        region[n] = if on_part_boundary[n] != 0 && P(coords[2 * n], coords[2 * n + 1]) { 1 } else { 0 };
        if region[n] != 0 { found = true; }
    }
    return found;
}

fn select_point(x: f64, y: f64) -> bool {
    guard region.len == node_count() && cell_count() > 0 else { abort("select_point: a mesh"); }
    region_part = 0;
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
