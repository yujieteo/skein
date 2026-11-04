import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .refine;
import .tri;

let MAX_LAYER_NODES = 5000;
let MAX_LAYER_COUNT = 1000;

var bl_wx: f64[>..] = [];
var bl_wy: f64[>..] = [];
var bl_nx: f64[>..] = [];
var bl_ny: f64[>..] = [];
var bl_lx: f64[>..] = [];
var bl_ly: f64[>..] = [];
var bl_ln = 0;
var bl_lk = 0;
var bl_ex_edge: i64[>..] = [];
var bl_ex_node: i64[>..] = [];
var bl_ex_at: f64[>..] = [];

fn layers_fail(what: u8[:]) -> bool {
    write_stderr(str("mesh_2d_layers: ", what, "\n"));
    return false;
}

fn area2(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> f64 {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

fn times(d: dictionary<i64, i64>&, key: i64) -> i64 {
    let have = d.get(key);
    if have { return have; }
    return 0;
}

fn read_mesh(xs: f64[>..]&, ys: f64[>..]&, ts: i64[>..]&) {
    xs.clear();
    ys.clear();
    ts.clear();
    for n in point_count() {
        let p = point(n);
        xs.push(p[0]);
        ys.push(p[1]);
    }
    for c in cell_count() { for k in 3 { ts.push(cell_node(c, k)); } }
}

fn count_edges(ts: i64[:], d: dictionary<i64, i64>&) {
    for c in ts.len / 3 {
        for k in 3 { d.update(edge_key(ts[3 * c + k], ts[3 * c + (k + 1) % 3]), 0) { it += 1; }; }
    }
}

fn wall_loop(wall: i64, size: f64) -> bool {
    guard mesh_2d(wall, size) else { return layers_fail("the wall did not mesh"); }
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    var ts: i64[>..] = [];
    read_mesh(xs, ys, ts);
    var count = dictionary<i64, i64> {};
    count_edges(ts, count);
    var succ: i64[>..] = [];
    for n in xs.len { succ.push(-1); }
    var edges = 0;
    var first = -1;
    for c in ts.len / 3 {
        for k in 3 {
            let a = ts[3 * c + k];
            let b = ts[3 * c + (k + 1) % 3];
            if times(count, edge_key(a, b)) != 1 { continue; }
            guard succ[a] < 0 else { return layers_fail("the wall boundary touches itself"); }
            succ[a] = b;
            first = a;
            edges++;
        }
    }
    guard edges >= 3 && edges <= MAX_LAYER_NODES else { return layers_fail("the wall boundary has no 3 to 5000 edges"); }
    bl_wx.clear();
    bl_wy.clear();
    var cur = first;
    for step in edges {
        bl_wx.push(xs[cur]);
        bl_wy.push(ys[cur]);
        cur = succ[cur];
        if cur == first || cur < 0 { break; }
    }
    guard cur == first && bl_wx.len == edges else { return layers_fail("the wall boundary is not one loop"); }
    bl_ln = edges;
    return true;
}

fn wall_normals() -> f64 {
    bl_nx.clear();
    bl_ny.clear();
    var least = 1.0e300;
    for i in bl_ln {
        let p = (i + bl_ln - 1) % bl_ln;
        let q = (i + 1) % bl_ln;
        let ax = bl_wx[i] - bl_wx[p];
        let ay = bl_wy[i] - bl_wy[p];
        let bx = bl_wx[q] - bl_wx[i];
        let by = bl_wy[q] - bl_wy[i];
        let la = sqrt(ax * ax + ay * ay);
        let lb = sqrt(bx * bx + by * by);
        let sx = ay / la + by / lb;
        let sy = -ax / la - bx / lb;
        let s = sqrt(sx * sx + sy * sy);
        bl_nx.push(if s > 1.0e-12 { sx / s } else { 0.0 });
        bl_ny.push(if s > 1.0e-12 { sy / s } else { 0.0 });
        least = min(least, 0.5 * (la + lb));
        if s <= 1.0e-12 { return -1.0; }
    }
    return least;
}

fn layer_positions(first: f64, ratio: f64, count: i64, least: f64) -> i64 {
    var h = first;
    var off = 0.0;
    bl_lx.clear();
    bl_ly.clear();
    for i in bl_ln {
        bl_lx.push(bl_wx[i]);
        bl_ly.push(bl_wy[i]);
    }
    var k = 0;
    for layer in count {
        if h > least { break; }
        off += h;
        for i in bl_ln {
            bl_lx.push(bl_wx[i] + off * bl_nx[i]);
            bl_ly.push(bl_wy[i] + off * bl_ny[i]);
        }
        h *= ratio;
        k++;
    }
    return k;
}

fn lid(k: i64, i: i64) -> i64 { return k * bl_ln + i % bl_ln; }

fn quad_split(k: i64, i: i64) -> i64 {
    let a = lid(k - 1, i);
    let b = lid(k - 1, i + 1);
    let c = lid(k, i + 1);
    let d = lid(k, i);
    let ac = (bl_lx[a] - bl_lx[c]) * (bl_lx[a] - bl_lx[c]) + (bl_ly[a] - bl_ly[c]) * (bl_ly[a] - bl_ly[c]);
    let bd = (bl_lx[b] - bl_lx[d]) * (bl_lx[b] - bl_lx[d]) + (bl_ly[b] - bl_ly[d]) * (bl_ly[b] - bl_ly[d]);
    let acok = area2(bl_lx[a], bl_ly[a], bl_lx[d], bl_ly[d], bl_lx[c], bl_ly[c]) > 0.0 && area2(bl_lx[a], bl_ly[a], bl_lx[c], bl_ly[c], bl_lx[b], bl_ly[b]) > 0.0;
    let bdok = area2(bl_lx[a], bl_ly[a], bl_lx[d], bl_ly[d], bl_lx[b], bl_ly[b]) > 0.0 && area2(bl_lx[d], bl_ly[d], bl_lx[c], bl_ly[c], bl_lx[b], bl_ly[b]) > 0.0;
    if acok && (ac <= bd || !bdok) { return 1; }
    if bdok { return 2; }
    return 0;
}

fn layer_ok(g: i64, k: i64) -> bool {
    for i in bl_ln {
        if phi(g, bl_lx[lid(k, i)], bl_ly[lid(k, i)], 0.0) >= 0.0 { return false; }
        if quad_split(k, i) == 0 { return false; }
    }
    return true;
}

fn loop_simple(k: i64) -> bool {
    for i in bl_ln {
        let a = lid(k, i);
        let b = lid(k, i + 1);
        for j in i + 2..bl_ln {
            if (j + 1) % bl_ln == i { continue; }
            let c = lid(k, j);
            let d = lid(k, j + 1);
            let s1 = area2(bl_lx[a], bl_ly[a], bl_lx[b], bl_ly[b], bl_lx[c], bl_ly[c]) * area2(bl_lx[a], bl_ly[a], bl_lx[b], bl_ly[b], bl_lx[d], bl_ly[d]);
            let s2 = area2(bl_lx[c], bl_ly[c], bl_lx[d], bl_ly[d], bl_lx[a], bl_ly[a]) * area2(bl_lx[c], bl_ly[c], bl_lx[d], bl_ly[d], bl_lx[b], bl_ly[b]);
            if s1 <= 0.0 && s2 <= 0.0 { return false; }
        }
    }
    return true;
}

fn valid_layers(g: i64, k0: i64) -> i64 {
    var k = 0;
    for layer in k0 {
        if !layer_ok(g, k + 1) { break; }
        k++;
    }
    for back in k0 {
        if k == 0 || loop_simple(k) { break; }
        k--;
    }
    return k;
}

fn core_mesh(g: i64, size: f64) -> bool {
    var xy: f64[>..] = [];
    for i in bl_ln {
        xy.push(bl_lx[lid(bl_lk, i)]);
        xy.push(bl_ly[lid(bl_lk, i)]);
    }
    let hole = polygon(xy[..]);
    let before = refines.len;
    for i in bl_ln {
        let a = lid(bl_lk, i);
        let b = lid(bl_lk, i + 1);
        let len = sqrt((bl_lx[a] - bl_lx[b]) * (bl_lx[a] - bl_lx[b]) + (bl_ly[a] - bl_ly[b]) * (bl_ly[a] - bl_ly[b]));
        refine_segment(bl_lx[a], bl_ly[a], bl_lx[b], bl_ly[b], 1.2 * len);
    }
    let done = mesh_2d(subtract(g, hole), size);
    for r in refines.len - before { refines.pop(); }
    return done;
}

fn cell_key(x: f64, y: f64, x0: f64, y0: f64, w: f64) -> i64 {
    return (floor((x - x0) / w) as! i64) * 4294967296 + (floor((y - y0) / w) as! i64);
}

fn outer_node(grid: dictionary<i64, i64>&, x: f64, y: f64, x0: f64, y0: f64, w: f64, tol: f64) -> i64 {
    for dx in 3 {
        for dy in 3 {
            let have = grid.get(cell_key(x + (dx - 1) * w, y + (dy - 1) * w, x0, y0, w));
            if !have { continue; }
            let v = lid(bl_lk, have - 1);
            if abs(bl_lx[v] - x) <= tol && abs(bl_ly[v] - y) <= tol { return v; }
        }
    }
    return -1;
}

fn find_hanging(cx: f64[:], cy: f64[:], onb: bool[:], nid: i64[:], tol: f64) {
    bl_ex_edge.clear();
    bl_ex_node.clear();
    bl_ex_at.clear();
    for j in cx.len {
        if nid[j] >= 0 || !onb[j] { continue; }
        for i in bl_ln {
            let a = lid(bl_lk, i);
            let b = lid(bl_lk, i + 1);
            let ux = bl_lx[b] - bl_lx[a];
            let uy = bl_ly[b] - bl_ly[a];
            let t = ((cx[j] - bl_lx[a]) * ux + (cy[j] - bl_ly[a]) * uy) / (ux * ux + uy * uy);
            let dx = cx[j] - bl_lx[a] - t * ux;
            let dy = cy[j] - bl_ly[a] - t * uy;
            if t > 0.0 && t < 1.0 && sqrt(dx * dx + dy * dy) <= tol {
                bl_ex_edge.push(i);
                bl_ex_node.push(j);
                bl_ex_at.push(t);
                break;
            }
        }
    }
}

fn boundary_nodes(cn: i64, ct: i64[:]) -> bool[>..] {
    var count = dictionary<i64, i64> {};
    count_edges(ct, count);
    var onb: bool[>..] = [];
    for n in cn { onb.push(false); }
    for c in ct.len / 3 {
        for k in 3 {
            let a = ct[3 * c + k];
            let b = ct[3 * c + (k + 1) % 3];
            if times(count, edge_key(a, b)) == 1 {
                onb[a] = true;
                onb[b] = true;
            }
        }
    }
    return onb;
}

fn add_tri(a: i64, b: i64, c: i64) {
    var v: i64[3] = [a, b, c];
    add_simplex(v);
}

fn chain(i: i64, nid: i64[:]) -> i64[>..] {
    var ids: i64[>..] = [];
    var at: f64[>..] = [];
    for e in bl_ex_edge.len {
        if bl_ex_edge[e] != i { continue; }
        var p = ids.len;
        ids.push(nid[bl_ex_node[e]]);
        at.push(bl_ex_at[e]);
        while p > 0 && at[p - 1] > at[p] {
            let t = at[p - 1];
            at[p - 1] = at[p];
            at[p] = t;
            let u = ids[p - 1];
            ids[p - 1] = ids[p];
            ids[p] = u;
            p--;
        }
    }
    return ids;
}

fn emit_quad(k: i64, i: i64, nid: i64[:]) {
    let a = lid(k - 1, i);
    let b = lid(k - 1, i + 1);
    let c = lid(k, i + 1);
    let d = lid(k, i);
    let diag = quad_split(k, i);
    if diag == 1 { add_tri(a, c, b); } else { add_tri(a, d, b); }
    let apex = if diag == 1 { a } else { b };
    if k < bl_lk {
        if diag == 1 { add_tri(a, d, c); } else { add_tri(d, c, b); }
        return;
    }
    var run: i64[>..] = [d];
    run.append(chain(i, nid));
    run.push(c);
    for s in run.len - 1 { add_tri(run[s], run[s + 1], apex); }
}

fn merge_layers(g: i64, tol: f64) {
    var cx: f64[>..] = [];
    var cy: f64[>..] = [];
    var ct: i64[>..] = [];
    read_mesh(cx, cy, ct);
    let b = bounds(g);
    let w = 2.0 * tol;
    var grid = dictionary<i64, i64> {};
    for i in bl_ln { grid.insert(cell_key(bl_lx[lid(bl_lk, i)], bl_ly[lid(bl_lk, i)], b[0], b[1], w), i + 1); }
    var nid: i64[>..] = [];
    for j in cx.len { nid.push(outer_node(grid, cx[j], cy[j], b[0], b[1], w, tol)); }
    let onb = boundary_nodes(cx.len, ct);
    find_hanging(cx, cy, onb, nid, tol);
    clear_mesh(2);
    for p in bl_lx.len { add_point(bl_lx[p], bl_ly[p], 0.0); }
    for j in cx.len { if nid[j] < 0 { nid[j] = add_point(cx[j], cy[j], 0.0); } }
    for c in ct.len / 3 { add_tri(nid[ct[3 * c]], nid[ct[3 * c + 1]], nid[ct[3 * c + 2]]); }
    for k in 1..bl_lk + 1 { for i in bl_ln { emit_quad(k, i, nid); } }
}

fn mesh_2d_layers(g: i64, wall: i64, size: f64, first: f64, ratio: f64, count: i64) -> bool {
    guard size > 0.0 && first > 0.0 && ratio >= 1.0 && count >= 1 && count <= MAX_LAYER_COUNT else {
        return layers_fail("give size > 0, first > 0, ratio >= 1 and 1 to 1000 layers");
    }
    guard wall_loop(wall, size) else { return false; }
    let least = wall_normals();
    guard least > 0.0 else { return layers_fail("the wall has a reversal"); }
    let k0 = layer_positions(first, ratio, count, least);
    bl_lk = valid_layers(g, k0);
    guard bl_lk >= 1 else { return layers_fail("not one layer fits: the wall nodes are too close or g is too narrow"); }
    guard core_mesh(g, size) else { return layers_fail("the core did not mesh"); }
    let b = bounds(g);
    merge_layers(g, 1.0e-9 * sqrt((b[3] - b[0]) * (b[3] - b[0]) + (b[4] - b[1]) * (b[4] - b[1])));
    return true;
}
