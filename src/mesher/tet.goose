import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .quality;
import .refine;
import .grid;
import .improve;
import .feature;

let ALPHA_LONG = 0.24999;
let ALPHA_SHORT = 0.41189;

var lattice_cubes: i64[3] = [0, 0, 0];
var lattice_lines: f64[>..] = [];
var line_start: i64[4] = [0, 0, 0, 0];

var lx: f64[>..] = [];
var lphi: f64[>..] = [];
var lnode: i64[>..] = [];
var on_surface: bool[>..] = [];
var cuts = dictionary<i64, i64> {};

fn corner_count() -> i64 {
    return (lattice_cubes[0] + 1) * (lattice_cubes[1] + 1) * (lattice_cubes[2] + 1);
}

fn corner_id(i: i64, j: i64, k: i64) -> i64 {
    return i + (lattice_cubes[0] + 1) * (j + (lattice_cubes[1] + 1) * k);
}

fn centre_id(i: i64, j: i64, k: i64) -> i64 {
    return corner_count() + i + lattice_cubes[0] * (j + lattice_cubes[1] * k);
}

fn lattice_point_count() -> i64 {
    return corner_count() + lattice_cubes[0] * lattice_cubes[1] * lattice_cubes[2];
}

fn half_position(p: i64) -> i64[3] {
    let n = corner_count();
    if p < n {
        let a = lattice_cubes[0] + 1;
        let b = lattice_cubes[1] + 1;
        return [2 * (p % a), 2 * (p / a % b), 2 * (p / (a * b))];
    }
    let q = p - n;
    return [2 * (q % lattice_cubes[0]) + 1, 2 * (q / lattice_cubes[0] % lattice_cubes[1]) + 1,
            2 * (q / (lattice_cubes[0] * lattice_cubes[1])) + 1];
}

fn box_leaves(g: i64, out: f64[>..]&) {
    var todo: i64[>..] = [g];
    while todo.len > 0 {
        let s = shapes[todo[todo.len - 1]];
        todo.pop();
        if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT {
            todo.push(s.b);
            todo.push(s.a);
        } else if s.kind == BOX {
            for i in 6 { out.push(s.p[i]); }
        }
    }
}

fn sign_of(v: f64) -> i64 {
    if v < 0.0 { return -1; }
    if v > 0.0 { return 1; }
    return 0;
}

fn cut_fraction(g: i64, a: i64, b: i64) -> f64 { return cut_fraction_steps(g, a, b, 60); }

fn cut_fraction_steps(g: i64, a: i64, b: i64, steps: i64) -> f64 {
    var lo = 0.0;
    var hi = 1.0;
    for k in steps {
        let t = 0.5 * (lo + hi);
        let x = lx[3 * a] + t * (lx[3 * b] - lx[3 * a]);
        let y = lx[3 * a + 1] + t * (lx[3 * b + 1] - lx[3 * a + 1]);
        let z = lx[3 * a + 2] + t * (lx[3 * b + 2] - lx[3 * a + 2]);
        if phi(g, x, y, z) < 0.0 { lo = t; } else { hi = t; }
    }
    return 0.5 * (lo + hi);
}

fn lattice_node(p: i64) -> i64 {
    if lnode[p] < 0 {
        lnode[p] = add_point(lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        on_surface.push(lphi[p] == 0.0);
    }
    return lnode[p];
}

fn cut_node(g: i64, a: i64, b: i64) -> i64 {
    let key = edge_key(a, b);
    let have = cuts.get(key);
    if have { return have; }
    let t = cut_fraction(g, a, b);
    let node = add_point(lx[3 * a] + t * (lx[3 * b] - lx[3 * a]),
                         lx[3 * a + 1] + t * (lx[3 * b + 1] - lx[3 * a + 1]),
                         lx[3 * a + 2] + t * (lx[3 * b + 2] - lx[3 * a + 2]));
    cuts.insert(key, node);
    on_surface.push(true);
    return node;
}

fn add_tet(a: i64, b: i64, c: i64, d: i64) {
    let pa = vertex(a);
    let u = sub3(vertex(b), pa);
    let v = sub3(vertex(c), pa);
    let w = sub3(vertex(d), pa);
    if dot3(u, cross3(v, w)) < 0.0 {
        add_simplex([a, c, b, d]);
    } else {
        add_simplex([a, b, c, d]);
    }
}

fn add_pyramid(z: i64, q0: i64, q1: i64, q2: i64, q3: i64) {
    let m = min(min(q0, q1), min(q2, q3));
    if m == q0 || m == q2 {
        add_tet(z, q0, q1, q2);
        add_tet(z, q0, q2, q3);
    } else {
        add_tet(z, q1, q2, q3);
        add_tet(z, q1, q3, q0);
    }
}

fn add_prism(a0: i64, a1: i64, a2: i64, b0: i64, b1: i64, b2: i64) {
    var a: i64[3] = [a0, a1, a2];
    var b: i64[3] = [b0, b1, b2];
    let low = min(min(min(a0, a1), min(a2, b0)), min(b1, b2));
    if low == b0 || low == b1 || low == b2 { swap(a, b); }
    var m = 0;
    for i in 3 {
        if a[i] == low { m = i; }
    }
    let r: i64[3] = [a[m], a[(m + 1) % 3], a[(m + 2) % 3]];
    let s: i64[3] = [b[m], b[(m + 1) % 3], b[(m + 2) % 3]];
    add_tet(r[0], s[0], s[1], s[2]);
    if min(r[1], s[2]) < min(r[2], s[1]) {
        add_tet(r[0], r[1], r[2], s[2]);
        add_tet(r[0], r[1], s[2], s[1]);
    } else {
        add_tet(r[0], r[1], r[2], s[1]);
        add_tet(r[0], r[2], s[2], s[1]);
    }
}

fn stuff(g: i64, p: i64[4]) {
    var minus: i64[..4] = [];
    var zero: i64[..4] = [];
    var plus: i64[..4] = [];
    for v in p {
        let s = sign_of(lphi[v]);
        if s < 0 { minus.push(v); } else if s == 0 { zero.push(v); } else { plus.push(v); }
    }
    if minus.len == 0 {
        guard zero.len == 4 else { return; }
        var c: f64[3] = [0.0, 0.0, 0.0];
        for v in p {
            for i in 3 { c[i] += 0.25 * lx[3 * v + i]; }
        }
        guard phi(g, c[0], c[1], c[2]) < 0.0 else { return; }
    }
    if plus.len == 0 {
        add_tet(lattice_node(p[0]), lattice_node(p[1]), lattice_node(p[2]), lattice_node(p[3]));
        return;
    }
    if minus.len == 3 {
        add_prism(lattice_node(minus[0]), lattice_node(minus[1]), lattice_node(minus[2]),
                  cut_node(g, minus[0], plus[0]), cut_node(g, minus[1], plus[0]), cut_node(g, minus[2], plus[0]));
    } else if minus.len == 2 && plus.len == 2 {
        add_prism(lattice_node(minus[0]), cut_node(g, minus[0], plus[0]), cut_node(g, minus[0], plus[1]),
                  lattice_node(minus[1]), cut_node(g, minus[1], plus[0]), cut_node(g, minus[1], plus[1]));
    } else if minus.len == 2 {
        add_pyramid(lattice_node(zero[0]), lattice_node(minus[0]), lattice_node(minus[1]),
                    cut_node(g, minus[1], plus[0]), cut_node(g, minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 1 {
        add_tet(lattice_node(minus[0]), lattice_node(zero[0]), lattice_node(zero[1]),
                cut_node(g, minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 2 {
        add_tet(lattice_node(minus[0]), lattice_node(zero[0]),
                cut_node(g, minus[0], plus[0]), cut_node(g, minus[0], plus[1]));
    } else if minus.len == 1 {
        add_tet(lattice_node(minus[0]), cut_node(g, minus[0], plus[0]),
                cut_node(g, minus[0], plus[1]), cut_node(g, minus[0], plus[2]));
    }
}

var btets: i64[>..] = [];

var edge_index = dictionary<i64, i64> {};
var ea: i64[>..] = [];
var eb: i64[>..] = [];
var elong: bool[>..] = [];
var edge_start: i64[>..] = [];
var edges_at: i64[>..] = [];

fn find_edges() {
    edge_index.clear();
    ea.clear();
    eb.clear();
    elong.clear();
    for t in btets.len / 4 {
        let v = btets[4 * t..4 * t + 4];
        let a, b, longest = longest_edge(lx, v);
        for i in 4 {
            for j in i + 1..4 {
                let e = edge_index.get_or_insert(edge_key(v[i], v[j]), ea.len);
                if e == ea.len {
                    ea.push(v[i]);
                    eb.push(v[j]);
                    elong.push(false);
                }
                if gap(lx, v[i], v[j]) >= longest * (1.0 - 1.0e-9) { elong[e] = true; }
            }
        }
    }
    let points = lx.len / 3;
    edge_start.clear();
    for p in points + 1 { edge_start.push(0); }
    for e in ea.len {
        edge_start[ea[e] + 1]++;
        edge_start[eb[e] + 1]++;
    }
    for p in points { edge_start[p + 1] += edge_start[p]; }
    edges_at.clear();
    for e in 2 * ea.len { edges_at.push(0); }
    var fill: i64[>..] = [];
    fill.append(edge_start[..points]);
    for e in ea.len {
        edges_at[fill[ea[e]]] = e;
        fill[ea[e]]++;
        edges_at[fill[eb[e]]] = e;
        fill[eb[e]]++;
    }
}

let TUBES: f64[4] = [0.3, 0.45, 0.6, 0.8];
let MIN_GAMMA = 0.03;
let SURFACE_BISECTIONS = 26;
let STEP_LOW = 0.25;
let STEP_HIGH = 1.35;
let MAX_TRIES = 8;
let FLAT_GAMMA = 0.01;
let RETRY_SCALES = [1.0, 0.97, 1.03, 0.94, 1.06];
let RELAX_DEPTH = 1;
let RELAX_ROUNDS = 6;
let RELAX = 0.6;
let SUPPORT = 1.05;
let END_REACH = 0.6;

var keep_edges = true;
var curves_found = 0;
var curves_missed = 0;
var chain_a: i64[>..] = [];
var chain_b: i64[>..] = [];
var chain_owner: i64[>..] = [];

fn node_edge_mean(p: i64) -> f64 {
    var sum = 0.0;
    var count = 0.0;
    for e in edges_at[edge_start[p]..edge_start[p + 1]] {
        sum += gap(lx[..], ea[e], eb[e]);
        count += 1.0;
    }
    return sum / count;
}

fn lattice_at(p: i64) -> f64[3] { return [lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]]; }

fn on_lattice_edge(x: f64[3], d: f64) -> bool {
    var planes = 0;
    for a in 3 {
        for i in line_start[a]..line_start[a + 1] {
            if abs(x[a] - lattice_lines[i]) <= 1.0e-9 * d {
                planes++;
                break;
            }
        }
    }
    return planes >= 2;
}

fn surface_slot(g: i64, key: i64, a: i64, b: i64, d: f64, slot_of: dictionary<i64, i64>&, slots: f64[>..]&) -> i64 {
    let have = slot_of.get(key);
    if have { return have; }
    var x = lattice_at(a);
    if a != b {
        let t = cut_fraction_steps(g, a, b, SURFACE_BISECTIONS);
        x = axpy3(x, t, sub3(lattice_at(b), x));
    }
    let eps = FD * d;
    let f0 = phi(g, x[0], x[1], x[2]);
    let n: f64[3] = [(phi(g, x[0] + eps, x[1], x[2]) - f0) / eps, (phi(g, x[0], x[1] + eps, x[2]) - f0) / eps,
                     (phi(g, x[0], x[1], x[2] + eps) - f0) / eps];
    let norm = norm3(n);
    let slot = slots.len / 6;
    slots.append(x);
    slots.append(if norm > 0.0 { [n[0] / norm, n[1] / norm, n[2] / norm] } else { [0.0, 0.0, 0.0] });
    slot_of.insert(key, slot);
    return slot;
}

fn find_curves(g: i64, size: f64) {
    clear_curves();
    let d = shape_scale(g);
    var slot_of = dictionary<i64, i64> {};
    var slots: f64[>..] = [];
    for t in btets.len / 4 {
        let v = btets[4 * t..4 * t + 4];
        var found: i64[10] = [0; 10];
        var n = 0;
        for i in 4 {
            if lphi[v[i]] == 0.0 {
                found[n] = surface_slot(g, -(v[i] + 1), v[i], v[i], d, slot_of, slots);
                n++;
            }
        }
        for i in 4 {
            for j in i + 1..4 {
                guard lphi[v[i]] * lphi[v[j]] < 0.0 else { continue; }
                let inside = if lphi[v[i]] < 0.0 { v[i] } else { v[j] };
                let outside = if lphi[v[i]] < 0.0 { v[j] } else { v[i] };
                found[n] = surface_slot(g, edge_key(v[i], v[j]), inside, outside, d, slot_of, slots);
                n++;
            }
        }
        var best = pair_cos;
        var bi = -1;
        var bj = -1;
        for i in n {
            for j in i + 1..n {
                let c = slots[6 * found[i] + 3] * slots[6 * found[j] + 3] + slots[6 * found[i] + 4] * slots[6 * found[j] + 4] +
                        slots[6 * found[i] + 5] * slots[6 * found[j] + 5];
                if c < best && c > -1.5 {
                    best = c;
                    bi = i;
                    bj = j;
                }
            }
        }
        guard bi >= 0 else { continue; }
        let u: f64[3] = [slots[6 * found[bi]], slots[6 * found[bi] + 1], slots[6 * found[bi] + 2]];
        let nu: f64[3] = [slots[6 * found[bi] + 3], slots[6 * found[bi] + 4], slots[6 * found[bi] + 5]];
        let w: f64[3] = [slots[6 * found[bj]], slots[6 * found[bj] + 1], slots[6 * found[bj] + 2]];
        let nw: f64[3] = [slots[6 * found[bj] + 3], slots[6 * found[bj] + 4], slots[6 * found[bj] + 5]];
        guard abs(dot3(nu, nu) - 1.0) < 1.0e-6 && abs(dot3(nw, nw) - 1.0) < 1.0e-6 else { continue; }
        let mid = axpy3(u, 0.5, sub3(w, u));
        let local = wanted_size(mid[0], mid[1], mid[2], 0.0, size);
        guard !near_curve(mid, local) else { continue; }
        let r = find_ridge(g, u, nu, w, nw, d);
        guard r.found && !on_lattice_edge(r.x, d) && !near_curve(r.x, 0.75 * local / CURVE_DIV) else { continue; }
        trace_curve(g, r, size, d);
    }
}

fn curve_arc(c: i64, arc: f64[>..]&) {
    arc.clear();
    arc.push(0.0);
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = if curve_closed[c] { n } else { n - 1 };
    for i in segs {
        arc.push(arc[i] + norm3(sub3(sample_at(first + (i + 1) % n), sample_at(first + i))));
    }
}

fn curve_foot(c: i64, arc: f64[:], x: f64[3]) -> f64, f64 {
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = arc.len - 1;
    var best = 1.0e300;
    var best_s = 0.0;
    for i in segs {
        let a = sample_at(first + i);
        let ab = sub3(sample_at(first + (i + 1) % n), a);
        let len2 = dot3(ab, ab);
        let lam = if len2 > 0.0 { clamp(dot3(sub3(x, a), ab) / len2, 0.0, 1.0) } else { 0.0 };
        let dist = norm3(sub3(x, axpy3(a, lam, ab)));
        if dist < best {
            best = dist;
            best_s = arc[i] + lam * sqrt(len2);
        }
    }
    return best, best_s;
}

fn curve_exact(g: i64, c: i64, arc: f64[:], s: f64, d: f64) -> f64[3] {
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = arc.len - 1;
    var i = 0;
    while i < segs - 1 && arc[i + 1] < s { i++; }
    let a = sample_at(first + i);
    let b = sample_at(first + (i + 1) % n);
    let len = arc[i + 1] - arc[i];
    let lam = if len > 0.0 { clamp((s - arc[i]) / len, 0.0, 1.0) } else { 0.0 };
    let y = axpy3(a, lam, sub3(b, a));
    guard len > 0.0 else { return y; }
    if !curve_closed[c] && (s <= 0.0 || s >= arc[segs]) { return y; }
    let na = sample_normal(first + i);
    let nb = sample_normal(first + (i + 1) % n);
    let nav = unit3([na[0] + nb[0], na[1] + nb[1], na[2] + nb[2]]);
    let t = unit3(sub3(b, a));
    let w = unit3(cross3(t, nav));
    let delta = 0.2 * len;
    let u, nu, ok_u = settle(g, axpy3(y, -delta, w), d);
    let v, nv, ok_v = settle(g, axpy3(y, delta, w), d);
    guard ok_u && ok_v else { return y; }
    let r = ridge_between(g, u, nu, v, nv, d);
    guard r.found && norm3(sub3(r.x, y)) <= 2.0 * delta else { return y; }
    return r.x;
}

fn node_cost(d: f64, l: f64) -> f64 { return (d / l) * (d / l); }

fn step_cost(advance: f64, le: f64, l: f64) -> f64 { return (advance - le) * (advance - le) / (l * l); }

fn find_candidates(c: i64, arc: f64[:], tube: f64, node: i64[>..]&, at_s: f64[>..]&, at_d: f64[>..]&, at_l: f64[>..]&) {
    node.clear();
    at_s.clear();
    at_d.clear();
    at_l.clear();
    let start = sample_at(curve_start[c]);
    var near = 0;
    var near_d = 1.0e300;
    for p in lx.len / 3 {
        let dist = norm3(sub3(lattice_at(p), start));
        if dist < near_d {
            near_d = dist;
            near = p;
        }
    }
    var seen: bool[>..] = [];
    for p in lx.len / 3 { seen.push(false); }
    var queue: i64[>..] = [near];
    seen[near] = true;
    var head = 0;
    while head < queue.len {
        let p = queue[head];
        head++;
        let l = node_edge_mean(p);
        let dist, s = curve_foot(c, arc, lattice_at(p));
        guard dist <= 1.2 * l else { continue; }
        if dist <= tube * l && chain_owner[p] == -1 {
            node.push(p);
            at_s.push(s);
            at_d.push(dist);
            at_l.push(l);
        }
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            if !seen[q] {
                seen[q] = true;
                queue.push(q);
            }
        }
    }
}

var tet_start: i64[>..] = [];
var tets_of: i64[>..] = [];

fn find_tets_of() {
    let points = lx.len / 3;
    tet_start.clear();
    for p in points + 1 { tet_start.push(0); }
    for t in btets.len / 4 {
        for i in 4 { tet_start[btets[4 * t + i] + 1]++; }
    }
    for p in points { tet_start[p + 1] += tet_start[p]; }
    tets_of.clear();
    for i in tet_start[points] { tets_of.push(0); }
    var fill: i64[>..] = [];
    fill.append(tet_start[..points]);
    for t in btets.len / 4 {
        for i in 4 {
            let p = btets[4 * t + i];
            tets_of[fill[p]] = t;
            fill[p]++;
        }
    }
}

fn stays_inside(m: i64) -> bool {
    guard lphi[m] < 0.0 else { return false; }
    for e in edges_at[edge_start[m]..edge_start[m + 1]] {
        let q = if ea[e] == m { eb[e] } else { ea[e] };
        guard lphi[q] > 0.0 else { continue; }
        let t = -lphi[m] / (lphi[q] - lphi[m]);
        if t < SUPPORT * (if elong[e] { ALPHA_LONG } else { ALPHA_SHORT }) { return false; }
    }
    return true;
}

fn edge_supported(a: i64, b: i64) -> bool {
    for t in tets_of[tet_start[a]..tet_start[a + 1]] {
        var has_b = false;
        for i in 4 {
            if btets[4 * t + i] == b { has_b = true; }
        }
        guard has_b else { continue; }
        for i in 4 {
            let m = btets[4 * t + i];
            if m != a && m != b && stays_inside(m) { return true; }
        }
    }
    return false;
}

var dp_pos: f64[>..] = [];
var dp_order: i64[>..] = [];
var dp_cost: f64[>..] = [];
var dp_prev: i64[>..] = [];

fn chain_start(count: i64, closed: bool, p0: i64, total: f64, at_s: f64[:], at_d: f64[:], at_l: f64[:]) {
    dp_pos.clear();
    dp_order.clear();
    dp_cost.clear();
    dp_prev.clear();
    var keys: i64[>..] = [];
    for k in count {
        let origin = if closed { at_s[p0] } else { 0.0 };
        var x = at_s[k] - origin;
        if x < 0.0 { x += total; }
        dp_pos.push(x);
        dp_cost.push(1.0e300);
        dp_prev.push(-1);
        keys.push(trunc(x / total * 1.0e9) as i64 * 4194304 + k);
    }
    sort(keys);
    for k in count { dp_order.push(keys[k] % 4194304); }
    if closed {
        dp_cost[p0] = node_cost(at_d[p0], at_l[p0]);
    } else {
        for k in count {
            if at_s[k] <= END_REACH * at_l[k] { dp_cost[k] = node_cost(at_d[k], at_l[k]); }
        }
    }
}

fn chain_forward(count: i64, closed: bool, p0: i64, node: i64[:], index_of: dictionary<i64, i64>&, at_d: f64[:], at_l: f64[:]) {
    for oi in count {
        let q = dp_order[oi];
        for e in edges_at[edge_start[node[q]]..edge_start[node[q] + 1]] {
            let pn = if ea[e] == node[q] { eb[e] } else { ea[e] };
            let have = index_of.get(pn);
            guard have else { continue; }
            let pk = have;
            guard dp_cost[pk] < 1.0e299 && dp_pos[pk] < dp_pos[q] else { continue; }
            let le = gap(lx[..], node[q], pn);
            let advance = dp_pos[q] - dp_pos[pk];
            guard advance >= STEP_LOW * le && advance <= STEP_HIGH * le else { continue; }
            if closed && q == p0 { continue; }
            guard edge_supported(node[pk], node[q]) else { continue; }
            let l = 0.5 * (at_l[q] + at_l[pk]);
            let sum = dp_cost[pk] + step_cost(advance, le, l) + node_cost(at_d[q], at_l[q]);
            if sum < dp_cost[q] {
                dp_cost[q] = sum;
                dp_prev[q] = pk;
            }
        }
    }
}

fn best_chain(c: i64, total: f64, node: i64[:], at_s: f64[:], at_d: f64[:], at_l: f64[:], path: i64[>..]&) -> bool {
    path.clear();
    let count = node.len;
    guard count >= 2 else { return false; }
    let closed = curve_closed[c];
    var index_of = dictionary<i64, i64> {};
    for k in count { index_of.insert(node[k], k); }
    var starts: i64[>..] = [];
    if closed {
        for k in count {
            if at_s[k] < 2.5 * at_l[k] || at_s[k] > total - 0.5 * at_l[k] { starts.push(k); }
        }
        guard starts.len > 0 else { return false; }
    } else {
        starts.push(-1);
    }
    var best_total = 1.0e300;
    var best_prev: i64[>..] = [];
    var best_last = -1;
    for p0 in starts {
        chain_start(count, closed, p0, total, at_s, at_d, at_l);
        chain_forward(count, closed, p0, node, index_of, at_d, at_l);
        for q in count {
            guard dp_cost[q] < 1.0e299 else { continue; }
            var sum = dp_cost[q];
            if closed {
                guard q != p0 else { continue; }
                var closes = false;
                for e in edges_at[edge_start[node[q]]..edge_start[node[q] + 1]] {
                    if ea[e] == node[p0] || eb[e] == node[p0] { closes = true; }
                }
                let le = gap(lx[..], node[q], node[p0]);
                let advance = total - dp_pos[q];
                guard closes && advance >= STEP_LOW * le && advance <= STEP_HIGH * le && edge_supported(node[q], node[p0]) else { continue; }
                sum += step_cost(advance, le, 0.5 * (at_l[q] + at_l[p0]));
            } else {
                guard total - at_s[q] <= END_REACH * at_l[q] else { continue; }
            }
            if sum < best_total {
                best_total = sum;
                best_last = q;
                best_prev.clear();
                best_prev.append(dp_prev);
            }
        }
    }
    guard best_last >= 0 else { return false; }
    var rev: i64[>..] = [];
    var k = best_last;
    while k >= 0 {
        rev.push(k);
        k = best_prev[k];
    }
    for i in rev.len { path.push(rev[rev.len - 1 - i]); }
    return true;
}

fn pinned(p: i64, planes: f64[:]) -> bool {
    if lphi[p] == 0.0 || chain_owner[p] >= 0 { return true; }
    for k in planes.len / 2 {
        if abs(lx[3 * p + planes[2 * k] as i64] - planes[2 * k + 1]) <= 1.0e-12 { return true; }
    }
    return false;
}

fn near_chain(node: i64[:], path: i64[:], planes: f64[:], free: i64[>..]&) {
    var queue: i64[>..] = [];
    var depth: i64[>..] = [];
    var seen = dictionary<i64, i64> {};
    for k in path.len {
        queue.push(node[path[k]]);
        depth.push(0);
        seen.insert(node[path[k]], 0);
    }
    var head = 0;
    while head < queue.len {
        let p = queue[head];
        let dp = depth[head];
        head++;
        guard dp < RELAX_DEPTH else { continue; }
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            guard !seen.get(q) else { continue; }
            seen.insert(q, dp + 1);
            queue.push(q);
            depth.push(dp + 1);
            if !pinned(q, planes) { free.push(q); }
        }
    }
}

fn relax_free(g: i64, free: i64[:], free_old: f64[:], free_phi: f64[:]) {
    for round in RELAX_ROUNDS {
        for p in free {
            var mean: f64[3] = [0.0, 0.0, 0.0];
            var count = 0.0;
            for e in edges_at[edge_start[p]..edge_start[p + 1]] {
                let q = if ea[e] == p { eb[e] } else { ea[e] };
                for a in 3 { mean[a] += lx[3 * q + a]; }
                count += 1.0;
            }
            for a in 3 { lx[3 * p + a] += RELAX * (mean[a] / count - lx[3 * p + a]); }
        }
    }
    for k in free.len {
        let p = free[k];
        let f = phi(g, lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        if f == 0.0 || (f < 0.0) != (free_phi[k] < 0.0) {
            for a in 3 { lx[3 * p + a] = free_old[3 * k + a]; }
        } else {
            lphi[p] = f;
        }
    }
}

fn cells_valid(cells: i64[:], sign: f64[:], marked: dictionary<i64, i64>&, target: f64[:], old: f64[:]) -> bool {
    var ok = true;
    for i in cells.len {
        let v = btets[4 * cells[i]..4 * cells[i] + 4];
        let q = if sign[i] < 0.0 { [v[1], v[0], v[2], v[3]] } else { [v[0], v[1], v[2], v[3]] };
        guard gamma_of(lattice_at(q[0]), lattice_at(q[1]), lattice_at(q[2]), lattice_at(q[3])) < MIN_GAMMA else { continue; }
        ok = false;
        var far = -1;
        var far_move = -1.0;
        for n in v {
            let k = marked.get(n);
            if k {
                let move = norm3(sub3([target[3 * k], target[3 * k + 1], target[3 * k + 2]],
                                      [old[3 * k], old[3 * k + 1], old[3 * k + 2]]));
                if move > far_move {
                    far_move = move;
                    far = n;
                }
            }
        }
        if far >= 0 { chain_owner[far] = -2; }
    }
    return ok;
}

fn snap_chain(g: i64, c: i64, arc: f64[:], d: f64, node: i64[:], at_s: f64[:], path: i64[:], planes: f64[:]) -> bool {
    let total = arc[arc.len - 1];
    var target: f64[>..] = [];
    for k in path.len {
        var s = at_s[path[k]];
        if !curve_closed[c] {
            if k == 0 { s = 0.0; }
            if k == path.len - 1 { s = total; }
        }
        target.append(curve_exact(g, c, arc, s, d));
    }
    var marked = dictionary<i64, i64> {};
    for k in path.len { marked.insert(node[path[k]], k); }
    var free: i64[>..] = [];
    near_chain(node, path, planes, free);
    var cells: i64[>..] = [];
    var sign: f64[>..] = [];
    var moves = dictionary<i64, i64> {};
    for p in free { moves.insert(p, 0); }
    for t in btets.len / 4 {
        var touched = false;
        for i in 4 {
            if marked.get(btets[4 * t + i]) || moves.get(btets[4 * t + i]) { touched = true; }
        }
        guard touched else { continue; }
        cells.push(t);
        sign.push(if lattice_volume(t) < 0.0 { -1.0 } else { 1.0 });
    }
    var old: f64[>..] = [];
    var old_phi: f64[>..] = [];
    for k in path.len {
        let p = node[path[k]];
        old.append(lattice_at(p));
        old_phi.push(lphi[p]);
        for a in 3 { lx[3 * p + a] = target[3 * k + a]; }
        lphi[p] = 0.0;
    }
    var free_old: f64[>..] = [];
    var free_phi: f64[>..] = [];
    for p in free {
        free_old.append(lattice_at(p));
        free_phi.push(lphi[p]);
    }
    relax_free(g, free[..], free_old[..], free_phi[..]);
    let ok = cells_valid(cells[..], sign[..], marked, target[..], old[..]);
    if !ok {
        for k in path.len {
            let p = node[path[k]];
            for a in 3 { lx[3 * p + a] = old[3 * k + a]; }
            lphi[p] = old_phi[k];
        }
        for k in free.len {
            for a in 3 { lx[3 * free[k] + a] = free_old[3 * k + a]; }
            lphi[free[k]] = free_phi[k];
        }
    }
    return ok;
}

fn lattice_volume(t: i64) -> f64 {
    let v = btets[4 * t..4 * t + 4];
    let a = lattice_at(v[0]);
    return dot3(sub3(lattice_at(v[1]), a), cross3(sub3(lattice_at(v[2]), a), sub3(lattice_at(v[3]), a)));
}

fn protect_edges(g: i64, size: f64) -> i64 {
    curves_found = 0;
    curves_missed = 0;
    chain_a.clear();
    chain_b.clear();
    find_curves(g, size);
    chain_owner.clear();
    for p in lx.len / 3 { chain_owner.push(-1); }
    guard curve_count() > 0 else { return 0; }
    find_tets_of();
    let d = shape_scale(g);
    var missed = 0;
    var arc: f64[>..] = [];
    var node: i64[>..] = [];
    var at_s: f64[>..] = [];
    var at_d: f64[>..] = [];
    var at_l: f64[>..] = [];
    var path: i64[>..] = [];
    var boxes: f64[>..] = [];
    box_leaves(g, boxes);
    var planes: f64[>..] = [];
    for k in boxes.len / 6 {
        for a in 3 {
            planes.append([a, boxes[6 * k + a], a, boxes[6 * k + 3 + a]]);
        }
    }
    for c in curve_count() {
        curve_arc(c, arc);
        var done = false;
        for tube in TUBES {
            for attempt in MAX_TRIES {
                find_candidates(c, arc[..], tube, node, at_s, at_d, at_l);
                guard best_chain(c, arc[arc.len - 1], node[..], at_s[..], at_d[..], at_l[..], path) else { break; }
                if snap_chain(g, c, arc[..], d, node[..], at_s[..], path[..], planes[..]) {
                    for k in path { chain_owner[node[k]] = c; }
                    for k in path.len - 1 {
                        chain_a.push(node[path[k]]);
                        chain_b.push(node[path[k + 1]]);
                    }
                    if curve_closed[c] {
                        chain_a.push(node[path[path.len - 1]]);
                        chain_b.push(node[path[0]]);
                    }
                    done = true;
                    break;
                }
            }
            if done { break; }
        }
        for p in chain_owner.len {
            if chain_owner[p] == -2 { chain_owner[p] = -1; }
        }
        if !done { missed++; }
    }
    curves_found = curve_count();
    curves_missed = missed;
    return missed;
}

fn warp(g: i64) {
    for p in lx.len / 3 {
        guard lphi[p] != 0.0 else { continue; }
        var best = 2.0;
        var best_x: f64[3] = [0.0, 0.0, 0.0];
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            guard lphi[q] != 0.0 && (lphi[q] < 0.0) != (lphi[p] < 0.0) else { continue; }
            let inside = if lphi[p] < 0.0 { p } else { q };
            let outside = if lphi[p] < 0.0 { q } else { p };
            let t = cut_fraction(g, inside, outside);
            let from_p = if inside == p { t } else { 1.0 - t };
            let alpha = if elong[e] { ALPHA_LONG } else { ALPHA_SHORT };
            if from_p < alpha && from_p < best {
                best = from_p;
                for c in 3 { best_x[c] = lx[3 * p + c] + from_p * (lx[3 * q + c] - lx[3 * p + c]); }
            }
        }
        if best < 2.0 {
            for c in 3 { lx[3 * p + c] = best_x[c]; }
            lphi[p] = 0.0;
        }
    }
}

fn lattice_tets() {
    btets.clear();
    let n = lattice_cubes;
    for k in n[2] {
        for j in n[1] {
            for i in n[0] {
                let c = centre_id(i, j, k);
                for axis in 3 {
                    var o: i64[3] = [i, j, k];
                    o[axis] += 1;
                    guard o[axis] < n[axis] else { continue; }
                    let c2 = centre_id(o[0], o[1], o[2]);
                    let u = (axis + 1) % 3;
                    let v = (axis + 2) % 3;
                    var f: i64[4] = [0, 0, 0, 0];
                    let du: i64[4] = [0, 1, 1, 0];
                    let dv: i64[4] = [0, 0, 1, 1];
                    for q in 4 {
                        var at: i64[3] = o;
                        at[u] += du[q];
                        at[v] += dv[q];
                        f[q] = corner_id(at[0], at[1], at[2]);
                    }
                    for q in 4 { btets.append([c, c2, f[q], f[(q + 1) % 4]]); }
                }
            }
        }
    }
}

fn build_lattice(g: i64, b: f64[6], size: f64) -> bool {
    var boxes: f64[>..] = [];
    box_leaves(g, boxes);
    lattice_lines.clear();
    var line: f64[>..] = [];
    for a in 3 {
        guard b[3 + a] > b[a] else { return false; }
        grid_lines(boxes[..], a, b[a], b[3 + a], size, line);
        line_start[a] = lattice_lines.len;
        lattice_lines.push(2.0 * line[0] - line[1]);
        lattice_lines.append(line);
        lattice_lines.push(2.0 * line[line.len - 1] - line[line.len - 2]);
        lattice_cubes[a] = line.len + 1;
    }
    line_start[3] = lattice_lines.len;
    lx.clear();
    lphi.clear();
    lnode.clear();
    cuts.clear();
    on_surface.clear();
    for p in lattice_point_count() {
        let at = half_position(p);
        for a in 3 {
            let lo = lattice_lines[line_start[a] + at[a] / 2];
            lx.push(if at[a] % 2 == 0 { lo } else { 0.5 * (lo + lattice_lines[line_start[a] + at[a] / 2 + 1]) });
        }
    }
    return true;
}

fn has_bad_cell() -> bool {
    for c in cell_count() {
        let a = vertex(cell_node(c, 0));
        let b = vertex(cell_node(c, 1));
        let d = vertex(cell_node(c, 2));
        let e = vertex(cell_node(c, 3));
        if gamma_of(a, b, d, e) < FLAT_GAMMA { return true; }
    }
    return false;
}

fn lattice_pass(g: i64, b: f64[6], scale: f64) -> bool {
    clear_mesh(3);
    guard build_lattice(g, b, scale) else { return false; }
    lattice_tets();
    guard refine_cells(lx, btets, 4, scale) else { return false; }
    var h = scale;
    for a in 3 {
        for i in lattice_cubes[a] { h = min(h, lattice_lines[line_start[a] + i + 1] - lattice_lines[line_start[a] + i]); }
    }
    for t in targets { h = min(h, t.size); }
    lphi.clear();
    lnode.clear();
    for p in lx.len / 3 {
        let f = phi(g, lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        lphi.push(if abs(f) <= 1.0e-9 * h { 0.0 } else { f });
        lnode.push(-1);
    }
    find_edges();
    return true;
}

fn chains_in_mesh() -> bool {
    let n = point_count();
    var want = dictionary<i64, i64> {};
    for k in chain_a.len {
        let a = lnode[chain_a[k]];
        let b = lnode[chain_b[k]];
        guard a >= 0 && b >= 0 else { return false; }
        want.insert(min(a, b) * n + max(a, b), 0);
    }
    for c in cell_count() {
        for i in 4 {
            for j in i + 1..4 {
                let a = cell_node(c, i);
                let b = cell_node(c, j);
                let key = min(a, b) * n + max(a, b);
                if want.get(key) { want.update(key, 0) { it = 1; }; }
            }
        }
    }
    var missing = 0;
    want.each() { key, v => if v == 0 { missing++; } };
    return missing == 0;
}

fn mesh_3d_body(g: i64, size: f64) -> bool {
    clear_mesh(3);
    curves_found = 0;
    curves_missed = 0;
    if shell_tets(g, size) { return true; }
    let b = bounds(g);
    if has_round(g) && targets_span_z(b[2], b[5]) {
        let s = section(g, b[2], b[5]);
        if s >= 0 { return sweep_tets(s, b[2], b[5], size); }
    }
    for shift in [1, 2] {
        guard has_round(g) && turnable() else { break; }
        let count = shapes.len;
        turn_shapes(count, shift);
        let t = bounds(g);
        let s = if targets_span_z(t[2], t[5]) { section(g, t[2], t[5]) } else { -1 };
        let done = s >= 0 && sweep_tets(s, t[2], t[5], size);
        turn_shapes(count, 3 - shift);
        if s >= 0 {
            turn_nodes(3 - shift);
            return done;
        }
    }
    let passes = RETRY_SCALES.len + 1;
    for pass in passes {
        let protect = keep_edges && pass < RETRY_SCALES.len;
        let scale = if protect { size * RETRY_SCALES[pass] } else { size };
        guard lattice_pass(g, b, scale) else { return false; }
        if protect {
            let missed = protect_edges(g, scale);
            if missed > 0 && (pass + 1 < RETRY_SCALES.len || missed == curves_found) { continue; }
        }
        warp(g);
        for t in btets.len / 4 { stuff(g, [btets[4 * t], btets[4 * t + 1], btets[4 * t + 2], btets[4 * t + 3]]); }
        if protect && curves_found > 0 && (has_bad_cell() || !chains_in_mesh()) {
            curves_missed = curves_found;
            continue;
        }
        if improve_cells { improve(on_surface[..]); }
        if protect && curves_found > 0 && !quality_passes(mesh_quality()) {
            curves_missed = curves_found;
            continue;
        }
        break;
    }
    return true;
}

let CUT_SNAP = 0.12;
var cut_px: f64[>..] = [];
var cut_cells: i64[>..] = [];
var cut_next: i64[>..] = [];
var cut_edge = dictionary<i64, i64> {};
var cut_sign: i64[>..] = [];
var cut_h = 0;

fn cut_add(a: i64, b: i64, c: i64, d: i64) {
    let pa: f64[3] = [cut_px[3 * a], cut_px[3 * a + 1], cut_px[3 * a + 2]];
    let u = sub3([cut_px[3 * b], cut_px[3 * b + 1], cut_px[3 * b + 2]], pa);
    let v = sub3([cut_px[3 * c], cut_px[3 * c + 1], cut_px[3 * c + 2]], pa);
    let w = sub3([cut_px[3 * d], cut_px[3 * d + 1], cut_px[3 * d + 2]], pa);
    if dot3(u, cross3(v, w)) < 0.0 {
        cut_next.append([a, c, b, d]);
    } else {
        cut_next.append([a, b, c, d]);
    }
}

fn cut_pyramid(z: i64, q0: i64, q1: i64, q2: i64, q3: i64) {
    let m = min(min(q0, q1), min(q2, q3));
    if m == q0 || m == q2 {
        cut_add(z, q0, q1, q2);
        cut_add(z, q0, q2, q3);
    } else {
        cut_add(z, q1, q2, q3);
        cut_add(z, q1, q3, q0);
    }
}

fn cut_prism(a0: i64, a1: i64, a2: i64, b0: i64, b1: i64, b2: i64) {
    var a: i64[3] = [a0, a1, a2];
    var b: i64[3] = [b0, b1, b2];
    let low = min(min(min(a0, a1), min(a2, b0)), min(b1, b2));
    if low == b0 || low == b1 || low == b2 { swap(a, b); }
    var m = 0;
    for i in 3 {
        if a[i] == low { m = i; }
    }
    let r: i64[3] = [a[m], a[(m + 1) % 3], a[(m + 2) % 3]];
    let s: i64[3] = [b[m], b[(m + 1) % 3], b[(m + 2) % 3]];
    cut_add(r[0], s[0], s[1], s[2]);
    if min(r[1], s[2]) < min(r[2], s[1]) {
        cut_add(r[0], r[1], r[2], s[2]);
        cut_add(r[0], r[1], s[2], s[1]);
    } else {
        cut_add(r[0], r[1], r[2], s[1]);
        cut_add(r[0], r[2], s[2], s[1]);
    }
}

fn cut_point(a: i64, b: i64) -> i64 {
    let key = min(a, b) * 1073741824 + max(a, b);
    let have = cut_edge.get(key);
    if have { return have; }
    var lo = 0.0;
    var hi = 1.0;
    for k in 50 {
        let t = 0.5 * (lo + hi);
        let x = cut_px[3 * a] + t * (cut_px[3 * b] - cut_px[3 * a]);
        let y = cut_px[3 * a + 1] + t * (cut_px[3 * b + 1] - cut_px[3 * a + 1]);
        let z = cut_px[3 * a + 2] + t * (cut_px[3 * b + 2] - cut_px[3 * a + 2]);
        if phi(cut_h, x, y, z) < 0.0 { lo = t; } else { hi = t; }
    }
    let t = 0.5 * (lo + hi);
    let n = cut_px.len / 3;
    for i in 3 { cut_px.push(cut_px[3 * a + i] + t * (cut_px[3 * b + i] - cut_px[3 * a + i])); }
    cut_sign.push(0);
    cut_edge.insert(key, n);
    return n;
}

fn cut_side(p: i64[4], sgn: i64[4]) {
    var minus: i64[..4] = [];
    var zero: i64[..4] = [];
    var plus: i64[..4] = [];
    for i in 4 {
        if sgn[i] < 0 { minus.push(p[i]); } else if sgn[i] == 0 { zero.push(p[i]); } else { plus.push(p[i]); }
    }
    guard minus.len > 0 else { return; }
    if plus.len == 0 {
        cut_add(p[0], p[1], p[2], p[3]);
    } else if minus.len == 3 {
        cut_prism(minus[0], minus[1], minus[2], cut_point(minus[0], plus[0]), cut_point(minus[1], plus[0]), cut_point(minus[2], plus[0]));
    } else if minus.len == 2 && plus.len == 2 {
        cut_prism(minus[0], cut_point(minus[0], plus[0]), cut_point(minus[0], plus[1]),
                  minus[1], cut_point(minus[1], plus[0]), cut_point(minus[1], plus[1]));
    } else if minus.len == 2 {
        cut_pyramid(zero[0], minus[0], minus[1], cut_point(minus[1], plus[0]), cut_point(minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 1 {
        cut_add(minus[0], zero[0], zero[1], cut_point(minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 2 {
        cut_add(minus[0], zero[0], cut_point(minus[0], plus[0]), cut_point(minus[0], plus[1]));
    } else if minus.len == 1 {
        cut_add(minus[0], cut_point(minus[0], plus[0]), cut_point(minus[0], plus[1]), cut_point(minus[0], plus[2]));
    }
}

fn cut_interface(h: i64) {
    cut_h = h;
    let n = cut_px.len / 3;
    var total: f64[>..] = [];
    var count: f64[>..] = [];
    var psi: f64[>..] = [];
    for i in n {
        total.push(0.0);
        count.push(0.0);
        psi.push(phi(h, cut_px[3 * i], cut_px[3 * i + 1], cut_px[3 * i + 2]));
    }
    for c in cut_cells.len / 4 {
        for i in 4 {
            for j in i + 1..4 {
                let a = cut_cells[4 * c + i];
                let b = cut_cells[4 * c + j];
                var l = 0.0;
                for k in 3 { l += (cut_px[3 * a + k] - cut_px[3 * b + k]) * (cut_px[3 * a + k] - cut_px[3 * b + k]); }
                l = sqrt(l);
                total[a] += l; count[a] += 1.0;
                total[b] += l; count[b] += 1.0;
            }
        }
    }
    cut_sign.clear();
    for i in n {
        let mean = if count[i] > 0.0 { total[i] / count[i] } else { 1.0 };
        cut_sign.push(if abs(psi[i]) <= CUT_SNAP * mean { 0 } else { sign_of(psi[i]) });
    }
    cut_edge.clear();
    cut_next.clear();
    for c in cut_cells.len / 4 {
        let p: i64[4] = [cut_cells[4 * c], cut_cells[4 * c + 1], cut_cells[4 * c + 2], cut_cells[4 * c + 3]];
        let s: i64[4] = [cut_sign[p[0]], cut_sign[p[1]], cut_sign[p[2]], cut_sign[p[3]]];
        if s[0] <= 0 && s[1] <= 0 && s[2] <= 0 && s[3] <= 0 || s[0] >= 0 && s[1] >= 0 && s[2] >= 0 && s[3] >= 0 {
            cut_next.append(p);
            continue;
        }
        cut_side(p, s);
        cut_side(p, [-s[0], -s[1], -s[2], -s[3]]);
    }
    cut_cells.clear();
    cut_cells.append(cut_next);
}

fn mesh_3d(g: i64, size: f64) -> bool {
    guard mesh_3d_body(g, size) else { return false; }
    guard interfaces.len > 0 else { return true; }
    cut_px.clear();
    cut_cells.clear();
    for n in point_count() {
        let x = point(n);
        for i in 3 { cut_px.push(x[i]); }
    }
    for c in cell_count() {
        for i in 4 { cut_cells.push(cell_node(c, i)); }
    }
    for h in interfaces { cut_interface(h); }
    clear_mesh(3);
    for n in cut_px.len / 3 { add_point(cut_px[3 * n], cut_px[3 * n + 1], cut_px[3 * n + 2]); }
    for c in cut_cells.len / 4 { add_simplex([cut_cells[4 * c], cut_cells[4 * c + 1], cut_cells[4 * c + 2], cut_cells[4 * c + 3]]); }
    return true;
}
