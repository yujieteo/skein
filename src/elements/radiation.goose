import std;
import math;
import os;
import dictionary;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import solvers from "../solvers";

let SUB = 4;
let RAD_CLOSURE = 1.0e-2;
let CLIP_TOLERANCE = 1.0e-12;
let SINKHORN_SWEEPS = 100;
let SINKHORN_TOLERANCE = 1.0e-14;
let PIVOT_FLOOR = 1.0e-300;
let RAD_STEPS = 200;
let RAD_TOLERANCE = 1.0e-10;

var rad_facet: i64[>..] = [];
var rad_eps: f64[>..] = [];
var rad_xy: f64[>..] = [];
var rad_q: f64[>..] = [];
var sigma_sb = 0.0;
var rad_steps = 0;

fn rad_count() -> i64 { return rad_eps.len; }

fn clear_radiation() {
    rad_facet.clear();
    rad_eps.clear();
    rad_xy.clear();
    rad_q.clear();
}

fn add_radiation(eps: f64, sb: f64) {
    guard physics == HEAT && dim == 2 && eps > 0.0 && eps <= 1.0 && sb > 0.0 else { abort("add_radiation: a 2d heat physics, 0 < eps <= 1 and sb > 0"); }
    sigma_sb = sb;
    var count = dictionary<i64, i64> {};
    for c in cell_count() {
        for f in 3 {
            let a = cell_node(c, facet_node(2, f, 0));
            let b = cell_node(c, facet_node(2, f, 1));
            count.update(min(a, b) * node_count() + max(a, b), 0) { it += 1; };
        }
    }
    for c in cell_count() {
        for f in 3 {
            let a = cell_node(c, facet_node(2, f, 0));
            let b = cell_node(c, facet_node(2, f, 1));
            if count.get_or_insert(min(a, b) * node_count() + max(a, b), 0) != 1 { continue; }
            var on = false;
            for k in region.len { if on_part(a, k) && on_part(b, k) { on = true; } }
            if on { add_facet(c, f, eps); }
        }
    }
}

fn add_facet(c: i64, f: i64, eps: f64) {
    let a = cell_node(c, facet_node(2, f, 0));
    let b = cell_node(c, facet_node(2, f, 1));
    let o = cell_node(c, (f + 2) % 3);
    rad_facet.push(a);
    rad_facet.push(b);
    rad_facet.push(if order == 2 { cell_node(c, facet_node(2, f, 2)) } else { -1 });
    rad_eps.push(eps);
    let tx = node_x(b, 0) - node_x(a, 0);
    let ty = node_x(b, 1) - node_x(a, 1);
    let out = ty * (node_x(a, 0) - node_x(o, 0)) - tx * (node_x(a, 1) - node_x(o, 1));
    let first = if out >= 0.0 { a } else { b };
    let second = if out >= 0.0 { b } else { a };
    rad_xy.push(node_x(first, 0));
    rad_xy.push(node_x(first, 1));
    rad_xy.push(node_x(second, 0));
    rad_xy.push(node_x(second, 1));
}

fn facet_length(i: i64) -> f64 {
    let dx = rad_xy[4 * i + 2] - rad_xy[4 * i];
    let dy = rad_xy[4 * i + 3] - rad_xy[4 * i + 1];
    return sqrt(dx * dx + dy * dy);
}

fn dist(a: f64[2], b: f64[2]) -> f64 { return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])); }

fn front(i: i64, p: f64[2]) -> f64 {
    let tx = rad_xy[4 * i + 2] - rad_xy[4 * i];
    let ty = rad_xy[4 * i + 3] - rad_xy[4 * i + 1];
    return ty * (p[0] - rad_xy[4 * i]) - tx * (p[1] - rad_xy[4 * i + 1]);
}

fn clip(i: i64, p: f64[2], q: f64[2]) -> f64[2], f64[2], bool {
    let tol = CLIP_TOLERANCE * facet_length(i) * facet_length(i);
    let fp = front(i, p);
    let fq = front(i, q);
    if fp <= tol && fq <= tol { return p, q, false; }
    if fp > tol && fq > tol { return p, q, true; }
    let s = fp / (fp - fq);
    let m: f64[2] = [p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])];
    return if fp > tol { p } else { m }, if fp > tol { m } else { q }, true;
}

fn hidden(p: f64[2], q: f64[2], i: i64, j: i64) -> bool {
    for k in rad_count() {
        if k == i || k == j { continue; }
        let a: f64[2] = [rad_xy[4 * k], rad_xy[4 * k + 1]];
        let b: f64[2] = [rad_xy[4 * k + 2], rad_xy[4 * k + 3]];
        let d1 = cross2(p, q, a);
        let d2 = cross2(p, q, b);
        let d3 = cross2(a, b, p);
        let d4 = cross2(a, b, q);
        if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)) { return true; }
    }
    return false;
}

fn cross2(o: f64[2], a: f64[2], b: f64[2]) -> f64 { return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); }

fn strings(p0: f64[2], p1: f64[2], q0: f64[2], q1: f64[2]) -> f64 {
    return 0.5 * abs(dist(p0, q0) + dist(p1, q1) - dist(p0, q1) - dist(p1, q0));
}

fn lerp(a: f64[2], b: f64[2], s: f64) -> f64[2] { return [a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])]; }

fn exchange(i: i64, j: i64) -> f64 {
    let pi0: f64[2] = [rad_xy[4 * i], rad_xy[4 * i + 1]];
    let pi1: f64[2] = [rad_xy[4 * i + 2], rad_xy[4 * i + 3]];
    let pj0: f64[2] = [rad_xy[4 * j], rad_xy[4 * j + 1]];
    let pj1: f64[2] = [rad_xy[4 * j + 2], rad_xy[4 * j + 3]];
    let a0, a1, seen_j = clip(i, pj0, pj1);
    guard seen_j else { return 0.0; }
    let b0, b1, seen_i = clip(j, pi0, pi1);
    guard seen_i else { return 0.0; }
    let whole = !hidden(b0, a0, i, j) && !hidden(b0, a1, i, j) && !hidden(b1, a0, i, j) && !hidden(b1, a1, i, j) &&
        !hidden(lerp(b0, b1, 0.5), lerp(a0, a1, 0.5), i, j);
    if whole { return strings(b0, b1, a0, a1); }
    var g = 0.0;
    let h = 1.0 / SUB;
    for s in SUB {
        let p0 = lerp(b0, b1, s * h);
        let p1 = lerp(b0, b1, (s + 1) * h);
        for t in SUB {
            let q0 = lerp(a0, a1, t * h);
            let q1 = lerp(a0, a1, (t + 1) * h);
            if !hidden(lerp(p0, p1, 0.5), lerp(q0, q1, 0.5), i, j) { g += strings(p0, p1, q0, q1); }
        }
    }
    return g;
}

fn view_factors(g: f64[:]) {
    let n = rad_count();
    guard g.len == n * n else { abort("view_factors: n × n entries"); }
    for i in n {
        for j in n { g[n * i + j] = if i == j { 0.0 } else { exchange(i, j) }; }
    }
    var worst = 0.0;
    for i in n {
        var s = 0.0;
        for j in i + 1..n {
            let m = 0.5 * (g[n * i + j] + g[n * j + i]);
            g[n * i + j] = m;
            g[n * j + i] = m;
        }
        for j in n { s += g[n * i + j]; }
        worst = max(worst, abs(s / facet_length(i) - 1.0));
    }
    if worst > RAD_CLOSURE {
        write_stderr(str("skein run: the faces of radiation(eps) do not close an enclosure: a view factor sum is off by ", worst, "\n"));
        exit(1);
    }
    close_rows(g);
}

fn close_rows(g: f64[:]) {
    let n = rad_count();
    var d: f64[>..] = [];
    for i in n { d.push(1.0); }
    for sweep in SINKHORN_SWEEPS {
        var worst = 0.0;
        for i in n {
            var s = 0.0;
            for j in n { s += g[n * i + j] * d[j]; }
            let target = facet_length(i) / (d[i] * s);
            worst = max(worst, abs(target - 1.0));
            d[i] *= sqrt(target);
        }
        if worst < SINKHORN_TOLERANCE { break; }
    }
    for i in n {
        for j in n { g[n * i + j] *= d[i] * d[j]; }
    }
}

fn exchange_matrix() {
    let n = rad_count();
    solvers::memory_check("radiation dense", n * n, n, 0, 0, ", or a coarser mesh on the faces of radiation(eps)");
    var f: f64[>..] = [];
    solvers::zeros(f, n * n);
    view_factors(f);
    for i in n {
        let l = facet_length(i);
        for j in n { f[n * i + j] /= l; }
    }
    var m: f64[>..] = [];
    solvers::zeros(m, n * n);
    for i in n {
        for j in n { m[n * i + j] = (if i == j { 1.0 } else { 0.0 }) - (1.0 - rad_eps[i]) * f[n * i + j]; }
    }
    var x: f64[>..] = [];
    solvers::zeros(x, n * n);
    for i in n { x[n * i + i] = rad_eps[i]; }
    guard dense_solve(m, x, n) else { abort("exchange_matrix: I - (1 - eps) F is singular"); }
    rad_q.clear();
    for i in n {
        for j in n {
            var fx = 0.0;
            for k in n { fx += f[n * i + k] * x[n * k + j]; }
            rad_q.push(x[n * i + j] - fx);
        }
    }
}

fn dense_solve(m: f64[:], b: f64[:], n: i64) -> bool {
    guard m.len == n * n && b.len == n * n else { abort("dense_solve: n × n matrices"); }
    for k in n {
        var p = k;
        for r in k + 1..n { if abs(m[n * r + k]) > abs(m[n * p + k]) { p = r; } }
        guard abs(m[n * p + k]) > PIVOT_FLOOR else { return false; }
        if p != k {
            for c in n {
                let t = m[n * k + c];
                m[n * k + c] = m[n * p + c];
                m[n * p + c] = t;
                let u = b[n * k + c];
                b[n * k + c] = b[n * p + c];
                b[n * p + c] = u;
            }
        }
        for r in k + 1..n {
            let s = m[n * r + k] / m[n * k + k];
            if s == 0.0 { continue; }
            for c in k..n { m[n * r + c] -= s * m[n * k + c]; }
            for c in n { b[n * r + c] -= s * b[n * k + c]; }
        }
    }
    for k in n {
        let r = n - 1 - k;
        for c in n {
            var s = b[n * r + c];
            for j in r + 1..n { s -= m[n * r + j] * b[n * j + c]; }
            b[n * r + c] = s / m[n * r + r];
        }
    }
    return true;
}

fn rad_weight(l: i64) -> f64 { return if order == 1 { 0.5 } else if l < 2 { 1.0 / 6.0 } else { 2.0 / 3.0 }; }

fn rad_nodes() -> i64 { return if order == 1 { 2 } else { 3 }; }

fn radiated(r: f64[:], dg: f64[:]) {
    let n = rad_count();
    var e: f64[>..] = [];
    var de: f64[>..] = [];
    for j in n {
        var s = 0.0;
        for l in rad_nodes() {
            let t = sol[rad_facet[3 * j + l]];
            s += rad_weight(l) * sigma_sb * t * t * t * t;
        }
        e.push(s);
    }
    for d in r.len {
        r[d] = 0.0;
        dg[d] = 0.0;
    }
    for i in n {
        var q = 0.0;
        for j in n { q += rad_q[n * i + j] * e[j]; }
        let lt = facet_length(i) * thickness;
        for l in rad_nodes() {
            let a = rad_facet[3 * i + l];
            r[a] += lt * rad_weight(l) * q;
            dg[a] += lt * rad_weight(l) * self_term(i, a);
        }
    }
}

fn self_term(i: i64, a: i64) -> f64 {
    let n = rad_count();
    var s = 0.0;
    let t = sol[a];
    for j in n {
        for l in rad_nodes() {
            if rad_facet[3 * j + l] == a { s += rad_q[n * i + j] * rad_weight(l) * 4.0 * sigma_sb * t * t * t; }
        }
    }
    return s;
}
