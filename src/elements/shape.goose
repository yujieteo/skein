import std;
import math;
namespace elements;

var hex_cells = false;

var cube_gauss = 2;

let HEX_CORNERS: i64[24] = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1];

let HEX_EDGES: i64[24] = [0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7];

let LINE_EDGES: i64[2] = [0, 1];
let TRI_EDGES: i64[6] = [0, 1, 1, 2, 2, 0];
let TET_EDGES: i64[12] = [0, 1, 1, 2, 0, 2, 0, 3, 1, 3, 2, 3];

fn edge_count(d: i64) -> i64 {
    guard d >= 1 && d <= 3 else { abort("edge_count: 1 <= d <= 3"); }
    if d == 1 { return 1; }
    if d == 2 { return 3; }
    return 6;
}

fn is_cube(d: i64) -> bool { return hex_cells && d >= 2; }

fn edge_node(d: i64, e: i64, side: i64) -> i64 {
    guard e >= 0 && e < edge_count(d) && (side == 0 || side == 1) else { abort("edge_node: e < edge_count(d) and side 0 or 1"); }
    if d == 1 { return LINE_EDGES[2 * e + side]; }
    if d == 2 { return TRI_EDGES[2 * e + side]; }
    return TET_EDGES[2 * e + side];
}

fn node_count_of(d: i64, order: i64) -> i64 {
    guard d >= 1 && d <= 3 && (order == 1 || order == 2) else { abort("node_count_of: 1 <= d <= 3 and order 1 or 2"); }
    if is_cube(d) { return if d == 3 { if order == 1 { 8 } else { 20 } } else { if order == 1 { 4 } else { 8 } }; }
    if order == 1 { return d + 1; }
    return d + 1 + edge_count(d);
}

fn bary_grad(k: i64, r: i64) -> f64 {
    guard k >= 0 && k <= 3 && r >= 0 && r < 3 else { abort("bary_grad: 0 <= k <= 3 and 0 <= r < 3"); }
    if k == 0 { return -1.0; }
    if k - 1 == r { return 1.0; }
    return 0.0;
}

fn cube_shape(d: i64, order: i64, xi: f64[3], n: f64[20]&, dn: f64[60]&) {
    if order == 2 {
        cube_shape_2(d, xi, n, dn);
        return;
    }
    for l in node_count_of(d, 1) {
        var f: f64[3] = [1.0; 3];
        var df: f64[3] = [0.0; 3];
        for r in d {
            let up = HEX_CORNERS[3 * l + r] == 1;
            f[r] = if up { xi[r] } else { 1.0 - xi[r] };
            df[r] = if up { 1.0 } else { -1.0 };
        }
        n[l] = f[0] * f[1] * f[2];
        for r in d {
            var g = df[r];
            for t in d { if t != r { g *= f[t]; } }
            dn[3 * l + r] = g;
        }
    }
}

fn natural(d: i64, l: i64, r: i64) -> f64 {
    if l < (1 << d) { return if HEX_CORNERS[3 * l + r] == 1 { 1.0 } else { -1.0 }; }
    let e = l - (1 << d);
    let a = HEX_CORNERS[3 * HEX_EDGES[2 * e] + r];
    let b = HEX_CORNERS[3 * HEX_EDGES[2 * e + 1] + r];
    if a != b { return 0.0; }
    return if a == 1 { 1.0 } else { -1.0 };
}

fn cube_shape_2(d: i64, xi: f64[3], n: f64[20]&, dn: f64[60]&) {
    var s: f64[3] = [0.0; 3];
    for r in d { s[r] = 2.0 * xi[r] - 1.0; }
    let scale = if d == 3 { 8.0 } else { 4.0 };
    for l in node_count_of(d, 2) {
        var a: f64[3] = [0.0; 3];
        var zeros = 0;
        for r in d {
            a[r] = natural(d, l, r);
            if a[r] == 0.0 { zeros++; }
        }
        var f: f64[3] = [1.0; 3];
        var df: f64[3] = [0.0; 3];
        var lin = -(d as f64 - 1.0);
        for r in d {
            if a[r] == 0.0 {
                f[r] = 1.0 - s[r] * s[r];
                df[r] = -2.0 * s[r];
            } else {
                f[r] = 1.0 + a[r] * s[r];
                df[r] = a[r];
                lin += a[r] * s[r];
            }
        }
        let corner = zeros == 0;
        var value = 1.0;
        for r in d { value *= f[r]; }
        n[l] = value * (if corner { lin } else { 1.0 }) / (if corner { scale } else { 0.5 * scale });
        for r in d {
            var g = df[r];
            for t in d { if t != r { g *= f[t]; } }
            var full = g * (if corner { lin } else { 1.0 });
            if corner { full += value * a[r]; }
            dn[3 * l + r] = 2.0 * full / (if corner { scale } else { 0.5 * scale });
        }
    }
}

fn shape(d: i64, order: i64, xi: f64[3], n: f64[20]&, dn: f64[60]&) {
    guard node_count_of(d, order) <= n.len && 3 * n.len <= dn.len else { abort("shape: n and dn hold the nodes of the simplex"); }
    if is_cube(d) {
        cube_shape(d, order, xi, n, dn);
        return;
    }
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

let TET_A = 0.1381966011250105;
let TET_B = 0.5854101966249685;
let LINE_G = 0.21132486540518713;
let GAUSS_3 = 0.7745966692414834;

fn quad_count(d: i64) -> i64 {
    guard d >= 1 && d <= 3 else { abort("quad_count: 1 <= d <= 3"); }
    if is_cube(d) { return if cube_gauss == 3 { if d == 3 { 27 } else if d == 2 { 9 } else { 3 } } else { 1 << d }; }
    return d + 1;
}

fn quad_point(d: i64, q: i64) -> f64[3], f64 {
    guard q >= 0 && q < quad_count(d) else { abort("quad_point: 0 <= q < quad_count(d)"); }
    if is_cube(d) {
        var p: f64[3] = [0.0; 3];
        var w = 1.0;
        var rest = q;
        for r in d {
            if cube_gauss == 3 {
                let k = rest % 3;
                rest /= 3;
                p[r] = if k == 0 { 0.5 * (1.0 - GAUSS_3) } else if k == 1 { 0.5 } else { 0.5 * (1.0 + GAUSS_3) };
                w *= if k == 1 { 8.0 / 18.0 } else { 5.0 / 18.0 };
            } else {
                p[r] = if rest % 2 == 0 { LINE_G } else { 1.0 - LINE_G };
                rest /= 2;
                w *= 0.5;
            }
        }
        return p, w;
    }
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

let TET_FACES: i64[24] = [
    0, 1, 2, 4, 5, 6,
    0, 1, 3, 4, 8, 7,
    1, 2, 3, 5, 9, 8,
    0, 2, 3, 6, 9, 7,
];

let HEX_FACES: i64[24] = [
    0, 3, 7, 4,
    1, 2, 6, 5,
    0, 1, 5, 4,
    3, 2, 6, 7,
    0, 1, 2, 3,
    4, 5, 6, 7,
];

fn facet_count(d: i64) -> i64 {
    guard d == 2 || d == 3 else { abort("facet_count: d is 2 or 3"); }
    if is_cube(d) && d == 3 { return 6; }
    return d + 1;
}

fn facet_vertices(d: i64) -> i64 {
    guard d == 2 || d == 3 else { abort("facet_vertices: d is 2 or 3"); }
    if is_cube(d) && d == 3 { return 4; }
    return d;
}

fn facet_node(d: i64, f: i64, l: i64) -> i64 {
    guard f >= 0 && f < facet_count(d) && l >= 0 && l < node_count_of(d - 1, 2) else { abort("facet_node: f < facet_count(d) and l < the nodes of a facet"); }
    if is_cube(d) && d == 3 { return if l < 4 { HEX_FACES[4 * f + l] } else { hex_mid_node(HEX_FACES[4 * f + l - 4], HEX_FACES[4 * f + (l - 3) % 4]) }; }
    if d == 3 { return TET_FACES[6 * f + l]; }
    if l < 2 { return TRI_EDGES[2 * f + l]; }
    return 3 + f;
}

fn hex_mid_node(a: i64, b: i64) -> i64 {
    for e in 12 {
        if (HEX_EDGES[2 * e] == a && HEX_EDGES[2 * e + 1] == b) || (HEX_EDGES[2 * e] == b && HEX_EDGES[2 * e + 1] == a) { return 8 + e; }
    }
    abort("hex_mid_node: a and b are the ends of an edge");
}
