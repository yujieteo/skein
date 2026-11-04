import std;
import math;
namespace mesher;
import .geometry;
import .refine;
import .quality;

var sharp_angle = 30.0;
var sharp_cos = 0.8660254037844387;
var pair_cos = 0.9135454576426009;

fn set_feature_angle(degrees: f64) {
    guard degrees > 0.0 && degrees < 90.0 else { abort("set_feature_angle: an angle from 0 to 90 degrees"); }
    sharp_angle = degrees;
    sharp_cos = cos(radians(degrees));
    pair_cos = cos(radians(0.8 * degrees));
}
let FD = 1.0e-8;
let BISECT_STOP = 30.0;
let CURVE_DIV = 2.0;
let MAX_SAMPLES = 100000;

struct Ridge {
    found: bool,
    x: f64[3],
    na: f64[3],
    nb: f64[3],
}

var curve_x: f64[>..] = [];
var curve_n: f64[>..] = [];
var curve_start: i64[>..] = [0];
var curve_closed: bool[>..] = [];

fn clear_curves() {
    curve_x.clear();
    curve_n.clear();
    curve_start.clear();
    curve_start.push(0);
    curve_closed.clear();
}

fn curve_count() -> i64 { return curve_closed.len; }

fn axpy3(a: f64[3], s: f64, b: f64[3]) -> f64[3] {
    return [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]];
}

fn unit3(a: f64[3]) -> f64[3] {
    let n = norm3(a);
    return [a[0] / n, a[1] / n, a[2] / n];
}

fn sample_at(k: i64) -> f64[3] { return [curve_x[3 * k], curve_x[3 * k + 1], curve_x[3 * k + 2]]; }

fn sample_normal(k: i64) -> f64[3] { return [curve_n[3 * k], curve_n[3 * k + 1], curve_n[3 * k + 2]]; }

fn shape_scale(g: i64) -> f64 {
    let b = bounds(g);
    return sqrt((b[3] - b[0]) * (b[3] - b[0]) + (b[4] - b[1]) * (b[4] - b[1]) + (b[5] - b[2]) * (b[5] - b[2]));
}

fn gradient3(g: i64, x: f64[3], eps: f64) -> f64[3] {
    var out: f64[3] = [0.0, 0.0, 0.0];
    for c in 3 {
        var lo = x;
        var hi = x;
        lo[c] -= eps;
        hi[c] += eps;
        out[c] = (phi(g, hi[0], hi[1], hi[2]) - phi(g, lo[0], lo[1], lo[2])) / (2.0 * eps);
    }
    return out;
}

fn settle(g: i64, x0: f64[3], d: f64) -> f64[3], f64[3], bool {
    let eps = FD * d;
    var x = x0;
    var n: f64[3] = [0.0, 0.0, 0.0];
    for step in 8 {
        let f = phi(g, x[0], x[1], x[2]);
        n = gradient3(g, x, eps);
        let n2 = dot3(n, n);
        guard n2 > 0.0 else { return x, n, false; }
        if abs(f) <= 1.0e-12 * d { return x, unit3(n), true; }
        x = axpy3(x, -f / n2, n);
    }
    return x, n, false;
}

fn ridge_between(g: i64, u: f64[3], nu: f64[3], v: f64[3], nv: f64[3], d: f64) -> Ridge {
    var r = Ridge { found: false, x: u, na: nu, nb: nv };
    var a = u;
    var na = nu;
    var b = v;
    var nb = nv;
    guard dot3(na, nb) < sharp_cos else { return r; }
    let stop = BISECT_STOP * FD * d;
    for k in 80 {
        guard norm3(sub3(a, b)) > stop else { break; }
        let m, nm, ok = settle(g, axpy3(a, 0.5, sub3(b, a)), d);
        guard ok else { return r; }
        if dot3(na, nm) < dot3(nm, nb) {
            b = m;
            nb = nm;
        } else {
            a = m;
            na = nm;
        }
    }
    guard dot3(na, nb) < sharp_cos else { return r; }
    let c = dot3(na, nb);
    let mid = axpy3(a, 0.5, sub3(b, a));
    let r1 = -dot3(na, sub3(mid, a));
    let r2 = -dot3(nb, sub3(mid, b));
    let det = 1.0 - c * c;
    let p = (r1 - c * r2) / det;
    let q = (r2 - c * r1) / det;
    r.x = axpy3(axpy3(mid, p, na), q, nb);
    r.na = na;
    r.nb = nb;
    r.found = true;
    return r;
}

fn angle_of(a: f64[3], b: f64[3]) -> f64 { return acos(clamp(dot3(a, b), -1.0, 1.0)); }

fn find_ridge(g: i64, u: f64[3], nu: f64[3], v: f64[3], nv: f64[3], d: f64) -> Ridge {
    var r = Ridge { found: false, x: u, na: nu, nb: nv };
    let m, nm, ok = settle(g, axpy3(u, 0.5, sub3(v, u)), d);
    guard ok else { return r; }
    let a1 = angle_of(nu, nm);
    let a2 = angle_of(nm, nv);
    let big = max(a1, a2);
    let small = min(a1, a2);
    guard big > radians(sharp_angle) else { return r; }
    guard big > 2.5 * small || angle_of(nu, nv) > radians(45.0) else { return r; }
    return ridge_between(g, u, nu, v, nv, d);
}

fn walk(g: i64, r0: Ridge, dir: f64, size: f64, d: f64, out: f64[>..]&, normals: f64[>..]&) -> bool {
    var x = r0.x;
    var t = unit3(cross3(r0.na, r0.nb));
    for c in 3 { t[c] *= dir; }
    var nav = unit3([r0.na[0] + r0.nb[0], r0.na[1] + r0.nb[1], r0.na[2] + r0.nb[2]]);
    var travelled = 0.0;
    for count in MAX_SAMPLES {
        let s0 = wanted_size(x[0], x[1], x[2], 0.0, size) / CURVE_DIV;
        var s = s0;
        var next = Ridge { found: false, x: x, na: r0.na, nb: r0.nb };
        for attempt in 8 {
            let y = axpy3(x, s, t);
            let w = unit3(cross3(t, nav));
            let u, nu, ok_u = settle(g, axpy3(y, -0.4 * s, w), d);
            let v, nv, ok_v = settle(g, axpy3(y, 0.4 * s, w), d);
            if ok_u && ok_v {
                let r = ridge_between(g, u, nu, v, nv, d);
                if r.found && norm3(sub3(r.x, y)) <= 0.5 * s {
                    next = r;
                    break;
                }
            }
            s *= 0.5;
        }
        guard next.found else { return false; }
        var t2 = unit3(cross3(next.na, next.nb));
        if dot3(t2, t) < 0.0 { t2 = [-t2[0], -t2[1], -t2[2]]; }
        guard dot3(t2, t) > 0.7 else { return false; }
        travelled += norm3(sub3(next.x, x));
        x = next.x;
        t = t2;
        nav = unit3([next.na[0] + next.nb[0], next.na[1] + next.nb[1], next.na[2] + next.nb[2]]);
        if travelled > 3.0 * s0 && norm3(sub3(x, r0.x)) < 1.2 * s0 { return true; }
        out.append(x);
        normals.append(nav);
    }
    return false;
}

fn trace_curve(g: i64, r0: Ridge, size: f64, d: f64) {
    var fwd: f64[>..] = [];
    var fwd_n: f64[>..] = [];
    let closed = walk(g, r0, 1.0, size, d, fwd, fwd_n);
    var back: f64[>..] = [];
    var back_n: f64[>..] = [];
    if !closed { walk(g, r0, -1.0, size, d, back, back_n); }
    guard back.len / 3 + 1 + fwd.len / 3 >= 3 else { return; }
    for k in back.len / 3 {
        let j = back.len / 3 - 1 - k;
        for c in 3 {
            curve_x.push(back[3 * j + c]);
            curve_n.push(back_n[3 * j + c]);
        }
    }
    curve_x.append(r0.x);
    curve_n.append(unit3([r0.na[0] + r0.nb[0], r0.na[1] + r0.nb[1], r0.na[2] + r0.nb[2]]));
    curve_x.append(fwd);
    curve_n.append(fwd_n);
    curve_start.push(curve_x.len / 3);
    curve_closed.push(closed);
}

fn near_curve(x: f64[3], tol: f64) -> bool {
    for k in curve_x.len / 3 {
        if norm3(sub3(sample_at(k), x)) <= tol { return true; }
    }
    return false;
}
