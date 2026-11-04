import std;
import math;
namespace shells_gen;
import syntax from "../syntax";

let S_BOX = 0;
let S_CYLINDER = 1;
let S_MOVE = 2;
let S_TURN = 3;
let S_UNION = 4;
let S_INTER = 5;
let S_SUB = 6;
let S_PATCH = 7;

struct Shape {
    kind: i64,
    a: i64,
    b: i64,
    p: f64[8],
}

var shapes: Shape[>..] = [];
var open_patch = false;
var patch_shape = -1;

fn add_shape(k: i64, a: i64, b: i64, p: f64[8]) -> i64 {
    shapes.push(Shape { kind: k, a: a, b: b, p: p });
    return shapes.len - 1;
}

fn dot3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

fn sub3(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }

recursive fn phi(g: i64, q: f64[3]) -> f64 {
    let s = shapes[g];
    if s.kind == S_BOX {
        let r = s.p[6];
        var outside = 0.0;
        var inner = -1.0e300;
        for k in 3 {
            let d = abs(q[k] - 0.5 * (s.p[k] + s.p[3 + k])) - (0.5 * (s.p[3 + k] - s.p[k]) - r);
            outside += sq(max(d, 0.0));
            inner = max(inner, d);
        }
        return sqrt(outside) + min(inner, 0.0) - r;
    }
    if s.kind == S_CYLINDER {
        let d = sub3(q, [s.p[0], s.p[1], s.p[2]]);
        let a: f64[3] = [s.p[3], s.p[4], s.p[5]];
        let along = dot3(d, a);
        let dr = sqrt(max(dot3(d, d) - along * along, 0.0)) - s.p[6];
        let dh = abs(along - 0.5 * s.p[7]) - 0.5 * s.p[7];
        return sqrt(sq(max(dr, 0.0)) + sq(max(dh, 0.0))) + min(max(dr, dh), 0.0);
    }
    if s.kind == S_MOVE { return phi(s.a, sub3(q, [s.p[0], s.p[1], s.p[2]])); }
    if s.kind == S_TURN { return phi(s.a, rodrigues(q, [s.p[0], s.p[1], s.p[2]], -s.p[3])); }
    if s.kind == S_UNION { return min(phi(s.a, q), phi(s.b, q)); }
    if s.kind == S_INTER { return max(phi(s.a, q), phi(s.b, q)); }
    if s.kind == S_SUB { return max(phi(s.a, q), -phi(s.b, q)); }
    if open_patch { return -1.0e30; }
    let t = theta_of(q);
    return max(max(x0 - q[0], q[0] - x1), radius * max(t0 - t, t - t1));
}

fn cut_phi(g: i64, q: f64[3]) -> f64 {
    open_patch = true;
    let f = phi(g, q);
    open_patch = false;
    return f;
}

var fn_names: u8[..64][>..] = [];
var fn_nodes: i64[>..] = [];
var bind_names: u8[..64][>..] = [];
var bind_vals: Val[>..] = [];

fn user_fn(v: i64) -> i64 {
    if kind(v) != syntax::N_CALL || kind(child(v, 0)) != syntax::N_NAME { return -1; }
    let i = find(fn_names, name(child(v, 0)));
    return if i >= 0 { fn_nodes[i] } else { -1 };
}

fn bind_call(v: i64, f: i64) -> i64 {
    var count = 0;
    var body = -1;
    var c = first(f);
    while c >= 0 {
        if kind(c) == syntax::N_FNPARAM { count++; }
        body = c;
        c = next(c);
    }
    if count != syntax::child_count(v) - 1 || count > 8 || depth > 100 {
        unsupported(tok(v), "the call does not give one value for each parameter", "give each parameter of the function");
        return -1;
    }
    var vals = default<Val[8]>();
    for k in count { vals[k] = value(arg(v, k)); }
    var k = 0;
    c = first(f);
    while c >= 0 {
        if kind(c) == syntax::N_FNPARAM {
            var w: u8[..64] = [];
            set64(w, name(c));
            bind_names.push(w);
            bind_vals.push(vals[k]);
            k++;
        }
        c = next(c);
    }
    depth++;
    return body;
}

fn unbind(f: i64) {
    var c = first(f);
    while c >= 0 {
        if kind(c) == syntax::N_FNPARAM {
            bind_names.pop();
            bind_vals.pop();
        }
        c = next(c);
    }
    depth--;
}

fn bound(word: u8[:]) -> i64 {
    for k in bind_names.len {
        let i = bind_names.len - 1 - k;
        if bind_names[i] == word { return i; }
    }
    return -1;
}

let SOLID_FIX = "write the solid with box, cylinder, move, rotate, fillet of a box, |, & and -, or a fn that gives one";

recursive fn shape_of(v: i64) -> i64 {
    if kind(v) == syntax::N_NAME {
        let w = name(v);
        if shell_tok >= 0 && w == shell_name { return patch_shape; }
        let su = find(surface_names, w);
        if su >= 0 { return surface_shapes[su]; }
        let si = find(solid_names, w);
        if si >= 0 && depth < 100 {
            if solid_shapes[si] < 0 {
                depth++;
                solid_shapes[si] = shape_of(solid_nodes[si]);
                depth--;
            }
            return solid_shapes[si];
        }
    }
    let f = user_fn(v);
    if f >= 0 {
        let body = bind_call(v, f);
        if body < 0 { return -1; }
        let g = shape_of(body);
        unbind(f);
        return g;
    }
    if kind(v) == syntax::N_BINARY && (op(v) == "|" || op(v) == "&" || op(v) == "-") {
        let a = shape_of(child(v, 0));
        let b = shape_of(child(v, 1));
        if a < 0 || b < 0 { return -1; }
        let k = if op(v) == "|" { S_UNION } else if op(v) == "&" { S_INTER } else { S_SUB };
        return add_shape(k, a, b, [0.0; 8]);
    }
    return primitive(v);
}

recursive fn primitive(v: i64) -> i64 {
    let count = if kind(v) == syntax::N_CALL { syntax::child_count(v) - 1 } else { -1 };
    if is_call(v, "box") && count == 2 {
        let a = vec3(value(arg(v, 0)));
        let b = vec3(value(arg(v, 1)));
        return add_shape(S_BOX, -1, -1, [min(a[0], b[0]), min(a[1], b[1]), min(a[2], b[2]),
                                         max(a[0], b[0]), max(a[1], b[1]), max(a[2], b[2]), 0.0, 0.0]);
    }
    if is_call(v, "cylinder") && count == 4 {
        let c = vec3(value(arg(v, 0)));
        let d = vec3(value(arg(v, 1)));
        let l = sqrt(dot3(d, d));
        return add_shape(S_CYLINDER, -1, -1, [c[0], c[1], c[2], d[0] / l, d[1] / l, d[2] / l,
                                              scalar_value(arg(v, 2)), scalar_value(arg(v, 3))]);
    }
    if is_call(v, "move") && count == 2 {
        let a = shape_of(arg(v, 0));
        let d = vec3(value(arg(v, 1)));
        return if a < 0 { -1 } else { add_shape(S_MOVE, a, -1, [d[0], d[1], d[2], 0.0, 0.0, 0.0, 0.0, 0.0]) };
    }
    if is_call(v, "rotate") && count == 3 {
        let a = shape_of(arg(v, 0));
        let d = vec3(value(arg(v, 1)));
        let l = sqrt(dot3(d, d));
        let angle = scalar_value(arg(v, 2));
        return if a < 0 { -1 } else { add_shape(S_TURN, a, -1, [d[0] / l, d[1] / l, d[2] / l, angle, 0.0, 0.0, 0.0, 0.0]) };
    }
    if (is_call(v, "pattern") || is_call(v, "pattern_polar")) && count == 3 {
        let a = shape_of(arg(v, 0));
        let n = round(scalar_value(arg(v, 2))) as! i64;
        if a < 0 || n < 1 { return -1; }
        let d = vec3(value(arg(v, 1)));
        let l = sqrt(dot3(d, d));
        var g = -1;
        for k in n {
            let c = if is_call(v, "pattern") { add_shape(S_MOVE, a, -1, [k * d[0], k * d[1], k * d[2], 0.0, 0.0, 0.0, 0.0, 0.0]) }
                    else { add_shape(S_TURN, a, -1, [d[0] / l, d[1] / l, d[2] / l, 2.0 * PI * k / n, 0.0, 0.0, 0.0, 0.0]) };
            g = if g < 0 { c } else { add_shape(S_UNION, g, c, [0.0; 8]) };
        }
        return g;
    }
    if is_call(v, "fillet") && count == 2 {
        let a = shape_of(arg(v, 0));
        return if a < 0 { -1 } else { rounded(a, scalar_value(arg(v, 1)), tok(v)) };
    }
    unsupported(tok(v), str("the solid '", source_of(v), "' is not in this slice"), SOLID_FIX);
    return -1;
}

recursive fn rounded(g: i64, r: f64, at: i64) -> i64 {
    let s = shapes[g];
    if s.kind == S_MOVE || s.kind == S_TURN {
        let a = rounded(s.a, r, at);
        return if a < 0 { -1 } else { add_shape(s.kind, a, -1, s.p) };
    }
    var half = 1.0e300;
    for k in 3 { half = min(half, 0.5 * (s.p[3 + k] - s.p[k])); }
    if s.kind != S_BOX || r < 0.0 || r > half {
        unsupported(at, "this slice has the fillet of a box only, with a radius up to half its smallest side", SOLID_FIX);
        return -1;
    }
    var p = s.p;
    p[6] = r;
    return add_shape(S_BOX, -1, -1, p);
}

fn is_condition(v: i64) -> bool {
    if kind(v) != syntax::N_BINARY { return false; }
    let o = op(v);
    return o == "==" || o == ">=" || o == "<=" || o == ">" || o == "<" || o == "&&";
}

fn holds(cond: i64, p: f64[3], tol: f64) -> bool {
    var parts: i64[>..] = [];
    if !terms(cond, parts) { return false; }
    for t in parts {
        let m, size = margin(t, p);
        if m < -tol { return false; }
    }
    return true;
}

recursive fn shell_curve(v: i64) -> bool {
    if kind(v) == syntax::N_NAME {
        let r = find_region(v);
        return r >= 0 && regions[r].node >= 0;
    }
    if is_call(v, "boundary") && syntax::child_count(v) == 3 && kind(arg(v, 0)) == syntax::N_NAME {
        let w = name(arg(v, 0));
        return (shell_tok >= 0 && w == shell_name) || find(surface_names, w) >= 0;
    }
    if kind(v) == syntax::N_BINARY && (op(v) == "&" || op(v) == "-" || op(v) == "|") { return shell_curve(child(v, 0)); }
    return false;
}

fn node_shape(v: i64) -> i64 {
    for n, i in shape_nodes { if n == v { return shape_ids[i]; } }
    let g = shape_of(v);
    shape_nodes.push(v);
    shape_ids.push(g);
    return g;
}

recursive fn on_curve(v: i64, p: f64[3], tol: f64) -> bool {
    if kind(v) == syntax::N_NAME { return on_curve(regions[find_region(v)].node, p, tol); }
    if kind(v) == syntax::N_BINARY && op(v) == "|" { return on_curve(child(v, 0), p, tol) || on_curve(child(v, 1), p, tol); }
    if kind(v) == syntax::N_BINARY {
        let f = phi(node_shape(child(v, 1)), p);
        return on_curve(child(v, 0), p, tol) && (if op(v) == "&" { f <= tol } else { f >= -tol });
    }
    let s = node_shape(arg(v, 0));
    let c = arg(v, 1);
    if is_condition(c) { return holds(c, p, tol) && abs(phi(s, p)) <= tol; }
    return abs(phi(node_shape(c), p)) <= tol && phi(s, p) <= tol;
}

let GROW = 0.2;
let SNAP = 0.45;

var order = 1;
var xs: f64[>..] = [];
var ts: f64[>..] = [];
var x_fixed: u8[>..] = [];
var t_fixed: u8[>..] = [];
var px_lat: f64[>..] = [];
var pt_lat: f64[>..] = [];
var cell_mat: i64[>..] = [];
var snaps: i64[>..] = [];

var refine_nodes: i64[>..] = [];
var refine_sizes: f64[>..] = [];
var refine_x: f64[>..] = [];
var refine_t: f64[>..] = [];
var refine_of: i64[>..] = [];

fn nx_cells() -> i64 { return xs.len - 1; }
fn nt_cells() -> i64 { return ts.len - 1; }
fn lat_nx() -> i64 { return order * nx_cells(); }
fn lat_nt() -> i64 { return order * nt_cells(); }
fn lat(i: i64, j: i64) -> i64 { return i + (lat_nx() + 1) * j; }
fn lat_point(k: i64) -> f64[3] { return shell_point(px_lat[k], pt_lat[k]); }

fn geo_tol() -> f64 { return 1.0e-6 * max(x1 - x0, radius * (t1 - t0)); }

fn extent(g: i64) -> f64[4] {
    let n = 400;
    var lo: f64[2] = [1.0e300, 1.0e300];
    var hi: f64[2] = [-1.0e300, -1.0e300];
    var at: f64[4] = [0.0; 4];
    for i in n {
        for j in n {
            let x = x0 + (x1 - x0) * (i + 0.5) / n;
            let t = t0 + (t1 - t0) * (j + 0.5) / n;
            if phi(g, shell_point(x, t)) >= 0.0 { continue; }
            if x < lo[0] { lo[0] = x; at[0] = t; }
            if x > hi[0] { hi[0] = x; at[1] = t; }
            if t < lo[1] { lo[1] = t; at[2] = x; }
            if t > hi[1] { hi[1] = t; at[3] = x; }
        }
    }
    if lo[0] > hi[0] { return [x0, x1, t0, t1]; }
    return [edge_bisect(g, lo[0], x0, at[0], true), edge_bisect(g, hi[0], x1, at[1], true),
            edge_bisect(g, lo[1], t0, at[2], false), edge_bisect(g, hi[1], t1, at[3], false)];
}

fn edge_bisect(g: i64, c: f64, side: f64, other: f64, along_x: bool) -> f64 {
    let pt = if along_x { shell_point(side, other) } else { shell_point(other, side) };
    if phi(g, pt) < 0.0 { return side; }
    var a = c;
    var b = side;
    for k in 100 {
        let m = 0.5 * (a + b);
        let q = if along_x { shell_point(m, other) } else { shell_point(other, m) };
        if phi(g, q) < 0.0 { a = m; } else { b = m; }
    }
    return 0.5 * (a + b);
}

fn spacing(c: f64, along_x: bool) -> f64 {
    var h = mesh_size;
    for k in refine_of.len {
        let d = if along_x { abs(c - refine_x[k]) } else { radius * abs(c - refine_t[k]) };
        h = min(h, refine_sizes[refine_of[k]] + GROW * d);
    }
    return h;
}

fn fill_lines(a: f64, b: f64, along_x: bool, out: f64[>..]&, fixed: u8[>..]&) {
    let steps = 256;
    let scale = if along_x { 1.0 } else { radius };
    var cum: f64[>..] = [0.0];
    for k in steps {
        let c = a + (b - a) * (k + 0.5) / steps;
        cum.push(cum[k] + scale * (b - a) / steps / spacing(c, along_x));
    }
    let count = max(ceil(cum[steps] - 1.0e-9) as! i64, 1);
    var k = 0;
    for m in 1..count {
        let target = cum[steps] * m / count;
        while cum[k + 1] < target { k++; }
        let f = (target - cum[k]) / (cum[k + 1] - cum[k]);
        out.push(a + (b - a) * (k + f) / steps);
        fixed.push(0);
    }
    out.push(b);
    fixed.push(1);
}

fn required(values: f64[:], lo: f64, hi: f64, tol: f64) -> f64[>..] {
    var out: f64[>..] = [lo, hi];
    for v in values {
        if v < lo + tol || v > hi - tol { continue; }
        var seen = false;
        for u in out { if abs(u - v) <= tol { seen = true; } }
        if !seen { out.push(v); }
    }
    sort(out);
    return out;
}

fn sample_refines(e: f64[4]) {
    for r in refine_nodes.len {
        let s = 0.25 * refine_sizes[r];
        let n = min(max(ceil((e[1] - e[0]) / s) as! i64, 1), 2000);
        let m = min(max(ceil(radius * (e[3] - e[2]) / s) as! i64, 1), 2000);
        for i in n + 1 {
            for j in m + 1 {
                let x = e[0] + (e[1] - e[0]) * i / n;
                let t = e[2] + (e[3] - e[2]) * j / m;
                if on_curve(refine_nodes[r], shell_point(x, t), 0.5 * s) {
                    refine_x.push(x);
                    refine_t.push(t);
                    refine_of.push(r);
                }
            }
        }
    }
}

fn make_lines(g: i64, ring_x: f64[:], line_t: f64[:]) {
    let e = extent(g);
    sample_refines(e);
    let tol = geo_tol();
    let rx = required(ring_x, e[0], e[1], tol);
    let rt = required(line_t, e[2], e[3], tol / radius);
    xs.clear();
    ts.clear();
    x_fixed.clear();
    t_fixed.clear();
    xs.push(rx[0]);
    x_fixed.push(1);
    for k in rx.len - 1 { fill_lines(rx[k], rx[k + 1], true, xs, x_fixed); }
    ts.push(rt[0]);
    t_fixed.push(1);
    for k in rt.len - 1 { fill_lines(rt[k], rt[k + 1], false, ts, t_fixed); }
}

fn corner(i: i64, j: i64) -> i64 { return lat(order * i, order * j); }

fn local_size(i: i64, j: i64) -> f64 {
    var h = 1.0e300;
    if i > 0 { h = min(h, xs[i] - xs[i - 1]); }
    if i < nx_cells() { h = min(h, xs[i + 1] - xs[i]); }
    if j > 0 { h = min(h, radius * (ts[j] - ts[j - 1])); }
    if j < nt_cells() { h = min(h, radius * (ts[j + 1] - ts[j])); }
    return h;
}

fn root_on(g: i64, xa: f64, ta: f64, xb: f64, tb: f64) -> f64, f64, bool {
    var fa = cut_phi(g, shell_point(xa, ta));
    let fb = cut_phi(g, shell_point(xb, tb));
    if (fa < 0.0) == (fb < 0.0) { return 0.0, 0.0, false; }
    var a = 0.0;
    var b = 1.0;
    for k in 100 {
        let m = 0.5 * (a + b);
        let f = cut_phi(g, shell_point(xa + m * (xb - xa), ta + m * (tb - ta)));
        if (f < 0.0) == (fa < 0.0) { a = m; } else { b = m; }
    }
    let m = 0.5 * (a + b);
    return xa + m * (xb - xa), ta + m * (tb - ta), true;
}

fn root_along(g: i64, x: f64, t: f64, reach: f64, along_x: bool) -> f64, f64, bool {
    let n = 16;
    let d = if along_x { reach / n } else { reach / (n * radius) };
    for k in n {
        for s in [1.0, -1.0] {
            let a = s * d * k;
            let b = s * d * (k + 1);
            if along_x {
                let xr, tr, found = root_on(g, x + a, t, x + b, t);
                if found { return xr, tr, true; }
            } else {
                let xr, tr, found = root_on(g, x, t + a, x, t + b);
                if found { return xr, tr, true; }
            }
        }
    }
    return x, t, false;
}

fn root_near(g: i64, x: f64, t: f64, h: f64) -> f64, f64, bool {
    var xc = x;
    var tc = t;
    let e = 1.0e-6 * h;
    for k in 30 {
        let f = cut_phi(g, shell_point(xc, tc));
        if abs(f) <= 1.0e-12 * h { return xc, tc, true; }
        let gx = (cut_phi(g, shell_point(xc + e, tc)) - cut_phi(g, shell_point(xc - e, tc))) / (2.0 * e);
        let gs = (cut_phi(g, shell_point(xc, tc + e / radius)) - cut_phi(g, shell_point(xc, tc - e / radius))) / (2.0 * e);
        let g2 = gx * gx + gs * gs;
        if g2 <= 1.0e-24 { return x, t, false; }
        xc -= f * gx / g2;
        tc -= f * gs / g2 / radius;
    }
    return xc, tc, abs(cut_phi(g, shell_point(xc, tc))) <= 1.0e-9 * h;
}

fn snap_corner(i: i64, j: i64) -> bool {
    if x_fixed[i] != 0 && t_fixed[j] != 0 { return false; }
    let k = corner(i, j);
    let h = local_size(i, j);
    let p = lat_point(k);
    var best = -1;
    var size = SNAP * h;
    for g in snaps {
        let f = abs(cut_phi(g, p));
        if f < size {
            best = g;
            size = f;
        }
    }
    if best < 0 { return false; }
    let x = px_lat[k];
    let t = pt_lat[k];
    var xr = x;
    var tr = t;
    var found = false;
    if x_fixed[i] != 0 {
        xr, tr, found = root_along(best, x, t, 0.6 * h, false);
    } else if t_fixed[j] != 0 {
        xr, tr, found = root_along(best, x, t, 0.6 * h, true);
    } else {
        xr, tr, found = root_near(best, x, t, h);
    }
    if !found || sqrt(sq(xr - x) + sq(radius * (tr - t))) > 0.6 * h { return false; }
    px_lat[k] = xr;
    pt_lat[k] = tr;
    return true;
}

fn default_place(a: i64, b: i64) -> f64, f64 {
    let i = a / order;
    let j = b / order;
    let x = if a % order == 0 { xs[i] } else { 0.5 * (xs[i] + xs[i + 1]) };
    let t = if b % order == 0 { ts[j] } else { 0.5 * (ts[j] + ts[j + 1]) };
    return x, t;
}

fn cell_quality(i: i64, j: i64) -> f64 {
    let ks = [corner(i, j), corner(i + 1, j), corner(i + 1, j + 1), corner(i, j + 1)];
    var q = 1.0;
    for c in 4 {
        let a = ks[c];
        let b = ks[(c + 1) % 4];
        let d = ks[(c + 3) % 4];
        let ex = px_lat[b] - px_lat[a];
        let es = radius * (pt_lat[b] - pt_lat[a]);
        let fx = px_lat[d] - px_lat[a];
        let fs = radius * (pt_lat[d] - pt_lat[a]);
        let l = sqrt((ex * ex + es * es) * (fx * fx + fs * fs));
        q = min(q, if l > 0.0 { (ex * fs - es * fx) / l } else { -1.0 });
    }
    return q;
}

fn cell_centre(i: i64, j: i64) -> f64[3] {
    var x = 0.0;
    var t = 0.0;
    for k in [corner(i, j), corner(i + 1, j), corner(i + 1, j + 1), corner(i, j + 1)] {
        x += 0.25 * px_lat[k];
        t += 0.25 * pt_lat[k];
    }
    return shell_point(x, t);
}

fn make_cells(mesh_g: i64, mats: i64[:]) {
    px_lat.clear();
    pt_lat.clear();
    for b in lat_nt() + 1 {
        for a in lat_nx() + 1 {
            let x, t = default_place(a, b);
            px_lat.push(x);
            pt_lat.push(t);
        }
    }
    snaps.clear();
    snaps.push(mesh_g);
    for g in mats { if g >= 0 && position(snaps[..], g) < 0 { snaps.push(g); } }
    var moved: u8[>..] = [];
    for j in nt_cells() + 1 { for i in nx_cells() + 1 { moved.push(if snap_corner(i, j) { 1 } else { 0 }); } }
    for pass in 10 {
        var undone = 0;
        for j in nt_cells() {
            for i in nx_cells() {
                if cell_quality(i, j) > 0.1 { continue; }
                for c in [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]] {
                    let m = c[0] + (nx_cells() + 1) * c[1];
                    if moved[m] == 0 { continue; }
                    let k = corner(c[0], c[1]);
                    let x, t = default_place(order * c[0], order * c[1]);
                    px_lat[k] = x;
                    pt_lat[k] = t;
                    moved[m] = 0;
                    undone++;
                }
            }
        }
        if undone == 0 { break; }
    }
    cell_mat.clear();
    for j in nt_cells() {
        for i in nx_cells() {
            let c = cell_centre(i, j);
            var m = -1;
            if phi(mesh_g, c) < 0.0 {
                for g, k in mats { if m < 0 && (g < 0 || phi(g, c) < 0.0) { m = k; } }
                if m < 0 { m = -2; }
            }
            cell_mat.push(m);
        }
    }
    if order == 2 { place_middles(); }
}

fn place_middles() {
    for b in lat_nt() + 1 {
        for a in lat_nx() + 1 {
            if a % 2 == 0 && b % 2 == 0 { continue; }
            let k = lat(a, b);
            if a % 2 == 1 && b % 2 == 1 { continue; }
            let ka = if a % 2 == 1 { lat(a - 1, b) } else { lat(a, b - 1) };
            let kb = if a % 2 == 1 { lat(a + 1, b) } else { lat(a, b + 1) };
            px_lat[k] = 0.5 * (px_lat[ka] + px_lat[kb]);
            pt_lat[k] = 0.5 * (pt_lat[ka] + pt_lat[kb]);
            for g in snaps {
                let tol = 1.0e-9 * mesh_size;
                if abs(cut_phi(g, lat_point(ka))) > tol || abs(cut_phi(g, lat_point(kb))) > tol { continue; }
                let xr, tr, found = root_near(g, px_lat[k], pt_lat[k], mesh_size);
                if found {
                    px_lat[k] = xr;
                    pt_lat[k] = tr;
                }
                break;
            }
        }
    }
    for b in lat_nt() {
        for a in lat_nx() {
            if a % 2 == 0 || b % 2 == 0 { continue; }
            var x = 0.0;
            var t = 0.0;
            for d in [[0, -1], [1, 0], [0, 1], [-1, 0]] {
                x += 0.5 * px_lat[lat(a + d[0], b + d[1])];
                t += 0.5 * pt_lat[lat(a + d[0], b + d[1])];
            }
            for d in [[-1, -1], [1, -1], [1, 1], [-1, 1]] {
                x -= 0.25 * px_lat[lat(a + d[0], b + d[1])];
                t -= 0.25 * pt_lat[lat(a + d[0], b + d[1])];
            }
            px_lat[lat(a, b)] = x;
            pt_lat[lat(a, b)] = t;
        }
    }
}

fn kept(i: i64, j: i64) -> bool {
    return i >= 0 && j >= 0 && i < nx_cells() && j < nt_cells() && cell_mat[i + nx_cells() * j] >= 0;
}

struct Edge {
    a: i64,
    m: i64,
    b: i64,
    nx: f64,
    ns: f64,
}

var edges: Edge[>..] = [];

fn make_edges() {
    edges.clear();
    for j in nt_cells() {
        for i in nx_cells() {
            if !kept(i, j) { continue; }
            let sides = [[i, j, i + 1, j, 0, -1], [i + 1, j, i + 1, j + 1, 1, 0], [i + 1, j + 1, i, j + 1, 0, 1], [i, j + 1, i, j, -1, 0]];
            for s in sides {
                if kept(i + s[4], j + s[5]) { continue; }
                let a = corner(s[0], s[1]);
                let b = corner(s[2], s[3]);
                let m = if order == 2 { lat(s[0] + s[2], s[1] + s[3]) } else { -1 };
                let dx = px_lat[b] - px_lat[a];
                let ds = radius * (pt_lat[b] - pt_lat[a]);
                let l = sqrt(dx * dx + ds * ds);
                edges.push(Edge { a: a, m: m, b: b, nx: ds / l, ns: -dx / l });
            }
        }
    }
}

fn edge_mid(e: Edge) -> f64[3] {
    if e.m >= 0 { return lat_point(e.m); }
    return shell_point(0.5 * (px_lat[e.a] + px_lat[e.b]), 0.5 * (pt_lat[e.a] + pt_lat[e.b]));
}

fn edge_length(e: Edge) -> f64 {
    return sqrt(sq(px_lat[e.b] - px_lat[e.a]) + sq(radius * (pt_lat[e.b] - pt_lat[e.a])));
}

fn edge_normal(e: Edge) -> f64[3] {
    let t = 0.5 * (pt_lat[e.a] + pt_lat[e.b]);
    return [e.nx, e.ns * cos(t), -e.ns * sin(t)];
}

recursive fn edge_on(v: i64, e: Edge) -> bool {
    let tol = geo_tol();
    if kind(v) == syntax::N_NAME { return edge_on(regions[find_region(v)].node, e); }
    if kind(v) == syntax::N_BINARY && op(v) == "|" { return edge_on(child(v, 0), e) || edge_on(child(v, 1), e); }
    if kind(v) == syntax::N_BINARY {
        let f = phi(node_shape(child(v, 1)), edge_mid(e));
        return edge_on(child(v, 0), e) && (if op(v) == "&" { f <= tol } else { f >= -tol });
    }
    let c = arg(v, 1);
    if is_condition(c) {
        return holds(c, lat_point(e.a), tol) && holds(c, lat_point(e.b), tol) && holds(c, edge_mid(e), tol);
    }
    let g = node_shape(c);
    let near = 1.0 * edge_length(e);
    return abs(phi(g, lat_point(e.a))) <= near && abs(phi(g, lat_point(e.b))) <= near && abs(phi(g, edge_mid(e))) <= near;
}
