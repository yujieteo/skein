// The 2D triangle mesher (spec §9.2): Delaunay refinement after Ruppert
// (1995) and Shewchuk (2002), on a boundary sampled from the CSG tree.
//
// 1. Each boundary curve of a primitive (a line of a rect or a polygon, a
//    circle, a span of the spline of a fill) is sampled. A sample is on the boundary of the whole shape
//    when φ changes sign across the curve there. The ends of each piece on
//    the boundary are found by bisection, then moved to the nearest true
//    intersection with another curve; they are the corners and stay
//    fixed. Each piece is cut into segments no longer than `size`. The
//    boundary of each interface shape inside the shape is sampled the same
//    way, as inner segments, so the triangles do not cross it.
// 2. The Bowyer-Watson algorithm inserts the points into a Delaunay
//    triangulation inside a large triangle.
// 3. A segment is encroached if a point lies inside its diametral circle.
//    Each encroached segment is split at the midpoint of its parameter on
//    the true curve, or on a concentric shell next to a corner below
//    60 degrees, until all segments are edges of the triangulation.
// 4. Each triangle inside the polygon of the outer segments with an angle below
//    MIN_ANGLE, or a circumradius above size / sqrt(3), gets a point at
//    its circumcentre. With the targets of refine.goose, the wanted size
//    at the circumcentre takes the place of size.
//    If that point encroaches a segment, the segment is split instead. A
//    triangle whose smallest angle is at a corner of 2 segments below
//    MIN_ANGLE cannot reach it, so only its circumradius counts.
//
// `size` is the size far from each `refine` (spec §9.1). Near a refined
// segment or circle, the size is the refine size plus (grade - 1) times
// the distance, so that the size of two neighbour cells changes by at most
// the factor `grade`.
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .refine;

let MIN_ANGLE = 25.0;     // degrees; Ruppert's proof needs at most 20.7

let LINE = 0;
let ARC = 1;
let SPAN = 2;

// A boundary curve of a primitive. A line goes from (p0, p1) to (p2, p3).
// An arc is a full circle of centre (p0, p1) and radius p2, from angle 0.
// A span is span p3 of fill node p2, moved by (p0, p1).
struct Curve {
    kind: i64,
    p: f64[4],
}

var curves: Curve[>..] = [];

// The points: x and y of each. Points 0, 1 and 2 are the corners of the
// large triangle.
var px: f64[>..] = [];
var py: f64[>..] = [];
var vtri: i64[>..] = [];        // a live triangle at each point
var corner: f64[>..] = [];      // the smallest angle of 2 segments at the point, degrees

// The triangles: 3 points, counter-clockwise, and the 3 neighbours. The
// neighbour k is across the edge opposite point k, or -1.
var tv: i64[>..] = [];
var tn: i64[>..] = [];
var talive: bool[>..] = [];
var tfree: i64[>..<] = [];
var tmark: i64[>..] = [];       // the cavity stamp
var tin: bool[>..] = [];        // inside the polygon of the segments, by classify
var stamp = 0;
var last_tri = 0;

// The segments of the boundary: the 2 points, the curve and its parameter
// at each point.
struct Segment {
    a: i64,
    b: i64,
    curve: i64,
    t0: f64,
    t1: f64,
    alive: bool,
    inner: bool,    // on the boundary of an interface, inside the shape
}

var segments: Segment[>..] = [];
var interfaces: i64[>..] = [];  // the shapes of add_interface
var fan: i64[>..] = [];         // scratch of insert, by point

// --- the size --------------------------------------------------------------------

// A `refine` on a segment from (p0, p1) to (p2, p3) (LINE), or on a circle
// of centre (p0, p1) and radius p2 (ARC), with the size h there.
struct Refine {
    kind: i64,
    p: f64[4],
    h: f64,
}

var refines: Refine[>..] = [];
var grade = 1.3;

// The calls before mesh_2d that set the refines and the grade. mesh_2d
// keeps them.
fn clear_refines() {
    refines.clear();
    grade = 1.3;
}

fn set_grade(g: f64) {
    guard g > 1.0 else { abort("set_grade: g > 1"); }
    grade = g;
}

fn refine_segment(x0: f64, y0: f64, x1: f64, y1: f64, h: f64) {
    guard h > 0.0 else { abort("refine_segment: h > 0"); }
    refines.push(Refine { kind: LINE, p: [x0, y0, x1, y1], h: h });
}

fn refine_circle(cx: f64, cy: f64, r: f64, h: f64) {
    guard h > 0.0 && r > 0.0 else { abort("refine_circle: h > 0 and r > 0"); }
    refines.push(Refine { kind: ARC, p: [cx, cy, r, 0.0], h: h });
}

// The distance from (x, y) to refine k.
fn refine_distance(k: i64, x: f64, y: f64) -> f64 {
    let p = refines[k].p;
    if refines[k].kind == ARC { return abs(sqrt((x - p[0]) * (x - p[0]) + (y - p[1]) * (y - p[1])) - p[2]); }
    let ux = p[2] - p[0];
    let uy = p[3] - p[1];
    let l2 = ux * ux + uy * uy;
    let t = if l2 > 0.0 { clamp(((x - p[0]) * ux + (y - p[1]) * uy) / l2, 0.0, 1.0) } else { 0.0 };
    let dx = x - p[0] - t * ux;
    let dy = y - p[1] - t * uy;
    return sqrt(dx * dx + dy * dy);
}

// The target edge length at (x, y).
fn local_size(size: f64, x: f64, y: f64) -> f64 {
    var h = size;
    for k in refines.len { h = min(h, refines[k].h + (grade - 1.0) * refine_distance(k, x, y)); }
    return h;
}

// The number of cells of edge local_size in the box b of g: the integral
// of 1 / h² on a lattice of 200 × 200 points.
fn cell_estimate(b: f64[6], size: f64) -> f64 {
    let w = b[3] - b[0];
    let h = b[4] - b[1];
    var sum = 0.0;
    for i in 200 {
        for j in 200 {
            let s = local_size(size, b[0] + (i + 0.5) * w / 200.0, b[1] + (j + 0.5) * h / 200.0);
            sum += 1.0 / (s * s);
        }
    }
    return sum * w * h / 40000.0;
}

// --- the boundary ----------------------------------------------------------------

// The boundary curves of g, each moved by (tx, ty): the sum of the moves
// of the patterns above g. The reader has no pattern_polar of a 2D shape
// (verify.goose).
fn collect_curves(g: i64, tx: f64, ty: f64) -> bool {
    // The nodes to walk, each with its move: a stack, with no recursion (rule 1).
    var todo: i64[>..] = [g];
    var moves: f64[>..] = [tx, ty];
    while todo.len > 0 {
        let h = todo[todo.len - 1];
        let mx = moves[moves.len - 2];
        let my = moves[moves.len - 1];
        todo.pop();
        moves.pop();
        moves.pop();
        let s = shapes[h];
        if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT {
            todo.append([s.b, s.a]);
            moves.append([mx, my, mx, my]);
        } else if s.kind == PATTERN {
            for i in s.b {
                let k = s.b - 1 - i;
                todo.push(s.a);
                moves.append([mx + k * s.p[0], my + k * s.p[1]]);
            }
        } else if !own_curves(h, mx, my) {
            return false;
        }
    }
    return true;
}

// The boundary curves of node h with no parts, moved by (tx, ty). False
// if h has no curves of this mesher.
fn own_curves(g: i64, tx: f64, ty: f64) -> bool {
    let s = shapes[g];
    let p = s.p;
    fn line(x0: f64, y0: f64, x1: f64, y1: f64) {
        curves.push(Curve { kind: LINE, p: [x0 + tx, y0 + ty, x1 + tx, y1 + ty] });
    }
    if s.kind == RECT {
        let xs: f64[5] = [p[0], p[2], p[2], p[0], p[0]];
        let ys: f64[5] = [p[1], p[1], p[3], p[3], p[1]];
        for i in 4 { line(xs[i], ys[i], xs[i + 1], ys[i + 1]); }
        return true;
    }
    if s.kind == CIRCLE {
        curves.push(Curve { kind: ARC, p: [p[0] + tx, p[1] + ty, p[2], 0.0] });
        return true;
    }
    if s.kind == POLYGON {
        for i in s.b {
            let j = (i + 1) % s.b;
            line(polygon_points[2 * (s.a + i)], polygon_points[2 * (s.a + i) + 1],
                 polygon_points[2 * (s.a + j)], polygon_points[2 * (s.a + j) + 1]);
        }
        return true;
    }
    if s.kind == FILL {
        for i in s.b - 1 { curves.push(Curve { kind: SPAN, p: [tx, ty, g, i] }); }
        let a = fill_point(g, s.b - 1);
        let b = fill_point(g, 0);
        if a[0] != b[0] || a[1] != b[1] { line(a[0], a[1], b[0], b[1]); }
        return true;
    }
    return false;
}

// The cubics of span curve c (geometry.goose).
fn span_of(c: i64) -> f64[12] { span(curves[c].p[2] as i64, curves[c].p[3] as i64) }

fn curve_point(c: i64, t: f64) -> f64, f64 {
    let p = curves[c].p;
    if curves[c].kind == LINE { return p[0] + t * (p[2] - p[0]), p[1] + t * (p[3] - p[1]); }
    if curves[c].kind == SPAN {
        let e = span_of(c);
        let x, dx, ddx = cubic(e, 0, t);
        let y, dy, ddy = cubic(e, 1, t);
        return p[0] + x, p[1] + y;
    }
    return p[0] + p[2] * cos(TAU * t), p[1] + p[2] * sin(TAU * t);
}

// The length of curve c. A span has the length of 16 chords.
fn curve_length(c: i64) -> f64 {
    let p = curves[c].p;
    if curves[c].kind == LINE { return sqrt((p[2] - p[0]) * (p[2] - p[0]) + (p[3] - p[1]) * (p[3] - p[1])); }
    if curves[c].kind == SPAN {
        var l = 0.0;
        var x0, y0 = curve_point(c, 0.0);
        for k in 16 {
            let x1, y1 = curve_point(c, (k + 1) / 16.0);
            l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
            x0 = x1;
            y0 = y1;
        }
        return l;
    }
    return TAU * p[2];
}

// True if an earlier curve holds the point of curve c at t: a line on the
// same line, the same circle, or the same span. That curve then has the
// boundary there.
fn covered(c: i64, t: f64, eps: f64) -> bool {
    let x, y = curve_point(c, t);
    let p = curves[c].p;
    for e in c {
        guard curves[e].kind == curves[c].kind else { continue; }
        let q = curves[e].p;
        if curves[c].kind == SPAN {
            if q[2] == p[2] && q[3] == p[3] && abs(q[0] - p[0]) <= eps && abs(q[1] - p[1]) <= eps { return true; }
            continue;
        }
        if curves[c].kind == ARC {
            if abs(q[0] - p[0]) <= eps && abs(q[1] - p[1]) <= eps && abs(q[2] - p[2]) <= eps { return true; }
            continue;
        }
        let l = curve_length(e);
        let dx = (q[2] - q[0]) / l;
        let dy = (q[3] - q[1]) / l;
        let off0 = (p[1] - q[1]) * dx - (p[0] - q[0]) * dy;
        let off1 = (p[3] - q[1]) * dx - (p[2] - q[0]) * dy;
        let along = (x - q[0]) * dx + (y - q[1]) * dy;
        if abs(off0) <= eps && abs(off1) <= eps && along >= 0.0 && along <= l { return true; }
    }
    return false;
}

fn clear_interfaces() { interfaces.clear(); }

// Keeps the boundary of shape h inside the shape of mesh_2d as edges of
// the triangles: a boundary between two materials, or of a region that a
// step kills or gives birth to (spec §14.3).
fn add_interface(h: i64) {
    guard h >= 0 && h < shapes.len else { abort("add_interface: a shape"); }
    interfaces.push(h);
}

// True if curve c has the boundary of shape g at t, and no earlier curve
// holds the point. With `inner`, true if c has the boundary of an
// interface there, and g is on both sides.
fn on_boundary(g: i64, c: i64, t: f64, eps: f64, inner: bool) -> bool {
    if covered(c, t, eps) { return false; }
    let a, b = sides(g, c, t, eps);
    guard inner else { return a != b; }
    guard a && b else { return false; }
    for h in interfaces {
        let p, q = sides(h, c, t, eps);
        if p != q { return true; }
    }
    return false;
}

// Whether the points at eps on the two sides of curve c at t are in g.
fn sides(g: i64, c: i64, t: f64, eps: f64) -> bool, bool {
    let x, y = curve_point(c, t);
    let p = curves[c].p;
    var nx = 0.0;
    var ny = 0.0;
    if curves[c].kind == LINE {
        let l = curve_length(c);
        nx = (p[3] - p[1]) / l;
        ny = (p[0] - p[2]) / l;
    } else if curves[c].kind == SPAN {
        let e = span_of(c);
        let x0, dx, ddx = cubic(e, 0, t);
        let y0, dy, ddy = cubic(e, 1, t);
        let l = sqrt(dx * dx + dy * dy);
        nx = dy / l;
        ny = -dx / l;
    } else {
        nx = cos(TAU * t);
        ny = sin(TAU * t);
    }
    return phi(g, x + eps * nx, y + eps * ny, 0.0) < 0.0, phi(g, x - eps * nx, y - eps * ny, 0.0) < 0.0;
}

// The points where the line through (x0, y0) along the unit vector
// (ux, uy) meets the circle of centre (cx, cy) and radius r. A line within
// eps of a tangent touches the circle at one point.
fn line_circle(x0: f64, y0: f64, ux: f64, uy: f64, cx: f64, cy: f64, r: f64, eps: f64, out: f64[>..]&) {
    let b = (x0 - cx) * ux + (y0 - cy) * uy;
    let h = abs((x0 - cx) * uy - (y0 - cy) * ux);
    guard h <= r + eps else { return; }
    let w = sqrt(max(r * r - h * h, 0.0));
    for sign in [-1.0, 1.0] {
        out.push(x0 + (sign * w - b) * ux);
        out.push(y0 + (sign * w - b) * uy);
    }
}

// The parameter of curve c at the intersection with another curve that is
// nearest to the point at t, or t if no intersection is within `reach`.
// The ends of a piece from the bisection of on_boundary are only near the
// true corner: up to eps / tan(angle), or sqrt(eps r) at a tangent. A span
// has no intersection here: the end from the bisection stays.
fn snap_end(c: i64, t: f64, eps: f64, reach: f64) -> f64 {
    let x, y = curve_point(c, t);
    let p = curves[c].p;
    var points: f64[>..] = [];
    for e in curves.len {
        guard e != c && curves[c].kind != SPAN && curves[e].kind != SPAN else { continue; }
        let q = curves[e].p;
        if curves[c].kind == LINE && curves[e].kind == LINE {
            let dx = p[2] - p[0];
            let dy = p[3] - p[1];
            let ex = q[2] - q[0];
            let ey = q[3] - q[1];
            let d = dx * ey - dy * ex;
            guard abs(d) > 1.0e-12 * (dx * dx + dy * dy + ex * ex + ey * ey) else { continue; }
            let s = ((q[0] - p[0]) * ey - (q[1] - p[1]) * ex) / d;
            points.push(p[0] + s * dx);
            points.push(p[1] + s * dy);
        } else if curves[c].kind == LINE || curves[e].kind == LINE {
            let l = if curves[c].kind == LINE { p } else { q };
            let a = if curves[c].kind == ARC { p } else { q };
            let n = sqrt((l[2] - l[0]) * (l[2] - l[0]) + (l[3] - l[1]) * (l[3] - l[1]));
            line_circle(l[0], l[1], (l[2] - l[0]) / n, (l[3] - l[1]) / n, a[0], a[1], a[2], eps, points);
        } else {
            let dx = q[0] - p[0];
            let dy = q[1] - p[1];
            let d = sqrt(dx * dx + dy * dy);
            guard d > eps else { continue; }
            let along = (d * d + p[2] * p[2] - q[2] * q[2]) / (2.0 * d);
            line_circle(p[0] + along * dx / d, p[1] + along * dy / d, -dy / d, dx / d, p[0], p[1], p[2], eps, points);
        }
    }
    var best = t;
    var nearest = reach;
    for i in points.len / 2 {
        let d = sqrt((points[2 * i] - x) * (points[2 * i] - x) + (points[2 * i + 1] - y) * (points[2 * i + 1] - y));
        guard d < nearest else { continue; }
        nearest = d;
        if curves[c].kind == LINE {
            let dx = p[2] - p[0];
            let dy = p[3] - p[1];
            best = ((points[2 * i] - p[0]) * dx + (points[2 * i + 1] - p[1]) * dy) / (dx * dx + dy * dy);
        } else {
            let a = atan2(points[2 * i + 1] - p[1], points[2 * i] - p[0]) / TAU;
            best = a + round(t - a);
        }
    }
    return best;
}

// The point at (x, y), or a new one. `tol` merges the ends of two pieces:
// the bisection puts an end up to eps / tan(angle) from the true corner.
fn boundary_point(x: f64, y: f64, tol: f64) -> i64 {
    for i in 3..px.len {
        if abs(px[i] - x) <= tol && abs(py[i] - y) <= tol { return i; }
    }
    px.push(x);
    py.push(y);
    vtri.push(-1);
    return px.len - 1;
}

// Cuts curve c from ta to tb into segments no longer than local_size:
// the cuts are at equal steps of the integral of 1 / local_size, on 64
// samples for each segment of `size`.
fn add_piece(c: i64, ta: f64, tb: f64, size: f64, tol: f64, inner: bool) {
    let length = curve_length(c) * (tb - ta);
    let m = max(64 * ceil(length / size) as i64, 64);
    var steps: f64[>..] = [0.0];     // the integral at sample k
    for k in m {
        let x, y = curve_point(c, ta + (tb - ta) * (k + 0.5) / m);
        steps.push(steps[k] + length / m / local_size(size, x, y));
    }
    let n = max(ceil(steps[m] - 1.0e-9) as i64, if curves[c].kind == ARC { 3 } else { 1 });
    var prev = -1;
    var tprev = ta;
    var k = 0;
    for i in n + 1 {
        // the parameter where the integral is steps[m] i / n
        let want = steps[m] * i / n;
        while k < m - 1 && steps[k + 1] < want { k++; }
        let part = if steps[k + 1] > steps[k] { clamp((want - steps[k]) / (steps[k + 1] - steps[k]), 0.0, 1.0) } else { 0.0 };
        let t = if i == n { tb } else { ta + (tb - ta) * (k + part) / m };
        let x, y = curve_point(c, t);
        let v = boundary_point(x, y, tol);
        if i > 0 && v != prev {
            segments.push(Segment { a: prev, b: v, curve: c, t0: tprev, t1: t, alive: true, inner: inner });
        }
        prev = v;
        tprev = t;
    }
}

// The state of curve c at t: 0 off the boundary (on_boundary), else 1,
// plus 2^(k + 1) for each interface k that holds a point at eps from the
// boundary of g there. Reason: an interface that ends on the boundary of
// g (the flag of examples/flow/turek_hron_fsi1_tri.skn on its cylinder)
// must end at a point of the boundary, or the refinement does not finish.
fn piece_state(g: i64, c: i64, t: f64, eps: f64, inner: bool) -> i64 {
    guard on_boundary(g, c, t, eps, inner) else { return 0; }
    guard !inner else { return 1; }
    var state = 1;
    for h, k in interfaces {
        let p, q = sides(h, c, t, eps);
        if p || q { state += 2 << k; }
    }
    return state;
}

// The pieces of each curve that are on the boundary of g, or with `inner`
// on the boundary of an interface inside g. A piece ends where the state
// of piece_state changes.
fn sample_boundary(g: i64, size: f64, diameter: f64, inner: bool) {
    let eps = 1.0e-7 * diameter;
    let tol = 100.0 * eps;
    let reach = 10.0 * sqrt(eps * diameter);
    for c in curves.len {
        let m = max(512, ceil(16.0 * curve_length(c) / size) as i64);
        var on: i64[>..] = [];
        for i in m { on.push(piece_state(g, c, (i + 0.5) / m, eps, inner)); }
        // The parameter between sample i and sample i + 1 where the state changes.
        fn change(i: i64) -> f64 {
            var lo = (i + 0.5) / m;
            var hi = (i + 1.5) / m;
            let state = on[i];
            for k in 60 {
                let mid = 0.5 * (lo + hi);
                if piece_state(g, c, mid, eps, inner) == state { lo = mid; } else { hi = mid; }
            }
            return snap_end(c, 0.5 * (lo + hi), eps, reach);
        }
        if curves[c].kind != ARC {
            var i = 0;
            while i < m {
                if on[i] == 0 { i++; continue; }
                var j = i;
                while j + 1 < m && on[j + 1] == on[i] { j++; }
                let ta = if i == 0 { 0.0 } else { change(i - 1) };
                let tb = if j == m - 1 { 1.0 } else { change(j) };
                add_piece(c, ta, tb, size, tol, inner);
                i = j + 1;
            }
            continue;
        }
        // A circle: start at a change, if there is one.
        var start = -1;
        for i in m {
            if on[i] != on[(i + 1) % m] { start = i; break; }
        }
        if start < 0 {
            if on[0] != 0 { add_piece(c, 0.0, 1.0, size, tol, inner); }
            continue;
        }
        var i = start + 1;
        while i < start + 1 + m {
            if on[i % m] == 0 { i++; continue; }
            var j = i;
            while j + 1 < start + 1 + m && on[(j + 1) % m] == on[i % m] { j++; }
            let ta = change((i - 1) % m) + ((i - 1) / m);
            let tb = change(j % m) + (j / m);
            add_piece(c, ta, tb, size, tol, inner);
            i = j + 1;
        }
    }
}

// Sets `corner` of each point of the sampled boundary.
fn find_corners() {
    for v in px.len { corner.push(180.0); }
    for s in segments.len {
        for e in s {
            for i in 2 {
                let v = if i == 0 { segments[s].a } else { segments[s].b };
                let a = if v == segments[s].a { segments[s].b } else { segments[s].a };
                let b = if v == segments[e].a { segments[e].b } else if v == segments[e].b { segments[e].a } else { -1 };
                guard b >= 0 && b != a else { continue; }
                let ax = px[a] - px[v];
                let ay = py[a] - py[v];
                let bx = px[b] - px[v];
                let by = py[b] - py[v];
                let c = (ax * bx + ay * by) / sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
                corner[v] = min(corner[v], degrees(acos(max(-1.0, min(1.0, c)))));
            }
        }
    }
}

// --- predicates ------------------------------------------------------------------

fn orient(a: i64, b: i64, x: f64, y: f64) -> f64 {
    return (px[b] - px[a]) * (y - py[a]) - (py[b] - py[a]) * (x - px[a]);
}

// > 0 if (x, y) is inside the circumcircle of the counter-clockwise
// triangle t.
fn incircle(t: i64, x: f64, y: f64) -> f64 {
    let a = tv[3 * t];
    let b = tv[3 * t + 1];
    let c = tv[3 * t + 2];
    let adx = px[a] - x;
    let ady = py[a] - y;
    let bdx = px[b] - x;
    let bdy = py[b] - y;
    let cdx = px[c] - x;
    let cdy = py[c] - y;
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
         - (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady)
         + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

// --- the triangulation -------------------------------------------------------------

fn new_triangle(a: i64, b: i64, c: i64) -> i64 {
    var t = 0;
    if tfree.len > 0 {
        t = tfree.pop();
    } else {
        t = talive.len;
        for k in 3 {
            tv.push(0);
            tn.push(-1);
        }
        talive.push(false);
        tmark.push(0);
        tin.push(false);
    }
    tv[3 * t] = a;
    tv[3 * t + 1] = b;
    tv[3 * t + 2] = c;
    for k in 3 { tn[3 * t + k] = -1; }
    talive[t] = true;
    tin[t] = false;
    vtri[a] = t;
    vtri[b] = t;
    vtri[c] = t;
    last_tri = t;
    return t;
}

fn kill_triangle(t: i64) {
    talive[t] = false;
    tfree.push(t);
}

// The triangle that holds (x, y), by a walk from the last new triangle.
fn locate(x: f64, y: f64) -> i64 {
    var t = last_tri;
    if !talive[t] {
        for i in talive.len {
            if talive[i] { t = i; break; }
        }
    }
    for step in 4 * talive.len + 16 {
        var moved = false;
        for k in 3 {
            let a = tv[3 * t + (k + 1) % 3];
            let b = tv[3 * t + (k + 2) % 3];
            if orient(a, b, x, y) < 0.0 && tn[3 * t + k] >= 0 {
                t = tn[3 * t + k];
                moved = true;
                break;
            }
        }
        if !moved { return t; }
    }
    // The walk cycled: search all triangles.
    for i in talive.len {
        guard talive[i] else { continue; }
        var inside = true;
        for k in 3 {
            if orient(tv[3 * i + (k + 1) % 3], tv[3 * i + (k + 2) % 3], x, y) < 0.0 { inside = false; }
        }
        if inside { return i; }
    }
    return t;
}

// Inserts the point (x, y) by Bowyer-Watson and returns its index.
fn insert_point(x: f64, y: f64) -> i64 {
    let t0 = locate(x, y);
    px.push(x);
    py.push(y);
    vtri.push(-1);
    let v = px.len - 1;
    insert_existing(v, t0);
    return v;
}

// The cavity of point v from triangle t0: the triangles whose circumcircle
// holds v, and that keep the cavity star-shaped from v. They have tmark
// stamp.
fn cavity_of(v: i64, t0: i64) -> i64[>..] {
    let x = px[v];
    let y = py[v];
    stamp++;
    var cavity: i64[>..] = [t0];
    tmark[t0] = stamp;
    var i = 0;
    while i < cavity.len {
        let t = cavity[i];
        for k in 3 {
            let n = tn[3 * t + k];
            if n >= 0 && tmark[n] != stamp && incircle(n, x, y) > 0.0 {
                tmark[n] = stamp;
                cavity.push(n);
            }
        }
        i++;
    }
    // The cavity must be star-shaped from the point: a triangle with a
    // boundary edge that does not face the point leaves the cavity.
    loop {
        var removed = false;
        for t in cavity {
            if t == t0 || tmark[t] != stamp { continue; }
            for k in 3 {
                let n = tn[3 * t + k];
                if n >= 0 && tmark[n] == stamp { continue; }
                if orient(tv[3 * t + (k + 1) % 3], tv[3 * t + (k + 2) % 3], x, y) <= 0.0 {
                    tmark[t] = stamp - 1;
                    removed = true;
                    break;
                }
            }
        }
        if !removed { break; }
    }
    return cavity;
}

fn insert_existing(v: i64, t0: i64) {
    let cavity = cavity_of(v, t0);
    // The new triangles: one for each boundary edge of the cavity.
    var made: i64[>..] = [];
    for t in cavity {
        if tmark[t] != stamp { continue; }
        for k in 3 {
            let n = tn[3 * t + k];
            if n >= 0 && tmark[n] == stamp { continue; }
            let a = tv[3 * t + (k + 1) % 3];
            let b = tv[3 * t + (k + 2) % 3];
            made.push(a);
            made.push(b);
            made.push(n);
        }
    }
    for t in cavity {
        if tmark[t] == stamp { kill_triangle(t); }
    }
    while fan.len < px.len { fan.push(-1); }
    var news: i64[>..] = [];
    for e in made.len / 3 {
        let a = made[3 * e];
        let b = made[3 * e + 1];
        let n = made[3 * e + 2];
        let nt = new_triangle(a, b, v);
        news.push(nt);
        tn[3 * nt + 2] = n;
        if n >= 0 {
            for k in 3 {
                let na = tv[3 * n + (k + 1) % 3];
                let nb = tv[3 * n + (k + 2) % 3];
                if na == b && nb == a { tn[3 * n + k] = nt; }
            }
        }
        fan[a] = nt;
    }
    // Triangle (a, b, v): across from a is (b, v), the triangle that starts
    // at b. Across from b is (v, a), the triangle that ends at a.
    for nt in news {
        let b = tv[3 * nt + 1];
        let next = fan[b];
        tn[3 * nt] = next;
        tn[3 * next + 1] = nt;
    }
    for nt in news { fan[tv[3 * nt]] = -1; }
}

// The triangle and the local index k such that the edge from a to b is
// the edge opposite k, or -1, -1. a is not a corner of the large triangle.
fn find_edge(a: i64, b: i64) -> i64, i64 {
    let first = vtri[a];
    var t = first;
    for step in px.len + 3 {
        var ia = 0;
        for k in 3 {
            if tv[3 * t + k] == a { ia = k; }
        }
        let next = tv[3 * t + (ia + 1) % 3];
        if next == b { return t, (ia + 2) % 3; }
        // Turn about a: the neighbour across the edge from a to `next`.
        t = tn[3 * t + (ia + 2) % 3];
        if t < 0 || t == first { break; }
    }
    return -1, -1;
}

fn encroaches(s: i64, x: f64, y: f64) -> bool {
    let a = segments[s].a;
    let b = segments[s].b;
    return (px[a] - x) * (px[b] - x) + (py[a] - y) * (py[b] - y) < 0.0;
}

// True if segment s is not an edge, or a point next to it is inside its
// diametral circle.
fn encroached(s: i64) -> bool {
    let a = segments[s].a;
    let b = segments[s].b;
    let t, k = find_edge(a, b);
    if t < 0 { return true; }
    let c = tv[3 * t + k];
    if encroaches(s, px[c], py[c]) { return true; }
    let n = tn[3 * t + k];
    if n >= 0 {
        for j in 3 {
            let d = tv[3 * n + j];
            if d != a && d != b && encroaches(s, px[d], py[d]) { return true; }
        }
    }
    return false;
}

// True if point v is a corner of 2 segments at less than `angle` degrees.
fn sharp(v: i64, angle: f64) -> bool { return v < corner.len && corner[v] < angle; }

// Splits segment s at its midpoint. A segment with one end at a corner
// below 60 degrees, where midpoint splits of the 2 segments can encroach
// each other without end, is split on a concentric shell: at a distance
// of a power of 2 from the corner. Then the 2 segments at the corner get
// the same length, and the split of one does not encroach the other.
fn split_segment(s: i64) {
    let seg = segments[s];
    let sharp_a = sharp(seg.a, 60.0);
    let sharp_b = sharp(seg.b, 60.0);
    var f = 0.5;
    if sharp_a != sharp_b {
        let l = sqrt((px[seg.b] - px[seg.a]) * (px[seg.b] - px[seg.a]) + (py[seg.b] - py[seg.a]) * (py[seg.b] - py[seg.a]));
        let shell = pow(2.0, round(log2(0.5 * l))) / l;
        f = if sharp_a { shell } else { 1.0 - shell };
    }
    let tm = seg.t0 + f * (seg.t1 - seg.t0);
    let x, y = curve_point(seg.curve, tm);
    let v = insert_point(x, y);
    segments[s].alive = false;
    segments.push(Segment { a: seg.a, b: v, curve: seg.curve, t0: seg.t0, t1: tm, alive: true, inner: seg.inner });
    segments.push(Segment { a: v, b: seg.b, curve: seg.curve, t0: tm, t1: seg.t1, alive: true, inner: seg.inner });
}

// Splits encroached segments until there are none. False if the point
// budget runs out.
fn recover_segments(budget: i64) -> bool {
    loop {
        var split = false;
        for s in segments.len {
            if segments[s].alive && encroached(s) {
                if px.len >= budget { return false; }
                split_segment(s);
                split = true;
            }
        }
        if !split { return true; }
    }
}

// --- refinement ------------------------------------------------------------------

// Sets tin of each live triangle from the polygon of the segments, not
// from φ: between a chord and its arc, φ has the wrong sign. A flood fill
// starts at the triangles on the corners of the large triangle, which are
// outside, and flips the side across each outer segment. The segments must be
// edges of the triangulation.
fn classify() {
    var seg: bool[>..] = [];
    var seen: bool[>..] = [];
    for e in tv.len { seg.push(false); }
    for t in talive.len {
        tin[t] = false;
        seen.push(false);
    }
    for s in segments.len {
        guard segments[s].alive && !segments[s].inner else { continue; }
        let t, k = find_edge(segments[s].a, segments[s].b);
        if t < 0 { continue; }
        seg[3 * t + k] = true;
        let n = tn[3 * t + k];
        if n >= 0 {
            for j in 3 {
                if tn[3 * n + j] == t { seg[3 * n + j] = true; }
            }
        }
    }
    var queue: i64[>..] = [];
    for t in talive.len {
        if talive[t] && (tv[3 * t] < 3 || tv[3 * t + 1] < 3 || tv[3 * t + 2] < 3) {
            seen[t] = true;
            queue.push(t);
        }
    }
    var i = 0;
    while i < queue.len {
        let t = queue[i];
        for k in 3 {
            let n = tn[3 * t + k];
            if n < 0 || seen[n] { continue; }
            seen[n] = true;
            tin[n] = tin[t] != seg[3 * t + k];
            queue.push(n);
        }
        i++;
    }
}

// The circumcentre of t, the ratio of its circumradius to its shortest
// edge, and the point opposite that edge, at the smallest angle.
fn circumcentre(t: i64) -> f64, f64, f64, f64, i64 {
    let a = tv[3 * t];
    let b = tv[3 * t + 1];
    let c = tv[3 * t + 2];
    let bx = px[b] - px[a];
    let by = py[b] - py[a];
    let cx = px[c] - px[a];
    let cy = py[c] - py[a];
    let d = 2.0 * (bx * cy - by * cx);
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    let ux = (cy * b2 - by * c2) / d;
    let uy = (bx * c2 - cx * b2) / d;
    let r = sqrt(ux * ux + uy * uy);
    let ex = px[c] - px[b];
    let ey = py[c] - py[b];
    let e2 = ex * ex + ey * ey;
    let apex = if e2 <= b2 && e2 <= c2 { a } else if b2 <= c2 { c } else { b };
    return px[a] + ux, py[a] + uy, r, r / sqrt(min(min(b2, c2), e2)), apex;
}

// Meshes the 2D shape g with triangles of edge length about `size`. False
// if g holds a 3D shape, or if the refinement does not finish in its
// budget of points; the mesh is then the last triangulation.
fn mesh_2d(g: i64, size: f64) -> bool {
    clear_triangulation();
    guard collect_curves(g, 0.0, 0.0) else { return false; }
    for h in interfaces { guard collect_curves(h, 0.0, 0.0) else { return false; }; }
    let b = bounds(g);
    let w = b[3] - b[0];
    let h = b[4] - b[1];
    let diameter = sqrt(w * w + h * h);
    let cx = 0.5 * (b[0] + b[3]);
    let cy = 0.5 * (b[1] + b[4]);
    let far = 100.0 * diameter;
    px.push(cx - far);
    py.push(cy - far);
    px.push(cx + far);
    py.push(cy - far);
    px.push(cx);
    py.push(cy + far);
    for i in 3 { vtri.push(-1); }
    new_triangle(0, 1, 2);
    sample_boundary(g, size, diameter, false);
    if interfaces.len > 0 { sample_boundary(g, size, diameter, true); }
    find_corners();
    for v in 3..px.len { insert_existing(v, locate(px[v], py[v])); }
    // A budget of points: 20 times the points of a uniform mesh of the box,
    // or of the cells of the box with the refines.
    let uniform = (ceil(w / size) + 2.0) * (ceil(h / size) + 2.0);
    let estimate = if refines.len > 0 { max(uniform, cell_estimate(b, size)) } else { uniform };
    let budget = 20 * (estimate as! i64) + 100 * px.len;
    var done = recover_segments(budget);
    let ratio = 1.0 / (2.0 * sin(radians(MIN_ANGLE)));
    while done {
        classify();
        var inserted = false;
        for t in talive.len {
            if !talive[t] || !tin[t] { continue; }
            let x, y, r, q, apex = circumcentre(t);
            let gx = (px[tv[3 * t]] + px[tv[3 * t + 1]] + px[tv[3 * t + 2]]) / 3.0;
            let gy = (py[tv[3 * t]] + py[tv[3 * t + 1]] + py[tv[3 * t + 2]]) / 3.0;
            let wanted = if targets.len > 0 { min(local_size(size, gx, gy), wanted_size(x, y, section_z, r, size)) } else { local_size(size, gx, gy) };
            let largest = wanted / sqrt(3.0) * (1.0 + 1.0e-9);
            if (q <= ratio || sharp(apex, MIN_ANGLE)) && r <= largest { continue; }
            if px.len >= budget { done = false; break; }
            var split = false;
            for s in segments.len {
                if segments[s].alive && encroaches(s, x, y) {
                    split_segment(s);
                    split = true;
                }
            }
            if !split { insert_point(x, y); }
            inserted = true;
            if !recover_segments(budget) { done = false; break; }
        }
        if !inserted { break; }
    }
    add_triangles();
    return done;
}

fn clear_triangulation() {
    clear_mesh(2);
    curves.clear();
    px.clear();
    py.clear();
    vtri.clear();
    corner.clear();
    tv.clear();
    tn.clear();
    talive.clear();
    tfree.clear();
    tmark.clear();
    tin.clear();
    segments.clear();
    fan.clear();
}

// The mesh: the triangles inside g, with their points renumbered.
fn add_triangles() {
    classify();
    var index: i64[>..] = [];
    for v in px.len { index.push(-1); }
    for t in talive.len {
        if !talive[t] || !tin[t] { continue; }
        var nodes: i64[3] = [0, 0, 0];
        for k in 3 {
            let v = tv[3 * t + k];
            if index[v] < 0 { index[v] = add_point(px[v], py[v], 0.0); }
            nodes[k] = index[v];
        }
        add_simplex(nodes);
    }
}
