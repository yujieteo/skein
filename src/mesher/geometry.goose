// The geometry of spec §8.1: a CSG tree of primitives, with the signed
// distance function φ of each node (φ < 0 inside, φ = 0 on the boundary,
// φ > 0 outside).
//
// A node is an index in `shapes`. The functions that make a node return
// its index. A 2D shape lies in the x-y plane and ignores z.
import std;
import math;
namespace mesher;

let BOX = 0;
let CYLINDER = 1;
let RECT = 2;
let CIRCLE = 3;
let POLYGON = 4;
let EXTRUDE = 5;
let REVOLVE = 6;
let UNION = 7;
let INTERSECT = 8;
let SUBTRACT = 9;
let LOFT = 10;
let PATTERN = 11;
let PATTERN_POLAR = 12;
let FILL = 13;

// `a` and `b` are child nodes, the first point and the point count of a
// polygon or a fill, or the child node and the copy count of a pattern. `p`
// holds the numbers of the primitive. A fill of a nurbs has p[0] = 1: its
// span i is the 12 numbers from polygon_points[2 a + 12 i], and b is the
// span count + 1.
struct Shape {
    kind: i64,
    a: i64,
    b: i64,
    p: f64[8],
}

var shapes: Shape[>..] = [];
var polygon_points: f64[>..] = [];    // x and y of each point of a polygon or a fill

fn add_shape(kind: i64, a: i64, b: i64, p: f64[8]) -> i64 {
    shapes.push(Shape { kind: kind, a: a, b: b, p: p });
    return shapes.len - 1;
}

// Removes all shapes.
fn clear_shapes() {
    shapes.clear();
    polygon_points.clear();
}

fn box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> i64 {
    return add_shape(BOX, -1, -1,
        [min(x0, x1), min(y0, y1), min(z0, z1), max(x0, x1), max(y0, y1), max(z0, z1), 0.0, 0.0]);
}

// Base centre (cx, cy, cz), axis direction (ax, ay, az), radius r and
// height h.
fn cylinder(cx: f64, cy: f64, cz: f64, ax: f64, ay: f64, az: f64, r: f64, h: f64) -> i64 {
    let l = sqrt(ax * ax + ay * ay + az * az);
    return add_shape(CYLINDER, -1, -1, [cx, cy, cz, ax / l, ay / l, az / l, r, h]);
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> i64 {
    return add_shape(RECT, -1, -1, [min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1), 0.0, 0.0, 0.0, 0.0]);
}

fn circle(cx: f64, cy: f64, r: f64) -> i64 {
    return add_shape(CIRCLE, -1, -1, [cx, cy, r, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// `xy` holds x and y of each point, counter-clockwise or clockwise.
fn polygon(xy: f64[:]) -> i64 {
    let first = polygon_points.len / 2;
    for v in xy { polygon_points.push(v); }
    return add_shape(POLYGON, first, xy.len / 2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// The region inside the spline through the points of `xy` and the line
// from the last point to the first (spec §8.1). Each point is apart from
// the next.
fn fill(xy: f64[:]) -> i64 {
    let first = polygon_points.len / 2;
    for v in xy { polygon_points.push(v); }
    return add_shape(FILL, first, xy.len / 2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// The region inside the NURBS curve of degree p (1 to 3) with the knots
// `knots` and the points x, y, w of `xyw`, and the line from its end to
// its start (spec §8.1). The knots do not fall, the first and the last
// are each p + 1 times, and an inner knot is p times or fewer. Knot insertion (Boehm) makes each inner knot
// p times: then each knot span is a rational Bézier curve of degree p,
// which is raised to degree 3.
fn fill_nurbs(p: i64, knots: f64[:], xyw: f64[:]) -> i64 {
    guard p >= 1 && p <= 3 && knots.len == xyw.len / 3 + p + 1 && xyw.len / 3 >= p + 1 else {
        abort("fill_nurbs: 1 <= p <= 3, p + 1 points or more, and a knot for each point and p + 1 more");
    }
    var k: f64[>..] = [];
    k.append(knots);
    // The homogeneous points: x w, y w, w.
    var h: f64[>..] = [];
    for i in xyw.len / 3 {
        let w = xyw[3 * i + 2];
        h.push(xyw[3 * i] * w);
        h.push(xyw[3 * i + 1] * w);
        h.push(w);
    }
    var i = p + 1;
    while i < k.len - p - 1 {
        var m = 1;
        while i + m < k.len - p - 1 && k[i + m] == k[i] { m++; }
        let u = k[i];
        for r in p - m { insert_knot(p, k, h, u); }
        i += max(m, p);
    }
    let first = polygon_points.len / 2;
    var spans = 0;
    for j in (h.len / 3 - 1) / p {
        let c = bezier_span(p, h[3 * j * p..]);
        for v in c { polygon_points.push(v); }
        spans++;
    }
    return add_shape(FILL, first, spans + 1, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// Inserts the knot u once into the curve of degree p with the knots k and
// the homogeneous points h (Boehm). k[p] <= u < k[k.len - p - 1].
fn insert_knot(p: i64, k: f64[>..]&, h: f64[>..]&, u: f64) {
    var s = p;
    while s + 1 < k.len - p - 1 && k[s + 1] <= u { s++; }
    var q: f64[>..] = [];
    for i in h.len / 3 + 1 {
        for c in 3 {
            if i <= s - p {
                q.push(h[3 * i + c]);
            } else if i > s {
                q.push(h[3 * (i - 1) + c]);
            } else {
                let a = (u - k[i]) / (k[i + p] - k[i]);
                q.push(a * h[3 * i + c] + (1.0 - a) * h[3 * (i - 1) + c]);
            }
        }
    }
    var kk: f64[>..] = [];
    kk.append(k[..s + 1]);
    kk.push(u);
    kk.append(k[s + 1..]);
    k = copy(kk);
    h = copy(q);
}

// The span of the rational Bézier curve of degree p with the homogeneous
// points h[0..3 (p + 1)], raised to degree 3, as the cubics of x w, y w
// and w (span).
fn bezier_span(p: i64, h: f64[:]) -> f64[12] {
    var c: f64[12] = default<f64[12]>();
    for d in 3 {
        fn at(i: i64) -> f64 { h[3 * i + d] }
        let b: f64[4] = if p == 1 { [at(0), (2.0 * at(0) + at(1)) / 3.0, (at(0) + 2.0 * at(1)) / 3.0, at(1)] }
                        else if p == 2 { [at(0), (at(0) + 2.0 * at(1)) / 3.0, (2.0 * at(1) + at(2)) / 3.0, at(2)] }
                        else { [at(0), at(1), at(2), at(3)] };
        c[4 * d] = b[0];
        c[4 * d + 1] = 3.0 * (b[1] - b[0]);
        c[4 * d + 2] = 3.0 * (b[2] - 2.0 * b[1] + b[0]);
        c[4 * d + 3] = b[3] - 3.0 * b[2] + 3.0 * b[1] - b[0];
    }
    return c;
}

// The 2D shape s, from z = 0 to z = h.
fn extrude(s: i64, h: f64) -> i64 {
    return add_shape(EXTRUDE, s, -1, [h, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// The 2D shape s (x is the radius) turned about the y axis by `angle`
// radians, from the +x axis towards -z.
fn revolve(s: i64, angle: f64) -> i64 {
    return add_shape(REVOLVE, s, -1, [angle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// From the 2D section s1 at z = 0 to the 2D section s2 at z = h. At a
// height z, φ is (1 - z/h) φ1 + (z/h) φ2.
fn loft(s1: i64, s2: i64, h: f64) -> i64 {
    return add_shape(LOFT, s1, s2, [h, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// n copies of shape a: copy k (0 <= k < n) is a moved by k (vx, vy, vz).
fn pattern(a: i64, vx: f64, vy: f64, vz: f64, n: i64) -> i64 {
    return add_shape(PATTERN, a, n, [vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// n copies of shape a: copy k (0 <= k < n) is a turned by k 2π/n about the
// axis (ax, ay, az) through the origin.
fn pattern_polar(a: i64, ax: f64, ay: f64, az: f64, n: i64) -> i64 {
    let l = sqrt(ax * ax + ay * ay + az * az);
    return add_shape(PATTERN_POLAR, a, n, [ax / l, ay / l, az / l, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// The point q turned by `angle` about the unit axis u through the origin
// (Rodrigues).
fn turn(q: f64[3], u: f64[3], angle: f64) -> f64[3] {
    let c = cos(angle);
    let s = sin(angle);
    let d = (u[0] * q[0] + u[1] * q[1] + u[2] * q[2]) * (1.0 - c);
    return [q[0] * c + (u[1] * q[2] - u[2] * q[1]) * s + u[0] * d,
            q[1] * c + (u[2] * q[0] - u[0] * q[2]) * s + u[1] * d,
            q[2] * c + (u[0] * q[1] - u[1] * q[0]) * s + u[2] * d];
}

fn union(a: i64, b: i64) -> i64 { return add_shape(UNION, a, b, default<f64[8]>()); }
fn intersect(a: i64, b: i64) -> i64 { return add_shape(INTERSECT, a, b, default<f64[8]>()); }
fn subtract(a: i64, b: i64) -> i64 { return add_shape(SUBTRACT, a, b, default<f64[8]>()); }

// The distance function of a slab |t - h/2| <= h/2 cut by a 2D or radial
// distance d: exact outside the edges, and the larger distance inside.
fn slab(d: f64, t: f64, h: f64) -> f64 {
    let e = abs(t - 0.5 * h) - 0.5 * h;
    let ox = max(d, 0.0);
    let oy = max(e, 0.0);
    return min(max(d, e), 0.0) + sqrt(ox * ox + oy * oy);
}

fn polygon_phi(first: i64, count: i64, x: f64, y: f64) -> f64 {
    var d2 = -1.0;
    var inside = false;
    for i in count {
        let j = (i + 1) % count;
        let ax = polygon_points[2 * (first + i)];
        let ay = polygon_points[2 * (first + i) + 1];
        let bx = polygon_points[2 * (first + j)];
        let by = polygon_points[2 * (first + j) + 1];
        let ex = bx - ax;
        let ey = by - ay;
        let t = clamp(((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey), 0.0, 1.0);
        let dx = x - ax - t * ex;
        let dy = y - ay - t * ey;
        let d = dx * dx + dy * dy;
        if d2 < 0.0 || d < d2 { d2 = d; }
        if (ay > y) != (by > y) && x < ax + (y - ay) * ex / ey { inside = !inside; }
    }
    return if inside { -sqrt(d2) } else { sqrt(d2) };
}

// Point k of fill node g. The points before the first and after the last
// are reflections: 2 P(0) - P(1) and 2 P(n-1) - P(n-2).
fn fill_point(g: i64, k: i64) -> f64[2] {
    let s = shapes[g];
    if s.p[0] == 1.0 {
        // A nurbs: the start of span 0 or the end of the last span.
        let c = stored_span(g, if k <= 0 { 0 } else { s.b - 2 });
        let u = if k <= 0 { 0.0 } else { 1.0 };
        let x, _, _ = cubic(c, 0, u);
        let y, _, _ = cubic(c, 1, u);
        return [x, y];
    }
    fn at(j: i64) -> f64[2] { [polygon_points[2 * (s.a + j)], polygon_points[2 * (s.a + j) + 1]] }
    if k < 0 {
        let a = at(0);
        let b = at(1);
        return [2.0 * a[0] - b[0], 2.0 * a[1] - b[1]];
    }
    if k >= s.b {
        let a = at(s.b - 1);
        let b = at(s.b - 2);
        return [2.0 * a[0] - b[0], 2.0 * a[1] - b[1]];
    }
    return at(k);
}

// Span i of the nurbs of fill node g.
fn stored_span(g: i64, i: i64) -> f64[12] {
    var c: f64[12] = default<f64[12]>();
    for j in 12 { c[j] = polygon_points[2 * shapes[g].a + 12 * i + j]; }
    return c;
}

// Span i of fill node g, from point i to point i + 1, as the cubics
// X = c0 + c1 u + c2 u^2 + c3 u^3 (c[0..4]), Y (c[4..8]) and W (c[8..12])
// for u in [0, 1]: the point is (X / W, Y / W). The spline is the
// centripetal Catmull-Rom spline (Yuksel et al. 2011) in its Hermite
// form, with W = 1: it has no cusp and no loop in a span. A nurbs keeps
// its spans in polygon_points.
fn span(g: i64, i: i64) -> f64[12] {
    if shapes[g].p[0] == 1.0 { return stored_span(g, i); }
    let p0 = fill_point(g, i - 1);
    let p1 = fill_point(g, i);
    let p2 = fill_point(g, i + 1);
    let p3 = fill_point(g, i + 2);
    fn knot(a: f64[2], b: f64[2]) -> f64 { sqrt(sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]))) }
    let d0 = knot(p0, p1);
    let d1 = knot(p1, p2);
    let d2 = knot(p2, p3);
    var c: f64[12] = default<f64[12]>();
    c[8] = 1.0;
    for k in 2 {
        let m1 = d1 * ((p1[k] - p0[k]) / d0 - (p2[k] - p0[k]) / (d0 + d1) + (p2[k] - p1[k]) / d1);
        let m2 = d1 * ((p2[k] - p1[k]) / d1 - (p3[k] - p1[k]) / (d1 + d2) + (p3[k] - p2[k]) / d2);
        c[4 * k] = p1[k];
        c[4 * k + 1] = m1;
        c[4 * k + 2] = 3.0 * (p2[k] - p1[k]) - 2.0 * m1 - m2;
        c[4 * k + 3] = 2.0 * (p1[k] - p2[k]) + m1 + m2;
    }
    return c;
}

// The value of coordinate j (0 is x, 1 is y) of span c at u, and its
// first and second derivatives: v = N / W, by the quotient rule.
fn cubic(c: f64[12], j: i64, u: f64) -> f64, f64, f64 {
    fn poly(a: i64) -> f64, f64, f64 {
        return c[a] + u * (c[a + 1] + u * (c[a + 2] + u * c[a + 3])),
               c[a + 1] + u * (2.0 * c[a + 2] + u * 3.0 * c[a + 3]),
               2.0 * c[a + 2] + 6.0 * u * c[a + 3];
    }
    let n, n1, n2 = poly(4 * j);
    let w, w1, w2 = poly(8);
    let v = n / w;
    let v1 = (n1 - w1 * v) / w;
    return v, v1, (n2 - 2.0 * w1 * v1 - w2 * v) / w;
}

// Whether span c has W = 1: a span of a spline.
fn polynomial(c: f64[12]) -> bool { c[8] == 1.0 && c[9] == 0.0 && c[10] == 0.0 && c[11] == 0.0 }

// The box x0, y0, x1, y1 of the Bézier points of span c. It holds the
// span: the weights are above 0.
fn span_box(c: f64[12]) -> f64[4] {
    var b: f64[4] = [1.0e300, 1.0e300, -1.0e300, -1.0e300];
    fn bezier(a: i64) -> f64[4] {
        [c[a], c[a] + c[a + 1] / 3.0, c[a] + (2.0 * c[a + 1] + c[a + 2]) / 3.0, c[a] + c[a + 1] + c[a + 2] + c[a + 3]]
    }
    let w = bezier(8);
    for j in 2 {
        let e = bezier(4 * j);
        for v, i in e {
            b[j] = min(b[j], v / w[i]);
            b[2 + j] = max(b[2 + j], v / w[i]);
        }
    }
    return b;
}

// The number of times span c crosses the ray from (x, y) along +x, with
// the half-open rule of polygon_phi. Each piece between 2 zeros of dy/du
// is monotone in y, so it crosses at most once.
fn span_crossings(c: f64[12], x: f64, y: f64) -> i64 {
    var cuts: f64[6] = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    var count = 1;
    if !polynomial(c) {
        // The zeros of dy/du of a rational span: its numerator has degree
        // 4 or less, so 64 samples and a bisection at each change of sign
        // find them, apart from a double zero, where y does not turn.
        var lo = 0.0;
        let _, d0, _ = cubic(c, 1, 0.0);
        var dlo = d0;
        for m in 64 {
            let hi = (m + 1) / 64.0;
            let _, dhi, _ = cubic(c, 1, hi);
            if (dlo > 0.0) != (dhi > 0.0) && count < 5 && hi < 1.0 {
                var a = lo;
                var b = hi;
                for r in 60 {
                    let mid = 0.5 * (a + b);
                    let _, dm, _ = cubic(c, 1, mid);
                    if (dm > 0.0) == (dlo > 0.0) { a = mid; } else { b = mid; }
                }
                cuts[count] = 0.5 * (a + b);
                count++;
            }
            lo = hi;
            dlo = dhi;
        }
        cuts[count] = 1.0;
        return piece_crossings(c, cuts, count, x, y);
    }
    // The zeros of c5 + 2 c6 u + 3 c7 u^2 in (0, 1).
    let qa = 3.0 * c[7];
    let qb = 2.0 * c[6];
    let disc = qb * qb - 4.0 * qa * c[5];
    if abs(qa) > 1.0e-300 && disc > 0.0 {
        let r0 = (-qb - sqrt(disc)) / (2.0 * qa);
        let r1 = (-qb + sqrt(disc)) / (2.0 * qa);
        for r in [min(r0, r1), max(r0, r1)] {
            if r > 0.0 && r < 1.0 { cuts[count] = r; count++; }
        }
    } else if abs(qa) <= 1.0e-300 && abs(qb) > 1.0e-300 {
        let r = -c[5] / qb;
        if r > 0.0 && r < 1.0 { cuts[count] = r; count++; }
    }
    cuts[count] = 1.0;
    return piece_crossings(c, cuts, count, x, y);
}

// The crossings of the ray from (x, y) along +x with the `count` pieces
// of span c between cuts[k] and cuts[k + 1]. y is monotone in each piece.
fn piece_crossings(c: f64[12], cuts: f64[6], count: i64, x: f64, y: f64) -> i64 {
    var n = 0;
    for k in count {
        var lo = cuts[k];
        var hi = cuts[k + 1];
        let ylo, _, _ = cubic(c, 1, lo);
        let yhi, _, _ = cubic(c, 1, hi);
        guard (ylo > y) != (yhi > y) else { continue; }
        for m in 60 {
            let mid = 0.5 * (lo + hi);
            let ym, _, _ = cubic(c, 1, mid);
            if (ym > y) == (ylo > y) { lo = mid; } else { hi = mid; }
        }
        let xr, _, _ = cubic(c, 0, 0.5 * (lo + hi));
        if x < xr { n++; }
    }
    return n;
}

// The square of the distance from (x, y) to span c: the nearest of 9
// samples, then Newton steps on (P(u) - q) . P'(u) = 0.
fn span_distance2(c: f64[12], x: f64, y: f64) -> f64 {
    fn d2(u: f64) -> f64 {
        let px, _, _ = cubic(c, 0, u);
        let py, _, _ = cubic(c, 1, u);
        return (px - x) * (px - x) + (py - y) * (py - y);
    }
    var best = 1.0e300;
    var u = 0.0;
    for k in 9 {
        let d = d2(k / 8.0);
        if d < best { best = d; u = k / 8.0; }
    }
    for k in 6 {
        let px, dx, ddx = cubic(c, 0, u);
        let py, dy, ddy = cubic(c, 1, u);
        let f = (px - x) * dx + (py - y) * dy;
        let df = dx * dx + dy * dy + (px - x) * ddx + (py - y) * ddy;
        guard df > 0.0 else { break; }
        u = clamp(u - f / df, 0.0, 1.0);
        best = min(best, d2(u));
    }
    return best;
}

// φ of fill node g: the distance to the nearest span or to the closing
// line, with the sign of the crossing count of the ray along +x.
fn fill_phi(g: i64, x: f64, y: f64) -> f64 {
    let s = shapes[g];
    let a = fill_point(g, s.b - 1);
    let b = fill_point(g, 0);
    let ex = b[0] - a[0];
    let ey = b[1] - a[1];
    let t = clamp(((x - a[0]) * ex + (y - a[1]) * ey) / max(ex * ex + ey * ey, 1.0e-300), 0.0, 1.0);
    var d2 = (x - a[0] - t * ex) * (x - a[0] - t * ex) + (y - a[1] - t * ey) * (y - a[1] - t * ey);
    var inside = (a[1] > y) != (b[1] > y) && x < a[0] + (y - a[1]) * ex / ey;
    for i in s.b - 1 {
        let c = span(g, i);
        let box = span_box(c);
        if y >= box[1] && y <= box[3] && x <= box[2] && span_crossings(c, x, y) % 2 == 1 { inside = !inside; }
        let ox = max(max(box[0] - x, x - box[2]), 0.0);
        let oy = max(max(box[1] - y, y - box[3]), 0.0);
        if ox * ox + oy * oy < d2 { d2 = min(d2, span_distance2(c, x, y)); }
    }
    return if inside { -sqrt(d2) } else { sqrt(d2) };
}

// The angle wedge of a revolve: 0 <= atan2(-z, x) <= angle.
fn wedge_phi(angle: f64, x: f64, z: f64) -> f64 {
    let low = z;    // the half-space atan2(-z, x) >= 0
    let high = cos(angle) * -z - sin(angle) * x;
    return if angle <= PI { max(low, high) } else { min(low, high) };
}

// A step of the walk of phi: shape g at the point q, the next child k,
// and the value of the children before k.
struct Frame {
    g: i64,
    q: f64[3],
    k: i64,
    acc: f64,
}

var frames: Frame[>..] = [];

fn is_leaf(kind: i64) -> bool {
    return kind == BOX || kind == RECT || kind == CYLINDER || kind == CIRCLE || kind == POLYGON || kind == FILL;
}

// The signed distance from (x, y, z) to shape g with no children.
fn leaf_phi(g: i64, x: f64, y: f64, z: f64) -> f64 {
    let s = shapes[g];
    let p = s.p;
    if s.kind == BOX || s.kind == RECT {
        let n = if s.kind == BOX { 3 } else { 2 };
        let q: f64[3] = [x, y, z];
        var outside = 0.0;
        var inside = -1.0e300;
        for i in n {
            let d = max(p[i] - q[i], q[i] - p[n + i]);
            outside += max(d, 0.0) * max(d, 0.0);
            inside = max(inside, d);
        }
        return min(inside, 0.0) + sqrt(outside);
    }
    if s.kind == CYLINDER {
        let vx = x - p[0];
        let vy = y - p[1];
        let vz = z - p[2];
        let t = vx * p[3] + vy * p[4] + vz * p[5];
        let rx = vx - t * p[3];
        let ry = vy - t * p[4];
        let rz = vz - t * p[5];
        return slab(sqrt(rx * rx + ry * ry + rz * rz) - p[6], t, p[7]);
    }
    if s.kind == CIRCLE {
        return sqrt((x - p[0]) * (x - p[0]) + (y - p[1]) * (y - p[1])) - p[2];
    }
    if s.kind == POLYGON { return polygon_phi(s.a, s.b, x, y); }
    if s.kind == FILL { return fill_phi(g, x, y); }
    return fill_phi(g, x, y);
}

fn child_count(s: Shape) -> i64 {
    if s.kind == EXTRUDE || s.kind == REVOLVE { return 1; }
    if s.kind == PATTERN || s.kind == PATTERN_POLAR { return s.b; }
    return 2;
}

// Child k of shape s, and the point q in the frame of that child.
fn child_at(s: Shape, k: i64, q: f64[3]) -> i64, f64[3] {
    let p = s.p;
    if s.kind == EXTRUDE || s.kind == LOFT { return if k == 0 { s.a } else { s.b }, [q[0], q[1], 0.0]; }
    if s.kind == REVOLVE { return s.a, [sqrt(q[0] * q[0] + q[2] * q[2]), q[1], 0.0]; }
    if s.kind == PATTERN { return s.a, [q[0] - k * p[0], q[1] - k * p[1], q[2] - k * p[2]]; }
    if s.kind == PATTERN_POLAR { return s.a, turn(q, [p[0], p[1], p[2]], -TAU * k / s.b); }
    return if k == 0 { s.a } else { s.b }, q;
}

// Adds the value v of child k of shape s at q to acc, the value of the
// children before k.
fn combine(s: Shape, k: i64, q: f64[3], acc: f64, v: f64) -> f64 {
    if s.kind == LOFT {
        let t = clamp(q[2] / s.p[0], 0.0, 1.0);
        return if k == 0 { (1.0 - t) * v } else { acc + t * v };
    }
    if s.kind == PATTERN || s.kind == PATTERN_POLAR { return min(acc, v); }
    if k == 0 { return v; }
    if s.kind == UNION { return min(acc, v); }
    if s.kind == INTERSECT { return max(acc, v); }
    return max(acc, -v);
}

// The value of shape s at q from acc, the value of its children.
fn finish(s: Shape, q: f64[3], acc: f64) -> f64 {
    if s.kind == EXTRUDE || s.kind == LOFT { return slab(acc, q[2], s.p[0]); }
    if s.kind == REVOLVE && s.p[0] < TAU - 1.0e-12 { return max(acc, wedge_phi(s.p[0], q[0], q[2])); }
    return acc;
}

// The signed distance from (x, y, z) to shape g, < 0 inside. A stack of
// frames walks the tree of g, with no recursion (rule 1 of
// docs/tradeoffs.md §6). Each child has a smaller index than its shape, so
// the stack holds at most g + 1 frames.
fn phi(g: i64, x: f64, y: f64, z: f64) -> f64 {
    if is_leaf(shapes[g].kind) { return leaf_phi(g, x, y, z); }
    frames.clear();
    frames.push(Frame { g: g, q: [x, y, z], k: 0, acc: 1.0e300 });
    var v = 0.0;    // the value of the last frame that finished
    while frames.len > 0 {
        let f = frames.len - 1;
        let s = shapes[frames[f].g];
        let q = frames[f].q;
        if is_leaf(s.kind) {
            v = leaf_phi(frames[f].g, q[0], q[1], q[2]);
            frames.pop();
            continue;
        }
        if frames[f].k > 0 { frames[f].acc = combine(s, frames[f].k - 1, q, frames[f].acc, v); }
        if frames[f].k == child_count(s) {
            v = finish(s, q, frames[f].acc);
            frames.pop();
            continue;
        }
        let c, cq = child_at(s, frames[f].k, q);
        frames[f].k++;
        frames.push(Frame { g: c, q: cq, k: 0, acc: 1.0e300 });
    }
    return v;
}

// The bounding box of shape g with no children.
fn leaf_bounds(g: i64) -> f64[6] {
    let s = shapes[g];
    let p = s.p;
    if s.kind == BOX { return [p[0], p[1], p[2], p[3], p[4], p[5]]; }
    if s.kind == RECT { return [p[0], p[1], 0.0, p[2], p[3], 0.0]; }
    if s.kind == CIRCLE { return [p[0] - p[2], p[1] - p[2], 0.0, p[0] + p[2], p[1] + p[2], 0.0]; }
    if s.kind == POLYGON {
        var b: f64[6] = [1.0e300, 1.0e300, 0.0, -1.0e300, -1.0e300, 0.0];
        for i in s.b {
            for c in 2 {
                let v = polygon_points[2 * (s.a + i) + c];
                b[c] = min(b[c], v);
                b[3 + c] = max(b[3 + c], v);
            }
        }
        return b;
    }
    if s.kind == FILL {
        var b: f64[6] = [1.0e300, 1.0e300, 0.0, -1.0e300, -1.0e300, 0.0];
        for i in s.b - 1 {
            let e = span_box(span(g, i));
            for c in 2 {
                b[c] = min(b[c], e[c]);
                b[3 + c] = max(b[3 + c], e[2 + c]);
            }
        }
        return b;
    }
    if s.kind == CYLINDER {
        // The two end discs, each in the box of its centre ± r sin(angle to
        // the axis) along each axis.
        var b: f64[6] = [1.0e300, 1.0e300, 1.0e300, -1.0e300, -1.0e300, -1.0e300];
        for c in 3 {
            let e = p[6] * sqrt(max(1.0 - p[3 + c] * p[3 + c], 0.0));
            let end = p[c] + p[7] * p[3 + c];
            b[c] = min(p[c], end) - e;
            b[3 + c] = max(p[c], end) + e;
        }
        return b;
    }
    abort("leaf_bounds: a shape with no children");
}

// The bounding box of shape g from a and b, the boxes of its children.
fn own_bounds(g: i64, a: f64[6], b: f64[6]) -> f64[6] {
    let s = shapes[g];
    let p = s.p;
    if s.kind == EXTRUDE {
        return [a[0], a[1], 0.0, a[3], a[4], p[0]];
    }
    if s.kind == REVOLVE {
        let r = max(abs(a[0]), abs(a[3]));
        return [-r, a[1], -r, r, a[4], r];
    }
    if s.kind == LOFT {
        return [min(a[0], b[0]), min(a[1], b[1]), 0.0, max(a[3], b[3]), max(a[4], b[4]), p[0]];
    }
    if s.kind == PATTERN {
        var u: f64[6] = a;
        for c in 3 {
            let shift = (s.b - 1) * p[c];
            u[c] = min(a[c], a[c] + shift);
            u[3 + c] = max(a[3 + c], a[3 + c] + shift);
        }
        return u;
    }
    if s.kind == PATTERN_POLAR {
        // The box of the 8 turned corners of the box of a, for each copy.
        var u: f64[6] = [1.0e300, 1.0e300, 1.0e300, -1.0e300, -1.0e300, -1.0e300];
        for k in s.b {
            for m in 8 {
                let q = turn([a[3 * (m % 2)], a[1 + 3 * (m / 2 % 2)], a[2 + 3 * (m / 4)]],
                             [p[0], p[1], p[2]], TAU * k / s.b);
                for c in 3 {
                    u[c] = min(u[c], q[c]);
                    u[3 + c] = max(u[3 + c], q[c]);
                }
            }
        }
        return u;
    }
    if s.kind == SUBTRACT { return a; }
    var u: f64[6] = a;
    for c in 3 {
        if s.kind == UNION {
            u[c] = min(a[c], b[c]);
            u[3 + c] = max(a[3 + c], b[3 + c]);
        } else {
            u[c] = max(a[c], b[c]);
            u[3 + c] = min(a[3 + c], b[3 + c]);
        }
    }
    return u;
}

var box_scratch: f64[>..] = [];    // scratch of bounds: 6 values for each shape

// The bounding box of shape g: x0, y0, z0, x1, y1, z1. A 2D shape has
// z0 = z1 = 0. Each child has a smaller index than its shape, so one pass
// over the shapes up to g gives the box of each child first.
fn bounds(g: i64) -> f64[6] {
    box_scratch.clear();
    for i in g + 1 {
        let s = shapes[i];
        var box: f64[6] = [0.0; 6];
        if is_leaf(s.kind) {
            box = leaf_bounds(i);
        } else {
            let two = s.kind == LOFT || s.kind == UNION || s.kind == INTERSECT;
            box = own_bounds(i, box_of(s.a), if two { box_of(s.b) } else { box_of(s.a) });
        }
        box_scratch.append(box);
    }
    return box_of(g);
}

fn box_of(i: i64) -> f64[6] {
    var box: f64[6] = [0.0; 6];
    for c in 6 { box[c] = box_scratch[6 * i + c]; }
    return box;
}

// The point x moved onto the boundary φ = 0 of g: Newton steps along the
// gradient of φ, from central differences.
fn project(g: i64, x0: f64[3]) -> f64[3] {
    let b = bounds(g);
    var d = 0.0;
    for c in 3 { d = max(d, b[3 + c] - b[c]); }
    let h = 1.0e-7 * d;
    var x = x0;
    for step in 8 {
        let f = phi(g, x[0], x[1], x[2]);
        guard abs(f) > 1.0e-12 * d else { break; };
        var grad: f64[3] = [0.0, 0.0, 0.0];
        var g2 = 0.0;
        for c in 3 {
            var lo = x;
            var hi = x;
            lo[c] -= h;
            hi[c] += h;
            grad[c] = (phi(g, hi[0], hi[1], hi[2]) - phi(g, lo[0], lo[1], lo[2])) / (2.0 * h);
            g2 += grad[c] * grad[c];
        }
        guard g2 > 0.0 else { break; };
        for c in 3 { x[c] -= f * grad[c] / g2; }
    }
    return x;
}
