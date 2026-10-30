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

// `a` and `b` are child nodes, or the first point and the point count of a
// polygon. `p` holds the numbers of the primitive.
struct Shape {
    kind: i64,
    a: i64,
    b: i64,
    p: f64[8],
}

var shapes: Shape[>..] = [];
var polygon_points: f64[>..] = [];    // x and y of each polygon point

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

// The 2D shape s, from z = 0 to z = h.
fn extrude(s: i64, h: f64) -> i64 {
    return add_shape(EXTRUDE, s, -1, [h, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
}

// The 2D shape s (x is the radius) turned about the y axis by `angle`
// radians, from the +x axis towards -z.
fn revolve(s: i64, angle: f64) -> i64 {
    return add_shape(REVOLVE, s, -1, [angle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
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

// The angle wedge of a revolve: 0 <= atan2(-z, x) <= angle.
fn wedge_phi(angle: f64, x: f64, z: f64) -> f64 {
    let low = z;    // the half-space atan2(-z, x) >= 0
    let high = cos(angle) * -z - sin(angle) * x;
    return if angle <= PI { max(low, high) } else { min(low, high) };
}

recursive fn phi(g: i64, x: f64, y: f64, z: f64) -> f64 {
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
    if s.kind == EXTRUDE { return slab(phi(s.a, x, y, 0.0), z, p[0]); }
    if s.kind == REVOLVE {
        let d = phi(s.a, sqrt(x * x + z * z), y, 0.0);
        if p[0] >= TAU - 1.0e-12 { return d; }
        return max(d, wedge_phi(p[0], x, z));
    }
    let fa = phi(s.a, x, y, z);
    let fb = phi(s.b, x, y, z);
    if s.kind == UNION { return min(fa, fb); }
    if s.kind == INTERSECT { return max(fa, fb); }
    return max(fa, -fb);
}

// The bounding box of node g: x0, y0, z0, x1, y1, z1. A 2D shape has
// z0 = z1 = 0.
recursive fn bounds(g: i64) -> f64[6] {
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
    if s.kind == EXTRUDE {
        let a = bounds(s.a);
        return [a[0], a[1], 0.0, a[3], a[4], p[0]];
    }
    if s.kind == REVOLVE {
        let a = bounds(s.a);
        let r = max(abs(a[0]), abs(a[3]));
        return [-r, a[1], -r, r, a[4], r];
    }
    let a = bounds(s.a);
    if s.kind == SUBTRACT { return a; }
    let b = bounds(s.b);
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
