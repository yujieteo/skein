// The point cloud: the points, the triangles between them, and the
// boundary edges.
//
// The triangles are not cells of the method. They give the nodal domain of
// each point for the integration of rkpm (one third of each triangle of
// the point, Chen et al. 2001), the local spacing h of each point (the
// largest middle edge of its triangles), the boundary edges for the
// conditions, and the cell of a point for value_at. Point p is at (xy[2p], xy[2p + 1]).
import std;
import math;
namespace meshless;

let FREE = -1;       // the condition of a boundary edge that no fix holds

var xy: f64[>..] = [];
var tri: i64[>..] = [];         // 3 points for each triangle, counterclockwise
var spacing: f64[>..] = [];     // h of each point
var edges: i64[>..] = [];       // 2 points for each boundary edge, with the body on the left
var edge_tri: i64[>..] = [];    // the triangle of each boundary edge
var edge_fix: i64[>..] = [];    // the condition k of each component of each edge (ncomp each), or FREE
var edge_load: f64[>..] = [];   // the traction or the outward flux of each component of each edge
var selected: u8[>..] = [];     // select_faces holds the edge
var point_tris: i64[>..] = [];  // the triangles of each point, from point_start
var point_start: i64[>..] = [];
var lo: f64[2] = [0.0, 0.0];    // the bounding box
var hi: f64[2] = [0.0, 0.0];

fn point_count() -> i64 { return xy.len / 2; }
fn triangle_count() -> i64 { return tri.len / 3; }
fn edge_count() -> i64 { return edges.len / 2; }

fn norm_2(x: f64, y: f64) -> f64 { return sqrt(x * x + y * y); }

fn px(p: i64) -> f64 { return xy[2 * p]; }
fn py(p: i64) -> f64 { return xy[2 * p + 1]; }

// The distance below which two coordinates are the same point.
fn tolerance() -> f64 { return 1.0e-9 * max(hi[0] - lo[0], hi[1] - lo[1]); }

// Twice the signed area of triangle (a, b, c).
fn twice_area(a: i64, b: i64, c: i64) -> f64 {
    return (px(b) - px(a)) * (py(c) - py(a)) - (py(b) - py(a)) * (px(c) - px(a));
}

// Takes the points `points` (x and y of each) and the triangles
// `triangles` (3 points each, in either orientation) as the cloud.
fn set_cloud(points: f64[:], triangles: i64[:]) {
    guard points.len >= 6 && points.len % 2 == 0 && triangles.len >= 3 && triangles.len % 3 == 0 &&
          min(triangles) >= 0 && max(triangles) < points.len / 2 else { abort("set_cloud: 3 points or more, and triangles of 3 of them"); }
    xy.clear();
    xy.append(points);
    tri.clear();
    for t in triangles.len / 3 {
        let a = triangles[3 * t];
        let b = triangles[3 * t + 1];
        let c = triangles[3 * t + 2];
        tri.append(if twice_area(a, b, c) > 0.0 { [a, b, c] } else { [a, c, b] });
    }
    for j in 2 {
        lo[j] = points[j];
        hi[j] = points[j];
        for p in point_count() {
            lo[j] = min(lo[j], points[2 * p + j]);
            hi[j] = max(hi[j], points[2 * p + j]);
        }
    }
    point_spacing();
    boundary_edges();
    start_physics();
}

// The spacing of each point, and the triangles of each point.
fn point_spacing() {
    spacing.clear();
    point_start.clear();
    for p in point_count() {
        spacing.push(0.0);
        point_start.push(0);
    }
    point_start.push(0);
    for t in triangle_count() {
        var d: f64[3] = [0.0, 0.0, 0.0];
        for k in 3 {
            let a = tri[3 * t + k];
            let b = tri[3 * t + (k + 1) % 3];
            d[k] = norm_2(px(b) - px(a), py(b) - py(a));
        }
        // The middle edge: the spacing of a lattice cut into triangles, and
        // the larger spacing of a stretched triangle.
        let middle = d[0] + d[1] + d[2] - min(d[0], min(d[1], d[2])) - max(d[0], max(d[1], d[2]));
        for k in 3 {
            let a = tri[3 * t + k];
            spacing[a] = max(spacing[a], middle);
            point_start[a + 1]++;
        }
    }
    for p in point_count() { point_start[p + 1] += point_start[p]; }
    point_tris.clear();
    for t in tri.len { point_tris.push(0); }
    var at: i64[>..] = [];
    at.append(point_start[..point_count()]);
    for t in triangle_count() {
        for k in 3 {
            let a = tri[3 * t + k];
            point_tris[at[a]] = t;
            at[a]++;
        }
    }
}

// The edges of one triangle, with the triangle on their left.
fn boundary_edges() {
    let n = point_count();
    var keys: i64[>..] = [];
    for t in triangle_count() {
        for k in 3 {
            let a = tri[3 * t + k];
            let b = tri[3 * t + (k + 1) % 3];
            keys.push(min(a, b) * n + max(a, b));
        }
    }
    sort(keys);
    edges.clear();
    edge_tri.clear();
    for t in triangle_count() {
        for k in 3 {
            let a = tri[3 * t + k];
            let b = tri[3 * t + (k + 1) % 3];
            let key = min(a, b) * n + max(a, b);
            let i = lower_bound(keys, key);
            let shared = (i > 0 && keys[i - 1] == key) || (i + 1 < keys.len && keys[i + 1] == key);
            if !shared {
                edges.append([a, b]);
                edge_tri.push(t);
            }
        }
    }
    selected.clear();
    for e in edge_count() { selected.push(0); }
}

// The first index of `key` in the sorted `keys`.
fn lower_bound(keys: i64[:], key: i64) -> i64 {
    var a = 0;
    var b = keys.len;
    while a < b {
        let m = (a + b) / 2;
        if keys[m] < key { a = m + 1; } else { b = m; }
    }
    return a;
}

// The outward normal of edge e, and its length.
fn edge_normal(e: i64) -> f64[2], f64 {
    let a = edges[2 * e];
    let b = edges[2 * e + 1];
    let dx = px(b) - px(a);
    let dy = py(b) - py(a);
    let l = norm_2(dx, dy);
    return [dy / l, -dx / l], l;
}

// The point on edge e at the parameter s from 0 to 1.
fn edge_point(e: i64, s: f64) -> f64[2] {
    let a = edges[2 * e];
    let b = edges[2 * e + 1];
    return [px(a) + s * (px(b) - px(a)), py(a) + s * (py(b) - py(a))];
}

// A lattice of (nx + 1) × (ny + 1) points on the rectangle from (xa, ya)
// to (xb, yb), with each point inside moved by up to JITTER of the
// spacing on each axis, so that the cloud is not a grid. The move of each
// point is the same on each run.
let JITTER = 0.15;

fn points_rect(xa: f64, ya: f64, xb: f64, yb: f64, nx: i64, ny: i64) {
    guard xb > xa && yb > ya && nx >= 1 && ny >= 1 else { abort("points_rect: xb > xa, yb > ya, and nx, ny >= 1"); }
    let hx = (xb - xa) / nx;
    let hy = (yb - ya) / ny;
    var points: f64[>..] = [];
    for j in ny + 1 {
        for i in nx + 1 {
            let inside = i > 0 && i < nx && j > 0 && j < ny;
            let p = i + (nx + 1) * j;
            points.push(xa + i * hx + (if inside { JITTER * hx * scatter(2 * p) } else { 0.0 }));
            points.push(ya + j * hy + (if inside { JITTER * hy * scatter(2 * p + 1) } else { 0.0 }));
        }
    }
    var triangles: i64[>..] = [];
    for j in ny {
        for i in nx {
            let a = i + (nx + 1) * j;
            triangles.append([a, a + 1, a + nx + 2, a, a + nx + 2, a + nx + 1]);
        }
    }
    set_cloud(points, triangles);
    // The spacing of the lattice, not of the moved points. Reason: with the
    // spacing of each moved point, the rkpm order of
    // examples/meshless/mms_rkpm.skn from 32 × 32 to 64 × 64 was 1.73;
    // with the lattice spacing it is 2.02.
    for &h in spacing { h = max(hx, hy); }
}

// A number from -1 to 1 for each k, from an integer hash.
fn scatter(k: i64) -> f64 {
    guard k >= 0 && k < 1 << 30 else { abort("scatter: 0 <= k < 2^30"); }
    var h = ((k + 1) * 0x9E3779B1) & 0x7FFFFFF;
    h = ((h ^ (h >> 13)) * 0x85EBCA77) & 0x7FFFFFF;
    h = (h ^ (h >> 11)) & 0xFFFFFF;
    return 2.0 * h / 16777215.0 - 1.0;
}
