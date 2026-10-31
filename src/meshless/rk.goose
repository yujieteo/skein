// The reproducing kernel shape functions of rkpm with the linear basis
// (Liu, Jun and Zhang 1995, Chen et al. 1996):
//
//     Ψ_I(x) = H(0)ᵀ M(x)⁻¹ H(x_I - x) φ_I(x),   M(x) = Σ_I H(x_I - x) H(x_I - x)ᵀ φ_I(x)
//
// with H(d) = (1, d_x/L, d_y/L) and φ_I the cubic B-spline of |x - x_I|/a_I.
// The support a_I is the `support` factor times the spacing of point I.
// The functions reproduce 1, x and y, so they pass the patch test of
// order 1. L is the largest support of the points near x: it scales H so
// that M is well conditioned in each unit system.
import std;
import math;
namespace meshless;
import .cloud;

var support: f64[>..] = [];     // a of each point
var bucket_size = 0.0;          // at least the largest support
var bucket_n: i64[2] = [0, 0];
var bucket_start: i64[>..] = [];
var bucket_items: i64[>..] = [];

// The points near the point of the last shape_at, and Ψ, ∂Ψ/∂x and ∂Ψ/∂y
// of each.
var near: i64[>..] = [];
var psi: f64[>..] = [];
var psi_x: f64[>..] = [];
var psi_y: f64[>..] = [];

// Gives each point the support `factor` times its spacing, and sorts the
// points into buckets of the largest support.
fn set_support(factor: f64) {
    guard factor > 0.0 && spacing.len == point_count() else { abort("set_support: factor > 0 and a cloud"); }
    support.clear();
    bucket_size = 0.0;
    for p in point_count() {
        support.push(factor * spacing[p]);
        bucket_size = max(bucket_size, factor * spacing[p]);
    }
    for j in 2 { bucket_n[j] = min(max(1, ceil((hi[j] - lo[j]) / bucket_size) as! i64), 4096); }
    let count = bucket_n[0] * bucket_n[1];
    bucket_start.clear();
    for b in count + 1 { bucket_start.push(0); }
    for p in point_count() { bucket_start[bucket_of(px(p), py(p)) + 1]++; }
    for b in count { bucket_start[b + 1] += bucket_start[b]; }
    bucket_items.clear();
    for p in point_count() { bucket_items.push(0); }
    var at: i64[>..] = [];
    at.append(bucket_start[..count]);
    for p in point_count() {
        let b = bucket_of(px(p), py(p));
        bucket_items[at[b]] = p;
        at[b]++;
    }
}

fn bucket_axis(j: i64, v: f64) -> i64 {
    return clamp(floor((v - lo[j]) / (bucket_size * (1.0 + 1.0e-12))) as! i64, 0, bucket_n[j] - 1);
}

fn bucket_of(x: f64, y: f64) -> i64 { return bucket_axis(0, x) + bucket_n[0] * bucket_axis(1, y); }

// The points within `radius` of (x, y), or within their own support for
// a radius <= 0, in `into`.
fn points_near(x: f64, y: f64, radius: f64, into: i64[>..]&) {
    guard support.len == point_count() && bucket_size > 0.0 else { abort("points_near: set_support"); }
    into.clear();
    let reach = max(1, ceil(radius / bucket_size) as! i64);
    let bi = bucket_axis(0, x);
    let bj = bucket_axis(1, y);
    for j in max(0, bj - reach)..min(bucket_n[1], bj + reach + 1) {
        for i in max(0, bi - reach)..min(bucket_n[0], bi + reach + 1) {
            let b = i + bucket_n[0] * j;
            for p in bucket_items[bucket_start[b]..bucket_start[b + 1]] {
                let r = norm_2(px(p) - x, py(p) - y);
                if r < (if radius > 0.0 { radius } else { support[p] }) { into.push(p); }
            }
        }
    }
}

// The cubic B-spline of z = r/a, and its derivative.
fn kernel(z: f64) -> f64, f64 {
    if z <= 0.5 { return 2.0 / 3.0 - 4.0 * z * z + 4.0 * z * z * z, -8.0 * z + 12.0 * z * z; }
    if z < 1.0 { return 4.0 / 3.0 - 4.0 * z + 4.0 * z * z - 4.0 / 3.0 * z * z * z, -4.0 + 8.0 * z - 4.0 * z * z; }
    return 0.0, 0.0;
}

// The inverse of the symmetric 3 × 3 matrix m, and false if it is
// singular.
fn inverse_3(m: f64[9]) -> f64[9], bool {
    var c: f64[9] = [0.0; 9];
    c[0] = m[4] * m[8] - m[5] * m[7];
    c[1] = m[2] * m[7] - m[1] * m[8];
    c[2] = m[1] * m[5] - m[2] * m[4];
    c[3] = m[5] * m[6] - m[3] * m[8];
    c[4] = m[0] * m[8] - m[2] * m[6];
    c[5] = m[2] * m[3] - m[0] * m[5];
    c[6] = m[3] * m[7] - m[4] * m[6];
    c[7] = m[1] * m[6] - m[0] * m[7];
    c[8] = m[0] * m[4] - m[1] * m[3];
    let det = m[0] * c[0] + m[1] * c[3] + m[2] * c[6];
    if !(abs(det) > 1.0e-12 * abs(m[0] * m[4] * m[8])) { return c, false; }
    for k in 9 { c[k] /= det; }
    return c, true;
}

fn times_3(m: f64[9], v: f64[3]) -> f64[3] {
    return [m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
            m[6] * v[0] + m[7] * v[1] + m[8] * v[2]];
}

fn dot_3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// Ψ and its gradient at (x, y) for each point near it, in near, psi,
// psi_x and psi_y. False if M is singular there: too few points cover it.
fn shape_at(x: f64, y: f64) -> bool {
    points_near(x, y, 0.0, near);
    psi.clear();
    psi_x.clear();
    psi_y.clear();
    var l = 0.0;
    for p in near { l = max(l, support[p]); }
    var m: f64[9] = [0.0; 9];
    var mx: f64[9] = [0.0; 9];
    var my: f64[9] = [0.0; 9];
    for p in near {
        let h = [1.0, (px(p) - x) / l, (py(p) - y) / l];
        let f, fx, fy = weight(p, x, y);
        for r in 3 {
            for c in 3 {
                m[3 * r + c] += f * h[r] * h[c];
                // ∂H/∂x = (0, -1/L, 0), ∂H/∂y = (0, 0, -1/L)
                mx[3 * r + c] += fx * h[r] * h[c] - f * ((if r == 1 { h[c] } else { 0.0 }) + (if c == 1 { h[r] } else { 0.0 })) / l;
                my[3 * r + c] += fy * h[r] * h[c] - f * ((if r == 2 { h[c] } else { 0.0 }) + (if c == 2 { h[r] } else { 0.0 })) / l;
            }
        }
    }
    let inv, regular = inverse_3(m);
    if !regular { return false; }
    let b = [inv[0], inv[3], inv[6]];
    let bx = times_3(inv, times_3(mx, b));
    let by = times_3(inv, times_3(my, b));
    for p in near {
        let h = [1.0, (px(p) - x) / l, (py(p) - y) / l];
        let f, fx, fy = weight(p, x, y);
        let bh = dot_3(b, h);
        psi.push(bh * f);
        psi_x.push(-dot_3(bx, h) * f - b[1] / l * f + bh * fx);
        psi_y.push(-dot_3(by, h) * f - b[2] / l * f + bh * fy);
    }
    return true;
}

// φ_p at (x, y), and its gradient.
fn weight(p: i64, x: f64, y: f64) -> f64, f64, f64 {
    let dx = x - px(p);
    let dy = y - py(p);
    let r = norm_2(dx, dy);
    let a = support[p];
    let w, dw = kernel(r / a);
    if r <= 0.0 { return w, 0.0, 0.0; }
    return w, dw * dx / (r * a), dw * dy / (r * a);
}
