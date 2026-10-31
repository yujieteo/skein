// gfd: the generalized finite differences of Benito, Ureña and Gavete
// (2001), collocated at each point (spec §17.5).
//
// The stencil of point I is the `stencil` points nearest to it. The
// derivatives D = (u_x, u_y, u_xx, u_yy, u_xy) at I minimize
//
//     Σ_j w_j² (u_I + p_jᵀ D - u_j)²,   p_j = (Δx, Δy, Δx²/2, Δy²/2, Δx Δy)
//
// with Δ = x_j - x_I and the weight w_j = 1/|Δ|³. So D = Σ_j c_j (u_j - u_I),
// with c_j = A⁻¹ p_j w_j² and A = Σ_j w_j² p_j p_jᵀ. Each point gives one
// row: the equation at a point inside, u = g at a point of a fixed edge,
// and (a u - d grad(u))·n = 0 at a point of a free edge, with n the mean
// normal of its edges.
import std;
import math;
namespace meshless;
import .cloud;
import .rk;
import .model;
import .solve;

var st_start: i64[>..] = [];    // the stencil of each point, from st_start[p]
var st_item: i64[>..] = [];
var st_c: f64[>..] = [];        // c_j of each entry: 5 values
var on_edge: i64[>..] = [];     // a boundary edge of each point, or -1
var normal: f64[>..] = [];      // the mean outward normal of each point on the boundary

// The stencils, the normals, and the pattern. False if a stencil cannot
// fit the 5 derivatives: its points are on one conic.
fn gfd_setup() -> bool {
    st_start.clear();
    st_item.clear();
    st_c.clear();
    st_start.push(0);
    for p in point_count() { if !stencil_of(p) { return false; } }
    boundary_normals();
    start_matrix();
    var cols: i64[>..] = [];
    for p in point_count() {
        cols.clear();
        cols.push(p);
        cols.append(st_item[st_start[p]..st_start[p + 1]]);
        sort(cols);
        add_rows(cols);
    }
    size_matrix();
    return true;
}

fn boundary_normals() {
    on_edge.clear();
    normal.clear();
    for p in point_count() {
        on_edge.push(-1);
        normal.append([0.0, 0.0]);
    }
    for e in edge_count() {
        let n, l = edge_normal(e);
        for s in 2 {
            let p = edges[2 * e + s];
            on_edge[p] = e;
            normal[2 * p] += l * n[0];
            normal[2 * p + 1] += l * n[1];
        }
    }
    for p in point_count() {
        let l = norm_2(normal[2 * p], normal[2 * p + 1]);
        if l > 0.0 {
            normal[2 * p] /= l;
            normal[2 * p + 1] /= l;
        }
    }
}

// The `stencil` points nearest to p, and their c_j.
fn stencil_of(p: i64) -> bool {
    var radius = 2.0 * spacing[p];
    var cand: i64[>..] = [];
    for tries in 8 {
        points_near(px(p), py(p), radius, cand);
        if cand.len > stencil { break; }
        radius *= 1.5;
    }
    guard cand.len > stencil else { return false; };
    var dist: f64[>..] = [];
    for j in cand { dist.push(if j == p { 1.0e300 } else { norm_2(px(j) - px(p), py(j) - py(p)) }); }
    var chosen: i64[>..] = [];
    for k in stencil {
        var best = 0;
        for i in dist.len { if dist[i] < dist[best] { best = i; } }
        chosen.push(cand[best]);
        dist[best] = 1.0e300;
    }
    let h = spacing[p];
    var a: f64[25] = [0.0; 25];
    for j in chosen {
        let q, w2 = taylor(p, j, h);
        for r in 5 { for c in 5 { a[5 * r + c] += w2 * q[r] * q[c]; } }
    }
    let inv, regular = inverse_5(a);
    if !regular { return false; }
    let scale = [1.0 / h, 1.0 / h, 1.0 / (h * h), 1.0 / (h * h), 1.0 / (h * h)];
    for j in chosen {
        let q, w2 = taylor(p, j, h);
        st_item.push(j);
        for r in 5 {
            var s = 0.0;
            for c in 5 { s += inv[5 * r + c] * q[c]; }
            st_c.push(s * w2 * scale[r]);
        }
    }
    st_start.push(st_item.len);
    return true;
}

// p_j of point j from point p, in units of h, and w_j².
fn taylor(p: i64, j: i64, h: f64) -> f64[5], f64 {
    let dx = (px(j) - px(p)) / h;
    let dy = (py(j) - py(p)) / h;
    let r = norm_2(dx, dy);
    let w = 1.0 / (r * r * r);
    return [dx, dy, 0.5 * dx * dx, 0.5 * dy * dy, dx * dy], w * w;
}

// The inverse of the 5 × 5 matrix a by Gauss-Jordan with partial
// pivoting, and false if it is singular.
fn inverse_5(a: f64[25]) -> f64[25], bool {
    var m = a;
    var inv: f64[25] = [0.0; 25];
    for i in 5 { inv[6 * i] = 1.0; }
    var largest = 0.0;
    for v in a { largest = max(largest, abs(v)); }
    for c in 5 {
        var piv = c;
        for r in c + 1..5 { if abs(m[5 * r + c]) > abs(m[5 * piv + c]) { piv = r; } }
        if !(abs(m[5 * piv + c]) > 1.0e-12 * largest) { return inv, false; }
        for k in 5 {
            let t = m[5 * c + k];
            m[5 * c + k] = m[5 * piv + k];
            m[5 * piv + k] = t;
            let u = inv[5 * c + k];
            inv[5 * c + k] = inv[5 * piv + k];
            inv[5 * piv + k] = u;
        }
        let d = m[6 * c];
        for k in 5 {
            m[5 * c + k] /= d;
            inv[5 * c + k] /= d;
        }
        for r in 5 {
            if r == c { continue; }
            let f = m[5 * r + c];
            for k in 5 {
                m[5 * r + k] -= f * m[5 * c + k];
                inv[5 * r + k] -= f * inv[5 * c + k];
            }
        }
    }
    return inv, true;
}

// The row of each point.
fn gfd_assemble<SRC, BC>() {
    for &v in val { v = 0.0; }
    var fix: i64[>..] = [];
    for p in point_count() { fix.push(FREE); }
    for e in edge_count() {
        if edge_fix[e] == FREE { continue; }
        fix[edges[2 * e]] = edge_fix[e];
        fix[edges[2 * e + 1]] = edge_fix[e];
    }
    for p in point_count() {
        let x = px(p);
        let y = py(p);
        if fix[p] != FREE {
            val[entry(p, p)] = 1.0;
            rhs[p] = BC(fix[p], x, y, time);
            continue;
        }
        let free = on_edge[p] >= 0;
        let n = [normal[2 * p], normal[2 * p + 1]];
        // The row: Σ_j coef_j (u_j - u_I).
        var centre = if free { ax * n[0] + ay * n[1] } else { 0.0 };
        for i in st_start[p]..st_start[p + 1] {
            let c = st_c[5 * i..5 * i + 5];
            let coef = if free { -diffusivity * (n[0] * c[0] + n[1] * c[1]) }
                       else { -diffusivity * (c[2] + c[3]) + ax * c[0] + ay * c[1] };
            val[entry(p, st_item[i])] += coef;
            centre -= coef;
        }
        val[entry(p, p)] += centre;
        rhs[p] = if free { 0.0 } else { SRC(x, y, time) };
    }
}
