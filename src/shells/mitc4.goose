// The MITC4 shell: a 4-node Reissner-Mindlin shell with the transverse
// shear strains of Dvorkin and Bathe (1984), "A continuum mechanics based
// four-node shell element for general non-linear analysis", Engineering
// Computations 1, 77-88.
//
// The shell is the mid-surface x = Σ Ni xi with a director Vi at each node.
// The displacement at the distance z from the mid-surface along V is
//     u + z (θ × V),
// with θ the rotation vector of the node. The strains are linear in z, with
// the covariant terms at the mid-surface (a thin shell):
//     E_ab = ½ (a_a · v_b + a_b · v_a)                      membrane
//     K_ab = ½ (a_a · w_b + b_a · v_b + a_b · w_a + b_b · v_a)  bending
//     G_a  = a_a · r + V · v_a                              shear
// where a_a = ∂x/∂s_a, b_a = ∂V/∂s_a, v_a = ∂u/∂s_a, w_a = ∂(θ × V)/∂s_a and
// r = θ × V. MITC4 takes G_1 from the points (0, ±1) and G_2 from (±1, 0).
// The section (laminate.goose) gives N, M and Q from the strains in the
// local axes (e1, e2, e3): e3 is the normal a_1 × a_2, e1 is the axis of
// the section on the surface.
//
// The rotation about the normal (the drilling rotation) has no stiffness
// in the shell. A small penalty ties it to the rotation of the membrane,
// ½ (∂u2/∂x1 - ∂u1/∂x2), as in Hughes and Brezzi (1989).
import std;
import math;
namespace shells;
import .laminate;
import .model;

let NODE_XI = [-1.0, 1.0, 1.0, -1.0];
let NODE_ETA = [-1.0, -1.0, 1.0, 1.0];
let GAUSS = 0.5773502691896258;     // 1/√3
let DRILL = 1.0e-3;                  // the drilling penalty relative to A66

var ke: f64[576];       // 24 × 24, row 6 l + d for dof d of local node l
var xn: f64[12];        // the nodes of the cell
var vn: f64[12];        // the directors of the cell

// The covariant strains at a point, one row of 24 for each: E11, E22, E12,
// then K11, K22, K12, then G1, G2.
var cov: f64[192];
var a1: f64[3];
var a2: f64[3];
var shape: f64[4];

// The strains in the local axes: eps 11, 22, 12; kappa 11, 22, 12; gamma
// 13, 23; the drilling rotation.
var rows: f64[216];
var e1: f64[3];
var e2: f64[3];
var e3: f64[3];
var area = 0.0;     // |a_1 × a_2|

fn load_cell(c: i64) {
    guard c >= 0 && c < shells.len && shells[c].count == 4 && directors.len == coords.len else { abort("load_cell: 0 <= c < shells.len, a 4-node cell, and a director for each node"); }
    for l in 4 {
        let n = shells[c].n[l];
        for k in 3 {
            xn[3 * l + k] = coords[3 * n + k];
            vn[3 * l + k] = directors[3 * n + k];
        }
    }
}

fn node_vec(v: f64[12], l: i64) -> f64[3] {
    guard l >= 0 && l < 4 else { abort("node_vec: 0 <= l < 4"); }
    return [v[3 * l], v[3 * l + 1], v[3 * l + 2]];
}

// The covariant strain rows `cov` at (xi, eta), with a1, a2 and the shape
// functions.
fn covariant(xi: f64, eta: f64) {
    guard abs(xi) <= 1.0 && abs(eta) <= 1.0 else { abort("covariant: -1 <= xi, eta <= 1"); }
    var dn: f64[8] = [0.0; 8];     // ∂Ni/∂ξ, then ∂Ni/∂η
    for l in 4 {
        shape[l] = 0.25 * (1.0 + xi * NODE_XI[l]) * (1.0 + eta * NODE_ETA[l]);
        dn[l] = 0.25 * NODE_XI[l] * (1.0 + eta * NODE_ETA[l]);
        dn[4 + l] = 0.25 * NODE_ETA[l] * (1.0 + xi * NODE_XI[l]);
    }
    var b1: f64[3] = [0.0, 0.0, 0.0];
    var b2: f64[3] = [0.0, 0.0, 0.0];
    var vbar: f64[3] = [0.0, 0.0, 0.0];
    a1 = [0.0, 0.0, 0.0];
    a2 = [0.0, 0.0, 0.0];
    for l in 4 {
        for k in 3 {
            a1[k] += dn[l] * xn[3 * l + k];
            a2[k] += dn[4 + l] * xn[3 * l + k];
            b1[k] += dn[l] * vn[3 * l + k];
            b2[k] += dn[4 + l] * vn[3 * l + k];
            vbar[k] += shape[l] * vn[3 * l + k];
        }
    }
    for &v in cov { v = 0.0; }
    for l in 4 {
        let v = node_vec(vn, l);
        let d1 = dn[l];
        let d2 = dn[4 + l];
        // (V × a)[c] = a · (e_c × V): the effect of rotation c on a · r.
        let va1 = cross(v, a1);
        let va2 = cross(v, a2);
        for c in 3 {
            let u = 6 * l + c;
            let r = 6 * l + 3 + c;
            cov[u] = a1[c] * d1;
            cov[24 + u] = a2[c] * d2;
            cov[48 + u] = 0.5 * (a1[c] * d2 + a2[c] * d1);
            cov[72 + u] = b1[c] * d1;
            cov[96 + u] = b2[c] * d2;
            cov[120 + u] = 0.5 * (b1[c] * d2 + b2[c] * d1);
            cov[72 + r] = va1[c] * d1;
            cov[96 + r] = va2[c] * d2;
            cov[120 + r] = 0.5 * (va1[c] * d2 + va2[c] * d1);
            cov[144 + u] = vbar[c] * d1;
            cov[168 + u] = vbar[c] * d2;
            cov[144 + r] = shape[l] * va1[c];
            cov[168 + r] = shape[l] * va2[c];
        }
    }
}

// The local strain rows at (xi, eta) of cell c, with the MITC4 shear.
fn local_rows(c: i64, xi: f64, eta: f64) {
    guard c >= 0 && c < shells.len && abs(xi) <= 1.0 && abs(eta) <= 1.0 else { abort("local_rows: 0 <= c < shells.len, and -1 <= xi, eta <= 1"); }
    // The tying points: G1 at A (0, -1) and C (0, 1), G2 at D (-1, 0) and B (1, 0).
    var g1: f64[24] = [0.0; 24];
    var g2: f64[24] = [0.0; 24];
    let tie_xi = [0.0, 0.0, -1.0, 1.0];
    let tie_eta = [-1.0, 1.0, 0.0, 0.0];
    let tie_w = [0.5 * (1.0 - eta), 0.5 * (1.0 + eta), 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)];
    for t in 4 {
        covariant(tie_xi[t], tie_eta[t]);
        for i in 24 {
            if t < 2 { g1[i] += tie_w[t] * cov[144 + i]; } else { g2[i] += tie_w[t] * cov[168 + i]; }
        }
    }
    covariant(xi, eta);
    let c00, c01, c10, c11 = local_frame(c);
    for i in 24 {
        for s in 2 {
            let e11 = cov[72 * s + i];
            let e22 = cov[72 * s + 24 + i];
            let e12 = cov[72 * s + 48 + i];
            rows[72 * s + i] = c00 * c00 * e11 + 2.0 * c00 * c10 * e12 + c10 * c10 * e22;
            rows[72 * s + 24 + i] = c01 * c01 * e11 + 2.0 * c01 * c11 * e12 + c11 * c11 * e22;
            rows[72 * s + 48 + i] = 2.0 * (c00 * c01 * e11 + (c00 * c11 + c10 * c01) * e12 + c10 * c11 * e22);
        }
        rows[144 + i] = c00 * g1[i] + c10 * g2[i];
        rows[168 + i] = c01 * g1[i] + c11 * g2[i];
    }
    drilling_rows(xi, eta, c00, c01, c10, c11);
}

// The local axes e1, e2, e3 and the area of cell c from a1 and a2, and
// c[2 a + k] = g^a · e_k from the inverse of the metric of a1 and a2.
fn local_frame(c: i64) -> f64, f64, f64, f64 {
    guard c >= 0 && c < shells.len else { abort("local_frame: 0 <= c < shells.len"); }
    let n = cross(a1, a2);
    area = norm(n);
    e3 = unit(n);
    let axis = shells[c].axis;
    var t = sub(axis, [dot3(axis, e3) * e3[0], dot3(axis, e3) * e3[1], dot3(axis, e3) * e3[2]]);
    if norm(t) < 1.0e-8 * norm(axis) || norm(axis) == 0.0 { t = a1; }
    e1 = unit(t);
    e2 = cross(e3, e1);
    let g11 = dot3(a1, a1);
    let g12 = dot3(a1, a2);
    let g22 = dot3(a2, a2);
    let det = g11 * g22 - g12 * g12;
    var up1: f64[3] = [0.0; 3];
    var up2: f64[3] = [0.0; 3];
    for k in 3 {
        up1[k] = (g22 * a1[k] - g12 * a2[k]) / det;
        up2[k] = (g11 * a2[k] - g12 * a1[k]) / det;
    }
    return dot3(up1, e1), dot3(up1, e2), dot3(up2, e1), dot3(up2, e2);
}

// The drilling row at (xi, eta): θ · e3 - ½ (∂u/∂x1 · e2 - ∂u/∂x2 · e1),
// with c[2 a + k] = g^a · e_k (local_rows).
fn drilling_rows(xi: f64, eta: f64, c00: f64, c01: f64, c10: f64, c11: f64) {
    guard abs(xi) <= 1.0 && abs(eta) <= 1.0 else { abort("drilling_rows: -1 <= xi, eta <= 1"); }
    for l in 4 {
        let d1 = 0.25 * NODE_XI[l] * (1.0 + eta * NODE_ETA[l]);
        let d2 = 0.25 * NODE_ETA[l] * (1.0 + xi * NODE_XI[l]);
        let dx1 = c00 * d1 + c10 * d2;
        let dx2 = c01 * d1 + c11 * d2;
        for k in 3 {
            rows[192 + 6 * l + k] = -0.5 * (dx1 * e2[k] - dx2 * e1[k]);
            rows[192 + 6 * l + 3 + k] = shape[l] * e3[k];
        }
    }
}

// The product of a symmetric 3 × 3 matrix m and the column j of rows r0.
fn mat3_row(m: f64[9], r0: i64, j: i64) -> f64[3] {
    guard (r0 == 0 || r0 == 72) && j >= 0 && j < 24 else { abort("mat3_row: r0 is 0 or 72, and 0 <= j < 24"); }
    let x = rows[r0 + j];
    let y = rows[r0 + 24 + j];
    let z = rows[r0 + 48 + j];
    return [m[0] * x + m[1] * y + m[2] * z, m[3] * x + m[4] * y + m[5] * z, m[6] * x + m[7] * y + m[8] * z];
}

// Adds a load per area (qx, qy, qz) on cell c. With `normal`, the load is
// a pressure qx along the top normal e3 instead.
fn shell_load(c: i64, qx: f64, qy: f64, qz: f64, normal: bool) {
    guard c >= 0 && c < shells.len && load.len == dof_count() else { abort("shell_load: 0 <= c < shells.len, and a load for each dof"); }
    if shells[c].count == 9 {
        shell_load9(c, qx, qy, qz, normal);
        return;
    }
    load_cell(c);
    for gp in 4 {
        let xi = if gp % 2 == 0 { -GAUSS } else { GAUSS };
        let eta = if gp < 2 { -GAUSS } else { GAUSS };
        covariant(xi, eta);
        let n = cross(a1, a2);
        let da = norm(n);
        let q: f64[3] = if normal { [qx * n[0] / da, qx * n[1] / da, qx * n[2] / da] } else { [qx, qy, qz] };
        for l in 4 {
            let node = shells[c].n[l];
            for k in 3 { load[6 * node + k] += shape[l] * q[k] * da; }
        }
    }
}

// The membrane strain and the curvature of cell c at (xi, eta) in the axes
// of the section: (11, 22, engineering 12).
fn shell_strains(c: i64, xi: f64, eta: f64) -> f64[3], f64[3] {
    guard c >= 0 && c < shells.len && disp.len == dof_count() else { abort("shell_strains: 0 <= c < shells.len, and a solution for each dof"); }
    if shells[c].count == 9 { return shell_strains9(c, xi, eta); }
    load_cell(c);
    local_rows(c, xi, eta);
    var eps: f64[3] = [0.0, 0.0, 0.0];
    var kap: f64[3] = [0.0, 0.0, 0.0];
    for l in 4 {
        let node = shells[c].n[l];
        for d in 6 {
            let u = disp[6 * node + d];
            for k in 3 {
                eps[k] += rows[24 * k + 6 * l + d] * u;
                kap[k] += rows[72 + 24 * k + 6 * l + d] * u;
            }
        }
    }
    return eps, kap;
}

// The membrane force N and the membrane strain of cell c at its centre
// (MITC4 or MITC9), as
// 3 × 3 tensors in global axes (row-major).
fn shell_resultants(c: i64) -> f64[9], f64[9] {
    return shell_resultants_at(c, 0.0, 0.0);
}

// shell_resultants at (xi, eta) of cell c.
fn shell_resultants_at(c: i64, xi: f64, eta: f64) -> f64[9], f64[9] {
    guard c >= 0 && c < shells.len && shells[c].section >= 0 && shells[c].section < sections.len else { abort("shell_resultants: 0 <= c < shells.len, and its section < sections.len"); }
    let s = sections[shells[c].section];
    let eps, kap = shell_strains(c, xi, eta);
    var nl: f64[3] = [0.0; 3];
    for k in 3 {
        nl[k] = s.a[3 * k] * eps[0] + s.a[3 * k + 1] * eps[1] + s.a[3 * k + 2] * eps[2]
            + s.b[3 * k] * kap[0] + s.b[3 * k + 1] * kap[1] + s.b[3 * k + 2] * kap[2];
    }
    var nt: f64[9] = [0.0; 9];
    var et: f64[9] = [0.0; 9];
    for i in 3 {
        for j in 3 {
            nt[3 * i + j] = nl[0] * e1[i] * e1[j] + nl[1] * e2[i] * e2[j] + nl[2] * (e1[i] * e2[j] + e2[i] * e1[j]);
            et[3 * i + j] = eps[0] * e1[i] * e1[j] + eps[1] * e2[i] * e2[j]
                + 0.5 * eps[2] * (e1[i] * e2[j] + e2[i] * e1[j]);
        }
    }
    return nt, et;
}

// The mean of shell_resultants over the shell cells at node n.
fn mean_resultants(n: i64) -> f64[9], f64[9] {
    guard n >= 0 && n < node_count() else { abort("mean_resultants: 0 <= n < node_count()"); }
    var nt: f64[9] = [0.0; 9];
    var et: f64[9] = [0.0; 9];
    var count = 0;
    for c in shells.len {
        if position(shells[c].n[..shells[c].count], n) < 0 { continue; }
        let a, b = shell_resultants(c);
        for k in 9 {
            nt[k] += a[k];
            et[k] += b[k];
        }
        count++;
    }
    for k in 9 {
        nt[k] /= max(count, 1);
        et[k] /= max(count, 1);
    }
    return nt, et;
}

// Sets the director of each node that has none to the mean normal of the
// shell cells at the node.
fn mean_directors() {
    guard directors.len == coords.len else { abort("mean_directors: a director for each node"); }
    var sum: f64[>..] = [];
    for i in coords.len { sum.push(0.0); }
    for c in shells.len {
        if shells[c].count == 9 {
            add_normals9(c, sum);
            continue;
        }
        load_cell(c);
        for l in 4 {
            covariant(NODE_XI[l], NODE_ETA[l]);
            let n = unit(cross(a1, a2));
            for k in 3 { sum[3 * shells[c].n[l] + k] += n[k]; }
        }
    }
    for node in node_count() {
        let d = node_vec3(directors, node);
        let s = node_vec3(sum, node);
        if norm(d) == 0.0 && norm(s) > 0.0 { set_director(node, s[0], s[1], s[2]); }
    }
}

fn node_vec3(v: f64[:], n: i64) -> f64[3] {
    guard n >= 0 && 3 * n + 2 < v.len else { abort("node_vec3: 0 <= n < v.len / 3"); }
    return [v[3 * n], v[3 * n + 1], v[3 * n + 2]];
}
