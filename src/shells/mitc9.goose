// The MITC9 shell: a 9-node Reissner-Mindlin shell with the mixed
// interpolation of Bucalem and Bathe (1993), "Higher-order MITC general
// shell elements", International Journal for Numerical Methods in
// Engineering 36, 3729-3754.
//
// The kinematics, the strains, the section and the drilling penalty are
// those of MITC4 (mitc4.goose), with the 9-node Lagrange shape functions.
// The node sequence is that of the comment on Shell (model.goose). The
// bending strains come from the displacements. The membrane strains and
// the transverse shear strains come from tying points, with a = 1/√3 and
// b = √(3/5):
//     E11 and G1 at (±a, -b), (±a, 0) and (±a, b),
//     E22 and G2 at (-b, ±a), (0, ±a) and (b, ±a),
//     E12 at (±a, ±a).
// The interpolation is linear between -a and a, and quadratic through -b,
// 0 and b. The stiffness and the loads use 3 × 3 Gauss points.
import std;
import math;
namespace shells;
import .laminate;
import .model;
import .mitc4;

let NODE9_I = [0, 2, 2, 0, 1, 2, 1, 0, 1];     // the place of each node along xi: 0, 1, 2 for -1, 0, 1
let NODE9_J = [0, 0, 2, 2, 0, 1, 2, 1, 1];     // the same along eta
let TIE_A = 0.5773502691896258;     // 1/√3
let TIE_B = 0.7745966692414834;     // √(3/5)
let GAUSS3 = [-0.7745966692414834, 0.0, 0.7745966692414834];
let GAUSS3_W = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

var ke9: f64[2916];     // 54 × 54, row 6 l + d for dof d of local node l
var xn9: f64[27];       // the nodes of the cell
var vn9: f64[27];       // the directors of the cell

// The covariant strain rows at a point, one row of 54 for each: E11, E22,
// E12, then K11, K22, K12, then G1, G2. covariant9 also sets a1 and a2 of
// mitc4.goose.
var cov9: f64[432];
var shape9: f64[9];
var dn9: f64[18];       // ∂Ni/∂ξ, then ∂Ni/∂η

// The covariant rows at the tying points of the cell (tie_cell9): E11 and
// G1 at the 6 points of the first family, E22 and G2 at the 6 points of
// the second, E12 at the 4 points of the third.
var tie9: f64[1512];
var mix9: f64[270];     // the tied E11, E22, E12, G1 and G2 at a point

// The strains in the local axes, as `rows` of mitc4.goose: eps 11, 22,
// 12; kappa 11, 22, 12; gamma 13, 23; the drilling rotation.
var rows9: f64[486];

// The quadratic Lagrange functions of the points -1, 0 and 1 at t, and
// their derivatives.
fn quad3(t: f64) -> f64[3] { return [0.5 * t * (t - 1.0), 1.0 - t * t, 0.5 * t * (t + 1.0)]; }
fn dquad3(t: f64) -> f64[3] { return [t - 0.5, -2.0 * t, t + 0.5]; }

// The linear functions of the points -a and a, and the quadratic functions
// of the points -b, 0 and b, at t.
fn lin_a(t: f64) -> f64[2] { return [0.5 * (1.0 - t / TIE_A), 0.5 * (1.0 + t / TIE_A)]; }
fn quad_b(t: f64) -> f64[3] {
    let b2 = TIE_B * TIE_B;
    return [0.5 * t * (t - TIE_B) / b2, 1.0 - t * t / b2, 0.5 * t * (t + TIE_B) / b2];
}

fn load_cell9(c: i64) {
    guard c >= 0 && c < shells.len && shells[c].count == 9 && directors.len == coords.len else { abort("load_cell9: 0 <= c < shells.len, a 9-node cell, and a director for each node"); }
    for l in 9 {
        let n = shells[c].n[l];
        for k in 3 {
            xn9[3 * l + k] = coords[3 * n + k];
            vn9[3 * l + k] = directors[3 * n + k];
        }
    }
}

// The shape functions, a1, a2 and the covariant rows cov9 at (xi, eta).
fn covariant9(xi: f64, eta: f64) {
    guard abs(xi) <= 1.0 && abs(eta) <= 1.0 else { abort("covariant9: -1 <= xi, eta <= 1"); }
    let fx = quad3(xi);
    let fe = quad3(eta);
    let gx = dquad3(xi);
    let ge = dquad3(eta);
    var b1: f64[3] = [0.0, 0.0, 0.0];
    var b2: f64[3] = [0.0, 0.0, 0.0];
    var vbar: f64[3] = [0.0, 0.0, 0.0];
    a1 = [0.0, 0.0, 0.0];
    a2 = [0.0, 0.0, 0.0];
    for l in 9 {
        shape9[l] = fx[NODE9_I[l]] * fe[NODE9_J[l]];
        dn9[l] = gx[NODE9_I[l]] * fe[NODE9_J[l]];
        dn9[9 + l] = fx[NODE9_I[l]] * ge[NODE9_J[l]];
        for k in 3 {
            a1[k] += dn9[l] * xn9[3 * l + k];
            a2[k] += dn9[9 + l] * xn9[3 * l + k];
            b1[k] += dn9[l] * vn9[3 * l + k];
            b2[k] += dn9[9 + l] * vn9[3 * l + k];
            vbar[k] += shape9[l] * vn9[3 * l + k];
        }
    }
    for &v in cov9 { v = 0.0; }
    for l in 9 {
        let v: f64[3] = [vn9[3 * l], vn9[3 * l + 1], vn9[3 * l + 2]];
        let d1 = dn9[l];
        let d2 = dn9[9 + l];
        // (V × a)[c] = a · (e_c × V): the effect of rotation c on a · r.
        let va1 = cross(v, a1);
        let va2 = cross(v, a2);
        for c in 3 {
            let u = 6 * l + c;
            let r = 6 * l + 3 + c;
            cov9[u] = a1[c] * d1;
            cov9[54 + u] = a2[c] * d2;
            cov9[108 + u] = 0.5 * (a1[c] * d2 + a2[c] * d1);
            cov9[162 + u] = b1[c] * d1;
            cov9[216 + u] = b2[c] * d2;
            cov9[270 + u] = 0.5 * (b1[c] * d2 + b2[c] * d1);
            cov9[162 + r] = va1[c] * d1;
            cov9[216 + r] = va2[c] * d2;
            cov9[270 + r] = 0.5 * (va1[c] * d2 + va2[c] * d1);
            cov9[324 + u] = vbar[c] * d1;
            cov9[378 + u] = vbar[c] * d2;
            cov9[324 + r] = shape9[l] * va1[c];
            cov9[378 + r] = shape9[l] * va2[c];
        }
    }
}

// Copies row `from` of cov9 (0 E11, 1 E22, 2 E12, 6 G1, 7 G2) to row `to`
// of tie9.
fn keep_tie(from: i64, to: i64) {
    guard from >= 0 && from < 8 && to >= 0 && to < 28 else { abort("keep_tie: 0 <= from < 8, and 0 <= to < 28"); }
    for i in 54 { tie9[54 * to + i] = cov9[54 * from + i]; }
}

// The rows tie9 of the cell in xn9 and vn9. Point p = i + 2 k of the first
// family is at (±a, (-b, 0, b)[k]), point p of the second family is at
// ((-b, 0, b)[k], ±a), and point i + 2 k of the third is at (±a, ±a).
fn tie_cell9() {
    let sa = [-TIE_A, TIE_A];
    let sb = [-TIE_B, 0.0, TIE_B];
    for k in 3 {
        for i in 2 {
            let p = i + 2 * k;
            covariant9(sa[i], sb[k]);
            keep_tie(0, 2 * p);
            keep_tie(6, 2 * p + 1);
            covariant9(sb[k], sa[i]);
            keep_tie(1, 12 + 2 * p);
            keep_tie(7, 12 + 2 * p + 1);
        }
    }
    for k in 2 {
        for i in 2 {
            covariant9(sa[i], sa[k]);
            keep_tie(2, 24 + i + 2 * k);
        }
    }
}

// The tied rows mix9 at (xi, eta): E11, E22, E12, G1 and G2, 54 each.
fn tied_strains9(xi: f64, eta: f64) {
    let lx = lin_a(xi);
    let le = lin_a(eta);
    let qx = quad_b(xi);
    let qe = quad_b(eta);
    for &v in mix9 { v = 0.0; }
    for k in 3 {
        for i in 2 {
            let p = i + 2 * k;
            let w1 = lx[i] * qe[k];
            let w2 = le[i] * qx[k];
            for j in 54 {
                mix9[j] += w1 * tie9[54 * (2 * p) + j];
                mix9[162 + j] += w1 * tie9[54 * (2 * p + 1) + j];
                mix9[54 + j] += w2 * tie9[54 * (12 + 2 * p) + j];
                mix9[216 + j] += w2 * tie9[54 * (13 + 2 * p) + j];
            }
        }
    }
    for k in 2 {
        for i in 2 {
            let w = lx[i] * le[k];
            for j in 54 { mix9[108 + j] += w * tie9[54 * (24 + i + 2 * k) + j]; }
        }
    }
}

// The local strain rows rows9 at (xi, eta) of cell c, from the tie9 of the
// cell: the tied membrane and shear strains, the bending strains of the
// point.
fn local_rows9(c: i64, xi: f64, eta: f64) {
    guard c >= 0 && c < shells.len && abs(xi) <= 1.0 && abs(eta) <= 1.0 else { abort("local_rows9: 0 <= c < shells.len, and -1 <= xi, eta <= 1"); }
    tied_strains9(xi, eta);
    covariant9(xi, eta);
    let c00, c01, c10, c11 = local_frame(c);
    for i in 54 {
        for s in 2 {
            // s = 0: the membrane rows of mix9; s = 1: the bending rows of cov9.
            let e11 = if s == 0 { mix9[i] } else { cov9[162 + i] };
            let e22 = if s == 0 { mix9[54 + i] } else { cov9[216 + i] };
            let e12 = if s == 0 { mix9[108 + i] } else { cov9[270 + i] };
            rows9[162 * s + i] = c00 * c00 * e11 + 2.0 * c00 * c10 * e12 + c10 * c10 * e22;
            rows9[162 * s + 54 + i] = c01 * c01 * e11 + 2.0 * c01 * c11 * e12 + c11 * c11 * e22;
            rows9[162 * s + 108 + i] = 2.0 * (c00 * c01 * e11 + (c00 * c11 + c10 * c01) * e12 + c10 * c11 * e22);
        }
        rows9[324 + i] = c00 * mix9[162 + i] + c10 * mix9[216 + i];
        rows9[378 + i] = c01 * mix9[162 + i] + c11 * mix9[216 + i];
    }
    // The drilling row: θ · e3 - ½ (∂u/∂x1 · e2 - ∂u/∂x2 · e1).
    for l in 9 {
        let dx1 = c00 * dn9[l] + c10 * dn9[9 + l];
        let dx2 = c01 * dn9[l] + c11 * dn9[9 + l];
        for k in 3 {
            rows9[432 + 6 * l + k] = -0.5 * (dx1 * e2[k] - dx2 * e1[k]);
            rows9[432 + 6 * l + 3 + k] = shape9[l] * e3[k];
        }
    }
}

// The product of a symmetric 3 × 3 matrix m and the column j of rows9 r0.
fn mat3_row9(m: f64[9], r0: i64, j: i64) -> f64[3] {
    guard (r0 == 0 || r0 == 162) && j >= 0 && j < 54 else { abort("mat3_row9: r0 is 0 or 162, and 0 <= j < 54"); }
    let x = rows9[r0 + j];
    let y = rows9[r0 + 54 + j];
    let z = rows9[r0 + 108 + j];
    return [m[0] * x + m[1] * y + m[2] * z, m[3] * x + m[4] * y + m[5] * z, m[6] * x + m[7] * y + m[8] * z];
}

// The stiffness ke9 of 9-node shell cell c, by 3 × 3 Gauss points.
fn mitc9_stiffness(c: i64) {
    guard c >= 0 && c < shells.len && shells[c].section >= 0 && shells[c].section < sections.len else { abort("mitc9_stiffness: 0 <= c < shells.len, and its section < sections.len"); }
    load_cell9(c);
    tie_cell9();
    let s = sections[shells[c].section];
    let kd = DRILL * s.a[8];
    for &v in ke9 { v = 0.0; }
    var gen: f64[432] = [0.0; 432];     // N, M and Q of each dof column j: 8 rows of 54
    for gp in 9 {
        local_rows9(c, GAUSS3[gp % 3], GAUSS3[gp / 3]);
        let w = GAUSS3_W[gp % 3] * GAUSS3_W[gp / 3] * area;
        for j in 54 {
            let am = mat3_row9(s.a, 0, j);
            let bk = mat3_row9(s.b, 162, j);
            let bm = mat3_row9(s.b, 0, j);
            let dk = mat3_row9(s.d, 162, j);
            for k in 3 {
                gen[54 * k + j] = w * (am[k] + bk[k]);
                gen[162 + 54 * k + j] = w * (bm[k] + dk[k]);
            }
            gen[324 + j] = w * (s.s[0] * rows9[324 + j] + s.s[1] * rows9[378 + j]);
            gen[378 + j] = w * (s.s[2] * rows9[324 + j] + s.s[3] * rows9[378 + j]);
        }
        for i in 54 {
            let di = w * kd * rows9[432 + i];
            for j in i..54 {
                var v = di * rows9[432 + j];
                for k in 8 { v += rows9[54 * k + i] * gen[54 * k + j]; }
                ke9[54 * i + j] += v;
            }
        }
    }
    for i in 54 {
        for j in i + 1..54 { ke9[54 * j + i] = ke9[54 * i + j]; }
    }
}

// Adds a load per area (qx, qy, qz) on 9-node cell c, or a pressure qx
// along the top normal with `normal` (shell_load).
fn shell_load9(c: i64, qx: f64, qy: f64, qz: f64, normal: bool) {
    guard c >= 0 && c < shells.len && load.len == dof_count() else { abort("shell_load9: 0 <= c < shells.len, and a load for each dof"); }
    load_cell9(c);
    for gp in 9 {
        covariant9(GAUSS3[gp % 3], GAUSS3[gp / 3]);
        let n = cross(a1, a2);
        let da = norm(n) * GAUSS3_W[gp % 3] * GAUSS3_W[gp / 3];
        let q: f64[3] = if normal { [qx * n[0] / norm(n), qx * n[1] / norm(n), qx * n[2] / norm(n)] } else { [qx, qy, qz] };
        for l in 9 {
            let node = shells[c].n[l];
            for k in 3 { load[6 * node + k] += shape9[l] * q[k] * da; }
        }
    }
}

// The membrane strain and the curvature of 9-node cell c at (xi, eta), as
// shell_strains.
fn shell_strains9(c: i64, xi: f64, eta: f64) -> f64[3], f64[3] {
    guard c >= 0 && c < shells.len && disp.len == dof_count() else { abort("shell_strains9: 0 <= c < shells.len, and a solution for each dof"); }
    load_cell9(c);
    tie_cell9();
    local_rows9(c, xi, eta);
    var eps: f64[3] = [0.0, 0.0, 0.0];
    var kap: f64[3] = [0.0, 0.0, 0.0];
    for l in 9 {
        let node = shells[c].n[l];
        for d in 6 {
            let u = disp[6 * node + d];
            for k in 3 {
                eps[k] += rows9[54 * k + 6 * l + d] * u;
                kap[k] += rows9[162 + 54 * k + 6 * l + d] * u;
            }
        }
    }
    return eps, kap;
}

// Adds the unit normal of 9-node cell c at each of its nodes to sum (3 for
// each node), for mean_directors.
fn add_normals9(c: i64, sum: f64[:]) {
    guard c >= 0 && c < shells.len && sum.len == coords.len else { abort("add_normals9: 0 <= c < shells.len, and sum has 3 values for each node"); }
    load_cell9(c);
    for l in 9 {
        covariant9(NODE9_I[l] - 1.0, NODE9_J[l] - 1.0);
        let n = unit(cross(a1, a2));
        for k in 3 { sum[3 * shells[c].n[l] + k] += n[k]; }
    }
}
