import std;
import math;
namespace runtime;
import .fem;

var young = 0.0;
var poisson = 0.0;

var cell_young: f64[>..] = [];
var cell_poisson: f64[>..] = [];
var cell_strain: f64[>..] = [];

fn cell_lame(c: i64) -> f64, f64 {
    guard c >= 0 && c < cell_young.len && cell_poisson.len == cell_young.len else { abort("cell_lame: a cell with a material"); }
    let e = cell_young[c];
    let nu = cell_poisson[c];
    return e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), e / (2.0 * (1.0 + nu));
}

fn shape_1d(a: i64, s: f64) -> f64 {
    guard a >= 0 && a < 3 else { abort("shape_1d: 0 <= a < 3"); }
    if a == 0 { return 0.5 * s * (s - 1.0); }
    if a == 1 { return 1.0 - s * s; }
    return 0.5 * s * (s + 1.0);
}

fn shape_1d_ds(a: i64, s: f64) -> f64 {
    guard a >= 0 && a < 3 else { abort("shape_1d_ds: 0 <= a < 3"); }
    if a == 0 { return s - 0.5; }
    if a == 1 { return -2.0 * s; }
    return s + 0.5;
}

fn hex27_shape(xi: f64[3], n: f64[>..]&, dn: f64[>..]&) {
    guard n.len == 0 && dn.len == 0 else { abort("hex27_shape: empty n and dn"); }
    for l in 27 {
        let a = l % 3;
        let b = l / 3 % 3;
        let c = l / 9;
        n.push(shape_1d(a, xi[0]) * shape_1d(b, xi[1]) * shape_1d(c, xi[2]));
        dn.push(shape_1d_ds(a, xi[0]) * shape_1d(b, xi[1]) * shape_1d(c, xi[2]));
        dn.push(shape_1d(a, xi[0]) * shape_1d_ds(b, xi[1]) * shape_1d(c, xi[2]));
        dn.push(shape_1d(a, xi[0]) * shape_1d(b, xi[1]) * shape_1d_ds(c, xi[2]));
    }
}

fn hex_face_weights(c: i64, f: i64, w: f64[>..]&) {
    guard c >= 0 && c < cell_count() && (cell_type[c] as i64 == HEX27 || cell_type[c] as i64 == HEX8I) && f >= 0 && f < 6 && w.len == 0 else { abort("hex_face_weights: a hexahedron, a local face and an empty w"); }
    let t = cell_type[c] as i64;
    let nn = CELL_NODES[t];
    let r = f / 2;
    let a = if r == 0 { 1 } else { 0 };
    let b = if r == 2 { 1 } else { 2 };
    var sign = if r == 1 { -1.0 } else { 1.0 };
    if f % 2 == 0 { sign = -sign; }
    for l in 4 * nn { w.push(0.0); }
    var xi: f64[3] = [0.0, 0.0, 0.0];
    xi[r] = if f % 2 == 0 { -1.0 } else { 1.0 };
    for gb in 3 {
        for ga in 3 {
            xi[a] = GAUSS3_X[ga];
            xi[b] = GAUSS3_X[gb];
            cell_shape(t, xi, shape_n, shape_dn);
            var ta: f64[3] = [0.0, 0.0, 0.0];
            var tb: f64[3] = [0.0, 0.0, 0.0];
            for l in nn {
                let n = cell_node(c, l);
                for s in 3 {
                    ta[s] += shape_dn[3 * l + a] * coords[3 * n + s];
                    tb[s] += shape_dn[3 * l + b] * coords[3 * n + s];
                }
            }
            let nx = sign * (ta[1] * tb[2] - ta[2] * tb[1]);
            let ny = sign * (ta[2] * tb[0] - ta[0] * tb[2]);
            let nz = sign * (ta[0] * tb[1] - ta[1] * tb[0]);
            let gw = GAUSS3_W[ga] * GAUSS3_W[gb];
            let da = gw * sqrt(nx * nx + ny * ny + nz * nz);
            for l in nn {
                w[4 * l] += shape_n[l] * da;
                w[4 * l + 1] += shape_n[l] * gw * nx;
                w[4 * l + 2] += shape_n[l] * gw * ny;
                w[4 * l + 3] += shape_n[l] * gw * nz;
            }
        }
    }
}

fn cell_jacobian(c: i64, dn: f64[:], jac: f64[9]&) {
    guard c >= 0 && c < cell_count() && dn.len == 3 * (cell_start[c + 1] - cell_start[c]) else { abort("cell_jacobian: a cell and 3 derivatives for each node"); }
    for i in 9 { jac[i] = 0.0; }
    for l in cell_start[c + 1] - cell_start[c] {
        let n = cell_node(c, l);
        for r in 3 {
            for s in 3 { jac[3 * r + s] += dn[3 * l + r] * coords[3 * n + s]; }
        }
    }
}

fn inverse3(jac: f64[9], inv: f64[9]&) -> f64 {
    let det = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7])
        - jac[1] * (jac[3] * jac[8] - jac[5] * jac[6])
        + jac[2] * (jac[3] * jac[7] - jac[4] * jac[6]);
    guard det != 0.0 else { return det; }
    inv[0] = (jac[4] * jac[8] - jac[5] * jac[7]) / det;
    inv[1] = (jac[2] * jac[7] - jac[1] * jac[8]) / det;
    inv[2] = (jac[1] * jac[5] - jac[2] * jac[4]) / det;
    inv[3] = (jac[5] * jac[6] - jac[3] * jac[8]) / det;
    inv[4] = (jac[0] * jac[8] - jac[2] * jac[6]) / det;
    inv[5] = (jac[2] * jac[3] - jac[0] * jac[5]) / det;
    inv[6] = (jac[3] * jac[7] - jac[4] * jac[6]) / det;
    inv[7] = (jac[1] * jac[6] - jac[0] * jac[7]) / det;
    inv[8] = (jac[0] * jac[4] - jac[1] * jac[3]) / det;
    return det;
}

var quad_xi: f64[>..] = [];
var quad_w: f64[>..] = [];
var shape_n: f64[>..] = [];
var shape_dn: f64[>..] = [];
var grad: f64[>..] = [];
