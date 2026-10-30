// The linear elastic solid: the material, the element routine of a 3D cell
// of any type in the registry, the Jacobian of a 3D cell, and the registry
// entry of the 27-node hexahedron (the triquadratic Lagrange cell).
//
// The local node a of a hex27 cell is at (a % 3, a / 3 % 3, a / 9) on the
// lattice of the reference points -1, 0 and 1 along each axis. Its local
// face 2 r + side is the face ξ_r = -1 (side 0) or ξ_r = 1 (side 1).
import std;
import math;
namespace runtime;
import .fem;

var young = 0.0;
var poisson = 0.0;

// The quadratic Lagrange function of node a (at -1, 0 and 1) at s.
fn shape_1d(a: i64, s: f64) -> f64 {
    if a == 0 { return 0.5 * s * (s - 1.0); }
    if a == 1 { return 1.0 - s * s; }
    return 0.5 * s * (s + 1.0);
}

fn shape_1d_ds(a: i64, s: f64) -> f64 {
    if a == 0 { return s - 0.5; }
    if a == 1 { return -2.0 * s; }
    return s + 0.5;
}

fn hex27_shape(xi: f64[3], n: f64[>..]&, dn: f64[>..]&) {
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

// The weights of local face f = 2 r + side of the hex27 cell c
// (cell_face_weights in fem.goose), by 3 × 3 Gauss points on the face.
fn hex27_face_weights(c: i64, f: i64, w: f64[>..]&) {
    let r = f / 2;
    let a = if r == 0 { 1 } else { 0 };
    let b = if r == 2 { 1 } else { 2 };
    for l in 27 { w.push(0.0); }
    var xi: f64[3] = [0.0, 0.0, 0.0];
    xi[r] = if f % 2 == 0 { -1.0 } else { 1.0 };
    for gb in 3 {
        for ga in 3 {
            xi[a] = GAUSS3_X[ga];
            xi[b] = GAUSS3_X[gb];
            shape_n.clear();
            shape_dn.clear();
            hex27_shape(xi, shape_n, shape_dn);
            // The tangents ∂x/∂ξ_a and ∂x/∂ξ_b of the face.
            var ta: f64[3] = [0.0, 0.0, 0.0];
            var tb: f64[3] = [0.0, 0.0, 0.0];
            for l in 27 {
                let n = cell_node(c, l);
                for s in 3 {
                    ta[s] += shape_dn[3 * l + a] * coords[3 * n + s];
                    tb[s] += shape_dn[3 * l + b] * coords[3 * n + s];
                }
            }
            let nx = ta[1] * tb[2] - ta[2] * tb[1];
            let ny = ta[2] * tb[0] - ta[0] * tb[2];
            let nz = ta[0] * tb[1] - ta[1] * tb[0];
            let da = GAUSS3_W[ga] * GAUSS3_W[gb] * sqrt(nx * nx + ny * ny + nz * nz);
            for l in 27 { w[l] += shape_n[l] * da; }
        }
    }
}

// jac[3 r + s] = ∂x_s / ∂ξ_r in cell c, from the shape function
// derivatives dn of cell_shape.
fn cell_jacobian(c: i64, dn: f64[:], jac: f64[9]&) {
    for i in 9 { jac[i] = 0.0; }
    for l in cell_start[c + 1] - cell_start[c] {
        let n = cell_node(c, l);
        for r in 3 {
            for s in 3 { jac[3 * r + s] += dn[3 * l + r] * coords[3 * n + s]; }
        }
    }
}

// inv = jac⁻¹, so ∂N/∂x_s = Σr inv[3 s + r] ∂N/∂ξ_r. It gives the
// determinant of jac, and leaves inv as it is if the determinant is 0.
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
var grad: f64[>..] = [];      // the gradient of each shape function at a quadrature point

// The stiffness of the 3D cell c, with 3 dofs for each node:
// K(i p, j q) = ∫ λ ∂p Ni ∂q Nj + μ (∂q Ni ∂p Nj + δpq ∇Ni · ∇Nj) dV.
fn solid_stiffness(c: i64, ke: f64[>..]&) -> bool {
    let t = cell_type[c] as i64;
    guard CELL_DIMS[t] == 3 && cell_quadrature(t, quad_xi, quad_w) else { return false; }
    let nn = CELL_NODES[t];
    let size = 3 * nn;
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    ke.clear();
    for i in size * size { ke.push(0.0); }
    grad.clear();
    for i in size { grad.push(0.0); }
    for w, g in quad_w {
        let xi: f64[3] = [quad_xi[3 * g], quad_xi[3 * g + 1], quad_xi[3 * g + 2]];
        guard cell_shape(t, xi, shape_n, shape_dn) else { return false; }
        var jac: f64[9] = [0.0; 9];
        cell_jacobian(c, shape_dn, jac);
        var inv: f64[9] = [0.0; 9];
        let det = inverse3(jac, inv);
        assert(det > 0.0);
        for l in nn {
            for s in 3 {
                grad[3 * l + s] = inv[3 * s] * shape_dn[3 * l] + inv[3 * s + 1] * shape_dn[3 * l + 1]
                    + inv[3 * s + 2] * shape_dn[3 * l + 2];
            }
        }
        let dv = w * det;
        for i in nn {
            for j in nn {
                let dot = grad[3 * i] * grad[3 * j] + grad[3 * i + 1] * grad[3 * j + 1]
                    + grad[3 * i + 2] * grad[3 * j + 2];
                for p in 3 {
                    let row = size * (3 * i + p) + 3 * j;
                    for q in 3 {
                        var v = lambda * grad[3 * i + p] * grad[3 * j + q]
                            + mu * grad[3 * i + q] * grad[3 * j + p];
                        if p == q { v += mu * dot; }
                        ke[row + q] += v * dv;
                    }
                }
            }
        }
    }
    return true;
}
