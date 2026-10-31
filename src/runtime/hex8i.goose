// The registry entry of the 8-node hexahedron with incompatible modes
// (Wilson et al. 1973, Taylor et al. 1976): the trilinear cell of C3D8I
// and of the femahoy grid solver (spec §13, `family: incompatible`).
//
// Each cell adds the modes P_m = 1 - ξ_m², m = 0, 1, 2, to each component
// of the displacement. The modes are internal to the cell: the element
// routine condenses them out, and hex8i_add_modes gets them again from the
// nodal displacements for the stress at a point (solve.goose).
//
// The local node a is at the corner HEX8_CORNERS[3 a .. 3 a + 3], in the
// VTK order of a hexahedron. Its local face 2 r + side is the face
// ξ_r = -1 (side 0) or ξ_r = 1 (side 1), as for hex27.
import std;
namespace runtime;
import .fem;
import .solid;
import .solve;

let HEX8_CORNERS: f64[24] = [-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0,
                             -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0];

let GAUSS2_X = 0.5773502691896258;   // 1/√3

// The 8 nodes and the 3 modes, with 3 dofs each.
let NODE_DOFS = 24;
let MODE_DOFS = 9;
let ALL_DOFS = 33;

fn hex8_shape(xi: f64[3], n: f64[>..]&, dn: f64[>..]&) {
    guard n.len == 0 && dn.len == 0 else { abort("hex8_shape: empty n and dn"); }
    for a in 8 {
        let ca: f64[3] = [HEX8_CORNERS[3 * a], HEX8_CORNERS[3 * a + 1], HEX8_CORNERS[3 * a + 2]];
        let f0 = 1.0 + ca[0] * xi[0];
        let f1 = 1.0 + ca[1] * xi[1];
        let f2 = 1.0 + ca[2] * xi[2];
        n.push(0.125 * f0 * f1 * f2);
        dn.push(0.125 * ca[0] * f1 * f2);
        dn.push(0.125 * f0 * ca[1] * f2);
        dn.push(0.125 * f0 * f1 * ca[2]);
    }
}

// The tensor product of 2-point Gauss rules on [-1, 1]^3.
fn gauss_hex2(xi: f64[>..]&, w: f64[>..]&) {
    guard xi.len == 0 && w.len == 0 else { abort("gauss_hex2: empty xi and w"); }
    for a in 8 {
        for r in 3 { xi.push(HEX8_CORNERS[3 * a + r] * GAUSS2_X); }
        w.push(1.0);
    }
}

var centre_n: f64[>..] = [];
var centre_dn: f64[>..] = [];

// Sets g[at + 3 m + s] = ∂P_m/∂x_s at the reference point xi of cell c,
// where the Jacobian determinant is det. It uses the Jacobian at the
// centre of the cell, times det0 / det (Taylor et al. 1976), so that the
// cell passes the patch test.
fn mode_gradient(c: i64, xi: f64[3], det: f64, g: f64[>..]&, at: i64) {
    guard det != 0.0 && at >= 0 && g.len == at + MODE_DOFS else { abort("mode_gradient: det != 0 and 9 entries of g from at"); }
    centre_n.clear();
    centre_dn.clear();
    hex8_shape([0.0, 0.0, 0.0], centre_n, centre_dn);
    var jac: f64[9] = [0.0; 9];
    cell_jacobian(c, centre_dn, jac);
    var inv: f64[9] = [0.0; 9];
    let det0 = inverse3(jac, inv);
    for m in 3 {
        let dp = -2.0 * xi[m] * det0 / det;
        for s in 3 { g[at + 3 * m + s] = inv[3 * s + m] * dp; }
    }
}

var full: f64[>..] = [];       // the 33 × 33 matrix of the nodes and the modes of a cell
var condensed: f64[>..] = [];  // X = Kaa⁻¹ Kau: 9 rows (the mode dofs), 24 columns (the node dofs)

// Fills `full` for cell c (hex8i_full in kernels.goose), then `condensed`.
fn hex8i_condense(c: i64) {
    guard c >= 0 && c < cell_count() && cell_type[c] as i64 == HEX8I && cell_start[c + 1] - cell_start[c] == 8 else { abort("hex8i_condense: a hex8i cell with 8 nodes"); }
    hex8i_full(c);
    // Gauss elimination of Kaa X = Kau. Kaa is symmetric positive definite,
    // so it needs no pivots.
    var a: f64[81] = [0.0; 81];
    condensed.clear();
    for r in MODE_DOFS {
        for q in MODE_DOFS { a[MODE_DOFS * r + q] = full[ALL_DOFS * (NODE_DOFS + r) + NODE_DOFS + q]; }
        for j in NODE_DOFS { condensed.push(full[ALL_DOFS * (NODE_DOFS + r) + j]); }
    }
    for k in MODE_DOFS {
        for r in k + 1..MODE_DOFS {
            let f = a[MODE_DOFS * r + k] / a[MODE_DOFS * k + k];
            for q in k..MODE_DOFS { a[MODE_DOFS * r + q] -= f * a[MODE_DOFS * k + q]; }
            for j in NODE_DOFS { condensed[NODE_DOFS * r + j] -= f * condensed[NODE_DOFS * k + j]; }
        }
    }
    var k = MODE_DOFS - 1;
    while k >= 0 {
        for j in NODE_DOFS {
            var v = condensed[NODE_DOFS * k + j];
            for q in k + 1..MODE_DOFS { v -= a[MODE_DOFS * k + q] * condensed[NODE_DOFS * q + j]; }
            condensed[NODE_DOFS * k + j] = v / a[MODE_DOFS * k + k];
        }
        k--;
    }
}

// The element routine: ke = Kuu - Kua Kaa⁻¹ Kau.
fn hex8i_matrix(c: i64, ke: f64[>..]&) -> bool {
    guard c >= 0 && c < cell_count() && cell_type[c] as i64 == HEX8I else { abort("hex8i_matrix: a hex8i cell"); }
    hex8i_condense(c);
    ke.clear();
    for i in NODE_DOFS {
        for j in NODE_DOFS {
            var v = full[ALL_DOFS * i + j];
            for r in MODE_DOFS { v -= full[ALL_DOFS * i + NODE_DOFS + r] * condensed[NODE_DOFS * r + j]; }
            ke.push(v);
        }
    }
    return true;
}

var mode_grad: f64[>..] = [];

// Adds the gradient of the modes of cell c at the reference point xi, where
// the Jacobian determinant is det, to g[3 k + s] = ∂u_k/∂x_s. The mode
// amplitudes are -Kaa⁻¹ Kau u of the nodal displacements u.
fn hex8i_add_modes(c: i64, xi: f64[3], det: f64, g: f64[9]&) {
    guard dofs_per_node == 3 && disp.len == dof_total() else { abort("hex8i_add_modes: 3 dofs for each node and the solution"); }
    hex8i_condense(c);
    var alpha: f64[9] = [0.0; 9];
    for r in MODE_DOFS {
        for j in NODE_DOFS { alpha[r] -= condensed[NODE_DOFS * r + j] * disp[dofs_per_node * cell_node(c, j / 3) + j % 3]; }
    }
    mode_grad.clear();
    for i in MODE_DOFS { mode_grad.push(0.0); }
    mode_gradient(c, xi, det, mode_grad, 0);
    for m in 3 {
        for k in 3 {
            for s in 3 { g[3 * k + s] += alpha[3 * m + k] * mode_grad[3 * m + s]; }
        }
    }
}
