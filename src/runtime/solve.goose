// Linear static elasticity: the cell stiffness, the sparse matrix (CSR in
// three arrays), CG with a Jacobi preconditioner and the reactions.
//
// The matrix holds all dofs, with the fixed dofs too, so that the
// reactions are K u - f on the fixed dofs. CG solves the system of the
// free dofs: it holds the residual and the search direction at 0 on the
// fixed dofs.
import std;
import math;
import os;
namespace runtime;
import .mesh;

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];

var last_iterations = 0;
var last_seconds = 0.0;

let CG_TOLERANCE = 1.0e-10;   // the residual relative to the load

fn clear_matrix() {
    row_start.clear();
    col.clear();
    val.clear();
}

// --- the sparsity pattern ------------------------------------------------------

// The lattice range of the nodes that share a cell with lattice node i,
// along an axis of n cells.
fn neighbour_range(i: i64, n: i64) -> i64, i64 {
    let first = if i % 2 == 0 { max(i / 2 - 1, 0) } else { (i - 1) / 2 };
    let last = if i % 2 == 0 { min(i / 2, n - 1) } else { (i - 1) / 2 };
    return 2 * first, 2 * last + 2;
}

// Two nodes share a cell if they share one along each axis. The 3 rows of
// a node have the same columns, in increasing order.
fn build_pattern() {
    clear_matrix();
    row_start.push(0);
    for k in lattice[2] {
        let k0, k1 = neighbour_range(k, cells[2]);
        for j in lattice[1] {
            let j0, j1 = neighbour_range(j, cells[1]);
            for i in lattice[0] {
                let i0, i1 = neighbour_range(i, cells[0]);
                for c in 3 {
                    for mk in k0..k1 + 1 {
                        for mj in j0..j1 + 1 {
                            for mi in i0..i1 + 1 {
                                let m = node_at(mi, mj, mk);
                                for q in 3 {
                                    col.push((3 * m + q) as i32);
                                    val.push(0.0);
                                }
                            }
                        }
                    }
                    row_start.push(col.len);
                }
            }
        }
    }
}

// The position of column c in row r.
fn entry(r: i64, c: i64) -> i64 {
    var lo = row_start[r];
    var hi = row_start[r + 1] - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if col[mid] < c { lo = mid + 1; } else { hi = mid; }
    }
    assert(col[lo] == c);
    return lo;
}

// --- the cell stiffness ------------------------------------------------------------

var ke: f64[6561];      // 81 × 81, row 3 l + p for component p of local node l
var grad: f64[81];      // the gradient of each shape function at a Gauss point
var dref: f64[81];      // the same in the reference cell

// The stiffness of cell (cx, cy, cz) by 3 × 3 × 3 Gauss points:
// K(i p, j q) = ∫ λ ∂p Ni ∂q Nj + μ (∂q Ni ∂p Nj + δpq ∇Ni · ∇Nj) dV.
fn cell_stiffness(cx: i64, cy: i64, cz: i64) {
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    for &v in ke { v = 0.0; }
    for gc in 3 {
        for gb in 3 {
            for ga in 3 {
                let s: f64[3] = [GAUSS_X[ga], GAUSS_X[gb], GAUSS_X[gc]];
                let w = GAUSS_W[ga] * GAUSS_W[gb] * GAUSS_W[gc];
                // jac[3 r + t] = ∂x_t / ∂s_r
                var jac: f64[9] = [0.0; 9];
                for l in 27 {
                    let a = l % 3;
                    let b = l / 3 % 3;
                    let c = l / 9;
                    dref[3 * l] = shape_1d_ds(a, s[0]) * shape_1d(b, s[1]) * shape_1d(c, s[2]);
                    dref[3 * l + 1] = shape_1d(a, s[0]) * shape_1d_ds(b, s[1]) * shape_1d(c, s[2]);
                    dref[3 * l + 2] = shape_1d(a, s[0]) * shape_1d(b, s[1]) * shape_1d_ds(c, s[2]);
                    let n = cell_node(cx, cy, cz, l);
                    for r in 3 {
                        for t in 3 { jac[3 * r + t] += dref[3 * l + r] * coords[3 * n + t]; }
                    }
                }
                let det = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7])
                    - jac[1] * (jac[3] * jac[8] - jac[5] * jac[6])
                    + jac[2] * (jac[3] * jac[7] - jac[4] * jac[6]);
                assert(det > 0.0);
                // inv = jac⁻¹, so ∂N/∂x_t = Σr inv[3 t + r] ∂N/∂s_r
                var inv: f64[9] = [0.0; 9];
                inv[0] = (jac[4] * jac[8] - jac[5] * jac[7]) / det;
                inv[1] = (jac[2] * jac[7] - jac[1] * jac[8]) / det;
                inv[2] = (jac[1] * jac[5] - jac[2] * jac[4]) / det;
                inv[3] = (jac[5] * jac[6] - jac[3] * jac[8]) / det;
                inv[4] = (jac[0] * jac[8] - jac[2] * jac[6]) / det;
                inv[5] = (jac[2] * jac[3] - jac[0] * jac[5]) / det;
                inv[6] = (jac[3] * jac[7] - jac[4] * jac[6]) / det;
                inv[7] = (jac[1] * jac[6] - jac[0] * jac[7]) / det;
                inv[8] = (jac[0] * jac[4] - jac[1] * jac[3]) / det;
                for l in 27 {
                    for t in 3 {
                        grad[3 * l + t] = inv[3 * t] * dref[3 * l] + inv[3 * t + 1] * dref[3 * l + 1]
                            + inv[3 * t + 2] * dref[3 * l + 2];
                    }
                }
                let dv = w * det;
                for i in 27 {
                    for j in 27 {
                        let dot = grad[3 * i] * grad[3 * j] + grad[3 * i + 1] * grad[3 * j + 1]
                            + grad[3 * i + 2] * grad[3 * j + 2];
                        for p in 3 {
                            let row = 81 * (3 * i + p) + 3 * j;
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
        }
    }
}

fn assemble() {
    build_pattern();
    for cz in cells[2] {
        for cy in cells[1] {
            for cx in cells[0] {
                cell_stiffness(cx, cy, cz);
                for i in 27 {
                    let ni = cell_node(cx, cy, cz, i);
                    for j in 27 {
                        let nj = cell_node(cx, cy, cz, j);
                        for p in 3 {
                            let e = entry(3 * ni + p, 3 * nj);
                            for q in 3 { val[e + q] += ke[81 * (3 * i + p) + 3 * j + q]; }
                        }
                    }
                }
            }
        }
    }
}

// --- CG ------------------------------------------------------------------------

// y = K x
fn multiply(x: f64[:], y: f64[:]) {
    for r in y.len {
        var s = 0.0;
        for e in row_start[r]..row_start[r + 1] { s += val[e] * x[col[e]]; }
        y[r] = s;
    }
}

fn dot(a: f64[:], b: f64[:]) -> f64 {
    var s = 0.0;
    for v, i in a { s += v * b[i]; }
    return s;
}

fn zeros(n: i64) -> f64[>..] {
    var v: f64[>..] = [];
    for i in n { v.push(0.0); }
    return v;
}

fn solve_linear_static() -> bool {
    assert(young > 0.0);
    let start = clock();
    assemble();
    let n = disp.len;
    var r = zeros(n);
    var z = zeros(n);
    var p = zeros(n);
    var q = zeros(n);
    var inv_diag = zeros(n);
    for d in n {
        disp[d] = 0.0;
        if fixed[d] == 0 {
            r[d] = load[d];
            inv_diag[d] = 1.0 / val[entry(d, d)];
        }
    }
    let target = CG_TOLERANCE * sqrt(dot(r, r));
    for d in n { z[d] = inv_diag[d] * r[d]; p[d] = z[d]; }
    var rz = dot(r, z);
    var converged = target == 0.0;
    var it = 0;
    while !converged && it < 10 * n {
        multiply(p, q);
        for d in n { if fixed[d] != 0 { q[d] = 0.0; } }
        let alpha = rz / dot(p, q);
        for d in n {
            disp[d] += alpha * p[d];
            r[d] -= alpha * q[d];
        }
        it++;
        if sqrt(dot(r, r)) <= target { converged = true; break; }
        for d in n { z[d] = inv_diag[d] * r[d]; }
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    last_iterations = it;
    last_seconds = clock() - start;
    return converged;
}

// --- reactions -------------------------------------------------------------------

// The force of the supports on the body: K u - f on the fixed dofs of the
// nodes of the face. With equilibrium, the reactions and the loads sum to 0.
fn face_reaction(axis: i64, value: f64) -> f64, f64, f64 {
    let f = face_index(axis, value);
    guard f >= 0 && row_start.len > 0 else { return 0.0, 0.0, 0.0; }
    let a, b = other_axes(axis);
    var sum: f64[3] = [0.0, 0.0, 0.0];
    for t in lattice[b] {
        for s in lattice[a] {
            let node = face_node(axis, f, s, t);
            for c in 3 {
                let d = 3 * node + c;
                if fixed[d] == 0 { continue; }
                var ku = 0.0;
                for e in row_start[d]..row_start[d + 1] { ku += val[e] * disp[col[e]]; }
                sum[c] += ku - load[d];
            }
        }
    }
    return sum[0], sum[1], sum[2];
}
