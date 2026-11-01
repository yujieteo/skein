// The matrix of the package optim and its direct solve: the Cholesky
// factor L of the band of the matrix, K = L Lᵀ. Reason: the matrix of a
// grid is symmetric positive definite with a narrow band, and a topology
// loop solves it again for each design and for each adjoint, so the
// factor is used two times for each solve. The row of each fixed dof is
// that of the identity.
//
// Row i of L holds the columns i - width..i, at i (width + 1) + (i - j).
import std;
import math;
namespace optim;
import .q4;
import solvers from "../solvers";

var width = 0;               // the half band width of the matrix
var factor: f64[>..] = [];
var rhs: f64[>..] = [];      // the loads, less K times the fixed values

// The half band width of the dofs of the grid: the largest difference of
// two dofs of a cell.
fn band_width() -> i64 {
    var w = 0;
    for c in cell_count() {
        let e = cell_dofs(c);
        for a in e { for b in e { w = max(w, a - b); } }
    }
    return w;
}

fn at(i: i64, j: i64) -> i64 { return i * (width + 1) + (i - j); }

// The loads of the step, and the band of K, with the identity on the
// fixed dofs and their values in rhs.
fn assemble() {
    let n = dof_count();
    for &v in rhs { v = 0.0; }
    for &v in factor { v = 0.0; }
    for l in loads { add_load(l, 1.0, rhs); }
    for d in n { if fixed[d] != 0 { rhs[d] = 0.0; } }
    for c in cell_count() {
        let e = cell_dofs(c);
        let m = cell_modulus(c);
        for r in 8 {
            let i = e[r];
            guard fixed[i] == 0 else { continue; }
            for q in 8 {
                let j = e[q];
                if fixed[j] != 0 {
                    rhs[i] -= m * k0[8 * q + r] * fixed_value[j];
                } else if j <= i {
                    factor[at(i, j)] += m * k0[8 * q + r];
                }
            }
        }
    }
    for d in n {
        if fixed[d] != 0 {
            factor[at(d, d)] = 1.0;
            rhs[d] = fixed_value[d];
        }
    }
}

// L in place of the band. False if a pivot is not positive: the
// conditions do not hold the body.
fn cholesky() -> bool {
    let n = dof_count();
    for i in n {
        let lo = max(0, i - width);
        for j in lo..i + 1 {
            var s = factor[at(i, j)];
            for k in max(lo, j - width)..j { s -= factor[at(i, k)] * factor[at(j, k)]; }
            if j < i {
                factor[at(i, j)] = s / factor[at(j, j)];
            } else {
                guard s > 1.0e-14 * abs(factor[at(i, i)]) + 1.0e-300 else { return false; }
                factor[at(i, i)] = sqrt(s);
            }
        }
    }
    return true;
}

// x = K⁻¹ b with the factor: L y = b, then Lᵀ x = y.
fn substitute(b: f64[:], x: f64[:]) {
    let n = dof_count();
    guard b.len == n && x.len == n && factor.len == n * (width + 1) else { abort("substitute: b and x of one value for each dof, and a factor"); }
    for i in n {
        var s = b[i];
        for k in max(0, i - width)..i { s -= factor[at(i, k)] * x[k]; }
        x[i] = s / factor[at(i, i)];
    }
    var i = n - 1;
    while i >= 0 {
        x[i] /= factor[at(i, i)];
        for k in max(0, i - width)..i { x[k] -= factor[at(i, k)] * x[i]; }
        i--;
    }
}

// The facts of the last solve, for the records of the program.
var solves = 0;              // the solves of the run: 1, or one for each design of a topology loop
var last_seconds = 0.0;

// Solves the step: u = K⁻¹ f with the fixed values. False if the
// conditions do not hold the body.
fn solve() -> bool {
    guard young > 0.0 && fixed.len == dof_count() else { abort("solve: elastic and clear_bcs"); }
    let start = clock();
    width = band_width();
    let n = dof_count();
    // The band, the loads, u, and the 2 vectors of the adjoint.
    solvers::memory_check("q4 band cholesky", n * (width + 1), n, 4, 0, "");
    if u.len != n || factor.len != n * (width + 1) {
        solvers::zeros(u, n);
        solvers::zeros(rhs, n);
        solvers::zeros(factor, n * (width + 1));
    }
    let sizes = state_sizes();
    assemble();
    let regular = cholesky();
    if regular { substitute(rhs, u); } else { for &v in u { v = 0.0; } }
    guard state_sizes() == sizes else { abort("solve: a solve grows no array of the state (rule 3)"); }
    solves++;
    last_seconds += clock() - start;
    return regular;
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6). The first
// solve of a grid sizes u, rhs and the factor.
fn state_sizes() -> i64[6] { return [density.len, fixed.len, u.len, factor.len, rhs.len, loads.len]; }
