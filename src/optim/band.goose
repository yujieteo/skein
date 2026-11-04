import std;
import math;
namespace optim;
import .q4;
import solvers from "../solvers";

var width = 0;
var factor: f64[>..] = [];
var rhs: f64[>..] = [];

fn band_width() -> i64 {
    var w = 0;
    for c in cell_count() {
        let e = cell_dofs(c);
        for a in e { for b in e { w = max(w, a - b); } }
    }
    return w;
}

fn at(i: i64, j: i64) -> i64 { return i * (width + 1) + (i - j); }

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

var solves = 0;
var last_seconds = 0.0;

fn solve() -> bool {
    guard young > 0.0 && fixed.len == dof_count() else { abort("solve: elastic and clear_bcs"); }
    let start = clock();
    width = band_width();
    let n = dof_count();
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

fn state_sizes() -> i64[6] { return [density.len, fixed.len, u.len, factor.len, rhs.len, loads.len]; }
