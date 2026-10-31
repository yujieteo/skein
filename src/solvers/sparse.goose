// The operations on CSR matrices and vectors that the solvers share.
import std;
import math;
namespace solvers;

// y = A x
fn multiply(row_start: i64[:], col: i32[:], val: f64[:], x: f64[:], y: f64[:]) {
    assert(row_start.len == y.len + 1 && col.len == val.len && x.len == y.len);
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

fn norm(a: f64[:]) -> f64 { return sqrt(dot(a, a)); }

// Sets `v` to n zeros.
fn zeros(v: f64[>..]&, n: i64) {
    v.clear();
    for i in n { v.push(0.0); }
}

// The position of column c in row r, or -1 if row r has no column c.
fn entry(row_start: i64[:], col: i32[:], r: i64, c: i64) -> i64 {
    var lo = row_start[r];
    var hi = row_start[r + 1] - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if col[mid] < c { lo = mid + 1; } else { hi = mid; }
    }
    if lo > hi || col[lo] != c { return -1; }
    return lo;
}

// Holds each dof d with fixed[d] != 0 at 0: its row and its column become
// those of the identity, and b[d] becomes 0. The matrix stays symmetric
// if it was, so that CG and GMRES solve the same system.
fn constrain(row_start: i64[:], col: i32[:], val: f64[:], fixed: u8[:], b: f64[:]) {
    assert(row_start.len == b.len + 1 && fixed.len == b.len && col.len == val.len);
    for r in b.len {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e];
            if fixed[r] != 0 || fixed[c] != 0 { val[e] = if c == r { 1.0 } else { 0.0 }; }
        }
        if fixed[r] != 0 { b[r] = 0.0; }
    }
}
