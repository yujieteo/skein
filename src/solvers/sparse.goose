// The operations on CSR matrices and vectors that the solvers share.
import std;
import math;
namespace solvers;

// y = A x
fn multiply(row_start: i64[:], col: i32[:], val: f64[:], x: f64[:], y: f64[:]) {
    guard row_start.len == y.len + 1 && col.len == val.len && x.len == y.len else { abort("multiply: row_start has y.len + 1 entries, col and val one size, x and y one size"); }
    multiply_rows(row_start, col, val, x, y);
}

// y = the rows of A in row_start, col and val, times x. row_start[0] is
// the start of the first of these rows: col and val hold only these rows.
// The workers of pool.goose and multiply share it, so a row is the same
// sum on each number of threads.
fn multiply_rows(row_start: i64[:], col: i32[:], val: f64[:], x: f64[:], y: f64[:]) {
    guard row_start.len == y.len + 1 && col.len == val.len && row_start[y.len] - row_start[0] == col.len else { abort("multiply_rows: row_start has y.len + 1 entries, and col and val hold these rows"); }
    let base = row_start[0];
    for r in y.len {
        var s = 0.0;
        for e in row_start[r] - base..row_start[r + 1] - base { s += val[e] * x[col[e]]; }
        y[r] = s;
    }
}

fn dot(a: f64[:], b: f64[:]) -> f64 {
    guard a.len == b.len else { abort("dot: a and b have one size"); }
    var s = 0.0;
    for v, i in a { s += v * b[i]; }
    return s;
}

fn norm(a: f64[:]) -> f64 { return sqrt(dot(a, a)); }

// Sets `v` to n zeros.
fn zeros(v: f64[>..]&, n: i64) {
    guard n >= 0 else { abort("zeros: n >= 0"); }
    v.clear();
    for i in n { v.push(0.0); }
}

// The position of column c in row r, or -1 if row r has no column c.
fn entry(row_start: i64[:], col: i32[:], r: i64, c: i64) -> i64 {
    guard r >= 0 && r + 1 < row_start.len && row_start[row_start.len - 1] == col.len else { abort("entry: 0 <= r < rows, and row_start ends at col.len"); }
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
    guard row_start.len == b.len + 1 && fixed.len == b.len && col.len == val.len else { abort("constrain: row_start has b.len + 1 entries, fixed has b.len, col and val one size"); }
    for r in b.len {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e];
            if fixed[r] != 0 || fixed[c] != 0 { val[e] = if c == r { 1.0 } else { 0.0 }; }
        }
        if fixed[r] != 0 { b[r] = 0.0; }
    }
}
