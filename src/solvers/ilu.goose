// ILU(0): the incomplete LU factors of A with the sparsity of A. L has a
// unit diagonal and U has the diagonal of the factor, both in one array
// `lu` with the pattern of A (Saad, Iterative Methods for Sparse Linear
// Systems, 2nd ed., Algorithm 10.4).
import std;
namespace solvers;
import .sparse;

// Factors A into `lu`, and puts the position of each diagonal entry in
// `diag`. Returns false if a row has no diagonal entry or a pivot is 0.
fn ilu0(row_start: i64[:], col: i32[:], val: f64[:], lu: f64[>..]&, diag: i64[>..]&) -> bool {
    guard row_start.len > 0 && col.len == val.len && row_start[row_start.len - 1] == col.len else { abort("ilu0: row_start ends at col.len, and col and val have one size"); }
    let n = row_start.len - 1;
    lu.clear();
    for v in val { lu.push(v); }
    diag.clear();
    for r in n {
        let d = entry(row_start, col, r, r);
        if d < 0 { return false; }
        diag.push(d);
    }
    // where[c] is the position of column c in row r, or -1
    var where: i64[>..] = [];
    for c in n { where.push(-1); }
    for r in n {
        for e in row_start[r]..row_start[r + 1] { where[col[e]] = e; }
        for e in row_start[r]..diag[r] {
            let k = col[e] as i64;
            let pivot = lu[diag[k]];
            if pivot == 0.0 { return false; }
            lu[e] /= pivot;
            let factor = lu[e];
            for f in diag[k] + 1..row_start[k + 1] {
                let w = where[col[f]];
                if w >= 0 { lu[w] -= factor * lu[f]; }
            }
        }
        for e in row_start[r]..row_start[r + 1] { where[col[e]] = -1; }
        if lu[diag[r]] == 0.0 { return false; }
    }
    return true;
}

// z = (LU)⁻¹ r: a forward and a back substitution.
fn ilu0_apply(row_start: i64[:], col: i32[:], lu: f64[:], diag: i64[:], r: f64[:], z: f64[:]) {
    guard row_start.len == r.len + 1 && z.len == r.len && diag.len == r.len && lu.len == col.len else { abort("ilu0_apply: row_start has r.len + 1 entries, z and diag r.len, lu col.len"); }
    let n = r.len;
    for i in n {
        var s = r[i];
        for e in row_start[i]..diag[i] { s -= lu[e] * z[col[e]]; }
        z[i] = s;
    }
    for k in n {
        let i = n - 1 - k;
        var s = z[i];
        for e in diag[i] + 1..row_start[i + 1] { s -= lu[e] * z[col[e]]; }
        z[i] = s / lu[diag[i]];
    }
}
