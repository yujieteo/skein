// The direct solve: a profile (skyline) LU factorization without
// pivoting, as COLSOL in K. J. Bathe, Finite Element Procedures, 1996,
// §8.2.3, with a lower part for a matrix that is not symmetric.
//
// Reason: the rkpm matrix has a large condition number (the Nitsche
// penalty and the wide supports), and on the plate of
// samples/05_meshless_plate.skn CG with AMG and GMRES with ILU(0) take
// thousands of iterations or do not converge. A direct solve does not
// depend on it.
//
// The points take the reverse Cuthill-McKee order of the pattern, which
// makes the profile small. Dof d of the pattern is row (and column)
// renumber(d) of the profile. Column j holds the rows top[j]..j of U, and
// row j the columns top[j]..j - 1 of L, at sky_start[j] + (i - top[j]).
import std;
import math;
namespace meshless;
import .model;
import .solve;
import solvers from "../solvers";

var position: i64[>..] = [];    // the position of each point in the order
var top: i64[>..] = [];
var sky_start: i64[>..] = [];
var upper: f64[>..] = [];
var lower: f64[>..] = [];

fn renumber(d: i64) -> i64 { return ncomp * position[d / ncomp] + d % ncomp; }

// The points that couple with point p in the pattern.
fn degree(p: i64) -> i64 { return (row_start[ncomp * p + 1] - row_start[ncomp * p]) / ncomp; }

// Appends to seq the points of the component of root, in the Cuthill-McKee
// order: by level from root, and in each level by degree. It marks each
// with stamp, and gives the last point.
fn cuthill(root: i64, stamp: i64, mark: i64[>..]&, seq: i64[>..]&) -> i64 {
    var head = seq.len;
    seq.push(root);
    mark[root] = stamp;
    var next: i64[>..] = [];
    while head < seq.len {
        let p = seq[head];
        head++;
        next.clear();
        var e = row_start[ncomp * p];
        while e < row_start[ncomp * p + 1] {
            let q = col[e] / ncomp;
            if mark[q] != stamp {
                mark[q] = stamp;
                next.push(q);
            }
            e += ncomp;
        }
        sort(next) { a, b => degree(a) < degree(b) };
        seq.append(next);
    }
    return seq[seq.len - 1];
}

// The reverse Cuthill-McKee order of the points. Each component starts
// from the last point of a first pass from its point of least degree: a
// point near the end of the component.
fn order_points() {
    let n = point_count();
    var mark: i64[>..] = [];
    for p in n { mark.push(-1); }
    var seq: i64[>..] = [];
    var trial: i64[>..] = [];
    var stamp = 0;
    while seq.len < n {
        var root = -1;
        for p in n { if mark[p] < 0 && (root < 0 || degree(p) < degree(root)) { root = p; } }
        trial.clear();
        let far = cuthill(root, n + stamp, mark, trial);
        for p in trial { mark[p] = -1; }
        cuthill(far, stamp, mark, seq);
        stamp++;
    }
    position.clear();
    for p in n { position.push(0); }
    for k in n { position[seq[k]] = n - 1 - k; }
}

// The order, the profile of the pattern and its storage.
fn size_profile() {
    order_points();
    let n = dof_count();
    top.clear();
    for j in n { top.push(j); }
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let i = renumber(r);
            let j = renumber(col[e]);
            top[max(i, j)] = min(top[max(i, j)], min(i, j));
        }
    }
    sky_start.clear();
    sky_start.push(0);
    for j in n { sky_start.push(sky_start[j] + j - top[j] + 1); }
    solvers::memory_check("meshless skyline lu", sky_start[n], n, 2, 0, "");
    solvers::zeros(upper, sky_start[n]);
    solvers::zeros(lower, sky_start[n]);
}

fn profile_size() -> i64 { return sky_start[sky_start.len - 1]; }

// The sum of a[k] b[k] for k from lo to hi - 1, with a and b at k from
// their starts.
fn span_dot(a: f64[:], a0: i64, b: f64[:], b0: i64, lo: i64, hi: i64) -> f64 {
    var s = 0.0;
    for k in lo..hi { s += a[a0 + k] * b[b0 + k]; }
    return s;
}

// A = L U in place (Crout). False if a pivot is 0 relative to its row: the
// matrix is singular, or needs a pivot order.
fn factor() -> bool {
    for &v in upper { v = 0.0; }
    for &v in lower { v = 0.0; }
    for r in dof_count() {
        for e in row_start[r]..row_start[r + 1] {
            let i = renumber(r);
            let j = renumber(col[e]);
            if i <= j { upper[sky_start[j] + i - top[j]] = val[e]; } else { lower[sky_start[i] + j - top[i]] = val[e]; }
        }
    }
    // A symmetric matrix has L(j, i) = U(i, j) / U(i, i): half the work.
    let symmetric = method == RKPM && (physics == ELASTIC || (ax == 0.0 && ay == 0.0));
    for j in dof_count() {
        // U(i, j) and L(j, i) for i < j, then U(j, j).
        let uj = sky_start[j] - top[j];
        for i in top[j]..j {
            let ui = sky_start[i] - top[i];
            let lo = max(top[i], top[j]);
            upper[uj + i] -= span_dot(lower, ui, upper, uj, lo, i);
            let l = if symmetric { upper[uj + i] } else { lower[uj + i] - span_dot(lower, uj, upper, ui, lo, i) };
            lower[uj + i] = l / upper[ui + i];
        }
        let size = abs(upper[uj + j]);
        upper[uj + j] -= span_dot(lower, uj, upper, uj, top[j], j);
        if !(abs(upper[uj + j]) > 1.0e-12 * size) { return false; }
    }
    return true;
}

// sol = A⁻¹ rhs with the factors: L y = b, then U x = y.
fn substitute() {
    let n = dof_count();
    var x: f64[>..] = [];
    solvers::zeros(x, n);
    for d in n { x[renumber(d)] = rhs[d]; }
    for j in n { x[j] -= span_dot(lower, sky_start[j] - top[j], x, 0, top[j], j); }
    var j = n - 1;
    while j >= 0 {
        let uj = sky_start[j] - top[j];
        x[j] /= upper[uj + j];
        for i in top[j]..j { x[i] -= upper[uj + i] * x[j]; }
        j--;
    }
    for d in n { sol[d] = x[renumber(d)]; }
}
