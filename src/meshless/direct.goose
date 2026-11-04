import std;
import math;
namespace meshless;
import .model;
import .solve;
import solvers from "../solvers";

var position: i64[>..] = [];
var top: i64[>..] = [];
var sky_start: i64[>..] = [];
var upper: f64[>..] = [];
var lower: f64[>..] = [];

fn renumber(d: i64) -> i64 { return ncomp * position[d / ncomp] + d % ncomp; }

fn degree(p: i64) -> i64 { return (row_start[ncomp * p + 1] - row_start[ncomp * p]) / ncomp; }

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

fn span_dot(a: f64[:], a0: i64, b: f64[:], b0: i64, lo: i64, hi: i64) -> f64 {
    var s = 0.0;
    for k in lo..hi { s += a[a0 + k] * b[b0 + k]; }
    return s;
}

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
    let symmetric = method == RKPM && (physics == ELASTIC || (ax == 0.0 && ay == 0.0));
    for j in dof_count() {
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
