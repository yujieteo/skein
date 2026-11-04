import std;
import math;
namespace elements;
import .mesh;
import .assemble;
import .solve;
import solvers from "../solvers";

let DIRECT_WORK = 5.0e8;

var sky_position: i64[>..] = [];
var sky_top: i64[>..] = [];
var sky_start: i64[>..] = [];
var sky_upper: f64[>..] = [];
var sky_lower: f64[>..] = [];

fn sky_row(d: i64) -> i64 { return ncomp * sky_position[d / ncomp] + d % ncomp; }

fn node_degree(p: i64) -> i64 { return row_start[ncomp * p + 1] - row_start[ncomp * p]; }

fn cuthill_nodes(root: i64, stamp: i64, mark: i64[>..]&, seq: i64[>..]&) -> i64 {
    var head = seq.len;
    seq.push(root);
    mark[root] = stamp;
    var next: i64[>..] = [];
    while head < seq.len {
        let p = seq[head];
        head++;
        next.clear();
        for e in row_start[ncomp * p]..row_start[ncomp * p + 1] {
            let q = col[e] / ncomp;
            if mark[q] != stamp {
                mark[q] = stamp;
                next.push(q);
            }
        }
        sort(next) { a, b => node_degree(a) < node_degree(b) };
        seq.append(next);
    }
    return seq[seq.len - 1];
}

fn order_nodes() {
    let n = node_count();
    var mark: i64[>..] = [];
    for p in n { mark.push(-1); }
    var seq: i64[>..] = [];
    var trial: i64[>..] = [];
    var stamp = 0;
    while seq.len < n {
        var root = -1;
        for p in n { if mark[p] < 0 && (root < 0 || node_degree(p) < node_degree(root)) { root = p; } }
        trial.clear();
        let far = cuthill_nodes(root, n + stamp, mark, trial);
        for p in trial { mark[p] = -1; }
        cuthill_nodes(far, stamp, mark, seq);
        stamp++;
    }
    sky_position.clear();
    for p in n { sky_position.push(0); }
    for k in n { sky_position[seq[k]] = n - 1 - k; }
}

fn direct_setup() -> bool {
    guard skyline_layout(DIRECT_WORK) else { return false; }
    let n = sol.len;
    solvers::memory_check("simplex newton skyline ldlt", sky_start[n], n, 2, 0, "");
    solvers::zeros(sky_upper, sky_start[n]);
    solvers::zeros(sky_lower, sky_start[n]);
    return true;
}

fn skyline_layout(limit: f64) -> bool {
    order_nodes();
    let n = sol.len;
    sky_top.clear();
    for j in n { sky_top.push(j); }
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let i = sky_row(r);
            let j = sky_row(col[e]);
            sky_top[max(i, j)] = min(sky_top[max(i, j)], min(i, j));
        }
    }
    var work = 0.0;
    for j in n { work += 0.5 * (j - sky_top[j]) * (j - sky_top[j]); }
    guard work <= limit else { return false; }
    sky_start.clear();
    sky_start.push(0);
    for j in n { sky_start.push(sky_start[j] + j - sky_top[j] + 1); }
    return true;
}

fn sky_dot(a: f64[:], a0: i64, b: f64[:], b0: i64, lo: i64, hi: i64) -> f64 {
    var s = 0.0;
    for k in lo..hi { s += a[a0 + k] * b[b0 + k]; }
    return s;
}

fn sky_factor() -> bool {
    for &v in sky_upper { v = 0.0; }
    for r in sol.len {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            guard c == r || fixed[r] == 0 && fixed[c] == 0 else { continue; };
            let i = sky_row(r);
            let j = sky_row(c);
            if i <= j { sky_upper[sky_start[j] + i - sky_top[j]] = if fixed[r] != 0 { 1.0 } else { val[e] }; }
        }
    }
    for j in sol.len {
        let uj = sky_start[j] - sky_top[j];
        for i in sky_top[j]..j {
            let ui = sky_start[i] - sky_top[i];
            sky_upper[uj + i] -= sky_dot(sky_lower, ui, sky_upper, uj, max(sky_top[i], sky_top[j]), i);
            sky_lower[uj + i] = sky_upper[uj + i] / sky_upper[ui + i];
        }
        let size = sky_upper[uj + j];
        sky_upper[uj + j] -= sky_dot(sky_lower, uj, sky_upper, uj, sky_top[j], j);
        if !(abs(sky_upper[uj + j]) > 1.0e-12 * abs(size)) { return false; }
    }
    return true;
}

fn direct_solve(start: f64) -> bool {
    let n = sol.len;
    guard sky_start.len == n + 1 else { abort("direct_solve: direct_setup of the pattern"); }
    let sizes = state_sizes();
    var x: f64[>..] = [];
    solvers::zeros(x, n);
    for r in n {
        guard fixed[r] == 0 else { continue; };
        var ku = 0.0;
        for e in row_start[r]..row_start[r + 1] { ku += val[e] * sol[col[e]]; }
        x[sky_row(r)] = load[r] + assembled_source[r] - ku;
    }
    let done = sky_factor();
    if done {
        for j in n { x[j] -= sky_dot(sky_lower, sky_start[j] - sky_top[j], x, 0, sky_top[j], j); }
        var j = n - 1;
        while j >= 0 {
            let uj = sky_start[j] - sky_top[j];
            x[j] /= sky_upper[uj + j];
            for i in sky_top[j]..j { x[i] -= sky_upper[uj + i] * x[j]; }
            j--;
        }
        for d in n { if fixed[d] == 0 { sol[d] += x[sky_row(d)]; } }
    }
    guard state_sizes() == sizes else { abort("direct_solve: a solve grows no array of the state (rule 3)"); }
    last_iterations = 1;
    last_seconds = clock() - start;
    return done;
}
