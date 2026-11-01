// The direct solve of a Newton iteration of the nonlinear solid
// (solid.goose): a skyline (profile) LDLᵀ factorization, as COLSOL in K. J.
// Bathe, Finite Element Procedures, 1996, §8.2.3, in the reverse
// Cuthill-McKee order of the nodes, as src/meshless/direct.goose.
//
// Reason: K of a nearly incompressible solid has a large condition number.
// On examples/femahoy/hyperelastic_tube.skn (K/μ = 9654, 2640 dofs), CG
// with AMG takes about 500 iterations for each Newton solve, and the model
// takes 34 s on 1 thread. With this solve it takes 7 s.
//
// A static step uses it when the factorization takes at most DIRECT_WORK
// multiplies, else CG with AMG. A fixed dof is a row and a column of the
// identity, with a right side of 0. Node n is at position
// sky_position[n] of the order, and dof d at row sky_row(d). Column j
// holds the rows sky_top[j]..j of U = D Lᵀ, and row j the columns
// sky_top[j]..j - 1 of L, at sky_start[j] + (i - sky_top[j]).
import std;
import math;
namespace elements;
import .mesh;
import .assemble;
import .solve;
import solvers from "../solvers";

let DIRECT_WORK = 5.0e8;    // about 0.5 s of multiplies

var sky_position: i64[>..] = [];
var sky_top: i64[>..] = [];
var sky_start: i64[>..] = [];
var sky_upper: f64[>..] = [];
var sky_lower: f64[>..] = [];

fn sky_row(d: i64) -> i64 { return ncomp * sky_position[d / ncomp] + d % ncomp; }

// The entries of the row of the first dof of node p: its degree.
fn node_degree(p: i64) -> i64 { return row_start[ncomp * p + 1] - row_start[ncomp * p]; }

// Appends to seq the nodes of the component of root, in the Cuthill-McKee
// order: by level from root, and in each level by degree. It marks each
// with stamp, and gives the last node.
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

// The reverse Cuthill-McKee order of the nodes. Each component starts from
// the last node of a first pass from its node of least degree.
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

// The order and the skyline of the pattern of the step. False if the
// factorization takes more than DIRECT_WORK multiplies: the step then
// solves with CG and AMG.
fn direct_setup() -> bool {
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
    guard work <= DIRECT_WORK else { return false; }
    sky_start.clear();
    sky_start.push(0);
    for j in n { sky_start.push(sky_start[j] + j - sky_top[j] + 1); }
    solvers::memory_check("simplex newton skyline ldlt", sky_start[n], n, 2, 0, "");
    solvers::zeros(sky_upper, sky_start[n]);
    solvers::zeros(sky_lower, sky_start[n]);
    return true;
}

// The sum of a[a0 + k] b[b0 + k] for k from lo to hi - 1.
fn sky_dot(a: f64[:], a0: i64, b: f64[:], b0: i64, lo: i64, hi: i64) -> f64 {
    var s = 0.0;
    for k in lo..hi { s += a[a0 + k] * b[b0 + k]; }
    return s;
}

// K = L U in place, with the fixed dofs as rows of the identity, without
// pivoting. A negative pivot is correct: the tangent of the first Newton
// iteration of a step is indefinite when the fixed dofs move (the tube has
// a pivot of -83 on a diagonal of 8299). False if a pivot is near 0
// relative to its diagonal.
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

// The Newton solve of solve_solid: sol of the free dofs from
// K sol = load, with sol of the fixed dofs as it is, as cg_solve. False
// if a pivot is near 0.
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
