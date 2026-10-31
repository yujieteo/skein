// The matrix and the solve: one direct solve (direct.goose). The matrix of
// gfd, and of rkpm with advection, is not symmetric.
import std;
import math;
import os;
namespace meshless;
import solvers from "../solvers";
import .cloud;
import .model;
import .rkpm;
import .gfd;
import .direct;

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];
var rhs: f64[>..] = [];

// The facts of the last solve, for the records of the program.
var iterations = 0;            // 1: each solve is linear
var last_seconds = 0.0;
var first_residual = 0.0;      // the norm of the right side
var last_residual = 0.0;       // the norm of rhs - A u

fn start_matrix() {
    row_start.clear();
    col.clear();
    row_start.push(0);
}

// Adds the ncomp rows of a point that couples with the points `cols`, in
// increasing order.
fn add_rows(cols: i64[:]) {
    for p in ncomp {
        for k in cols {
            for q in ncomp { col.push((ncomp * k + q) as i32); }
        }
        row_start.push(col.len);
    }
}

fn size_matrix() {
    guard row_start.len == dof_count() + 1 else { abort("size_matrix: a row for each dof"); }
    solvers::zeros(val, col.len);
    solvers::zeros(rhs, dof_count());
    size_profile();
}

// The position of column c in row r.
fn entry(r: i64, c: i64) -> i64 { return solvers::entry(row_start, col, r, c); }

// Solves the static step with the source SRC(x, y, t) and the fixed values
// BC(k, x, y, t). False if the matrix is singular (factor).
fn solve_static<SRC, BC>() -> bool {
    guard method >= 0 && sol.len == dof_count() && val.len == col.len && upper.len == profile_size() else { abort("solve_static: a method and clear_bcs"); }
    let start = clock();
    let sizes = state_sizes();
    // size_profile checked the factors, which the RSS now holds.
    solvers::memory_check("meshless skyline lu", 0, dof_count(), 2, 0, "");
    if method == RKPM { rkpm_assemble(SRC, BC); } else { gfd_assemble(SRC, BC); }
    let regular = factor();
    if regular { substitute(); } else { for &v in sol { v = 0.0; } }
    iterations = 1;
    first_residual = solvers::norm(rhs);
    var r: f64[>..] = [];
    solvers::zeros(r, sol.len);
    solvers::multiply(row_start, col, val, sol, r);
    for d in r.len { r[d] = rhs[d] - r[d]; }
    last_residual = solvers::norm(r);
    guard state_sizes() == sizes else { abort("solve_static: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - start;
    return regular;
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[12] {
    return [xy.len, tri.len, edge_fix.len, edge_load.len, sol.len, row_start.len, col.len, val.len, rhs.len, upper.len,
            lower.len, item.len];
}
