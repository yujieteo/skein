// The line mesh of quadratic cells and the assembly of M and K for heat.
//
// Cell c has the nodes 2c, 2c + 1 and 2c + 2, at its left end, its middle
// and its right end. Each node has one dof, the temperature.
import std;
import math;
namespace transient;
import .step;
import solvers from "../solvers";

var line_lo = 0.0;
var line_hi = 0.0;
var cells = 0;
var spacing = 0.0;     // the length of a cell

var conductivity = 0.0;
var capacity = 0.0;    // rho cp

fn node_count() -> i64 { return 2 * cells + 1; }

fn build_line(x0: f64, x1: f64, n: i64) {
    assert(n > 0 && x1 > x0);
    line_lo = x0;
    line_hi = x1;
    cells = n;
    spacing = (x1 - x0) / n;
    conductivity = 0.0;
    capacity = 0.0;
    temp.clear();
    fixed.clear();
    fixed_value.clear();
    for d in node_count() {
        temp.push(0.0);
        fixed.push(0);
        fixed_value.push(0.0);
    }
    now = 0.0;
    clear_matrices();
}

fn set_initial(T0: f64) {
    for d in temp.len {
        temp[d] = if fixed[d] != 0 { fixed_value[d] } else { T0 };
    }
    now = 0.0;
}

fn mark_fixed_end(x: f64, T: f64) -> bool {
    let tol = 1.0e-9 * (line_hi - line_lo);
    var d = -1;
    if abs(x - line_lo) <= tol { d = 0; }
    if abs(x - line_hi) <= tol { d = node_count() - 1; }
    guard d >= 0 else { return false; }
    fixed[d] = 1;
    fixed_value[d] = T;
    temp[d] = T;
    return true;
}

// The exact matrices of a quadratic cell of length h, for k = 1 and
// rho cp = 1: K = 1/(3h) KE and M = h/30 ME.
let KE: f64[9] = [7.0, -8.0, 1.0, -8.0, 16.0, -8.0, 1.0, -8.0, 7.0];
let ME: f64[9] = [4.0, 2.0, -1.0, 2.0, 16.0, 2.0, -1.0, 2.0, 4.0];

// The row of node i holds the nodes of the cells that hold i, in
// increasing order: i - 2 to i + 2 for an end node of a cell, and i - 1 to
// i + 1 for a middle node.
fn assemble_line() {
    assert(conductivity > 0.0 && capacity > 0.0);
    clear_matrices();
    let n = node_count();
    row_start.push(0);
    for i in n {
        let reach = if i % 2 == 0 { 2 } else { 1 };
        for j in max(i - reach, 0)..min(i + reach, n - 1) + 1 {
            col.push(j as i32);
            mass.push(0.0);
            stiff.push(0.0);
        }
        row_start.push(col.len);
    }
    let kf = conductivity / (3.0 * spacing);
    let mf = capacity * spacing / 30.0;
    for c in cells {
        for a in 3 {
            let r = 2 * c + a;
            for bb in 3 {
                let e = solvers::entry(row_start, col, r, 2 * c + bb);
                stiff[e] += kf * KE[3 * a + bb];
                mass[e] += mf * ME[3 * a + bb];
            }
        }
    }
}

fn interpolate_T(x: f64) -> f64 {
    guard temp.len > 0 else { return 0.0; }
    let c = clamp(floor((x - line_lo) / spacing) as i64, 0, cells - 1);
    let s = 2.0 * (x - line_lo - c * spacing) / spacing - 1.0;
    let n0 = 0.5 * s * (s - 1.0);
    let n1 = 1.0 - s * s;
    let n2 = 0.5 * s * (s + 1.0);
    return n0 * temp[2 * c] + n1 * temp[2 * c + 1] + n2 * temp[2 * c + 2];
}
