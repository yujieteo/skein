import std;
import math;
import os;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .solve;
import .eigen;
import .stages;
import solvers from "../solvers";

let STEP_COUNT_SLACK = 1.0e-9;

var capacity = 0.0;
var th_k: f64[>..] = [];
var th_b: f64[>..] = [];
var th_load: f64[>..] = [];
var th_work: f64[>..] = [];
var th_start = 0.0;
var th_index = 0;
var th_steps = 0;
var th_dt = 0.0;

fn heat_capacity(c: f64) {
    guard c > 0.0 else { abort("heat_capacity: rho cp > 0"); }
    capacity = c;
}

fn capacity_count() -> i64 { return if hex_cells { quad_count(dim) } else { rich_count() }; }

fn capacity_point(q: i64) -> f64[3], f64 {
    if hex_cells { return quad_point(dim, q); }
    return rich_point(q);
}

fn capacity_matrix(m: f64[:]) {
    guard m.len == val.len else { abort("capacity_matrix: the entries of K"); }
    var me: f64[100] = [0.0; 100];
    for c in cell_count() {
        if !is_live(c) { continue; }
        for &v in me { v = 0.0; }
        for q in capacity_count() {
            let xi, w = capacity_point(q);
            let dv = w * cell_geometry(c, xi) * thickness;
            for i in per_cell {
                for j in per_cell { me[per_cell * i + j] += capacity * shp[i] * shp[j] * dv; }
            }
        }
        for i in per_cell {
            for j in per_cell { m[entry(cell_node(c, i), cell_node(c, j))] += me[per_cell * i + j]; }
        }
    }
}

fn transient_heat_begin(end: f64, dt: f64, theta: f64) -> i64 {
    guard physics == HEAT && dofs_sized() && capacity > 0.0 && end > 0.0 && dt > 0.0 && theta >= 0.5 && theta <= 1.0 && curve_k.len == 0 && amb_cell.len == 0 else { abort("transient_heat_begin: a heat physics with a capacity, end > 0, dt > 0 and 1/2 <= theta <= 1, k constant, and no radiation to the surroundings"); }
    th_start = clock();
    assemble("simplex cg amg", 6 + 2 * mode_count());
    hold_dead_nodes();
    th_steps = max(ceil(end / dt - STEP_COUNT_SLACK) as! i64, 1);
    th_dt = end / th_steps;
    let n = sol.len;
    let from_last = state_valid && held.len == n;
    for d in n { sol[d] = if fixed[d] != 0 { prescribed[d] } else if from_last { held[d] } else { initial_value }; }
    th_k.clear();
    th_k.append(val);
    th_load.clear();
    th_load.append(load);
    var m = zeros(val.len);
    capacity_matrix(m);
    th_b = zeros(val.len);
    for e in val.len {
        th_b[e] = m[e] / th_dt - (1.0 - theta) * th_k[e];
        val[e] = m[e] / th_dt + theta * th_k[e];
    }
    th_work = zeros(n);
    th_index = 0;
    keep_factor = true;
    factor_ready = false;
    return th_steps;
}

fn transient_heat_advance() -> bool {
    guard th_index < th_steps && th_work.len == sol.len else { abort("transient_heat_advance: a step started by transient_heat_begin, and a step left"); }
    solvers::multiply(row_start, col, th_b, sol, th_work);
    for d in sol.len { load[d] = th_load[d] + th_work[d]; }
    let done = linear_solve(th_start, if th_index == 0 { solvers::warm_key() } else { -1 }, false);
    th_index++;
    return done;
}

fn transient_heat_end() {
    for e in val.len { val[e] = th_k[e]; }
    for d in load.len { load[d] = th_load[d]; }
    keep_factor = false;
    factor_ready = false;
    state_valid = true;
    last_seconds = clock() - th_start;
}
