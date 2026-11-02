// The transient step of a heat physics (spec §14.2, `transient` with
// scheme backward_euler or crank_nicolson): C ddt(T) + K T = f with the
// capacity C_ij = ∫ rho cp Ni Nj dV, in the θ-scheme
//
//     (C/dt + θ K) T' = (C/dt - (1 - θ) K) T + f
//
// θ = 1 is backward Euler and θ = 1/2 Crank-Nicolson. The fixed dofs hold
// their prescribed value from the start of the step. The step starts from
// the temperature at the end of the last solve or step, or from the initial
// value of the unknown, and runs one time step at a time, so that the
// program writes the fields each k steps.
//
// The solve of a time step is the solve of A du = f + B T - A T of
// solve.goose from the state T in `sol`, with A = C/dt + θ K and B = C/dt
// - (1 - θ) K. So the matrix is A and the load is f + B T. Both A stay the
// same in all time steps, so the AMG levels or the direct factor of the
// first time step stay.
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

var capacity = 0.0;             // rho cp of the material
var th_k: f64[>..] = [];        // K, on the pattern
var th_b: f64[>..] = [];        // B = C/dt - (1 - theta) K
var th_load: f64[>..] = [];     // the load f of the step
var th_work: f64[>..] = [];     // B T
var th_start = 0.0;             // the clock at the start of the step
var th_index = 0;               // the time steps done
var th_steps = 0;               // the time steps of the step
var th_dt = 0.0;                // the time step

fn heat_capacity(c: f64) {
    guard c > 0.0 else { abort("heat_capacity: rho cp > 0"); }
    capacity = c;
}

// The quadrature of the capacity: the rule of a cell for a hexahedron, else
// the conical product rule of the mass of eigen.goose, exact for a P2 mass.
fn capacity_count() -> i64 { return if hex_cells { quad_count(dim) } else { rich_count() }; }

fn capacity_point(q: i64) -> f64[3], f64 {
    if hex_cells { return quad_point(dim, q); }
    return rich_point(q);
}

// ∫ rho cp Ni Nj dV of all live cells, on the pattern of K, into m.
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

// Assembles K and f, the capacity and the scheme, and puts the state at the
// start of the step. It gives the number of time steps, the least that has
// a step of at most dt.
fn transient_heat_begin(end: f64, dt: f64, theta: f64) -> i64 {
    guard physics == HEAT && dofs_sized() && capacity > 0.0 && end > 0.0 && dt > 0.0 && theta >= 0.5 && theta <= 1.0 else { abort("transient_heat_begin: a heat physics with a capacity, end > 0, dt > 0 and 1/2 <= theta <= 1"); }
    th_start = clock();
    assemble("simplex cg amg", 6 + 2 * mode_count());
    hold_dead_nodes();
    th_steps = max(ceil(end / dt - 1.0e-9) as! i64, 1);
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

// One time step. False if its solve does not converge.
fn transient_heat_advance() -> bool {
    guard th_index < th_steps && th_work.len == sol.len else { abort("transient_heat_advance: a step started by transient_heat_begin, and a step left"); }
    solvers::multiply(row_start, col, th_b, sol, th_work);
    for d in sol.len { load[d] = th_load[d] + th_work[d]; }
    let done = linear_solve(th_start, if th_index == 0 { solvers::warm_key() } else { -1 }, false);
    th_index++;
    return done;
}

// Puts the matrix K and the load f back, for the reactions and the fields,
// and keeps the temperature as the start of the next step.
fn transient_heat_end() {
    for e in val.len { val[e] = th_k[e]; }
    for d in load.len { load[d] = th_load[d]; }
    keep_factor = false;
    factor_ready = false;
    state_valid = true;
    last_seconds = clock() - th_start;
}
