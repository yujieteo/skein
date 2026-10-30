// The API of the transient heat solver. A program calls these functions in
// this sequence:
//
//     transient::line(0.0, 0.02, 100);
//     transient::material(1.0, 1.0e6);          // k, rho cp
//     transient::initial(100.0);
//     transient::fix_end(0.0, 0.0);
//     transient::fix_end(0.02, 0.0);
//     let ok = transient::run(20.0, 0.1, transient::BACKWARD_EULER);
//     let T = transient::temperature_at(0.01);
//
// The package has no units. The caller converts each value to one
// consistent system, for example SI, before the call.
//
// The mesh is a line of quadratic cells (3 nodes each), for heat that flows
// along one axis of a slab or a bar. The stepper of step.goose takes the
// matrices M and K and does not know the mesh, so a 3D mesh can give it
// its own M and K.
//
// The state is global: one mesh, one material, one set of fixed ends and
// one temperature field. `line` clears all of the state.
import std;
namespace transient;
import .line;
import .step;

let BACKWARD_EULER = 0;
let CRANK_NICOLSON = 1;
let BDF2 = 2;

// Meshes the line from x0 to x1 with n quadratic cells. Clears the
// material, the fixed ends and the temperature.
fn line(x0: f64, x1: f64, n: i64) {
    build_line(x0, x1, n);
}

// The conductivity k and the heat capacity per volume rho cp of all cells.
fn material(k: f64, rho_cp: f64) {
    conductivity = k;
    capacity = rho_cp;
}

// Sets the temperature of all nodes to T0, and the time to 0.
fn initial(T0: f64) {
    set_initial(T0);
}

// Holds the temperature of the end at x at T from time 0. Returns false if
// x is not an end of the line.
fn fix_end(x: f64, T: f64) -> bool {
    return mark_fixed_end(x, T);
}

// Steps from the current time to `end` with steps of `dt` and the scheme
// BACKWARD_EULER, CRANK_NICOLSON or BDF2. If `dt` does not divide the
// interval, the steps are made a little shorter, all of the same size.
// Returns false if a linear solve does not converge.
fn run(end: f64, dt: f64, scheme: i64) -> bool {
    assemble_line();
    return advance(end, dt, scheme);
}

// The temperature at x, interpolated in the cell that holds x.
fn temperature_at(x: f64) -> f64 {
    return interpolate_T(x);
}

// Facts about the last run, for the report.
fn time() -> f64 { return now; }
fn step_count() -> i64 { return last_steps; }
fn cg_iterations() -> i64 { return last_iterations; }
