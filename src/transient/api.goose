import std;
namespace transient;
import .line;
import .step;

let BACKWARD_EULER = 0;
let CRANK_NICOLSON = 1;
let BDF2 = 2;

fn line(x0: f64, x1: f64, n: i64) {
    build_line(x0, x1, n);
}

fn material(k: f64, rho_cp: f64) {
    conductivity = k;
    capacity = rho_cp;
}

fn initial(T0: f64) {
    set_initial(T0);
}

fn fix_end(x: f64, T: f64) -> bool {
    return mark_fixed_end(x, T);
}

fn run(end: f64, dt: f64, scheme: i64) -> bool {
    assemble_line();
    return advance(end, dt, scheme);
}

fn temperature_at(x: f64) -> f64 {
    return interpolate_T(x);
}

fn time() -> f64 { return now; }
fn step_count() -> i64 { return last_steps; }
fn cg_iterations() -> i64 { return last_iterations; }
