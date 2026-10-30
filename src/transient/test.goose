// Tests of the transient package. `./build.sh test` runs this program from
// the root of the repository, and a non-zero exit status is a failure.
//
// The acceptance case is femahoy validation/cases/grid-cooling-slab
// (examples/femahoy/grid_cooling_slab.skn): a slab 20 mm thick at 100 degC,
// k = 1 W/(m K), rho cp = 10⁶ J/(m³ K). At time 0 both faces go to 0 degC.
// The reference is the series of Carslaw and Jaeger (1959), section 3.3:
//
//     T(x, t) = Σ 4 T0 / (m π) sin(m π x / H) exp(-(m π / H)² α t),  m odd
//
// with α = k / (rho cp). Units: m, s, K.
import std;
import math;
import transient from "../transient";

var failures = 0;

let H = 0.02;
let K = 1.0;
let RHO_CP = 1.0e6;
let T0 = 100.0;
let END = 20.0;

// The value of femahoy and of the `expect` line of grid_cooling_slab.skn:
// the centre at t = 20 s (Fo = α t / (H/2)² = 0.2).
let SERIES_CENTRE = 77.231160686;

fn series(x: f64, t: f64) -> f64 {
    let alpha = K / RHO_CP;
    var s = 0.0;
    for i in 200 {
        let m = (2 * i + 1) as f64;
        let w = m * PI / H;
        s += 4.0 * T0 / (m * PI) * sin(w * x) * exp(-w * w * alpha * t);
    }
    return s;
}

fn slab(cells: i64) {
    transient::line(0.0, H, cells);
    transient::material(K, RHO_CP);
    assert(transient::fix_end(0.0, 0.0));
    assert(transient::fix_end(H, 0.0));
    transient::initial(T0);
}

// The centre temperature at END, or -1 if a solve does not converge.
fn centre(cells: i64, dt: f64, scheme: i64) -> f64 {
    slab(cells);
    guard transient::run(END, dt, scheme) else {
        failures++;
        print("FAIL run: CG does not converge, dt ", dt, ", scheme ", scheme);
        return -1.0;
    }
    return transient::temperature_at(0.5 * H);
}

// The femahoy case: backward Euler, dt = 0.1 s, within the 0.06 K of its
// `expect` line. The test uses quadratic cells of the same 0.5 mm size as
// the femahoy model, which uses linear cells.
fn test_femahoy() {
    let series_centre = series(0.5 * H, END);
    if abs(series_centre - SERIES_CENTRE) > 1.0e-8 {
        failures++;
        print("FAIL series: ", series_centre, ", want ", SERIES_CENTRE);
    }
    let T = centre(40, 0.1, transient::BACKWARD_EULER);
    print("cooling slab: 40 quadratic cells, backward Euler, dt 0.1 s, ",
        transient::step_count(), " steps, ", transient::cg_iterations(), " CG iterations");
    print("  centre T at 20 s ", T, " (series ", SERIES_CENTRE, ", error ", T - SERIES_CENTRE, ")");
    if abs(T - SERIES_CENTRE) > 0.06 {
        failures++;
        print("FAIL femahoy centre: error ", T - SERIES_CENTRE, " K, limit 0.06 K");
    }
    if abs(transient::time() - END) > 1.0e-12 {
        failures++;
        print("FAIL time: ", transient::time(), ", want ", END);
    }
}

// The error at the centre for dt = 2, 1, 0.5, 0.25 and 0.125 s on 200
// cells, where the error of the mesh is far below the error of the time
// steps. The observed order of the two smallest steps must be at most 0.3
// below the order of the scheme (docs/implementation.md, convergence).
fn test_order(name: u8[:], scheme: i64, order: f64) {
    print(name, ":");
    var last_err = 0.0;
    var observed = 0.0;
    var dt = 2.0;
    for k in 5 {
        let err = abs(centre(200, dt, scheme) - SERIES_CENTRE);
        if k == 0 {
            print("  dt ", dt, " s, error ", err, " K");
        } else {
            observed = log2(last_err / err);
            print("  dt ", dt, " s, error ", err, " K, order ", observed);
        }
        last_err = err;
        dt *= 0.5;
    }
    if observed < order - 0.3 {
        failures++;
        print("FAIL ", name, " order: ", observed, ", want ", order);
    }
}

fn main() {
    test_femahoy();
    test_order("backward Euler", transient::BACKWARD_EULER, 1.0);
    test_order("Crank-Nicolson", transient::CRANK_NICOLSON, 2.0);
    test_order("BDF2", transient::BDF2, 2.0);
    if failures > 0 { exit(1); }
    print("transient: all tests pass");
}
