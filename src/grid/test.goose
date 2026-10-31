// Tests of the fv and dg methods. `./build.sh test` runs this program from
// the root of the repository, and a non-zero exit status is a failure.
//
// Each manufactured solution is solved on 3 grids, each with half the
// cell size of the one before, and the observed order of the L2 error on
// the two finest grids must be within 0.3 of the order of the method
// (docs/checks.md §4.3): 2 for fv with `linear_upwind` or with no
// advection, 1 for fv with `upwind`, and p + 1 for dg of degree p. The
// models of examples/grid/ check the same orders through
// `skein verify convergence`.
import std;
import math;
import grid from "../grid";
import solvers from "../solvers";   // the grid imports it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

// The smooth field exp(x + y/2), and its source for a = (1, 0.5), d.
let AX = 1.0;
let AY = 0.5;
fn smooth(x: f64, y: f64) -> f64 { return exp(x + 0.5 * y); }
fn smooth_source(d: f64, x: f64, y: f64) -> f64 { return (AX + 0.5 * AY - d * 1.25) * smooth(x, y); }

// sin(πx) sin(πy), and its source for -d div(grad(u)) = f.
fn bump(x: f64, y: f64) -> f64 { return sin(PI * x) * sin(PI * y); }

// The terms of the problem that the functions below give to the solves:
// the factor of a, and d.
var a_now = 1.0;
var d_now = 0.0;

fn bump_source(x: f64, y: f64, t: f64) -> f64 { return 2.0 * PI * PI * d_now * bump(x, y); }
fn zero(k: i64, x: f64, y: f64, t: f64) -> f64 { return 0.0; }
fn zero_source(x: f64, y: f64, t: f64) -> f64 { return 0.0; }
fn smooth_static_source(x: f64, y: f64, t: f64) -> f64 { return smooth_source(d_now, x, y) + (a_now - 1.0) * (AX + 0.5 * AY) * smooth(x, y); }
fn smooth_value(k: i64, x: f64, y: f64, t: f64) -> f64 { return smooth(x, y); }

// exp(x - t + y/2), the transient manufactured solution.
fn moving(x: f64, y: f64, t: f64) -> f64 { return smooth(x, y) * exp(-t); }
fn moving_source(x: f64, y: f64, t: f64) -> f64 { return (smooth_source(d_now, x, y) - smooth(x, y)) * exp(-t); }
fn moving_value(k: i64, x: f64, y: f64, t: f64) -> f64 { return moving(x, y, t); }

// Solves the static problem on n × n cells of the unit square, with the
// method of `setup`, and gives the L2 error against `smooth` or `bump`.
fn static_error(n: i64, kind: i64, a: f64, d: f64, use_bump: bool, flux: i64, limiter: i64, gradient: i64, p: i64) -> f64 {
    grid::mesh_rect(0.0, 0.0, 1.0, 1.0, n, n);
    grid::equation(1.0, a * AX, a * AY, d);
    if kind == grid::FV { grid::fv(flux, gradient, limiter); } else { grid::dg(p, flux, 2.0); }
    grid::clear_bcs();
    assert(grid::select_faces() { x, y => true });
    grid::fix_selected(0);
    a_now = a;
    d_now = d;
    let converged = if use_bump { grid::solve_static(1.0e-10, bump_source, zero) } else { grid::solve_static(1.0e-10, smooth_static_source, smooth_value) };
    if !converged { fail(str("the static solve on ", n, " x ", n, " cells does not converge")); }
    if use_bump {
        let e, norm = grid::l2_error() { x, y, t => bump(x, y) };
        return e / norm;
    }
    let e, norm = grid::l2_error() { x, y, t => smooth(x, y) };
    return e / norm;
}

// The order on the two finest of the grids of 8, 16 and 32 cells, against
// `want`.
fn check_order(what: u8[:], want: f64, kind: i64, a: f64, d: f64, use_bump: bool, flux: i64, limiter: i64, gradient: i64, p: i64) {
    var e: f64[3] = [0.0; 3];
    for l in 3 { e[l] = static_error(8 << l, kind, a, d, use_bump, flux, limiter, gradient, p); }
    let got = log(e[1] / e[2]) / log(2.0);
    print(what, ": L2 errors ", e[0], ", ", e[1], ", ", e[2], ", order ", got, " (want ", want, ")");
    if got < want - 0.3 { fail(str(what, ": the order is ", got, ", less than ", want, " - 0.3")); }
}

// A square wave advected along x for one half of its width, with the
// limiter: the values must stay in [0, 1] (TVD). Gives the overshoot.
fn square_wave(limiter: i64) -> f64 {
    grid::mesh_rect(0.0, 0.0, 1.0, 0.1, 50, 1);
    grid::equation(1.0, 1.0, 0.0, 0.0);
    grid::fv(grid::LINEAR_UPWIND, grid::GREEN_GAUSS, limiter);
    grid::clear_bcs();
    assert(grid::select_faces() { x, y => abs(x) <= 1.0e-9 || abs(x - 1.0) <= 1.0e-9 });
    grid::fix_selected(0);
    grid::start() { x, y => if x > 0.2 && x < 0.4 { 1.0 } else { 0.0 } };
    if !grid::solve_transient(0.1, 0.0, zero_source, zero) { fail("the square wave is not finite"); }
    var over = 0.0;
    for v in grid::sol { over = max(over, max(v - 1.0, -v)); }
    return over;
}

// The transient manufactured solution u = exp(x - t + y/2) of
// ddt(u) + div(a u) - d div(grad(u)) = f, from t = 0 to 0.1 with the
// stable step.
fn transient_error(n: i64, kind: i64, limiter: i64, p: i64) -> f64 {
    let d = 0.01;
    grid::mesh_rect(0.0, 0.0, 1.0, 1.0, n, n);
    grid::equation(1.0, AX, AY, d);
    if kind == grid::FV { grid::fv(grid::LINEAR_UPWIND, grid::LEAST_SQUARES, limiter); } else { grid::dg(p, grid::UPWIND, 2.0); }
    grid::clear_bcs();
    assert(grid::select_faces() { x, y => true });
    grid::fix_selected(0);
    grid::start() { x, y => smooth(x, y) };
    d_now = d;
    let ok = grid::solve_transient(0.1, 0.0, moving_source, moving_value);
    if !ok { fail("the transient solve is not finite"); }
    let e, norm = grid::l2_error(moving);
    return e / norm;
}

fn check_transient(what: u8[:], want: f64, kind: i64, limiter: i64, p: i64) {
    var e: f64[3] = [0.0; 3];
    for l in 3 { e[l] = transient_error(8 << l, kind, limiter, p); }
    let got = log(e[1] / e[2]) / log(2.0);
    print(what, ": L2 errors ", e[0], ", ", e[1], ", ", e[2], ", order ", got, " (want ", want, ")");
    if got < want - 0.3 { fail(str(what, ": the order is ", got, ", less than ", want, " - 0.3")); }
}

fn main() {
    let fv = grid::FV;
    let dg = grid::DG;
    check_order("fv diffusion, green_gauss", 2.0, fv, 0.0, 1.0, true, grid::UPWIND, grid::NONE, grid::GREEN_GAUSS, 0);
    check_order("fv advection-diffusion, upwind", 1.0, fv, 1.0, 0.01, false, grid::UPWIND, grid::NONE, grid::GREEN_GAUSS, 0);
    check_order("fv advection-diffusion, central", 2.0, fv, 1.0, 0.01, false, grid::CENTRAL, grid::NONE, grid::GREEN_GAUSS, 0);
    for lim, i in [grid::NONE, grid::MINMOD, grid::VAN_LEER, grid::VENKATAKRISHNAN] {
        let what = str("fv advection-diffusion, linear_upwind, green_gauss, limiter ", ["none", "minmod", "van_leer", "venkatakrishnan"][i]);
        check_order(what, 2.0, fv, 1.0, 0.01, false, grid::LINEAR_UPWIND, lim, grid::GREEN_GAUSS, 0);
    }
    check_order("fv advection-diffusion, linear_upwind, least_squares, limiter van_leer", 2.0, fv, 1.0, 0.01, false, grid::LINEAR_UPWIND,
                grid::VAN_LEER, grid::LEAST_SQUARES, 0);
    for p in 1..4 {
        check_order(str("dg p = ", p, " diffusion"), p + 1.0, dg, 0.0, 1.0, true, grid::UPWIND, grid::NONE, 0, p);
        check_order(str("dg p = ", p, " advection-diffusion, upwind"), p + 1.0, dg, 1.0, 0.01, false, grid::UPWIND, grid::NONE, 0, p);
    }
    check_order("dg p = 2 advection-diffusion, lax_friedrichs", 3.0, dg, 1.0, 0.01, false, grid::LAX_FRIEDRICHS, grid::NONE, 0, 2);
    check_order("dg p = 2 advection", 3.0, dg, 1.0, 0.0, false, grid::UPWIND, grid::NONE, 0, 2);
    check_transient("fv ssprk3, linear_upwind, van_leer", 2.0, fv, grid::VAN_LEER, 0);
    check_transient("dg ssprk3, p = 2", 3.0, dg, grid::NONE, 2);
    let free = square_wave(grid::NONE);
    print("square wave, limiter none: overshoot ", free);
    if free < 1.0e-3 { fail("square wave, limiter none: no overshoot, so the test does not see a limiter"); }
    for lim, i in [grid::MINMOD, grid::VAN_LEER] {
        let over = square_wave(lim);
        let what = ["minmod", "van_leer"][i];
        print("square wave, limiter ", what, ": overshoot ", over);
        if over > 1.0e-12 { fail(str("square wave, limiter ", what, ": the overshoot is ", over, ", not 0")); }
    }
    // Venkatakrishnan is not TVD: ε lets a small overshoot through, so
    // that the limiter does not clip a smooth field.
    let venkat = square_wave(grid::VENKATAKRISHNAN);
    print("square wave, limiter venkatakrishnan: overshoot ", venkat);
    if venkat > 0.5 * free { fail(str("square wave, limiter venkatakrishnan: the overshoot is ", venkat, ", more than half of ", free)); }
    if failures > 0 { exit(1); }
    print("grid: all tests pass");
}
