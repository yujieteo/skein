import std;
import math;
import vof from "../vof";
import solvers from "../solvers";

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

var rho_one = 1000.0;
var rho_two = 1.0;

fn properties(a: f64) -> f64[2] { return [a * rho_one + (1.0 - a) * rho_two, a * 1.0e-3 + (1.0 - a) * 1.0e-5]; }
fn no_slip(k: i64, x: f64, y: f64, t: f64) -> f64[2] { return [0.0, 0.0]; }
fn stream(k: i64, x: f64, y: f64, t: f64) -> f64[2] { return [1.0, 0.0]; }
fn after(k: i64, t: f64) {}
fn stream_field(x: f64, y: f64) -> f64[2] { return [1.0, 0.0]; }

fn half_filled(x: f64, y: f64) -> f64 { return if y < 0.5 { 1.0 } else { 0.0 }; }

fn volume() -> f64 {
    var v = 0.0;
    for c in vof::cell_count() { v += vof::alpha[c] * vof::hx * vof::hy; }
    return v;
}

fn largest_velocity() -> f64 {
    var m = 0.0;
    for v in vof::vel { m = max(m, abs(v)); }
    return m;
}

fn test_hydrostatic() {
    vof::mesh_rect(0.0, 0.0, 1.0, 1.0, 20, 20);
    vof::set_gravity(0.0, -9.81);
    vof::set_scheme(vof::LINEAR_UPWIND, vof::VAN_LEER);
    vof::clear_bcs();
    let tol = vof::tolerance();
    let found = vof::select_faces() { x, y => y < 1.0 - tol };
    if !found { fail("hydrostatic: no wall"); return; }
    vof::fix_selected(0);
    vof::start(half_filled);
    if !vof::solve_transient(0.05, 0.001, properties, no_slip, after) { fail("hydrostatic: the solve stops"); return; }
    if largest_velocity() > 1.0e-6 { fail(str("hydrostatic: the velocity at rest is ", largest_velocity())); }
    let bottom = vof::pres[10];
    let want = 1000.0 * 9.81 * (0.5 - 0.025);
    if abs(bottom - want) > 0.01 * want { fail(str("hydrostatic: the pressure at the bottom is ", bottom, ", the weight of the water ", want)); }
}

fn drop(x: f64, y: f64) -> f64 {
    let r = sqrt((x - 0.4) * (x - 0.4) + (y - 0.5) * (y - 0.5));
    return clamp(0.5 - (r - 0.15) / 0.025, 0.0, 1.0);
}

fn test_advection() {
    vof::mesh_rect(0.0, 0.0, 2.0, 1.0, 80, 40);
    vof::set_gravity(0.0, 0.0);
    vof::set_tension(0.0);
    vof::set_scheme(vof::LINEAR_UPWIND, vof::VAN_LEER);
    vof::clear_bcs();
    let tol = vof::tolerance();
    let found = vof::select_faces() { x, y => x < tol || y < tol || y > 1.0 - tol };
    if !found { fail("advection: no wall"); return; }
    vof::fix_selected(0);
    vof::start(drop);
    vof::start_velocity(stream_field);
    let before = volume();
    var cx = 0.0;
    for c in vof::cell_count() { cx += vof::alpha[c] * vof::centre(c)[0] * vof::hx * vof::hy; }
    cx /= before;
    if !vof::solve_transient(0.4, 0.005, properties, stream, after) { fail("advection: the solve stops"); return; }
    var cx2 = 0.0;
    var smeared = 0.0;
    for c in vof::cell_count() {
        cx2 += vof::alpha[c] * vof::centre(c)[0] * vof::hx * vof::hy;
        if vof::alpha[c] > 0.05 && vof::alpha[c] < 0.95 { smeared += vof::hx * vof::hy; }
        if vof::alpha[c] < 0.0 || vof::alpha[c] > 1.0 { fail("advection: alpha leaves 0..1"); return; }
    }
    cx2 /= volume();
    for c in vof::cell_count() {
        if abs(vof::vel[2 * c] - 1.0) > 1.0e-6 || abs(vof::vel[2 * c + 1]) > 1.0e-6 {
            fail(str("advection: the velocity of cell ", c, " is (", vof::vel[2 * c], ", ", vof::vel[2 * c + 1], "), not (1, 0)"));
            return;
        }
    }
    if abs(volume() - before) > 0.01 * before { fail(str("advection: the volume goes from ", before, " to ", volume())); }
    if abs(cx2 - cx - 0.4) > 0.01 { fail(str("advection: the drop moves ", cx2 - cx, ", not 0.4")); }
    if smeared > 0.6 * 2.0 * PI * 0.15 * 0.0125 * 10.0 { fail(str("advection: the interface is smeared over an area of ", smeared)); }
}

fn bubble(x: f64, y: f64) -> f64 {
    let r = sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
    return clamp(0.5 - (r - 0.25) / 0.025, 0.0, 1.0);
}

fn test_laplace() {
    rho_one = 1.0;
    rho_two = 1.0;
    vof::mesh_rect(0.0, 0.0, 1.0, 1.0, 40, 40);
    vof::set_gravity(0.0, 0.0);
    vof::set_tension(0.1);
    vof::set_scheme(vof::LINEAR_UPWIND, vof::VAN_LEER);
    vof::clear_bcs();
    let found = vof::select_faces() { x, y => true };
    if !found { fail("laplace: no wall"); return; }
    vof::fix_selected(0);
    vof::start(bubble);
    if !vof::solve_transient(0.2, 0.002, properties, no_slip, after) { fail("laplace: the solve stops"); return; }
    let inside = vof::pres[vof::locate(0.5, 0.5)];
    let outside = vof::pres[vof::locate(0.06, 0.06)];
    let jump = inside - outside;
    if abs(jump - 0.1 / 0.25) > 0.1 * 0.1 / 0.25 { fail(str("laplace: the pressure jump is ", jump, ", σ/R is ", 0.1 / 0.25)); }
    if largest_velocity() > 0.05 { fail(str("laplace: the parasitic velocity is ", largest_velocity())); }
    rho_one = 1000.0;
    rho_two = 1.0;
}

fn main() {
    test_hydrostatic();
    test_advection();
    test_laplace();
    if failures > 0 { exit(1); }
    print("vof: ok");
}
