// Tests of the meshless methods. `./build.sh test` runs this program from
// the root of the repository, and a non-zero exit status is a failure.
//
// Each manufactured solution is solved on 3 clouds, each with half the
// spacing of the one before, and the observed order of the L2 error on
// the two finest clouds must be within 0.3 of 2 (docs/checks.md §4.3):
// rkpm with the linear basis, and gfd with its fit of order 2. The models
// of examples/meshless/ check the same orders through
// `skein verify convergence`.
import std;
import math;
import meshless from "../meshless";
import solvers from "../solvers";   // the meshless import it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

// The smooth field exp(x + y/2), and its source for a = (1, 0.5), d.
let AX = 1.0;
let AY = 0.5;
let D = 0.5;
fn smooth(x: f64, y: f64) -> f64 { return exp(x + 0.5 * y); }
fn smooth_source(x: f64, y: f64, t: f64) -> f64 { return (AX + 0.5 * AY - D * 1.25) * smooth(x, y); }
fn smooth_value(k: i64, x: f64, y: f64, t: f64) -> f64 { return smooth(x, y); }

// Solves the problem on the jittered n × n lattice of the unit square,
// fixed on each edge, and gives the relative L2 error.
fn error_on(n: i64, use_gfd: bool) -> f64 {
    meshless::points_rect(0.0, 0.0, 1.0, 1.0, n, n);
    meshless::equation(AX, AY, D);
    let ready = if use_gfd { meshless::gfd(12) } else { meshless::rkpm(2.5) };
    if !ready { fail(str("the method has no shape functions on ", n, " x ", n, " points")); }
    meshless::clear_bcs();
    assert(meshless::select_faces() { x, y => true });
    meshless::fix_selected(0);
    if !meshless::solve_static(smooth_source, smooth_value) { fail(str("the solve on ", n, " x ", n, " is singular")); }
    let e, norm = meshless::l2_error() { x, y, t => smooth(x, y) };
    return e / norm;
}

fn check_order(what: u8[:], use_gfd: bool) {
    var e: f64[3] = [0.0; 3];
    for l in 3 { e[l] = error_on(8 << l, use_gfd); }
    let got = log(e[1] / e[2]) / log(2.0);
    print(what, ": L2 errors ", e[0], ", ", e[1], ", ", e[2], ", order ", got, " (want 2)");
    if got < 2.0 - 0.3 { fail(str(what, ": the order is ", got, ", less than 2 - 0.3")); }
}

// The linear field u = (1e-3 + 2e-3 x - 1e-3 y, -2e-3 + 1e-3 x + 3e-3 y) is
// exact for elasticity with no load, so rkpm with SCNI and Nitsche must
// give it on a cloud that is not a grid (the patch test).
fn linear_u(k: i64, x: f64, y: f64, t: f64) -> f64 {
    return if k == 0 { 1.0e-3 + 2.0e-3 * x - 1.0e-3 * y } else { -2.0e-3 + 1.0e-3 * x + 3.0e-3 * y };
}
fn no_source(x: f64, y: f64, t: f64) -> f64 { return 0.0; }

fn patch_test() {
    meshless::points_rect(0.0, 0.0, 2.0, 1.0, 10, 5);
    meshless::elastic(1000.0, 0.3, true, 1.0);
    if !meshless::rkpm(2.5) { fail("patch test: no shape functions"); }
    meshless::clear_bcs();
    assert(meshless::select_faces() { x, y => true });
    meshless::fix_component(0, 0);
    meshless::fix_component(1, 1);
    if !meshless::solve_static(no_source, linear_u) { fail("patch test: the solve is singular"); }
    var worst = 0.0;
    for p in [[0.3, 0.2], [1.1, 0.5], [1.9, 0.9], [0.77, 0.61]] {
        for c in 2 {
            let v, found = meshless::value_at(c, p[0], p[1]);
            if !found { fail("patch test: a point inside is not found"); }
            worst = max(worst, abs(v - linear_u(c, p[0], p[1], 0.0)));
        }
        let sxx, found = meshless::stress_at(0, 0, p[0], p[1]);
        // plane stress: sxx = E/(1 - nu²) (exx + nu eyy)
        worst = max(worst, abs(sxx - 1000.0 / 0.91 * (2.0e-3 + 0.3 * 3.0e-3)) * 1.0e-3);
    }
    let r, size, f = meshless::balance();
    print("patch test: largest error ", worst, ", reaction (", r[0], ", ", r[1], ")");
    if worst > 1.0e-8 { fail(str("patch test: the error is ", worst)); }
    if abs(r[0]) + abs(r[1]) > 1.0e-8 * size { fail("patch test: the reaction is not in balance"); }
}

fn main() {
    patch_test();
    check_order("rkpm, advection-diffusion", false);
    check_order("gfd, advection-diffusion", true);
    if failures > 0 { exit(1); }
    print("meshless: all tests pass");
}
