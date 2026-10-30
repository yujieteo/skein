// Tests of the runtime. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
//
// The acceptance case is femahoy validation/cases/grid-cantilever
// (examples/femahoy/grid_cantilever.skn): a steel beam 100 × 10 × 10 mm,
// E = 210000 MPa, nu = 0.3, fixed at x = 0, with a traction of -1 MPa in z
// on the tip x = 100 mm, so P = 100 N. Units: mm, N, MPa.
import std;
import math;
import os;
import runtime from "../runtime";

var failures = 0;

fn expect_near(what: u8[:], got: f64, want: f64, rel: f64) {
    let err = abs(got - want) / abs(want);
    if err > rel {
        failures++;
        print("FAIL ", what, ": got ", got, ", want ", want, " within ", 100.0 * rel, " %");
    }
}

let L = 100.0;
let B = 10.0;
let H = 10.0;
let E = 210000.0;
let NU = 0.3;
let P = 100.0;

// The femahoy grid result at max_cell 2.5 mm (derivation.md, mesh study).
let FEMAHOY_TIP = -0.190113;

fn euler_bernoulli_tip() -> f64 {
    let inertia = B * H * H * H / 12.0;
    return -P * L * L * L / (3.0 * E * inertia);
}

// The cantilever on 2.5 mm cells (40 × 4 × 4 27-node hexahedra), as the
// femahoy grid case.
fn test_cantilever() {
    runtime::mesh_box(0.0, 0.0, 0.0, L, B, H, 40, 4, 4);
    runtime::material(E, NU);
    assert(runtime::fix_face(runtime::X, 0.0));
    assert(runtime::traction_face(runtime::X, L, 0.0, 0.0, -P / (B * H)));
    if !runtime::solve_static() {
        failures++;
        print("FAIL cantilever: CG does not converge in ", runtime::cg_iterations(), " iterations");
    }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    let rx, ry, rz = runtime::reaction_face(runtime::X, 0.0);
    let fx, fy, fz = runtime::applied_load();
    print("cantilever: 40 x 4 x 4 hex27, ", runtime::dof_count(), " dofs, ",
        runtime::cg_iterations(), " CG iterations, ", runtime::solve_seconds(), " s");
    print("  tip uz ", uz, " mm (femahoy ", FEMAHOY_TIP, ", Euler-Bernoulli ",
        euler_bernoulli_tip(), ")");
    print("  reaction ", rx, " ", ry, " ", rz, " N, applied ", fx, " ", fy, " ", fz, " N");
    expect_near("tip uz against femahoy", uz, FEMAHOY_TIP, 0.02);
    expect_near("tip uz against Euler-Bernoulli", uz, euler_bernoulli_tip(), 0.02);
    expect_near("applied load", fz, -P, 1.0e-12);
    // The equilibrium check of docs/checks.md: 1e-6 of the largest load.
    if abs(rz + fz) > 1.0e-6 * P || abs(rx + fx) > 1.0e-6 * P || abs(ry + fy) > 1.0e-6 * P {
        failures++;
        print("FAIL equilibrium: reaction ", rx, " ", ry, " ", rz, ", load ", fx, " ", fy, " ", fz);
    }
    let path = "build/test_runtime_cantilever.vtu";
    if !runtime::write_vtu(path) {
        failures++;
        print("FAIL write_vtu: cannot write ", path);
    }
}

// A second load step on the same mesh: a point force at the tip, after
// `clear_loads`. A plane that is not a face of the box selects nothing.
fn test_second_step() {
    runtime::clear_loads();
    runtime::force_point(L, 0.5 * B, 0.5 * H, 0.0, 0.0, -2.0 * P);
    if runtime::traction_face(runtime::X, 0.5 * L, 0.0, 0.0, 1.0) {
        failures++;
        print("FAIL traction_face: x = L/2 is not a face of the box");
    }
    if !runtime::solve_static() {
        failures++;
        print("FAIL second step: CG does not converge");
    }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, H);
    let rx, ry, rz = runtime::reaction_face(runtime::X, 0.0);
    let fx, fy, fz = runtime::applied_load();
    print("second step: point force ", fz, " N, tip uz at the top ", uz, " mm, reaction z ", rz, " N");
    expect_near("second step applied load", fz, -2.0 * P, 1.0e-12);
    expect_near("second step reaction", rz, 2.0 * P, 1.0e-6);
    // The beam is linear: twice the load gives about twice the deflection.
    expect_near("second step tip uz", uz, 2.0 * FEMAHOY_TIP, 0.02);
}

// A quarter of a closed-end steel tube, a = 100 mm, b = 105 mm, 10 mm
// long, with an internal pressure of 1 MPa. Symmetry on x = 0, y = 0 and
// z = 0. The closed end is the uniform axial stress A = p a² / (b² - a²)
// on z = 10 mm. The Lamé solution is exact for these loads:
// σθ = A (1 + b²/r²), σr = A (1 - b²/r²) and σz = A. The 0.1% limit is
// about 4 times the error of this mesh.
fn test_annulus() {
    let a = 100.0;
    let b = 105.0;
    let p = 1.0;
    let A = p * a * a / (b * b - a * a);
    runtime::mesh_annulus(0.0, 0.0, a, b, 0.0, 0.5 * PI, 0.0, 10.0, 2, 48, 2);
    runtime::material(E, NU);
    assert(runtime::fix_face_component(runtime::THETA, 0.5 * PI, runtime::X));
    assert(runtime::fix_face_component(runtime::THETA, 0.0, runtime::Y));
    assert(runtime::fix_face_component(runtime::Z, 0.0, runtime::Z));
    assert(runtime::normal_traction_face(runtime::R, a, -p));
    assert(runtime::traction_face(runtime::Z, 10.0, 0.0, 0.0, A));
    if !runtime::solve_static() {
        failures++;
        print("FAIL annulus: CG does not converge in ", runtime::cg_iterations(), " iterations");
    }
    let r = 102.5;
    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(r, 0.0, 5.0);
    let rx, ry, rz = runtime::reaction_face(runtime::Z, 0.0);
    let fx, fy, fz = runtime::applied_load();
    print("annulus: 2 x 48 x 2 hex27, ", runtime::dof_count(), " dofs, ",
        runtime::cg_iterations(), " CG iterations, ", runtime::solve_seconds(), " s");
    print("  at r = ", r, ": hoop ", syy, ", axial ", szz, " MPa");
    print("  reaction z ", rz, " N, applied ", fx, " ", fy, " ", fz, " N");
    expect_near("annulus hoop against Lame", syy, A * (1.0 + b * b / (r * r)), 1.0e-3);
    expect_near("annulus axial against Lame", szz, A, 1.0e-3);
    // The axial load is the pressure on the end cap, p π a² / 4.
    expect_near("annulus axial reaction", rz, -0.25 * PI * a * a * p, 1.0e-6);
    // The hoop stress at 45°: (σθ + σr) / 2 on xx and yy, (σθ - σr) / 2 on xy.
    let s45 = r * sqrt(0.5);
    let txx, tyy, tzz, tyz, txz, txy = runtime::stress_at(s45, s45, 5.0);
    expect_near("annulus shear at 45 degrees", -txy, A * b * b / (r * r), 1.0e-2);
    let ox, oy, oz = runtime::displacement_at(0.0, 0.0, 5.0);
    if ox != 0.0 || oy != 0.0 || oz != 0.0 {
        failures++;
        print("FAIL annulus: the axis is not in the mesh, but u there is ", ox, " ", oy, " ", oz);
    }
}

fn main() {
    test_cantilever();
    test_second_step();
    test_annulus();
    if failures > 0 { exit(1); }
    print("runtime: all tests pass");
}
