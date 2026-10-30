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

// The linear problem through the Newton interface: R(u) = K u - f, with
// the Jacobian K. One Newton iteration solves it.
fn linear_residual(u: f64[:], r: f64[:]) {
    runtime::multiply(u, r);
    for d in r.len { r[d] -= runtime::load[d]; }
}

fn linear_jacobian(u: f64[:]) { assert(runtime::assemble()); }

fn test_newton() {
    runtime::clear_loads();
    assert(runtime::traction_face(runtime::X, L, 0.0, 0.0, -P / (B * H)));
    assert(runtime::assemble());
    let converged, it = runtime::newton(runtime::disp, 1.0e-8, 5, linear_residual, linear_jacobian);
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    print("newton: ", it, " iterations, tip uz ", uz, " mm");
    if !converged || it != 1 {
        failures++;
        print("FAIL newton: converged ", converged, " in ", it, " iterations, want 1");
    }
    expect_near("newton tip uz against the linear solve", uz, -0.1905018167560175, 1.0e-6);
}

// The face sets of the box, and a cell type with no registry entry.
fn test_registry() {
    runtime::mesh_box(0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);
    let s = runtime::find_set("x1", true);
    if s < 0 || runtime::set_members(s).len != 1 || runtime::set_members(s)[0] != 6 * 1 + 1 {
        failures++;
        print("FAIL face set x1: want cell 1, local face 1");
    }
    let n = runtime::find_set("z0", false);
    if n < 0 || runtime::set_members(n).len != 5 * 3 {
        failures++;
        print("FAIL node set z0: want 15 nodes");
    }
    var nodes: i64[4] = [0, 1, 5, 25];
    runtime::add_cell(runtime::TET4, nodes);
    runtime::material(E, NU);
    assert(runtime::fix_face(runtime::X, 0.0));
    if runtime::solve_static() {
        failures++;
        print("FAIL registry: a tet4 cell with no element routine solves");
    }
}

fn main() {
    test_cantilever();
    test_second_step();
    test_newton();
    test_registry();
    if failures > 0 { exit(1); }
    print("runtime: all tests pass");
}
