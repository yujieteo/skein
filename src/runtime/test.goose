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

fn main() {
    test_cantilever();
    if failures > 0 { exit(1); }
    print("runtime: all tests pass");
}
