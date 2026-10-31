// Tests of the flow solver. `./build.sh test` runs this program from the
// root of the repository, and a non-zero exit status is a failure.
//
// Plane Poiseuille flow in a channel of length L and height H, driven by a
// body force G along x, with v = 0 on the walls and v.y = 0 at the ends:
// u = G/(2μ) y (H - y), v.y = 0, p = 0. The convection term (∇v) v is 0,
// so Stokes and Navier-Stokes have the same solution. P2 velocity holds it
// exactly. The force of each wall on the flow is -G H L / 2 along x.
//
// The gate of phase 3 is the lid-driven cavity of samples/02_cavity.skn
// against Ghia et al. (1982): scripts/generated_c.sh runs it.
import std;
import math;
import flow from "../flow";
import mesher from "../mesher";
import runtime from "../runtime";   // the mesher imports it
import solvers from "../solvers";   // the runtime imports it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

let L = 4.0;
let H = 1.0;
let G = 2.0;
let MU = 0.5;

fn exact_u(y: f64) -> f64 { return G / (2.0 * MU) * y * (H - y); }

// Solves the channel and gives the largest error of u at 9 points on the
// line x = L/2, relative to the peak G H²/(8 μ).
fn poiseuille(what: u8[:], order: i64, navier_stokes: bool, stab: i64, tolerance: f64) {
    mesher::clear_shapes();
    mesher::clear_refines();
    if !mesher::mesh_2d(mesher::rect(0.0, 0.0, L, H), 0.125) { fail(str(what, ": the mesher did not finish")); }
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    for n in mesher::point_count() {
        xy.push(mesher::vertex(n)[0]);
        xy.push(mesher::vertex(n)[1]);
    }
    for c in mesher::cell_count() {
        for k in 3 { tri.push(mesher::cell_node(c, k)); }
    }
    flow::mesh_from_arrays(xy, tri, order);
    flow::fluid(1.0, MU, navier_stokes, stab);
    flow::body_force(G, 0.0);
    flow::clear_bcs();
    let tol = flow::tolerance();
    assert(flow::select_boundary() { x, y => abs(x) <= tol || abs(x - L) <= tol });
    flow::fix_velocity(1, 0.0);
    assert(flow::select_boundary() { x, y => abs(y) <= tol || abs(y - H) <= tol });
    flow::fix_velocity(0, 0.0);
    flow::fix_velocity(1, 0.0);
    let converged = flow::solve_steady(1.0e-10);
    if !converged { fail(str(what, ": Newton does not converge")); }
    let peak = G * H * H / (8.0 * MU);
    var error = 0.0;
    for i in 1..10 {
        let y = i * H / 10.0;
        let u, found = flow::value_at(flow::VALUE_VX, 0.5 * L, y);
        if !found { fail(str(what, ": no cell holds (", 0.5 * L, ", ", y, ")")); }
        error = max(error, abs(u - exact_u(y)) / peak);
    }
    assert(flow::select_boundary() { x, y => abs(y - H) <= tol });
    let fx = flow::region_force(0);
    let want = -0.5 * G * H * L;
    print(what, ": ", flow::cell_count(), " cells, ", flow::dof_count(), " dofs, ", flow::newton_steps, " Newton steps, ",
          flow::linear_iterations, " GMRES iterations; error of u ", error, ", force of the top wall ", fx, " (want ", want, ")");
    if error > tolerance { fail(str(what, ": the error of u is ", error, ", more than ", tolerance)); }
    if abs(fx - want) > tolerance * abs(want) { fail(str(what, ": the force of the top wall is ", fx, ", not ", want)); }
}

fn main() {
    poiseuille("Stokes, P2/P1", 2, false, flow::NONE, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg", 2, true, flow::SUPG, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg_pspg_lsic", 2, true, flow::SUPG_PSPG_LSIC, 1.0e-6);
    poiseuille("Navier-Stokes, P1/P1, supg_pspg", 1, true, flow::SUPG_PSPG, 0.03);
    if failures > 0 { exit(1); }
    print("flow: all tests pass");
}
