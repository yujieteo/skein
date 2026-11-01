// Tests of the flow solver. `./build.sh test` runs this program from the
// root of the repository, and a non-zero exit status is a failure.
//
// Plane Poiseuille flow in a channel of length L and height H, driven by a
// body force G along x, with v = 0 on the walls and v.y = 0 at the ends:
// u = G/(2μ) y (H - y), v.y = 0, p = 0. The convection term (∇v) v is 0,
// so Stokes and Navier-Stokes have the same solution. P2 velocity holds it
// exactly. The force of each wall on the flow is -G H L / 2 along x.
//
// A cantilever of length 10 and height 1 in plane strain, held at x = 0,
// with a shear load P at x = L on the solid of fsi.goose: the tip
// deflection of Timoshenko's beam, P L³/(3 E' I) + P L/(κ G A), with
// E' = E/(1 - ν²) and κ = 5/6.
//
// The gate of phase 3 is samples/02_cavity.skn, with an expect for each
// row of Ghia et al. (1982), and samples/06_cylinder_sweep.skn.
// scripts/generated_c.sh runs the two samples.
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

// A fully developed turbulent channel of half height 1 at Re_τ = 395,
// driven by the body force ρ u_τ² / h with u_τ = 1, ρ = 1, ν = 1/395,
// with the k-ω SST model. The ends are periodic: the mesh is half the
// channel and its mirror, so that the nodes of the ends match. The DNS of
// Moser, Kim and Mansour (1999) gives U_b+ = 17.55 and U_c+ = 20.1. SST
// gives about 6% more on this mesh, and 5% more on a mesh of half the size.
fn channel() {
    let len = 0.1;
    let hw = 1.0 / 395.0;
    mesher::clear_shapes();
    mesher::clear_refines();
    mesher::refine_segment(0.0, 0.0, len, 0.0, hw);
    mesher::refine_segment(0.0, 2.0, len, 2.0, hw);
    if !mesher::mesh_2d(mesher::rect(0.0, 0.0, 0.5 * len, 2.0), 0.1) { fail("channel: the mesher did not finish"); }
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    var image: i64[>..] = [];
    let half = mesher::point_count();
    for n in half {
        xy.push(mesher::vertex(n)[0]);
        xy.push(mesher::vertex(n)[1]);
    }
    for n in half {
        let x = mesher::vertex(n)[0];
        if abs(x - 0.5 * len) <= 1.0e-9 { image.push(n); continue; }
        image.push(xy.len / 2);
        xy.push(len - x);
        xy.push(mesher::vertex(n)[1]);
    }
    for c in mesher::cell_count() {
        for k in 3 { tri.push(mesher::cell_node(c, k)); }
        for k in 3 { tri.push(image[mesher::cell_node(c, k)]); }
    }
    flow::set_turbulent(true);
    flow::mesh_from_arrays(xy, tri, 1);
    flow::fluid(1.0, 1.0 / 395.0, true, flow::SUPG_PSPG);
    flow::body_force(1.0, 0.0);
    let tol = flow::tolerance();
    assert(flow::select_boundary() { x, y => abs(x - len) <= tol });
    if !flow::periodic(len, 0.0) { fail("channel: the ends do not match"); }
    flow::clear_bcs();
    assert(flow::select_boundary() { x, y => abs(y) <= tol || abs(y - 2.0) <= tol });
    flow::fix_velocity(0, 0.0);
    flow::fix_velocity(1, 0.0);
    if !flow::solve_steady(1.0e-8) { fail("channel: the solve does not converge"); }
    var bulk = 0.0;
    for i in 200 {
        let u, found = flow::value_at(flow::VALUE_VX, 0.5 * len, (i + 0.5) / 100.0);
        bulk += u / 200.0;
    }
    let centre, found = flow::value_at(flow::VALUE_VX, 0.5 * len, 1.0);
    print("channel, k-omega SST: ", flow::cell_count(), " cells, ", flow::dof_count(), " dofs, ", flow::newton_steps, " steps, ",
          flow::linear_iterations, " GMRES iterations; U_b+ ", bulk, " (DNS 17.55), U_c+ ", centre, " (DNS 20.1)");
    if abs(bulk - 17.55) > 0.07 * 17.55 { fail(str("channel: U_b+ is ", bulk, ", not within 7% of 17.55")); }
    if abs(centre - 20.1) > 0.07 * 20.1 { fail(str("channel: U_c+ is ", centre, ", not within 7% of 20.1")); }
    flow::set_turbulent(false);
}

fn cantilever() {
    let len = 10.0;
    let e = 1000.0;
    let nu = 0.3;
    let p = 0.1;
    mesher::clear_shapes();
    mesher::clear_refines();
    if !mesher::mesh_2d(mesher::rect(0.0, 0.0, len, 1.0), 0.25) { fail("cantilever: the mesher did not finish"); }
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    for n in mesher::point_count() {
        xy.push(mesher::vertex(n)[0]);
        xy.push(mesher::vertex(n)[1]);
    }
    for c in mesher::cell_count() {
        for k in 3 { tri.push(mesher::cell_node(c, k)); }
    }
    flow::solid_from_arrays(xy, tri);
    flow::solid(e, nu);
    assert(flow::select_solid_boundary() { x, y => abs(x) <= 1.0e-9 });
    flow::fix_solid();
    if !flow::factor_solid() { fail("cantilever: the stiffness is not positive definite"); }
    // the consistent load of a uniform shear on the edges at x = L
    var load: f64[>..] = [];
    var u: f64[>..] = [];
    solvers::zeros(load, 2 * flow::solid_nodes());
    solvers::zeros(u, 2 * flow::solid_nodes());
    for c in flow::solid_cells() {
        for k in 3 {
            let a = flow::s_conn[6 * c + k];
            let b = flow::s_conn[6 * c + (k + 1) % 3];
            if abs(flow::s_coords[2 * a] - len) > 1.0e-9 || abs(flow::s_coords[2 * b] - len) > 1.0e-9 { continue; }
            let t = p * abs(flow::s_coords[2 * a + 1] - flow::s_coords[2 * b + 1]);
            load[2 * a + 1] -= t / 6.0;
            load[2 * b + 1] -= t / 6.0;
            load[2 * flow::s_conn[6 * c + 3 + k] + 1] -= 4.0 * t / 6.0;
        }
    }
    flow::solve_solid(load, u);
    for d in u.len { flow::s_u[d] = u[d]; }
    let tip, found = flow::solid_value_at(1, len, 0.5);
    let ep = e / (1.0 - nu * nu);
    let exact = -(p * len * len * len / (3.0 * ep / 12.0) + p * len / (5.0 / 6.0 * e / (2.0 * (1.0 + nu))));
    print("cantilever, P2 plane strain: ", flow::solid_cells(), " cells, tip ", tip, " (Timoshenko ", exact, ")");
    if !found || abs(tip - exact) > 0.02 * abs(exact) { fail(str("cantilever: the tip deflection is ", tip, ", not within 2% of ", exact)); }
}

fn main() {
    poiseuille("Stokes, P2/P1", 2, false, flow::NONE, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg", 2, true, flow::SUPG, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg_pspg_lsic", 2, true, flow::SUPG_PSPG_LSIC, 1.0e-6);
    poiseuille("Navier-Stokes, P1/P1, supg_pspg", 1, true, flow::SUPG_PSPG, 0.03);
    channel();
    cantilever();
    if failures > 0 { exit(1); }
    print("flow: all tests pass");
}
