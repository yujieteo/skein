import std;
import math;
import flow from "../flow";
import mesher from "../mesher";
import runtime from "../runtime";
import solvers from "../solvers";

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

fn channel(sa: bool) {
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
    if sa { flow::set_spalart(true); } else { flow::set_turbulent(true); }
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
    if sa {
        flow::wall_distance_region();
        flow::fix_spalart(0.0);
        flow::spalart_start(0.0, 0.0, 3.0 / 395.0);
    }
    if !flow::solve_steady(1.0e-8) { fail("channel: the solve does not converge"); }
    var bulk = 0.0;
    for i in 200 {
        let u, found = flow::value_at(flow::VALUE_VX, 0.5 * len, (i + 0.5) / 100.0);
        bulk += u / 200.0;
    }
    let centre, found = flow::value_at(flow::VALUE_VX, 0.5 * len, 1.0);
    print("channel, ", if sa { "Spalart-Allmaras" } else { "k-omega SST" }, ": ", flow::cell_count(), " cells, ", flow::dof_count(), " dofs, ", flow::newton_steps, " steps, ",
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

fn quadrature_weights_sum_to_one_and_points_are_barycentric() {
    var total = 0.0;
    for q in 7 {
        let l, w = flow::quad_point(q);
        total += w;
        if abs(l[0] + l[1] + l[2] - 1.0) > 1.0e-12 { fail(str("quad_point ", q, ": the barycentric coordinates do not sum to 1")); }
    }
    if abs(total - 1.0) > 1.0e-12 { fail(str("quad_point: the weights sum to ", total, ", not 1")); }
}

fn main() {
    quadrature_weights_sum_to_one_and_points_are_barycentric();
    poiseuille("Stokes, P2/P1", 2, false, flow::NONE, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg", 2, true, flow::SUPG, 1.0e-6);
    poiseuille("Navier-Stokes, P2/P1, supg_pspg_lsic", 2, true, flow::SUPG_PSPG_LSIC, 1.0e-6);
    poiseuille("Navier-Stokes, P1/P1, supg_pspg", 1, true, flow::SUPG_PSPG, 0.03);
    channel(false);
    channel(true);
    cantilever();
    if failures > 0 { exit(1); }
    print("flow: all tests pass");
}
