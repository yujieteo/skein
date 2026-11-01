// Tests of the package optim. `./build.sh test` runs this program from
// the root of the repository, and a non-zero exit status is a failure.
//
// - the patch test: a uniform tension gives the exact linear field;
// - the adjoint: the derivatives of a point objective by the density of
//   cells and by E, nu, t and the loads equal central finite differences;
// - the half MBB beam of Andreassen et al. (2011), 60 x 20 cells: the
//   loop gives the compliance c = 216.81 of their figure 3 within 0.2 %
//   (femahoy validation/cases/topology-mbb-beam).
import std;
import math;
import optim from "../optim";
import solvers from "../solvers";   // optim imports it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

fn near(what: u8[:], got: f64, want: f64, rel: f64) {
    if !(abs(got - want) <= rel * max(abs(want), 1.0e-300)) { fail(str(what, ": ", got, ", not ", want, " within ", rel)); }
}

// A plate of 4 x 2 cells of 2 x 1, fixed on x = 0 along x and at the
// origin along y, with the traction s along x on x = 8: u = (s x / E,
// -nu s y / E) in plane stress, at each node.
fn test_patch() {
    optim::mesh_rect(0.0, 0.0, 8.0, 2.0, 4, 2);
    optim::elastic(200.0, 0.25, 1.0, optim::PLANE_STRESS);
    optim::clear_bcs();
    assert(optim::fix_line(optim::X, 0.0, optim::X, 0.0));
    assert(optim::fix_point(0.0, 0.0, optim::Y, 0.0));
    assert(optim::edge_load(optim::X, 8.0, 3.0, 0.0));
    if !optim::solve() { fail("the patch test is singular"); }
    var worst = 0.0;
    for i in 5 {
        for j in 3 {
            let x = 2.0 * i;
            let y = 1.0 * j;
            let ux, _ = optim::value_at(optim::X, x, y);
            let uy, _ = optim::value_at(optim::Y, x, y);
            worst = max(worst, max(abs(ux - 3.0 * x / 200.0), abs(uy + 0.25 * 3.0 * y / 200.0)));
        }
    }
    if worst > 1.0e-12 { fail(str("the patch test is not exact: ", worst)); }
}

// A small design with densities that are not all equal, held as a
// cantilever, a point load at the free end, and the objective u.x at an
// inside point: not a compliance, so that λ is not u.
fn small_design(E: f64, nu: f64, t: f64, fy: f64) {
    optim::mesh_rect(0.0, 0.0, 6.0, 3.0, 6, 3);
    optim::elastic(E, nu, t, optim::PLANE_STRESS);
    optim::penalty = 3.0;
    for c in optim::cell_count() { optim::density[c] = 0.3 + 0.6 * ((7 * c) % 11) / 10.0; }
    optim::clear_bcs();
    assert(optim::fix_line(optim::X, 0.0, optim::X, 0.0));
    assert(optim::fix_line(optim::X, 0.0, optim::Y, 0.0));
    assert(optim::point_load(6.0, 0.0, 0.5, fy));
    assert(optim::edge_load(optim::Y, 3.0, 0.0, -0.2));
    assert(optim::objective(optim::X, 3.5, 1.5, 2.0));
}

fn objective_of(E: f64, nu: f64, t: f64, fy: f64, cell: i64, dx: f64) -> f64 {
    small_design(E, nu, t, fy);
    if cell >= 0 { optim::density[cell] += dx; }
    if !optim::solve() { fail("the small design is singular"); }
    return optim::objective_value();
}

fn test_adjoint() {
    small_design(1.0, 0.3, 2.0, -1.0);
    if !optim::solve() { fail("the small design is singular"); }
    optim::solve_adjoint();
    let g = optim::input_gradient();
    var dc: f64[>..] = [];
    for c in optim::cell_count() { dc.push(optim::cell_sensitivity(c)); }
    let h = 1.0e-5;
    // The error of a difference is 1e-6 of the largest sensitivity: a
    // cell far from the load and the point has a small one.
    var largest = 0.0;
    for v in dc { largest = max(largest, abs(v)); }
    for c in [0, 5, 11, 17] {
        let fd = (objective_of(1.0, 0.3, 2.0, -1.0, c, h) - objective_of(1.0, 0.3, 2.0, -1.0, c, -h)) / (2.0 * h);
        if abs(dc[c] - fd) > 1.0e-6 * largest { fail(str("dJ/dx of cell ", c, ": ", dc[c], ", not ", fd)); }
    }
    near("dJ/dE", g[0], (objective_of(1.0 + h, 0.3, 2.0, -1.0, -1, 0.0) - objective_of(1.0 - h, 0.3, 2.0, -1.0, -1, 0.0)) / (2.0 * h), 1.0e-6);
    near("dJ/dnu", g[1], (objective_of(1.0, 0.3 + h, 2.0, -1.0, -1, 0.0) - objective_of(1.0, 0.3 - h, 2.0, -1.0, -1, 0.0)) / (2.0 * h), 1.0e-6);
    near("dJ/dfy", g[4], (objective_of(1.0, 0.3, 2.0, -1.0 + h, -1, 0.0) - objective_of(1.0, 0.3, 2.0, -1.0 - h, -1, 0.0)) / (2.0 * h), 1.0e-6);
    near("dJ/dt", g[2], (objective_of(1.0, 0.3, 2.0 + h, -1.0, -1, 0.0) - objective_of(1.0, 0.3, 2.0 - h, -1.0, -1, 0.0)) / (2.0 * h), 1.0e-6);
}

// The half MBB beam: ux = 0 on the left edge (a plane of symmetry), uy =
// 0 at the bottom right corner, and the force 1 down at the top left
// corner. E = 1, nu = 0.3, plane stress, a volume of 0.5, p = 3 and the
// filter radius 2.4. The objective is the compliance -u_y at the load.
fn test_mbb() {
    optim::mesh_rect(0.0, 0.0, 60.0, 20.0, 60, 20);
    optim::elastic(1.0, 0.3, 1.0, optim::PLANE_STRESS);
    optim::clear_bcs();
    assert(optim::fix_line(optim::X, 0.0, optim::X, 0.0));
    assert(optim::fix_point(60.0, 0.0, optim::Y, 0.0));
    assert(optim::point_load(0.0, 20.0, 0.0, -1.0));
    assert(optim::objective(optim::Y, 0.0, 20.0, -1.0));
    if !optim::topology(0.5, 3.0, 2.4, 300) { fail("the MBB beam is singular"); }
    let k = optim::history_objective.len;
    let c = optim::history_objective[k - 1];
    print("optim: the MBB beam gives c = ", c, " in ", k, " iterations, with the volume ", optim::history_volume[k - 1]);
    near("the compliance of the MBB beam", c, 216.81, 0.002);
    near("the volume of the MBB beam", optim::history_volume[k - 1], 0.5, 0.002);
}

fn main() {
    test_patch();
    test_adjoint();
    test_mbb();
    if failures > 0 {
        print("optim: ", failures, " failures");
        exit(1);
    }
    print("optim: all tests pass");
}
