// Tests of the simplex elements. `./build.sh test` runs this program from
// the root of the repository, and a non-zero exit status is a failure.
//
// The gates of phase 2 (docs/roadmap.md):
// - the patch test: a linear field on a mesh of distorted cells is exact,
//   for each element (P1 and P2 triangles and tetrahedra), in elasticity
//   and in heat conduction;
// - the plate with a hole under tension: the stress at the top of the hole
//   against the Kirsch factor 3;
// - a square bar with a uniform heat source: the centre temperature against
//   the series solution (Timoshenko and Goodier 1970, §109).
import std;
import math;
import os;
import elements from "../elements";
import mesher from "../mesher";
import runtime from "../runtime";   // the mesher imports it
import solvers from "../solvers";   // the runtime imports it

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

var failures = 0;

fn expect_near(what: u8[:], got: f64, want: f64, rel: f64) {
    let err = abs(got - want) / abs(want);
    if !(err <= rel) {
        failures++;
        print("FAIL ", what, ": got ", got, ", want ", want, " within ", 100.0 * rel, " %");
    }
}

fn kind_name(d: i64, ord: i64) -> u8[>..] {
    return str(if d == 2 { "tri" } else { "tet" }, elements::node_count_of(d, ord));
}

// --- the patch test ------------------------------------------------------------

// Moves each vertex node inside the unit box by up to 0.2 of a lattice
// cell, then puts each midpoint node back on the midpoint of its edge, so
// the cells are distorted but straight.
fn distort(cells: i64) {
    var seed: u64 = 12345;
    let h = 1.0 / cells;
    let tol = 1.0e-9;
    for n in elements::node_count() {
        for t in elements::dim {
            seed = seed * 6364136223846793005 + 1442695040888963407;
            let x = elements::coords[3 * n + t];
            if x > tol && x < 1.0 - tol {
                let r = ((seed >> 11) as f64) / 9007199254740992.0 - 0.5;
                elements::coords[3 * n + t] = x + 0.4 * h * r;
            }
        }
    }
    guard elements::order == 2 else { return; }
    let d = elements::dim;
    for c in elements::cell_count() {
        for e in elements::edge_count(d) {
            let a = elements::cell_node(c, elements::edge_node(d, e, 0));
            let b = elements::cell_node(c, elements::edge_node(d, e, 1));
            let m = elements::cell_node(c, d + 1 + e);
            for t in 3 {
                elements::coords[3 * m + t] = 0.5 * (elements::coords[3 * a + t] + elements::coords[3 * b + t]);
            }
        }
    }
}

fn on_boundary(n: i64) -> bool {
    for t in elements::dim {
        let x = elements::coords[3 * n + t];
        if x < 1.0e-9 || x > 1.0 - 1.0e-9 { return true; }
    }
    return false;
}

// The linear field of the patch test: u_a = Σt G[a][t] x_t + 0.01 (a + 1).
let G: f64[9] = [0.010, 0.004, -0.003, 0.002, -0.006, 0.005, -0.001, 0.003, 0.008];

fn exact(n: i64, a: i64) -> f64 {
    var v = 0.01 * (a + 1);
    for t in elements::dim { v += G[3 * a + t] * elements::coords[3 * n + t]; }
    return v;
}

// Fixes the boundary nodes at the linear field, solves, and gives the
// largest error of the nodal values and of the cell stresses (or fluxes),
// relative to the size of the field.
fn patch(d: i64, ord: i64, heat: bool) -> f64, f64 {
    let cells = 3;
    elements::mesh_box(d, ord, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, cells, cells, cells);
    distort(cells);
    let E = 1000.0;
    let nu = 0.3;
    if heat { elements::heat(2.0, 0.0, 1.0); } else { elements::elastic_plane_strain(E, nu, 1.0); }
    if d == 3 && !heat { elements::elastic(E, nu); }
    let nc = elements::ncomp;
    for n in elements::node_count() {
        if !on_boundary(n) { continue; }
        for a in nc {
            elements::fixed[nc * n + a] = 1;
            elements::prescribed[nc * n + a] = exact(n, a);
        }
    }
    if !elements::solve() {
        failures++;
        print("FAIL patch ", kind_name(d, ord), ": CG does not converge");
    }
    var err_u = 0.0;
    var size_u = 0.0;
    for n in elements::node_count() {
        for a in nc {
            err_u = max(err_u, abs(elements::sol[nc * n + a] - exact(n, a)));
            size_u = max(size_u, abs(exact(n, a)));
        }
    }
    // The exact stress (or flux) of the linear field.
    let lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = E / (2.0 * (1.0 + nu));
    var want: f64[9] = [0.0; 9];
    var tr = 0.0;
    for a in d { tr += G[4 * a]; }
    for i in 3 {
        for j in 3 {
            if heat {
                if i == 0 && j < d { want[j] = -2.0 * G[j]; }
            } else if i < d && j < d {
                want[3 * i + j] = mu * (G[3 * i + j] + G[3 * j + i]);
                if i == j { want[3 * i + j] += lambda * tr; }
            } else if i == 2 && j == 2 {
                want[8] = lambda * tr;
            }
        }
    }
    var size_s = 0.0;
    for v in want { size_s = max(size_s, abs(v)); }
    var err_s = 0.0;
    let what = if heat { elements::FIELD_FLUX } else { elements::FIELD_STRESS };
    let comps = if heat { 3 } else { 9 };
    for c in elements::cell_count() {
        // A point inside the cell that is not a quadrature point.
        let xi: f64[3] = [0.15, 0.2, if d == 3 { 0.25 } else { 0.0 }];
        for k in comps {
            let s = elements::cell_field(c, xi, what, k);
            err_s = max(err_s, abs(s - want[k]));
        }
    }
    return err_u / size_u, err_s / size_s;
}

fn test_patches() {
    for heat in [false, true] {
        for d in [2, 3] {
            for ord in [1, 2] {
                let eu, es = patch(d, ord, heat);
                let name = str(kind_name(d, ord), if heat { " heat" } else { " elastic" });
                print("patch ", name, ": ", elements::dof_count(), " dofs, ", elements::cg_iterations(),
                      " CG iterations, error of the nodal values ", eu, ", of the ",
                      if heat { "flux " } else { "stress " }, es);
                // Round-off of a solve with CG to 1e-12 on a few hundred dofs.
                if !(eu <= 1.0e-10) || !(es <= 1.0e-9) {
                    failures++;
                    print("FAIL patch ", name, ": the linear field is not exact");
                }
            }
        }
    }
}

// --- the plate with a hole -------------------------------------------------------

// A quarter of a plate 2W × 2W with a hole of radius A, under a tension S
// along x on x = W. The Kirsch solution for an infinite plate gives
// sigma_xx = 3 S at the top of the hole (0, A) and sigma_yy = -S at its side
// (A, 0). The plate has the width 2W = 80 A, so the finite width changes
// the factor by less than 0.1 % (Heywood: Kt = 3.002 for d/W = 0.025).
let A = 1.0;
let W = 40.0;
let S = 100.0;

// Solves the quarter plate of half width w on the graded structured mesh
// and gives the factor sigma_xx / S at the top of the hole.
fn kirsch(d: i64, ord: i64, w: f64, size: f64, fine: f64) -> f64 {
    let th = 1.0;
    elements::mesh_hole(d, ord, 0.0, 0.0, 0.0, w, w, th, A, size, fine);
    if d == 2 { elements::elastic_plane_stress(200000.0, 0.3, th); } else { elements::elastic(200000.0, 0.0); }
    assert(elements::plane(elements::X, 0.0));
    elements::fix(elements::X, 0.0);
    assert(elements::plane(elements::Y, 0.0));
    elements::fix(elements::Y, 0.0);
    if d == 3 {
        assert(elements::plane(elements::Z, 0.0));
        elements::fix(elements::Z, 0.0);
    }
    assert(elements::plane(elements::X, w));
    elements::traction(S, 0.0, 0.0);
    let name = str("kirsch ", kind_name(d, ord), " w ", w, " size ", size, " fine ", fine);
    if !elements::solve() {
        failures++;
        print("FAIL ", name, ": CG does not converge");
    }
    let z = 0.5 * th;
    let top, found = elements::stress_at(elements::X, elements::X, 0.0, A, z);
    let side, found2 = elements::stress_at(elements::Y, elements::Y, A, 0.0, z);
    assert(elements::plane(elements::X, 0.0));
    let rx, ry, rz = elements::reaction();
    let fx, fy, fz = elements::applied_load();
    print(name, ": ", elements::element_count(), " cells, ", elements::dof_count(), " dofs, ",
          elements::cg_iterations(), " CG iterations, ", elements::solve_seconds(), " s");
    print("  sigma_xx at the top ", top, " (factor ", top / S, "), sigma_yy at the side ", side,
          ", reaction x ", rx, ", load x ", fx);
    expect_near(str(name, " load"), fx, S * w * th, 1.0e-12);
    expect_near(str(name, " equilibrium"), -rx, fx, 1.0e-8);
    return top / S;
}

// The quarter plate of half width w on the unstructured triangles of
// src/mesher/, all of the size h, and the factor at the top of the hole.
// For order 2, the midpoint node of each edge on the hole goes onto the
// arc, so the cell edge there is curved.
fn kirsch_mesher(ord: i64, w: f64, h: f64) -> f64 {
    mesher::clear_shapes();
    let g = mesher::subtract(mesher::rect(0.0, 0.0, w, w), mesher::circle(0.0, 0.0, A));
    if !mesher::mesh_2d(g, h) {
        failures++;
        print("FAIL kirsch mesher: the mesher does not finish");
        return 0.0;
    }
    var vertices: i64[>..] = [];
    vertices.append(runtime::cell_nodes);
    elements::mesh_simplices(2, ord, runtime::coords, vertices);
    let tol_r = 1.0e-9 * w;
    if ord == 2 {
        for c in elements::cell_count() {
            for e in 3 {
                let a = elements::cell_node(c, elements::edge_node(2, e, 0));
                let b = elements::cell_node(c, elements::edge_node(2, e, 1));
                guard abs(radius(a) - A) <= tol_r && abs(radius(b) - A) <= tol_r else { continue; }
                let m = elements::cell_node(c, 3 + e);
                let r = radius(m);
                for t in 2 { elements::coords[3 * m + t] *= A / r; }
            }
        }
    }
    let th = 1.0;
    elements::elastic_plane_stress(200000.0, 0.3, th);
    assert(elements::plane(elements::X, 0.0));
    elements::fix(elements::X, 0.0);
    assert(elements::plane(elements::Y, 0.0));
    elements::fix(elements::Y, 0.0);
    assert(elements::plane(elements::X, w));
    elements::traction(S, 0.0, 0.0);
    let name = str("kirsch mesher ", kind_name(2, ord), " w ", w, " size ", h);
    if !elements::solve() {
        failures++;
        print("FAIL ", name, ": CG does not converge");
    }
    let top, found = elements::stress_at(elements::X, elements::X, 0.0, A, 0.0);
    assert(elements::plane(elements::X, 0.0));
    let rx, ry, rz = elements::reaction();
    let fx, fy, fz = elements::applied_load();
    print(name, ": ", elements::element_count(), " cells, ", elements::dof_count(), " dofs, ",
          elements::cg_iterations(), " CG iterations, ", elements::solve_seconds(), " s");
    print("  sigma_xx at the top ", top, " (factor ", top / S, "), reaction x ", rx, ", load x ", fx);
    if !found {
        failures++;
        print("FAIL ", name, ": no cell holds the top of the hole");
    }
    expect_near(str(name, " load"), fx, S * w * th, 1.0e-12);
    expect_near(str(name, " equilibrium"), -rx, fx, 1.0e-8);
    return top / S;
}

fn radius(n: i64) -> f64 {
    let x = elements::coords[3 * n];
    let y = elements::coords[3 * n + 1];
    return sqrt(x * x + y * y);
}

// The gate: the factor against Kirsch on the plate of half width W = 40 A.
// The tolerance of P1 is 5 % because its stress is constant in each cell.
// Then the mesh of src/mesher/ on a plate of half width 10 A, which is
// small enough for cells of one size: its factor against the structured
// mesh of the same plate. The finite plate raises the factor to about
// 3.08, so Kirsch is not the reference there. The tolerance is 2 %: the
// stress at the top is the mean of the 2 cells at that node, which differ
// by 2.3 % at h = 0.2 A.
fn test_kirsch() {
    expect_near("kirsch tri3 factor at the top", kirsch(2, 1, W, 4.0, 0.05), 3.0, 0.05);
    expect_near("kirsch tri6 factor at the top", kirsch(2, 2, W, 4.0, 0.1), 3.0, 0.01);
    expect_near("kirsch tet10 factor at the top", kirsch(3, 2, W, 4.0, 0.1), 3.0, 0.01);
    let small = 10.0 * A;
    let structured = kirsch(2, 2, small, 1.0, 0.05);
    expect_near("kirsch mesher tri6 against the structured mesh", kirsch_mesher(2, small, 0.2), structured, 0.02);
}

// --- the heated bar --------------------------------------------------------------

// The section a × a of a long bar with a uniform source q and its sides at
// T = 0. The series of Timoshenko and Goodier (1970), §109, for the centre:
// T = (q b² / 2k) (1 - (32/π³) Σ_n odd (-1)^((n-1)/2) / (n³ cosh(nπ/2))),
// b = a/2. femahoy (validation/cases/grid-heated-bar) gives 14.734270656 K
// for a = 0.02 m, q = 1e7 W/m³, k = 20 W/(m K).
let BAR_A = 0.02;
let BAR_Q = 1.0e7;
let BAR_K = 20.0;
let BAR_CENTRE = 14.734270656;

fn series_centre() -> f64 {
    let b = 0.5 * BAR_A;
    var sum = 0.0;
    var n = 1;
    while n < 41 {
        let sgn = if (n / 2) % 2 == 0 { 1.0 } else { -1.0 };
        sum += sgn / (n * n * n * (0.5 * (exp(0.5 * n * PI) + exp(-0.5 * n * PI))));
        n += 2;
    }
    return BAR_Q * b * b / (2.0 * BAR_K) * (1.0 - 32.0 / (PI * PI * PI) * sum);
}

fn heated_bar(d: i64, ord: i64, cells: i64, tol: f64) {
    elements::mesh_box(d, ord, 0.0, 0.0, 0.0, BAR_A, BAR_A, 0.25 * BAR_A, cells, cells, max(cells / 4, 1));
    elements::heat(BAR_K, BAR_Q, 1.0);
    for axis in 2 {
        for v in [0.0, BAR_A] {
            assert(elements::plane(axis, v));
            elements::fix(0, 0.0);
        }
    }
    let name = str("heated bar ", kind_name(d, ord), " ", cells, " cells");
    if !elements::solve() {
        failures++;
        print("FAIL ", name, ": CG does not converge");
    }
    let t, found = elements::value_at(0, 0.5 * BAR_A, 0.5 * BAR_A, 0.125 * BAR_A);
    let fx, fy, fz = elements::applied_load();
    assert(elements::plane(elements::X, 0.0));
    assert(elements::also_plane(elements::X, BAR_A));
    for v in [0.0, BAR_A] { assert(elements::also_plane(elements::Y, v)); }
    let heat_out, ry, rz = elements::reaction();
    print(name, ": ", elements::dof_count(), " dofs, ", elements::cg_iterations(), " CG iterations, centre T ", t,
          " K (series ", series_centre(), "), source ", fx, ", reaction ", heat_out);
    expect_near(str(name, " centre"), t, BAR_CENTRE, tol);
    expect_near(str(name, " heat balance"), -heat_out, fx, 1.0e-8);
}

fn test_heated_bar() {
    expect_near("the series against femahoy", series_centre(), BAR_CENTRE, 1.0e-9);
    heated_bar(2, 1, 32, 0.01);
    heated_bar(2, 2, 16, 0.001);
    heated_bar(3, 2, 8, 0.002);
}

// --- a body force ------------------------------------------------------------------

// A column of length L under its own weight rho*g, held at z = 0 with
// rollers on x = 0 and y = 0, with nu = 0: u_z = -(rho g / E)(L z - z²/2).
// P2 cells hold the quadratic, so the tip value is exact. The base takes
// the weight rho g V.
fn test_body_force() {
    let L = 4.0;
    let w = 1.0;
    let E = 1000.0;
    let rg = 2.0;
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, w, w, L, 2, 2, 4);
    elements::elastic(E, 0.0);
    elements::body_force(0.0, 0.0, -rg);
    for axis in 3 {
        assert(elements::plane(axis, 0.0));
        elements::fix(axis, 0.0);
    }
    if !elements::solve() {
        failures++;
        print("FAIL body force: CG does not converge");
    }
    let uz, found = elements::value_at(elements::Z, 0.5 * w, 0.5 * w, L);
    let fx, fy, fz = elements::applied_load();
    assert(elements::plane(elements::Z, 0.0));
    let rx, ry, rz = elements::reaction();
    print("body force: tip u_z ", uz, ", weight ", fz, ", base reaction ", rz);
    expect_near("body force tip", uz, -rg * L * L / (2.0 * E), 1.0e-8);
    expect_near("body force weight", fz, -rg * w * w * L, 1.0e-10);
    expect_near("body force base", rz, rg * w * w * L, 1.0e-8);
}

// --- the scale of the check `equilibrium` -----------------------------------------

// One nonzero temperature at all fixed nodes and no source: the reactions
// are round-off of K u, and they are inside the limit. A solution that is
// not the solve of K u = f (one free node moved by 10 K) is outside it.
fn test_equilibrium_scale() {
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, BAR_A, BAR_A, 0.25 * BAR_A, 8, 8, 2);
    elements::heat(BAR_K, 0.0, 1.0);
    for axis in 2 {
        for v in [0.0, BAR_A] {
            assert(elements::plane(axis, v));
            elements::fix(0, 273.15);
        }
    }
    if !elements::solve() {
        failures++;
        print("FAIL equilibrium scale: CG does not converge");
    }
    for moved in [false, true] {
        if moved {
            var n = 0;
            while elements::fixed[n] != 0 { n++; }
            elements::sol[n] += 10.0;
        }
        let rx, ry, rz = elements::total_reaction();
        let fx, fy, fz = elements::applied_load();
        let r = abs(rx + fx);
        let limit = 1.0e-6 * elements::load_scale();
        print("equilibrium scale, ", if moved { "one node moved" } else { "solved" }, ": ", r, " W, limit ", limit, " W");
        if (r <= limit) == moved {
            failures++;
            print("FAIL equilibrium scale: ", if moved { "a real imbalance passes" } else { "round-off fails" });
        }
    }
}

// --- skein run ------------------------------------------------------------------

// `skein run` of the gate models: each model holds its own `expect`, so a
// run that does not pass exits with a status other than 0. `./build.sh
// test` builds build/skein first. The runs write .skein under build/.
let MODELS = [
    "examples/femahoy/grid_kirsch_hole.skn",
    "examples/elements/kirsch_plate_tri.skn",
    "examples/elements/heated_bar_tet.skn",
    "examples/elements/heated_bar_template_tet.skn",
    "examples/elements/conduction_slab_tet.skn",
    "examples/elements/uniform_temperature_tet.skn",
    "examples/elements/column_weight_tet.skn",
];

fn test_runs() {
    for m in MODELS {
        print("skein run ", m, ":");
        var c: u8[>..] = [];
        format(c, "cd build && ./skein run ../", m);
        c.push(0);
        if system(c[0]) != 0 {
            failures++;
            print("FAIL skein run ", m);
        }
    }
}

// --- the submodel cut ------------------------------------------------------------

// A box with the linear field of the patch test is the source of a smaller
// box with other cells. The smaller box takes the field on its faces, so it
// has the same field inside, and the stress resultants on its faces match.
fn test_submodel() {
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    elements::elastic(1000.0, 0.3);
    for n in elements::node_count() {
        if !on_boundary(n) { continue; }
        for a in 3 {
            elements::fixed[3 * n + a] = 1;
            elements::prescribed[3 * n + a] = exact(n, a);
        }
    }
    if !elements::solve() || !elements::write_state("build", "submodel_test") {
        failures++;
        print("FAIL submodel: the source does not solve or save");
        return;
    }
    let lo: f64[3] = [0.2, 0.1, 0.3];
    let hi: f64[3] = [0.7, 0.8, 0.9];
    elements::mesh_box(3, 2, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 3, 3, 3);
    elements::elastic(1000.0, 0.3);
    guard elements::read_state("build", "submodel_test", 1.0, 1.0, 1.0) else {
        failures++;
        print("FAIL submodel: read_state");
        return;
    };
    elements::plane(elements::X, lo[0]);
    for t in 3 {
        if t > 0 { elements::also_plane(t, lo[t]); }
        elements::also_plane(t, hi[t]);
    }
    let outside = elements::fix_state();
    if outside != 0 || !elements::solve() {
        failures++;
        print("FAIL submodel: ", outside, " cut nodes outside the source, or no solve");
        return;
    }
    var err = 0.0;
    for n in elements::node_count() {
        for a in 3 { err = max(err, abs(elements::sol[3 * n + a] - exact(n, a))); }
    }
    let d, scale = elements::cut_difference();
    print("submodel tet10: ", elements::dof_count(), " dofs, error of the nodal values ", err,
          ", of the cut resultants ", d, " of ", scale);
    if !(err <= 1.0e-10) || !(scale > 0.0) || !(d <= 1.0e-9 * scale) {
        failures++;
        print("FAIL submodel: the field or the cut resultants differ");
    }
    if elements::read_state("build", "no_such_step", 1.0, 1.0, 1.0) {
        failures++;
        print("FAIL submodel: read_state of a step with no state");
    }
}

fn main() {
    test_patches();
    test_kirsch();
    test_heated_bar();
    test_body_force();
    test_equilibrium_scale();
    test_submodel();
    test_runs();
    if failures > 0 { exit(1); }
    print("elements: all tests pass");
}
