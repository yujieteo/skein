import std;
import math;
import os;
import elements from "../elements";
import mesher from "../mesher";
import runtime from "../runtime";
import solvers from "../solvers";

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
    fn setenv(name: const u8&, value: const u8&, overwrite: i32) -> i32;
    fn unsetenv(name: const u8&) -> i32;
}

var failures = 0;

let PATCH_NODAL_TOLERANCE = 1.0e-10;
let PATCH_STRESS_TOLERANCE = 1.0e-9;
let KIRSCH_TRI3_TOLERANCE = 0.05;
let KIRSCH_P2_TOLERANCE = 0.01;
let KIRSCH_MESHER_TOLERANCE = 0.02;
let HARMONIC_STATIC_TOLERANCE = 1.0e-4;
let HARMONIC_PEAK_TOLERANCE = 0.01;

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

let G: f64[9] = [0.010, 0.004, -0.003, 0.002, -0.006, 0.005, -0.001, 0.003, 0.008];

fn exact(n: i64, a: i64) -> f64 {
    var v = 0.01 * (a + 1);
    for t in elements::dim { v += G[3 * a + t] * elements::coords[3 * n + t]; }
    return v;
}

var written_out = false;

fn isotropic_stiffness(E: f64, nu: f64) -> f64[36] {
    let lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = E / (2.0 * (1.0 + nu));
    var c: f64[36] = [0.0; 36];
    for i in 3 {
        for j in 3 { c[6 * i + j] = lambda + (if i == j { 2.0 * mu } else { 0.0 }); }
    }
    for i in 3..6 { c[7 * i] = 2.0 * mu; }
    return c;
}

fn patch(d: i64, ord: i64, heat: bool) -> f64, f64 {
    let cells = 3;
    elements::mesh_box(d, ord, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, cells, cells, cells);
    distort(cells);
    let E = 1000.0;
    let nu = 0.3;
    if heat { elements::heat(2.0, 0.0, 1.0); } else { elements::elastic_plane_strain(E, nu, 1.0); }
    if d == 3 && !heat { elements::elastic(E, nu); }
    if written_out && !heat { elements::anisotropic(isotropic_stiffness(E, nu), [0.0; 6]); }
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
                if !(eu <= PATCH_NODAL_TOLERANCE) || !(es <= PATCH_STRESS_TOLERANCE) {
                    failures++;
                    print("FAIL patch ", name, ": the linear field is not exact");
                }
            }
        }
    }
}

fn test_written_out() {
    written_out = true;
    for d in [2, 3] {
        for ord in [1, 2] {
            let eu, es = patch(d, ord, false);
            print("patch ", kind_name(d, ord), " written out: error of the nodal values ", eu, ", of the stress ", es);
            if !(eu <= PATCH_NODAL_TOLERANCE) || !(es <= PATCH_STRESS_TOLERANCE) {
                failures++;
                print("FAIL patch ", kind_name(d, ord), " written out: the linear field is not exact");
            }
        }
    }
    written_out = false;
}

let A = 1.0;
let W = 40.0;
let S = 100.0;

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

fn test_kirsch() {
    expect_near("kirsch tri3 factor at the top", kirsch(2, 1, W, 4.0, 0.05), 3.0, KIRSCH_TRI3_TOLERANCE);
    expect_near("kirsch tri6 factor at the top", kirsch(2, 2, W, 4.0, 0.1), 3.0, KIRSCH_P2_TOLERANCE);
    expect_near("kirsch tet10 factor at the top", kirsch(3, 2, W, 4.0, 0.1), 3.0, KIRSCH_P2_TOLERANCE);
    let small = 10.0 * A;
    let structured = kirsch(2, 2, small, 1.0, 0.05);
    expect_near("kirsch mesher tri6 against the structured mesh", kirsch_mesher(2, small, 0.2), structured, KIRSCH_MESHER_TOLERANCE);
}

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

fn test_direct_fallback() {
    elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
    elements::elastic(1000.0, 0.49999999999999);
    assert(elements::plane(elements::X, 0.0));
    elements::fix(0, 0.0);
    elements::fix(1, 0.0);
    elements::fix(2, 0.0);
    assert(elements::plane(elements::X, 1.0));
    elements::flux(1.0);
    if !elements::solve() {
        failures++;
        print("FAIL direct fallback: the solve does not converge");
    }
    if elements::direct {
        failures++;
        print("FAIL direct fallback: the factor passes");
    }
}

fn warm_bar() -> bool {
    solvers::warm_solves = 0;
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, BAR_A, BAR_A, 0.25 * BAR_A, 8, 8, 2);
    elements::heat(BAR_K, 0.0, 1.0);
    assert(elements::plane(elements::X, 0.0));
    elements::fix(0, 273.15);
    assert(elements::plane(elements::X, BAR_A));
    elements::fix(0, 373.15);
    return elements::solve();
}

fn c_string(s: u8[:]) -> u8[>..] {
    var c: u8[>..] = [];
    c.append(s);
    c.push(0);
    return c;
}

fn test_warm_start() {
    let name = c_string("SKEIN_WARM");
    let dir = c_string("build/warm_elements");
    let clean = c_string("rm -rf build/warm_elements");
    let make = c_string("mkdir -p build/warm_elements");
    system(clean[0]);
    system(make[0]);
    setenv(name[0], dir[0], 1);
    let cold = warm_bar();
    let cold_its = elements::cg_iterations();
    let warm = warm_bar();
    let warm_its = elements::cg_iterations();
    unsetenv(name[0]);
    system(clean[0]);
    print("warm start: ", cold_its, " CG iterations cold, ", warm_its, " warm");
    if !cold || !warm || warm_its > 2 {
        failures++;
        print("FAIL warm start: the warm solve of a model with no free load does not converge at once");
    }
}

let MODELS = [
    "examples/femahoy/grid_kirsch_hole.skn",
    "examples/femahoy/birth_death_column.skn",
    "examples/femahoy/buckling_column.skn",
    "examples/femahoy/fracture_centre_crack.skn",
    "examples/femahoy/submodel_kirsch_hole_global.skn",
    "examples/femahoy/submodel_kirsch_hole.skn",
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

fn test_stages() {
    let E = 1000.0;
    let E2 = 3000.0;
    mesher::clear_shapes();
    let right = mesher::rect(1.0, 0.0, 2.0, 1.0);
    let left = mesher::rect(0.0, 0.0, 1.0, 1.0);
    fn stretch() {
        elements::clear_bcs();
        elements::plane(elements::Y, 0.0);
        elements::fix(elements::ALL, 0.0);
        elements::plane(elements::Y, 1.0);
        elements::fix(elements::Y, 0.1);
    }
    fn top_force(what: u8[:], want: f64) {
        if !elements::solve() {
            failures++;
            print("FAIL stages ", what, ": no solve");
            return;
        }
        elements::plane(elements::Y, 1.0);
        let rx, ry, rz = elements::reaction();
        let tx, ty, tz = elements::total_reaction();
        let fx, fy, fz = elements::applied_load();
        let sl, found_l = elements::stress_at(elements::Y, elements::Y, 0.5, 0.5, 0.0);
        let sr, found_r = elements::stress_at(elements::Y, elements::Y, 1.5, 0.5, 0.0);
        print("stages ", what, ": force at y = 1 ", ry, ", sigma_yy left ", sl, ", right ", sr, ", out of balance ", ty + fy);
        if abs(ry - want) > 1.0e-8 * E || abs(ty + fy) > 1.0e-8 * E {
            failures++;
            print("FAIL stages ", what, ": the force at y = 1 is ", ry, ", not ", want, ", or the reactions and loads do not balance");
        }
    }
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 4, 2, 1);
    elements::elastic_plane_strain(E, 0.0, 1.0);
    if !elements::kill(right) { failures++; print("FAIL stages: kill finds no cell"); }
    stretch();
    top_force("step 1, the left block", 0.1 * E);
    stretch();
    if !elements::birth(right) { failures++; print("FAIL stages: birth finds no cell"); }
    top_force("step 2, the right block born strain-free", 0.1 * E);
    let sr, found = elements::stress_at(elements::Y, elements::Y, 1.5, 0.5, 0.0);
    if abs(sr) > 1.0e-8 * E { failures++; print("FAIL stages: the born block has the stress ", sr); }
    stretch();
    elements::kill(left);
    top_force("step 3, the left block dead", 0.0);
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 4, 2, 1);
    elements::elastic_plane_strain(E, 0.0, 1.0);
    if !elements::material_in(right, E2, 0.0, 0.0, 0.0, 0.0) { failures++; print("FAIL stages: material_in finds no cell"); }
    stretch();
    top_force("two materials", 0.1 * (E + E2));
}

fn test_birth_origin() {
    mesher::clear_shapes();
    let right = mesher::rect(1.0, 0.0, 2.0, 1.0);
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 4, 2, 1);
    elements::elastic_plane_strain(1000.0, 0.0, 1.0);
    elements::kill(right);
    elements::plane(elements::Y, 0.0);
    elements::fix(elements::ALL, 0.0);
    elements::plane(elements::Y, 1.0);
    elements::fix(elements::Y, 0.1);
    if !elements::solve() { failures++; print("FAIL birth origin: no solve"); return; }
    elements::clear_bcs();
    elements::plane(elements::Y, 0.0);
    elements::fix(elements::ALL, 0.0);
    elements::plane(elements::Y, 1.0);
    elements::fix(elements::Y, 0.1);
    elements::birth(right);
    if !elements::solve() { failures++; print("FAIL birth origin: no solve after the birth"); return; }
    let near, found_n = elements::value_at(elements::Y, 1.0, 1.0, 0.0);
    let far, found_f = elements::value_at(elements::Y, 2.0, 1.0, 0.0);
    print("birth origin: u.y ", near, " at the old face, ", far, " at the new face");
    if abs(near - 0.1) > 1.0e-9 || abs(far) > 1.0e-9 {
        failures++;
        print("FAIL birth origin: a born node shows the change since its birth: 0.1 and 0 expected");
    }
}

fn test_quarter_point() {
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 1);
    elements::elastic_plane_strain(1000.0, 0.3, 1.0);
    let moved = elements::quarter_point(0.0, 0.0);
    if moved < 2 { failures++; print("FAIL quarter point: ", moved, " edges moved"); return; }
    for n in elements::node_count() {
        let r = sqrt(elements::node_x(n, 0) * elements::node_x(n, 0) + elements::node_x(n, 1) * elements::node_x(n, 1));
        elements::sol[2 * n] = sqrt(r);
        elements::sol[2 * n + 1] = 0.0;
    }
    var worst = 0.0;
    for n in elements::node_count() {
        let x = elements::node_x(n, 0);
        let y = elements::node_x(n, 1);
        let r = sqrt(x * x + y * y);
        let v, found = elements::value_at(elements::X, x, y, 0.0);
        worst = max(worst, abs(v - sqrt(r)));
        if !found { failures++; print("FAIL quarter point: node ", n, " is not found"); }
    }
    let on_edge, found_e = elements::value_at(elements::X, 0.125, 0.0, 0.0);
    print("quarter point: ", moved, " edges moved, largest error at the nodes ", worst, ", at r = 0.125: ", on_edge, " (", sqrt(0.125), ")");
    if worst > 1.0e-9 || abs(on_edge - sqrt(0.125)) > 1.0e-9 {
        failures++;
        print("FAIL quarter point: a node does not show its value");
    }
}

fn test_eigen() {
    let E = 210000.0;
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, 100.0, 10.0, 5.0, 20, 2, 1);
    elements::elastic(E, 0.3);
    elements::set_density(7.85e-9);
    elements::clear_bcs();
    elements::plane(elements::X, 0.0);
    elements::fix(elements::ALL, 0.0);
    let modal = elements::modal(3);
    let f1 = elements::eigenvalue(1);
    print("modal: ", elements::dof_count(), " dofs, ", elements::cg_iterations(), " LOBPCG iterations, ",
          elements::solve_seconds(), " s, f1 ", f1, " Hz, f2 ", elements::eigenvalue(2), " Hz, f3 ", elements::eigenvalue(3), " Hz");
    if !modal || abs(f1 - 417.758) > 0.01 * 417.758 {
        failures++;
        print("FAIL modal: f1 is ", f1, " Hz, not 417.758 Hz within 1 %, or LOBPCG did not converge");
    }
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, 10.0, 10.0, 200.0, 2, 2, 20);
    elements::elastic(E, 0.3);
    elements::clear_bcs();
    elements::plane(elements::Z, 0.0);
    elements::fix(elements::ALL, 0.0);
    elements::plane(elements::Z, 200.0);
    elements::traction(0.0, 0.0, -0.01);
    if !elements::solve() { failures++; print("FAIL buckling: the preload does not solve"); }
    elements::clear_bcs();
    elements::plane(elements::Z, 0.0);
    elements::fix(elements::ALL, 0.0);
    let buckled = elements::buckling(2);
    let p1 = elements::eigenvalue(1);
    print("buckling: ", elements::dof_count(), " dofs, ", elements::cg_iterations(), " LOBPCG iterations, ",
          elements::solve_seconds(), " s, load factors ", p1, ", ", elements::eigenvalue(2));
    if !buckled || abs(p1 - 10794.88) > 0.01 * 10794.88 {
        failures++;
        print("FAIL buckling: the load factor is ", p1, ", not 10794.88 within 1 %, or LOBPCG did not converge");
    }
}

fn hex_bar(curve: bool) {
    elements::mesh_hex_box(0.0, 0.0, 0.0, BAR_A, BAR_A, 0.25 * BAR_A, 16, 16, 4);
    elements::heat(BAR_K, BAR_Q, 1.0);
    if curve { elements::conductivity_curve(0.0, 100.0, [BAR_K, BAR_K]); }
    for axis in 2 {
        for v in [0.0, BAR_A] {
            assert(elements::plane(axis, v));
            elements::fix(0, 0.0);
        }
    }
    let name = str("hex heated bar", if curve { " with k(T)" } else { "" });
    if !elements::solve() {
        failures++;
        print("FAIL ", name, ": the solve does not converge");
    }
    let t, found = elements::value_at(0, 0.5 * BAR_A, 0.5 * BAR_A, 0.125 * BAR_A);
    print(name, ": ", elements::dof_count(), " dofs, centre T ", t, " K (series ", BAR_CENTRE, ")");
    if !found { failures++; print("FAIL ", name, ": the centre is not found"); }
    expect_near(name, t, BAR_CENTRE, 0.01);
}

fn hex_cooling() {
    let h = 0.02;
    elements::mesh_hex_box(0.0, 0.0, 0.0, 0.001, 0.001, h, 1, 1, 20);
    elements::heat(1.0, 0.0, 1.0);
    elements::set_capacity(1.0e6);
    elements::initial_temperature(100.0);
    for v in [0.0, h] {
        assert(elements::plane(elements::Z, v));
        elements::fix(0, 0.0);
    }
    let n = elements::heat_begin(20.0, 0.1, 1.0);
    for i in n { assert(elements::heat_advance()); }
    elements::heat_end();
    let t, found = elements::value_at(0, 0.0005, 0.0005, 0.5 * h);
    print("hex cooling slab: ", n, " time steps, centre T ", t, " K");
    if !found { failures++; print("FAIL hex cooling slab: the centre is not found"); }
    expect_near("hex cooling slab", t, 77.231160686, 0.01);
}

fn hex_modal() {
    elements::mesh_hex_box(0.0, 0.0, 0.0, 100.0, 10.0, 5.0, 20, 2, 2);
    elements::elastic(210000.0, 0.3);
    elements::incompatible(true);
    elements::set_density(7.85e-9);
    elements::clear_bcs();
    elements::plane(elements::X, 0.0);
    elements::fix(elements::ALL, 0.0);
    let done = elements::modal(1);
    let f1 = elements::eigenvalue(1);
    elements::incompatible(false);
    print("hex modal: ", elements::dof_count(), " dofs, f1 ", f1, " Hz");
    if !done { failures++; print("FAIL hex modal: LOBPCG does not converge"); }
    expect_near("hex modal", f1, 417.758, 0.03);
}

fn test_hex() {
    hex_bar(false);
    hex_bar(true);
    hex_cooling();
    hex_modal();
}

fn wave_at(x: f64) -> f64[2] {
    let re, f1 = elements::value_at(0, x, 0.05, 0.0);
    let im, f2 = elements::value_at(1, x, 0.05, 0.0);
    return [re, im];
}

fn expect_wave(what: u8[:], got: f64[2], want: f64[2], tol: f64) {
    let err = sqrt((got[0] - want[0]) * (got[0] - want[0]) + (got[1] - want[1]) * (got[1] - want[1]));
    if !(err <= tol) {
        failures++;
        print("FAIL ", what, ": got ", got[0], " + ", got[1], "i, want ", want[0], " + ", want[1], "i within ", tol);
    }
}

fn test_solid_harmonic() {
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, 100.0, 10.0, 5.0, 20, 2, 1);
    elements::elastic(210000.0, 0.3);
    elements::set_density(7.85e-9);
    elements::body_force(0.0, 0.0, -7.85e-9 * 1000.0);
    elements::clear_bcs();
    elements::plane(elements::X, 0.0);
    elements::fix(elements::ALL, 0.0);
    if !elements::solve() { failures++; print("FAIL solid harmonic: the static step does not solve"); }
    elements::plane(elements::X, 100.0);
    let still = elements::mean_value(elements::Z);
    elements::set_harmonic(2.0 * PI, 0.02);
    if !elements::solve() { failures++; print("FAIL solid harmonic: the step at 1 Hz does not solve"); }
    elements::plane(elements::X, 100.0);
    expect_near("solid harmonic, static limit", elements::harmonic_mean(0, elements::Z), still / (1.0 + 4.0 * 0.02 * 0.02), HARMONIC_STATIC_TOLERANCE);
    var best = 0.0;
    var best_f = 0.0;
    var f = 380.0;
    while f <= 460.0 {
        elements::set_harmonic(2.0 * PI * f, 0.02);
        if !elements::solve() { failures++; print("FAIL solid harmonic: a step does not solve"); }
        elements::plane(elements::X, 100.0);
        let re = elements::harmonic_mean(0, elements::Z);
        let im = elements::harmonic_mean(1, elements::Z);
        let size = sqrt(re * re + im * im);
        if size > best {
            best = size;
            best_f = f;
        }
        f += 4.0;
    }
    print("solid harmonic: static ", still, " mm, peak ", best, " mm at ", best_f, " Hz");
    expect_near("solid harmonic, peak frequency", best_f, 417.758, HARMONIC_PEAK_TOLERANCE);
}

fn test_harmonic() {
    let L = 1.1;
    let k = 2.0 * PI / 0.25;
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, L, 0.1, 0.0, 88, 4, 1);
    elements::wave(1.0);
    elements::wave_terms(k * k, 0.0, 0.0, 0.0, 0.0);
    elements::plane(elements::X, 0.0);
    elements::fix(0, 1.0);
    elements::fix(1, 0.0);
    elements::plane(elements::X, L);
    elements::fix(0, 0.0);
    elements::fix(1, 1.0);
    if !elements::solve() { failures++; print("FAIL harmonic: GMRES does not converge"); }
    let x = 0.37;
    let s = sin(k * L);
    expect_wave("harmonic strip", wave_at(x), [sin(k * (L - x)) / s, sin(k * x) / s], 1.0e-3);
}

fn test_pml() {
    let k = 2.0 * PI / 0.25;
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, 1.5, 0.1, 0.0, 120, 2, 1);
    elements::wave(1.0);
    elements::wave_terms(k * k, 0.0, 0.0, 0.0, 0.0);
    elements::pml(mesher::rect(1.0, -1.0, 2.0, 1.0), -10.0, -10.0, -10.0, 1.0, 10.0, 10.0, 0.5);
    elements::plane(elements::X, 0.0);
    elements::fix(0, 1.0);
    elements::fix(1, 0.0);
    if !elements::solve() { failures++; print("FAIL pml: GMRES does not converge"); }
    for i in 9 {
        let x = 0.1 + 0.1 * i as f64;
        expect_wave("pml reflection", wave_at(x), [cos(k * x), -sin(k * x)], 2.0e-3);
    }
}

fn test_bloch() {
    let L = 1.0;
    let kap = 7.0;
    let theta = 0.5 * PI;
    elements::mesh_box(2, 2, 0.0, 0.0, 0.0, L, 0.1, 0.0, 40, 2, 1);
    elements::wave(1.0);
    elements::wave_terms(kap * kap, 0.0, 1.0, 0.0, 0.0);
    elements::plane(elements::X, 0.0);
    elements::bloch_low();
    elements::plane(elements::X, L);
    if !elements::bloch(L, 0.0, 0.0, theta) { failures++; print("FAIL bloch: the high face does not tie to the low face"); }
    if !elements::solve() { failures++; print("FAIL bloch: GMRES does not converge"); }
    expect_wave("bloch cell", wave_at(0.3), [0.0247168, -0.0093576], 1.0e-5);
}

fn edge_exact(p: f64[3]) -> f64[6] {
    if elements::dim == 2 { return [sin(PI * p[1]), sin(PI * p[0]), 0.0, 0.0, 0.0, 0.0]; }
    return [sin(PI * p[1]) * sin(PI * p[2]), sin(PI * p[2]) * sin(PI * p[0]), sin(PI * p[0]) * sin(PI * p[1]), 0.0, 0.0, 0.0];
}

fn edge_exact_curl(p: f64[3]) -> f64[3] {
    if elements::dim == 2 { return [0.0, 0.0, PI * cos(PI * p[0]) - PI * cos(PI * p[1])]; }
    let s: f64[3] = [sin(PI * p[0]), sin(PI * p[1]), sin(PI * p[2])];
    let c: f64[3] = [cos(PI * p[0]), cos(PI * p[1]), cos(PI * p[2])];
    return [PI * s[0] * (c[1] - c[2]), PI * s[1] * (c[2] - c[0]), PI * s[2] * (c[0] - c[1])];
}

fn edge_errors() -> f64, f64 {
    var eu = 0.0;
    var ec = 0.0;
    let d = elements::dim;
    for c in elements::cell_count() {
        for q in elements::rich_count() {
            let xi, w = elements::rich_point(q);
            var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
            for r in d {
                l[r + 1] = xi[r];
                l[0] -= xi[r];
            }
            var p: f64[3] = [0.0; 3];
            for v in d + 1 {
                for t in 3 { p[t] += l[v] * elements::node_x(elements::cell_node(c, v), t); }
            }
            let u = elements::edge_cell_value(c, l, false);
            let cu = elements::edge_cell_value(c, l, true);
            let dv = w * elements::cell_geometry(c, [0.0; 3]);
            let ue = edge_exact(p);
            let ce = edge_exact_curl(p);
            for t in 3 {
                eu += (u[t] - ue[t]) * (u[t] - ue[t]) * dv;
                ec += (cu[t] - ce[t]) * (cu[t] - ce[t]) * dv;
            }
        }
    }
    return sqrt(eu), sqrt(ec);
}

fn edge_box(d: i64, n: i64, c: f64, p: i64) -> f64, f64 {
    elements::mesh_box(d, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, if d == 3 { n } else { 1 });
    distort(n);
    elements::edge_order(p);
    elements::edge(1.0, c, 0.0, 1.0);
    let k = if d == 2 { PI * PI } else { 2.0 * PI * PI };
    elements::current() { p =>
        let u = edge_exact(p);
        [(k + c) * u[0], (k + c) * u[1], (k + c) * u[2], 0.0, 0.0, 0.0]
    };
    elements::plane(elements::X, 0.0);
    elements::also_plane(elements::X, 1.0);
    elements::also_plane(elements::Y, 0.0);
    elements::also_plane(elements::Y, 1.0);
    if d == 3 {
        elements::also_plane(elements::Z, 0.0);
        elements::also_plane(elements::Z, 1.0);
    }
    elements::tangent() { p => edge_exact(p) };
    if c == 0.0 && elements::gauge() <= 0 { failures++; print("FAIL edge: the tree gauge holds no edge"); }
    if !elements::solve() { failures++; print("FAIL edge: GMRES does not converge on ", n, " cells"); }
    return edge_errors();
}

fn test_edge() {
    for d in [2, 3] {
        for c in [1.0, 0.0] {
            let n0 = if d == 2 { 8 } else { 4 };
            let eu0, ec0 = edge_box(d, n0, c, 1);
            let eu1, ec1 = edge_box(d, 2 * n0, c, 1);
            let ou = log(eu0 / eu1) / log(2.0);
            let oc = log(ec0 / ec1) / log(2.0);
            let what = str("nedelec ", d, "d c=", c);
            if c != 0.0 && !(ou >= 0.7) { failures++; print("FAIL ", what, ": the L2 order of u is ", ou, ", want 1 - 0.3"); }
            if !(oc >= 0.7) { failures++; print("FAIL ", what, ": the L2 order of curl u is ", oc, ", want 1 - 0.3"); }
            print(what, ": L2 orders ", ou, " (u) and ", oc, " (curl u), errors ", eu1, " and ", ec1);
        }
    }
}

fn test_edge_order2() {
    for d in [2, 3] {
        let n0 = if d == 2 { 4 } else { 2 };
        let eu0, ec0 = edge_box(d, n0, 1.0, 2);
        let eu1, ec1 = edge_box(d, 2 * n0, 1.0, 2);
        let ou = log(eu0 / eu1) / log(2.0);
        let oc = log(ec0 / ec1) / log(2.0);
        let what = str("nedelec order 2 ", d, "d");
        if !(ou >= 1.6) { failures++; print("FAIL ", what, ": the L2 order of u is ", ou, ", want 2 - 0.4"); }
        if !(oc >= 1.6) { failures++; print("FAIL ", what, ": the L2 order of curl u is ", oc, ", want 2 - 0.4"); }
        print(what, ": L2 orders ", ou, " (u) and ", oc, " (curl u), errors ", eu1, " and ", ec1);
    }
    elements::edge_order(1);
}

fn test_edge_modal() {
    for d in [2, 3] {
        let a = 1.0;
        let c = 1.3;
        elements::mesh_box(d, 1, 0.0, 0.0, 0.0, a, 0.8, c, if d == 2 { 10 } else { 4 }, if d == 2 { 7 } else { 3 }, if d == 3 { 5 } else { 1 });
        elements::edge_order(2);
        elements::edge(1.0, 0.0, 0.0, 1.0);
        elements::edge_mass_in(-1, 1.0, 0.0);
        elements::plane(elements::X, 0.0);
        elements::also_plane(elements::X, a);
        elements::also_plane(elements::Y, 0.0);
        elements::also_plane(elements::Y, 0.8);
        if d == 3 {
            elements::also_plane(elements::Z, 0.0);
            elements::also_plane(elements::Z, c);
        }
        elements::tangent() { p => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] };
        let want = if d == 2 { PI } else { PI * sqrt(1.0 / (a * a) + 1.0 / (c * c)) };
        let ok = elements::resonance(3, 0.9 * want) { g, w => [0.0, 0.0] };
        if !ok { failures++; print("FAIL edge modal ", d, "d: no convergence"); continue; }
        let f, fi = elements::modal_frequency(1);
        expect_near(str("edge modal ", d, "d frequency"), f * 2.0 * PI, want, 1.0e-3);
        if abs(fi) > 1.0e-8 * f { failures++; print("FAIL edge modal ", d, "d: a lossless mode has Im f = ", fi); }
        let f2, f2i = elements::modal_frequency(2);
        let f3, f3i = elements::modal_frequency(3);
        if !(f2 > f * (1.0 + 1.0e-6) && f3 > f2 * (1.0 + 1.0e-6)) { failures++; print("FAIL edge modal ", d, "d: the modes are not different and in the sequence of f"); }
        if d == 2 { expect_near("edge modal 2d second frequency", f2 * 2.0 * PI, PI / 0.8, 2.0e-3); }
        print("edge modal ", d, "d: omega ", f * 2.0 * PI, ", want ", want, ", second ", f2 * 2.0 * PI);
    }
    test_edge_lossy();
    test_edge_lossy_fixed();
}

fn test_edge_lossy() {
    var w_found = 0.0;
    var f_found = 0.0;
    for pass in 2 {
        elements::mesh_box(2, 1, 0.0, 0.0, 0.0, 1.0, 0.8, 0.0, 10, 7, 1);
        elements::edge_order(2);
        elements::edge(1.0, 0.0, 0.0, 1.0);
        elements::edge_mass_in(-1, 1.0, 0.0);
        elements::plane(elements::X, 1.0);
        elements::also_plane(elements::Y, 0.0);
        elements::also_plane(elements::Y, 0.8);
        elements::tangent() { p => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] };
        elements::plane(elements::X, 0.0);
        elements::impedance_group(0);
        let ok = if pass == 0 {
            elements::resonance(1, 3.0) { g, w => [0.0, 0.2 * sqrt(w)] }
        } else {
            elements::resonance(1, 3.0) { g, w => [0.0, 0.2 * sqrt(w_found)] }
        };
        if !ok { failures++; print("FAIL edge lossy: no convergence"); return; }
        let f, fi = elements::modal_frequency(1);
        if pass == 0 {
            w_found = 2.0 * PI * f;
            f_found = f;
            if !(abs(fi) > 1.0e-3 * f) { failures++; print("FAIL edge lossy: Im f = ", fi, " for a wall with loss"); }
            print("edge lossy: f ", f, " + ", fi, " i");
        } else {
            expect_near("edge lossy fixed point", f, f_found, 1.0e-9);
        }
    }
}

fn test_edge_lossy_fixed() {
    elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 0.8, 1.3, 4, 3, 5);
    elements::edge_order(2);
    elements::edge(1.0, 0.0, 0.0, 1.0);
    elements::edge_mass_in(-1, 1.0, 0.0);
    elements::plane(elements::X, 1.0);
    elements::also_plane(elements::Y, 0.0);
    elements::also_plane(elements::Y, 0.8);
    elements::also_plane(elements::Z, 0.0);
    elements::also_plane(elements::Z, 1.3);
    elements::tangent() { p => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] };
    elements::plane(elements::X, 0.0);
    elements::impedance_group(0);
    let ok = elements::resonance(1, 3.0) { g, w => [0.0, 0.2 * sqrt(w)] };
    if !ok { failures++; print("FAIL edge lossy fixed: no convergence"); return; }
    var top = 0.0;
    var at_fixed = 0.0;
    for d in elements::fixed.len {
        let v = abs(elements::edge_modes[d]);
        top = max(top, v);
        if elements::fixed[d] != 0 { at_fixed = max(at_fixed, v); }
    }
    if !(at_fixed <= 1.0e-12 * top) { failures++; print("FAIL edge lossy fixed: the mode is ", at_fixed, " on a fixed dof, max ", top); }
    print("edge lossy fixed: max on a fixed dof ", at_fixed, ", max ", top);
}

fn dipole_field(p: f64[3], k: f64, eta: f64) -> f64[6], f64[6] {
    let r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    let ct = p[2] / r;
    let st = sqrt(p[0] * p[0] + p[1] * p[1]) / r;
    let cp = if st > 0.0 { p[0] / (r * st) } else { 1.0 };
    let sp = if st > 0.0 { p[1] / (r * st) } else { 0.0 };
    let kr = k * r;
    let ph: f64[2] = [cos(kr), -sin(kr)];
    let b1: f64[2] = [1.0, -1.0 / kr];
    let b2: f64[2] = [1.0 - 1.0 / (kr * kr), -1.0 / kr];
    let er = elements::cmul(elements::cmul([eta * ct / (2.0 * PI * r * r), 0.0], b1), ph);
    let eth = elements::cmul(elements::cmul([0.0, eta * k * st / (4.0 * PI * r)], b2), ph);
    let hph = elements::cmul(elements::cmul([0.0, k * st / (4.0 * PI * r)], b1), ph);
    let rh: f64[3] = [st * cp, st * sp, ct];
    let th: f64[3] = [ct * cp, ct * sp, -st];
    let phv: f64[3] = [-sp, cp, 0.0];
    var e: f64[6] = [0.0; 6];
    var h: f64[6] = [0.0; 6];
    for t in 3 {
        for part in 2 {
            e[3 * part + t] = er[part] * rh[t] + eth[part] * th[t];
            h[3 * part + t] = hph[part] * phv[t];
        }
    }
    return e, h;
}

fn test_far_field() {
    let k = 2.0 * PI;
    let eta = 376.730313;
    let m = 40;
    elements::clear_samples();
    for axis in 3 {
        for side in [-1.0, 1.0] {
            for i in m {
                for j in m {
                    var p: f64[3] = [0.0; 3];
                    p[axis] = 0.5 * side;
                    p[(axis + 1) % 3] = (i + 0.5) / m - 0.5;
                    p[(axis + 2) % 3] = (j + 0.5) / m - 0.5;
                    var n: f64[3] = [0.0; 3];
                    n[axis] = side;
                    let e, h = dipole_field(p, k, eta);
                    elements::far_p.append(p);
                    elements::far_n.append(n);
                    elements::far_w.push(1.0 / (m * m));
                    elements::far_e.append(e);
                    elements::far_h.append(h);
                }
            }
        }
    }
    let u = elements::far_intensity([1.0, 1.0, 0.0], k, eta);
    let power = elements::surface_power();
    expect_near("far field intensity", u, eta * k * k / (32.0 * PI * PI), 1.0e-2);
    expect_near("far field directivity", 4.0 * PI * u / power, 1.5, 1.0e-2);
    let u45 = elements::far_intensity([1.0, 0.0, 1.0], k, eta);
    expect_near("far field at 45 degrees", u45, 0.5 * eta * k * k / (32.0 * PI * PI), 1.0e-2);
}

fn darcy_errors() -> f64, f64 {
    var ep = 0.0;
    var eq = 0.0;
    for c in elements::cell_count() {
        for qp in elements::rich_count() {
            let xi, w = elements::rich_point(qp);
            let p = elements::point_of(c, xi);
            let v = elements::mixed_cell_value(c, p);
            let dv = w * elements::cell_geometry(c, [0.0; 3]);
            let pe = sin(PI * p[0]) * sin(PI * p[1]);
            let qe: f64[2] = [-PI * cos(PI * p[0]) * sin(PI * p[1]), -PI * sin(PI * p[0]) * cos(PI * p[1])];
            ep += (v[3] - pe) * (v[3] - pe) * dv;
            eq += ((v[0] - qe[0]) * (v[0] - qe[0]) + (v[1] - qe[1]) * (v[1] - qe[1])) * dv;
        }
    }
    return sqrt(ep), sqrt(eq);
}

fn darcy_square(n: i64) -> f64, f64 {
    elements::mesh_box(2, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, n, n, 1);
    distort(n);
    elements::darcy(1.0, 0.0, 1.0);
    elements::darcy_source() { p => 2.0 * PI * PI * sin(PI * p[0]) * sin(PI * p[1]) };
    elements::plane(elements::X, 0.0);
    elements::also_plane(elements::X, 1.0);
    elements::also_plane(elements::Y, 0.0);
    elements::also_plane(elements::Y, 1.0);
    if !elements::pressure(0.0) { failures++; print("FAIL darcy: the boundary holds no facet"); }
    if !elements::solve() { failures++; print("FAIL darcy: CG does not converge on ", n, " cells"); }
    return darcy_errors();
}

fn test_darcy() {
    let ep0, eq0 = darcy_square(8);
    let ep1, eq1 = darcy_square(16);
    let op = log(ep0 / ep1) / log(2.0);
    let oq = log(eq0 / eq1) / log(2.0);
    if !(op >= 0.7) { failures++; print("FAIL raviart_thomas: the L2 order of p is ", op, ", want 1 - 0.3"); }
    if !(oq >= 0.7) { failures++; print("FAIL raviart_thomas: the L2 order of q is ", oq, ", want 1 - 0.3"); }
    print("raviart_thomas 2d: L2 orders ", op, " (p) and ", oq, " (q), errors ", ep1, " and ", eq1);
    elements::mesh_box(2, 1, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 8, 4, 1);
    elements::darcy(1.0, 0.0, 1.0);
    elements::plane(elements::X, 0.0);
    elements::normal_flux(-3.0);
    elements::plane(elements::X, 2.0);
    elements::pressure(0.0);
    if !elements::solve() { failures++; print("FAIL darcy strip: CG does not converge"); }
    let p, found = elements::darcy_at(3, 0.5, 0.5, 0.0);
    expect_near("darcy strip pressure", p, 3.0 * 1.5, 0.05);
    expect_near("darcy strip flux", elements::boundary_flux(), 3.0, 1.0e-9);
}

fn rel_diff(a: f64[:], b: f64[:]) -> f64 {
    var d = 0.0;
    var m = 0.0;
    for x, i in a {
        d = max(d, abs(x - b[i]));
        m = max(m, abs(x));
    }
    return d / max(m, 1.0);
}

fn neo_hooke(c: elements::Hd[3]) -> elements::Hd {
    let lnj = elements::hd_ln(elements::hd_mul(elements::hd_mul(c[0], c[1]), c[2])) * 0.5;
    return (c[0] + c[1] + c[2] - elements::hd(3.0)) * 40.0 - lnj * 80.0 + elements::hd_mul(lnj, lnj) * 60.0;
}

fn neo_hooke_s(cm: f64[9]) -> f64[6] {
    let ci = elements::mat_inv(cm);
    let lnj = 0.5 * log(elements::mat_det(cm));
    return elements::to_voigt((elements::EYE - ci) * 80.0 + ci * (120.0 * lnj));
}

fn spectral_s(cm: f64[9], s: f64[6]&, t: f64[36]&) {
    let c, e = elements::sym_eigen(cm);
    var w1: f64[3] = [0.0; 3];
    var w2: f64[9] = [0.0; 9];
    elements::energy_derivatives(c, 3, w1, w2, neo_hooke);
    elements::energy_stress(c, e, 3, w1, w2, s, t);
}

fn tangent_column(cm: f64[9], w: i64) -> f64[6] {
    let h = 1.0e-6;
    var up = cm;
    var down = cm;
    let i = elements::VOIGT_I[w];
    let j = elements::VOIGT_J[w];
    let dc = if i == j { 2.0 * h } else { h };
    up[3 * i + j] += dc;
    down[3 * i + j] -= dc;
    if i != j {
        up[3 * j + i] += dc;
        down[3 * j + i] -= dc;
    }
    return (neo_hooke_s(up) - neo_hooke_s(down)) * (0.5 / h);
}

fn test_solid_points() {
    let a: f64[9] = [4.0, 1.0, 0.5, 1.0, 3.0, -0.2, 0.5, -0.2, 1.0];
    let l, e = elements::sym_eigen(a);
    var worst = 0.0;
    for k in 3 {
        for i in 3 {
            var ae = 0.0;
            for j in 3 { ae += a[3 * i + j] * e[3 * j + k]; }
            worst = max(worst, abs(ae - l[k] * e[3 * i + k]));
        }
    }
    print("solid: Jacobi residual ", worst);
    if !(worst < 1.0e-12) { failures++; print("FAIL solid Jacobi: residual ", worst); }
    let cs: f64[9][2] = [[1.3, 0.2, -0.1, 0.2, 0.9, 0.15, -0.1, 0.15, 1.1], elements::EYE];
    for cm, k in cs {
        var s: f64[6] = [0.0; 6];
        var t: f64[36] = [0.0; 36];
        spectral_s(cm, s, t);
        let exact = neo_hooke_s(cm);
        let es = rel_diff(exact, s);
        var et = 0.0;
        for w in 6 {
            let col = tangent_column(cm, w);
            var tw: f64[6] = [0.0; 6];
            for v in 6 { tw[v] = t[6 * v + w]; }
            et = max(et, rel_diff(col, tw));
        }
        print("solid: neo-Hooke at C ", k, ": error of S ", es, ", of the tangent ", et);
        if !(es < 1.0e-12 && et < 1.0e-6) { failures++; print("FAIL solid neo-Hooke at C ", k); }
    }
    test_j2_points();
}

fn test_j2_points() {
    let m = elements::J2 { lame: 115384.6, shear: 76923.1, yield0: 250.0, hard: 1000.0 };
    let eps: f64[9] = [0.004, 0.001, 0.0, 0.001, -0.002, 0.0005, 0.0, 0.0005, 0.0];
    var ep: f64[6] = [0.0; 6];
    var alpha = 0.0;
    var s0: f64[6] = [0.0; 6];
    var t: f64[36] = [0.0; 36];
    elements::j2_small(m, eps, ep, alpha, s0, t);
    var et = 0.0;
    for w in 6 {
        let h = 1.0e-8;
        let i = elements::VOIGT_I[w];
        let j = elements::VOIGT_J[w];
        var up = eps;
        var down = eps;
        let de = if i == j { h } else { 0.5 * h };
        up[3 * i + j] += de;
        down[3 * i + j] -= de;
        if i != j {
            up[3 * j + i] += de;
            down[3 * j + i] -= de;
        }
        var sp: f64[6] = [0.0; 6];
        var sm: f64[6] = [0.0; 6];
        var tt: f64[36] = [0.0; 36];
        var ep2: f64[6] = [0.0; 6];
        var a2 = 0.0;
        elements::j2_small(m, up, ep2, a2, sp, tt);
        ep2 = [0.0; 6];
        a2 = 0.0;
        elements::j2_small(m, down, ep2, a2, sm, tt);
        var tw: f64[6] = [0.0; 6];
        for v in 6 { tw[v] = t[6 * v + w]; }
        let col = (sp - sm) * (0.5 / h);
        et = max(et, rel_diff(col, tw));
    }
    var cpi: f64[6] = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
    var af = 0.0;
    var tau: f64[6] = [0.0; 6];
    let f: f64[9] = [1.3, 0.1, 0.0, 0.0, 0.85, 0.05, 0.0, 0.0, 0.9];
    let p = elements::j2_finite(m, f, cpi, af, tau);
    let det = elements::mat_det(elements::from_voigt(cpi));
    print("solid: J2 alpha ", alpha, ", error of the tangent ", et, "; finite strain alpha ", af, ", det Cp⁻¹ ", det);
    if !(alpha > 0.0 && et < 1.0e-5) { failures++; print("FAIL solid J2 tangent: ", et); }
    if !(af > 0.0 && abs(det - 1.0) < 1.0e-12) { failures++; print("FAIL solid J2 finite: det Cp⁻¹ ", det); }
}

fn g_x(p: f64[3]) -> f64 { return 0.3 * p[0] + 0.1 * p[1]; }
fn g_y(p: f64[3]) -> f64 { return -0.1 * p[1] + 0.05 * p[2]; }
fn g_z(p: f64[3]) -> f64 { return 0.02 * p[0] - 0.08 * p[2]; }

fn test_contact_points() {
    elements::clear_contact();
    elements::ct_n.append([0.0, 0.0, 1.0]);
    elements::ct_g0.push(0.001);
    elements::ct_law.append([1000.0, 500.0, 0.3]);
    for k in 3 {
        elements::ct_t.push(0.0);
        elements::ct_r.push(0.0);
        elements::ct_trial_t.push(0.0);
        elements::ct_trial_r.push(0.0);
    }
    elements::ct_p.push(0.0);
    var d: f64[9] = [0.0; 9];
    let r: f64[3] = [1.0e-4, 0.0, -0.003];
    let f = elements::point_traction(0, r, d);
    expect_near("contact stick t", f[0], 0.05, 1.0e-12);
    expect_near("contact stick p", f[2], -2.0, 1.0e-12);
    var worst = 0.0;
    var dd: f64[9] = [0.0; 9];
    for b in 3 {
        var up = r;
        var down = r;
        up[b] += 1.0e-7;
        down[b] -= 1.0e-7;
        let fu = elements::point_traction(0, up, dd);
        let fd = elements::point_traction(0, down, dd);
        for a in 3 { worst = max(worst, abs((fu[a] - fd[a]) / 2.0e-7 - d[3 * a + b]) / 1000.0); }
    }
    print("contact: stick tangent error ", worst);
    if !(worst < 1.0e-9) { failures++; print("FAIL contact stick tangent: ", worst); }
    let s = elements::point_traction(0, [0.01, 0.005, -0.003], d);
    expect_near("contact slip |t|", sqrt(s[0] * s[0] + s[1] * s[1]), 0.6, 1.0e-12);
    expect_near("contact slip direction", s[1] / s[0], 0.5, 1.0e-12);
    let o = elements::point_traction(0, [0.01, 0.0, 0.0], d);
    if !(o[0] == 0.0 && o[2] == 0.0 && elements::ct_p[0] == 0.0) { failures++; print("FAIL contact open: ", o[0], " ", o[2]); }
    elements::clear_contact();
}

fn test_solid_energy() {
    elements::mesh_box(3, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
    elements::hyperelastic(1.0);
    for axis in 3 {
        for side in 2 {
            assert(elements::plane(axis, side * 1.0));
            elements::fix_with(elements::X, g_x);
            elements::fix_with(elements::Y, g_y);
            elements::fix_with(elements::Z, g_z);
        }
    }
    let done = elements::solve_energy(4, neo_hooke);
    let f: f64[9] = [1.3, 0.1, 0.0, 0.0, 0.9, 0.05, 0.02, 0.0, 0.92];
    let p = elements::mat_mul(f, elements::from_voigt(neo_hooke_s(elements::mat_mul(elements::mat_t(f), f))));
    assert(elements::plane(elements::X, 1.0));
    let rx, ry, rz = elements::reaction();
    let r: f64[3] = [rx, ry, rz];
    let pe: f64[3] = [p[0], p[3], p[6]];
    let err = rel_diff(pe, r);
    print("solid: neo-Hooke box ", elements::newton_iterations(), " Newton iterations, reaction (", rx, ", ", ry, ", ", rz,
          "), P e_x (", p[0], ", ", p[3], ", ", p[6], ")");
    if !(done && err < 1.0e-8) { failures++; print("FAIL solid neo-Hooke box: converged ", done, ", error ", err); }
}

fn test_solid_plastic() {
    let E = 200000.0;
    let sy = 250.0;
    let H = 1000.0;
    elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
    elements::elastic(E, 0.3);
    elements::plastic(sy, H, false);
    for axis in 3 {
        assert(elements::plane(axis, 0.0));
        elements::fix(axis, 0.0);
    }
    assert(elements::plane(elements::X, 1.0));
    elements::fix(elements::X, 0.01);
    let done = elements::solve_plastic(5);
    let rx, ry, rz = elements::reaction();
    let want = sy + E * H / (E + H) * (0.01 - sy / E);
    print("solid: J2 bar ", elements::newton_iterations(), " Newton iterations, stress ", rx, ", want ", want);
    if !done { failures++; print("FAIL solid J2 bar: no convergence"); }
    expect_near("solid J2 bar", rx, want, 1.0e-8);
}

fn test_quadrature_weights_sum_to_reference_size() {
    let sizes: f64[3] = [1.0, 0.5, 1.0 / 6.0];
    let second_moments: f64[3] = [1.0 / 3.0, 1.0 / 12.0, 1.0 / 60.0];
    for d in [1, 2, 3] {
        var total = 0.0;
        var moment = 0.0;
        for q in elements::quad_count(d) {
            let xi, w = elements::quad_point(d, q);
            total += w;
            moment += w * xi[0] * xi[0];
        }
        expect_near(str("quadrature weights d=", d), total, sizes[d - 1], 1.0e-14);
        expect_near(str("quadrature degree 2 exact d=", d), moment, second_moments[d - 1], 1.0e-14);
    }
}

fn main() {
    test_quadrature_weights_sum_to_reference_size();
    test_solid_points();
    test_solid_energy();
    test_solid_plastic();
    test_contact_points();
    test_patches();
    test_written_out();
    test_kirsch();
    test_heated_bar();
    test_body_force();
    test_equilibrium_scale();
    test_direct_fallback();
    test_submodel();
    test_stages();
    test_birth_origin();
    test_quarter_point();
    test_eigen();
    test_hex();
    test_harmonic();
    test_solid_harmonic();
    test_pml();
    test_bloch();
    test_edge();
    test_edge_order2();
    test_edge_modal();
    test_far_field();
    test_darcy();
    test_warm_start();
    test_runs();
    failures += elements::el2_failures();
    if failures > 0 { exit(1); }
    print("elements: all tests pass");
}
