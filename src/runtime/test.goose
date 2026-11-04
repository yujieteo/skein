import std;
import math;
import os;
import runtime from "../runtime";
import solvers from "../solvers";

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

let EQUILIBRIUM_TOLERANCE = 1.0e-6;
let FEMAHOY_TIP = -0.190113;

fn euler_bernoulli_tip() -> f64 {
    let inertia = B * H * H * H / 12.0;
    return -P * L * L * L / (3.0 * E * inertia);
}

fn test_cantilever() {
    runtime::mesh_box(0.0, 0.0, 0.0, L, B, H, 40, 4, 4);
    runtime::material(E, NU);
    assert(runtime::fix("x0"));
    assert(runtime::traction("x1", 0.0, 0.0, -P / (B * H)));
    if !runtime::solve_static() {
        failures++;
        print("FAIL cantilever: CG does not converge in ", runtime::cg_iterations(), " iterations");
    }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    let rx, ry, rz = runtime::reaction("x0");
    let fx, fy, fz = runtime::applied_load();
    print("cantilever: 40 x 4 x 4 hex27, ", runtime::dof_count(), " dofs, ",
        runtime::cg_iterations(), " CG iterations, ", runtime::solve_seconds(), " s");
    print("  tip uz ", uz, " mm (femahoy ", FEMAHOY_TIP, ", Euler-Bernoulli ",
        euler_bernoulli_tip(), ")");
    print("  reaction ", rx, " ", ry, " ", rz, " N, applied ", fx, " ", fy, " ", fz, " N");
    expect_near("tip uz against femahoy", uz, FEMAHOY_TIP, 0.02);
    expect_near("tip uz against Euler-Bernoulli", uz, euler_bernoulli_tip(), 0.02);
    expect_near("applied load", fz, -P, 1.0e-12);
    expect_near("tip face area", runtime::face_area("x1"), B * H, 1.0e-12);
    if abs(rz + fz) > EQUILIBRIUM_TOLERANCE * P || abs(rx + fx) > EQUILIBRIUM_TOLERANCE * P || abs(ry + fy) > EQUILIBRIUM_TOLERANCE * P {
        failures++;
        print("FAIL equilibrium: reaction ", rx, " ", ry, " ", rz, ", load ", fx, " ", fy, " ", fz);
    }
    let path = "build/test_runtime_cantilever.vtu";
    if !runtime::write_vtu(path) {
        failures++;
        print("FAIL write_vtu: cannot write ", path);
    }
}

fn test_cantilever_incompatible() {
    runtime::mesh_box_incompatible(0.0, 0.0, 0.0, L, B, H, 40, 4, 4);
    runtime::material(E, NU);
    assert(runtime::fix("x0"));
    assert(runtime::traction("x1", 0.0, 0.0, -P / (B * H)));
    if !runtime::solve_static() {
        failures++;
        print("FAIL hex8i cantilever: CG does not converge in ", runtime::cg_iterations(), " iterations");
    }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, 0.5 * H);
    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(51.25, 0.5 * B, H);
    let bending = P * (L - 51.25) * 0.5 * H / (B * H * H * H / 12.0);
    print("cantilever: 40 x 4 x 4 hex8i, ", runtime::dof_count(), " dofs, ",
        runtime::cg_iterations(), " CG iterations, ", runtime::solve_seconds(), " s");
    print("  tip uz ", uz, " mm (femahoy ", FEMAHOY_TIP, "), top sxx ", sxx, " MPa (Saint-Venant ", bending, ")");
    expect_near("hex8i tip uz against femahoy", uz, FEMAHOY_TIP, 0.01);
    expect_near("hex8i top sxx against Saint-Venant", sxx, bending, 0.015);
}

fn test_patch_incompatible() {
    runtime::mesh_box_incompatible(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2, 2, 2);
    let centre = 1 + 3 * (1 + 3 * 1);
    runtime::coords[3 * centre] += 0.2;
    runtime::coords[3 * centre + 1] -= 0.15;
    runtime::coords[3 * centre + 2] += 0.1;
    runtime::material(E, NU);
    assert(runtime::fix_component("x0", runtime::X));
    assert(runtime::fix_component("y0", runtime::Y));
    assert(runtime::fix_component("z0", runtime::Z));
    assert(runtime::traction("x1", 1.0, 0.0, 0.0));
    if !runtime::solve_static() {
        failures++;
        print("FAIL hex8i patch: CG does not converge");
    }
    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(0.7, 1.3, 0.6);
    let ux, uy, uz = runtime::displacement_at(2.0, 2.0, 2.0);
    let other = abs(syy) + abs(szz) + abs(syz) + abs(sxz) + abs(sxy);
    print("patch: hex8i, sxx ", sxx, " MPa, the other components ", other, " MPa");
    expect_near("hex8i patch sxx", sxx, 1.0, 1.0e-8);
    expect_near("hex8i patch ux", ux, 2.0 / E, 1.0e-8);
    if other > 1.0e-8 {
        failures++;
        print("FAIL hex8i patch: the other stress components are ", other, " MPa, want 0");
    }
}

fn test_second_step() {
    runtime::clear_loads();
    runtime::force_point(L, 0.5 * B, 0.5 * H, 0.0, 0.0, -2.0 * P);
    if runtime::traction("mid", 0.0, 0.0, 1.0) || runtime::fix("mid") {
        failures++;
        print("FAIL traction and fix: mid is not a set of the mesh");
    }
    if !runtime::solve_static() {
        failures++;
        print("FAIL second step: CG does not converge");
    }
    let ux, uy, uz = runtime::displacement_at(L, 0.5 * B, H);
    let rx, ry, rz = runtime::reaction("x0");
    let fx, fy, fz = runtime::applied_load();
    print("second step: point force ", fz, " N, tip uz at the top ", uz, " mm, reaction z ", rz, " N");
    expect_near("second step applied load", fz, -2.0 * P, 1.0e-12);
    expect_near("second step reaction", rz, 2.0 * P, 1.0e-6);
    expect_near("second step tip uz", uz, 2.0 * FEMAHOY_TIP, 0.02);
}

fn linear_residual(u: f64[:], r: f64[:]) {
    runtime::multiply(u, r);
    for d in r.len { r[d] -= runtime::load[d]; }
}

fn linear_jacobian(u: f64[:]) { assert(runtime::assemble()); }

fn test_newton() {
    runtime::clear_loads();
    assert(runtime::traction("x1", 0.0, 0.0, -P / (B * H)));
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

fn test_shared_edge() {
    runtime::mesh_box(0.0, 0.0, 0.0, L, B, H, 10, 1, 1);
    runtime::material(E, NU);
    assert(runtime::fix("x0") && runtime::fix("z0"));
    assert(runtime::traction("x1", 0.0, 0.0, -P / (B * H)));
    if !runtime::solve_static() {
        failures++;
        print("FAIL shared edge: CG does not converge");
    }
    if !runtime::sets_share_nodes("x0", "z0") || runtime::sets_share_nodes("x0", "x1") {
        failures++;
        print("FAIL shared nodes: x0 and z0 share an edge, and x0 and x1 share no node");
    }
    let rx, ry, rz = runtime::total_reaction();
    let fx, fy, fz = runtime::applied_load();
    if abs(rx + fx) > 1.0e-6 * P || abs(ry + fy) > 1.0e-6 * P || abs(rz + fz) > 1.0e-6 * P {
        failures++;
        print("FAIL shared edge: total reaction ", rx, " ", ry, " ", rz, ", load ", fx, " ", fy, " ", fz);
    }
}

fn test_annulus() {
    let a = 100.0;
    let b = 105.0;
    let p = 1.0;
    let A = p * a * a / (b * b - a * a);
    runtime::mesh_annulus(0.0, 0.0, a, b, 0.0, 0.5 * PI, 0.0, 10.0, 2, 48, 2);
    runtime::material(E, NU);
    assert(runtime::fix_component("t1", runtime::X));
    assert(runtime::fix_component("t0", runtime::Y));
    assert(runtime::fix_component("z0", runtime::Z));
    assert(runtime::normal_traction("r0", -p));
    assert(runtime::traction("z1", 0.0, 0.0, A));
    if !runtime::solve_static() {
        failures++;
        print("FAIL annulus: CG does not converge in ", runtime::cg_iterations(), " iterations");
    }
    let r = 102.5;
    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(r, 0.0, 5.0);
    let rx, ry, rz = runtime::reaction("z0");
    let fx, fy, fz = runtime::applied_load();
    print("annulus: 2 x 48 x 2 hex27, ", runtime::dof_count(), " dofs, ",
        runtime::cg_iterations(), " CG iterations, ", runtime::solve_seconds(), " s");
    print("  at r = ", r, ": hoop ", syy, ", axial ", szz, " MPa");
    print("  reaction z ", rz, " N, applied ", fx, " ", fy, " ", fz, " N");
    expect_near("annulus hoop against Lame", syy, A * (1.0 + b * b / (r * r)), 1.0e-3);
    expect_near("annulus axial against Lame", szz, A, 1.0e-3);
    expect_near("annulus axial reaction", rz, -0.25 * PI * a * a * p, 1.0e-6);
    let s45 = r * sqrt(0.5);
    let txx, tyy, tzz, tyz, txz, txy = runtime::stress_at(s45, s45, 5.0);
    expect_near("annulus shear at 45 degrees", -txy, A * b * b / (r * r), 1.0e-2);
    let ox, oy, oz = runtime::displacement_at(0.0, 0.0, 5.0);
    if ox != 0.0 || oy != 0.0 || oz != 0.0 {
        failures++;
        print("FAIL annulus: the axis is not in the mesh, but u there is ", ox, " ", oy, " ", oz);
    }
}

fn test_scalar_field() {
    runtime::mesh_box(0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);
    runtime::dofs_per_node = 1;
    runtime::point_fields.clear();
    runtime::add_point_field("T", 0, 1);
    runtime::clear_dofs();
    for d in runtime::disp.len { runtime::disp[d] = d; }
    let path = "build/test_runtime_scalar.vtu";
    var text: u8[>..] = [];
    if !runtime::write_vtu(path) || !read_file(path, text) {
        failures++;
        print("FAIL scalar field: cannot write or read ", path);
        return;
    }
    let head = "<DataArray type=\"Float64\" Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    let at = find(text, head);
    var values: u8[>..] = [];
    if at >= 0 {
        let rest = text[at + head.len..];
        values.append(rest[..find(rest, "</DataArray>")]);
    }
    let lines = split(values, "\n");
    var good = at >= 0 && lines.len == runtime::node_count() + 1;
    for i in runtime::node_count() {
        if good && lines[i] != str(i as f64) { good = false; }
    }
    if !good {
        failures++;
        print("FAIL scalar field: want one value of T for each of ", runtime::node_count(), " nodes in ", path);
    }
}

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
    assert(runtime::fix("x0"));
    if runtime::solve_static() {
        failures++;
        print("FAIL registry: a tet4 cell with no element routine solves");
    }
}

fn test_thermal_grid() {
    let s = 1.0e-3;
    runtime::mesh_union([0.0, 0.0, 0.0, 2.0, 1.0, 0.5, 0.0, 0.0, 0.5, 2.0, 1.0, 0.6], [], 1.0);
    runtime::material(E, NU);
    let bottom = runtime::material_box(E, NU, s, 0.0, 0.0, 0.0, 2.0, 1.0, 0.5);
    let top = runtime::material_box(E / 3.0, 0.25, s, 0.0, 0.0, 0.5, 2.0, 1.0, 0.6);
    assert(runtime::fix_component("x0", runtime::X));
    assert(runtime::fix_component("y0", runtime::Y));
    assert(runtime::fix_component("z0", runtime::Z));
    if runtime::cell_count() != 4 || bottom != 2 || top != 2 || !runtime::point_set("corner", 2.0, 1.0, 0.6) ||
       runtime::point_set("off", 0.3, 0.0, 0.0) {
        failures++;
        print("FAIL thermal grid: want 4 cells, 2 in each box, and a node at the corner only");
    }
    runtime::thermal_load();
    if !runtime::solve_static() {
        failures++;
        print("FAIL thermal grid: CG does not converge");
    }
    let ux, uy, uz = runtime::displacement_at(2.0, 1.0, 0.6);
    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(1.5, 0.5, 0.55);
    let rx, ry, rz = runtime::total_reaction();
    print("thermal grid: u at the corner ", ux, " ", uy, " ", uz, " mm, sxx in the top ", sxx, " MPa");
    expect_near("thermal grid ux", ux, 2.0 * s, 1.0e-8);
    expect_near("thermal grid uy", uy, s, 1.0e-8);
    expect_near("thermal grid uz", uz, 0.6 * s, 1.0e-8);
    let scale = E * s;
    if abs(sxx) + abs(syy) + abs(szz) > 1.0e-6 * scale || abs(rx) + abs(ry) + abs(rz) > 1.0e-6 * runtime::load_spread() {
        failures++;
        print("FAIL thermal grid: want sigma = 0 and no reaction, got sxx ", sxx, " syy ", syy, " szz ", szz, " reaction ", rx, " ", ry, " ", rz);
    }
}

fn main() {
    test_cantilever();
    test_second_step();
    test_newton();
    test_shared_edge();
    test_annulus();
    test_cantilever_incompatible();
    test_patch_incompatible();
    test_scalar_field();
    test_registry();
    test_thermal_grid();
    if failures > 0 { exit(1); }
    print("runtime: all tests pass");
}
