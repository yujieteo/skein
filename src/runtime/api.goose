import std;
namespace runtime;
import .fem;
import .mesh;
import .solid;
import .solve;
import .vtu;

let X = 0;
let Y = 1;
let Z = 2;

fn mesh_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
            nx: i64, ny: i64, nz: i64) {
    build_box(HEX27, x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

fn mesh_box_incompatible(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
                         nx: i64, ny: i64, nz: i64) {
    build_box(HEX8I, x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

fn mesh_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                nr: i64, nt: i64, nz: i64) {
    build_annulus(cx, cy, r0, r1, t0, t1, z0, z1, nr, nt, nz);
}

fn mesh_union(boxes: f64[:], points: f64[:], size: f64) {
    build_grid(grid_lines(boxes, points, 0, size), grid_lines(boxes, points, 1, size), grid_lines(boxes, points, 2, size), boxes);
}

fn write_mesh(dir: const u8[:]) -> bool { return save_mesh(dir); }

fn read_mesh(dir: const u8[:]) -> bool {
    guard load_mesh(dir) else { return false; };
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
    return true;
}

fn write_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return save_array(dir, str(step, ".fixed"), fixed) && save_array(dir, str(step, ".load"), load) &&
           save_array(dir, str(step, ".disp"), disp) && save_array(dir, str(step, ".row_start"), row_start) &&
           save_array(dir, str(step, ".col"), col) && save_array(dir, str(step, ".val"), val);
}

fn read_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return load_array(dir, str(step, ".fixed"), fixed) && load_array(dir, str(step, ".load"), load) &&
           load_array(dir, str(step, ".disp"), disp) && load_array(dir, str(step, ".row_start"), row_start) &&
           load_array(dir, str(step, ".col"), col) && load_array(dir, str(step, ".val"), val);
}

fn cells_for(length: f64, size: f64) -> i64 {
    guard length > 0.0 && size > 0.0 else { abort("cells_for: length > 0 and size > 0"); }
    return max(ceil(length / size - 1.0e-9) as i64, 1);
}

fn material(E: f64, nu: f64) {
    guard E > 0.0 && nu > -1.0 && nu < 0.5 else { abort("material: E > 0 and -1 < nu < 0.5"); }
    young = E;
    poisson = nu;
    cell_young.clear();
    cell_poisson.clear();
    cell_strain.clear();
    for c in cell_count() {
        cell_young.push(E);
        cell_poisson.push(nu);
        cell_strain.push(0.0);
    }
}

fn material_box(E: f64, nu: f64, strain: f64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> i64 {
    guard E > 0.0 && nu > -1.0 && nu < 0.5 && cell_young.len == cell_count() else { abort("material_box: E > 0, -1 < nu < 0.5, after material"); }
    return set_box_material(E, nu, strain, [x0, y0, z0], [x1, y1, z1]);
}

fn thermal_load() {
    guard cell_strain.len == cell_count() && load.len == dof_total() && dofs_per_node == 3 else { abort("thermal_load: a material and the dofs of a solid"); }
    add_thermal_load();
}

fn load_spread() -> f64 {
    guard load.len == dof_total() else { abort("load_spread: the loads of the dofs"); }
    return positive_load();
}

fn point_set(name: const u8[:], x: f64, y: f64, z: f64) -> bool {
    guard node_count() > 0 else { abort("point_set: a mesh"); }
    return add_point_set(name, [x, y, z]);
}

fn fix(set: const u8[:]) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fix: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 else { return false; }
    fix_set(s, 0, dofs_per_node);
    return true;
}

fn fix_component(set: const u8[:], c: i64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fix_component: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 && c >= 0 && c < dofs_per_node else { return false; }
    fix_set(s, c, c + 1);
    return true;
}

fn traction(set: const u8[:], tx: f64, ty: f64, tz: f64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("traction: a mesh and its dofs"); }
    let s = find_set(set, true);
    guard s >= 0 else { return false; }
    return add_set_traction(s, tx, ty, tz, 0.0);
}

fn normal_traction(set: const u8[:], s: f64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("normal_traction: a mesh and its dofs"); }
    let fs = find_set(set, true);
    guard fs >= 0 else { return false; }
    return add_set_traction(fs, 0.0, 0.0, 0.0, s);
}

fn face_area(set: const u8[:]) -> f64 {
    guard node_count() > 0 else { abort("face_area: a mesh"); }
    let s = find_set(set, true);
    guard s >= 0 else { return 0.0; }
    return set_area(s);
}

fn force_point(x: f64, y: f64, z: f64, fx: f64, fy: f64, fz: f64) {
    add_point_force(x, y, z, fx, fy, fz);
}

fn clear_loads() {
    clear_load_vector();
}

fn linear_solver(solver: const u8[:], preconditioner: const u8[:]) {
    solver_name = [];
    solver_name.append(solver[..min(solver.len, 64)]);
    preconditioner_name = [];
    preconditioner_name.append(preconditioner[..min(preconditioner.len, 64)]);
}

fn solve_static() -> bool {
    return solve_linear_static();
}

fn displacement_at(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    return interpolate_u(x, y, z);
}

fn stress_at(x: f64, y: f64, z: f64) -> f64, f64, f64, f64, f64, f64 {
    return stress_at_point(x, y, z);
}

fn reaction(set: const u8[:]) -> f64, f64, f64 {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("reaction: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 else { return 0.0, 0.0, 0.0; }
    return set_reaction(s);
}

fn sets_share_nodes(a: const u8[:], b: const u8[:]) -> bool {
    guard node_count() > 0 else { abort("sets_share_nodes: a mesh"); }
    let sa = find_set(a, false);
    let sb = find_set(b, false);
    guard sa >= 0 && sb >= 0 else { return false; }
    return sets_meet(sa, sb);
}

fn total_reaction() -> f64, f64, f64 {
    return fixed_reaction();
}

fn applied_load() -> f64, f64, f64 {
    return load_sum();
}

fn write_vtu(path: const u8[:]) -> bool {
    return write_vtu_file(path);
}

fn dof_count() -> i64 { return dof_total(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
