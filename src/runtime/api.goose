// The API that generated code calls (docs/implementation.md §6). A
// generated program imports the package with
// `import runtime from "<path to src>/runtime";` and calls these functions
// in this sequence:
//
//     runtime::mesh_box(0.0, 0.0, 0.0, 100.0, 10.0, 10.0, 40, 4, 4);
//     runtime::material(210000.0, 0.3);
//     runtime::fix("x0");
//     runtime::fix_component("y0", runtime::Y);
//     runtime::traction("x1", 0.0, 0.0, -1.0);
//     runtime::linear_solver("cg", "jacobi");
//     let ok = runtime::solve_static();
//     let ux, uy, uz = runtime::displacement_at(100.0, 5.0, 5.0);
//     let rx, ry, rz = runtime::reaction("x0");
//     runtime::write_vtu("out.vtu");
//
// The runtime has no units. The generated code converts each value to one
// consistent system, for example mm, N and MPa, before the call.
//
// A support, a load and a reaction name a set of the mesh: a node set for
// `fix` and `reaction`, a face set for `traction` and `face_area`. A
// mesher gives a node set and a face set of the same name for each
// surface. A name that the mesh does not have selects nothing, and the
// call returns false. The box mesher names the faces of the box x0, x1,
// y0, y1, z0 and z1: x0 is the face x = x0, as in
// `boundary(beam, x == 0 mm)`. The annulus mesher names the faces of the
// sector r0, r1, t0, t1, z0 and z1 (r, θ and z).
//
// The mesh is the unstructured mesh of fem.goose. mesh_box and
// mesh_annulus are its first meshers, and mesh_union meshes a union of boxes. The state is global: one mesh, a material for each cell, one set of supports and one
// set of loads. A second load step calls `clear_loads`, applies its loads
// and solves again. A second mesh (a sweep over a dimension) starts again
// at `mesh_box`, which clears all of the state.
import std;
namespace runtime;
import .fem;
import .mesh;
import .solid;
import .solve;
import .vtu;

// The components of a vector, for fix_component.
let X = 0;
let Y = 1;
let Z = 2;

// Meshes the box from (x0, y0, z0) to (x1, y1, z1) with nx × ny × nz
// 27-node hexahedra (triquadratic Lagrange cells). Clears the material,
// the supports, the loads and the solution.
fn mesh_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
            nx: i64, ny: i64, nz: i64) {
    build_box(HEX27, x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

// Meshes the box as mesh_box does, with 8-node hexahedra with incompatible
// modes (C3D8I, hex8i.goose).
fn mesh_box_incompatible(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
                         nx: i64, ny: i64, nz: i64) {
    build_box(HEX8I, x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

// Meshes the sector r0 <= r <= r1, t0 <= θ <= t1, z0 <= z <= z1 of an
// annulus about the axis x = cx, y = cy with nr × nt × nz 27-node
// hexahedra. Clears the state as `mesh_box` does.
fn mesh_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                nr: i64, nt: i64, nz: i64) {
    build_annulus(cx, cy, r0, r1, t0, t1, z0, z1, nr, nt, nz);
}

// Meshes the union of the boxes of `boxes` (x0, y0, z0, x1, y1, z1 for
// each) with 27-node hexahedra of at most `size`, on the grid lines of
// `grid_lines`, with a plane through each point of `points`. A cell of the
// grid is in the mesh if its centre is in a box. The node set and the face
// set x0, x1, y0, y1, z0 and z1 hold the boundary on each bounding plane of
// the grid. The program computes the grid, so that a new value of
// a param changes the values of the call and not its form (the browser
// solve, docs/implementation.md §8). Clears the state as `mesh_box` does.
fn mesh_union(boxes: f64[:], points: f64[:], size: f64) {
    build_grid(grid_lines(boxes, points, 0, size), grid_lines(boxes, points, 1, size), grid_lines(boxes, points, 2, size), boxes);
}

// Writes the mesh and its named sets in the directory `dir`, for the
// other cases of a study that share it (spec §20.4). False if it cannot
// write a file.
fn write_mesh(dir: const u8[:]) -> bool { return save_mesh(dir); }

// Reads the mesh that `write_mesh` wrote in `dir`, in place of a mesh call.
// Clears the state as `mesh_box` does. False if it cannot read it.
fn read_mesh(dir: const u8[:]) -> bool {
    guard load_mesh(dir) else { return false; };
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
    return true;
}

// Writes the solved state of step `step` in the directory `dir`: the
// supports, the loads, the matrix and the solution, for the other cases of
// a study that run only the post-processing (spec §20.4). False if it
// cannot write a file.
fn write_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return save_array(dir, str(step, ".fixed"), fixed) && save_array(dir, str(step, ".load"), load) &&
           save_array(dir, str(step, ".disp"), disp) && save_array(dir, str(step, ".row_start"), row_start) &&
           save_array(dir, str(step, ".col"), col) && save_array(dir, str(step, ".val"), val);
}

// Reads the solved state that `write_solution` wrote for step `step`, in
// place of the supports, the loads and the solve of the step. False if it
// cannot read it.
fn read_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return load_array(dir, str(step, ".fixed"), fixed) && load_array(dir, str(step, ".load"), load) &&
           load_array(dir, str(step, ".disp"), disp) && load_array(dir, str(step, ".row_start"), row_start) &&
           load_array(dir, str(step, ".col"), col) && load_array(dir, str(step, ".val"), val);
}

// The number of cells along `length` for a cell size of at most `size`.
fn cells_for(length: f64, size: f64) -> i64 {
    guard length > 0.0 && size > 0.0 else { abort("cells_for: length > 0 and size > 0"); }
    return max(ceil(length / size - 1.0e-9) as i64, 1);
}

// A linear isotropic material for all cells: Young's modulus E and
// Poisson's ratio nu.
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

// Gives the material E, nu and the isotropic thermal strain `strain`
// (alpha (T - T0)) to each cell whose centre is in the box from
// (x0, y0, z0) to (x1, y1, z1). Call it after `material`. It gives the
// number of these cells.
fn material_box(E: f64, nu: f64, strain: f64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> i64 {
    guard E > 0.0 && nu > -1.0 && nu < 0.5 && cell_young.len == cell_count() else { abort("material_box: E > 0, -1 < nu < 0.5, after material"); }
    return set_box_material(E, nu, strain, [x0, y0, z0], [x1, y1, z1]);
}

// Adds the loads of the thermal strain of each cell: f = ∫ Bᵀ D ε_th dV.
fn thermal_load() {
    guard cell_strain.len == cell_count() && load.len == dof_total() && dofs_per_node == 3 else { abort("thermal_load: a material and the dofs of a solid"); }
    add_thermal_load();
}

// The largest sum of the positive loads along x, y or z: the scale of a
// load that sums to 0, as a thermal load does, for the check equilibrium.
fn load_spread() -> f64 {
    guard load.len == dof_total() else { abort("load_spread: the loads of the dofs"); }
    return positive_load();
}

// Adds the node set `name` of the one node at (x, y, z), for a support at
// a point. False if no node is there.
fn point_set(name: const u8[:], x: f64, y: f64, z: f64) -> bool {
    guard node_count() > 0 else { abort("point_set: a mesh"); }
    return add_point_set(name, [x, y, z]);
}

// Holds all dofs at 0 on the nodes of a node set.
fn fix(set: const u8[:]) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fix: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 else { return false; }
    fix_set(s, 0, dofs_per_node);
    return true;
}

// Holds dof c (X, Y or Z for a solid) at 0 on the nodes of a node set, as
// a plane of symmetry does.
fn fix_component(set: const u8[:], c: i64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fix_component: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 && c >= 0 && c < dofs_per_node else { return false; }
    fix_set(s, c, c + 1);
    return true;
}

// Adds a uniform traction (force per area) (tx, ty, tz) on a face set.
fn traction(set: const u8[:], tx: f64, ty: f64, tz: f64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("traction: a mesh and its dofs"); }
    let s = find_set(set, true);
    guard s >= 0 else { return false; }
    return add_set_traction(s, tx, ty, tz, 0.0);
}

// Adds the traction s n on a face set, with n the unit normal out of the
// body. An internal pressure p is s = -p.
fn normal_traction(set: const u8[:], s: f64) -> bool {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("normal_traction: a mesh and its dofs"); }
    let fs = find_set(set, true);
    guard fs >= 0 else { return false; }
    return add_set_traction(fs, 0.0, 0.0, 0.0, s);
}

// The area of a face set, or 0 if the mesh has no face set of that name.
fn face_area(set: const u8[:]) -> f64 {
    guard node_count() > 0 else { abort("face_area: a mesh"); }
    let s = find_set(set, true);
    guard s >= 0 else { return 0.0; }
    return set_area(s);
}

// Adds the force (fx, fy, fz) at the node nearest to (x, y, z).
fn force_point(x: f64, y: f64, z: f64, fx: f64, fy: f64, fz: f64) {
    add_point_force(x, y, z, fx, fy, fz);
}

// Removes all loads and the solution. The mesh, the material and the
// supports stay.
fn clear_loads() {
    clear_load_vector();
}

// Sets the linear solver and the preconditioner of the next solves by
// name. "cg" with "jacobi" is the default. The solvers lane registers the
// others in lane_solvers.goose.
fn linear_solver(solver: const u8[:], preconditioner: const u8[:]) {
    solver_name = [];
    solver_name.append(solver[..min(solver.len, 64)]);
    preconditioner_name = [];
    preconditioner_name.append(preconditioner[..min(preconditioner.len, 64)]);
}

// Assembles the stiffness matrix and solves K u = f by the linear solver
// (CG with a Jacobi preconditioner by default). Returns false if the solver
// does not converge.
fn solve_static() -> bool {
    return solve_linear_static();
}

// The displacement at (x, y, z), interpolated in the cell that holds the
// point. A point outside the 3D cells gives 0, 0, 0.
fn displacement_at(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    return interpolate_u(x, y, z);
}

// The stress at (x, y, z) in the cell that holds the point: xx, yy, zz,
// yz, xz and xy. A point outside the 3D cells gives 0.
fn stress_at(x: f64, y: f64, z: f64) -> f64, f64, f64, f64, f64, f64 {
    return stress_at_point(x, y, z);
}

// The sum of the reaction forces on the fixed nodes of a node set.
fn reaction(set: const u8[:]) -> f64, f64, f64 {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("reaction: a mesh and its dofs"); }
    let s = find_set(set, false);
    guard s >= 0 else { return 0.0, 0.0, 0.0; }
    return set_reaction(s);
}

// Whether the node sets a and b have a node in common. False if the mesh
// has no node set of one of the names.
fn sets_share_nodes(a: const u8[:], b: const u8[:]) -> bool {
    guard node_count() > 0 else { abort("sets_share_nodes: a mesh"); }
    let sa = find_set(a, false);
    let sb = find_set(b, false);
    guard sa >= 0 && sb >= 0 else { return false; }
    return sets_meet(sa, sb);
}

// The sum of the reaction forces on all fixed nodes, each node once.
fn total_reaction() -> f64, f64, f64 {
    return fixed_reaction();
}

// The sum of all applied loads.
fn applied_load() -> f64, f64, f64 {
    return load_sum();
}

// Writes the mesh and the fields of the dofs to a VTU file for ParaView.
fn write_vtu(path: const u8[:]) -> bool {
    return write_vtu_file(path);
}

// Facts about the last solve, for the report.
fn dof_count() -> i64 { return dof_total(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
