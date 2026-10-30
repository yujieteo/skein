// The API that generated code calls (docs/implementation.md §6). A
// generated program imports the package with
// `import runtime from "<path to src>/runtime";` and calls these functions
// in this sequence:
//
//     runtime::mesh_box(0.0, 0.0, 0.0, 100.0, 10.0, 10.0, 40, 4, 4);
//     runtime::material(210000.0, 0.3);
//     runtime::fix_face(runtime::X, 0.0);
//     runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -1.0);
//     let ok = runtime::solve_static();
//     let ux, uy, uz = runtime::displacement_at(100.0, 5.0, 5.0);
//     let rx, ry, rz = runtime::reaction_face(runtime::X, 0.0);
//     runtime::write_vtu("out.vtu");
//
// The runtime has no units. The generated code converts each value to one
// consistent system, for example mm, N and MPa, before the call.
//
// The mesh is a box (`mesh_box`) or a sector of an annulus about an axis
// parallel to z (`mesh_annulus`). A face is the part of the boundary where
// the lattice coordinate `axis` is equal to `value`. For a box the lattice
// coordinates are x, y and z, as in `boundary(beam, x == 0 mm)`. For an
// annulus they are r, θ (in radians) and z. A value that is not on a face
// selects no nodes, and the call returns false.
//
// The state is global: one mesh, one material, one set of supports and one
// set of loads. A second load step calls `clear_loads`, applies its loads
// and solves again. A second mesh (a sweep over a dimension) starts again
// at `mesh_box`, which clears all of the state.
import std;
namespace runtime;
import .mesh;
import .solve;
import .vtu;

let X = 0;
let Y = 1;
let Z = 2;
let R = 0;
let THETA = 1;

// Meshes the box from (x0, y0, z0) to (x1, y1, z1) with nx × ny × nz
// 27-node hexahedra (triquadratic Lagrange cells). Clears the material,
// the supports, the loads and the solution.
fn mesh_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
            nx: i64, ny: i64, nz: i64) {
    build_box(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

// Meshes the sector r0 <= r <= r1, t0 <= θ <= t1, z0 <= z <= z1 of an
// annulus about the axis x = cx, y = cy with nr × nt × nz 27-node
// hexahedra. Clears the state as `mesh_box` does.
fn mesh_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                nr: i64, nt: i64, nz: i64) {
    build_annulus(cx, cy, r0, r1, t0, t1, z0, z1, nr, nt, nz);
}

// The number of cells along `length` for a cell size of at most `size`.
fn cells_for(length: f64, size: f64) -> i64 {
    return max(ceil(length / size - 1.0e-9) as i64, 1);
}

// A linear isotropic material for all cells: Young's modulus E and
// Poisson's ratio nu.
fn material(E: f64, nu: f64) {
    young = E;
    poisson = nu;
}

// Holds all 3 components of the displacement at 0 on a face.
fn fix_face(axis: i64, value: f64) -> bool {
    return mark_fixed(axis, value, 0, 3);
}

// Holds component c (X, Y or Z) of the displacement at 0 on a face, as a
// plane of symmetry does.
fn fix_face_component(axis: i64, value: f64, c: i64) -> bool {
    return mark_fixed(axis, value, c, c + 1);
}

// Adds a uniform traction (force per area) (tx, ty, tz) on a face.
fn traction_face(axis: i64, value: f64, tx: f64, ty: f64, tz: f64) -> bool {
    return add_face_load(axis, value, tx, ty, tz, 0.0);
}

// Adds the traction s n on a face, with n the unit normal out of the body.
// An internal pressure p is s = -p.
fn normal_traction_face(axis: i64, value: f64, s: f64) -> bool {
    return add_face_load(axis, value, 0.0, 0.0, 0.0, s);
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

// Assembles the stiffness matrix and solves K u = f by CG with a Jacobi
// preconditioner. Returns false if CG does not converge.
fn solve_static() -> bool {
    return solve_linear_static();
}

// The displacement at (x, y, z), interpolated in the cell that holds the
// point. A point outside the mesh gives 0, 0, 0.
fn displacement_at(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    return interpolate_u(x, y, z);
}

// The stress at (x, y, z) in the cell that holds the point: xx, yy, zz,
// yz, xz and xy. A point outside the mesh gives 0.
fn stress_at(x: f64, y: f64, z: f64) -> f64, f64, f64, f64, f64, f64 {
    return stress_at_point(x, y, z);
}

// The sum of the reaction forces on the fixed dofs of the nodes of a face.
// A node on an edge counts in each face that holds it.
fn reaction_face(axis: i64, value: f64) -> f64, f64, f64 {
    return face_reaction(axis, value);
}

// The sum of all applied loads.
fn applied_load() -> f64, f64, f64 {
    return load_sum();
}

// Writes the mesh and the displacement to a VTU file for ParaView.
fn write_vtu(path: const u8[:]) -> bool {
    return write_vtu_file(path);
}

// Facts about the last solve, for the report.
fn dof_count() -> i64 { return 3 * node_count(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
