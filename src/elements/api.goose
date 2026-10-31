// The API of the simplex elements that generated code calls. A generated
// program imports the package with
// `import elements from "<path to src>/elements";` and calls these
// functions in this sequence:
//
//     elements::mesh_hole(3, 2, 0.0, 0.0, 0.0, 40.0, 40.0, 1.0, 1.0, 4.0, 0.1);
//     elements::elastic(200000.0, 0.0);
//     elements::plane(elements::X, 0.0);
//     elements::fix(elements::X, 0.0);
//     elements::plane(elements::X, 40.0);
//     elements::traction(100.0, 0.0, 0.0);
//     let ok = elements::solve();
//     let sxx, found = elements::stress_at(elements::X, elements::X, 0.0, 1.0, 0.5);
//
// The runtime has no units, like src/runtime: the generated code converts
// each value to one consistent system before the call.
//
// The state is global: one mesh, one physics, one set of fixed dofs and one
// set of loads. A boundary call (`fix`, `traction`, `flux`, `reaction`)
// acts on the region of the last `plane`, `also_plane`, `hole`, `surface`
// or `also_surface` call.
// These give false when the region holds no node.
import std;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .solve;
import .post;
import .csg;
import .source;

let X = 0;
let Y = 1;
let Z = 2;
let ALL = -1;

// Meshes the box from (x0, y0, z0) to (x1, y1, z1) in `d` = 3, or the
// rectangle from (x0, y0) to (x1, y1) in `d` = 2, with nx × ny (× nz)
// lattice cells, each cut into 6 tetrahedra or 2 triangles of order `ord`.
fn mesh_box(d: i64, ord: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
            nx: i64, ny: i64, nz: i64) {
    build_box(d, ord, [x0, y0, z0], [x1, y1, z1], nx, ny, nz);
    physics = -1;
}

// Meshes the box (or the rectangle) less a hole of radius r whose axis is
// z through the corner (x0, y0): the quarter model of a plate with a hole.
// The cells are of size `fine` at the hole and grow to `size`.
fn mesh_hole(d: i64, ord: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             r: f64, size: f64, fine: f64) {
    build_hole(d, ord, [x0, y0, z0], [x1, y1, z1], r, size, fine);
    physics = -1;
}

// Meshes the shape g of src/mesher/ (built with its functions, for
// example mesher::subtract(mesher::box(...), mesher::cylinder(...))) with
// cells of order `ord` and size `size` in dimension d, and smaller at the
// targets of mesher::add_refine. False if the mesher does not finish.
fn mesh_shape(d: i64, ord: i64, g: i64, size: f64) -> bool {
    physics = -1;
    return build_shape(d, ord, g, size);
}

// Meshes with the cells `vertices` of dimension `d` (d + 1 nodes each, in
// positive orientation) on the nodes `xyz`, for example the mesh of
// src/mesher/ in src/runtime/fem.goose. Order 2 adds a node at the
// midpoint of each edge.
fn mesh_simplices(d: i64, ord: i64, xyz: f64[:], vertices: i64[:]) {
    build_simplices(d, ord, xyz, vertices);
    physics = -1;
}

// Writes the mesh in the directory `dir`, for the other cases of a study
// that share it (spec §20.4). False if it cannot write a file.
fn write_mesh(dir: const u8[:]) -> bool { return save_mesh(dir); }

// Reads the mesh that `write_mesh` wrote in `dir`, in place of a mesh call.
// False if it cannot read it.
fn read_mesh(dir: const u8[:]) -> bool {
    physics = -1;
    return load_mesh(dir);
}

// Writes the solved state of step `step` in the directory `dir`: the
// fixed values, the loads, the matrix and the solution, for the other cases
// of a study that run only the post-processing (spec §20.4). False if it
// cannot write a file.
fn write_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return save_array(dir, str(step, ".fixed"), fixed) && save_array(dir, str(step, ".prescribed"), prescribed) &&
           save_array(dir, str(step, ".load"), load) && save_array(dir, str(step, ".sol"), sol) &&
           save_array(dir, str(step, ".row_start"), row_start) && save_array(dir, str(step, ".col"), col) &&
           save_array(dir, str(step, ".val"), val) && save_array(dir, str(step, ".source"), assembled_source);
}

// Reads the solved state that `write_solution` wrote for step `step`, in
// place of the conditions and the solve of the step. False if it cannot
// read it.
fn read_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return load_array(dir, str(step, ".fixed"), fixed) && load_array(dir, str(step, ".prescribed"), prescribed) &&
           load_array(dir, str(step, ".load"), load) && load_array(dir, str(step, ".sol"), sol) &&
           load_array(dir, str(step, ".row_start"), row_start) && load_array(dir, str(step, ".col"), col) &&
           load_array(dir, str(step, ".val"), val) && load_array(dir, str(step, ".source"), assembled_source);
}

// Writes the mesh, the material and the solution of step `step` in the
// directory `dir`, for a `results` of another model (spec §15.1). False if
// it cannot write a file.
fn write_state(dir: const u8[:], step: const u8[:]) -> bool { return save_state(dir, step); }

// Reads the state that `write_state` wrote for step `step` in `dir` as the
// field of a submodel cut. The scales convert a length, a value of the
// unknown and a modulus of the other model into the units of this model.
// False if it cannot read it.
fn read_state(dir: const u8[:], step: const u8[:], length_scale: f64, value_scale: f64, modulus_scale: f64) -> bool {
    return load_source(dir, step, length_scale, value_scale, modulus_scale);
}

// The number of cells along `length` for a cell size of at most `size`.
fn cells_for(length: f64, size: f64) -> i64 {
    guard length > 0.0 && size > 0.0 else { abort("cells_for: length > 0 and size > 0"); }
    return max(ceil(length / size - 1.0e-9) as i64, 1);
}

// Linear isotropic elasticity with Young's modulus E and Poisson's ratio
// nu: 3D, or in 2D plane stress or plane strain with thickness t.
fn elastic(E: f64, nu: f64) { set_elastic(E, nu, SOLID_3D, 1.0); }
fn elastic_plane_stress(E: f64, nu: f64, t: f64) { set_elastic(E, nu, PLANE_STRESS, t); }
fn elastic_plane_strain(E: f64, nu: f64, t: f64) { set_elastic(E, nu, PLANE_STRAIN, t); }

// Sets a uniform body force per volume (bx, by, bz) on the elastic
// physics of the last `elastic` call, for example rho*g.
fn body_force(bx: f64, by: f64, bz: f64) {
    guard physics == ELASTIC else { abort("body_force: after elastic"); }
    volume_force = [bx, by, bz];
}

// Steady heat conduction with conductivity k and a uniform source q per
// volume. In 2D, t is the thickness.
fn heat(k: f64, q: f64, t: f64) { set_heat(k, q, t); }

// Selects the boundary plane where coordinate `axis` is `value`.
fn plane(axis: i64, value: f64) -> bool {
    region.clear();
    return also_plane(axis, value);
}

// Adds a plane to the region of the last `plane` or `hole` call, so the
// region is their union. It gives false when the plane holds no node.
fn also_plane(axis: i64, value: f64) -> bool {
    guard axis >= 0 && axis < 3 && region.len < 8 else { abort("also_plane: 0 <= axis < 3 and at most 8 parts"); }
    region.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1 });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

// Selects the surface of the hole of radius r about the z axis through
// (cx, cy).
fn hole(cx: f64, cy: f64, r: f64) -> bool {
    guard r > 0.0 else { abort("hole: r > 0"); }
    region.clear();
    region.push(Part { hole: true, axis: 0, value: 0.0, c: [cx, cy], r: r, shape: -1 });
    return region_node_count() > 0;
}

// Selects the surface of the shape g of src/mesher/: the nodes where φ of
// g is 0.
fn surface(g: i64) -> bool {
    region.clear();
    return also_surface(g);
}

// Adds the surface of the shape g to the region of the last call that
// selects a region. It gives false when the surface holds no node.
fn also_surface(g: i64) -> bool {
    guard g >= 0 && region.len < 8 else { abort("also_surface: a shape and at most 8 parts"); }
    region.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: g });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

// Holds component `comp` (X, Y, Z, or ALL) of u, or T (comp 0), at `v` on
// the region.
fn fix(comp: i64, v: f64) { fix_region(comp, v); }

// Holds the unknown on the region at the field of the last `read_state`,
// interpolated at each node with the shape functions of the other mesh. It
// gives the number of nodes of the region outside the other mesh.
fn fix_state() -> i64 { return fix_region_source(); }

// Adds a uniform traction (force per area) on the region.
fn traction(tx: f64, ty: f64, tz: f64) { add_region_load([tx, ty, tz]); }

// Adds a uniform outward heat flux q·n = g (power per area) on the region.
fn flux(g: f64) { add_region_load([-g, 0.0, 0.0]); }

// Removes the fixed values, the loads and the solution, for the next
// step. The mesh and the physics stay.
fn clear_bcs() { clear_conditions(); }

// Assembles and solves by CG with an AMG preconditioner. Returns false
// if CG does not converge.
fn solve() -> bool { return solve_system(); }

// Component `comp` of u, or T (comp 0), at a point, and whether a cell
// holds the point.
fn value_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_VALUE, comp, x, y, z);
}

// Component (i, j) of the stress at a point.
fn stress_at(i: i64, j: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_STRESS, 3 * i + j, x, y, z);
}

// Component `comp` of the heat flux at a point.
fn flux_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_FLUX, comp, x, y, z);
}

// The sum of the reactions on the fixed dofs of the region.
fn reaction() -> f64, f64, f64 {
    let r = region_reaction();
    return r[0], r[1], r[2];
}

// The check `submodel_cut`: on each part of the region, the resultant force
// of the stress of this model and of the field of the last `read_state`.
// It gives the largest difference of the two, and the largest resultant of
// the other model.
fn cut_difference() -> f64, f64 {
    var most = 0.0;
    var scale = 0.0;
    for k in region.len {
        let own, other = part_resultants(k);
        var d = 0.0;
        var s = 0.0;
        for t in 3 {
            d += (own[t] - other[t]) * (own[t] - other[t]);
            s += other[t] * other[t];
        }
        most = max(most, sqrt(d));
        scale = max(scale, sqrt(s));
    }
    return most, scale;
}

// The largest von Mises stress at the nodes of the region.
fn max_von_mises() -> f64 { return region_max_von_mises(); }

// The sum of the reactions on all fixed dofs, for the check `equilibrium`.
fn total_reaction() -> f64, f64, f64 {
    let r = all_reaction();
    return r[0], r[1], r[2];
}

// The sum of all applied loads and sources.
fn applied_load() -> f64, f64, f64 {
    let f = total_load();
    return f[0], f[1], f[2];
}

// The scale of the check `equilibrium`: the larger of the sum of the sizes
// of the nodal loads, sources and reactions, and of the terms of K u on the
// fixed rows.
fn load_scale() -> f64 { return load_magnitude(); }

// Writes the mesh and the solution, with the name `field`, to a VTU file.
fn write_vtu(path: const u8[:], field: const u8[:]) -> bool { return write_vtu_file(path, field); }

// Facts about the mesh and the last solve, for the report.
fn dof_count() -> i64 { return dof_total(); }
fn element_count() -> i64 { return cell_count(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
