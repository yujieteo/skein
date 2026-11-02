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
// or `also_surface` call, less the nodes of each `except_plane` and
// `except_surface` call after it. A load acts on each facet of the region.
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
import .stages;
import .eigen;
import .wave;
import .radiation;
import .explicit;
import .edge;
import .mixed;
import .material;
import .solid;
import .contact;

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

// Meshes the box from (x0, y0, z0) to (x1, y1, z1) with nx × ny × nz
// hexahedra of order 1.
fn mesh_hex_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64, nx: i64, ny: i64, nz: i64) {
    build_hex_box([x0, y0, z0], [x1, y1, z1], nx, ny, nz);
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

// Meshes each shape of `shapes` alone, as mesh_shape, for the bodies of a
// contact: two bodies share no node. g is their union. False if the
// mesher does not finish.
fn mesh_bodies(d: i64, ord: i64, shapes: i64[:], g: i64, size: f64) -> bool {
    physics = -1;
    return build_bodies(d, ord, shapes, g, size);
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

// Gives the cells in the shape g of src/mesher/ the elastic material of E
// and nu, with the body force per volume (bx, by, bz), after `elastic` and
// `body_force` set the material of the other cells. False if g holds no
// cell.
fn material_in(g: i64, E: f64, nu: f64, bx: f64, by: f64, bz: f64) -> bool {
    return add_law(g, E, nu, [bx, by, bz]) > 0;
}

// Removes the cells in the shape g from the live set, for the next solves
// (spec §14.3). False if g holds no cell.
fn kill(g: i64) -> bool { return kill_cells(g) > 0; }

// Adds the dead cells in the shape g to the live set, strain-free on the
// solution of the last step. Call it after `clear_bcs`. False if g holds
// no cell.
fn birth(g: i64) -> bool { return birth_cells(g) > 0; }

// Steady heat conduction with conductivity k and a uniform source q per
// volume. In 2D, t is the thickness.
fn heat(k: f64, q: f64, t: f64) { set_heat(k, q, t); }

// The complex scalar wave div(grad(u)) + c u = s of a harmonic step
// (wave.goose), with u as components 0 (the real part) and 1 (the
// imaginary part). In 2D, t is the thickness.
fn wave(t: f64) { set_wave(t); }

// The terms of the next solve: c = c_re + i c_im, s = s_re + i s_im, and
// the damping ratio ζ that multiplies the stiffness by 1 + 2iζ.
fn wave_terms(c_re: f64, c_im: f64, s_re: f64, s_im: f64, damping: f64) {
    set_wave_terms([c_re, c_im], [s_re, s_im], damping);
}

// The PML of the cells of shape g (spec §13): past the box from
// (x0, y0, z0) to (x1, y1, z1), with the thickness h. Call it after
// wave_terms, which gives the wave number.
fn pml(g: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64, h: f64) {
    set_pml(g, [x0, y0, z0], [x1, y1, z1], h);
}

// The low region of a bloch group: the region of the last call that
// selects one. False if it holds no node.
fn bloch_low() -> bool { return keep_low(); }

// Ties each node of the region of the last call that selects one to the
// node of the low region at its point less (px, py, pz), with the phase
// theta (spec §12.6). False if a node has no partner.
fn bloch(px: f64, py: f64, pz: f64, theta: f64) -> bool { return tie_high([px, py, pz], theta); }

// The edge physics of `family: nedelec` (edge.goose): curl(ν curl(u)) +
// c u = s for a vector u on the edges of a P1 mesh, with ν = nu and
// c = c_re + i c_im on each cell. In 2D, t is the thickness.
fn edge(nu: f64, c_re: f64, c_im: f64, t: f64) { set_edge(nu, [c_re, c_im], t); }

// Gives the cells in the shape g of src/mesher/ the coefficients ν and c
// of the edge physics. False if g holds no cell.
fn edge_in(g: i64, nu: f64, c_re: f64, c_im: f64) -> bool { return edge_material(g, nu, [c_re, c_im]) > 0; }

// Adds the source s of the edge physics, with s = S(p) at each point p
// (f64[3]): the real parts of its 3 components, then the imaginary parts
// (f64[6]).
fn current<S>() { edge_source() { p => S(p) }; }

// The essential condition cross(n, u) = cross(n, g) of the edge physics
// on the region, with g = G(p) as in `current`. False if the region holds
// no edge.
fn tangent<G>() -> bool { return fix_tangent() { p => G(p) } > 0; }

// The natural condition cross(n, h) = z u of the edge physics on the
// region, with h = ν curl(u) and z = z_re + i z_im. False if the region
// holds no facet.
fn impedance(z_re: f64, z_im: f64) -> bool { return edge_impedance([z_re, z_im]) > 0; }

// The tree gauge of the edge physics (`gauge: tree`): holds at 0 the
// edges of a spanning tree of the cells with c = 0. Call it after the
// essential conditions. It gives the edges of the tree.
fn gauge() -> i64 { return tree_gauge(); }

// Part `part` (0 the real part, 1 the imaginary part) of component `comp`
// of u, or of curl(u) with `curl`, at a point of the edge physics, and
// whether a cell holds the point.
fn edge_at(curl: bool, part: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard (part == 0 || part == 1) && comp >= 0 && comp < 3 else { abort("edge_at: part 0 or 1 and 0 <= comp < 3"); }
    let v, found = edge_value_at([x, y, z], curl);
    return v[3 * part + comp], found;
}

// The integral of component `comp` of u, or of curl(u) with `curl`, along
// the line from (x0, y0, z0) to (x1, y1, z1), by the midpoint rule on n
// pieces: its real and imaginary parts, and whether each point is near a
// cell.
fn edge_line_integral(curl: bool, comp: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64, n: i64) -> f64, f64, bool {
    let s, found = edge_line([x0, y0, z0], [x1, y1, z1], comp, curl, n);
    return s[0], s[1], found;
}

// The radiation intensity (power per solid angle) in the direction
// (dx, dy, dz) of the field on the closed surface of the region, by the
// near-to-far-field transform into a medium of wave number k and
// impedance eta, with E = u and H = (hf_re + i hf_im) curl(u)
// (edge.goose, gap G13). 3D only.
fn far_field(dx: f64, dy: f64, dz: f64, k: f64, eta: f64, hf_re: f64, hf_im: f64) -> f64 {
    surface_samples([hf_re, hf_im]);
    return far_intensity([dx, dy, dz], k, eta);
}

// The power out of the closed surface of the region, ∮ ½ Re(E × H*) · n dA,
// with E and H as in far_field.
fn radiated_power(hf_re: f64, hf_im: f64) -> f64 {
    surface_samples([hf_re, hf_im]);
    return surface_power();
}

// The mixed Darcy flow of `family: raviart_thomas` (mixed.goose):
// q + k grad(p) = 0 and div(q) = f, with q on the facets and p on the
// cells of a P1 mesh. In 2D, t is the thickness.
fn darcy(k: f64, f: f64, t: f64) { set_mixed(k, f, t); }

// Gives the cells in the shape g of src/mesher/ the coefficients k and f
// of the Darcy flow. False if g holds no cell.
fn darcy_in(g: i64, k: f64, f: f64) -> bool { return mixed_material(g, k, f) > 0; }

// Sets f of each cell to its mean of F(p) (p: f64[3], F gives an f64).
fn darcy_source<F>() { mixed_source() { p => F(p) }; }

// The pressure p = v on the boundary facets of the region: a natural
// condition of the mixed form. False if the region holds no facet.
fn pressure(v: f64) -> bool { return mixed_pressure(v) > 0; }

// The flux dot(q, n) = g on the boundary facets of the region: an
// essential condition. False if the region holds no facet.
fn normal_flux(g: f64) -> bool { return mixed_normal_flux(g) > 0; }

// Component `comp` of q (0, 1 or 2), or p (comp 3), at a point of the
// Darcy flow, and whether a cell holds the point.
fn darcy_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard comp >= 0 && comp < 4 else { abort("darcy_at: 0 <= comp < 4"); }
    let v, found = mixed_value_at([x, y, z]);
    return v[comp], found;
}

// The flux out of the boundary facets of the region, ∮ dot(q, n) dA.
fn boundary_flux() -> f64 { return mixed_flux(); }

// Writes the mesh with p and q of the Darcy flow, named `p_name` and
// `q_name`, to a VTU file.
fn write_darcy(path: const u8[:], p_name: const u8[:], q_name: const u8[:]) -> bool { return write_mixed_vtu(path, p_name, q_name); }

// Selects the boundary plane where coordinate `axis` is `value`.
fn plane(axis: i64, value: f64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_plane(axis, value);
}

// Adds a plane to the region of the last `plane` or `hole` call, so the
// region is their union. It gives false when the plane holds no node.
fn also_plane(axis: i64, value: f64) -> bool {
    guard axis >= 0 && axis < 3 && region.len < 16 else { abort("also_plane: 0 <= axis < 3 and at most 16 parts"); }
    region.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1 });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

// Selects the surface of the hole of radius r about the z axis through
// (cx, cy).
fn hole(cx: f64, cy: f64, r: f64) -> bool {
    guard r > 0.0 else { abort("hole: r > 0"); }
    region.clear();
    region_out.clear();
    region_body = -1;
    region.push(Part { hole: true, axis: 0, value: 0.0, c: [cx, cy], r: r, shape: -1, node: -1 });
    return region_node_count() > 0;
}

// Selects the surface of the shape g of src/mesher/: the nodes where φ of
// g is 0.
fn surface(g: i64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_surface(g);
}

// Adds the surface of the shape g to the region of the last call that
// selects a region. It gives false when the surface holds no node.
fn also_surface(g: i64) -> bool {
    guard g >= 0 && region.len < 16 else { abort("also_surface: a shape and at most 16 parts"); }
    region.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: g, node: -1 });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

// Selects the node nearest to the point (x, y, z).
fn point(x: f64, y: f64, z: f64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_point(x, y, z);
}

// Adds the node nearest to the point (x, y, z) to the region of the last
// call that selects a region.
fn also_point(x: f64, y: f64, z: f64) -> bool {
    guard region.len < 16 && node_count() > 0 else { abort("also_point: a mesh and at most 16 parts"); }
    var best = 0;
    var least = 1.0e300;
    for n in node_count() {
        let d = (node_x(n, 0) - x) * (node_x(n, 0) - x) + (node_x(n, 1) - y) * (node_x(n, 1) - y) + (node_x(n, 2) - z) * (node_x(n, 2) - z);
        if d < least {
            least = d;
            best = n;
        }
    }
    region.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: -1, node: best });
    return true;
}

// Keeps only the nodes of body k of mesh_bodies in the region of the last
// call that selects a region. It gives false when the region holds no node.
fn in_body(k: i64) -> bool {
    guard k >= 0 && k < body_shapes.len else { abort("in_body: a body of mesh_bodies"); }
    region_body = k;
    return region_node_count() > 0;
}

// The force (fx, fy, fz) on each node of the region.
fn node_force(fx: f64, fy: f64, fz: f64) { add_node_force([fx, fy, fz]); }

// The mean of the stress component (i, j) over the body.
fn mean_stress(i: i64, j: i64) -> f64 { return body_mean(FIELD_STRESS, 3 * i + j); }
// Removes the nodes of the plane from the region: for example the node
// of a face on a symmetry plane, whose reaction also holds the force of
// that plane. It gives false when the region holds no node.
fn except_plane(axis: i64, value: f64) -> bool {
    guard axis >= 0 && axis < 3 && region_out.len < 8 else { abort("except_plane: 0 <= axis < 3 and at most 8 parts"); }
    region_out.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1 });
    return region_node_count() > 0;
}

// Removes the nodes of the surface of shape g from the region.
fn except_surface(g: i64) -> bool {
    guard g >= 0 && region_out.len < 8 else { abort("except_surface: a shape and at most 8 parts"); }
    region_out.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: g, node: -1 });
    return region_node_count() > 0;
}

// Holds component `comp` (X, Y, Z, or ALL) of u, or T (comp 0), at `v` on
// the region.
fn fix(comp: i64, v: f64) { fix_region(comp, v); }

// Holds component `comp` of each node of the region at VALUE([x, y, z]), a
// function of the program: a fixed value that varies on the region.
fn fix_with<VALUE>(comp: i64) { fix_region_with(comp, VALUE); }

// Holds the unknown on the region at the field of the last `read_state`,
// interpolated at each node with the shape functions of the other mesh. It
// gives the number of nodes of the region outside the other mesh.
fn fix_state() -> i64 { return fix_region_source(); }

// Adds a uniform traction (force per area) on the region.
fn traction(tx: f64, ty: f64, tz: f64) { add_region_load([tx, ty, tz]); }

// Adds the pressure p on the region: the traction -p n on its outward
// normal n.
fn pressure_load(p: f64) { add_region_pressure(p); }

// Adds a uniform outward heat flux q·n = g (power per area) on the region.
fn flux(g: f64) { add_region_load([-g, 0.0, 0.0]); }

// The condition dot(q, n) = h (T - t_inf) on the region.
fn film(h: f64, t_inf: f64) { add_film(h, t_inf); }

// The radiation to the surroundings at t_inf on the region: dot(q, n) = c
// (T^4 - t_inf^4) with c = eps sigma in the units of the model.
fn radiate_ambient(c: f64, t_inf: f64) { add_ambient(c, t_inf); }

// The conductivity k(T) of a heat physics at the temperatures lo, lo + step,
// .. hi, linear in between and beyond the ends. A step with it iterates.
fn conductivity_curve(lo: f64, hi: f64, v: f64[:]) { set_conductivity_curve(lo, hi, v); }

// The temperature of the free nodes before the first step, which is the
// start of the iteration of a nonlinear step.
fn initial_temperature(t: f64) { initial_value = t; }

// The tolerance of the nonlinear steps: the change of T relative to its
// largest value.
fn nonlinear_tolerance(rtol: f64) {
    guard rtol > 0.0 else { abort("nonlinear_tolerance: rtol > 0"); }
    nl_tolerance = rtol;
}

// Adds the facets of the region to the enclosure of the step, with the
// emissivity eps: the natural condition dot(q, n) = radiation(eps).
// `sb` is the Stefan-Boltzmann constant in the units of the model.
fn radiate(eps: f64, sb: f64) { add_radiation(eps, sb); }

// Removes the fixed values, the loads and the solution, for the next
// step. The mesh, the physics and the live set stay. The solution stays as
// the value of the nodes that no live cell holds.
fn clear_bcs() {
    clear_conditions();
    clear_rigid();
    clear_films();
    clear_ambient();
    clear_radiation();
    clear_impedance();
    if physics == MIXED { clear_pressure(); }
}

// Assembles and solves by CG with an AMG preconditioner, or the wave
// physics by GMRES with ILU(0). Returns false if the solver does not
// converge.
fn solve() -> bool {
    if physics == EDGE { return edge_solve(); }
    if physics == MIXED { return mixed_solve(); }
    return if physics == WAVE { wave_solve() } else { solve_system() };
}

// --- contact and the rigid face (contact.goose, spec §12.4, §12.1) -------------

// The facets of the region are the master of the next contact_slave. It
// gives false when the region holds no facet.
fn contact_master() -> bool { return set_master() > 0; }

// Pairs each quadrature point of the facets of the region, the slave, with
// the master of the last contact_master: penalty contact with the normal
// stiffness kn, the stick stiffness kt and the friction coefficient mu (kn
// and kt are a pressure per length). It gives the points paired.
fn contact_slave(kn: f64, kt: f64, mu: f64) -> i64 { return add_slave(kn, kt, mu); }

// The region is the rigid face of the step (`u = rigid(p)`): its nodes
// move as a rigid body. It gives false when the region holds no node.
fn rigid_face() -> bool { return set_rigid() > 0; }

// The moment (mx, my, mz) on the rigid face of the step: about z in 2D.
fn rigid_moment(mx: f64, my: f64, mz: f64) { add_rigid_moment([mx, my, mz]); }

// A static step of the linear law with the contact and the rigid face in
// `ramp` increments of Newton's method (solid.goose). False if an
// increment does not converge.
fn solve_contact(ramp: i64) -> bool { return solve_solid(ramp, no_energy); }

// The largest contact pressure of the last solve.
fn contact_pressure() -> f64 { return max_contact_pressure(); }

// --- the nonlinear solid (solid.goose, spec §10.6) ----------------------------

// J2 plasticity of the von Mises yield sy + hk α, with the E and nu of the
// last `elastic` call: at small strain, or at finite strain with the
// multiplicative split (`finite`).
fn plastic(sy: f64, hk: f64, finite: bool) { set_solid(if finite { SOLID_J2_FINITE } else { SOLID_J2 }, sy, hk); }

// An isotropic hyperelastic energy of the program (solve_energy), in 3d,
// or in 2d in plane strain with the thickness t. It needs no E and nu.
fn hyperelastic(t: f64) {
    set_elastic(1.0, 0.0, if dim == 2 { PLANE_STRAIN } else { SOLID_3D }, t);
    set_solid(SOLID_ENERGY, 0.0, 0.0);
}

// A static step of J2 at small strain in `ramp` increments. False if an
// increment does not converge.
fn solve_plastic(ramp: i64) -> bool { return solve_solid(ramp, no_energy); }

// A static step of the energy ENERGY in `ramp` increments. ENERGY takes
// the principal values of C as hyper-dual numbers (hd_seed) and gives W.
fn solve_energy<ENERGY>(ramp: i64) -> bool { return solve_solid(ramp, ENERGY); }

// The explicit step of J2 at finite strain to `end` with the step dt (dt
// <= 0 for the stable step), from the velocity (vx, vy, vz) of each free
// node. False if the state is not finite.
fn explicit_solid(end: f64, dt: f64, vx: f64, vy: f64, vz: f64) -> bool { return solid_explicit(end, dt, [vx, vy, vz]); }

// The Newton iterations of the last static step of the nonlinear solid.
fn newton_iterations() -> i64 { return newton_steps; }

// The transient step central_difference (explicit.goose) to the time
// `end`, with the step dt, or the stable step for dt <= 0, from rest.
// Returns false if the state is not finite.
fn central_difference(end: f64, dt: f64) -> bool { return central_difference_cpu(end, dt); }

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

// The mean of component comp of the solution on the region, weighted by
// area. NaN if the region holds no facet.
fn mean_value(comp: i64) -> f64 { return region_mean(comp); }

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

// The largest sign σ_ij at the nodes of the region.
fn max_stress(i: i64, j: i64, sign: f64) -> f64 { return region_max_stress(i, j, sign); }

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

// The density of the mass of a modal step.
fn set_density(rho: f64) {
    guard rho > 0.0 else { abort("set_density: rho > 0"); }
    density = rho;
}

// The `count` lowest modes of K x = ω² M x on the free dofs (spec §14.2).
// False if LOBPCG does not converge.
fn modal(count: i64) -> bool { return eigen_solve(MASS, count); }

// The `count` lowest load factors λ of (K + λ K_G) x = 0, with K_G of the
// stress of the solution of the last step. False if LOBPCG does not
// converge.
fn buckling(count: i64) -> bool { return eigen_solve(GEOMETRIC, count); }

// eigenvalue(k) of the last modal or buckling step: a frequency in 1/s of
// the time of the units, or a load factor.
fn eigenvalue(k: i64) -> f64 { return mode_value(k); }

// Writes mode k of the last modal or buckling step as the unknown.
fn write_mode(path: const u8[:], field: const u8[:], k: i64) -> bool {
    guard k >= 1 && modes.len >= k * sol.len else { abort("write_mode: 1 <= k <= the count of the last eigen step"); }
    var keep: f64[>..] = [];
    keep.append(sol);
    let n = sol.len;
    for d in n { sol[d] = modes[n * (k - 1) + d]; }
    let done = write_vtu_file(path, field);
    for d in n { sol[d] = keep[d]; }
    return done;
}

// Facts about the mesh and the last solve, for the report.
fn dof_count() -> i64 { return dof_total(); }
fn element_count() -> i64 { return cell_count(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
