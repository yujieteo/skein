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
import .curve;
import .harmonic;
import .radiation;
import .explicit;
import .edge;
import .edge_modal;
import .mixed;
import .material;
import .solid;
import .contact;
import .heat_transient;
import .field;

let X = 0;
let Y = 1;
let Z = 2;
let ALL = -1;

fn mesh_box(d: i64, ord: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
            nx: i64, ny: i64, nz: i64) {
    build_box(d, ord, [x0, y0, z0], [x1, y1, z1], nx, ny, nz);
    physics = -1;
}

fn set_cube_gauss(n: i64) {
    guard n == 2 || n == 3 else { abort("set_cube_gauss: 2 or 3 points"); }
    cube_gauss = n;
}

fn promote_hexahedra() {
    promote_hex();
    physics = -1;
}

fn mesh_hex_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64, nx: i64, ny: i64, nz: i64) {
    build_hex_box([x0, y0, z0], [x1, y1, z1], nx, ny, nz);
    physics = -1;
}

fn mesh_hole(d: i64, ord: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             r: f64, size: f64, fine: f64) {
    build_hole(d, ord, [x0, y0, z0], [x1, y1, z1], r, size, fine);
    physics = -1;
}

fn mesh_shape(d: i64, ord: i64, g: i64, size: f64) -> bool {
    physics = -1;
    return build_shape(d, ord, g, size);
}

fn mesh_hex_shape(g: i64, size: f64, planes: i64[:]) -> bool {
    physics = -1;
    return build_hex_shape(g, size, planes);
}

fn mesh_bodies(d: i64, ord: i64, shapes: i64[:], g: i64, size: f64) -> bool {
    physics = -1;
    return build_bodies(d, ord, shapes, g, size);
}

fn mesh_simplices(d: i64, ord: i64, xyz: f64[:], vertices: i64[:]) {
    build_simplices(d, ord, xyz, vertices);
    physics = -1;
}

fn write_mesh(dir: const u8[:]) -> bool { return save_mesh(dir); }

fn read_mesh(dir: const u8[:]) -> bool {
    physics = -1;
    return load_mesh(dir);
}

fn write_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return save_array(dir, str(step, ".fixed"), fixed) && save_array(dir, str(step, ".prescribed"), prescribed) &&
           save_array(dir, str(step, ".load"), load) && save_array(dir, str(step, ".sol"), sol) &&
           save_array(dir, str(step, ".row_start"), row_start) && save_array(dir, str(step, ".col"), col) &&
           save_array(dir, str(step, ".val"), val) && save_array(dir, str(step, ".source"), assembled_source);
}

fn read_solution(dir: const u8[:], step: const u8[:]) -> bool {
    return load_array(dir, str(step, ".fixed"), fixed) && load_array(dir, str(step, ".prescribed"), prescribed) &&
           load_array(dir, str(step, ".load"), load) && load_array(dir, str(step, ".sol"), sol) &&
           load_array(dir, str(step, ".row_start"), row_start) && load_array(dir, str(step, ".col"), col) &&
           load_array(dir, str(step, ".val"), val) && load_array(dir, str(step, ".source"), assembled_source);
}

fn write_state(dir: const u8[:], step: const u8[:]) -> bool { return save_state(dir, step); }

fn read_state(dir: const u8[:], step: const u8[:], length_scale: f64, value_scale: f64, modulus_scale: f64) -> bool {
    return load_source(dir, step, length_scale, value_scale, modulus_scale);
}

fn cells_for(length: f64, size: f64) -> i64 {
    guard length > 0.0 && size > 0.0 else { abort("cells_for: length > 0 and size > 0"); }
    return max(ceil(length / size - 1.0e-9) as i64, 1);
}

fn elastic(E: f64, nu: f64) { set_elastic(E, nu, SOLID_3D, 1.0); }
fn anisotropic(c: f64[36], s0: f64[6]) { set_aniso(c, s0); }

fn anisotropic_in(g: i64, c: f64[36], s0: f64[6], bx: f64, by: f64, bz: f64) -> bool { return add_aniso_law(g, c, s0, [bx, by, bz]) > 0; }

fn elastic_plane_stress(E: f64, nu: f64, t: f64) { set_elastic(E, nu, PLANE_STRESS, t); }
fn elastic_plane_strain(E: f64, nu: f64, t: f64) { set_elastic(E, nu, PLANE_STRAIN, t); }

fn incompatible(on: bool) { incompatible_modes = on; }

fn body_force(bx: f64, by: f64, bz: f64) {
    guard physics == ELASTIC else { abort("body_force: after elastic"); }
    volume_force = [bx, by, bz];
}

fn material_in(g: i64, E: f64, nu: f64, bx: f64, by: f64, bz: f64) -> bool {
    return add_law(g, E, nu, [bx, by, bz]) > 0;
}

fn heat_in(g: i64, k: f64, q: f64) -> bool { return add_heat_law(g, k, q) > 0; }

fn kill(g: i64) -> bool { return kill_cells(g) > 0; }

fn birth(g: i64) -> bool { return birth_cells(g) > 0; }

fn keep(slot: i64) { keep_field(slot); }

fn field_term(kind: i64, code: i64[:], consts: f64[:]) { add_field_term(kind, code, consts); }

fn field_term_in(law: i64, kind: i64, code: i64[:], consts: f64[:]) { add_field_term_in(law, kind, code, consts); }

fn heat(k: f64, q: f64, t: f64) { set_heat(k, q, t); }

fn wave(t: f64) { set_wave(t); }

fn wave_terms(c_re: f64, c_im: f64, s_re: f64, s_im: f64, damping: f64) {
    set_wave_terms([c_re, c_im], [s_re, s_im], damping);
}

fn wave_source_in<S>(g: i64) { wave_source(g) { p => S(p) }; }

fn wave_nu(nu: f64) { wave_material_all(nu); }

fn wave_region_integral(g: i64) -> f64[3] { return wave_cells_integral(g); }

fn wave_in(g: i64, nu: f64, c_re: f64, c_im: f64, s_re: f64, s_im: f64) -> bool {
    return wave_material(g, nu, [c_re, c_im], [s_re, s_im]) > 0;
}

fn pml(g: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64, h: f64) {
    set_pml(g, [x0, y0, z0], [x1, y1, z1], h);
}

fn bloch_low() -> bool { return keep_low(); }

fn bloch(px: f64, py: f64, pz: f64, theta: f64) -> bool { return tie_high([px, py, pz], theta); }

fn edge(nu: f64, c_re: f64, c_im: f64, t: f64) { set_edge(nu, [c_re, c_im], t); }

fn edge_order(p: i64) { set_edge_order(p); }

fn edge_mass_in(g: i64, m_re: f64, m_im: f64) -> bool { return edge_mass(g, [m_re, m_im]) > 0; }

fn impedance_group(group: i64) -> bool { return edge_impedance_in([0.0, 0.0], group) > 0; }

fn resonance<Z>(count: i64, shift: f64) -> bool { return edge_modal(count, shift) { g, w => Z(g, w) }; }

fn modal_frequency(k: i64) -> f64, f64 {
    let f = edge_frequency(k);
    return f[0], f[1];
}

fn write_modal_mode(path: const u8[:], field: const u8[:], k: i64) -> bool { return edge_write_mode(path, field, k); }

fn edge_in(g: i64, nu: f64, c_re: f64, c_im: f64) -> bool { return edge_material(g, nu, [c_re, c_im]) > 0; }

fn current<S>() { edge_source(-2) { p => S(p) }; }

fn current_in<S>(g: i64) { edge_source(g) { p => S(p) }; }

fn tangent<G>() -> bool { return fix_tangent() { p => G(p) } > 0; }

fn impedance(z_re: f64, z_im: f64) -> bool { return edge_impedance([z_re, z_im]) > 0; }

fn traction<G>() -> bool { return edge_traction() { p => G(p) } > 0; }

fn facet_integral<G>(field: bool) -> f64[2] { return edge_facet_integral(field) { p => G(p) }; }

fn gauge() -> i64 { return tree_gauge(); }

fn edge_at(curl: bool, part: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard (part == 0 || part == 1) && comp >= 0 && comp < 3 else { abort("edge_at: part 0 or 1 and 0 <= comp < 3"); }
    let v, found = edge_value_at([x, y, z], curl);
    return v[3 * part + comp], found;
}


fn far_field(dx: f64, dy: f64, dz: f64, k: f64, eta: f64, hf_re: f64, hf_im: f64) -> f64 {
    surface_samples([hf_re, hf_im]);
    return far_intensity([dx, dy, dz], k, eta);
}

fn darcy(k: f64, f: f64, t: f64) { set_mixed(k, f, t); }

fn darcy_source<F>() { mixed_source() { p => F(p) }; }

fn pressure(v: f64) -> bool { return mixed_pressure(v) > 0; }

fn normal_flux(g: f64) -> bool { return mixed_normal_flux(g) > 0; }

fn darcy_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard comp >= 0 && comp < 4 else { abort("darcy_at: 0 <= comp < 4"); }
    let v, found = mixed_value_at([x, y, z]);
    return v[comp], found;
}

fn boundary_flux() -> f64 { return mixed_flux(); }

fn write_darcy(path: const u8[:], p_name: const u8[:], q_name: const u8[:]) -> bool { return write_mixed_vtu(path, p_name, q_name); }

fn plane(axis: i64, value: f64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_plane(axis, value);
}

fn also_plane_cut(axis: i64, value: f64, la: i64, lv: f64, lr: i64) -> bool {
    guard axis >= 0 && axis < 3 && la >= 0 && la < 3 && region.len < 16 && (lr == 2 || lr == 1 || lr == -2 || lr == -1) else {
        abort("also_plane_cut: 0 <= axis < 3, a limit axis, a relation and at most 16 parts");
    }
    region.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1, rel: 0, la: la, lv: lv, lr: lr });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

fn plane_cut(axis: i64, value: f64, la: i64, lv: f64, lr: i64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_plane_cut(axis, value, la, lv, lr);
}

fn also_plane(axis: i64, value: f64) -> bool {
    guard axis >= 0 && axis < 3 && region.len < 16 else { abort("also_plane: 0 <= axis < 3 and at most 16 parts"); }
    region.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1, rel: 0, la: -1, lv: 0.0, lr: 0 });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

fn hole(cx: f64, cy: f64, r: f64) -> bool {
    guard r > 0.0 else { abort("hole: r > 0"); }
    region.clear();
    region_out.clear();
    region_body = -1;
    region.push(Part { hole: true, axis: 0, value: 0.0, c: [cx, cy], r: r, shape: -1, node: -1, rel: 0, la: -1, lv: 0.0, lr: 0 });
    return region_node_count() > 0;
}

fn surface(g: i64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_surface(g);
}

fn also_surface(g: i64) -> bool {
    guard g >= 0 && region.len < 16 else { abort("also_surface: a shape and at most 16 parts"); }
    region.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: g, node: -1, rel: 0, la: -1, lv: 0.0, lr: 0 });
    let k = region.len - 1;
    for n in node_count() { if on_part(n, k) { return true; } }
    return false;
}

fn point(x: f64, y: f64, z: f64) -> bool {
    region.clear();
    region_out.clear();
    region_body = -1;
    return also_point(x, y, z);
}

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
    region.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: -1, node: best, rel: 0, la: -1, lv: 0.0, lr: 0 });
    return true;
}

fn in_body(k: i64) -> bool {
    guard k >= 0 && k < body_shapes.len else { abort("in_body: a body of mesh_bodies"); }
    region_body = k;
    return region_node_count() > 0;
}

fn node_force(fx: f64, fy: f64, fz: f64) { add_node_force([fx, fy, fz]); }

fn mean_stress(i: i64, j: i64) -> f64 { return body_mean(FIELD_STRESS, 3 * i + j); }
fn mean_stress_box(i: i64, j: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> f64 {
    return box_mean(FIELD_STRESS, 3 * i + j, [x0, y0, z0], [x1, y1, z1]);
}
fn except_plane(axis: i64, value: f64) -> bool {
    guard axis >= 0 && axis < 3 && region_out.len < 8 else { abort("except_plane: 0 <= axis < 3 and at most 8 parts"); }
    region_out.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1, rel: 0, la: -1, lv: 0.0, lr: 0 });
    return region_node_count() > 0;
}

fn except_coord(axis: i64, value: f64, rel: i64) -> bool {
    guard axis >= 0 && axis < 3 && region_out.len < 8 && (rel == 2 || rel == 1 || rel == -2 || rel == -1) else { abort("except_coord: 0 <= axis < 3, a relation and at most 8 parts"); }
    region_out.push(Part { hole: false, axis: axis, value: value, c: [0.0, 0.0], r: 0.0, shape: -1, node: -1, rel: rel, la: -1, lv: 0.0, lr: 0 });
    return region_node_count() > 0;
}

fn except_radius(cx: f64, cy: f64, r: f64, rel: i64) -> bool {
    guard region_out.len < 8 && r > 0.0 && (rel == 2 || rel == 1 || rel == -2 || rel == -1) else { abort("except_radius: r > 0, a relation and at most 8 parts"); }
    region_out.push(Part { hole: true, axis: 0, value: 0.0, c: [cx, cy], r: r, shape: -1, node: -1, rel: rel, la: -1, lv: 0.0, lr: 0 });
    return region_node_count() > 0;
}

fn except_surface(g: i64) -> bool {
    guard g >= 0 && region_out.len < 8 else { abort("except_surface: a shape and at most 8 parts"); }
    region_out.push(Part { hole: false, axis: 0, value: 0.0, c: [0.0, 0.0], r: 0.0, shape: g, node: -1, rel: 0, la: -1, lv: 0.0, lr: 0 });
    return region_node_count() > 0;
}

fn fix(comp: i64, v: f64) { fix_region(comp, v); }

fn fix_with<VALUE>(comp: i64) { fix_region_with(comp, VALUE); }

fn fix_state() -> i64 { return fix_region_source(); }

fn traction(tx: f64, ty: f64, tz: f64) { add_region_load([tx, ty, tz]); }

fn traction_with<G>() { add_region_traction_with(G); }

fn pressure_load(p: f64) { add_region_pressure(p); }

fn flux(g: f64) { add_region_load([-g, 0.0, 0.0]); }

fn film(h: f64, t_inf: f64) { add_film(h, t_inf); }

fn robin(z_re: f64, z_im: f64) { add_robin([z_re, z_im]); }

fn radiate_ambient(c: f64, t_inf: f64) { add_ambient(c, t_inf); }

fn conductivity_curve(lo: f64, hi: f64, v: f64[:]) { set_conductivity_curve(lo, hi, v); }

fn initial_temperature(t: f64) { initial_value = t; }

fn nonlinear_tolerance(rtol: f64) {
    guard rtol > 0.0 else { abort("nonlinear_tolerance: rtol > 0"); }
    nl_tolerance = rtol;
}

fn radiate(eps: f64, sb: f64) { add_radiation(eps, sb); }

fn clear_bcs() {
    clear_conditions();
    clear_rigid();
    clear_films();
    clear_ambient();
    clear_radiation();
    clear_impedance();
    clear_robin();
    if physics == MIXED { clear_pressure(); }
}

fn solve() -> bool {
    if physics == EDGE { return edge_solve(); }
    if physics == MIXED { return mixed_solve(); }
    if physics == ELASTIC && harmonic_omega > 0.0 { return harmonic_solve(); }
    return if physics == WAVE { wave_solve() } else { solve_system() };
}

fn contact_master() -> bool { return set_master() > 0; }

fn contact_slave(kn: f64, kt: f64, mu: f64) -> i64 { return add_slave(kn, kt, mu); }

fn rigid_face() -> bool { return set_rigid() > 0; }

fn rigid_moment(mx: f64, my: f64, mz: f64) { add_rigid_moment([mx, my, mz]); }

fn solve_contact(ramp: i64) -> bool { return solve_solid(ramp, no_energy); }

fn plastic(sy: f64, hk: f64, finite: bool) { set_solid(if finite { SOLID_J2_FINITE } else { SOLID_J2 }, sy, hk); }

fn hyperelastic(t: f64) {
    set_elastic(1.0, 0.0, if dim == 2 { PLANE_STRAIN } else { SOLID_3D }, t);
    set_solid(SOLID_ENERGY, 0.0, 0.0);
}

fn solve_plastic(ramp: i64) -> bool { return solve_solid(ramp, no_energy); }

fn solve_energy<ENERGY>(ramp: i64) -> bool { return solve_solid(ramp, ENERGY); }

fn explicit_solid(end: f64, dt: f64, vx: f64, vy: f64, vz: f64) -> bool { return solid_explicit(end, dt, [vx, vy, vz]); }

fn newton_iterations() -> i64 { return newton_steps; }

fn central_difference(end: f64, dt: f64) -> bool { return central_difference_cpu(end, dt); }

fn value_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_VALUE, comp, x, y, z);
}

fn stress_at(i: i64, j: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_STRESS, 3 * i + j, x, y, z);
}

fn flux_at(comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    return field_at(FIELD_FLUX, comp, x, y, z);
}

fn reaction() -> f64, f64, f64 {
    let r = region_reaction();
    return r[0], r[1], r[2];
}

fn mean_value(comp: i64) -> f64 { return region_mean(comp); }

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

fn integral_expr(code: i64[:], consts: f64[:]) -> f64 { return region_integral(code, consts); }

fn extreme_value(comp: i64, sign: f64) -> f64 { return region_extreme(comp, sign); }

fn max_von_mises() -> f64 { return region_max_von_mises(); }

fn max_stress(i: i64, j: i64, sign: f64) -> f64 { return region_max_stress(i, j, sign); }

fn total_reaction() -> f64, f64, f64 {
    let r = all_reaction();
    return r[0], r[1], r[2];
}

fn applied_load() -> f64, f64, f64 {
    let f = total_load();
    return f[0], f[1], f[2];
}

fn load_scale() -> f64 { return load_magnitude(); }

fn write_vtu(path: const u8[:], field: const u8[:]) -> bool { return write_vtu_file(path, field); }

fn set_density(rho: f64) {
    guard rho > 0.0 else { abort("set_density: rho > 0"); }
    density = rho;
}

fn set_capacity(c: f64) { heat_capacity(c); }
fn heat_begin(end: f64, dt: f64, theta: f64) -> i64 { return transient_heat_begin(end, dt, theta); }
fn heat_advance() -> bool { return transient_heat_advance(); }
fn heat_end() { transient_heat_end(); }

fn modal(count: i64) -> bool { return eigen_solve(MASS, count); }

fn buckling(count: i64) -> bool { return eigen_solve(GEOMETRIC, count); }

fn eigenvalue(k: i64) -> f64 { return mode_value(k); }

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

fn dof_count() -> i64 { return dof_total(); }
fn element_count() -> i64 { return cell_count(); }
fn cg_iterations() -> i64 { return last_iterations; }
fn solve_seconds() -> f64 { return last_seconds; }
