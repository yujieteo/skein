import std;
import math;
namespace shells;
import .laminate;
import .mesh;

let UX = 0;
let UY = 1;
let UZ = 2;
let RX = 3;
let RY = 4;
let RZ = 5;

struct Shell {
    n: i64[9],
    count: i64,
    section: i64,
    axis: f64[3],
}

struct Beam {
    n: i64[2],
    section: i64,
    y_axis: f64[3],
    offset: f64[6],
}

struct BeamSection {
    e: f64,
    g: f64,
    area: f64,
    iy: f64,
    iz: f64,
    j: f64,
    rho: f64,
}

var coords: f64[>..] = [];
var directors: f64[>..] = [];
var shells: Shell[>..] = [];
var beams: Beam[>..] = [];
var beam_sections: BeamSection[>..] = [];
var fixed: u8[>..] = [];
var load: f64[>..] = [];
var disp: f64[>..] = [];
var disp_is_mode = false;

fn node_count() -> i64 { return coords.len / 3; }
fn dof_count() -> i64 { return 6 * node_count(); }

fn clear() {
    coords.clear();
    directors.clear();
    shells.clear();
    beams.clear();
    beam_sections.clear();
    sections.clear();
    grids.clear();
    grid_lines.clear();
    grid_cells.clear();
    fixed.clear();
    load.clear();
    disp.clear();
    disp_is_mode = false;
}

fn add_node(x: f64, y: f64, z: f64) -> i64 {
    coords.push(x);
    coords.push(y);
    coords.push(z);
    for c in 3 { directors.push(0.0); }
    for d in 6 {
        fixed.push(0);
        load.push(0.0);
        disp.push(0.0);
    }
    return node_count() - 1;
}

fn node_point(n: i64) -> f64[3] {
    guard n >= 0 && n < node_count() else { abort("node_point: 0 <= n < node_count()"); }
    return [coords[3 * n], coords[3 * n + 1], coords[3 * n + 2]];
}

fn set_director(n: i64, vx: f64, vy: f64, vz: f64) {
    guard n >= 0 && n < node_count() && directors.len == coords.len else { abort("set_director: 0 <= n < node_count(), and a director for each node"); }
    let v = unit([vx, vy, vz]);
    for c in 3 { directors[3 * n + c] = v[c]; }
}

fn add_shell(n0: i64, n1: i64, n2: i64, n3: i64, section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    guard min(min(n0, n1), min(n2, n3)) >= 0 && max(max(n0, n1), max(n2, n3)) < node_count() && section >= 0 && section < sections.len else { abort("add_shell: each node < node_count(), and 0 <= section < sections.len"); }
    shells.push(Shell { n: [n0, n1, n2, n3, -1, -1, -1, -1, -1], count: 4, section: section, axis: [ax, ay, az] });
    return shells.len - 1;
}

fn add_shell9(n: i64[9], section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    var lo = n[0];
    var hi = n[0];
    for k in n { lo = min(lo, k); hi = max(hi, k); }
    guard lo >= 0 && hi < node_count() && section >= 0 && section < sections.len else { abort("add_shell9: each node < node_count(), and 0 <= section < sections.len"); }
    shells.push(Shell { n: n, count: 9, section: section, axis: [ax, ay, az] });
    return shells.len - 1;
}

fn beam_section(e: f64, nu: f64, area: f64, iy: f64, iz: f64, j: f64) -> i64 {
    guard e > 0.0 && nu > -1.0 && area > 0.0 && iy > 0.0 && iz > 0.0 && j > 0.0 else { abort("beam_section: e, area, iy, iz and j > 0, and nu > -1"); }
    beam_sections.push(BeamSection { e: e, g: e / (2.0 * (1.0 + nu)), area: area, iy: iy, iz: iz, j: j, rho: 0.0 });
    return beam_sections.len - 1;
}

fn add_beam(n0: i64, n1: i64, section: i64, y: f64[3], o0: f64[3], o1: f64[3]) -> i64 {
    guard min(n0, n1) >= 0 && max(n0, n1) < node_count() && n0 != n1 && section >= 0 && section < beam_sections.len else { abort("add_beam: two nodes < node_count(), and 0 <= section < beam_sections.len"); }
    beams.push(Beam { n: [n0, n1], section: section, y_axis: y,
                      offset: [o0[0], o0[1], o0[2], o1[0], o1[1], o1[2]] });
    return beams.len - 1;
}

fn fix(n: i64, dof: i64) {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && fixed.len == dof_count() else { abort("fix: 0 <= n < node_count(), and 0 <= dof < 6"); }
    fixed[6 * n + dof] = 1;
}

fn force(n: i64, dof: i64, value: f64) {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && load.len == dof_count() else { abort("force: 0 <= n < node_count(), and 0 <= dof < 6"); }
    load[6 * n + dof] += value;
}

fn edge_load(a: i64, b: i64, fx: f64, fy: f64, fz: f64) {
    guard min(a, b) >= 0 && max(a, b) < node_count() && load.len == dof_count() else { abort("edge_load: 0 <= a, b < node_count()"); }
    let half = 0.5 * norm(sub(node_point(b), node_point(a)));
    let f: f64[3] = [fx, fy, fz];
    for c in 3 {
        load[6 * a + c] += half * f[c];
        load[6 * b + c] += half * f[c];
    }
}

fn edge_load3(a: i64, m: i64, b: i64, fx: f64, fy: f64, fz: f64) {
    guard min(min(a, m), b) >= 0 && max(max(a, m), b) < node_count() && load.len == dof_count() else { abort("edge_load3: 0 <= a, m, b < node_count()"); }
    let pa = node_point(a);
    let pm = node_point(m);
    let pb = node_point(b);
    let f: f64[3] = [fx, fy, fz];
    let gs = [-0.7745966692414834, 0.0, 0.7745966692414834];
    let gw = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
    for g in 3 {
        let s = gs[g];
        let w: f64[3] = [0.5 * s * (s - 1.0), 1.0 - s * s, 0.5 * s * (s + 1.0)];
        var t: f64[3] = [0.0; 3];
        for k in 3 { t[k] = (s - 0.5) * pa[k] - 2.0 * s * pm[k] + (s + 0.5) * pb[k]; }
        let ds = gw[g] * norm(t);
        let nodes = [a, m, b];
        for l in 3 {
            for k in 3 { load[6 * nodes[l] + k] += w[l] * ds * f[k]; }
        }
    }
}

fn clear_supports() { for &v in fixed { v = 0; } }

fn clear_loads() {
    for &v in load { v = 0.0; }
    for &v in disp { v = 0.0; }
}

fn displacement(n: i64, dof: i64) -> f64 {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && disp.len == dof_count() else { abort("displacement: 0 <= n < node_count(), and 0 <= dof < 6"); }
    return disp[6 * n + dof];
}

fn nearest_node(x: f64, y: f64, z: f64) -> i64, f64 {
    var best = -1;
    var d = 0.0;
    for n in node_count() {
        let e = norm(sub(node_point(n), [x, y, z]));
        if best < 0 || e < d {
            best = n;
            d = e;
        }
    }
    return best, d;
}

fn load_sum() -> f64, f64, f64 {
    var s: f64[3] = [0.0, 0.0, 0.0];
    for d in load.len { if d % 6 < 3 { s[d % 6] += load[d]; } }
    return s[0], s[1], s[2];
}

fn dot3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

fn cross(a: f64[3], b: f64[3]) -> f64[3] {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

fn sub(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }

fn norm(a: f64[3]) -> f64 { return sqrt(dot3(a, a)); }

fn unit(a: f64[3]) -> f64[3] {
    let l = norm(a);
    return [a[0] / l, a[1] / l, a[2] / l];
}
