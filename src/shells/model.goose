// The model of the shell runtime: the nodes, the shell and beam cells, the
// supports, the loads and the solution.
//
// Each node has 6 dofs in global axes: the displacement (ux, uy, uz) and the
// rotation vector (rx, ry, rz). Dof d is component d % 6 of node d / 6. A
// shell node and a beam node at the same point are one node, so a beam
// shares the dofs of the skin on its line (spec §12.1, a tie).
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
    n: i64[4],      // the nodes, anticlockwise about the top normal
    section: i64,
    axis: f64[3],   // the 1 axis of the section, projected on the surface
}

struct Beam {
    n: i64[2],
    section: i64,
    y_axis: f64[3],    // the y axis of the section, made normal to the beam
    offset: f64[6],    // from each node to the axis of the beam
}

struct BeamSection {
    e: f64,
    g: f64,
    area: f64,
    iy: f64,     // ∫ z² dA, for bending in the x-z plane
    iz: f64,     // ∫ y² dA, for bending in the x-y plane
    j: f64,      // the torsion constant
}

var coords: f64[>..] = [];          // x, y and z of each node
var directors: f64[>..] = [];       // the unit normal at each node, or 0
var shells: Shell[>..] = [];
var beams: Beam[>..] = [];
var beam_sections: BeamSection[>..] = [];
var fixed: u8[>..] = [];            // 1 for each dof that is held at 0
var load: f64[>..] = [];            // the external force on each dof
var disp: f64[>..] = [];            // the solution

fn node_count() -> i64 { return coords.len / 3; }
fn dof_count() -> i64 { return 6 * node_count(); }

// Removes all of the model, with the sections.
fn clear() {
    coords.clear();
    directors.clear();
    shells.clear();
    beams.clear();
    beam_sections.clear();
    sections.clear();
    grids.clear();
    fixed.clear();
    load.clear();
    disp.clear();
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

// Sets the normal of the shell surface at node n, for example the exact
// normal of a cylinder. A node with no normal gets the mean of the normals
// of its cells.
fn set_director(n: i64, vx: f64, vy: f64, vz: f64) {
    guard n >= 0 && n < node_count() && directors.len == coords.len else { abort("set_director: 0 <= n < node_count(), and a director for each node"); }
    let v = unit([vx, vy, vz]);
    for c in 3 { directors[3 * n + c] = v[c]; }
}

fn add_shell(n0: i64, n1: i64, n2: i64, n3: i64, section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    guard min(min(n0, n1), min(n2, n3)) >= 0 && max(max(n0, n1), max(n2, n3)) < node_count() && section >= 0 && section < sections.len else { abort("add_shell: each node < node_count(), and 0 <= section < sections.len"); }
    shells.push(Shell { n: [n0, n1, n2, n3], section: section, axis: [ax, ay, az] });
    return shells.len - 1;
}

fn beam_section(e: f64, nu: f64, area: f64, iy: f64, iz: f64, j: f64) -> i64 {
    beam_sections.push(BeamSection { e: e, g: e / (2.0 * (1.0 + nu)), area: area, iy: iy, iz: iz, j: j });
    return beam_sections.len - 1;
}

// A beam from node n0 to node n1. y is the y axis of the section, o0 and
// o1 the offsets from n0 and n1 to the axis of the beam.
fn add_beam(n0: i64, n1: i64, section: i64, y: f64[3], o0: f64[3], o1: f64[3]) -> i64 {
    guard min(n0, n1) >= 0 && max(n0, n1) < node_count() && n0 != n1 && section >= 0 && section < beam_sections.len else { abort("add_beam: two nodes < node_count(), and 0 <= section < beam_sections.len"); }
    beams.push(Beam { n: [n0, n1], section: section, y_axis: y,
                      offset: [o0[0], o0[1], o0[2], o1[0], o1[1], o1[2]] });
    return beams.len - 1;
}

// --- supports and loads ------------------------------------------------------

fn fix(n: i64, dof: i64) {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && fixed.len == dof_count() else { abort("fix: 0 <= n < node_count(), and 0 <= dof < 6"); }
    fixed[6 * n + dof] = 1;
}

fn force(n: i64, dof: i64, value: f64) {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && load.len == dof_count() else { abort("force: 0 <= n < node_count(), and 0 <= dof < 6"); }
    load[6 * n + dof] += value;
}

// A force per length (fx, fy, fz) on the straight edge from node a to node b.
fn edge_load(a: i64, b: i64, fx: f64, fy: f64, fz: f64) {
    guard min(a, b) >= 0 && max(a, b) < node_count() && load.len == dof_count() else { abort("edge_load: 0 <= a, b < node_count()"); }
    let half = 0.5 * norm(sub(node_point(b), node_point(a)));
    let f: f64[3] = [fx, fy, fz];
    for c in 3 {
        load[6 * a + c] += half * f[c];
        load[6 * b + c] += half * f[c];
    }
}

// Frees each dof again, for the supports of a new step.
fn clear_supports() { for &v in fixed { v = 0; } }

fn clear_loads() {
    for &v in load { v = 0.0; }
    for &v in disp { v = 0.0; }
}

fn displacement(n: i64, dof: i64) -> f64 {
    guard n >= 0 && n < node_count() && dof >= 0 && dof < 6 && disp.len == dof_count() else { abort("displacement: 0 <= n < node_count(), and 0 <= dof < 6"); }
    return disp[6 * n + dof];
}

// The node nearest to (x, y, z), and its distance.
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

// The sum of the external loads on each dof kind (fx, fy, fz).
fn load_sum() -> f64, f64, f64 {
    var s: f64[3] = [0.0, 0.0, 0.0];
    for d in load.len { if d % 6 < 3 { s[d % 6] += load[d]; } }
    return s[0], s[1], s[2];
}

// --- small vectors -------------------------------------------------------------

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
