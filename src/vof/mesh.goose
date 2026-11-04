import std;
import math;
namespace vof;

let INTERIOR = -2;
let OPEN = -1;

let NORMALS: f64[8] = [-1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0];

var nx = 0;
var ny = 0;
var x0 = 0.0;
var y0 = 0.0;
var hx = 0.0;
var hy = 0.0;

var gravity: f64[2] = [0.0, 0.0];
var tension = 0.0;

var side_kind: i64[>..] = [];
var selected: u8[>..] = [];
var side_vx: f64[>..] = [];
var side_vy: f64[>..] = [];

var alpha: f64[>..] = [];
var vel: f64[>..] = [];
var pres: f64[>..] = [];
var uf: f64[>..] = [];
var time = 0.0;

fn cell_count() -> i64 { return nx * ny; }
fn x_face_count() -> i64 { return (nx + 1) * ny; }
fn face_count() -> i64 { return (nx + 1) * ny + nx * (ny + 1); }

fn mesh_rect(xa: f64, ya: f64, xb: f64, yb: f64, cells_x: i64, cells_y: i64) {
    guard xb > xa && yb > ya && cells_x >= 2 && cells_y >= 2 else { abort("mesh_rect: xb > xa, yb > ya, and cells_x, cells_y >= 2"); }
    nx = cells_x;
    ny = cells_y;
    x0 = xa;
    y0 = ya;
    hx = (xb - xa) / cells_x;
    hy = (yb - ya) / cells_y;
    side_kind.clear();
    selected.clear();
    side_vx.clear();
    side_vy.clear();
    for c in cell_count() {
        for s in 4 {
            side_kind.push(if neighbour(c, s) >= 0 { INTERIOR } else { OPEN });
            selected.push(0);
            side_vx.push(0.0);
            side_vy.push(0.0);
        }
    }
    time = 0.0;
    gravity = [0.0, 0.0];
    tension = 0.0;
    size_state();
}

fn size_state() {
    alpha.clear();
    pres.clear();
    for c in cell_count() {
        alpha.push(0.0);
        pres.push(0.0);
    }
    vel.clear();
    for c in 2 * cell_count() { vel.push(0.0); }
    uf.clear();
    for f in face_count() { uf.push(0.0); }
    size_solver();
}

fn neighbour(c: i64, s: i64) -> i64 {
    let i = c % nx;
    let j = c / nx;
    if s == 0 { return if i > 0 { c - 1 } else { -1 }; }
    if s == 1 { return if i + 1 < nx { c + 1 } else { -1 }; }
    if s == 2 { return if j > 0 { c - nx } else { -1 }; }
    return if j + 1 < ny { c + nx } else { -1 };
}

fn face_of(c: i64, s: i64) -> i64 {
    let i = c % nx;
    let j = c / nx;
    if s == 0 { return i + (nx + 1) * j; }
    if s == 1 { return i + 1 + (nx + 1) * j; }
    if s == 2 { return x_face_count() + i + nx * j; }
    return x_face_count() + i + nx * (j + 1);
}

fn outward(s: i64) -> f64 { return if s % 2 == 1 { 1.0 } else { -1.0 }; }

fn centre(c: i64) -> f64[2] { return [x0 + (c % nx + 0.5) * hx, y0 + (c / nx + 0.5) * hy]; }

fn face_centre(c: i64, s: i64) -> f64[2] {
    let p = centre(c);
    return [p[0] + 0.5 * hx * NORMALS[2 * s], p[1] + 0.5 * hy * NORMALS[2 * s + 1]];
}

fn side_length(s: i64) -> f64 { return if s < 2 { hy } else { hx }; }
fn side_distance(s: i64) -> f64 { return if s < 2 { hx } else { hy }; }

fn tolerance() -> f64 { return 1.0e-9 * max(nx * hx, ny * hy); }

fn set_gravity(gx: f64, gy: f64) { gravity = [gx, gy]; }

fn set_tension(sigma: f64) {
    guard sigma >= 0.0 else { abort("set_tension: sigma >= 0"); }
    tension = sigma;
}

fn clear_bcs() {
    guard side_kind.len == 4 * cell_count() else { abort("clear_bcs: a mesh"); }
    for &b in side_kind { if b != INTERIOR { b = OPEN; } }
}

fn select_faces<P>() -> bool {
    guard side_kind.len == 4 * cell_count() && selected.len == side_kind.len else { abort("select_faces: a mesh"); }
    var found = false;
    for c in cell_count() {
        for s in 4 {
            let p = face_centre(c, s);
            let on = side_kind[4 * c + s] != INTERIOR && P(p[0], p[1]);
            selected[4 * c + s] = if on { 1 } else { 0 };
            if on { found = true; }
        }
    }
    return found;
}

fn fix_selected(k: i64) {
    guard k >= 0 && selected.len == side_kind.len else { abort("fix_selected: k >= 0 and select_faces"); }
    for f in side_kind.len { if selected[f] != 0 { side_kind[f] = k; } }
}

fn set_boundary_values<VALUE>(t: f64) {
    guard side_vx.len == side_kind.len else { abort("set_boundary_values: a mesh"); }
    for c in cell_count() {
        for s in 4 {
            let k = side_kind[4 * c + s];
            if k < 0 { continue; }
            let p = face_centre(c, s);
            let v = VALUE(k, p[0], p[1], t);
            side_vx[4 * c + s] = v[0];
            side_vy[4 * c + s] = v[1];
        }
    }
}

fn side_normal_velocity(c: i64, s: i64) -> f64 {
    return if s < 2 { side_vx[4 * c + s] } else { side_vy[4 * c + s] };
}

fn start<INIT>() {
    guard alpha.len == cell_count() else { abort("start: a mesh"); }
    time = 0.0;
    for c in cell_count() {
        let p = centre(c);
        alpha[c] = clamp(INIT(p[0], p[1]), 0.0, 1.0);
        pres[c] = 0.0;
        vel[2 * c] = 0.0;
        vel[2 * c + 1] = 0.0;
    }
    for f in uf.len { uf[f] = 0.0; }
}

fn start_velocity<FIELD>() {
    guard alpha.len == cell_count() && uf.len == face_count() else { abort("start_velocity: a mesh"); }
    for c in cell_count() {
        let p = centre(c);
        let v = FIELD(p[0], p[1]);
        vel[2 * c] = v[0];
        vel[2 * c + 1] = v[1];
        for s in 4 {
            let q = face_centre(c, s);
            let w = FIELD(q[0], q[1]);
            uf[face_of(c, s)] = if s < 2 { w[0] } else { w[1] };
        }
    }
}
