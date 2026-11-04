import std;
import math;
namespace grid;
import .dg;

let FV = 0;
let DG = 1;

let UPWIND = 0;
let LINEAR_UPWIND = 1;
let CENTRAL = 2;
let LAX_FRIEDRICHS = 3;
let GREEN_GAUSS = 0;
let LEAST_SQUARES = 1;
let NONE = 0;
let MINMOD = 1;
let VAN_LEER = 2;
let VENKATAKRISHNAN = 3;

let INTERIOR = -2;
let FREE = -1;

let NORMALS: f64[8] = [-1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0];
let MAX_ORDER = 3;

var nx = 0;
var ny = 0;
var x0 = 0.0;
var y0 = 0.0;
var hx = 0.0;
var hy = 0.0;

var ct = 0.0;
var ax = 0.0;
var ay = 0.0;
var diffusivity = 0.0;

var method = -1;
var flux = UPWIND;
var gradient = GREEN_GAUSS;
var limiter = NONE;
var order = 0;
var penalty = 0.0;

var face_bc: i64[>..] = [];
var selected: u8[>..] = [];
var sol: f64[>..] = [];
var time = 0.0;

fn mesh_rect(xa: f64, ya: f64, xb: f64, yb: f64, cells_x: i64, cells_y: i64) {
    guard xb > xa && yb > ya && cells_x >= 1 && cells_y >= 1 else { abort("mesh_rect: xb > xa, yb > ya, and cells_x, cells_y >= 1"); }
    nx = cells_x;
    ny = cells_y;
    x0 = xa;
    y0 = ya;
    hx = (xb - xa) / cells_x;
    hy = (yb - ya) / cells_y;
    face_bc.clear();
    selected.clear();
    for c in cell_count() {
        for s in 4 {
            face_bc.push(if neighbour(c, s) >= 0 { INTERIOR } else { FREE });
            selected.push(0);
        }
    }
    method = -1;
    sol.clear();
    time = 0.0;
}

fn cell_count() -> i64 { return nx * ny; }

fn neighbour(c: i64, s: i64) -> i64 {
    let i = c % nx;
    let j = c / nx;
    if s == 0 { return if i > 0 { c - 1 } else { -1 }; }
    if s == 1 { return if i + 1 < nx { c + 1 } else { -1 }; }
    if s == 2 { return if j > 0 { c - nx } else { -1 }; }
    return if j + 1 < ny { c + nx } else { -1 };
}

fn centre(c: i64) -> f64[2] { return [x0 + (c % nx + 0.5) * hx, y0 + (c / nx + 0.5) * hy]; }

fn face_centre(c: i64, s: i64) -> f64[2] {
    let p = centre(c);
    return [p[0] + 0.5 * hx * NORMALS[2 * s], p[1] + 0.5 * hy * NORMALS[2 * s + 1]];
}

fn side_length(s: i64) -> f64 { return if s < 2 { hy } else { hx }; }
fn side_distance(s: i64) -> f64 { return if s < 2 { hx } else { hy }; }

fn normal_velocity(s: i64) -> f64 { return ax * NORMALS[2 * s] + ay * NORMALS[2 * s + 1]; }

fn tolerance() -> f64 { return 1.0e-9 * max(nx * hx, ny * hy); }

fn equation(c_t: f64, a_x: f64, a_y: f64, d: f64) {
    guard cell_count() > 0 && c_t >= 0.0 && d >= 0.0 else { abort("equation: a mesh, c_t >= 0 and d >= 0"); }
    ct = c_t;
    ax = a_x;
    ay = a_y;
    diffusivity = d;
}

fn fv(face_flux: i64, face_gradient: i64, face_limiter: i64) {
    guard cell_count() > 0 && face_flux >= UPWIND && face_flux <= CENTRAL && face_gradient >= GREEN_GAUSS && face_gradient <= LEAST_SQUARES &&
          face_limiter >= NONE && face_limiter <= VENKATAKRISHNAN else { abort("fv: a mesh, and a flux, a gradient and a limiter of fv"); }
    method = FV;
    flux = face_flux;
    gradient = face_gradient;
    limiter = face_limiter;
    order = 0;
    size_state();
}

fn dg(p: i64, face_flux: i64, sigma: f64) {
    guard cell_count() > 0 && p >= 1 && p <= MAX_ORDER && (face_flux == UPWIND || face_flux == LAX_FRIEDRICHS) && sigma > 0.0 else {
        abort("dg: a mesh, 1 <= p <= 3, the flux upwind or lax_friedrichs, and sigma > 0");
    }
    method = DG;
    flux = face_flux;
    order = p;
    penalty = sigma;
    basis_tables();
    size_state();
}

fn per_cell() -> i64 { return if method == DG { (order + 1) * (order + 1) } else { 1 }; }

fn dof_count() -> i64 { return cell_count() * per_cell(); }

fn size_state() {
    sol.clear();
    for d in dof_count() { sol.push(0.0); }
    time = 0.0;
    size_solver();
}

fn clear_bcs() {
    guard face_bc.len == 4 * cell_count() && method >= 0 else { abort("clear_bcs: a mesh and a method"); }
    for &b in face_bc { if b != INTERIOR { b = FREE; } }
}

fn select_faces<P>() -> bool {
    guard face_bc.len == 4 * cell_count() && selected.len == face_bc.len else { abort("select_faces: a mesh"); }
    var found = false;
    for c in cell_count() {
        for s in 4 {
            let p = face_centre(c, s);
            let on = face_bc[4 * c + s] != INTERIOR && P(p[0], p[1]);
            selected[4 * c + s] = if on { 1 } else { 0 };
            if on { found = true; }
        }
    }
    return found;
}

fn fix_selected(k: i64) {
    guard k >= 0 && selected.len == face_bc.len && method >= 0 else { abort("fix_selected: k >= 0, a method and select_faces"); }
    for f in face_bc.len { if selected[f] != 0 { face_bc[f] = k; } }
}

fn start<INIT>() {
    guard method >= 0 && sol.len == dof_count() else { abort("start: a method"); }
    time = 0.0;
    if method == FV {
        for c in cell_count() {
            let p = centre(c);
            sol[c] = INIT(p[0], p[1]);
        }
        return;
    }
    project_field(sol, INIT);
}
