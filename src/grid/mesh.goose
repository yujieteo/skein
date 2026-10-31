// The grid of nx × ny rectangles, the equation, the method and the
// boundary faces.
//
// Cell c = i + nx j has the centre (x0 + (i + 1/2) hx, y0 + (j + 1/2) hy).
// Side s of a cell is 0 (-x), 1 (+x), 2 (-y) or 3 (+y), with the outward
// normal NORMALS[s]. A side on the boundary of the rectangle is a face of
// the boundary: `face_bc` gives it the essential condition k >= 0 that
// holds on it, or FREE.
import std;
import math;
namespace grid;
import .dg;

let FV = 0;
let DG = 1;

// The values of the `flux` setting (spec §13).
let UPWIND = 0;
let LINEAR_UPWIND = 1;
let CENTRAL = 2;
let LAX_FRIEDRICHS = 3;
// The values of `gradient`.
let GREEN_GAUSS = 0;
let LEAST_SQUARES = 1;
// The values of `limiter`.
let NONE = 0;
let MINMOD = 1;
let VAN_LEER = 2;
let VENKATAKRISHNAN = 3;

let INTERIOR = -2;   // the face_bc of a side between two cells
let FREE = -1;       // the face_bc of a boundary face with no condition

let NORMALS: f64[8] = [-1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0];
let MAX_ORDER = 3;

var nx = 0;
var ny = 0;
var x0 = 0.0;
var y0 = 0.0;
var hx = 0.0;
var hy = 0.0;

// ct ddt(u) + div(a u) - d div(grad(u)) = f
var ct = 0.0;
var ax = 0.0;
var ay = 0.0;
var diffusivity = 0.0;

var method = -1;
var flux = UPWIND;
var gradient = GREEN_GAUSS;
var limiter = NONE;
var order = 0;        // the polynomial degree of dg on each axis
var penalty = 0.0;    // the factor of the interior penalty of dg

var face_bc: i64[>..] = [];    // 4 for each cell
var selected: u8[>..] = [];    // 4 for each cell: select_faces holds the face
var sol: f64[>..] = [];        // the dofs: a value for each cell (fv), or (order + 1)² for each cell (dg)
var time = 0.0;                // the time of sol

// Meshes the rectangle from (xa, ya) to (xb, yb) with nx × ny cells.
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

// The cell beside side s of cell c, or -1 on the boundary.
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

// The length of side s, and the distance between the centres of two cells
// across it.
fn side_length(s: i64) -> f64 { return if s < 2 { hy } else { hx }; }
fn side_distance(s: i64) -> f64 { return if s < 2 { hx } else { hy }; }

// a·n on side s.
fn normal_velocity(s: i64) -> f64 { return ax * NORMALS[2 * s] + ay * NORMALS[2 * s + 1]; }

// The distance below which two coordinates are the same point.
fn tolerance() -> f64 { return 1.0e-9 * max(nx * hx, ny * hy); }

// Sets the terms of the equation.
fn equation(c_t: f64, a_x: f64, a_y: f64, d: f64) {
    guard cell_count() > 0 && c_t >= 0.0 && d >= 0.0 else { abort("equation: a mesh, c_t >= 0 and d >= 0"); }
    ct = c_t;
    ax = a_x;
    ay = a_y;
    diffusivity = d;
}

// The finite volume method: one value for each cell.
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

// The discontinuous Galerkin method of degree p on each axis.
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

// The dofs of each cell.
fn per_cell() -> i64 { return if method == DG { (order + 1) * (order + 1) } else { 1 }; }

fn dof_count() -> i64 { return cell_count() * per_cell(); }

fn size_state() {
    sol.clear();
    for d in dof_count() { sol.push(0.0); }
    time = 0.0;
    size_solver();
}

// Frees each boundary face.
fn clear_bcs() {
    guard face_bc.len == 4 * cell_count() && method >= 0 else { abort("clear_bcs: a mesh and a method"); }
    for &b in face_bc { if b != INTERIOR { b = FREE; } }
}

// Selects the boundary faces whose centre (x, y) has P(x, y). False if it
// selects none.
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

// Puts the essential condition k on each selected face.
fn fix_selected(k: i64) {
    guard k >= 0 && selected.len == face_bc.len && method >= 0 else { abort("fix_selected: k >= 0, a method and select_faces"); }
    for f in face_bc.len { if selected[f] != 0 { face_bc[f] = k; } }
}

// Sets the state to the field initial(x, y) at time 0: the value at the
// centre of each cell (fv), or the L2 projection on each cell (dg).
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
