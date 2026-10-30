// The structured mesh of a box or of a sector of an annulus, the supports
// and the loads.
//
// The cells are 27-node hexahedra. Their nodes are the points of a lattice
// with 2n + 1 points along an edge of n cells, so node (i, j, k) of the
// lattice has the index i + lx (j + ly k). The local node l of a cell is
// at lattice offset (l % 3, l / 3 % 3, l / 9) from the first corner of the
// cell. Each dof d is component d % 3 of node d / 3.
//
// The lattice is uniform in its own coordinates (q0, q1, q2). For a box
// they are x, y and z. For an annulus sector they are r, θ and z, and the
// node of (r, θ, z) is at (cx + r cos θ, cy + r sin θ, z). The cells
// are isoparametric, so a cell edge on an arc is a parabola through 3
// nodes on the arc.
import std;
import math;
namespace runtime;

// --- the mesh ----------------------------------------------------------------

var annulus = false;
var centre: f64[2] = [0.0, 0.0];   // the axis of an annulus
var lattice_lo: f64[3] = [0.0, 0.0, 0.0];
var lattice_hi: f64[3] = [0.0, 0.0, 0.0];
var cells: i64[3] = [0, 0, 0];      // cells along q0, q1 and q2
var lattice: i64[3] = [0, 0, 0];    // nodes along q0, q1 and q2: 2 cells + 1
var spacing: f64[3] = [0.0, 0.0, 0.0];   // the size of a cell in q
var coords: f64[>..] = [];          // x, y and z of each node

// --- the material, the supports, the loads and the solution ------------------

var young = 0.0;
var poisson = 0.0;
var fixed: u8[>..] = [];            // 1 for each dof that is held at 0
var load: f64[>..] = [];            // the external force on each dof
var disp: f64[>..] = [];            // the displacement of each dof

fn node_count() -> i64 { return lattice[0] * lattice[1] * lattice[2]; }
fn cell_count() -> i64 { return cells[0] * cells[1] * cells[2]; }

fn node_at(i: i64, j: i64, k: i64) -> i64 {
    return i + lattice[0] * (j + lattice[1] * k);
}

// The node of local index l in cell (cx, cy, cz).
fn cell_node(cx: i64, cy: i64, cz: i64, l: i64) -> i64 {
    return node_at(2 * cx + l % 3, 2 * cy + l / 3 % 3, 2 * cz + l / 9);
}

fn build_lattice(q0: f64, q1: f64, q2: f64, p0: f64, p1: f64, p2: f64,
                 n0: i64, n1: i64, n2: i64) {
    assert(n0 > 0 && n1 > 0 && n2 > 0);
    lattice_lo = [q0, q1, q2];
    lattice_hi = [p0, p1, p2];
    cells = [n0, n1, n2];
    for a in 3 {
        lattice[a] = 2 * cells[a] + 1;
        spacing[a] = (lattice_hi[a] - lattice_lo[a]) / cells[a];
    }
    coords.clear();
    for k in lattice[2] {
        for j in lattice[1] {
            for i in lattice[0] {
                let a = q0 + 0.5 * i * spacing[0];
                let b = q1 + 0.5 * j * spacing[1];
                let c = q2 + 0.5 * k * spacing[2];
                if annulus {
                    coords.push(centre[0] + a * cos(b));
                    coords.push(centre[1] + a * sin(b));
                } else {
                    coords.push(a);
                    coords.push(b);
                }
                coords.push(c);
            }
        }
    }
    young = 0.0;
    poisson = 0.0;
    fixed.clear();
    load.clear();
    disp.clear();
    for d in 3 * node_count() {
        fixed.push(0);
        load.push(0.0);
        disp.push(0.0);
    }
    clear_matrix();
}

fn build_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             nx: i64, ny: i64, nz: i64) {
    annulus = false;
    build_lattice(x0, y0, z0, x1, y1, z1, nx, ny, nz);
}

fn build_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                 nr: i64, nt: i64, nz: i64) {
    assert(r0 > 0.0);
    annulus = true;
    centre = [cx, cy];
    build_lattice(r0, t0, z0, r1, t1, z1, nr, nt, nz);
}

// --- faces -------------------------------------------------------------------

// The lattice index of the face where lattice coordinate `axis` is
// `value`, or -1 if it is not a face of the lattice.
fn face_index(axis: i64, value: f64) -> i64 {
    guard axis >= 0 && axis < 3 && cells[0] > 0 else { return -1; }
    let tol = 1.0e-9 * max(abs(lattice_hi[axis] - lattice_lo[axis]), 1.0e-300);
    if abs(value - lattice_lo[axis]) <= tol { return 0; }
    if abs(value - lattice_hi[axis]) <= tol { return lattice[axis] - 1; }
    return -1;
}

// The node at position (s, t) of the lattice of the face that is normal to
// `axis` at lattice index `f`. s and t run over the two other axes in the
// sequence q0, q1, q2.
fn face_node(axis: i64, f: i64, s: i64, t: i64) -> i64 {
    if axis == 0 { return node_at(f, s, t); }
    if axis == 1 { return node_at(s, f, t); }
    return node_at(s, t, f);
}

fn other_axes(axis: i64) -> i64, i64 {
    if axis == 0 { return 1, 2; }
    if axis == 1 { return 0, 2; }
    return 0, 1;
}

// Holds the components first to last - 1 of each node of a face at 0.
fn mark_fixed(axis: i64, value: f64, first: i64, last: i64) -> bool {
    let f = face_index(axis, value);
    guard f >= 0 else { return false; }
    let a, b = other_axes(axis);
    for t in lattice[b] {
        for s in lattice[a] {
            let n = face_node(axis, f, s, t);
            for c in first..last { fixed[3 * n + c] = 1; }
        }
    }
    clear_solution();
    return true;
}

// The consistent nodal forces of the traction (tx, ty, tz) + s n on a
// face, with n the unit normal out of the body: the integral of each shape
// function over each cell face, by 3 × 3 Gauss points.
fn add_face_load(axis: i64, value: f64, tx: f64, ty: f64, tz: f64, s: f64) -> bool {
    let f = face_index(axis, value);
    guard f >= 0 else { return false; }
    let a, b = other_axes(axis);
    // The tangents along a and b give the normal along +axis for the axes
    // 0 and 2, and along -axis for axis 1. The face at f = 0 looks to -axis.
    var sign = if axis == 1 { -1.0 } else { 1.0 };
    if f == 0 { sign = -sign; }
    for cb in cells[b] {
        for ca in cells[a] {
            for gb in 3 {
                for ga in 3 {
                    var da: f64[3] = [0.0, 0.0, 0.0];
                    var db: f64[3] = [0.0, 0.0, 0.0];
                    for l in 9 {
                        let n = face_node(axis, f, 2 * ca + l % 3, 2 * cb + l / 3);
                        let wa = shape_1d_ds(l % 3, GAUSS_X[ga]) * shape_1d(l / 3, GAUSS_X[gb]);
                        let wb = shape_1d(l % 3, GAUSS_X[ga]) * shape_1d_ds(l / 3, GAUSS_X[gb]);
                        for c in 3 {
                            da[c] += wa * coords[3 * n + c];
                            db[c] += wb * coords[3 * n + c];
                        }
                    }
                    // The normal times the area of the face per unit area
                    // of the reference square.
                    let nx = sign * (da[1] * db[2] - da[2] * db[1]);
                    let ny = sign * (da[2] * db[0] - da[0] * db[2]);
                    let nz = sign * (da[0] * db[1] - da[1] * db[0]);
                    let area = sqrt(nx * nx + ny * ny + nz * nz);
                    let w = GAUSS_W[ga] * GAUSS_W[gb];
                    for l in 9 {
                        let n = face_node(axis, f, 2 * ca + l % 3, 2 * cb + l / 3);
                        let wl = w * shape_1d(l % 3, GAUSS_X[ga]) * shape_1d(l / 3, GAUSS_X[gb]);
                        load[3 * n] += wl * (tx * area + s * nx);
                        load[3 * n + 1] += wl * (ty * area + s * ny);
                        load[3 * n + 2] += wl * (tz * area + s * nz);
                    }
                }
            }
        }
    }
    clear_solution();
    return true;
}

fn add_point_force(x: f64, y: f64, z: f64, fx: f64, fy: f64, fz: f64) {
    var best = 0;
    var best_d = -1.0;
    for n in node_count() {
        let dx = coords[3 * n] - x;
        let dy = coords[3 * n + 1] - y;
        let dz = coords[3 * n + 2] - z;
        let d = dx * dx + dy * dy + dz * dz;
        if best_d < 0.0 || d < best_d { best = n; best_d = d; }
    }
    load[3 * best] += fx;
    load[3 * best + 1] += fy;
    load[3 * best + 2] += fz;
    clear_solution();
}

fn clear_load_vector() {
    for &v in load { v = 0.0; }
    clear_solution();
}

fn clear_solution() {
    for &v in disp { v = 0.0; }
}

fn load_sum() -> f64, f64, f64 {
    var s: f64[3] = [0.0, 0.0, 0.0];
    for v, d in load { s[d % 3] += v; }
    return s[0], s[1], s[2];
}

// --- shape functions -----------------------------------------------------------

let GAUSS_X: f64[3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
let GAUSS_W: f64[3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

// The quadratic Lagrange function of node a (at -1, 0 and 1) at s.
fn shape_1d(a: i64, s: f64) -> f64 {
    if a == 0 { return 0.5 * s * (s - 1.0); }
    if a == 1 { return 1.0 - s * s; }
    return 0.5 * s * (s + 1.0);
}

fn shape_1d_ds(a: i64, s: f64) -> f64 {
    if a == 0 { return s - 0.5; }
    if a == 1 { return -2.0 * s; }
    return s + 0.5;
}

// --- points ------------------------------------------------------------------

// The cell that holds the last point that `find_point` found, and the local
// coordinates (-1 to 1) of the point in it.
var point_cell: i64[3] = [0, 0, 0];
var point_local: f64[3] = [0.0, 0.0, 0.0];

// Finds the cell that holds (x, y, z). False if the point is outside the
// mesh. In an annulus the local coordinates come from r and θ of the
// point, not from the parabolic edges of the cell. The two differ by less
// than 1e-6 of the radius for a cell of 2° or less.
fn find_point(x: f64, y: f64, z: f64) -> bool {
    guard cells[0] > 0 else { return false; }
    var q: f64[3] = [x, y, z];
    if annulus {
        let dx = x - centre[0];
        let dy = y - centre[1];
        q[0] = sqrt(dx * dx + dy * dy);
        q[1] = atan2(dy, dx);
    }
    for a in 3 {
        let tol = 1.0e-9 * abs(lattice_hi[a] - lattice_lo[a]);
        guard q[a] >= lattice_lo[a] - tol && q[a] <= lattice_hi[a] + tol else { return false; }
        let c = clamp(floor((q[a] - lattice_lo[a]) / spacing[a]) as i64, 0, cells[a] - 1);
        point_cell[a] = c;
        point_local[a] = 2.0 * (q[a] - lattice_lo[a] - c * spacing[a]) / spacing[a] - 1.0;
    }
    return true;
}

fn interpolate_u(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    guard disp.len > 0 && find_point(x, y, z) else { return 0.0, 0.0, 0.0; }
    var u: f64[3] = [0.0, 0.0, 0.0];
    for l in 27 {
        let w = shape_1d(l % 3, point_local[0]) * shape_1d(l / 3 % 3, point_local[1])
            * shape_1d(l / 9, point_local[2]);
        let n = cell_node(point_cell[0], point_cell[1], point_cell[2], l);
        for c in 3 { u[c] += w * disp[3 * n + c]; }
    }
    return u[0], u[1], u[2];
}
