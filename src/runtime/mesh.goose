// The box mesher: a structured mesh of 27-node hexahedra in the
// unstructured mesh of fem.goose, the supports and the loads on the faces
// of the box.
//
// The nodes are the points of a lattice with 2n + 1 points along an edge of
// n cells, so node (i, j, k) of the lattice has the index i + lx (j + ly k).
// Cell (cx, cy, cz) has the index cx + nx (cy + ny cz). Each face of the
// box has a node set and a face set, named x0, x1, y0, y1, z0 and z1 (the
// low and the high face along each axis).
import std;
import math;
namespace runtime;
import .fem;
import .solid;
import .solve;

var box_lo: f64[3] = [0.0, 0.0, 0.0];
var box_hi: f64[3] = [0.0, 0.0, 0.0];
var cells: i64[3] = [0, 0, 0];      // cells along x, y and z
var lattice: i64[3] = [0, 0, 0];    // nodes along x, y and z: 2 cells + 1
var spacing: f64[3] = [0.0, 0.0, 0.0];   // the size of a cell

let FACE_NAMES = ["x0", "x1", "y0", "y1", "z0", "z1"];

fn node_at(i: i64, j: i64, k: i64) -> i64 {
    return i + lattice[0] * (j + lattice[1] * k);
}

fn box_cell(cx: i64, cy: i64, cz: i64) -> i64 {
    return cx + cells[0] * (cy + cells[1] * cz);
}

fn build_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             nx: i64, ny: i64, nz: i64) {
    assert(nx > 0 && ny > 0 && nz > 0);
    mesh_clear();
    box_lo = [x0, y0, z0];
    box_hi = [x1, y1, z1];
    cells = [nx, ny, nz];
    for a in 3 {
        lattice[a] = 2 * cells[a] + 1;
        spacing[a] = (box_hi[a] - box_lo[a]) / cells[a];
    }
    for k in lattice[2] {
        for j in lattice[1] {
            for i in lattice[0] {
                add_node(x0 + 0.5 * i * spacing[0], y0 + 0.5 * j * spacing[1], z0 + 0.5 * k * spacing[2]);
            }
        }
    }
    var nodes: i64[27] = [0; 27];
    for cz in nz {
        for cy in ny {
            for cx in nx {
                for l in 27 { nodes[l] = node_at(2 * cx + l % 3, 2 * cy + l / 3 % 3, 2 * cz + l / 9); }
                add_cell(HEX27, nodes);
            }
        }
    }
    var members: i64[>..] = [];
    for axis in 3 {
        let a, b = other_axes(axis);
        for side in 2 {
            let f = side * (lattice[axis] - 1);
            members.clear();
            for t in lattice[b] {
                for s in lattice[a] { members.push(face_node(axis, f, s, t)); }
            }
            add_node_set(FACE_NAMES[2 * axis + side], members);
            members.clear();
            let cf = side * (cells[axis] - 1);
            for tb in cells[b] {
                for ta in cells[a] {
                    var c: i64[3] = [0, 0, 0];
                    c[axis] = cf;
                    c[a] = ta;
                    c[b] = tb;
                    members.push(FACE_STRIDE * box_cell(c[0], c[1], c[2]) + 2 * axis + side);
                }
            }
            add_face_set(FACE_NAMES[2 * axis + side], members);
        }
    }
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
}

// --- faces -------------------------------------------------------------------

// The lattice index of the face where coordinate `axis` is `value`, or -1
// if the plane is not a face of the box.
fn face_index(axis: i64, value: f64) -> i64 {
    guard axis >= 0 && axis < 3 && cells[0] > 0 else { return -1; }
    let tol = 1.0e-9 * max(abs(box_hi[axis] - box_lo[axis]), 1.0e-300);
    if abs(value - box_lo[axis]) <= tol { return 0; }
    if abs(value - box_hi[axis]) <= tol { return lattice[axis] - 1; }
    return -1;
}

// The node at position (s, t) of the lattice of the face that is normal to
// `axis` at lattice index `f`. s and t run over the two other axes in the
// sequence x, y, z.
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

// The node set of the face where coordinate `axis` is `value`, or -1.
fn face_set(axis: i64, value: f64) -> i64 {
    let f = face_index(axis, value);
    guard f >= 0 else { return -1; }
    return find_set(FACE_NAMES[2 * axis + (if f == 0 { 0 } else { 1 })], false);
}

fn mark_fixed(axis: i64, value: f64) -> bool {
    let s = face_set(axis, value);
    guard s >= 0 else { return false; }
    fix_set(s);
    return true;
}

fn face_reaction(axis: i64, value: f64) -> f64, f64, f64 {
    let s = face_set(axis, value);
    guard s >= 0 else { return 0.0, 0.0, 0.0; }
    return set_reaction(s);
}

// The consistent nodal forces of a uniform traction: the integral of each
// shape function over each cell face, by 3 × 3 Gauss points.
fn add_traction(axis: i64, value: f64, tx: f64, ty: f64, tz: f64) -> bool {
    let f = face_index(axis, value);
    guard f >= 0 else { return false; }
    let a, b = other_axes(axis);
    let area = 0.25 * spacing[a] * spacing[b];   // the Jacobian of a face
    for cb in cells[b] {
        for ca in cells[a] {
            for l in 9 {
                let la = l % 3;
                let lb = l / 3;
                var w = 0.0;
                for gb in 3 {
                    for ga in 3 {
                        w += GAUSS3_W[ga] * GAUSS3_W[gb]
                            * shape_1d(la, GAUSS3_X[ga]) * shape_1d(lb, GAUSS3_X[gb]);
                    }
                }
                let n = face_node(axis, f, 2 * ca + la, 2 * cb + lb);
                load[3 * n] += w * area * tx;
                load[3 * n + 1] += w * area * ty;
                load[3 * n + 2] += w * area * tz;
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

// --- points ---------------------------------------------------------------------

// The cell that holds coordinate v along axis a, and the local coordinate
// of v in it (-1 to 1).
fn locate(a: i64, v: f64) -> i64, f64 {
    let c = clamp(floor((v - box_lo[a]) / spacing[a]) as i64, 0, cells[a] - 1);
    let s = 2.0 * (v - box_lo[a] - c * spacing[a]) / spacing[a] - 1.0;
    return c, s;
}

fn interpolate_u(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    guard disp.len > 0 else { return 0.0, 0.0, 0.0; }
    let p: f64[3] = [x, y, z];
    for a in 3 {
        let tol = 1.0e-9 * abs(box_hi[a] - box_lo[a]);
        guard p[a] >= box_lo[a] - tol && p[a] <= box_hi[a] + tol else { return 0.0, 0.0, 0.0; }
    }
    let cx, sx = locate(0, x);
    let cy, sy = locate(1, y);
    let cz, sz = locate(2, z);
    var u: f64[3] = [0.0, 0.0, 0.0];
    for l in 27 {
        let w = shape_1d(l % 3, sx) * shape_1d(l / 3 % 3, sy) * shape_1d(l / 9, sz);
        let n = cell_node(box_cell(cx, cy, cz), l);
        for c in 3 { u[c] += w * disp[3 * n + c]; }
    }
    return u[0], u[1], u[2];
}
