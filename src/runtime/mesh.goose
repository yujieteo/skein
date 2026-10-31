// The lattice mesher: a structured mesh of 27-node hexahedra of a box or
// of a sector of an annulus, or of 8-node hexahedra with incompatible
// modes of a box, in the unstructured mesh of fem.goose. It is one producer
// of that mesh. The supports, the loads and the values at a point use only
// the mesh and its named sets (solve.goose).
//
// The nodes are the points of a lattice with 2n + 1 points (hex27) or
// n + 1 points (hex8i) along an edge of n cells, so node (i, j, k) of the
// lattice has the index i + lx (j + ly k).
// Cell (cx, cy, cz) has the index cx + nx (cy + ny cz).
//
// The lattice is uniform in its own coordinates (q0, q1, q2). For a box
// they are x, y and z. For an annulus sector they are r, θ and z, and the
// node of (r, θ, z) is at (cx + r cos θ, cy + r sin θ, z). The cells
// are isoparametric, so a cell edge on an arc is a parabola through 3
// nodes on the arc.
//
// Each face of the lattice has a node set and a face set: x0, x1, y0, y1,
// z0 and z1 on a box, r0, r1, t0, t1, z0 and z1 on an annulus (the low and
// the high face along each lattice axis).
import std;
import math;
namespace runtime;
import .fem;
import .solid;
import .solve;

var cells: i64[3] = [0, 0, 0];      // cells along x, y and z
var lattice: i64[3] = [0, 0, 0];    // nodes along x, y and z: 2 cells + 1 (hex27) or cells + 1

let FACE_NAMES = ["x0", "x1", "y0", "y1", "z0", "z1"];
let ANNULUS_FACE_NAMES = ["r0", "r1", "t0", "t1", "z0", "z1"];

fn node_at(i: i64, j: i64, k: i64) -> i64 {
    guard i >= 0 && j >= 0 && k >= 0 && i < lattice[0] && j < lattice[1] && k < lattice[2] else { abort("node_at: a point of the lattice"); }
    return i + lattice[0] * (j + lattice[1] * k);
}

fn box_cell(cx: i64, cy: i64, cz: i64) -> i64 {
    guard cx >= 0 && cy >= 0 && cz >= 0 && cx < cells[0] && cy < cells[1] && cz < cells[2] else { abort("box_cell: a cell of the lattice"); }
    return cx + cells[0] * (cy + cells[1] * cz);
}

// t is HEX27 or HEX8I.
fn build_box(t: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             nx: i64, ny: i64, nz: i64) {
    guard t == HEX27 || t == HEX8I else { abort("build_box: t is HEX27 or HEX8I"); }
    build_lattice(t, false, 0.0, 0.0, [x0, y0, z0], [x1, y1, z1], nx, ny, nz);
}

fn build_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                 nr: i64, nt: i64, nz: i64) {
    guard r0 > 0.0 else { abort("build_annulus: r0 > 0"); }
    build_lattice(HEX27, true, cx, cy, [r0, t0, z0], [r1, t1, z1], nr, nt, nz);
}

fn build_lattice(t: i64, annulus: bool, cx: f64, cy: f64, lo: f64[3], hi: f64[3], nx: i64, ny: i64, nz: i64) {
    guard nx > 0 && ny > 0 && nz > 0 && hi[0] > lo[0] && hi[1] > lo[1] && hi[2] > lo[2] else { abort("build_lattice: counts > 0 and hi > lo on each axis"); }
    let per = if t == HEX27 { 2 } else { 1 };   // lattice steps along a cell edge
    mesh_clear();
    var spacing: f64[3] = [0.0, 0.0, 0.0];   // the size of a cell in q
    cells = [nx, ny, nz];
    for a in 3 {
        lattice[a] = per * cells[a] + 1;
        spacing[a] = (hi[a] - lo[a]) / cells[a];
    }
    for k in lattice[2] {
        for j in lattice[1] {
            for i in lattice[0] {
                let a = lo[0] + i * spacing[0] / per;
                let b = lo[1] + j * spacing[1] / per;
                let c = lo[2] + k * spacing[2] / per;
                if annulus {
                    add_node(cx + a * cos(b), cy + a * sin(b), c);
                } else {
                    add_node(a, b, c);
                }
            }
        }
    }
    var nodes: i64[27] = [0; 27];
    for cz in nz {
        for cy in ny {
            for cx in nx {
                if t == HEX27 {
                    for l in 27 { nodes[l] = node_at(2 * cx + l % 3, 2 * cy + l / 3 % 3, 2 * cz + l / 9); }
                    add_cell(HEX27, nodes);
                } else {
                    // The corners in the order of HEX8_CORNERS (hex8i.goose).
                    var corners: i64[8] = [0; 8];
                    for l in 8 { corners[l] = node_at(cx + (l + 1) / 2 % 2, cy + l / 2 % 2, cz + l / 4); }
                    add_cell(HEX8I, corners);
                }
            }
        }
    }
    add_lattice_face_sets(annulus);
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
}

// The node set and the face set of each of the 6 faces of the lattice.
fn add_lattice_face_sets(annulus: bool) {
    guard cells[0] > 0 && cells[1] > 0 && cells[2] > 0 else { abort("add_lattice_face_sets: a lattice"); }
    let names = if annulus { ANNULUS_FACE_NAMES } else { FACE_NAMES };
    var members: i64[>..] = [];
    for axis in 3 {
        let a, b = other_axes(axis);
        for side in 2 {
            let f = side * (lattice[axis] - 1);
            members.clear();
            for t in lattice[b] {
                for s in lattice[a] { members.push(face_node(axis, f, s, t)); }
            }
            add_node_set(names[2 * axis + side], members);
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
            add_face_set(names[2 * axis + side], members);
        }
    }
}

// The node at position (s, t) of the lattice of the face that is normal to
// `axis` at lattice index `f`. s and t run over the two other axes in the
// sequence x, y, z.
fn face_node(axis: i64, f: i64, s: i64, t: i64) -> i64 {
    guard axis >= 0 && axis < 3 else { abort("face_node: 0 <= axis < 3"); }
    if axis == 0 { return node_at(f, s, t); }
    if axis == 1 { return node_at(s, f, t); }
    return node_at(s, t, f);
}

fn other_axes(axis: i64) -> i64, i64 {
    guard axis >= 0 && axis < 3 else { abort("other_axes: 0 <= axis < 3"); }
    if axis == 0 { return 1, 2; }
    if axis == 1 { return 0, 2; }
    return 0, 1;
}

// --- the grid of a union of boxes ----------------------------------------------

// The grid mesher: 27-node hexahedra on the grid lines lx, ly and lz, which
// need not be uniform. Cell (i, j, k) spans lx[i] to lx[i + 1], and so on.
// It is in the mesh if its centre is in a box of `boxes`. The nodes are the
// grid lines and the midpoints between them, and only the nodes of a cell
// of the mesh are made.
fn build_grid(lx: f64[:], ly: f64[:], lz: f64[:], boxes: f64[:]) {
    guard lx.len >= 2 && ly.len >= 2 && lz.len >= 2 && boxes.len > 0 && boxes.len % 6 == 0 else { abort("build_grid: 2 grid lines on each axis and a box"); }
    mesh_clear();
    let px = with_midpoints(lx);
    let py = with_midpoints(ly);
    let pz = with_midpoints(lz);
    cells = [lx.len - 1, ly.len - 1, lz.len - 1];
    lattice = [px.len, py.len, pz.len];
    var ids: i64[>..] = [];
    for i in px.len * py.len * pz.len { ids.push(-1); }
    var grid: i64[>..] = [];   // the cell of each cell of the grid, or -1
    var nodes: i64[27] = [0; 27];
    for k in cells[2] {
        for j in cells[1] {
            for i in cells[0] {
                let centre: f64[3] = [0.5 * (lx[i] + lx[i + 1]), 0.5 * (ly[j] + ly[j + 1]), 0.5 * (lz[k] + lz[k + 1])];
                if !in_boxes(boxes, centre) {
                    grid.push(-1);
                    continue;
                }
                for l in 27 {
                    let at = node_at(2 * i + l % 3, 2 * j + l / 3 % 3, 2 * k + l / 9);
                    if ids[at] < 0 { ids[at] = add_node(px[2 * i + l % 3], py[2 * j + l / 3 % 3], pz[2 * k + l / 9]); }
                    nodes[l] = ids[at];
                }
                grid.push(add_cell(HEX27, nodes));
            }
        }
    }
    add_grid_face_sets(grid);
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
}

// The grid lines `l` and the midpoint between each two of them.
fn with_midpoints(l: f64[:]) -> f64[>..] {
    guard l.len >= 2 else { abort("with_midpoints: 2 grid lines"); }
    var out: f64[>..] = [l[0]];
    for i in 1..l.len {
        guard l[i] > l[i - 1] else { abort("with_midpoints: grid lines in increasing order"); }
        out.push(0.5 * (l[i - 1] + l[i]));
        out.push(l[i]);
    }
    return out;
}

// Whether the point p is in one of the boxes (x0, y0, z0, x1, y1, z1).
fn in_boxes(boxes: f64[:], p: f64[3]) -> bool {
    guard boxes.len % 6 == 0 else { abort("in_boxes: 6 numbers for each box"); }
    for b in boxes.len / 6 {
        var all = true;
        for a in 3 {
            if p[a] < boxes[6 * b + a] || p[a] > boxes[6 * b + 3 + a] { all = false; }
        }
        if all { return true; }
    }
    return false;
}

// The node set and the face set of each bounding plane of the grid: the
// faces there of the cells of the mesh. `grid` holds the cell of each cell
// of the grid, or -1.
fn add_grid_face_sets(grid: i64[:]) {
    guard grid.len == cells[0] * cells[1] * cells[2] else { abort("add_grid_face_sets: a cell or -1 for each cell of the grid"); }
    var faces: i64[>..] = [];
    var nodes: i64[>..] = [];
    var mark: bool[>..] = [];
    for axis in 3 {
        for side in 2 {
            faces.clear();
            nodes.clear();
            mark.clear();
            for n in node_count() { mark.push(false); }
            for g, at in grid {
                let index: i64[3] = [at % cells[0], at / cells[0] % cells[1], at / (cells[0] * cells[1])];
                if g < 0 || index[axis] != side * (cells[axis] - 1) { continue; }
                faces.push(FACE_STRIDE * g + 2 * axis + side);
                for l in 27 {
                    let local: i64[3] = [l % 3, l / 3 % 3, l / 9];
                    let n = cell_node(g, l);
                    if local[axis] == 2 * side && !mark[n] {
                        mark[n] = true;
                        nodes.push(n);
                    }
                }
            }
            add_node_set(FACE_NAMES[2 * axis + side], nodes);
            add_face_set(FACE_NAMES[2 * axis + side], faces);
        }
    }
}

// The centre of cell c: the mean of its nodes.
fn cell_centre(c: i64) -> f64[3] {
    guard c >= 0 && c < cell_count() else { abort("cell_centre: a cell"); }
    var p: f64[3] = [0.0, 0.0, 0.0];
    let nn = cell_start[c + 1] - cell_start[c];
    for l in nn {
        let n = cell_node(c, l);
        for s in 3 { p[s] += coords[3 * n + s] / nn; }
    }
    return p;
}

// Gives E, nu and the thermal strain to each cell whose centre is in the
// box lo to hi, and gives the number of these cells.
fn set_box_material(E: f64, nu: f64, strain: f64, lo: f64[3], hi: f64[3]) -> i64 {
    guard cell_young.len == cell_count() && cell_poisson.len == cell_count() && cell_strain.len == cell_count() else { abort("set_box_material: a material for each cell"); }
    var count = 0;
    for c in cell_count() {
        if !in_boxes([lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]], cell_centre(c)) { continue; }
        cell_young[c] = E;
        cell_poisson[c] = nu;
        cell_strain[c] = strain;
        count++;
    }
    return count;
}

// Adds the node set `name` of the node at p, to 1e-9 of the size of the
// mesh. False if no node is there.
fn add_point_set(name: const u8[:], p: f64[3]) -> bool {
    guard node_count() > 0 else { abort("add_point_set: a mesh"); }
    var lo: f64[3] = [coords[0], coords[1], coords[2]];
    var hi: f64[3] = lo;
    for n in node_count() {
        for s in 3 {
            lo[s] = min(lo[s], coords[3 * n + s]);
            hi[s] = max(hi[s], coords[3 * n + s]);
        }
    }
    let tol = 1.0e-9 * max(hi[0] - lo[0], max(hi[1] - lo[1], hi[2] - lo[2]));
    for n in node_count() {
        if abs(coords[3 * n] - p[0]) <= tol && abs(coords[3 * n + 1] - p[1]) <= tol && abs(coords[3 * n + 2] - p[2]) <= tol {
            add_node_set(name, [n]);
            return true;
        }
    }
    return false;
}
