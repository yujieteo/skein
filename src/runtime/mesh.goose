// The lattice mesher: a structured mesh of 27-node hexahedra of a box or
// of a sector of an annulus in the unstructured mesh of fem.goose. It is
// one producer of that mesh. The supports, the loads and the values at a
// point use only the mesh and its named sets (solve.goose).
//
// The nodes are the points of a lattice with 2n + 1 points along an edge of
// n cells, so node (i, j, k) of the lattice has the index i + lx (j + ly k).
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
var lattice: i64[3] = [0, 0, 0];    // nodes along x, y and z: 2 cells + 1

let FACE_NAMES = ["x0", "x1", "y0", "y1", "z0", "z1"];
let ANNULUS_FACE_NAMES = ["r0", "r1", "t0", "t1", "z0", "z1"];

fn node_at(i: i64, j: i64, k: i64) -> i64 {
    return i + lattice[0] * (j + lattice[1] * k);
}

fn box_cell(cx: i64, cy: i64, cz: i64) -> i64 {
    return cx + cells[0] * (cy + cells[1] * cz);
}

fn build_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64,
             nx: i64, ny: i64, nz: i64) {
    build_lattice(false, 0.0, 0.0, [x0, y0, z0], [x1, y1, z1], nx, ny, nz);
}

fn build_annulus(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64, z0: f64, z1: f64,
                 nr: i64, nt: i64, nz: i64) {
    assert(r0 > 0.0);
    build_lattice(true, cx, cy, [r0, t0, z0], [r1, t1, z1], nr, nt, nz);
}

fn build_lattice(annulus: bool, cx: f64, cy: f64, lo: f64[3], hi: f64[3], nx: i64, ny: i64, nz: i64) {
    assert(nx > 0 && ny > 0 && nz > 0);
    mesh_clear();
    var spacing: f64[3] = [0.0, 0.0, 0.0];   // the size of a cell in q
    cells = [nx, ny, nz];
    for a in 3 {
        lattice[a] = 2 * cells[a] + 1;
        spacing[a] = (hi[a] - lo[a]) / cells[a];
    }
    for k in lattice[2] {
        for j in lattice[1] {
            for i in lattice[0] {
                let a = lo[0] + 0.5 * i * spacing[0];
                let b = lo[1] + 0.5 * j * spacing[1];
                let c = lo[2] + 0.5 * k * spacing[2];
                if annulus {
                    add_node(cx + a * cos(b), cy + a * sin(b), c);
                } else {
                    add_node(a, b, c);
                }
            }
        }
    }
    let names = if annulus { ANNULUS_FACE_NAMES } else { FACE_NAMES };
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
    young = 0.0;
    poisson = 0.0;
    clear_dofs();
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
