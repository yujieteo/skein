// The box mesher: a structured mesh of 27-node hexahedra in the
// unstructured mesh of fem.goose. It is one producer of that mesh. The
// supports, the loads and the values at a point use only the mesh and its
// named sets (solve.goose).
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

var cells: i64[3] = [0, 0, 0];      // cells along x, y and z
var lattice: i64[3] = [0, 0, 0];    // nodes along x, y and z: 2 cells + 1

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
    let lo: f64[3] = [x0, y0, z0];
    let hi: f64[3] = [x1, y1, z1];
    var spacing: f64[3] = [0.0, 0.0, 0.0];   // the size of a cell
    cells = [nx, ny, nz];
    for a in 3 {
        lattice[a] = 2 * cells[a] + 1;
        spacing[a] = (hi[a] - lo[a]) / cells[a];
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
