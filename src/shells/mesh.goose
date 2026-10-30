// A structured mesh of a cylinder patch: (nx + 1) × (nt + 1) nodes at
// (x, R sin θ, R cos θ), so θ turns from +z to +y about +x. The cells have
// their top normal outward and the exact normal of the cylinder as the
// director of each node.
import std;
import math;
namespace shells;
import .model;

struct Grid {
    first_node: i64,
    first_cell: i64,
    nx: i64,
    nt: i64,
    x_fast: bool,       // the node number runs along x first
}

var grids: Grid[>..] = [];

// The node (i, j) of grid g: i along x, j along θ.
fn grid_node(g: i64, i: i64, j: i64) -> i64 {
    let gr = grids[g];
    let k = if gr.x_fast { i + (gr.nx + 1) * j } else { j + (gr.nt + 1) * i };
    return gr.first_node + k;
}

// The cell (i, j) of grid g, or -1 if the cell was left out.
fn grid_cell(g: i64, i: i64, j: i64) -> i64 { return grids[g].first_cell + i + grids[g].nx * j; }

// Makes the nodes and the cells of the patch x0 ≤ x ≤ x1, t0 ≤ θ ≤ t1 (rad),
// with the section and the section axis (ax, ay, az) on each cell. The
// node numbers run along the shorter side first, for a small profile.
fn mesh_cylinder(r: f64, x0: f64, x1: f64, t0: f64, t1: f64, nx: i64, nt: i64,
                 section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    grids.push(Grid { first_node: node_count(), first_cell: shells.len, nx: nx, nt: nt, x_fast: nx <= nt });
    let g = grids.len - 1;
    let outer = if nx <= nt { nt + 1 } else { nx + 1 };
    let inner = if nx <= nt { nx + 1 } else { nt + 1 };
    for a in outer {
        for b in inner {
            let i = if nx <= nt { b } else { a };
            let j = if nx <= nt { a } else { b };
            let x = x0 + (x1 - x0) * i / nx;
            let t = t0 + (t1 - t0) * j / nt;
            let n = add_node(x, r * sin(t), r * cos(t));
            set_director(n, 0.0, sin(t), cos(t));
        }
    }
    for j in nt {
        for i in nx {
            add_shell(grid_node(g, i, j), grid_node(g, i + 1, j), grid_node(g, i + 1, j + 1),
                      grid_node(g, i, j + 1), section, ax, ay, az);
        }
    }
    return g;
}
