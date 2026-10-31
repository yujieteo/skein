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
    r: f64,
    x0: f64,
    x1: f64,
    t0: f64,
    t1: f64,
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

// The nodes of side `side` of grid g in sequence: side 0 is x = x0, 1 is
// x = x1, 2 is θ = t0 and 3 is θ = t1.
fn grid_side(g: i64, side: i64) -> i64[>..] {
    let gr = grids[g];
    var out: i64[>..] = [];
    for k in (if side < 2 { gr.nt + 1 } else { gr.nx + 1 }) {
        out.push(if side == 0 { grid_node(g, 0, k) } else if side == 1 { grid_node(g, gr.nx, k) }
                 else if side == 2 { grid_node(g, k, 0) } else { grid_node(g, k, gr.nt) });
    }
    return out;
}

// The cell (i, j) of grid g that holds the point (x, y, z), and the place of
// the point in it, from 0 to 1 along x and along θ. False if the point is
// farther than tol from the patch.
fn grid_locate(g: i64, x: f64, y: f64, z: f64, tol: f64) -> i64, i64, f64, f64, bool {
    let gr = grids[g];
    let span = gr.t1 - gr.t0;
    // θ from t0, in [0, 2π), or a little below 0 near t0.
    var d = atan2(y, z) - gr.t0;
    d -= 2.0 * PI * floor(d / (2.0 * PI));
    if d > span && 2.0 * PI - d < d - span { d -= 2.0 * PI; }
    let u = (x - gr.x0) / (gr.x1 - gr.x0) * gr.nx;
    let v = d / span * gr.nt;
    let i = min(max(floor(u) as! i64, 0), gr.nx - 1);
    let j = min(max(floor(v) as! i64, 0), gr.nt - 1);
    let fx = u - i;
    let ft = v - j;
    let off_x = max(max(gr.x0 - x, x - gr.x1), 0.0);
    let off_t = gr.r * max(max(-d, d - span), 0.0);
    let off_r = abs(sqrt(y * y + z * z) - gr.r);
    return i, j, fx, ft, off_x <= tol && off_t <= tol && off_r <= tol;
}

// Makes the nodes and the cells of the patch x0 ≤ x ≤ x1, t0 ≤ θ ≤ t1 (rad),
// with the section and the section axis (ax, ay, az) on each cell. The
// node numbers run along the shorter side first, for a small profile.
fn mesh_cylinder(r: f64, x0: f64, x1: f64, t0: f64, t1: f64, nx: i64, nt: i64,
                 section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    assert(r > 0.0 && x1 > x0 && t1 > t0 && nx > 0 && nt > 0 && section >= 0 && section < sections.len);
    grids.push(Grid { first_node: node_count(), first_cell: shells.len, nx: nx, nt: nt, x_fast: nx <= nt,
                      r: r, x0: x0, x1: x1, t0: t0, t1: t1 });
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
