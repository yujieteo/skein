// A structured mesh of a cylinder patch, on given lines of x and of θ. A
// point of the patch is (x, R sin θ, R cos θ), so θ turns from +z to +y
// about +x. The cell (i, j) is between the x lines i and i + 1 and the θ
// lines j and j + 1. The cells have their top normal outward and the exact
// normal of the cylinder as the director of each node.
//
// A grid of order 1 has MITC4 cells, and a grid of order 2 has MITC9 cells
// with the mid-side and centre nodes at the middle of x and of θ. The nodes
// are on a lattice of (order nx + 1) × (order nt + 1) points: the lattice
// point (a, b) of cell (i, j) has order i <= a <= order (i + 1) and the
// same for b. Each lattice point has a node, also a point that no cell
// uses (fix_free_nodes holds it).
import std;
import math;
namespace shells;
import .model;

struct Grid {
    first_node: i64,
    nx: i64,            // the cells along x
    nt: i64,            // the cells along θ
    order: i64,         // 1 or 2
    x_fast: bool,       // the node number runs along x first
    r: f64,
    lines: i64,         // the x lines at grid_lines[lines..lines + nx + 1], then the θ lines
    cells: i64,         // the cell (i, j) is grid_cells[cells + i + nx j]
}

var grids: Grid[>..] = [];
var grid_lines: f64[>..] = [];
var grid_cells: i64[>..] = [];     // a shell cell, or -1 for a cell that was left out

// The node of lattice point (a, b) of grid g: a along x, b along θ.
fn grid_node(g: i64, a: i64, b: i64) -> i64 {
    guard g >= 0 && g < grids.len && a >= 0 && a <= grids[g].order * grids[g].nx && b >= 0 && b <= grids[g].order * grids[g].nt else { abort("grid_node: 0 <= g < grids.len, 0 <= a <= order nx and 0 <= b <= order nt"); }
    let gr = grids[g];
    let k = if gr.x_fast { a + (gr.order * gr.nx + 1) * b } else { b + (gr.order * gr.nt + 1) * a };
    return gr.first_node + k;
}

// The cell (i, j) of grid g, or -1 if the cell was left out.
fn grid_cell(g: i64, i: i64, j: i64) -> i64 {
    guard g >= 0 && g < grids.len && i >= 0 && i < grids[g].nx && j >= 0 && j < grids[g].nt else { abort("grid_cell: 0 <= g < grids.len, 0 <= i < nx and 0 <= j < nt"); }
    return grid_cells[grids[g].cells + i + grids[g].nx * j];
}

// The x line i (axis 0) or the θ line i (axis 1) of grid g.
fn grid_line(g: i64, axis: i64, i: i64) -> f64 {
    guard g >= 0 && g < grids.len && (axis == 0 || axis == 1) && i >= 0 && i <= (if axis == 0 { grids[g].nx } else { grids[g].nt }) else { abort("grid_line: 0 <= g < grids.len, axis 0 or 1, and 0 <= i <= nx or nt"); }
    return grid_lines[grids[g].lines + (if axis == 0 { 0 } else { grids[g].nx + 1 }) + i];
}

// The nodes of side `side` of grid g in sequence: side 0 is the first x
// line, 1 the last x line, 2 the first θ line and 3 the last θ line.
fn grid_side(g: i64, side: i64) -> i64[>..] {
    guard g >= 0 && g < grids.len && side >= 0 && side < 4 else { abort("grid_side: 0 <= g < grids.len, and 0 <= side < 4"); }
    let gr = grids[g];
    let la = gr.order * gr.nx;
    let lb = gr.order * gr.nt;
    var out: i64[>..] = [];
    for k in (if side < 2 { lb + 1 } else { la + 1 }) {
        out.push(if side == 0 { grid_node(g, 0, k) } else if side == 1 { grid_node(g, la, k) }
                 else if side == 2 { grid_node(g, k, 0) } else { grid_node(g, k, lb) });
    }
    return out;
}

// The cell i of the lines of grid g on axis `axis` that holds t: the last
// line i with line i <= t, from 0 to the cell count - 1. A binary search.
fn line_cell(g: i64, axis: i64, t: f64) -> i64 {
    var lo = 0;
    var hi = if axis == 0 { grids[g].nx } else { grids[g].nt };
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if grid_line(g, axis, mid) <= t { lo = mid; } else { hi = mid; }
    }
    return lo;
}

// The cell (i, j) of grid g that holds the point (x, y, z), and the place of
// the point in it, from 0 to 1 along x and along θ. False if the point is
// farther than tol from the patch. For a grid of order 2 the place is also
// that of the lattice: 0.5 is the middle node.
fn grid_locate(g: i64, x: f64, y: f64, z: f64, tol: f64) -> i64, i64, f64, f64, bool {
    guard g >= 0 && g < grids.len else { abort("grid_locate: 0 <= g < grids.len"); }
    let gr = grids[g];
    let x0 = grid_line(g, 0, 0);
    let x1 = grid_line(g, 0, gr.nx);
    let t0 = grid_line(g, 1, 0);
    let span = grid_line(g, 1, gr.nt) - t0;
    // θ from t0, in [0, 2π), or a little below 0 near t0.
    var d = atan2(y, z) - t0;
    d -= 2.0 * PI * floor(d / (2.0 * PI));
    if d > span && 2.0 * PI - d < d - span { d -= 2.0 * PI; }
    let i = line_cell(g, 0, x);
    let j = line_cell(g, 1, t0 + d);
    let fx = (x - grid_line(g, 0, i)) / (grid_line(g, 0, i + 1) - grid_line(g, 0, i));
    let ft = (t0 + d - grid_line(g, 1, j)) / (grid_line(g, 1, j + 1) - grid_line(g, 1, j));
    let off_x = max(max(x0 - x, x - x1), 0.0);
    let off_t = gr.r * max(max(-d, d - span), 0.0);
    let off_r = abs(sqrt(y * y + z * z) - gr.r);
    return i, j, fx, ft, off_x <= tol && off_t <= tol && off_r <= tol;
}

// The value of lattice point k on the lines v of a grid of order `order`:
// a line, or the middle of two lines.
fn lattice_value(v: f64[:], k: i64, order: i64) -> f64 {
    guard k >= 0 && k <= order * (v.len - 1) && (order == 1 || order == 2) else { abort("lattice_value: 0 <= k <= order (v.len - 1), and order 1 or 2"); }
    if k % order == 0 { return v[k / order]; }
    return 0.5 * (v[k / 2] + v[k / 2 + 1]);
}

// Makes the cylinder patch of radius r on the x lines xs and the θ lines ts
// (rad), each ascending, with cells of order 1 (MITC4) or 2 (MITC9). Cell
// (i, j) has the section cell_sections[i + nx j], or no cell for -1 (a
// cutout), and the section axis (ax, ay, az). The node numbers run along
// the shorter side of the lattice first, for a small profile.
fn mesh_cylinder_lines(r: f64, xs: f64[:], ts: f64[:], order: i64, cell_sections: i64[:],
                       ax: f64, ay: f64, az: f64) -> i64 {
    let nx = xs.len - 1;
    let nt = ts.len - 1;
    guard r > 0.0 && nx > 0 && nt > 0 && (order == 1 || order == 2) && cell_sections.len == nx * nt else { abort("mesh_cylinder_lines: r > 0, 2 lines or more on each axis, order 1 or 2, and a section for each cell"); }
    for i in nx { guard xs[i + 1] > xs[i] else { abort("mesh_cylinder_lines: xs ascending"); }; }
    for j in nt { guard ts[j + 1] > ts[j] else { abort("mesh_cylinder_lines: ts ascending"); }; }
    for s in cell_sections { guard s >= -1 && s < sections.len else { abort("mesh_cylinder_lines: each section -1 or < sections.len"); }; }
    let la = order * nx + 1;
    let lb = order * nt + 1;
    grids.push(Grid { first_node: node_count(), nx: nx, nt: nt, order: order, x_fast: la <= lb, r: r,
                      lines: grid_lines.len, cells: grid_cells.len });
    let g = grids.len - 1;
    for x in xs { grid_lines.push(x); }
    for t in ts { grid_lines.push(t); }
    for outer in (if la <= lb { lb } else { la }) {
        for inner in (if la <= lb { la } else { lb }) {
            let a = if la <= lb { inner } else { outer };
            let b = if la <= lb { outer } else { inner };
            let t = lattice_value(ts, b, order);
            let n = add_node(lattice_value(xs, a, order), r * sin(t), r * cos(t));
            set_director(n, 0.0, sin(t), cos(t));
        }
    }
    for j in nt {
        for i in nx {
            let s = cell_sections[i + nx * j];
            grid_cells.push(if s < 0 { -1 } else { grid_shell(g, i, j, s, [ax, ay, az]) });
        }
    }
    return g;
}

// The shell cell (i, j) of grid g, with section s and section axis axis.
fn grid_shell(g: i64, i: i64, j: i64, s: i64, axis: f64[3]) -> i64 {
    let o = grids[g].order;
    let a = o * i;
    let b = o * j;
    if o == 1 {
        return add_shell(grid_node(g, a, b), grid_node(g, a + 1, b), grid_node(g, a + 1, b + 1),
                         grid_node(g, a, b + 1), s, axis[0], axis[1], axis[2]);
    }
    return add_shell9([grid_node(g, a, b), grid_node(g, a + 2, b), grid_node(g, a + 2, b + 2), grid_node(g, a, b + 2),
                       grid_node(g, a + 1, b), grid_node(g, a + 2, b + 1), grid_node(g, a + 1, b + 2),
                       grid_node(g, a, b + 1), grid_node(g, a + 1, b + 1)], s, axis[0], axis[1], axis[2]);
}

// Makes the nodes and the MITC4 cells of the patch x0 ≤ x ≤ x1, t0 ≤ θ ≤ t1
// (rad) with nx × nt equal cells, with the section and the section axis
// (ax, ay, az) on each cell.
fn mesh_cylinder(r: f64, x0: f64, x1: f64, t0: f64, t1: f64, nx: i64, nt: i64,
                 section: i64, ax: f64, ay: f64, az: f64) -> i64 {
    guard r > 0.0 && x1 > x0 && t1 > t0 && nx > 0 && nt > 0 && section >= 0 && section < sections.len else { abort("mesh_cylinder: r > 0, x1 > x0, t1 > t0, nx, nt > 0, and 0 <= section < sections.len"); }
    var xs: f64[>..] = [];
    var ts: f64[>..] = [];
    var cs: i64[>..] = [];
    for i in nx + 1 { xs.push(x0 + (x1 - x0) * i / nx); }
    for j in nt + 1 { ts.push(t0 + (t1 - t0) * j / nt); }
    for k in nx * nt { cs.push(section); }
    return mesh_cylinder_lines(r, xs, ts, 1, cs, ax, ay, az);
}
