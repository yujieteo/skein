// The grid, the element and the conditions of the package optim. The grid
// has nx x ny equal rectangles from (x0, y0). Node (i, j) is at
// (x0 + i dx, y0 + j dy), and cell (i, j) has the nodes (i, j), (i + 1, j),
// (i + 1, j + 1) and (i, j + 1), counter-clockwise. Cell (i, j) is
// i ny + j, as in the 88 line code. The nodes run along the shorter side
// first, so that the band of the matrix is narrow. Dof 2 n + c is
// component c of node n.
import std;
import math;
namespace optim;

let X = 0;
let Y = 1;
let PLANE_STRESS = 0;
let PLANE_STRAIN = 1;

var x0 = 0.0;
var y0 = 0.0;
var dx = 1.0;
var dy = 1.0;
var nx = 0;
var ny = 0;
var along_y = true;          // node (i, j) is i (ny + 1) + j, else j (nx + 1) + i
var young = 0.0;
var poisson = 0.0;
var thickness = 0.0;
var plane = PLANE_STRESS;
var k0: f64[64];             // the matrix of a cell for E = 1
var density: f64[>..] = [];  // the density of each cell, 1 without a topology loop
var penalty = 1.0;           // the SIMP power p
let E_MIN = 1.0e-9;          // the modulus of a void cell, a fraction of E (Andreassen et al. 2011)

// A load of a step: a force (fx, fy) at a point, or a traction (fx, fy)
// on the boundary edges on the grid line `axis == value`, which acts on
// the thickness.
struct Load {
    point: bool,
    axis: i64,
    at: f64[2],
    f: f64[2],
}

var fixed: u8[>..] = [];
var fixed_value: f64[>..] = [];
var loads: Load[>..] = [];
var u: f64[>..] = [];        // the solution, with the fixed values

fn node_count() -> i64 { return (nx + 1) * (ny + 1); }
fn dof_count() -> i64 { return 2 * node_count(); }
fn cell_count() -> i64 { return nx * ny; }
fn node(i: i64, j: i64) -> i64 { return if along_y { i * (ny + 1) + j } else { j * (nx + 1) + i }; }

// The nodes of cell c, counter-clockwise from its corner at the lowest x and y.
fn cell_nodes(c: i64) -> i64[4] {
    let i = c / ny;
    let j = c % ny;
    return [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)];
}

// The grid of cx x cy rectangles from (xa, ya) to (xb, yb).
fn mesh_rect(xa: f64, ya: f64, xb: f64, yb: f64, cx: i64, cy: i64) {
    guard xb > xa && yb > ya && cx >= 1 && cy >= 1 else { abort("mesh_rect: xb > xa, yb > ya, and cx, cy >= 1"); }
    x0 = xa;
    y0 = ya;
    nx = cx;
    ny = cy;
    dx = (xb - xa) / cx;
    dy = (yb - ya) / cy;
    along_y = ny <= nx;
    density.clear();
    for c in cell_count() { density.push(1.0); }
    penalty = 1.0;
    clear_bcs();
}

// The material and the plane assumption: E, nu and the thickness t.
fn elastic(E: f64, nu: f64, t: f64, which: i64) {
    guard nx > 0 && E > 0.0 && nu > -1.0 && nu < 0.5 && t > 0.0 && (which == PLANE_STRESS || which == PLANE_STRAIN) else { abort("elastic: a mesh, E > 0, -1 < nu < 0.5, t > 0 and a plane assumption"); }
    young = E;
    poisson = nu;
    thickness = t;
    plane = which;
    k0 = cell_matrix(d_matrix(1.0, nu));
}

// The matrix D of the stresses (xx, yy, xy) from the strains (xx, yy, 2 xy).
fn d_matrix(E: f64, nu: f64) -> f64[9] {
    if plane == PLANE_STRAIN {
        let f = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        return [f * (1.0 - nu), f * nu, 0.0, f * nu, f * (1.0 - nu), 0.0, 0.0, 0.0, f * (0.5 - nu)];
    }
    let f = E / (1.0 - nu * nu);
    return [f, f * nu, 0.0, f * nu, f, 0.0, 0.0, 0.0, f * 0.5 * (1.0 - nu)];
}

// The derivative of D by nu, for the gradient (adjoint.goose).
fn d_by_nu(E: f64, nu: f64) -> f64[9] {
    let h = 1.0e-6;
    let a = d_matrix(E, nu + h);
    let b = d_matrix(E, nu - h);
    var d: f64[9] = default<f64[9]>();
    for k in 9 { d[k] = (a[k] - b[k]) / (2.0 * h); }
    return d;
}

// The matrix B (3 x 8) of a cell at (s, t) in [-1, 1]^2.
fn b_matrix(s: f64, t: f64) -> f64[24] {
    let ds: f64[4] = [-(1.0 - t), 1.0 - t, 1.0 + t, -(1.0 + t)];
    let dt: f64[4] = [-(1.0 - s), -(1.0 + s), 1.0 + s, 1.0 - s];
    var b: f64[24] = default<f64[24]>();
    for a in 4 {
        let gx = 0.5 * ds[a] / dx;
        let gy = 0.5 * dt[a] / dy;
        b[2 * a] = gx;
        b[8 + 2 * a + 1] = gy;
        b[16 + 2 * a] = gy;
        b[16 + 2 * a + 1] = gx;
    }
    return b;
}

// The matrix t ∫ Bᵀ D B of a cell, by 2 x 2 Gauss points.
fn cell_matrix(d: f64[9]) -> f64[64] {
    var k: f64[64] = default<f64[64]>();
    let g = 1.0 / sqrt(3.0);
    let w = 0.25 * dx * dy * thickness;
    for q in 4 {
        let b = b_matrix(if q == 1 || q == 2 { g } else { -g }, if q >= 2 { g } else { -g });
        for r in 8 {
            let db: f64[3] = [d[0] * b[r] + d[1] * b[8 + r] + d[2] * b[16 + r],
                              d[3] * b[r] + d[4] * b[8 + r] + d[5] * b[16 + r],
                              d[6] * b[r] + d[7] * b[8 + r] + d[8] * b[16 + r]];
            for c in 8 { k[8 * c + r] += w * (b[c] * db[0] + b[8 + c] * db[1] + b[16 + c] * db[2]); }
        }
    }
    return k;
}

// The modulus of cell c: E (E_MIN + (1 - E_MIN) x^p).
fn cell_modulus(c: i64) -> f64 { return young * (E_MIN + (1.0 - E_MIN) * pow(density[c], penalty)); }

// The 8 dofs of cell c.
fn cell_dofs(c: i64) -> i64[8] {
    let n = cell_nodes(c);
    return [2 * n[0], 2 * n[0] + 1, 2 * n[1], 2 * n[1] + 1, 2 * n[2], 2 * n[2] + 1, 2 * n[3], 2 * n[3] + 1];
}

// aᵀ k b on the dofs of cell c, for the matrix k of a cell.
fn cell_product(k: f64[64], c: i64, a: f64[:], b: f64[:]) -> f64 {
    let e = cell_dofs(c);
    var s = 0.0;
    for r in 8 {
        var kb = 0.0;
        for q in 8 { kb += k[8 * q + r] * b[e[q]]; }
        s += a[e[r]] * kb;
    }
    return s;
}

// Removes the conditions and the loads of the last step.
fn clear_bcs() {
    fixed.clear();
    fixed_value.clear();
    for d in dof_count() {
        fixed.push(0);
        fixed_value.push(0.0);
    }
    loads.clear();
}

// The tolerance of a point on a grid line: 1e-6 of the size of the grid.
fn tolerance() -> f64 { return 1.0e-6 * max(nx * dx, ny * dy); }

// Fixes component comp of each node on the grid line axis == value at v.
// False if no node is on it.
fn fix_line(axis: i64, value: f64, comp: i64, v: f64) -> bool {
    guard nx > 0 && (axis == X || axis == Y) && (comp == X || comp == Y) && fixed.len == dof_count() else { abort("fix_line: a mesh, clear_bcs, and an axis and a component X or Y"); }
    var found = false;
    for i in nx + 1 {
        for j in ny + 1 {
            let at = if axis == X { x0 + i * dx } else { y0 + j * dy };
            guard abs(at - value) <= tolerance() else { continue; }
            fixed[2 * node(i, j) + comp] = 1;
            fixed_value[2 * node(i, j) + comp] = v;
            found = true;
        }
    }
    return found;
}

// The node at (x, y), or -1 if no node is there.
fn node_at(x: f64, y: f64) -> i64 {
    let i = round((x - x0) / dx) as! i64;
    let j = round((y - y0) / dy) as! i64;
    guard i >= 0 && i <= nx && j >= 0 && j <= ny && abs(x0 + i * dx - x) <= tolerance() && abs(y0 + j * dy - y) <= tolerance() else { return -1; }
    return node(i, j);
}

// Fixes component comp of the node at (x, y) at v. False if no node is there.
fn fix_point(x: f64, y: f64, comp: i64, v: f64) -> bool {
    guard nx > 0 && (comp == X || comp == Y) && fixed.len == dof_count() else { abort("fix_point: a mesh, clear_bcs, and a component X or Y"); }
    let n = node_at(x, y);
    guard n >= 0 else { return false; }
    fixed[2 * n + comp] = 1;
    fixed_value[2 * n + comp] = v;
    return true;
}

// The force (fx, fy) at the point (x, y) of the grid. False if the point
// is outside the grid.
fn point_load(x: f64, y: f64, fx: f64, fy: f64) -> bool {
    guard nx > 0 && fixed.len == dof_count() else { abort("point_load: a mesh and clear_bcs"); }
    let _, _, found = point_weights(x, y);
    loads.push(Load { point: true, axis: -1, at: [x, y], f: [fx, fy] });
    return found;
}

// The traction (tx, ty) on the boundary edges on the line axis == value.
// False if the line is not a side of the grid.
fn edge_load(axis: i64, value: f64, tx: f64, ty: f64) -> bool {
    guard nx > 0 && (axis == X || axis == Y) && fixed.len == dof_count() else { abort("edge_load: a mesh, clear_bcs and an axis X or Y"); }
    let lo = if axis == X { x0 } else { y0 };
    let hi = lo + (if axis == X { nx * dx } else { ny * dy });
    loads.push(Load { point: false, axis: axis, at: [value, 0.0], f: [tx, ty] });
    return abs(value - lo) <= tolerance() || abs(value - hi) <= tolerance();
}

// The cell of the point (x, y), its nodes and their bilinear weights. A
// point on the side of two cells takes the cell of lower index. found is
// false outside the grid.
fn point_weights(x: f64, y: f64) -> i64[4], f64[4], bool {
    let fx = (x - x0) / dx;
    let fy = (y - y0) / dy;
    let tol = tolerance() / min(dx, dy);
    guard fx >= -tol && fx <= nx + tol && fy >= -tol && fy <= ny + tol else { return [0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0], false; }
    let i = clamp(floor(fx) as! i64, 0, nx - 1);
    let j = clamp(floor(fy) as! i64, 0, ny - 1);
    let s = clamp(fx - i, 0.0, 1.0);
    let t = clamp(fy - j, 0.0, 1.0);
    return cell_nodes(i * ny + j), [(1.0 - s) * (1.0 - t), s * (1.0 - t), s * t, (1.0 - s) * t], true;
}

// Adds the nodal forces of load l, times `scale`, to f.
fn add_load(l: Load, scale: f64, f: f64[>..]&) {
    if l.point {
        let n, w, _ = point_weights(l.at[0], l.at[1]);
        for k in 4 { for c in 2 { f[2 * n[k] + c] += scale * w[k] * l.f[c]; } }
        return;
    }
    // Each edge of length h on the line gives h t / 2 of the traction to each of its 2 nodes.
    let count = if l.axis == X { ny } else { nx };
    let h = if l.axis == X { dy } else { dx };
    let i = round((l.at[0] - (if l.axis == X { x0 } else { y0 })) / (if l.axis == X { dx } else { dy })) as! i64;
    for e in count {
        for end in 2 {
            let n = if l.axis == X { node(i, e + end) } else { node(e + end, i) };
            for c in 2 { f[2 * n + c] += scale * 0.5 * h * thickness * l.f[c]; }
        }
    }
}

// The value of component comp of u at (x, y), bilinear in its cell.
fn value_at(comp: i64, x: f64, y: f64) -> f64, bool {
    guard u.len == dof_count() && (comp == X || comp == Y) else { abort("value_at: a solve, and a component X or Y"); }
    let n, w, found = point_weights(x, y);
    var v = 0.0;
    for k in 4 { v += w[k] * u[2 * n[k] + comp]; }
    return v, found;
}
