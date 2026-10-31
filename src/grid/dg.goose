// The discontinuous Galerkin method (spec §17.4) on the grid of
// rectangles: the matrix A and the right side b(t) of A u = b.
//
// The basis of a cell is P_i(ξ) P_j(η), 0 <= i, j <= p, with the Legendre
// polynomials P_i on the reference square [-1, 1]². It is orthogonal, so
// the mass matrix is diagonal. Dof b = i + (p + 1) j of cell c is
// c (p + 1)² + b.
//
// With n the normal from the cell - to the cell +, [w] = w⁻ - w⁺ and
// {w} = (w⁻ + w⁺)/2, the form of ct ddt(u) + div(a u) - d div(grad(u)) = f is
//
//   Σ_K ∫ d ∇u·∇v - (a·∇v) u
//   + Σ_F ∫ F̂ [v] - {d ∇u·n}[v] - {d ∇v·n}[u] + σ [u][v]  =  Σ_K ∫ f v
//
// with the flux F̂ = a·n {u} + α/2 [u]: α = |a·n| (`upwind`, Roe) or |a|
// (`lax_friedrichs`). The diffusion is the symmetric interior penalty form
// (Arnold et al. 2002), with σ = penalty d (p + 1)²/h. On a face with the
// value g, u⁺ = g, and the gradient of the cell is the mean: the condition
// is weak (spec §17.4). A free face has no term.
//
// Each cell of the grid has the same size, so the matrix of each cell and
// of each face of a side is the same: dg_blocks computes them once.
import std;
import math;
namespace grid;
import .mesh;
import .solve;

let MAX_BASIS = 16;    // (MAX_ORDER + 1)²
let MAX_POINTS = 5;    // the Gauss points on each axis: p + 2

// The Gauss-Legendre points and weights on [-1, 1], n = 1 to 5 points:
// those of n points start at GAUSS_FIRST[n - 1].
let GAUSS_FIRST = [0, 1, 3, 6, 10];
let GAUSS_X: f64[15] = [
    0.0,
    -0.5773502691896257, 0.5773502691896257,
    -0.7745966692414834, 0.0, 0.7745966692414834,
    -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
    -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640,
];
let GAUSS_W: f64[15] = [
    2.0,
    1.0, 1.0,
    0.5555555555555556, 0.8888888888888888, 0.5555555555555556,
    0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538,
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891,
];

fn points() -> i64 { return order + 2; }
fn gauss_x(q: i64) -> f64 { return GAUSS_X[GAUSS_FIRST[points() - 1] + q]; }
fn gauss_w(q: i64) -> f64 { return GAUSS_W[GAUSS_FIRST[points() - 1] + q]; }

// P_i(x) and P_i'(x).
fn legendre(i: i64, x: f64) -> f64, f64 {
    if i == 0 { return 1.0, 0.0; }
    if i == 1 { return x, 1.0; }
    if i == 2 { return 0.5 * (3.0 * x * x - 1.0), 3.0 * x; }
    return 0.5 * (5.0 * x * x * x - 3.0 * x), 0.5 * (15.0 * x * x - 3.0);
}

// The basis functions and their x and y derivatives at (ξ, η).
struct Basis {
    v: f64[16],
    gx: f64[16],
    gy: f64[16],
}

fn basis(xi: f64, eta: f64) -> Basis {
    var b = Basis { v: [0.0; 16], gx: [0.0; 16], gy: [0.0; 16] };
    let m = order + 1;
    for j in m {
        let pj, dj = legendre(j, eta);
        for i in m {
            let pi, di = legendre(i, xi);
            b.v[i + m * j] = pi * pj;
            b.gx[i + m * j] = di * pj * 2.0 / hx;
            b.gy[i + m * j] = pi * dj * 2.0 / hy;
        }
    }
    return b;
}

// The basis at the Gauss points of a cell, qi + points() qj, and of each
// side, points() s + q, which basis_tables computes once for the method.
var vol_basis: Basis[>..] = [];
var side_basis: Basis[>..] = [];

fn basis_tables() {
    vol_basis.clear();
    side_basis.clear();
    for qj in points() {
        for qi in points() { vol_basis.push(basis(gauss_x(qi), gauss_x(qj))); }
    }
    for s in 4 {
        for q in points() {
            let p = side_point(s, gauss_x(q));
            side_basis.push(basis(p[0], p[1]));
        }
    }
}

// The diagonal of the mass matrix of dof b of a cell.
fn mass(b: i64) -> f64 {
    let m = order + 1;
    return 0.25 * hx * hy * (2.0 / (2 * (b % m) + 1)) * (2.0 / (2 * (b / m) + 1));
}

// The reference point (ξ, η) of side s at the face coordinate τ.
fn side_point(s: i64, tau: f64) -> f64[2] {
    if s < 2 { return [NORMALS[2 * s], tau]; }
    return [tau, NORMALS[2 * s + 1]];
}

// The value of the flux coefficient α on side s.
fn alpha(s: i64) -> f64 { return if flux == LAX_FRIEDRICHS { sqrt(ax * ax + ay * ay) } else { abs(normal_velocity(s)) }; }

fn sigma(s: i64) -> f64 { return penalty * diffusivity * (order + 1) * (order + 1) / side_distance(s); }

// The matrix of a cell: ∫ d ∇φ_b·∇φ_a - (a·∇φ_a) φ_b, row a and column b.
fn volume_block(out: f64[:]) {
    let n = per_cell();
    for &v in out { v = 0.0; }
    for qj in points() {
        for qi in points() {
            let f = vol_basis[qi + points() * qj];
            let w = gauss_w(qi) * gauss_w(qj) * 0.25 * hx * hy;
            for a in n {
                for b in n {
                    let k = diffusivity * (f.gx[a] * f.gx[b] + f.gy[a] * f.gy[b]) - (ax * f.gx[a] + ay * f.gy[a]) * f.v[b];
                    out[n * a + b] += w * k;
                }
            }
        }
    }
}

// The 4 matrices of an interior face of side s (1 or 3) of a cell: rows
// of the cell - or +, columns of the cell - or +, in the sequence --, -+,
// +-, ++.
fn face_blocks(s: i64, out: f64[:]) {
    let n = per_cell();
    for &v in out { v = 0.0; }
    let an = normal_velocity(s);
    let wm = 0.5 * (an + alpha(s));
    let wp = 0.5 * (an - alpha(s));
    let sg = sigma(s);
    let nx_ = NORMALS[2 * s];
    let ny_ = NORMALS[2 * s + 1];
    let d = 0.5 * diffusivity;
    for q in points() {
        let w = gauss_w(q) * 0.5 * side_length(s);
        let fm = side_basis[points() * s + q];
        let fp = side_basis[points() * opposite_side(s) + q];
        for a in n {
            let vm = fm.v[a];
            let vp = fp.v[a];
            let gvm = fm.gx[a] * nx_ + fm.gy[a] * ny_;
            let gvp = fp.gx[a] * nx_ + fp.gy[a] * ny_;
            for b in n {
                let um = fm.v[b];
                let up = fp.v[b];
                let gum = fm.gx[b] * nx_ + fm.gy[b] * ny_;
                let gup = fp.gx[b] * nx_ + fp.gy[b] * ny_;
                out[n * a + b] += w * (wm * vm * um - d * gum * vm - d * gvm * um + sg * um * vm);
                out[n * n + n * a + b] += w * (wp * vm * up - d * gup * vm + d * gvm * up - sg * up * vm);
                out[2 * n * n + n * a + b] += w * (-wm * vp * um + d * gum * vp - d * gvp * um - sg * um * vp);
                out[3 * n * n + n * a + b] += w * (-wp * vp * up + d * gup * vp + d * gvp * up + sg * up * vp);
            }
        }
    }
}

fn opposite_side(s: i64) -> i64 { return s ^ 1; }

// The matrix of a boundary face of side s with a value.
fn boundary_block(s: i64, out: f64[:]) {
    let n = per_cell();
    for &v in out { v = 0.0; }
    let wm = 0.5 * (normal_velocity(s) + alpha(s));
    let sg = sigma(s);
    for q in points() {
        let w = gauss_w(q) * 0.5 * side_length(s);
        let f = side_basis[points() * s + q];
        for a in n {
            let gv = f.gx[a] * NORMALS[2 * s] + f.gy[a] * NORMALS[2 * s + 1];
            for b in n {
                let gu = f.gx[b] * NORMALS[2 * s] + f.gy[b] * NORMALS[2 * s + 1];
                out[n * a + b] += w * (wm * f.v[a] * f.v[b] - diffusivity * (gu * f.v[a] + gv * f.v[b]) + sg * f.v[a] * f.v[b]);
            }
        }
    }
}

// Adds the n × n block at `at` of `blocks` to the cells row and column.
fn add_block(row: i64, column: i64, blocks: f64[:], at: i64) {
    let n = per_cell();
    for a in n {
        for b in n { add_entry(row * n + a, column * n + b, blocks[at + n * a + b]); }
    }
}

// The matrix A in val.
fn dg_matrix() {
    guard val.len == col.len && row_start.len == dof_count() + 1 else { abort("dg_matrix: the pattern of the dofs"); }
    for &v in val { v = 0.0; }
    let nn = per_cell() * per_cell();
    var vol: f64[256] = [0.0; 256];
    var faces: f64[2048] = [0.0; 2048];      // the 4 blocks of sides 1 and 3
    var bounds: f64[1024] = [0.0; 1024];     // the block of each side on the boundary
    volume_block(vol[..nn]);
    face_blocks(1, faces[..4 * nn]);
    face_blocks(3, faces[4 * nn..8 * nn]);
    for s in 4 { boundary_block(s, bounds[s * nn..(s + 1) * nn]); }
    for c in cell_count() {
        add_block(c, c, vol, 0);
        for s in 4 {
            let nb = neighbour(c, s);
            if nb < 0 {
                if face_bc[4 * c + s] >= 0 { add_block(c, c, bounds, s * nn); }
                continue;
            }
            if s % 2 == 0 { continue; }
            let at = (s / 2) * 4 * nn;
            add_block(c, c, faces, at);
            add_block(c, nb, faces, at + nn);
            add_block(nb, c, faces, at + 2 * nn);
            add_block(nb, nb, faces, at + 3 * nn);
        }
    }
}

// The point of cell c at the reference point (ξ, η).
fn cell_point(c: i64, xi: f64, eta: f64) -> f64[2] {
    let p = centre(c);
    return [p[0] + 0.5 * hx * xi, p[1] + 0.5 * hy * eta];
}

// The right side b(t): ∫ f v, and on each face with a value g the terms
// of g: -∫ d g ∇v·n + σ g v - (a·n - α)/2 g v.
fn dg_rhs<SRC, BC>(t: f64, out: f64[:]) {
    guard out.len == dof_count() else { abort("dg_rhs: one entry for each dof"); }
    let n = per_cell();
    for &v in out { v = 0.0; }
    for c in cell_count() {
        for qj in points() {
            for qi in points() {
                let f = vol_basis[qi + points() * qj];
                let p = cell_point(c, gauss_x(qi), gauss_x(qj));
                let w = gauss_w(qi) * gauss_w(qj) * 0.25 * hx * hy * SRC(p[0], p[1], t);
                for a in n { out[c * n + a] += w * f.v[a]; }
            }
        }
        for s in 4 {
            let k = face_bc[4 * c + s];
            if k < 0 { continue; }
            let wp = 0.5 * (normal_velocity(s) - alpha(s));
            for q in points() {
                let r = side_point(s, gauss_x(q));
                let f = side_basis[points() * s + q];
                let p = cell_point(c, r[0], r[1]);
                let g = BC(k, p[0], p[1], t) * gauss_w(q) * 0.5 * side_length(s);
                for a in n {
                    let gv = f.gx[a] * NORMALS[2 * s] + f.gy[a] * NORMALS[2 * s + 1];
                    out[c * n + a] += g * (-diffusivity * gv + sigma(s) * f.v[a] - wp * f.v[a]);
                }
            }
        }
    }
}

// The L2 projection of the field F(x, y) on each cell, in out.
fn project_field<F>(out: f64[:]) {
    guard out.len == dof_count() && method == DG else { abort("project_field: dg, and one entry for each dof"); }
    let n = per_cell();
    for &v in out { v = 0.0; }
    for c in cell_count() {
        for qj in points() {
            for qi in points() {
                let f = vol_basis[qi + points() * qj];
                let p = cell_point(c, gauss_x(qi), gauss_x(qj));
                let w = gauss_w(qi) * gauss_w(qj) * 0.25 * hx * hy * F(p[0], p[1]);
                for a in n { out[c * n + a] += w * f.v[a] / mass(a); }
            }
        }
    }
}
