// The cell-centred finite volume method (spec §17.3): the residual of each
// cell and the matrix of the defect correction.
//
// The residual of cell c is the outward flux through its sides less the
// source, ∮ (a u - d ∇u)·n - f(x_c, t) V. On a face between two cells:
//
//   - the diffusive flux is -d (u_N - u_C)/h, which is second order on the
//     grid of rectangles, and -d (g - u_C)/(h/2) on a face with the value g;
//   - the advective flux is a·n u_f, with u_f from the upwind cell C:
//     `upwind` u_C, `central` (u_C + u_D)/2, and `linear_upwind`
//     u_C + ψ(r) (u_D - u_C)/2 with r = 2 ∇u_C·d_CD/(u_D - u_C) - 1
//     (Darwish and Moukalled 2003). The limiter `none` is ψ = (1 + r)/2,
//     which gives u_C + ∇u_C·(x_f - x_C); `minmod` is max(0, min(1, r))
//     and `van_leer` (r + |r|)/(1 + |r|), which are TVD (Sweby 1984).
//     `venkatakrishnan` gives u_C + φ_C ∇u_C·(x_f - x_C) with the limiter
//     φ_C of the cell (Venkatakrishnan 1995).
//
// A boundary face with the value g is a ghost cell at the mirror of C,
// with the value 2g - u_C: u_f is g at an inflow face, and the formula
// above at an outflow face. A free face has no flux.
//
// The gradient of a cell is `green_gauss`, (1/V) Σ u_f n A with the face
// value the mean of the two cells or g, or `least_squares` over the
// neighbours and the faces with a value, weighted by 1/|d|².
import std;
import math;
namespace grid;
import .mesh;
import .solve;

let VENKAT_K = 1.0;   // ε² = (K h/ℓ)³ (u_max - u_min)² of the limiter of Venkatakrishnan

// The value of the condition of side s of cell c at time t.
fn face_value<BC>(c: i64, s: i64, t: f64) -> f64 {
    let p = face_centre(c, s);
    return BC(face_bc[4 * c + s], p[0], p[1], t);
}

// The gradient of u in each cell, in grad_u.
fn fv_gradients<BC>(u: f64[:], t: f64) {
    guard u.len == cell_count() && grad_u.len == 2 * u.len else { abort("fv_gradients: u and grad_u of the cells"); }
    for c in cell_count() {
        var a: f64[3] = [0.0; 3];    // the normal matrix of least squares: xx, xy, yy
        var g: f64[2] = [0.0; 2];
        for s in 4 {
            let nb = neighbour(c, s);
            let b = face_bc[4 * c + s];
            guard nb >= 0 || b >= 0 else { if gradient == GREEN_GAUSS { for k in 2 { g[k] += u[c] * NORMALS[2 * s + k] * side_length(s); } } continue; };
            let uf = if nb >= 0 { u[nb] } else { face_value(c, s, t, BC) };
            let h = if nb >= 0 { side_distance(s) } else { 0.5 * side_distance(s) };
            if gradient == GREEN_GAUSS {
                let um = if nb >= 0 { 0.5 * (u[c] + uf) } else { uf };
                for k in 2 { g[k] += um * NORMALS[2 * s + k] * side_length(s); }
                continue;
            }
            let dx = h * NORMALS[2 * s];
            let dy = h * NORMALS[2 * s + 1];
            let w = 1.0 / (h * h);
            a[0] += w * dx * dx;
            a[1] += w * dx * dy;
            a[2] += w * dy * dy;
            g[0] += w * dx * (uf - u[c]);
            g[1] += w * dy * (uf - u[c]);
        }
        grad_u[2 * c] = 0.0;
        grad_u[2 * c + 1] = 0.0;
        if gradient == GREEN_GAUSS {
            grad_u[2 * c] = g[0] / (hx * hy);
            grad_u[2 * c + 1] = g[1] / (hx * hy);
            continue;
        }
        let det = a[0] * a[2] - a[1] * a[1];
        if det > 1.0e-12 * (a[0] * a[0] + a[2] * a[2]) {
            grad_u[2 * c] = (a[2] * g[0] - a[1] * g[1]) / det;
            grad_u[2 * c + 1] = (a[0] * g[1] - a[1] * g[0]) / det;
        } else {
            if a[0] > 0.0 { grad_u[2 * c] = g[0] / a[0]; }
            if a[2] > 0.0 { grad_u[2 * c + 1] = g[1] / a[2]; }
        }
    }
}

// The limiter of Venkatakrishnan of each cell, in phi.
fn fv_venkat<BC>(u: f64[:], t: f64) {
    guard u.len == cell_count() && phi.len == u.len && grad_u.len == 2 * u.len else { abort("fv_venkat: u, phi and grad_u of the cells"); }
    var lo = u[0];
    var hi = u[0];
    for v in u {
        lo = min(lo, v);
        hi = max(hi, v);
    }
    let l = sqrt(nx * hx * ny * hy);
    let kh = VENKAT_K * sqrt(hx * hy) / l;
    let e2 = max(kh * kh * kh * (hi - lo) * (hi - lo), 1.0e-300);
    for c in cell_count() {
        var umax = u[c];
        var umin = u[c];
        for s in 4 {
            let nb = neighbour(c, s);
            guard nb >= 0 || face_bc[4 * c + s] >= 0 else { continue; };
            // A face with a value is the ghost cell of the flux.
            let v = if nb >= 0 { u[nb] } else { 2.0 * face_value(c, s, t, BC) - u[c] };
            umax = max(umax, v);
            umin = min(umin, v);
        }
        var f = 1.0;
        for s in 4 {
            let d2 = 0.5 * side_distance(s) * (grad_u[2 * c] * NORMALS[2 * s] + grad_u[2 * c + 1] * NORMALS[2 * s + 1]);
            if d2 == 0.0 { continue; }
            let d1 = if d2 > 0.0 { umax - u[c] } else { umin - u[c] };
            f = min(f, (d1 * d1 + e2 + 2.0 * d1 * d2) / (d1 * d1 + 2.0 * d2 * d2 + d1 * d2 + e2));
        }
        phi[c] = f;
    }
}

// ψ(r) of the limiter.
fn psi(r: f64) -> f64 {
    if limiter == MINMOD { return max(0.0, min(1.0, r)); }
    if limiter == VAN_LEER { return (r + abs(r)) / (1.0 + abs(r)); }
    return 0.5 * (1.0 + r);
}

// The advected value on the face of the upwind cell c, with the value ud
// of the cell (or the ghost) downwind at the distance h along normal s.
fn upwind_value(u: f64[:], c: i64, s: i64, ud: f64, h: f64) -> f64 {
    if flux == UPWIND { return u[c]; }
    if flux == CENTRAL { return 0.5 * (u[c] + ud); }
    let g = h * (grad_u[2 * c] * NORMALS[2 * s] + grad_u[2 * c + 1] * NORMALS[2 * s + 1]);
    if limiter == VENKATAKRISHNAN { return u[c] + 0.5 * phi[c] * g; }
    if limiter == NONE { return u[c] + 0.5 * g; }
    let du = ud - u[c];
    if abs(du) <= 1.0e-14 * (abs(u[c]) + abs(ud)) { return u[c]; }
    return u[c] + 0.5 * psi(2.0 * g / du - 1.0) * du;
}

// The opposite side of side s.
fn opposite(s: i64) -> i64 { return s ^ 1; }

// The outward flux of cell c through its side s to the cell nb.
fn interior_flux(u: f64[:], c: i64, s: i64, nb: i64) -> f64 {
    let an = normal_velocity(s);
    let h = side_distance(s);
    let uf = if an >= 0.0 { upwind_value(u, c, s, u[nb], h) } else { upwind_value(u, nb, opposite(s), u[c], h) };
    return (an * uf - diffusivity * (u[nb] - u[c]) / h) * side_length(s);
}

// The outward flux of cell c through its boundary side s, with the value g.
fn boundary_flux(u: f64[:], c: i64, s: i64, g: f64) -> f64 {
    let an = normal_velocity(s);
    let h = side_distance(s);
    let uf = if an >= 0.0 { upwind_value(u, c, s, 2.0 * g - u[c], h) } else { g };
    return (an * uf - diffusivity * (g - u[c]) / (0.5 * h)) * side_length(s);
}

// The residual of each cell at time t, in r.
fn fv_residual<SRC, BC>(u: f64[:], t: f64, r: f64[:]) {
    guard u.len == cell_count() && r.len == u.len else { abort("fv_residual: u and r of the cells"); }
    if flux == LINEAR_UPWIND { fv_gradients(u, t, BC); }
    if flux == LINEAR_UPWIND && limiter == VENKATAKRISHNAN { fv_venkat(u, t, BC); }
    for c in cell_count() {
        let p = centre(c);
        r[c] = -SRC(p[0], p[1], t) * hx * hy;
    }
    for c in cell_count() {
        for s in 4 {
            let nb = neighbour(c, s);
            if nb >= 0 {
                if s % 2 == 0 { continue; }
                let f = interior_flux(u, c, s, nb);
                r[c] += f;
                r[nb] -= f;
            } else if face_bc[4 * c + s] >= 0 {
                r[c] += boundary_flux(u, c, s, face_value(c, s, t, BC));
            }
        }
    }
}

// Adds v at row a and column b of val.
fn add_entry(a: i64, b: i64, v: f64) { val[solvers::entry(row_start, col, a, b)] += v; }

// The matrix of the defect correction in val: the diffusion, and the
// advection with the flux `upwind` (or `central`). Its rows are those of
// the residual. It is the exact Jacobian for `upwind` and `central`.
fn fv_matrix() {
    guard val.len == col.len && row_start.len == cell_count() + 1 else { abort("fv_matrix: the pattern of the cells"); }
    for &v in val { v = 0.0; }
    for c in cell_count() {
        for s in 4 {
            let nb = neighbour(c, s);
            let an = normal_velocity(s) * side_length(s);
            let k = diffusivity * side_length(s) / side_distance(s);
            if nb >= 0 {
                if s % 2 == 0 { continue; }
                let wc = if flux == CENTRAL { 0.5 * an } else { max(an, 0.0) };
                let wn = an - wc;
                add_entry(c, c, wc + k);
                add_entry(c, nb, wn - k);
                add_entry(nb, c, -wc - k);
                add_entry(nb, nb, -wn + k);
            } else if face_bc[4 * c + s] >= 0 {
                add_entry(c, c, (if flux == CENTRAL { 0.0 } else { max(an, 0.0) }) + 2.0 * k);
            }
        }
    }
}
