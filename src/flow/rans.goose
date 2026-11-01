// The k-ω SST model of turbulence (Menter 1994) on P1 cells: the cell
// residual of the mean flow, k and ω, and the treatment of the walls.
//
// The unknowns of a vertex are vx, vy, p, k and w (mesh.goose). With the
// kinematic viscosity ν = μ/ρ, the equations are
//
//     ρ (∇v) v - div(σ) = f,   div v = 0
//     σ = -p I + ρ (ν + ν_t) (∇v + ∇vᵀ)
//     v·∇k = P_k - β* k ω + div((ν + σ_k ν_t) ∇k)
//     v·∇ω = γ S² - β ω² + div((ν + σ_ω ν_t) ∇ω) + 2 (1 - F1) σ_ω2 ∇k·∇ω / ω
//
// with S² = 2 sym(∇v):sym(∇v), ν_t = a1 k / max(a1 ω, S F2) and
// P_k = min(ν_t S², 10 β* k ω). Each constant φ of the k and ω equations
// is F1 φ1 + (1 - F1) φ2. The residual has SUPG on each row and PSPG
// on the pressure rows (cell.goose), with the effective viscosity. The
// sources of k and ω have the vertex rule (lumped). Reason: with the
// 7-point rule, the term β ω² of the cell at a wall holds the large wall
// value of ω, and the ω row of the first vertex off the wall had no root.
// The ω rows are divided by omega_scale, so that the norm of the residual
// does not hold only the ω rows. Reason: in the channel at Re_τ = 395 the
// wall ω was 1.2e6 in scaled units, and the ω rows held 1e5 times the
// norm of the momentum rows.
//
// A wall is a boundary vertex with v = 0 fixed. The solve fixes k = 0 and
// ω = 60 ν / (β1 d1²) there, with d1 the wall distance of the nearest
// vertex off the wall (Menter 1994). The wall distance d of each vertex
// is the distance to the nearest boundary edge between two wall vertices.
import std;
import math;
namespace flow;
import .mesh;
import .cell;

let BETA_STAR = 0.09;
let A1 = 0.31;
let SIGMA_K1 = 0.85;
let SIGMA_W1 = 0.5;
let BETA1 = 0.075;
let GAMMA1 = 0.5532;
let SIGMA_K2 = 1.0;
let SIGMA_W2 = 0.856;
let BETA2 = 0.0828;
let GAMMA2 = 0.4403;
let OMEGA_FLOOR = 1.0e-8;   // scaled; the clip of ω in the coefficients

var omega_scale = 1.0;      // the ω rows are divided by it: the largest ω of a wall
var wall_d: f64[>..] = [];  // the wall distance of each vertex
var on_wall: u8[>..] = [];  // 1 for a wall vertex

// The sizes of wall_d and on_wall, before a solve.
fn clear_walls() {
    wall_d.clear();
    on_wall.clear();
    for n in vertex_count {
        wall_d.push(0.0);
        on_wall.push(0);
    }
}

// Marks the walls from the fixed velocity dofs, sets wall_d, and gives
// the value of ω at each wall vertex in omega_wall. False if the flow has
// no wall.
fn find_walls(fixed: u8[:], prescribed: f64[:], omega_wall: f64[:]) -> bool {
    guard wall_d.len == vertex_count && fixed.len == dof_count() && omega_wall.len == vertex_count else { abort("find_walls: clear_bcs on a turbulent flow"); }
    var any = false;
    for n in vertex_count {
        let wall = on_boundary[n] != 0 && fixed[vdof(n, 0)] != 0 && fixed[vdof(n, 1)] != 0
                   && prescribed[vdof(n, 0)] == 0.0 && prescribed[vdof(n, 1)] == 0.0;
        on_wall[n] = if wall { 1 } else { 0 };
        if wall { any = true; }
    }
    guard any else { return false; }
    for n in vertex_count { wall_d[n] = if on_wall[n] != 0 { 0.0 } else { 1.0e300 }; }
    for c in cell_count() {
        for e in 3 {
            let a = cell_node(c, EDGES[2 * e]);
            let b = cell_node(c, EDGES[2 * e + 1]);
            if on_wall[a] == 0 || on_wall[b] == 0 || a > b { continue; }
            for n in vertex_count {
                if on_wall[n] == 0 { wall_d[n] = min(wall_d[n], segment_distance(n, a, b)); }
            }
        }
    }
    let nu = mu / rho;
    for n in vertex_count { omega_wall[n] = 0.0; }
    for c in cell_count() {
        for a in 3 {
            let w = cell_node(c, a);
            if on_wall[w] == 0 { continue; }
            for b in 3 {
                let o = cell_node(c, b);
                if on_wall[o] != 0 { continue; }
                omega_wall[w] = max(omega_wall[w], 60.0 * nu / (BETA1 * wall_d[o] * wall_d[o]));
            }
        }
    }
    omega_scale = 1.0;
    for n in vertex_count { omega_scale = max(omega_scale, omega_wall[n]); }
    return true;
}

// The distance from vertex n to the segment from vertex a to vertex b.
fn segment_distance(n: i64, a: i64, b: i64) -> f64 {
    let ax = coords[2 * a];
    let ay = coords[2 * a + 1];
    let ux = coords[2 * b] - ax;
    let uy = coords[2 * b + 1] - ay;
    let px = coords[2 * n] - ax;
    let py = coords[2 * n + 1] - ay;
    let t = max(0.0, min(1.0, (px * ux + py * uy) / (ux * ux + uy * uy)));
    let dx = px - t * ux;
    let dy = py - t * uy;
    return sqrt(dx * dx + dy * dy);
}

// The SST values at a point.
struct Turbulence {
    nut: f64,   // ν_t
    f1: f64,
    pk: f64,    // P_k
    s2: f64,    // S²
}

fn sst(k: f64, w: f64, d: f64, s2: f64, cross: f64) -> Turbulence {
    let nu = mu / rho;
    let kp = max(k, 0.0);
    let wp = max(w, OMEGA_FLOOR);
    let dp = max(d, 1.0e-12 * diameter);
    let cd = max(2.0 * SIGMA_W2 * cross / wp, 1.0e-20);
    let root = sqrt(kp);
    let viscous = 500.0 * nu / (dp * dp * wp);
    let arg1 = min(max(root / (BETA_STAR * wp * dp), viscous), 4.0 * SIGMA_W2 * kp / (cd * dp * dp));
    let arg2 = max(2.0 * root / (BETA_STAR * wp * dp), viscous);
    let f2 = tanh_of(arg2 * arg2);
    let s = sqrt(s2);
    let nut = A1 * kp / max(A1 * wp, s * f2);
    return Turbulence { nut: nut, f1: tanh_of(arg1 * arg1 * arg1 * arg1), pk: min(nut * s2, 10.0 * BETA_STAR * kp * wp), s2: s2 };
}

// tanh(x) for x >= 0.
fn tanh_of(x: f64) -> f64 {
    let e = exp(-2.0 * x);
    return (1.0 - e) / (1.0 + e);
}

fn blend(f1: f64, a: f64, b: f64) -> f64 { return f1 * a + (1.0 - f1) * b; }

// The parameter of SUPG for a speed, a diffusivity and a reaction rate.
fn supg_tau(h: f64, speed2: f64, nu: f64, react: f64) -> f64 {
    return 1.0 / sqrt(4.0 * speed2 / (h * h) + 16.0 * nu * nu / (h * h * h * h) + react * react);
}

// The local unknowns of a turbulent cell: vx, vy of 3 vertices at 0..6,
// p at 6..9, k at 9..12 and w at 12..15.
let K0 = 9;
let W0 = 12;

// The residual r of the turbulent cell c at the local dofs u.
fn rans_residual(c: i64, u: f64[15], r: f64[15]&) {
    guard turbulent && order == 1 && c >= 0 && c < cell_count() && wall_d.len == vertex_count else { abort("rans_residual: a turbulent P1 flow after clear_bcs"); }
    let g = geometry(c);
    for d in 15 { r[d] = 0.0; }
    var gv: f64[4] = [0.0; 4];
    var gk: f64[2] = [0.0; 2];
    var gw: f64[2] = [0.0; 2];
    var gp: f64[2] = [0.0; 2];
    for a in 3 {
        for j in 2 {
            gv[j] += g.gl[2 * a + j] * u[2 * a];
            gv[2 + j] += g.gl[2 * a + j] * u[2 * a + 1];
            gp[j] += g.gl[2 * a + j] * u[6 + a];
            gk[j] += g.gl[2 * a + j] * u[K0 + a];
            gw[j] += g.gl[2 * a + j] * u[W0 + a];
        }
    }
    let s2 = 2.0 * (gv[0] * gv[0] + gv[3] * gv[3]) + (gv[1] + gv[2]) * (gv[1] + gv[2]);
    let div = gv[0] + gv[3];
    let cross = gk[0] * gw[0] + gk[1] * gw[1];
    for q in 7 {
        let l, wq = quad_point(q);
        rans_point(c, g, l, wq * g.area, u, gv, gp, gk, gw, s2, div, cross, r);
    }
    for a in 3 {
        let k = u[K0 + a];
        let om = u[W0 + a];
        let t = sst(k, om, wall_d[cell_node(c, a)], s2, cross);
        let wp = max(om, OMEGA_FLOOR);
        let m = g.area / 3.0;
        r[K0 + a] += m * (BETA_STAR * max(k, 0.0) * wp - t.pk);
        r[W0 + a] += m / omega_scale * (blend(t.f1, BETA1, BETA2) * wp * om - blend(t.f1, GAMMA1, GAMMA2) * s2 - 2.0 * (1.0 - t.f1) * SIGMA_W2 * cross / wp);
    }
}

// The terms of one quadrature point with barycentric l and weight w.
fn rans_point(c: i64, g: Geometry, l: f64[3], w: f64, u: f64[15], gv: f64[4], gp: f64[2], gk: f64[2], gw: f64[2],
              s2: f64, div: f64, cross: f64, r: f64[15]&) {
    var v: f64[2] = [0.0; 2];
    var p = 0.0;
    var k = 0.0;
    var om = 0.0;
    var d = 0.0;
    for a in 3 {
        v[0] += l[a] * u[2 * a];
        v[1] += l[a] * u[2 * a + 1];
        p += l[a] * u[6 + a];
        k += l[a] * u[K0 + a];
        om += l[a] * u[W0 + a];
        d += l[a] * wall_d[cell_node(c, a)];
    }
    let t = sst(k, om, d, s2, cross);
    let nu = mu / rho;
    let nu_eff = nu + t.nut;
    let speed2 = v[0] * v[0] + v[1] * v[1];
    let tm = supg_tau(g.h, speed2, nu_eff, 0.0);
    var res: f64[2] = [0.0; 2];
    for i in 2 { res[i] = rho * (gv[2 * i] * v[0] + gv[2 * i + 1] * v[1]) + gp[i] - force[i]; }
    let wp = max(om, OMEGA_FLOOR);
    let f1 = t.f1;
    let rk = v[0] * gk[0] + v[1] * gk[1] - t.pk + BETA_STAR * max(k, 0.0) * wp;
    let beta = blend(f1, BETA1, BETA2);
    let rw = v[0] * gw[0] + v[1] * gw[1] - blend(f1, GAMMA1, GAMMA2) * s2 + beta * wp * om - 2.0 * (1.0 - f1) * SIGMA_W2 * cross / wp;
    let nk = nu + blend(f1, SIGMA_K1, SIGMA_K2) * t.nut;
    let nw = nu + blend(f1, SIGMA_W1, SIGMA_W2) * t.nut;
    let tk = supg_tau(g.h, speed2, nk, BETA_STAR * wp);
    let tw = supg_tau(g.h, speed2, nw, 2.0 * beta * wp);
    let tc = 1.0 / (tm * g.trace_g);
    for a in 3 {
        let gx = g.gl[2 * a];
        let gy = g.gl[2 * a + 1];
        let streamline = v[0] * gx + v[1] * gy;
        for i in 2 {
            var row = l[a] * (res[i] - gp[i]) - p * g.gl[2 * a + i];
            row += rho * nu_eff * (gx * (gv[2 * i] + gv[i]) + gy * (gv[2 * i + 1] + gv[2 + i]));
            row += tm * streamline * res[i] + tc * rho * g.gl[2 * a + i] * div;
            r[2 * a + i] += w * row;
        }
        r[6 + a] += w * (l[a] * div + tm * (gx * res[0] + gy * res[1]) / rho);
        r[K0 + a] += w * (l[a] * (v[0] * gk[0] + v[1] * gk[1]) + nk * (gx * gk[0] + gy * gk[1]) + tk * streamline * rk);
        r[W0 + a] += w / omega_scale * (l[a] * (v[0] * gw[0] + v[1] * gw[1]) + nw * (gx * gw[0] + gy * gw[1]) + tw * streamline * rw);
    }
}

// The Jacobian k[15 a + b] = ∂r_a/∂u_b of the turbulent cell c at u, by
// forward differences as in cell.goose, and its residual r.
fn rans_jacobian(c: i64, u: f64[15], r: f64[15]&, k: f64[225]&) {
    rans_residual(c, u, r);
    var shifted = u;
    var rs: f64[15] = [0.0; 15];
    for b in 15 {
        let step = 1.0e-7 * max(1.0, abs(u[b]));
        shifted[b] = u[b] + step;
        rans_residual(c, shifted, rs);
        shifted[b] = u[b];
        for a in 15 { k[15 * a + b] = (rs[a] - r[a]) / step; }
    }
}
