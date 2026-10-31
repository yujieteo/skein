// The cell residual of the steady incompressible Navier-Stokes or Stokes
// equations, with the stabilization of spec §17.2, and its Jacobian.
//
// The unknowns are scaled (spec §21.2): the generated program gives rho,
// mu and the body force f in the scaled units of the scales line, so the
// values of the unknowns are near 1.
//
// For the test functions w = N_a e_i and q = L_k, the residual of a cell is
//
//     r(a i) = ∫ N_a (ρ (∇v) v - f)_i + ∂_j N_a σ_ij
//              + τ_M (v·∇N_a) R_i                     (supg)
//              + τ_C ρ ∂_i N_a div v                  (lsic)
//     r(k)   = ∫ L_k div v + τ_M ∂_i L_k R_i / ρ      (pspg)
//
// with σ = -p I + μ (∇v + ∇vᵀ) and the strong residual of the momentum
// equation R = ρ (∇v) v + ∇p - div(μ (∇v + ∇vᵀ)) - f. For Stokes the
// convection term (∇v) v is 0. The stabilization parameters are
// τ_M = ((2|v|/h)² + (4ν/h²)²)^(-1/2) and τ_C = 1 / (τ_M tr(G)), after
// Tezduyar (1992), with h = 2 sqrt(A/π) / order and the metric G of the
// cell. A static step has no 2/Δt term.
//
// The Jacobian of a cell is the forward difference of its residual, one
// dof at a time. Reason: τ_M depends on v, and the difference includes
// that derivative, so Newton's method keeps its quadratic rate.
import std;
import math;
namespace flow;
import .mesh;

let NONE = 0;
let SUPG = 1;
let SUPG_PSPG = 2;
let SUPG_PSPG_LSIC = 3;

var rho = 1.0;
var mu = 1.0;
var force: f64[2] = [0.0, 0.0];
var convective = true;      // Navier-Stokes; false for Stokes
var stabilize = NONE;

// The 7-point rule of degree 5 on a triangle (Radon 1948): barycentric
// coordinates and weights that sum to 1.
let QA1 = 0.059715871789770;
let QB1 = 0.470142064105115;
let QA2 = 0.797426985353087;
let QB2 = 0.101286507323456;
let QW0 = 0.225;
let QW1 = 0.132394152788506;
let QW2 = 0.125939180544827;

fn quad_point(q: i64) -> f64[3], f64 {
    guard q >= 0 && q < 7 else { abort("quad_point: 0 <= q < 7"); }
    if q == 0 { return [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], QW0; }
    let k = (q - 1) % 3;
    let a = if q <= 3 { QA1 } else { QA2 };
    let b = if q <= 3 { QB1 } else { QB2 };
    var l: f64[3] = [b, b, b];
    l[k] = a;
    return l, if q <= 3 { QW1 } else { QW2 };
}

// The geometry of a cell: the gradients of the barycentric coordinates,
// the area and the size h.
struct Geometry {
    gl: f64[6],     // ∂L_k/∂x, ∂L_k/∂y at 2 k, 2 k + 1
    area: f64,
    h: f64,
    trace_g: f64,   // tr(G), G = (∂ξ/∂x)ᵀ(∂ξ/∂x) with ξ = (L1, L2)
}

fn geometry(c: i64) -> Geometry {
    guard c >= 0 && c < cell_count() else { abort("geometry: c < cell_count()"); }
    var x: f64[3] = [0.0; 3];
    var y: f64[3] = [0.0; 3];
    for k in 3 {
        x[k] = coords[2 * cell_node(c, k)];
        y[k] = coords[2 * cell_node(c, k) + 1];
    }
    let det = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
    var g = Geometry { gl: [0.0; 6], area: 0.5 * det, h: 0.0, trace_g: 0.0 };
    for k in 3 {
        let i = (k + 1) % 3;
        let j = (k + 2) % 3;
        g.gl[2 * k] = (y[i] - y[j]) / det;
        g.gl[2 * k + 1] = (x[j] - x[i]) / det;
    }
    g.h = 2.0 * sqrt(g.area / PI) / order;
    for k in 1..3 { g.trace_g += g.gl[2 * k] * g.gl[2 * k] + g.gl[2 * k + 1] * g.gl[2 * k + 1]; }
    return g;
}

// The shape functions of the velocity at the barycentric point l: values
// n, gradients dn (2 per node) and the Laplacians lap and the Hessians hs
// (xx, xy, yy per node). For P1 the second derivatives are 0.
struct Shape {
    n: f64[6],
    dn: f64[12],
    hs: f64[18],
}

fn shape(g: Geometry, l: f64[3]) -> Shape {
    var s = Shape { n: [0.0; 6], dn: [0.0; 12], hs: [0.0; 18] };
    for k in 3 {
        let gx = g.gl[2 * k];
        let gy = g.gl[2 * k + 1];
        if order == 1 {
            s.n[k] = l[k];
            s.dn[2 * k] = gx;
            s.dn[2 * k + 1] = gy;
            continue;
        }
        s.n[k] = l[k] * (2.0 * l[k] - 1.0);
        s.dn[2 * k] = (4.0 * l[k] - 1.0) * gx;
        s.dn[2 * k + 1] = (4.0 * l[k] - 1.0) * gy;
        s.hs[3 * k] = 4.0 * gx * gx;
        s.hs[3 * k + 1] = 4.0 * gx * gy;
        s.hs[3 * k + 2] = 4.0 * gy * gy;
    }
    guard order == 2 else { return s; }
    for e in 3 {
        let a = EDGES[2 * e];
        let b = EDGES[2 * e + 1];
        let m = 3 + e;
        s.n[m] = 4.0 * l[a] * l[b];
        for t in 2 { s.dn[2 * m + t] = 4.0 * (l[b] * g.gl[2 * a + t] + l[a] * g.gl[2 * b + t]); }
        s.hs[3 * m] = 8.0 * g.gl[2 * a] * g.gl[2 * b];
        s.hs[3 * m + 1] = 4.0 * (g.gl[2 * a] * g.gl[2 * b + 1] + g.gl[2 * a + 1] * g.gl[2 * b]);
        s.hs[3 * m + 2] = 8.0 * g.gl[2 * a + 1] * g.gl[2 * b + 1];
    }
    return s;
}

// The flow at a point: v, ∇v (∂_j v_i at 2 i + j), p, ∇p and the strong
// residual R of the momentum equation.
struct Point {
    v: f64[2],
    gv: f64[4],
    p: f64,
    gp: f64[2],
    r: f64[2],
}

fn flow_point(g: Geometry, s: Shape, l: f64[3], u: f64[15]) -> Point {
    var f = Point { v: [0.0; 2], gv: [0.0; 4], p: 0.0, gp: [0.0; 2], r: [0.0; 2] };
    var visc: f64[2] = [0.0; 2];
    for a in per_cell {
        for i in 2 {
            let ui = u[2 * a + i];
            f.v[i] += s.n[a] * ui;
            for j in 2 { f.gv[2 * i + j] += s.dn[2 * a + j] * ui; }
        }
        // div(∇v + ∇vᵀ)_i = Δv_i + ∂_i div v
        let hxx = s.hs[3 * a];
        let hxy = s.hs[3 * a + 1];
        let hyy = s.hs[3 * a + 2];
        visc[0] += (hxx + hyy) * u[2 * a] + hxx * u[2 * a] + hxy * u[2 * a + 1];
        visc[1] += (hxx + hyy) * u[2 * a + 1] + hxy * u[2 * a] + hyy * u[2 * a + 1];
    }
    for k in 3 {
        let pk = u[2 * per_cell + k];
        f.p += l[k] * pk;
        f.gp[0] += g.gl[2 * k] * pk;
        f.gp[1] += g.gl[2 * k + 1] * pk;
    }
    for i in 2 {
        let conv = if convective { rho * (f.gv[2 * i] * f.v[0] + f.gv[2 * i + 1] * f.v[1]) } else { 0.0 };
        f.r[i] = conv + f.gp[i] - mu * visc[i] - force[i];
    }
    return f;
}

// τ_M and τ_C at a point with velocity v.
fn taus(g: Geometry, v: f64[2]) -> f64, f64 {
    let speed2 = v[0] * v[0] + v[1] * v[1];
    let nu = mu / rho;
    let advect = if convective { 4.0 * speed2 / (g.h * g.h) } else { 0.0 };
    let diffuse = 16.0 * nu * nu / (g.h * g.h * g.h * g.h);
    let tm = 1.0 / sqrt(advect + diffuse);
    return tm, 1.0 / (tm * g.trace_g);
}

// The residual r of cell c at the local dofs u (cell_dofs).
fn cell_residual(c: i64, u: f64[15], r: f64[15]&) {
    guard c >= 0 && c < cell_count() && rho > 0.0 && mu > 0.0 else { abort("cell_residual: c < cell_count(), rho > 0 and mu > 0"); }
    let g = geometry(c);
    for d in 15 { r[d] = 0.0; }
    for q in 7 {
        let l, wq = quad_point(q);
        let w = wq * g.area;
        let s = shape(g, l);
        let f = flow_point(g, s, l, u);
        let tm, tc = taus(g, f.v);
        let div = f.gv[0] + f.gv[3];
        var sigma: f64[4] = [0.0; 4];
        for i in 2 {
            for j in 2 { sigma[2 * i + j] = mu * (f.gv[2 * i + j] + f.gv[2 * j + i]) - (if i == j { f.p } else { 0.0 }); }
        }
        momentum_rows(s, f, sigma, tm, tc, div, w, r);
        for k in 3 {
            var row = l[k] * div;
            if stabilize >= SUPG_PSPG { row += tm * (g.gl[2 * k] * f.r[0] + g.gl[2 * k + 1] * f.r[1]) / rho; }
            r[2 * per_cell + k] += w * row;
        }
    }
}

fn momentum_rows(s: Shape, f: Point, sigma: f64[4], tm: f64, tc: f64, div: f64, w: f64, r: f64[15]&) {
    for a in per_cell {
        let streamline = f.v[0] * s.dn[2 * a] + f.v[1] * s.dn[2 * a + 1];
        for i in 2 {
            let conv = if convective { rho * (f.gv[2 * i] * f.v[0] + f.gv[2 * i + 1] * f.v[1]) } else { 0.0 };
            var row = s.n[a] * (conv - force[i]) + s.dn[2 * a] * sigma[2 * i] + s.dn[2 * a + 1] * sigma[2 * i + 1];
            if stabilize >= SUPG { row += tm * streamline * f.r[i]; }
            if stabilize == SUPG_PSPG_LSIC { row += tc * rho * s.dn[2 * a + i] * div; }
            r[2 * a + i] += w * row;
        }
    }
}

// The Jacobian k[15 a + b] = ∂r_a/∂u_b of cell c at u, and its residual r.
fn cell_jacobian(c: i64, u: f64[15], r: f64[15]&, k: f64[225]&) {
    guard c >= 0 && c < cell_count() else { abort("cell_jacobian: c < cell_count()"); }
    cell_residual(c, u, r);
    let n = cell_dof_count();
    var shifted = u;
    var rs: f64[15] = [0.0; 15];
    for b in n {
        let step = 1.0e-7 * max(1.0, abs(u[b]));
        shifted[b] = u[b] + step;
        cell_residual(c, shifted, rs);
        shifted[b] = u[b];
        for a in n { k[15 * a + b] = (rs[a] - r[a]) / step; }
    }
}

// The PSPG term of the pressure rows, τ_M/ρ ∫ ∇L_k·∇L_l, of cell c at u,
// with τ_M at the centroid: k[3 k + l]. The preconditioner adds it when
// the equations have no PSPG (solve.goose).
fn pressure_term(c: i64, u: f64[15], k: f64[9]&) {
    guard c >= 0 && c < cell_count() else { abort("pressure_term: c < cell_count()"); }
    let g = geometry(c);
    let l: f64[3] = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    let f = flow_point(g, shape(g, l), l, u);
    let tm, tc = taus(g, f.v);
    for a in 3 {
        for b in 3 { k[3 * a + b] = tm / rho * g.area * (g.gl[2 * a] * g.gl[2 * b] + g.gl[2 * a + 1] * g.gl[2 * b + 1]); }
    }
}
