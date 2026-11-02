// The cell-centred finite volume method for a system of conservation laws,
//
//     ddt(U) + div(F(U)) = S(U, x, y, t)
//
// on triangles, with the explicit SSP Runge-Kutta scheme of order 3 (Shu
// and Osher 1988). U holds the N components of the unknowns of all
// equations, and F(U)·n is a function that the generated program gives
// (`flux`), written from the flux of each `div` of the model. Nothing in
// this file is the Euler equations: a model with any flux that is linear
// in n runs here.
//
// The value of U on a side comes from the cell by a limited linear
// reconstruction (MUSCL, van Leer 1979): the gradient of each component is
// the Green-Gauss gradient, and the face value is
// U_C + ψ(r) (U_D - U_C)/2 with r = 2 ∇U_C·d_CD/(U_D - U_C) - 1 and ψ the
// limiter (the same as the fv method of src/grid/fv.goose). The numerical
// flux is the Rusanov flux
//
//     F̂ = ((F(U_L) + F(U_R))·n - λ (U_R - U_L)) / 2
//
// with λ the larger spectral radius of dF·n/dU of the two cells. The
// program gives F only, so the Jacobian of each cell is a forward
// difference of F, and its spectral radius is the power iteration on its
// square. Reason: a model of any physics can use it, and a table of the
// wave speeds of each system would run the Euler equations only. The
// iteration gives the radius from below, thus λ takes the factor 1.02.
//
// A side on the boundary has the mask of its conditions (mesh.goose). With
// none, F̂ is F(U_C)·n (the flow leaves as it is). A condition on an unknown
// gives the ghost state U_G with that component set, and F̂ is the Rusanov
// flux of (U_C, U_G); a condition on a flux sets the component of F̂.
import std;
import math;
import os;
namespace cons;
import solvers from "../solvers";
import .mesh;

let CFL_ABORT = 2.0;            // dt × the rate of a cell: the Rusanov step is unstable above 2
let STABLE_FRACTION = 0.75;     // dt: stable takes this fraction of the step above which the Rusanov step is unstable
let RADIUS_SAFETY = 1.02;
let POWER_STEPS = 8;
let FD_RELATIVE = 1.0e-6;

var grad: f64[>..] = [];        // 2 N for each cell
var jac_x: f64[>..] = [];       // N × N for each cell: dF·(1, 0)/dU
var jac_y: f64[>..] = [];
var radius_side: f64[>..] = []; // the spectral radius at side k of cell c, 3 for each cell
var comp_scale: f64[>..] = [];
var r: f64[>..] = [];
var u1: f64[>..] = [];
var u2: f64[>..] = [];
var ul: f64[>..] = [];
var ur: f64[>..] = [];
var fa: f64[>..] = [];
var fb: f64[>..] = [];
var ff: f64[>..] = [];          // the numerical flux of a side
var pt: f64[>..] = [];
var ug: f64[>..] = [];
var pw: f64[>..] = [];
var pz: f64[>..] = [];
var jm: f64[>..] = [];

var iterations = 0;
var last_seconds = 0.0;
var max_rate = 0.0;             // the largest of the sum of the lengths times the radius over the area

fn size_solver() {
    let n = comps;
    let m = cell_total;
    solvers::zeros(grad, 2 * n * m);
    solvers::zeros(jac_x, n * n * m);
    solvers::zeros(jac_y, n * n * m);
    solvers::zeros(radius_side, 3 * m);
    solvers::zeros(comp_scale, n);
    solvers::zeros(r, n * m);
    solvers::zeros(u1, n * m);
    solvers::zeros(u2, n * m);
    solvers::zeros(ul, n);
    solvers::zeros(ur, n);
    solvers::zeros(fa, n);
    solvers::zeros(fb, n);
    solvers::zeros(ff, n);
    solvers::zeros(pt, n);
    solvers::zeros(ug, n);
    solvers::zeros(pw, n);
    solvers::zeros(pz, n);
    solvers::zeros(jm, n * n);
}

// ψ(r) of the limiter.
fn psi(q: f64) -> f64 {
    if limiter == MINMOD { return max(0.0, min(1.0, q)); }
    if limiter == VAN_LEER { return (q + abs(q)) / (1.0 + abs(q)); }
    return 0.5 * (1.0 + q);
}

// The Green-Gauss gradient of each component in each cell, in grad. A side
// on the boundary has the value of the cell.
fn gradients(u: f64[:]) {
    let n = comps;
    for c in cell_total {
        for k in n { grad[(c * n + k) * 2] = 0.0; grad[(c * n + k) * 2 + 1] = 0.0; }
        for s in 3 {
            let nb = neighbour[3 * c + s];
            for k in n {
                let um = if nb >= 0 { 0.5 * (u[c * n + k] + u[nb * n + k]) } else { u[c * n + k] };
                grad[(c * n + k) * 2] += um * side_nx[3 * c + s];
                grad[(c * n + k) * 2 + 1] += um * side_ny[3 * c + s];
            }
        }
        for k in n * 2 { grad[c * n * 2 + k] /= area[c]; }
    }
}

// The value of component k of the cell c on the side to the cell d.
fn face_value(u: f64[:], c: i64, d: i64, k: i64) -> f64 {
    let n = comps;
    let uc = u[c * n + k];
    let du = u[d * n + k] - uc;
    let g = grad[(c * n + k) * 2] * (centre_x[d] - centre_x[c]) + grad[(c * n + k) * 2 + 1] * (centre_y[d] - centre_y[c]);
    if limiter == NONE { return uc + 0.5 * g; }
    if abs(du) <= 1.0e-14 * (abs(uc) + abs(u[d * n + k])) { return uc; }
    return uc + 0.5 * psi(2.0 * g / du - 1.0) * du;
}

// The Jacobian of F·(1, 0) and F·(0, 1) of each cell, by forward differences.
fn jacobians<FLUX>(u: f64[:]) {
    let n = comps;
    for k in n {
        var big = 0.0;
        for c in cell_total { big = max(big, abs(u[c * n + k])); }
        comp_scale[k] = big;
    }
    for c in cell_total {
        for k in n { ul[k] = u[c * n + k]; }
        for dir in 2 {
            FLUX(ul, if dir == 0 { 1.0 } else { 0.0 }, if dir == 0 { 0.0 } else { 1.0 }, fa);
            for j in n {
                let h = FD_RELATIVE * (abs(ul[j]) + 1.0e-3 * comp_scale[j]) + 1.0e-300;
                let saved = ul[j];
                ul[j] = saved + h;
                FLUX(ul, if dir == 0 { 1.0 } else { 0.0 }, if dir == 0 { 0.0 } else { 1.0 }, fb);
                ul[j] = saved;
                for i in n {
                    let v = (fb[i] - fa[i]) / h;
                    if dir == 0 { jac_x[(c * n + i) * n + j] = v; } else { jac_y[(c * n + i) * n + j] = v; }
                }
            }
        }
    }
}

// The spectral radius of nx Jx + ny Jy of cell c, by the power iteration on
// its square.
fn spectral(c: i64, nx: f64, ny: f64) -> f64 {
    let n = comps;
    for i in n * n { jm[i] = nx * jac_x[c * n * n + i] + ny * jac_y[c * n * n + i]; }
    for i in n { pw[i] = 1.0 + 0.1 * i; }
    var lam2 = 0.0;
    for it in POWER_STEPS {
        var nv = 0.0;
        for i in n { nv += pw[i] * pw[i]; }
        for i in n {
            var s = 0.0;
            for j in n { s += jm[i * n + j] * pw[j]; }
            pz[i] = s;
        }
        for i in n {
            var s = 0.0;
            for j in n { s += jm[i * n + j] * pz[j]; }
            pt[i] = s;
        }
        var nz = 0.0;
        for i in n { nz += pt[i] * pt[i]; }
        if nz <= 0.0 || nv <= 0.0 { return 0.0; }
        lam2 = sqrt(nz / nv);
        let scale = 1.0 / sqrt(nz);
        for i in n { pw[i] = pt[i] * scale; }
    }
    return RADIUS_SAFETY * sqrt(lam2);
}

// The radius at each side of each cell, and the largest rate.
fn radii() {
    max_rate = 0.0;
    for c in cell_total {
        var rate = 0.0;
        for s in 3 {
            let lx = side_nx[3 * c + s];
            let ly = side_ny[3 * c + s];
            let l = sqrt(lx * lx + ly * ly);
            let q = spectral(c, lx / l, ly / l);
            radius_side[3 * c + s] = q;
            rate += q * l;
        }
        max_rate = max(max_rate, rate / area[c]);
    }
}

// The Rusanov flux of the states ul and ur on a side with the unit normal
// (nx, ny) and the radius lam, per unit of length, in out.
fn rusanov<FLUX>(nx: f64, ny: f64, lam: f64, out: f64[:]) {
    let n = comps;
    FLUX(ul, nx, ny, fa);
    FLUX(ur, nx, ny, fb);
    for k in n { out[k] = 0.5 * (fa[k] + fb[k]) - 0.5 * lam * (ur[k] - ul[k]); }
}

// The flux through side s of cell c, times its length, added to the
// residual of the cell and taken from the residual of the neighbour.
fn side_flux<FLUX, GHOST, WALL>(u: f64[:], c: i64, s: i64, t: f64) {
    let n = comps;
    let lx = side_nx[3 * c + s];
    let ly = side_ny[3 * c + s];
    let l = sqrt(lx * lx + ly * ly);
    let nx = lx / l;
    let ny = ly / l;
    let nb = neighbour[3 * c + s];
    if nb >= 0 {
        if nb < c { return; }
        for k in n {
            ul[k] = face_value(u, c, nb, k);
            ur[k] = face_value(u, nb, c, k);
        }
        let opposite = if neighbour[3 * nb] == c { 0 } else if neighbour[3 * nb + 1] == c { 1 } else { 2 };
        let lam = max(radius_side[3 * c + s], radius_side[3 * nb + opposite]);
        rusanov(nx, ny, lam, ff, FLUX);
        for k in n {
            r[c * n + k] += ff[k] * l;
            r[nb * n + k] -= ff[k] * l;
        }
        return;
    }
    for k in n { ul[k] = u[c * n + k]; }
    let mask = side_mask[3 * c + s];
    if mask == FREE {
        FLUX(ul, nx, ny, ff);
    } else {
        GHOST(mask, ul, side_mx[3 * c + s], side_my[3 * c + s], t, nx, ny, ug);
        for k in n { ur[k] = ug[k]; }
        rusanov(nx, ny, radius_side[3 * c + s], ff, FLUX);
        WALL(mask, ul, side_mx[3 * c + s], side_my[3 * c + s], t, nx, ny, ff);
    }
    for k in n { r[c * n + k] += ff[k] * l; }
}

// The residual of each cell at time t, in r: the flux out less the source
// times the area.
fn residual<FLUX, SRC, GHOST, WALL>(u: f64[:], t: f64, refresh: bool) {
    let n = comps;
    if refresh {
        jacobians(u, FLUX);
        radii();
    }
    gradients(u);
    for c in cell_total {
        for k in n { ug[k] = 0.0; }
        for k in n { ul[k] = u[c * n + k]; }
        SRC(ul, centre_x[c], centre_y[c], t, ug);
        for k in n { r[c * n + k] = -ug[k] * area[c]; }
    }
    for c in cell_total {
        for s in 3 { side_flux(u, c, s, t, FLUX, GHOST, WALL); }
    }
}

// out = c0 u0 + c1 (v - dt R(v, t)/A), the stage of SSP-RK3.
fn stage<FLUX, SRC, GHOST, WALL>(u0: f64[:], v: f64[:], t: f64, dt: f64, c0: f64, c1: f64, refresh: bool, out: f64[:]) {
    let n = comps;
    residual(v, t, refresh, FLUX, SRC, GHOST, WALL);
    for c in cell_total {
        let w = dt / area[c];
        for k in n { out[c * n + k] = c0 * u0[c * n + k] + c1 * (v[c * n + k] - w * r[c * n + k]); }
    }
}

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[12] {
    return [sol.len, grad.len, jac_x.len, jac_y.len, radius_side.len, r.len, u1.len, u2.len, ul.len, ur.len, ff.len, jm.len];
}

// Whether each value of the state is finite.
fn finite() -> bool {
    for v in sol { if !(abs(v) < 1.0e300) { return false; } }
    return true;
}

// Integrates from the time of the state to `end` with the steps of dt. It
// stops with a message, and gives false, if dt is above the step that the
// radius allows (dt <= 0 takes the stable step), or if the state is not finite. It tells the time that the
// whole run takes after the first steps.
fn solve_transient<FLUX, SRC, GHOST, WALL>(end: f64, dt: f64) -> bool {
    guard comps >= 1 && end > time && sol.len == dof_count() else { abort("solve_transient: setup and an end after the time of the state"); }
    let sizes = state_sizes();
    solvers::memory_check("cons ssprk3", 0, dof_count(), 6, 0, "");
    var h = dt;
    if h <= 0.0 {
        jacobians(sol, FLUX);
        radii();
        h = STABLE_FRACTION * CFL_ABORT / max(max_rate, 1.0e-300);
    }
    let steps = max(1, ceil((end - time) / h - 1.0e-9) as! i64);
    let k = (end - time) / steps;
    let t0 = time;
    let start = clock();
    var good = true;
    var i = 0;
    while good && i < steps {
        let t = t0 + i * k;
        stage(sol, sol, t, k, 0.0, 1.0, true, u1, FLUX, SRC, GHOST, WALL);
        if i == 0 && k * max_rate > CFL_ABORT {
            write_stderr(str("skein run: the step ", k, " s is above the stable step ", CFL_ABORT / max_rate,
                             " s of the mesh and the waves of the initial state. Write a dt of ", 0.8 * CFL_ABORT / max_rate, " s or less.\n"));
            return false;
        }
        stage(sol, u1, t + k, k, 0.75, 0.25, false, u2, FLUX, SRC, GHOST, WALL);
        stage(sol, u2, t + 0.5 * k, k, 1.0 / 3.0, 2.0 / 3.0, false, sol, FLUX, SRC, GHOST, WALL);
        i++;
        if i == 10 { write_stderr(str("skein run: ", steps, " steps of ", cell_total, " cells take about ", (clock() - start) / 10.0 * steps, " s\n")); }
        if k * max_rate > CFL_ABORT && i % 100 == 0 {
            write_stderr(str("skein run: at the time ", t, " s the step ", k, " s is above the stable step ", CFL_ABORT / max_rate, " s. Write a dt of ", 0.8 * CFL_ABORT / max_rate, " s or less.\n"));
            good = false;
        }
        if i % 100 == 0 && !finite() {
            write_stderr(str("skein run: the state is not finite at the time ", t + k, " s. Write a smaller dt.\n"));
            good = false;
        }
    }
    time = t0 + i * k;
    iterations = i;
    last_seconds = clock() - start;
    guard state_sizes() == sizes else { abort("solve_transient: a solve grows no array of the state (rule 3)"); }
    return good && finite();
}
