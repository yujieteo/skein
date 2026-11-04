import std;
import math;
import os;
namespace cons;
import solvers from "../solvers";
import .mesh;

let CFL_ABORT = 2.0;
let STABLE_FRACTION = 0.75;
let ADVISED_FRACTION = 0.8;
let RADIUS_SAFETY = 1.02;
let POWER_STEPS = 8;
let FD_RELATIVE = 1.0e-6;

var grad: f64[>..] = [];
var jac_x: f64[>..] = [];
var jac_y: f64[>..] = [];
var radius_side: f64[>..] = [];
var comp_scale: f64[>..] = [];
var r: f64[>..] = [];
var u1: f64[>..] = [];
var u2: f64[>..] = [];
var ul: f64[>..] = [];
var ur: f64[>..] = [];
var fa: f64[>..] = [];
var fb: f64[>..] = [];
var ff: f64[>..] = [];
var pt: f64[>..] = [];
var ug: f64[>..] = [];
var pw: f64[>..] = [];
var pz: f64[>..] = [];
var jm: f64[>..] = [];

var iterations = 0;
var last_seconds = 0.0;
var max_rate = 0.0;

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

fn psi(q: f64) -> f64 {
    if limiter == MINMOD { return max(0.0, min(1.0, q)); }
    if limiter == VAN_LEER { return (q + abs(q)) / (1.0 + abs(q)); }
    return 0.5 * (1.0 + q);
}

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

fn face_value(u: f64[:], c: i64, d: i64, k: i64) -> f64 {
    let n = comps;
    let uc = u[c * n + k];
    let du = u[d * n + k] - uc;
    let g = grad[(c * n + k) * 2] * (centre_x[d] - centre_x[c]) + grad[(c * n + k) * 2 + 1] * (centre_y[d] - centre_y[c]);
    if limiter == NONE { return uc + 0.5 * g; }
    if abs(du) <= 1.0e-14 * (abs(uc) + abs(u[d * n + k])) { return uc; }
    return uc + 0.5 * psi(2.0 * g / du - 1.0) * du;
}

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

fn rusanov<FLUX>(nx: f64, ny: f64, lam: f64, out: f64[:]) {
    let n = comps;
    FLUX(ul, nx, ny, fa);
    FLUX(ur, nx, ny, fb);
    for k in n { out[k] = 0.5 * (fa[k] + fb[k]) - 0.5 * lam * (ur[k] - ul[k]); }
}

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

fn stage<FLUX, SRC, GHOST, WALL>(u0: f64[:], v: f64[:], t: f64, dt: f64, c0: f64, c1: f64, refresh: bool, out: f64[:]) {
    let n = comps;
    residual(v, t, refresh, FLUX, SRC, GHOST, WALL);
    for c in cell_total {
        let w = dt / area[c];
        for k in n { out[c * n + k] = c0 * u0[c * n + k] + c1 * (v[c * n + k] - w * r[c * n + k]); }
    }
}

fn state_sizes() -> i64[12] {
    return [sol.len, grad.len, jac_x.len, jac_y.len, radius_side.len, r.len, u1.len, u2.len, ul.len, ur.len, ff.len, jm.len];
}

fn finite() -> bool {
    for v in sol { if !(abs(v) < 1.0e300) { return false; } }
    return true;
}

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
                             " s of the mesh and the waves of the initial state. Write a dt of ", ADVISED_FRACTION * CFL_ABORT / max_rate, " s or less.\n"));
            return false;
        }
        stage(sol, u1, t + k, k, 0.75, 0.25, false, u2, FLUX, SRC, GHOST, WALL);
        stage(sol, u2, t + 0.5 * k, k, 1.0 / 3.0, 2.0 / 3.0, false, sol, FLUX, SRC, GHOST, WALL);
        i++;
        if i == 10 { write_stderr(str("skein run: ", steps, " steps of ", cell_total, " cells take about ", (clock() - start) / 10.0 * steps, " s\n")); }
        if k * max_rate > CFL_ABORT && i % 100 == 0 {
            write_stderr(str("skein run: at the time ", t, " s the step ", k, " s is above the stable step ", CFL_ABORT / max_rate, " s. Write a dt of ", ADVISED_FRACTION * CFL_ABORT / max_rate, " s or less.\n"));
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
