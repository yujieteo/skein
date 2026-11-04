import std;
import math;
import solvers from "../solvers";
namespace vof;
import .mesh;
import .advect;

let UPWIND = 0;
let LINEAR_UPWIND = 1;
let CENTRAL = 2;
let NONE = 0;
let MINMOD = 1;
let VAN_LEER = 2;

let CG_RTOL = 1.0e-10;
let CG_MAX = 500;

var flux_scheme = LINEAR_UPWIND;
var limiter = VAN_LEER;

var rho_c: f64[>..] = [];
var mu_c: f64[>..] = [];
var rho_adv: f64[>..] = [];
var vstar: f64[>..] = [];
var gradv: f64[>..] = [];
var ustar: f64[>..] = [];
var ui: f64[>..] = [];
var coef: f64[>..] = [];
var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];
var rhs: f64[>..] = [];
var dp: f64[>..] = [];
var cg_r: f64[>..] = [];
var cg_z: f64[>..] = [];
var cg_p: f64[>..] = [];
var cg_q: f64[>..] = [];

var steps_done = 0;
var linear_iterations = 0;
var last_seconds = 0.0;
var max_courant = 0.0;

fn set_scheme(face_flux: i64, face_limiter: i64) {
    guard face_flux >= UPWIND && face_flux <= CENTRAL && face_limiter >= NONE && face_limiter <= VAN_LEER else { abort("set_scheme: a flux and a limiter"); }
    flux_scheme = face_flux;
    limiter = face_limiter;
}

fn pattern() {
    row_start.clear();
    col.clear();
    row_start.push(0);
    for c in cell_count() {
        for m in [neighbour(c, 2), neighbour(c, 0), c, neighbour(c, 1), neighbour(c, 3)] {
            if m >= 0 { col.push(m as i32); }
        }
        row_start.push(col.len);
    }
}

fn size_solver() {
    pattern();
    let n = cell_count();
    let faces = face_count();
    solvers::zeros(val, col.len);
    solvers::zeros(af, faces);
    solvers::zeros(ustar, faces);
    solvers::zeros(ui, faces);
    solvers::zeros(coef, faces);
    solvers::zeros(rho_adv, faces);
    solvers::zeros(cour, n);
    solvers::zeros(alpha_new, n);
    solvers::zeros(alpha_s, n);
    solvers::zeros(kappa, n);
    solvers::zeros(normal, 2 * n);
    solvers::zeros(rho_c, n);
    solvers::zeros(mu_c, n);
    solvers::zeros(vstar, 2 * n);
    solvers::zeros(gradv, 4 * n);
    solvers::zeros(rhs, n);
    solvers::zeros(dp, n);
    solvers::zeros(cg_r, n);
    solvers::zeros(cg_z, n);
    solvers::zeros(cg_p, n);
    solvers::zeros(cg_q, n);
}

fn side_value(c: i64, s: i64, i: i64) -> f64 {
    let n = neighbour(c, s);
    if n >= 0 { return vel[2 * n + i]; }
    guard side_kind[4 * c + s] >= 0 else { return vel[2 * c + i]; };
    return 2.0 * (if i == 0 { side_vx[4 * c + s] } else { side_vy[4 * c + s] }) - vel[2 * c + i];
}

fn velocity_gradients() {
    for c in cell_count() {
        for i in 2 {
            gradv[4 * c + 2 * i] = (side_value(c, 1, i) - side_value(c, 0, i)) / (2.0 * hx);
            gradv[4 * c + 2 * i + 1] = (side_value(c, 3, i) - side_value(c, 2, i)) / (2.0 * hy);
        }
    }
}

fn properties<PROPS>() {
    for c in cell_count() {
        let p = PROPS(alpha[c]);
        rho_c[c] = p[0];
        mu_c[c] = p[1];
    }
    for f in face_count() { rho_adv[f] = PROPS(af[f])[0]; }
}

fn psi(num: f64, den: f64) -> f64 {
    if limiter == NONE { return if num * den > 0.0 { 0.5 * (1.0 + num / den) } else { 0.0 }; }
    guard num * den > 0.0 else { return 0.0; };
    let r = num / den;
    return if limiter == MINMOD { min(1.0, r) } else { 2.0 * r / (1.0 + r) };
}

fn face_velocity(c: i64, s: i64, i: i64) -> f64 {
    let n = neighbour(c, s);
    if n < 0 { return if side_kind[4 * c + s] >= 0 { (side_value(c, s, i) + vel[2 * c + i]) / 2.0 } else { vel[2 * c + i] }; }
    let out = outward(s) * uf[face_of(c, s)] >= 0.0;
    let d = if out { c } else { n };
    let a = if out { n } else { c };
    let far = neighbour(d, if out { s ^ 1 } else { s });
    if flux_scheme == CENTRAL { return 0.5 * (vel[2 * d + i] + vel[2 * a + i]); }
    guard flux_scheme == LINEAR_UPWIND && far >= 0 else { return vel[2 * d + i]; };
    let den = vel[2 * a + i] - vel[2 * d + i];
    return vel[2 * d + i] + 0.5 * psi(vel[2 * d + i] - vel[2 * far + i], den) * den;
}

fn viscous_side(c: i64, s: i64, i: i64) -> f64 {
    let n = neighbour(c, s);
    let k = side_kind[4 * c + s];
    guard n >= 0 || k >= 0 else { return 0.0; };
    let a = s / 2;
    let d = side_distance(s);
    let mu = if n >= 0 { 0.5 * (mu_c[c] + mu_c[n]) } else { mu_c[c] };
    let dn = if n >= 0 { (vel[2 * n + i] - vel[2 * c + i]) / d } else { (0.5 * (side_value(c, s, i) + vel[2 * c + i]) - vel[2 * c + i]) / (0.5 * d) };
    guard i != a else { return mu * 2.0 * dn * side_length(s); };
    let tangent = if n >= 0 { 0.5 * (gradv[4 * c + 2 * a + i] + gradv[4 * n + 2 * a + i]) } else { gradv[4 * c + 2 * a + i] };
    return mu * (dn + outward(s) * tangent) * side_length(s);
}

fn predict(dt: f64) {
    velocity_gradients();
    for c in cell_count() {
        for i in 2 {
            var force = 0.0;
            for s in 4 {
                let f = face_of(c, s);
                let flux = outward(s) * uf[f] * face_area(f) * rho_adv[f];
                force += viscous_side(c, s, i) - flux * (face_velocity(c, s, i) - vel[2 * c + i]);
            }
            vstar[2 * c + i] = vel[2 * c + i] + dt * force / (rho_c[c] * cell_volume());
        }
    }
}
