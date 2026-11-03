// The `modal` step of the edge physics (spec §14.2): the resonances of
//
//     curl(ν curl(u)) + C u + B(ω) u = ω² M u
//
// with K the curl term, C the part of c that does not depend on ω, M the
// matrix of the term m u of edge_m (the coefficient of ddt2(u)), and B(ω)
// the matrix of the impedance conditions, whose z can depend on ω (a wall
// of copper has z ~ sqrt(ω)). K, C, M and B are complex symmetric, so the
// modes are M-orthogonal in the bilinear form x^T M y with no conjugate.
//
// The eigenvalue λ = ω² is complex when the walls or the material lose
// energy. For the mode k, the iteration is nonlinear in two levels:
//
//   inner  with B frozen at ω_r (a real number), the shift-invert iteration
//          x <- (A - τ M)^-1 M x, with τ = ω_r², by the direct factor of
//          src/solvers/slu.goose. A found mode is projected out of x in the
//          bilinear form, so that the next mode is the next one nearest to
//          the shift. The Rayleigh quotient λ = x^T A x / x^T M x gives λ.
//   outer  ω_r <- Re sqrt(λ), until it changes by less than 1e-10 relative.
//          B changes with ω by about B/(2 ω), thus the outer iteration
//          converges in 3 or 4 passes, and the z of a lossy wall at Re ω
//          instead of ω changes λ at the second order of the loss.
//
// The first mode is the nearest to the shift, and each next one the nearest
// to the last found, from above (a shift of 1.01 ω). The modes are sorted by
// Re ω. The dofs that an essential condition holds
// are at 0, and the null space of the curl (gradients, λ = 0) is far from a
// shift above it: the shift must be > 0.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .solve;
import .wave;
import .post;
import .edge;
import solvers from "../solvers";

let MODAL_OUTER = 60;               // the passes of the outer iteration. Reason: a pass that ends before its inner iteration converges moves the shift to the estimate (a restart of the Rayleigh quotient iteration), and a mode in a crowd of null-space modes needs several
let MODAL_INNER = 12;               // the passes of the inner iteration
let MODAL_TOLERANCE = 1.0e-12;      // the change of λ in a pass, relative
let MODAL_RESIDUAL = 1.0e-6;        // the residual of a mode, relative to |λ| |M x|
let MODAL_OMEGA = 1.0e-10;          // the change of Re ω in a pass of the outer iteration, relative

var edge_omega: f64[>..] = [];      // ω of each mode of the last edge_modal, as re and im
var edge_modes: f64[>..] = [];      // mode k at dof_total() k
var val_m: f64[>..] = [];           // M on the pattern
var val_k: f64[>..] = [];           // K + C on the pattern
var val_s: f64[>..] = [];           // K + C + B - τ M on the pattern
var modal_mass: f64[>..] = [];      // the bilinear x^T M x of each mode found, as re and im

fn csqrt(a: f64[2]) -> f64[2] {
    let r = sqrt(a[0] * a[0] + a[1] * a[1]);
    let re = sqrt(max((r + a[0]) / 2.0, 0.0));
    let im = sqrt(max((r - a[0]) / 2.0, 0.0));
    return [re, if a[1] < 0.0 { -im } else { im }];
}

// y = V x, with V on the pattern of the matrix.
fn csr_multiply(v: f64[:], x: f64[:], y: f64[:]) {
    for r in x.len {
        var s = 0.0;
        for e in row_start[r]..row_start[r + 1] { s += v[e] * x[col[e] as i64]; }
        y[r] = s;
    }
}

// The bilinear form (no conjugate) of two complex vectors with the re and
// im of dof d at 2 d and 2 d + 1.
fn bilinear(a: f64[:], b: f64[:]) -> f64[2] {
    var re = 0.0;
    var im = 0.0;
    for d in a.len / 2 {
        re += a[2 * d] * b[2 * d] - a[2 * d + 1] * b[2 * d + 1];
        im += a[2 * d] * b[2 * d + 1] + a[2 * d + 1] * b[2 * d];
    }
    return [re, im];
}

// The entries of a fixed row or column of v at 0, and with `diagonal` the
// diagonal of a fixed row at 1.
fn mask_fixed(v: f64[:], diagonal: bool) {
    for r in dof_total() {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            if fixed[r] != 0 || fixed[c] != 0 { v[e] = if diagonal && r == c { 1.0 } else { 0.0 }; }
        }
    }
}

// A start vector of the free dofs that no mode is orthogonal to.
fn start_vector(x: f64[>..]&) {
    x.clear();
    for d in dof_total() {
        let s = sin(12.9898 * (d + 1)) * 43758.5453;
        x.push(if fixed[d] != 0 { 0.0 } else { s - floor(s) - 0.5 });
    }
}

// The matrices M and K + C on the pattern, with the fixed dofs held.
fn modal_matrices() {
    edge_pattern();
    for &v in val { v = 0.0; }
    add_cells(true);
    val_m.clear();
    val_m.append(val);
    for &v in val { v = 0.0; }
    add_cells(false);
    val_k.clear();
    val_k.append(val);
    mask_fixed(val_m, false);
    mask_fixed(val_k, true);
}

// The matrix K + C + B(ω_r) - τ M into val_s.
fn modal_shifted<Z>(w: f64, tau: f64) {
    for e in val.len { val[e] = val_k[e]; }
    for i in imp_group.len {
        let z: f64[2] = Z(imp_group[i], w);
        imp_z[2 * i] = z[0];
        imp_z[2 * i + 1] = z[1];
    }
    add_impedance();
    val_s.clear();
    for e in val.len { val_s.push(val[e] - tau * val_m[e]); }
}

// One inner iteration run from x: the eigenvalue nearest to τ of the
// problem with B frozen, with the modes found projected out. It gives
// false if it does not converge in MODAL_INNER passes (λ is then the last
// estimate, or NaN if x is 0).
fn modal_inner(x: f64[:], tau: f64, lambda: f64[2]&) -> bool {
    let n = x.len;
    var u: f64[>..] = [];
    var y: f64[>..] = [];
    var sx: f64[>..] = [];
    for d in n {
        u.push(0.0);
        y.push(0.0);
        sx.push(0.0);
    }
    var previous: f64[2] = [0.0, 0.0];
    for it in MODAL_INNER {
        csr_multiply(val_m, x, u);
        for i in modal_mass.len / 2 {
            let v = edge_modes[i * n..(i + 1) * n];
            let b = cdiv(bilinear(v, u), [modal_mass[2 * i], modal_mass[2 * i + 1]]);
            for d in n / 2 {
                x[2 * d] -= b[0] * v[2 * d] - b[1] * v[2 * d + 1];
                x[2 * d + 1] -= b[0] * v[2 * d + 1] + b[1] * v[2 * d];
            }
            csr_multiply(val_m, x, u);
        }
        solvers::slu_apply(u, y);
        var top = 0.0;
        for d in n { top = max(top, abs(y[d])); }
        if !(top > 0.0) || top != top {
            lambda[0] = 0.0 / 0.0;
            return false;
        }
        for d in n { x[d] = y[d] / top; }
        csr_multiply(val_m, x, u);
        csr_multiply(val_s, x, sx);
        let q = cdiv(bilinear(x, sx), bilinear(x, u));
        lambda[0] = tau + q[0];
        lambda[1] = q[1];
        let change = sqrt((lambda[0] - previous[0]) * (lambda[0] - previous[0]) + (lambda[1] - previous[1]) * (lambda[1] - previous[1]));
        previous = [lambda[0], lambda[1]];
        if it > 0 && change <= MODAL_TOLERANCE * sqrt(lambda[0] * lambda[0] + lambda[1] * lambda[1]) {
            // The residual of S x - (λ - τ) M x, relative to |λ| |M x|.
            var r2 = 0.0;
            var m2 = 0.0;
            for d in n {
                let rr = sx[d] - (q[0] * u[d] - (if d % 2 == 0 { q[1] * u[d + 1] } else { -q[1] * u[d - 1] }));
                r2 += rr * rr;
                m2 += u[d] * u[d];
            }
            let size = sqrt(lambda[0] * lambda[0] + lambda[1] * lambda[1]) * sqrt(m2);
            if sqrt(r2) <= MODAL_RESIDUAL * size { return true; }
        }
    }
    return false;
}

// Store the mode x with the frequency om, and its bilinear mass.
fn accept_mode(x: f64[:], om: f64[2], n: i64) {
    edge_omega.push(om[0]);
    edge_omega.push(om[1]);
    edge_modes.append(x);
    var mx: f64[>..] = [];
    for d in n { mx.push(0.0); }
    csr_multiply(val_m, x, mx);
    modal_mass.append(bilinear(x, mx));
}

// The `count` resonances nearest to the shift ω_s (> 0) of the problem
// above, with the z of the conditions of the group g at ω as Z(g, ω). False
// if a mode does not converge.
fn edge_modal<Z>(count: i64, shift: f64) -> bool {
    guard physics == EDGE && dofs_sized() && edge_eps_set && count >= 1 && shift > 0.0 else {
        abort("edge_modal: an edge physics with a mass term, count >= 1 and a shift > 0");
    }
    let start = clock();
    let n = dof_total();
    for d in n { if prescribed[d] != 0.0 { abort("edge_modal: the essential conditions are 0"); } }
    modal_matrices();
    edge_omega.clear();
    edge_modes.clear();
    modal_mass.clear();
    var x: f64[>..] = [];
    for j in count {
        start_vector(x);
        // The first mode is the nearest to the shift, and each next one is the nearest to 1% above the last: the null space of the curl
        // (λ = 0) is then not nearer to the shift than the next mode, if that is below 2 times the last.
        var w = if j == 0 { shift } else { 1.01 * edge_omega[2 * j - 2] };
        var lambda: f64[2] = [w * w, 0.0];
        var done = false;
        for round in MODAL_OUTER {
            let tau = w * w;
            modal_shifted(w, tau) { g, om => Z(g, om) };
            wave_choose_path();
            guard solvers::slu_factor(row_start, col, val_s) else {
                solvers::slu_free();
                write_stderr("skein run: the factor of the modal step has a zero row\n");
                return false;
            };
            let ok = modal_inner(x, tau, lambda);
            solvers::slu_free();
            if lambda[0] != lambda[0] || !(lambda[0] > 0.0) {
                write_stderr(str("skein run: mode ", j + 1, " of the modal step has no resonance above 0 near ", w / (2.0 * PI), "\n"));
                return false;
            }
            let om = csqrt(lambda);
            let change = abs(om[0] - w) / w;
            w = om[0];
            if ok && change <= MODAL_OMEGA {
                accept_mode(x, om, n);
                done = true;
                break;
            }
        }
        if !done {
            write_stderr(str("skein run: mode ", j + 1, " of the modal step: the nonlinear iteration does not converge near ", w / (2.0 * PI), "\n"));
            return false;
        }
    }
    sort_modes(n);
    val_m.clear();
    val_k.clear();
    val_s.clear();
    clear_matrix();
    last_seconds = clock() - start;
    return true;
}

// The modes in the sequence of Re ω.
fn sort_modes(n: i64) {
    let count = edge_omega.len / 2;
    for i in count {
        var m = i;
        for j in i + 1..count { if edge_omega[2 * j] < edge_omega[2 * m] { m = j; } }
        if m == i { continue; }
        for p in 2 {
            let t = edge_omega[2 * i + p];
            edge_omega[2 * i + p] = edge_omega[2 * m + p];
            edge_omega[2 * m + p] = t;
        }
        for d in n {
            let t = edge_modes[i * n + d];
            edge_modes[i * n + d] = edge_modes[m * n + d];
            edge_modes[m * n + d] = t;
        }
        for p in 2 {
            let t = modal_mass[2 * i + p];
            modal_mass[2 * i + p] = modal_mass[2 * m + p];
            modal_mass[2 * m + p] = t;
        }
    }
}

// The frequency ω / (2 π) of mode k (from 1), as re and im.
fn edge_frequency(k: i64) -> f64[2] {
    guard k >= 1 && 2 * k <= edge_omega.len else { abort("edge_frequency: 1 <= k <= the count of the last edge_modal"); }
    return [edge_omega[2 * k - 2] / (2.0 * PI), edge_omega[2 * k - 1] / (2.0 * PI)];
}

// Writes mode k as the unknown.
fn edge_write_mode(path: const u8[:], field: const u8[:], k: i64) -> bool {
    let n = dof_total();
    guard k >= 1 && edge_modes.len >= k * n && sol.len == n else { abort("edge_write_mode: 1 <= k <= the count of the last edge_modal"); }
    for d in n { sol[d] = edge_modes[n * (k - 1) + d]; }
    return write_edge_vtu(path, field);
}
