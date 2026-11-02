// The harmonic step of a real physics (spec §14.2): the steady response of
// K u + ddt2 terms to a load of the form F exp(i ω t), with u a `Complex`
// vector. This is the complex lift of the real assembly of the physics:
// the matrix of the frequency is
//
//     A = (1 + 2 i ζ) K - ω² M,
//
// with K the stiffness (assemble.goose) and M the mass (eigen.goose). The
// load is real. ζ is the damping ratio of each mode (the structural damping
// of spec §14.2, as in wave.goose).
//
// K and M are symmetric, so A is complex symmetric. A skyline LDLᵀ with
// complex entries (no conjugate, no pivoting) solves it, in the order of
// direct.goose. A is indefinite for ω above the first mode, and the factor
// of ω near a resonance has a small pivot: a pivot below 1e-12 of its
// diagonal is a failure, as in direct.goose.
//
// Reason: the random vibration model of femahoy needs the response at 130
// frequencies of a Complex vector u, and the modal sum of a LOBPCG run
// would be a truncated answer to a question that the solver can answer.
//
// The real part of the response is in `sol` after the solve; harmonic_part
// puts the imaginary part (1) or the real part (0) in `sol`, so that the
// functions of post.goose (value_at, region_mean, write_vtu) read either.
import std;
import math;
namespace elements;
import .mesh;
import .assemble;
import .eigen;
import .direct;
import .solve;
import .post;
import solvers from "../solvers";

let HARMONIC_WORK = 4.0e9;    // multiplies of one factorization: about 8 s

var harmonic_damping = 0.0;
var sol_re: f64[>..] = [];
var sol_im: f64[>..] = [];
var up_re: f64[>..] = [];     // the upper part of the profile, column by column as sky_upper
var up_im: f64[>..] = [];
var lo_re: f64[>..] = [];     // L, as sky_lower
var lo_im: f64[>..] = [];

// The next solve of the elastic physics is the harmonic one, at ω and ζ.
fn set_harmonic(omega: f64, damping: f64) {
    guard physics == ELASTIC && omega > 0.0 && damping >= 0.0 && density > 0.0 else {
        abort("set_harmonic: an elastic physics with a density, omega > 0 and damping >= 0");
    }
    harmonic_omega = omega;
    harmonic_damping = damping;
}

// Puts part p (0 real, 1 imaginary) of the response in sol.
fn harmonic_part(p: i64) {
    guard (p == 0 || p == 1) && sol_re.len == sol.len && sol_im.len == sol.len else { abort("harmonic_part: part 0 or 1 of a solved harmonic step"); }
    for d in sol.len { sol[d] = if p == 0 { sol_re[d] } else { sol_im[d] }; }
}

// Part p of component comp at a point, and whether a cell holds it.
fn harmonic_at(p: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    harmonic_part(p);
    return field_at(FIELD_VALUE, comp, x, y, z);
}

// Part p of the mean of component comp on the region of the last call that
// selects one.
fn harmonic_mean(p: i64, comp: i64) -> f64 {
    harmonic_part(p);
    return region_mean(comp);
}

// The sum of L[ai + k] U[bj + k], for k from lo to hi - 1, as re and im.
fn sky_cdot(lr: f64[:], li: f64[:], a0: i64, ur: f64[:], ui: f64[:], b0: i64, lo: i64, hi: i64) -> f64, f64 {
    var re = 0.0;
    var im = 0.0;
    for k in lo..hi {
        re += lr[a0 + k] * ur[b0 + k] - li[a0 + k] * ui[b0 + k];
        im += lr[a0 + k] * ui[b0 + k] + li[a0 + k] * ur[b0 + k];
    }
    return re, im;
}

// The quotient a / b of two complex numbers.
fn quotient(ar: f64, ai: f64, br: f64, bi: f64) -> f64, f64 {
    let d = br * br + bi * bi;
    return (ar * br + ai * bi) / d, (ai * br - ar * bi) / d;
}

// A = L D Lᵀ in place, with the fixed dofs as rows of the identity. False
// if a pivot is near 0 relative to its diagonal.
fn harmonic_factor(omega: f64, damping: f64) -> bool {
    for &v in up_re { v = 0.0; }
    for &v in up_im { v = 0.0; }
    let w2 = omega * omega;
    for r in sol.len {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            guard c == r || fixed[r] == 0 && fixed[c] == 0 else { continue; };
            let i = sky_row(r);
            let j = sky_row(c);
            guard i <= j else { continue; };
            let at = sky_start[j] + i - sky_top[j];
            if fixed[r] != 0 {
                up_re[at] = 1.0;
            } else {
                up_re[at] = val[e] - w2 * val_b[e];
                up_im[at] = 2.0 * damping * val[e];
            }
        }
    }
    for j in sol.len {
        let uj = sky_start[j] - sky_top[j];
        for i in sky_top[j]..j {
            let ui = sky_start[i] - sky_top[i];
            let sr, si = sky_cdot(lo_re, lo_im, ui, up_re, up_im, uj, max(sky_top[i], sky_top[j]), i);
            up_re[uj + i] -= sr;
            up_im[uj + i] -= si;
            let qr, qi = quotient(up_re[uj + i], up_im[uj + i], up_re[ui + i], up_im[ui + i]);
            lo_re[uj + i] = qr;
            lo_im[uj + i] = qi;
        }
        let size = sqrt(up_re[uj + j] * up_re[uj + j] + up_im[uj + j] * up_im[uj + j]);
        let sr, si = sky_cdot(lo_re, lo_im, uj, up_re, up_im, uj, sky_top[j], j);
        up_re[uj + j] -= sr;
        up_im[uj + j] -= si;
        let pivot = sqrt(up_re[uj + j] * up_re[uj + j] + up_im[uj + j] * up_im[uj + j]);
        if !(pivot > 1.0e-12 * size) { return false; }
    }
    return true;
}

// The right side in the order of the skyline: the load less A times the
// fixed values, on the free rows.
fn harmonic_rhs(br: f64[>..]&, bi: f64[>..]&) {
    let n = sol.len;
    for d in n { sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 }; }
    let w2 = harmonic_omega * harmonic_omega;
    solvers::zeros(br, n);
    solvers::zeros(bi, n);
    for r in n {
        guard fixed[r] == 0 else { continue; };
        var kr = 0.0;
        var ki = 0.0;
        for e in row_start[r]..row_start[r + 1] {
            let v = sol[col[e]];
            kr += (val[e] - w2 * val_b[e]) * v;
            ki += 2.0 * harmonic_damping * val[e] * v;
        }
        br[sky_row(r)] = load[r] + assembled_source[r] - kr;
        bi[sky_row(r)] = -ki;
    }
}

// The forward and backward substitution of the factor on (br, bi), then
// sol_re, sol_im and sol (the real part).
fn harmonic_back(br: f64[>..]&, bi: f64[>..]&) {
    let n = sol.len;
    for j in n {
        let sr, si = sky_cdot(lo_re, lo_im, sky_start[j] - sky_top[j], br, bi, 0, sky_top[j], j);
        br[j] -= sr;
        bi[j] -= si;
    }
    var j = n - 1;
    while j >= 0 {
        let uj = sky_start[j] - sky_top[j];
        let qr, qi = quotient(br[j], bi[j], up_re[uj + j], up_im[uj + j]);
        br[j] = qr;
        bi[j] = qi;
        for i in sky_top[j]..j {
            br[i] -= up_re[uj + i] * qr - up_im[uj + i] * qi;
            bi[i] -= up_re[uj + i] * qi + up_im[uj + i] * qr;
        }
        j--;
    }
    solvers::zeros(sol_re, n);
    solvers::zeros(sol_im, n);
    for d in n {
        sol_re[d] = if fixed[d] != 0 { prescribed[d] } else { br[sky_row(d)] };
        sol_im[d] = if fixed[d] != 0 { 0.0 } else { bi[sky_row(d)] };
    }
    harmonic_part(0);
    state_valid = true;
}

// Assembles K, M and the load, and solves A u = f at the omega and the
// damping of set_harmonic. False if a pivot is near 0.
fn harmonic_solve() -> bool {
    guard physics == ELASTIC && harmonic_omega > 0.0 && dofs_sized() else { abort("harmonic_solve: set_harmonic of an elastic physics and the arrays of the dofs"); }
    let start = clock();
    assemble("simplex elastic cg amg", 6);
    assemble_second(MASS);
    guard skyline_layout(HARMONIC_WORK) else {
        write_stderr("skein run: the harmonic step needs a factorization of too many multiplies: give a larger size\n");
        return false;
    };
    let n = sol.len;
    solvers::memory_check("simplex harmonic skyline ldlt", sky_start[n], n, 2, 0, "");
    solvers::zeros(up_re, sky_start[n]);
    solvers::zeros(up_im, sky_start[n]);
    solvers::zeros(lo_re, sky_start[n]);
    solvers::zeros(lo_im, sky_start[n]);
    var br: f64[>..] = [];
    var bi: f64[>..] = [];
    harmonic_rhs(br, bi);
    let done = harmonic_factor(harmonic_omega, harmonic_damping);
    if done { harmonic_back(br, bi); }
    last_iterations = 1;
    last_seconds = clock() - start;
    return done;
}
