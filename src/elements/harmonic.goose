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

let HARMONIC_WORK = 4.0e9;
let PIVOT_TOLERANCE = 1.0e-12;

var harmonic_damping = 0.0;
var sol_re: f64[>..] = [];
var sol_im: f64[>..] = [];
var up_re: f64[>..] = [];
var up_im: f64[>..] = [];
var lo_re: f64[>..] = [];
var lo_im: f64[>..] = [];

fn set_harmonic(omega: f64, damping: f64) {
    guard physics == ELASTIC && omega > 0.0 && damping >= 0.0 && density > 0.0 else {
        abort("set_harmonic: an elastic physics with a density, omega > 0 and damping >= 0");
    }
    harmonic_omega = omega;
    harmonic_damping = damping;
}

fn harmonic_part(p: i64) {
    guard (p == 0 || p == 1) && sol_re.len == sol.len && sol_im.len == sol.len else { abort("harmonic_part: part 0 or 1 of a solved harmonic step"); }
    for d in sol.len { sol[d] = if p == 0 { sol_re[d] } else { sol_im[d] }; }
}

fn harmonic_at(p: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    harmonic_part(p);
    return field_at(FIELD_VALUE, comp, x, y, z);
}

fn harmonic_mean(p: i64, comp: i64) -> f64 {
    harmonic_part(p);
    return region_mean(comp);
}

fn sky_cdot(lr: f64[:], li: f64[:], a0: i64, ur: f64[:], ui: f64[:], b0: i64, lo: i64, hi: i64) -> f64, f64 {
    var re = 0.0;
    var im = 0.0;
    for k in lo..hi {
        re += lr[a0 + k] * ur[b0 + k] - li[a0 + k] * ui[b0 + k];
        im += lr[a0 + k] * ui[b0 + k] + li[a0 + k] * ur[b0 + k];
    }
    return re, im;
}

fn quotient(ar: f64, ai: f64, br: f64, bi: f64) -> f64, f64 {
    let d = br * br + bi * bi;
    return (ar * br + ai * bi) / d, (ai * br - ar * bi) / d;
}

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
        if !(pivot > PIVOT_TOLERANCE * size) { return false; }
    }
    return true;
}

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
