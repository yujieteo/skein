// The laws of the nonlinear solid (solid.goose) at one integration point,
// and their tensor tools: 3 × 3 matrices in rows (m[3 i + j]), the
// symmetric ones also in Voigt order xx, yy, zz, xy, yz, xz, and the
// eigenvalues and eigenvectors of a symmetric matrix by Jacobi rotations.
//
//   J2 at small strain    the radial return of the von Mises yield
//                         q - (sy + Hk α) with linear isotropic hardening,
//                         and its consistent tangent (Simo and Hughes 1998,
//                         Box 3.2).
//   J2 at finite strain   the same return in the principal logarithmic
//                         strains of the elastic part Fe = F Fp⁻¹ (Hencky
//                         elasticity), which integrates the flow with the
//                         exponential map (Simo 1992; de Souza Neto, Perić
//                         and Owen 2008, Box 14.3). The state is
//                         Cp⁻¹ = Fp⁻¹ Fp⁻ᵀ, the trial elastic part is
//                         be = F Cp⁻¹ Fᵀ.
//   an energy             an isotropic hyperelastic energy W(c1, c2, c3) of
//                         the principal values of C = Fᵀ F. The program
//                         gives W in hyper-dual numbers (Fike and Alonso
//                         2011), which give its exact first and second
//                         derivatives. S and its tangent come from the
//                         spectral formulas of Holzapfel (2000, §6.6).
import std;
import math;
namespace elements;

let VOIGT_I: i64[6] = [0, 1, 2, 0, 1, 0];
let VOIGT_J: i64[6] = [0, 1, 2, 1, 2, 2];
let EYE: f64[9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
let JACOBI_SWEEPS = 50;

fn mat_mul(a: f64[9], b: f64[9]) -> f64[9] {
    var c: f64[9] = [0.0; 9];
    for i in 3 {
        for j in 3 {
            for k in 3 { c[3 * i + j] += a[3 * i + k] * b[3 * k + j]; }
        }
    }
    return c;
}

fn mat_t(a: f64[9]) -> f64[9] {
    return [a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8]];
}

fn mat_det(a: f64[9]) -> f64 {
    return a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) + a[2] * (a[3] * a[7] - a[4] * a[6]);
}

fn mat_inv(a: f64[9]) -> f64[9] {
    let det = mat_det(a);
    guard det != 0.0 else { abort("mat_inv: a matrix with det != 0"); }
    return [(a[4] * a[8] - a[5] * a[7]) / det, (a[2] * a[7] - a[1] * a[8]) / det, (a[1] * a[5] - a[2] * a[4]) / det,
            (a[5] * a[6] - a[3] * a[8]) / det, (a[0] * a[8] - a[2] * a[6]) / det, (a[2] * a[3] - a[0] * a[5]) / det,
            (a[3] * a[7] - a[4] * a[6]) / det, (a[1] * a[6] - a[0] * a[7]) / det, (a[0] * a[4] - a[1] * a[3]) / det];
}

fn to_voigt(a: f64[9]) -> f64[6] { return [a[0], a[4], a[8], a[1], a[5], a[2]]; }

fn from_voigt(v: f64[6]) -> f64[9] { return [v[0], v[3], v[5], v[3], v[1], v[4], v[5], v[4], v[2]]; }

// Σ_k d[k] n_k ⊗ n_k, with n_k column k of the rows of e.
fn spectral(d: f64[3], e: f64[9]) -> f64[9] {
    var a: f64[9] = [0.0; 9];
    for k in 3 {
        for i in 3 {
            for j in 3 { a[3 * i + j] += d[k] * e[3 * i + k] * e[3 * j + k]; }
        }
    }
    return a;
}

// The eigenvalues of the symmetric matrix a0 and its eigenvectors, as the
// columns of the second result, by cyclic Jacobi rotations (Press et al.
// 2007, §11.1).
fn sym_eigen(a0: f64[9]) -> f64[3], f64[9] {
    var a = a0;
    var v = EYE;
    var size = 0.0;
    for x in a0 { size += x * x; }
    for sweep in JACOBI_SWEEPS {
        let off = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
        if off <= 1.0e-30 * size { break; }
        for k in 3 {
            let p = if k == 2 { 1 } else { 0 };
            let q = if k == 0 { 1 } else { 2 };
            if a[3 * p + q] == 0.0 { continue; }
            let theta = (a[4 * q] - a[4 * p]) / (2.0 * a[3 * p + q]);
            let t = (if theta >= 0.0 { 1.0 } else { -1.0 }) / (abs(theta) + sqrt(theta * theta + 1.0));
            let c = 1.0 / sqrt(t * t + 1.0);
            var r = EYE;
            r[4 * p] = c;
            r[4 * q] = c;
            r[3 * p + q] = t * c;
            r[3 * q + p] = -t * c;
            a = mat_mul(mat_t(r), mat_mul(a, r));
            v = mat_mul(v, r);
        }
    }
    return [a[0], a[4], a[8]], v;
}

// --- J2 plasticity ------------------------------------------------------------

// The material of the J2 laws: the Lamé constants, the yield stress and the
// hardening modulus.
struct J2 {
    lame: f64,
    shear: f64,
    yield0: f64,
    hard: f64,
}

// The radial return of the deviatoric stress s of a trial state with the
// hardening variable alpha: the plastic multiplier dg (the increment of
// alpha), and the factor 1 - 3 μ dg / q of s. dg is 0 in the elastic range.
fn j2_return(m: J2, s: f64[9], alpha: f64) -> f64, f64, f64 {
    var ss = 0.0;
    for x in s { ss += x * x; }
    let q = sqrt(1.5 * ss);
    let f = q - (m.yield0 + m.hard * alpha);
    guard f > 1.0e-12 * m.yield0 && q > 0.0 else { return 0.0, 1.0, q; };
    let dg = f / (3.0 * m.shear + m.hard);
    return dg, 1.0 - 3.0 * m.shear * dg / q, q;
}

fn deviator(a: f64[9]) -> f64[9], f64 {
    let p = (a[0] + a[4] + a[8]) / 3.0;
    var s = a;
    for i in 3 { s[4 * i] -= p; }
    return s, p;
}

// The stress (Voigt) and the tangent of J2 at small strain for the strain
// eps, from the plastic strain ep (Voigt, tensor shears) and alpha. The
// plastic strain and alpha of the end of the increment replace ep and
// alpha.
fn j2_small(m: J2, eps: f64[9], ep: f64[6]&, alpha: f64&, stress: f64[6]&, tangent: f64[36]&) {
    var ee = eps - from_voigt(ep);
    let tr = ee[0] + ee[4] + ee[8];
    var sig = ee * (2.0 * m.shear);
    for i in 3 { sig[4 * i] += m.lame * tr; }
    let s, p = deviator(sig);
    let dg, theta, q = j2_return(m, s, alpha);
    let bulk = m.lame + 2.0 * m.shear / 3.0;
    var th_bar = 0.0;
    var n: f64[6] = [0.0; 6];
    if dg > 0.0 {
        ep = ep + to_voigt(s * (1.5 * dg / q));
        alpha += dg;
        th_bar = 1.0 / (1.0 + m.hard / (3.0 * m.shear)) - (1.0 - theta);
        n = to_voigt(s * (1.0 / (q * sqrt(2.0 / 3.0))));
    }
    stress = to_voigt(s * theta);
    for i in 3 { stress[i] += p; }
    // K 1⊗1 + 2μθ (I - 1/3 1⊗1) - 2μ θ̄ n⊗n, with ½ on the shear diagonal
    for a in 6 {
        for b in 6 {
            var dev = if a == b { if a < 3 { 1.0 } else { 0.5 } } else { 0.0 };
            if a < 3 && b < 3 { dev -= 1.0 / 3.0; }
            let vol = if a < 3 && b < 3 { bulk } else { 0.0 };
            tangent[6 * a + b] = vol + 2.0 * m.shear * theta * dev - 2.0 * m.shear * th_bar * n[a] * n[b];
        }
    }
}

// The first Piola-Kirchhoff stress P = τ F⁻ᵀ of J2 at finite strain for
// the deformation gradient f, from cpi = Cp⁻¹ (Voigt) and alpha, which the
// values of the end of the step replace. tau gets the Kirchhoff stress.
fn j2_finite(m: J2, f: f64[9], cpi: f64[6]&, alpha: f64&, tau: f64[6]&) -> f64[9] {
    let be = mat_mul(f, mat_mul(from_voigt(cpi), mat_t(f)));
    let l2, e = sym_eigen(be);
    var eps: f64[3] = [0.0; 3];
    for k in 3 { eps[k] = 0.5 * log(l2[k]); }
    let tr = eps[0] + eps[1] + eps[2];
    var t: f64[3] = [0.0; 3];
    for k in 3 { t[k] = m.lame * tr + 2.0 * m.shear * eps[k]; }
    let p = (t[0] + t[1] + t[2]) / 3.0;
    let s = [t[0] - p, 0.0, 0.0, 0.0, t[1] - p, 0.0, 0.0, 0.0, t[2] - p];
    let dg, theta, q = j2_return(m, s, alpha);
    if dg > 0.0 {
        for k in 3 { eps[k] -= 1.5 * dg * s[4 * k] / q; }
        alpha += dg;
    }
    var b2: f64[3] = [0.0; 3];
    for k in 3 {
        t[k] = p + theta * s[4 * k];
        b2[k] = exp(2.0 * eps[k]);
    }
    let finv = mat_inv(f);
    let tk = spectral(t, e);
    cpi = to_voigt(mat_mul(finv, mat_mul(spectral(b2, e), mat_t(finv))));
    tau = to_voigt(tk);
    return mat_mul(tk, mat_t(finv));
}

// --- hyper-dual numbers and the energy ----------------------------------------

// A hyper-dual number f + a ε1 + b ε2 + ab ε1 ε2 with ε1² = ε2² = 0. A
// function of it carries the exact first derivatives along the two seeds
// in a and b, and the mixed second derivative in ab. + and - and the
// product with an f64 work on each part, which is right for them.
struct Hd {
    f: f64,
    a: f64,
    b: f64,
    ab: f64,
}

fn hd(x: f64) -> Hd { Hd { f: x, a: 0.0, b: 0.0, ab: 0.0 } }

fn hd_mul(x: Hd, y: Hd) -> Hd {
    return Hd { f: x.f * y.f, a: x.a * y.f + x.f * y.a, b: x.b * y.f + x.f * y.b,
                ab: x.ab * y.f + x.a * y.b + x.b * y.a + x.f * y.ab };
}

// g(x) from g, g' and g'' at x.f.
fn hd_chain(x: Hd, g0: f64, g1: f64, g2: f64) -> Hd {
    return Hd { f: g0, a: g1 * x.a, b: g1 * x.b, ab: g1 * x.ab + g2 * x.a * x.b };
}

fn hd_div(x: Hd, y: Hd) -> Hd { return hd_mul(x, hd_chain(y, 1.0 / y.f, -1.0 / (y.f * y.f), 2.0 / (y.f * y.f * y.f))); }

fn hd_pow(x: Hd, p: f64) -> Hd {
    if p == 2.0 { return hd_mul(x, x); }
    return hd_chain(x, pow(x.f, p), p * pow(x.f, p - 1.0), p * (p - 1.0) * pow(x.f, p - 2.0));
}

fn hd_ln(x: Hd) -> Hd { return hd_chain(x, log(x.f), 1.0 / x.f, -1.0 / (x.f * x.f)); }

fn hd_exp(x: Hd) -> Hd {
    let e = exp(x.f);
    return hd_chain(x, e, e, e);
}

fn hd_sqrt(x: Hd) -> Hd {
    let r = sqrt(x.f);
    return hd_chain(x, r, 0.5 / r, -0.25 / (r * x.f));
}

// x^y: exp(y ln x), or pow for a constant y.
fn hd_pow_hd(x: Hd, y: Hd) -> Hd {
    // A constant exponent: x can be 0 or negative.
    if y.a == 0.0 && y.b == 0.0 && y.ab == 0.0 { return hd_pow(x, y.f); }
    return hd_exp(hd_mul(y, hd_ln(x)));
}

// The principal values c with the seed ε1 on c[i] and ε2 on c[j].
fn hd_seed(c: f64[3], i: i64, j: i64) -> Hd[3] {
    var h: Hd[3] = [hd(c[0]), hd(c[1]), hd(c[2])];
    h[i].a = 1.0;
    h[j].b = 1.0;
    return h;
}

// The derivatives dW/dc_i in w1 and d²W/dc_i dc_j in w2 of the energy
// ENERGY of the program at the principal values c, for i, j < n (the
// dimension: a 2d model has the principal values of the plane).
fn energy_derivatives<ENERGY>(c: f64[3], n: i64, w1: f64[3]&, w2: f64[9]&) {
    for i in n {
        for j in i..n {
            let w = ENERGY(hd_seed(c, i, j));
            w1[i] = w.a;
            w1[j] = w.b;
            w2[3 * i + j] = w.ab;
            w2[3 * j + i] = w.ab;
        }
    }
}

// The second Piola-Kirchhoff stress S (Voigt) and the tangent dS/dE of an
// isotropic energy at the principal values c and the principal directions
// e of C, from w1 and w2 of energy_derivatives, for i, j < n (Holzapfel
// 2000, eq. 6.196 and 6.197).
fn energy_stress(c: f64[3], e: f64[9], n: i64, w1: f64[3], w2: f64[9], stress: f64[6]&, tangent: f64[36]&) {
    var s: f64[3] = [0.0; 3];
    for i in n { s[i] = 2.0 * w1[i]; }
    stress = to_voigt(spectral(s, e));
    for &x in tangent { x = 0.0; }
    for i in n {
        for j in n {
            add_dyad(tangent, e, i, i, j, j, 4.0 * w2[3 * i + j]);
            guard i != j else { continue; };
            // the limit of equal values: ½ (4 w2_jj - 4 w2_ij)
            let gap = c[j] - c[i];
            let k = if abs(gap) > 1.0e-8 * max(abs(c[i]), abs(c[j])) { (s[j] - s[i]) / gap } else { 2.0 * (w2[3 * j + j] - w2[3 * i + j]) };
            add_dyad(tangent, e, i, j, i, j, k);
            add_dyad(tangent, e, i, j, j, i, k);
        }
    }
}

// Adds k (n_a ⊗ n_b) ⊗ (n_c ⊗ n_d), symmetric in each pair, to the Voigt
// tangent t.
fn add_dyad(t: f64[36]&, e: f64[9], a: i64, b: i64, c: i64, d: i64, k: f64) {
    for v in 6 {
        let i = VOIGT_I[v];
        let j = VOIGT_J[v];
        let left = 0.5 * (e[3 * i + a] * e[3 * j + b] + e[3 * j + a] * e[3 * i + b]);
        for w in 6 {
            let p = VOIGT_I[w];
            let q = VOIGT_J[w];
            t[6 * v + w] += k * left * 0.5 * (e[3 * p + c] * e[3 * q + d] + e[3 * q + c] * e[3 * p + d]);
        }
    }
}
