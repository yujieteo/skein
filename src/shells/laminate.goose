// Shell sections: the ABD matrix of a laminate (classical laminate theory)
// and the transverse shear stiffness, in the axes of the section.
//
// A section holds, in the local axes (1, 2) of the shell surface,
//     N = A eps + B kappa,   M = B eps + D kappa,   Q = S gamma
// with eps = (eps_11, eps_22, gamma_12), kappa the same for the curvature,
// and gamma = (gamma_13, gamma_23). Each 3 × 3 matrix is row-major.
//
// The plies stack from the top face (z = +h/2) down. A ply angle turns
// from the 1 axis of the section to the fibre, right-handed about the top
// normal (spec §4.2, `reissner_mindlin`). Reference: R. M. Jones,
// Mechanics of Composite Materials, 2nd ed., §2.6 (Q-bar) and §4.2 (ABD).
import std;
import math;
namespace shells;

let SHEAR_FACTOR = 5.0 / 6.0;   // the transverse shear correction factor

struct Section {
    a: f64[9],
    b: f64[9],
    d: f64[9],
    s: f64[4],
    h: f64,
}

struct Ply {
    angle: f64,     // rad
    t: f64,
    e1: f64,
    e2: f64,
    g12: f64,
    nu12: f64,
    g13: f64,
    g23: f64,
}

var sections: Section[>..] = [];
var plies: Ply[>..] = [];

// The reduced stiffness of a ply in the axes of the section (Jones eq. 2.85):
// q[0..9] for (11, 22, 12) and s[0..4] for (13, 23).
fn ply_qbar(p: Ply) -> f64[9], f64[4] {
    let nu21 = p.nu12 * p.e2 / p.e1;
    let den = 1.0 - p.nu12 * nu21;
    let q11 = p.e1 / den;
    let q22 = p.e2 / den;
    let q12 = p.nu12 * p.e2 / den;
    let q66 = p.g12;
    let m = cos(p.angle);
    let n = sin(p.angle);
    let m2 = m * m;
    let n2 = n * n;
    let b11 = q11 * m2 * m2 + 2.0 * (q12 + 2.0 * q66) * m2 * n2 + q22 * n2 * n2;
    let b22 = q11 * n2 * n2 + 2.0 * (q12 + 2.0 * q66) * m2 * n2 + q22 * m2 * m2;
    let b12 = (q11 + q22 - 4.0 * q66) * m2 * n2 + q12 * (m2 * m2 + n2 * n2);
    let b66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * m2 * n2 + q66 * (m2 * m2 + n2 * n2);
    let b16 = (q11 - q12 - 2.0 * q66) * m2 * m * n + (q12 - q22 + 2.0 * q66) * m * n2 * n;
    let b26 = (q11 - q12 - 2.0 * q66) * m * n2 * n + (q12 - q22 + 2.0 * q66) * m2 * m * n;
    let s55 = p.g13 * m2 + p.g23 * n2;
    let s44 = p.g23 * m2 + p.g13 * n2;
    let s45 = (p.g13 - p.g23) * m * n;
    return [b11, b12, b16, b12, b22, b26, b16, b26, b66], [s55, s45, s45, s44];
}

// --- the laminate builder ------------------------------------------------------

fn laminate_clear() { plies.clear(); }

// Adds a ply below the plies so far. The angle is in rad. G13 and G23 are
// the transverse shear moduli.
fn laminate_ply(angle: f64, t: f64, e1: f64, e2: f64, g12: f64, nu12: f64, g13: f64, g23: f64) {
    assert(t > 0.0 && e1 > 0.0 && e2 > 0.0 && g12 > 0.0 && g13 > 0.0 && g23 > 0.0);
    plies.push(Ply { angle: angle, t: t, e1: e1, e2: e2, g12: g12, nu12: nu12, g13: g13, g23: g23 });
}

// Adds the plies so far again in the reverse sequence: a symmetric laminate.
fn laminate_mirror() {
    let n = plies.len;
    for k in n { plies.push(plies[n - 1 - k]); }
}

// The section of the plies so far. Returns its index.
fn laminate_section() -> i64 {
    var s = Section { a: [0.0; 9], b: [0.0; 9], d: [0.0; 9], s: [0.0; 4], h: 0.0 };
    for p in plies { s.h += p.t; }
    var top = 0.5 * s.h;
    for p in plies {
        let q, qs = ply_qbar(p);
        let bot = top - p.t;
        for i in 9 {
            s.a[i] += q[i] * (top - bot);
            s.b[i] += q[i] * (top * top - bot * bot) / 2.0;
            s.d[i] += q[i] * (top * top * top - bot * bot * bot) / 3.0;
        }
        for i in 4 { s.s[i] += SHEAR_FACTOR * qs[i] * (top - bot); }
        top = bot;
    }
    sections.push(s);
    return sections.len - 1;
}

// The section of an isotropic shell of thickness h.
fn isotropic_section(e: f64, nu: f64, h: f64) -> i64 {
    laminate_clear();
    let g = e / (2.0 * (1.0 + nu));
    laminate_ply(0.0, h, e, e, g, nu, g, g);
    return laminate_section();
}

// The membrane compliance a = A⁻¹ of section k (row-major 3 × 3).
fn membrane_compliance(k: i64) -> f64[9] {
    let a = sections[k].a;
    let det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6]);
    return [
        (a[4] * a[8] - a[5] * a[7]) / det, (a[2] * a[7] - a[1] * a[8]) / det, (a[1] * a[5] - a[2] * a[4]) / det,
        (a[5] * a[6] - a[3] * a[8]) / det, (a[0] * a[8] - a[2] * a[6]) / det, (a[2] * a[3] - a[0] * a[5]) / det,
        (a[3] * a[7] - a[4] * a[6]) / det, (a[1] * a[6] - a[0] * a[7]) / det, (a[0] * a[4] - a[1] * a[3]) / det,
    ];
}
