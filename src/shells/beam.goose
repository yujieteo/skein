// The 2-node Timoshenko beam with a rigid offset at each end.
//
// The axis of the beam is at node + offset. The offset is rigid, so the
// axis point moves by u + θ × offset. The stiffness in local axes is the
// exact stiffness of a straight Timoshenko beam with the shear factor 5/6
// (J. S. Przemieniecki, Theory of Matrix Structural Analysis, 1968, §5.6).
// The local x axis runs from the axis point of node 0 to that of node 1,
// y is the y_axis of the beam made normal to x, and z = x × y.
import std;
import math;
namespace shells;
import .laminate;
import .model;

var kb: f64[144];       // 12 × 12, row 6 l + d
var tb: f64[144];       // global dofs to local dofs, with the offsets
var beam_len = 0.0;

// The local stiffness in kb and the transformation in tb of beam b.
fn beam_local(b: i64) {
    guard b >= 0 && b < beams.len && beams[b].section >= 0 && beams[b].section < beam_sections.len else { abort("beam_local: 0 <= b < beams.len, and its section < beam_sections.len"); }
    let bm = beams[b];
    let s = beam_sections[bm.section];
    let o = bm.offset;
    let p0 = add3(node_point(bm.n[0]), [o[0], o[1], o[2]]);
    let p1 = add3(node_point(bm.n[1]), [o[3], o[4], o[5]]);
    let ex = unit(sub(p1, p0));
    let y = bm.y_axis;
    let ey = unit(sub(y, scale(ex, dot3(ex, y))));
    let ez = cross(ex, ey);
    let l = norm(sub(p1, p0));
    beam_len = l;
    for &v in kb { v = 0.0; }
    let ea = s.e * s.area / l;
    let gj = s.g * s.j / l;
    set_pair(0, 6, ea);
    set_pair(3, 9, gj);
    // Bending in the x-y plane (v, θz) with Iz, and in the x-z plane
    // (w, θy) with Iy, where θy = -w'.
    bend(1, 5, s.e * s.iz, 12.0 * s.e * s.iz / (SHEAR_FACTOR * s.g * s.area * l * l), l, 1.0);
    bend(2, 4, s.e * s.iy, 12.0 * s.e * s.iy / (SHEAR_FACTOR * s.g * s.area * l * l), l, -1.0);
    // tb = blockdiag(Λ, Λ, Λ, Λ) [[I, -skew(e)], [0, I]] for each node.
    for &v in tb { v = 0.0; }
    let lam: f64[9] = [ex[0], ex[1], ex[2], ey[0], ey[1], ey[2], ez[0], ez[1], ez[2]];
    for nd in 2 {
        let e: f64[3] = [o[3 * nd], o[3 * nd + 1], o[3 * nd + 2]];
        // -skew(e) θ = θ × e
        let off: f64[9] = [0.0, e[2], -e[1], -e[2], 0.0, e[0], e[1], -e[0], 0.0];
        let k0 = 6 * nd;
        for i in 3 {
            for j in 3 {
                tb[12 * (k0 + i) + k0 + j] = lam[3 * i + j];
                tb[12 * (k0 + 3 + i) + k0 + 3 + j] = lam[3 * i + j];
                var v = 0.0;
                for k in 3 { v += lam[3 * i + k] * off[3 * k + j]; }
                tb[12 * (k0 + i) + k0 + 3 + j] = v;
            }
        }
    }
}

fn set_pair(i: i64, j: i64, k: f64) {
    guard min(i, j) >= 0 && max(i, j) < 12 else { abort("set_pair: 0 <= i, j < 12"); }
    kb[12 * i + i] += k;
    kb[12 * j + j] += k;
    kb[12 * i + j] -= k;
    kb[12 * j + i] -= k;
}

// The bending terms of dofs (w, θ) at both nodes. sign is 1 when θ = w'.
fn bend(w: i64, t: i64, ei: f64, phi: f64, l: f64, sign: f64) {
    guard min(w, t) >= 0 && max(w, t) < 6 && l > 0.0 else { abort("bend: 0 <= w, t < 6, and l > 0"); }
    let f = ei / (l * (1.0 + phi));
    let dofs = [w, t, w + 6, t + 6];
    let a = 12.0 * f / (l * l);
    let b = sign * 6.0 * f / l;
    let c = (4.0 + phi) * f;
    let d = (2.0 - phi) * f;
    let m: f64[16] = [a, b, -a, b, b, c, -b, d, -a, -b, a, -b, b, d, -b, c];
    for i in 4 {
        for j in 4 { kb[12 * dofs[i] + dofs[j]] += m[4 * i + j]; }
    }
}

// The local end forces of beam b: kb tb q. Index 6 is the axial force Nx
// at node 1 (positive in tension).
fn beam_forces(b: i64) -> f64[12] {
    guard b >= 0 && b < beams.len && disp.len == dof_count() else { abort("beam_forces: 0 <= b < beams.len, and a solution for each dof"); }
    beam_local(b);
    var q: f64[12] = [0.0; 12];
    for l in 2 {
        for d in 6 { q[6 * l + d] = disp[6 * beams[b].n[l] + d]; }
    }
    var ql: f64[12] = [0.0; 12];
    for i in 12 {
        var v = 0.0;
        for k in 12 { v += tb[12 * i + k] * q[k]; }
        ql[i] = v;
    }
    var f: f64[12] = [0.0; 12];
    for i in 12 {
        var v = 0.0;
        for k in 12 { v += kb[12 * i + k] * ql[k]; }
        f[i] = v;
    }
    return f;
}

fn beam_axial_force(b: i64) -> f64 { return beam_forces(b)[6]; }

fn add3(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }

fn scale(a: f64[3], s: f64) -> f64[3] { return [a[0] * s, a[1] * s, a[2] * s]; }

// --- the mass and the geometric stiffness ------------------------------------
//
// Both use the shape functions of the exact Timoshenko beam (Z. Friedman and
// J. B. Kosmatka, "An improved two-node Timoshenko beam finite element",
// Computers & Structures 47, 1993, 473-481): w is cubic and the rotation ψ
// is quadratic in ξ = x/l, with φ = 12 EI/(k G A l²). With φ = 0 they are
// the Hermite cubics, and the matrices are those of Przemieniecki (1968),
// §11.4 and §15.4. A 4-point Gauss rule integrates them exactly.

// The points and the weights of the 4-point Gauss rule on [0, 1].
let BEAM_GAUSS_X: f64[4] = [0.0694318442029737, 0.33000947820757187, 0.6699905217924281, 0.9305681557970263];
let BEAM_GAUSS_W: f64[4] = [0.17392742256872692, 0.3260725774312731, 0.3260725774312731, 0.17392742256872692];

// Sets the density of a beam section, for the mass of a modal step.
fn beam_density(section: i64, rho: f64) {
    guard section >= 0 && section < beam_sections.len && rho > 0.0 else { abort("beam_density: 0 <= section < beam_sections.len, and rho > 0"); }
    beam_sections[section].rho = rho;
}

// The terms of one bending plane, dofs (w, θ) at both nodes, into the
// local matrix m: ∫ (ca w w + ci ψ ψ + cg w' w') dx. sign is 1 when θ = ψ
// (the x-y plane) and -1 when θ = -ψ (the x-z plane).
fn beam_plane(m: f64[:], w: i64, t: i64, phi: f64, l: f64, sign: f64, ca: f64, ci: f64, cg: f64) {
    guard m.len == 144 && min(w, t) >= 0 && max(w, t) < 6 && l > 0.0 && phi >= 0.0 else { abort("beam_plane: m has 144 entries, 0 <= w, t < 6, l > 0 and phi >= 0"); }
    let dofs = [w, t, w + 6, t + 6];
    let s: f64[4] = [1.0, sign, 1.0, sign];
    let f = 1.0 / (1.0 + phi);
    for g in 4 {
        let x = BEAM_GAUSS_X[g];
        let dx = BEAM_GAUSS_W[g] * l;
        let nw: f64[4] = [f * (2.0 * x * x * x - 3.0 * x * x - phi * x + 1.0 + phi),
                          f * l * (x * x * x - (2.0 + 0.5 * phi) * x * x + (1.0 + 0.5 * phi) * x),
                          -f * (2.0 * x * x * x - 3.0 * x * x - phi * x),
                          f * l * (x * x * x - (1.0 - 0.5 * phi) * x * x - 0.5 * phi * x)];
        let nt: f64[4] = [6.0 * f / l * (x * x - x),
                          f * (3.0 * x * x - (4.0 + phi) * x + 1.0 + phi),
                          -6.0 * f / l * (x * x - x),
                          f * (3.0 * x * x - (2.0 - phi) * x)];
        let nd: f64[4] = [f / l * (6.0 * x * x - 6.0 * x - phi),
                          f * (3.0 * x * x - (4.0 + phi) * x + 1.0 + 0.5 * phi),
                          -f / l * (6.0 * x * x - 6.0 * x - phi),
                          f * (3.0 * x * x - (2.0 - phi) * x - 0.5 * phi)];
        for i in 4 {
            for j in 4 {
                let v = ca * nw[i] * nw[j] + ci * nt[i] * nt[j] + cg * nd[i] * nd[j];
                m[12 * dofs[i] + dofs[j]] += s[i] * s[j] * v * dx;
            }
        }
    }
}

// The consistent mass of linear N on dofs i and j: c (2, 1; 1, 2)/6, with
// c the mass of the beam (ρ A l) or its inertia about x (ρ I l).
fn beam_linear_mass(m: f64[:], i: i64, j: i64, c: f64) {
    guard m.len == 144 && min(i, j) >= 0 && max(i, j) < 12 else { abort("beam_linear_mass: m has 144 entries, and 0 <= i, j < 12"); }
    m[12 * i + i] += c / 3.0;
    m[12 * j + j] += c / 3.0;
    m[12 * i + j] += c / 6.0;
    m[12 * j + i] += c / 6.0;
}

// tbᵀ m tb into out: a local matrix of the last beam_local in global dofs.
fn beam_to_global(m: f64[:], out: f64[:]) {
    guard m.len == 144 && out.len == 144 else { abort("beam_to_global: m and out have 144 entries"); }
    var mt: f64[144] = [0.0; 144];
    for i in 12 {
        for j in 12 {
            var v = 0.0;
            for k in 12 { v += m[12 * i + k] * tb[12 * k + j]; }
            mt[12 * i + j] = v;
        }
    }
    for i in 12 {
        for j in 12 {
            var v = 0.0;
            for k in 12 { v += tb[12 * k + i] * mt[12 * k + j]; }
            out[12 * i + j] = v;
        }
    }
}

// The shear terms φ = 12 EI/(k G A l²) of the x-y plane (Iz) and of the
// x-z plane (Iy) of section s on length l.
fn beam_phi(s: BeamSection, l: f64) -> f64, f64 {
    let c = 12.0 * s.e / (SHEAR_FACTOR * s.g * s.area * l * l);
    return c * s.iz, c * s.iy;
}

// The consistent mass of beam b in global dofs, into out: the translation
// with ρA, the rotation about x with ρ(Iy + Iz), and the rotary inertia of
// bending with ρIz and ρIy. The mass is on the axis of the beam, so tb
// moves it to the nodes with the offsets.
fn beam_mass(b: i64, out: f64[:]) {
    guard b >= 0 && b < beams.len && out.len == 144 else { abort("beam_mass: 0 <= b < beams.len, and out has 144 entries"); }
    beam_local(b);
    let s = beam_sections[beams[b].section];
    let l = beam_len;
    let phi_z, phi_y = beam_phi(s, l);
    var m: f64[144] = [0.0; 144];
    beam_linear_mass(m, 0, 6, s.rho * s.area * l);
    beam_linear_mass(m, 3, 9, s.rho * (s.iy + s.iz) * l);
    beam_plane(m, 1, 5, phi_z, l, 1.0, s.rho * s.area, s.rho * s.iz, 0.0);
    beam_plane(m, 2, 4, phi_y, l, -1.0, s.rho * s.area, s.rho * s.iy, 0.0);
    beam_to_global(m, out);
}

// The geometric stiffness of beam b with the axial force nx (positive in
// tension) in global dofs, into out: nx ∫ (v'² + w'²) dx for the bending,
// and nx (Iy + Iz)/A ∫ θx'² dx for the torsion (W. McGuire, R. H. Gallagher
// and R. D. Ziemian, Matrix Structural Analysis, 2000, §9.2). It has no
// terms of the end moments, and no axial term.
fn beam_geometric(b: i64, nx: f64, out: f64[:]) {
    guard b >= 0 && b < beams.len && out.len == 144 else { abort("beam_geometric: 0 <= b < beams.len, and out has 144 entries"); }
    beam_local(b);
    let s = beam_sections[beams[b].section];
    let l = beam_len;
    let phi_z, phi_y = beam_phi(s, l);
    var m: f64[144] = [0.0; 144];
    let t = nx * (s.iy + s.iz) / (s.area * l);
    m[12 * 3 + 3] += t;
    m[12 * 9 + 9] += t;
    m[12 * 3 + 9] -= t;
    m[12 * 9 + 3] -= t;
    beam_plane(m, 1, 5, phi_z, l, 1.0, 0.0, 0.0, nx);
    beam_plane(m, 2, 4, phi_y, l, -1.0, 0.0, 0.0, nx);
    beam_to_global(m, out);
}
