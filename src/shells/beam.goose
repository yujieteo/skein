import std;
import math;
namespace shells;
import .laminate;
import .model;

var kb: f64[144];
var tb: f64[144];
var beam_len = 0.0;

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
    bend(1, 5, s.e * s.iz, 12.0 * s.e * s.iz / (SHEAR_FACTOR * s.g * s.area * l * l), l, 1.0);
    bend(2, 4, s.e * s.iy, 12.0 * s.e * s.iy / (SHEAR_FACTOR * s.g * s.area * l * l), l, -1.0);
    for &v in tb { v = 0.0; }
    let lam: f64[9] = [ex[0], ex[1], ex[2], ey[0], ey[1], ey[2], ez[0], ez[1], ez[2]];
    for nd in 2 {
        let e: f64[3] = [o[3 * nd], o[3 * nd + 1], o[3 * nd + 2]];
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

let BEAM_GAUSS_X: f64[4] = [0.0694318442029737, 0.33000947820757187, 0.6699905217924281, 0.9305681557970263];
let BEAM_GAUSS_W: f64[4] = [0.17392742256872692, 0.3260725774312731, 0.3260725774312731, 0.17392742256872692];

fn beam_density(section: i64, rho: f64) {
    guard section >= 0 && section < beam_sections.len && rho > 0.0 else { abort("beam_density: 0 <= section < beam_sections.len, and rho > 0"); }
    beam_sections[section].rho = rho;
}

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

fn beam_linear_mass(m: f64[:], i: i64, j: i64, c: f64) {
    guard m.len == 144 && min(i, j) >= 0 && max(i, j) < 12 else { abort("beam_linear_mass: m has 144 entries, and 0 <= i, j < 12"); }
    m[12 * i + i] += c / 3.0;
    m[12 * j + j] += c / 3.0;
    m[12 * i + j] += c / 6.0;
    m[12 * j + i] += c / 6.0;
}

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

fn beam_phi(s: BeamSection, l: f64) -> f64, f64 {
    let c = 12.0 * s.e / (SHEAR_FACTOR * s.g * s.area * l * l);
    return c * s.iz, c * s.iy;
}

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
