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
    kb[12 * i + i] += k;
    kb[12 * j + j] += k;
    kb[12 * i + j] -= k;
    kb[12 * j + i] -= k;
}

// The bending terms of dofs (w, θ) at both nodes. sign is 1 when θ = w'.
fn bend(w: i64, t: i64, ei: f64, phi: f64, l: f64, sign: f64) {
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

// The stiffness of beam b in global dofs: tbᵀ kb tb, into ke[0..144].
fn beam_stiffness(b: i64, out: f64[:]) {
    beam_local(b);
    var kt: f64[144] = [0.0; 144];
    for i in 12 {
        for j in 12 {
            var v = 0.0;
            for k in 12 { v += kb[12 * i + k] * tb[12 * k + j]; }
            kt[12 * i + j] = v;
        }
    }
    for i in 12 {
        for j in 12 {
            var v = 0.0;
            for k in 12 { v += tb[12 * k + i] * kt[12 * k + j]; }
            out[12 * i + j] = v;
        }
    }
}

// The local end forces of beam b: kb tb q. Index 6 is the axial force Nx
// at node 1 (positive in tension).
fn beam_forces(b: i64) -> f64[12] {
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
