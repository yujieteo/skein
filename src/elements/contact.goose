// Contact with penalty friction (spec §12.4), and the rigid face of
// `u = rigid(p)` (spec §12.1), as penalty terms of the solid of solid.goose.
//
// The contact is surface to surface in the small sliding form: each
// quadrature point of a facet of the slave region is paired once, on the
// mesh at u = 0, with the point of the nearest facet of the master region
// below it along the master normal n. The point keeps the nodes and the
// shape values of both facets, so the relative displacement of the pair
// is r = Σ w_k u_k, with w_k = N_k on the slave and -N_k on the master. The
// normal gap is g = g0 + n·r, and the pressure is p = kn <-g>. The
// tangential traction t, with the stick stiffness kt (femahoy gap G23),
// is the return map of Coulomb friction:
//
//     t_trial = t_c + kt T (r - r_c),  t = t_trial  if |t_trial| <= mu p
//                                      t = mu p t_trial / |t_trial|  else
//
// with T = I - n nᵀ and t_c, r_c of the end of the last increment. The
// internal force of node k is w_k A (t - p n) for the area A of the point.
// The tangent is w_k w_j A (kn n nᵀ + D), with D = kt T in stick and the
// symmetric part kt mu p / |t_trial| (T - t̂ t̂ᵀ) in slip, so that CG can
// solve it. Reason: the pairs of the femahoy contact cases slide less than
// a facet, and a fixed pairing keeps the terms on a fixed pattern of the
// matrix, as on the GPU (src/gpu/).
//
// A rigid face holds the displacement of its nodes to a rigid motion by
// the penalty k (I - Q Qᵀ), with Q the orthonormal rigid motions of the
// nodes of the face: the translations and the rotations about the
// centroid. k is RIGID_SCALE times the largest diagonal term of K on the
// face, so the face deforms by about 1/RIGID_SCALE of an elastic face.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;

let CT_SLOTS = 12;          // the nodes of the slave facet and of the master facet of a point
let RIGID_SCALE = 1.0e3;
let RIGID_MODES = 6;

var ct_master: i64[>..] = [];   // 4 c + f of each facet of the master of the next contact_slave
var ct_node: i64[>..] = [];     // CT_SLOTS nodes for each point, -1 for an empty slot
var ct_w: f64[>..] = [];        // the weight w of each slot
var ct_n: f64[>..] = [];        // 3 for each point: the outward normal of the master
var ct_area: f64[>..] = [];     // the area of each point
var ct_g0: f64[>..] = [];       // the gap of each point at u = 0
var ct_law: f64[>..] = [];      // 3 for each point: kn, kt and mu
var ct_t: f64[>..] = [];        // 3 for each point: t at the end of the last increment
var ct_r: f64[>..] = [];        // 3 for each point: r there
var ct_trial_t: f64[>..] = [];  // t and r of the last evaluation
var ct_trial_r: f64[>..] = [];
var ct_p: f64[>..] = [];        // the pressure of each point at the last evaluation

var rigid_nodes: i64[>..] = []; // the nodes of the rigid face of the step, or empty
var rigid_q: f64[>..] = [];     // RIGID_MODES columns of ncomp × rigid_nodes.len values
var rigid_modes = 0;
var rigid_k = 0.0;              // the penalty, or 0 before the first assembly of the step

fn contact_points() -> i64 { return ct_area.len; }

// Whether the solve of a static step needs the Newton iterations of
// solid.goose for the contact or the rigid face.
fn has_ties() -> bool { return contact_points() > 0 || rigid_nodes.len > 0; }

// --- the pairs -----------------------------------------------------------------

// The facets of the region are the master of the next contact_slave.
fn set_master() -> i64 {
    guard physics == ELASTIC else { abort("contact_master: an elastic physics"); }
    ct_master.clear();
    for c in cell_count() {
        for f in facet_count(dim) { if region_facet(c, f) { ct_master.push(4 * c + f); } }
    }
    return ct_master.len;
}

// The quadrature of a slave facet: 4 Gauss points on a line, and the 6
// points of degree 4 of Dunavant (1985) on a triangle. Reason: the
// pressure of a P2 facet is a curve, and its edge has a kink in it.
let CT_LINE: f64[8] = [0.0694318442029737, 0.1739274225687269, 0.3300094782075719, 0.3260725774312731,
                       0.6699905217924281, 0.3260725774312731, 0.9305681557970263, 0.1739274225687269];
let CT_TRI: f64[18] = [0.445948490915965, 0.445948490915965, 0.111690794839005,
                       0.445948490915965, 0.108103018168070, 0.111690794839005,
                       0.108103018168070, 0.445948490915965, 0.111690794839005,
                       0.091576213509771, 0.091576213509771, 0.054975871827661,
                       0.091576213509771, 0.816847572980459, 0.054975871827661,
                       0.816847572980459, 0.091576213509771, 0.054975871827661];

fn ct_quad_count() -> i64 { return if dim == 2 { 4 } else { 6 }; }

fn ct_quad_point(q: i64) -> f64[3], f64 {
    if dim == 2 { return [CT_LINE[2 * q], 0.0, 0.0], CT_LINE[2 * q + 1]; }
    return [CT_TRI[3 * q], CT_TRI[3 * q + 1], 0.0], CT_TRI[3 * q + 2];
}

// The point x, the shape values sv of the facet nodes, the tangents and
// the area factor of facet f of cell c at xi.
fn facet_at(c: i64, f: i64, xi: f64[3], sv: f64[10]&, tan: f64[6]&) -> f64[3], f64 {
    let fd = dim - 1;
    var dref: f64[30] = [0.0; 30];
    shape(fd, order, xi, sv, dref);
    var x: f64[3] = [0.0; 3];
    for &t in tan { t = 0.0; }
    for l in node_count_of(fd, order) {
        let n = cell_node(c, facet_node(dim, f, l));
        for t in 3 {
            x[t] += sv[l] * node_x(n, t);
            for r in fd { tan[3 * r + t] += dref[3 * l + r] * node_x(n, t); }
        }
    }
    if fd == 1 { return x, sqrt(tan[0] * tan[0] + tan[1] * tan[1]) * thickness; }
    let nx = tan[1] * tan[5] - tan[2] * tan[4];
    let ny = tan[2] * tan[3] - tan[0] * tan[5];
    let nz = tan[0] * tan[4] - tan[1] * tan[3];
    return x, sqrt(nx * nx + ny * ny + nz * nz);
}

// The projection of x on the plane of the vertices of master facet m: the
// reference point xi of the facet, the gap along its outward normal n, and
// whether the point is in the facet.
fn project(m: i64, x: f64[3]) -> f64[3], f64, f64[3], bool {
    let c = m / 4;
    let f = m % 4;
    var v: f64[9] = [0.0; 9];
    for l in dim {
        let n = cell_node(c, facet_node(dim, f, l));
        for t in 3 { v[3 * l + t] = node_x(n, t); }
    }
    var e: f64[6] = [0.0; 6];
    for r in dim - 1 {
        for t in 3 { e[3 * r + t] = v[3 * (r + 1) + t] - v[t]; }
    }
    let n = facet_normal(c, [v[0], v[1], v[2]], e);
    var d: f64[3] = [0.0; 3];
    for t in 3 { d[t] = x[t] - v[t]; }
    let gap = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
    var xi: f64[3] = [0.0; 3];
    // the least squares of d - gap n = Σ xi_r e_r
    var a: f64[4] = [0.0; 4];
    var b: f64[2] = [0.0; 2];
    for r in dim - 1 {
        for t in 3 { b[r] += e[3 * r + t] * (d[t] - gap * n[t]); }
        for s in dim - 1 {
            for t in 3 { a[2 * r + s] += e[3 * r + t] * e[3 * s + t]; }
        }
    }
    if dim == 2 {
        xi[0] = b[0] / a[0];
    } else {
        let det = a[0] * a[3] - a[1] * a[2];
        xi[0] = (a[3] * b[0] - a[1] * b[1]) / det;
        xi[1] = (a[0] * b[1] - a[2] * b[0]) / det;
    }
    let tol = 1.0e-6;
    let inside = xi[0] >= -tol && xi[1] >= -tol && xi[0] + xi[1] <= 1.0 + tol;
    return xi, gap, n, inside;
}

// Pairs each quadrature point of the facets of the region, the slave, with
// the master of the last contact_master, with the penalty kn, the stick
// stiffness kt and the friction coefficient mu. It gives the points
// paired.
fn add_slave(kn: f64, kt: f64, mu: f64) -> i64 {
    guard physics == ELASTIC && dofs_sized() && kn > 0.0 && kt >= 0.0 && mu >= 0.0 else { abort("contact_slave: an elastic physics, kn > 0, kt >= 0 and mu >= 0"); }
    let first = contact_points();
    let fd = dim - 1;
    let fn_count = node_count_of(fd, order);
    var sv: f64[10] = [0.0; 10];
    var mv: f64[10] = [0.0; 10];
    var dref: f64[30] = [0.0; 30];
    var tan: f64[6] = [0.0; 6];
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            for q in ct_quad_count() {
                let xi, w = ct_quad_point(q);
                let x, da = facet_at(c, f, xi, sv, tan);
                var best = -1;
                var best_xi: f64[3] = [0.0; 3];
                var best_gap = 0.0;
                var best_n: f64[3] = [0.0; 3];
                for m in ct_master {
                    let mxi, gap, n, inside = project(m, x);
                    if inside && (best < 0 || abs(gap) < abs(best_gap)) {
                        best = m;
                        best_xi = mxi;
                        best_gap = gap;
                        best_n = n;
                    }
                }
                if best < 0 { continue; }
                shape(fd, order, best_xi, mv, dref);
                for k in CT_SLOTS {
                    let on_slave = k < fn_count;
                    let l = if on_slave { k } else { k - fn_count };
                    let used = l < fn_count;
                    let cell = if on_slave { c } else { best / 4 };
                    let facet = if on_slave { f } else { best % 4 };
                    ct_node.push(if used { cell_node(cell, facet_node(dim, facet, l)) } else { -1 });
                    ct_w.push(if !used { 0.0 } else if on_slave { sv[l] } else { -mv[l] });
                }
                ct_n.append(best_n);
                ct_area.push(w * da);
                ct_g0.push(best_gap);
                ct_law.append([kn, kt, mu]);
                ct_t.append([0.0; 3]);
                ct_r.append([0.0; 3]);
                ct_trial_t.append([0.0; 3]);
                ct_trial_r.append([0.0; 3]);
                ct_p.push(0.0);
            }
        }
    }
    return contact_points() - first;
}

// The node pairs of the terms of the contact and the rigid face, for the
// pattern of the matrix (build_pattern).
fn tie_keys(keys: i64[>..]&) {
    let nn = node_count();
    for i in contact_points() {
        for a in CT_SLOTS {
            let na = ct_node[CT_SLOTS * i + a];
            if na < 0 { continue; }
            for b in CT_SLOTS {
                let nb = ct_node[CT_SLOTS * i + b];
                if nb >= 0 { keys.push(na * nn + nb); }
            }
        }
    }
    for a in rigid_nodes {
        for b in rigid_nodes { keys.push(a * nn + b); }
    }
}

// --- the terms -----------------------------------------------------------------

// The relative displacement r of point i at u.
fn relative(i: i64, u: f64[:]) -> f64[3] {
    var r: f64[3] = [0.0; 3];
    for k in CT_SLOTS {
        let n = ct_node[CT_SLOTS * i + k];
        if n < 0 { continue; }
        for t in ncomp { r[t] += ct_w[CT_SLOTS * i + k] * u[ncomp * n + t]; }
    }
    return r;
}

// The traction t - p n of point i at the relative displacement r into the
// trial state, and the 3 × 3 tangent of the traction in d.
fn point_traction(i: i64, r: f64[3], d: f64[9]&) -> f64[3] {
    var n: f64[3] = [ct_n[3 * i], ct_n[3 * i + 1], ct_n[3 * i + 2]];
    let kn = ct_law[3 * i];
    let kt = ct_law[3 * i + 1];
    let mu = ct_law[3 * i + 2];
    for &x in d { x = 0.0; }
    let rn = n[0] * r[0] + n[1] * r[1] + n[2] * r[2];
    let g = ct_g0[i] + rn;
    var dr: f64[3] = [0.0; 3];
    for t in 3 {
        ct_trial_r[3 * i + t] = r[t];
        dr[t] = r[t] - ct_r[3 * i + t];
    }
    let drn = n[0] * dr[0] + n[1] * dr[1] + n[2] * dr[2];
    var tr: f64[3] = [0.0; 3];
    for t in 3 { tr[t] = ct_t[3 * i + t] + kt * (dr[t] - drn * n[t]); }
    let p = if g < 0.0 { -kn * g } else { 0.0 };
    ct_p[i] = p;
    let size = sqrt(tr[0] * tr[0] + tr[1] * tr[1] + tr[2] * tr[2]);
    let slip = size > mu * p;
    var tt: f64[3] = [0.0; 3];
    if p > 0.0 && mu > 0.0 {
        for t in 3 { tt[t] = if slip { mu * p * tr[t] / size } else { tr[t] }; }
    }
    for t in 3 { ct_trial_t[3 * i + t] = tt[t]; }
    guard p > 0.0 else { return [0.0; 3]; }
    let ks = if mu <= 0.0 { 0.0 } else if slip { kt * mu * p / size } else { kt };
    for a in 3 {
        for b in 3 {
            let tb = (if a == b { 1.0 } else { 0.0 }) - n[a] * n[b];
            let hat = if slip && mu > 0.0 { tr[a] * tr[b] / (size * size) } else { 0.0 };
            d[3 * a + b] = kn * n[a] * n[b] + ks * (tb - hat);
        }
    }
    return [tt[0] - p * n[0], tt[1] - p * n[1], tt[2] - p * n[2]];
}

// Adds the internal forces of the contact at u to out, and with `tangent`
// their tangent to val.
fn contact_terms(u: f64[:], out: f64[:], tangent: bool) {
    var d: f64[9] = [0.0; 9];
    for i in contact_points() {
        let f = point_traction(i, relative(i, u), d);
        let a = ct_area[i];
        for k in CT_SLOTS {
            let nk = ct_node[CT_SLOTS * i + k];
            if nk < 0 { continue; }
            let wk = ct_w[CT_SLOTS * i + k];
            for t in ncomp { out[ncomp * nk + t] += wk * a * f[t]; }
            guard tangent && ct_p[i] > 0.0 else { continue; }
            for j in CT_SLOTS {
                let nj = ct_node[CT_SLOTS * i + j];
                if nj < 0 { continue; }
                let s = wk * ct_w[CT_SLOTS * i + j] * a;
                for p in ncomp {
                    let e = entry(ncomp * nk + p, ncomp * nj);
                    for q in ncomp { val[e + q] += s * d[3 * p + q]; }
                }
            }
        }
    }
}

// The state of the contact at the end of an increment or a time step.
fn commit_contact() {
    for k in ct_t.len {
        ct_t[k] = ct_trial_t[k];
        ct_r[k] = ct_trial_r[k];
    }
}

// The bound Σ_e |K_de| of the contact terms of each dof into rate, at the
// largest stiffness of each point, for the stable step of
// central_difference.
fn contact_rows(rate: f64[:]) {
    for i in contact_points() {
        let kn = ct_law[3 * i];
        let kt = if ct_law[3 * i + 2] > 0.0 { ct_law[3 * i + 1] } else { 0.0 };
        var sum = 0.0;
        for k in CT_SLOTS { sum += abs(ct_w[CT_SLOTS * i + k]); }
        for k in CT_SLOTS {
            let nk = ct_node[CT_SLOTS * i + k];
            if nk < 0 { continue; }
            for p in ncomp {
                var row = 0.0;
                for q in ncomp {
                    let npq = ct_n[3 * i + p] * ct_n[3 * i + q];
                    row += abs(kn * npq + kt * ((if p == q { 1.0 } else { 0.0 }) - npq));
                }
                rate[ncomp * nk + p] += abs(ct_w[CT_SLOTS * i + k]) * sum * ct_area[i] * row;
            }
        }
    }
}

// The largest contact pressure of the points.
fn max_contact_pressure() -> f64 {
    var p = 0.0;
    for x in ct_p { p = max(p, x); }
    return p;
}

// --- the rigid face ---------------------------------------------------------------

// The nodes of the region are the rigid face of the step. Q gets the
// rigid motions of these nodes, orthonormal by Gram-Schmidt.
fn set_rigid() -> i64 {
    guard physics == ELASTIC && dofs_sized() else { abort("rigid_face: an elastic physics"); }
    rigid_nodes.clear();
    for n in node_count() { if on_region(n) { rigid_nodes.push(n); } }
    rigid_k = 0.0;
    rigid_q.clear();
    rigid_modes = 0;
    let m = ncomp * rigid_nodes.len;
    let centre = rigid_centre();
    let rotations = if dim == 2 { 1 } else { 3 };
    var col: f64[>..] = [];
    for mode in dim + rotations {
        col.clear();
        for n in rigid_nodes {
            var r: f64[3] = [0.0; 3];
            for t in 3 { r[t] = node_x(n, t) - centre[t]; }
            for t in ncomp {
                if mode < dim {
                    col.push(if t == mode { 1.0 } else { 0.0 });
                } else {
                    // the rotation about axis a: a × r
                    let a = if dim == 2 { 2 } else { mode - dim };
                    var w: f64[3] = [0.0; 3];
                    w[a] = 1.0;
                    let x = [w[1] * r[2] - w[2] * r[1], w[2] * r[0] - w[0] * r[2], w[0] * r[1] - w[1] * r[0]];
                    col.push(x[t]);
                }
            }
        }
        for j in rigid_modes {
            var s = 0.0;
            for i in m { s += rigid_q[j * m + i] * col[i]; }
            for i in m { col[i] -= s * rigid_q[j * m + i]; }
        }
        var size = 0.0;
        for x in col { size += x * x; }
        // A face on a line has no rotation about the line.
        if size <= 1.0e-12 * m { continue; }
        for x in col { rigid_q.push(x / sqrt(size)); }
        rigid_modes++;
    }
    return rigid_nodes.len;
}

fn rigid_centre() -> f64[3] {
    var c: f64[3] = [0.0; 3];
    for n in rigid_nodes {
        for t in 3 { c[t] += node_x(n, t) / rigid_nodes.len; }
    }
    return c;
}

// Adds the moment mv about the rigid face as forces on its nodes: f_i =
// w × r_i about the centroid, with J w = mv for J = Σ (|r|² I - r rᵀ).
// The forces sum to 0, and their moment is mv about each point.
fn add_rigid_moment(mv: f64[3]) {
    guard rigid_nodes.len > 1 && dofs_sized() else { abort("rigid_moment: a rigid face of the step"); }
    let centre = rigid_centre();
    var j: f64[9] = [0.0; 9];
    for n in rigid_nodes {
        var r: f64[3] = [0.0; 3];
        for t in 3 { r[t] = node_x(n, t) - centre[t]; }
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        for a in 3 {
            for b in 3 { j[3 * a + b] += (if a == b { r2 } else { 0.0 }) - r[a] * r[b]; }
        }
    }
    var w: f64[3] = [0.0; 3];
    if dim == 2 {
        w[2] = mv[2] / j[8];
    } else {
        w = mat_solve3(j, mv);
    }
    for n in rigid_nodes {
        var r: f64[3] = [0.0; 3];
        for t in 3 { r[t] = node_x(n, t) - centre[t]; }
        let f = [w[1] * r[2] - w[2] * r[1], w[2] * r[0] - w[0] * r[2], w[0] * r[1] - w[1] * r[0]];
        for t in ncomp { load[ncomp * n + t] += f[t]; }
    }
}

// x of a x = b for a 3 × 3 matrix a, by Cramer's rule.
fn mat_solve3(a: f64[9], b: f64[3]) -> f64[3] {
    let det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) + a[2] * (a[3] * a[7] - a[4] * a[6]);
    guard abs(det) > 0.0 else { abort("rigid_moment: a face that is not on one line"); }
    var x: f64[3] = [0.0; 3];
    for k in 3 {
        var m = a;
        for r in 3 { m[3 * r + k] = b[r]; }
        x[k] = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
    }
    return x;
}

// Adds the internal forces k (I - Q Qᵀ) u of the rigid face to out, and
// with `tangent` the penalty to val. The first call of a step sets k from
// the diagonal of val, which then holds the terms of the cells.
fn rigid_terms(u: f64[:], out: f64[:], tangent: bool) {
    guard rigid_nodes.len > 0 else { return; }
    let m = ncomp * rigid_nodes.len;
    if rigid_k == 0.0 {
        guard tangent else { abort("rigid_terms: the first call of a step assembles the tangent"); }
        for n in rigid_nodes {
            for t in ncomp { rigid_k = max(rigid_k, RIGID_SCALE * val[entry(ncomp * n + t, ncomp * n + t)]); }
        }
    }
    var coef: f64[RIGID_MODES] = [0.0; RIGID_MODES];
    for j in rigid_modes {
        for i in m { coef[j] += rigid_q[j * m + i] * u[ncomp * rigid_nodes[i / ncomp] + i % ncomp]; }
    }
    for i in m {
        let di = ncomp * rigid_nodes[i / ncomp] + i % ncomp;
        var s = u[di];
        for j in rigid_modes { s -= rigid_q[j * m + i] * coef[j]; }
        out[di] += rigid_k * s;
        guard tangent else { continue; }
        for k in m {
            let dk = ncomp * rigid_nodes[k / ncomp] + k % ncomp;
            var q = if i == k { 1.0 } else { 0.0 };
            for j in rigid_modes { q -= rigid_q[j * m + i] * rigid_q[j * m + k]; }
            val[entry(di, dk)] += rigid_k * q;
        }
    }
}

// Removes the contact points and the rigid face, for a new physics.
fn clear_contact() {
    ct_master.clear();
    ct_node.clear();
    ct_w.clear();
    ct_n.clear();
    ct_area.clear();
    ct_g0.clear();
    ct_law.clear();
    ct_t.clear();
    ct_r.clear();
    ct_trial_t.clear();
    ct_trial_r.clear();
    ct_p.clear();
    clear_rigid();
}

// Removes the rigid face, for the next step.
fn clear_rigid() {
    rigid_nodes.clear();
    rigid_q.clear();
    rigid_modes = 0;
    rigid_k = 0.0;
}
