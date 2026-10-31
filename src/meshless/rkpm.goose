// rkpm: the Galerkin form with the reproducing kernel shape functions of
// rk.goose (spec §17.5).
//
// The stiffness takes stabilized conforming nodal integration (SCNI, Chen
// et al. 2001): the gradient of each shape function at point I is its
// mean over the nodal domain of I, the integral of Ψ n on the boundary of
// that domain divided by its area. The domain of I in each of its
// triangles is the quadrilateral of I, the midpoints of its two edges and
// the centroid; each side takes 2 Gauss points. The smoothed gradients
// pass the patch test of order 1 exactly. The source and the advection
// take the Gauss rule of 3 points on each triangle.
//
// The shape functions do not interpolate, so a fixed value is weak, by the
// symmetric method of Nitsche (1971): on an edge where component c of u is
// g, the form adds
//
//     - ∫ (σ(v) n)_c (u_c - g) - ∫ (σ(u) n)_c v_c + β ∫ (u_c - g) v_c
//
// with β = NITSCHE × (λ + 2μ)/h_e (or d/h_e for the scalar physics), h_e
// the length of the edge. Each of these integrals takes 2 Gauss points on
// each half of the edge.
//
// Each integration point is a record: its kind, its point, its weight,
// and the points whose functions are not 0 there, with Ψ and its gradient
// (or the smoothed gradient). The records do not change with the
// conditions, so rkpm() builds them once, with the matrix pattern.
import std;
import math;
namespace meshless;
import .cloud;
import .rk;
import .model;

// The factor of the Nitsche penalty. Reason: with 4, the patch test and
// the orders of test.goose pass, and a larger factor makes the matrix
// worse conditioned (CG with AMG took 831 iterations with 4 and 2155 with
// 100 on the 32 × 32 cloud).
let NITSCHE = 4.0;

let SMOOTH = 0;   // the nodal domain of a point: the smoothed gradient
let DOMAIN = 1;   // a Gauss point of a triangle
let EDGE = 2;     // a Gauss point of a boundary edge

let GAUSS_2: f64[2] = [0.21132486540518713, 0.7886751345948129];   // on [0, 1], weight 1/2 each

var rec_start: i64[>..] = [];   // the entries of each record, from rec_start[r]
var rec_kind: i64[>..] = [];
var rec_x: f64[>..] = [];
var rec_y: f64[>..] = [];
var rec_w: f64[>..] = [];       // the area or the length that the point stands for
var rec_edge: i64[>..] = [];    // the edge of an EDGE record, else -1
var rec_owner: i64[>..] = [];   // the SMOOTH record whose side holds an EDGE record, else -1
var rec_g: f64[>..] = [];       // ncomp for each record: the fixed value of the last solve
var item: i64[>..] = [];        // the point of each entry
var item_psi: f64[>..] = [];
var item_dx: f64[>..] = [];
var item_dy: f64[>..] = [];

// The scratch of the smoothed gradients: the sum of each point, and the
// points with a sum.
var acc: f64[>..] = [];
var touched: i64[>..] = [];

fn record_count() -> i64 { return rec_kind.len; }

fn open_record(kind: i64, x: f64, y: f64, w: f64, e: i64, owner: i64) {
    rec_kind.push(kind);
    rec_owner.push(owner);
    rec_x.push(x);
    rec_y.push(y);
    rec_w.push(w);
    rec_edge.push(e);
    for c in ncomp { rec_g.push(0.0); }
}

// Adds a record at (x, y) with the shape functions there. False if M is
// singular there.
fn shape_record(kind: i64, x: f64, y: f64, w: f64, e: i64, owner: i64) -> bool {
    if !shape_at(x, y) { return false; }
    open_record(kind, x, y, w, e, owner);
    item.append(near);
    item_psi.append(psi);
    item_dx.append(psi_x);
    item_dy.append(psi_y);
    rec_start.push(item.len);
    return true;
}

// Adds the integral of Ψ n along the segment from p to q to acc.
fn segment(p: f64[2], q: f64[2]) -> bool {
    for s in GAUSS_2 {
        if !shape_at(p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])) { return false; }
        for j, i in near {
            if acc[2 * j] == 0.0 && acc[2 * j + 1] == 0.0 { touched.push(j); }
            // n ds = (dy, -dx) dt on a side of a counterclockwise domain
            acc[2 * j] += 0.5 * psi[i] * (q[1] - p[1]);
            acc[2 * j + 1] -= 0.5 * psi[i] * (q[0] - p[0]);
        }
    }
    return true;
}

fn mid(a: i64, b: i64) -> f64[2] { return [0.5 * (px(a) + px(b)), 0.5 * (py(a) + py(b))]; }

// The SMOOTH record of the part of the nodal domain of corner k of
// triangle t: the quadrilateral of the corner, the midpoints of its two
// edges and the centroid. It is record 3 t + k.
fn smooth_record(t: i64, k: i64) -> bool {
    touched.clear();
    let p = tri[3 * t + k];
    let a = tri[3 * t + (k + 1) % 3];
    let b = tri[3 * t + (k + 2) % 3];
    let c = [(px(p) + px(a) + px(b)) / 3.0, (py(p) + py(a) + py(b)) / 3.0];
    let o = [px(p), py(p)];
    if !(segment(o, mid(p, a)) && segment(mid(p, a), c) && segment(c, mid(p, b)) && segment(mid(p, b), o)) { return false; }
    let w = twice_area(p, a, b) / 6.0;
    open_record(SMOOTH, px(p), py(p), w, -1, -1);
    for j in touched {
        item.push(j);
        item_psi.push(0.0);
        item_dx.push(acc[2 * j] / w);
        item_dy.push(acc[2 * j + 1] / w);
        acc[2 * j] = 0.0;
        acc[2 * j + 1] = 0.0;
    }
    rec_start.push(item.len);
    return true;
}

// The SMOOTH record of the end `half` of boundary edge e: the record of
// that corner of the triangle of the edge.
fn owner_record(e: i64, half: i64) -> i64 {
    let t = edge_tri[e];
    let p = edges[2 * e + half];
    var k = 0;
    while tri[3 * t + k] != p { k++; }
    return 3 * t + k;
}

// The records of the cloud: SMOOTH for each corner of each triangle, DOMAIN for the scalar
// physics, and EDGE for each boundary edge. False if a shape function is
// singular.
fn build_records() -> bool {
    rec_start.clear();
    rec_start.push(0);
    rec_kind.clear();
    rec_x.clear();
    rec_y.clear();
    rec_w.clear();
    rec_edge.clear();
    rec_owner.clear();
    rec_g.clear();
    item.clear();
    item_psi.clear();
    item_dx.clear();
    item_dy.clear();
    acc.clear();
    for v in 2 * point_count() { acc.push(0.0); }
    for t in triangle_count() {
        for k in 3 { if !smooth_record(t, k) { return false; } }
    }
    if physics == SCALAR {
        for t in triangle_count() {
            let a = tri[3 * t];
            let b = tri[3 * t + 1];
            let c = tri[3 * t + 2];
            let w = twice_area(a, b, c) / 6.0;
            for k in 3 {
                let x = (3.0 * xy[2 * tri[3 * t + k]] + px(a) + px(b) + px(c)) / 6.0;
                let y = (3.0 * xy[2 * tri[3 * t + k] + 1] + py(a) + py(b) + py(c)) / 6.0;
                if !shape_record(DOMAIN, x, y, w, -1, -1) { return false; }
            }
        }
    }
    // The points of each half of the edge, as on the side of the nodal
    // domain of each end: the boundary integrals of the two forms are then
    // the same sum, and the patch test holds.
    for e in edge_count() {
        let n, l = edge_normal(e);
        for half in 2 {
            for s in GAUSS_2 {
                let q = edge_point(e, 0.5 * (half + s));
                if !shape_record(EDGE, q[0], q[1], 0.25 * l, e, owner_record(e, half)) { return false; }
            }
        }
    }
    return true;
}

// The matrix couples two points if one record holds both. A DOMAIN record
// couples only with advection.
fn rkpm_pattern() {
    let n = point_count();
    let advection = physics == SCALAR && (ax != 0.0 || ay != 0.0);
    var start: i64[>..] = [];
    for p in n + 1 { start.push(0); }
    for r in record_count() {
        if rec_kind[r] == DOMAIN && !advection { continue; }
        for p in item[rec_start[r]..rec_start[r + 1]] { start[p + 1]++; }
    }
    for p in n { start[p + 1] += start[p]; }
    var recs: i64[>..] = [];
    for k in start[n] { recs.push(0); }
    var at: i64[>..] = [];
    at.append(start[..n]);
    for r in record_count() {
        if rec_kind[r] == DOMAIN && !advection { continue; }
        for p in item[rec_start[r]..rec_start[r + 1]] {
            recs[at[p]] = r;
            at[p]++;
        }
    }
    var mark: i64[>..] = [];
    for p in n { mark.push(-1); }
    start_matrix();
    var cols: i64[>..] = [];
    for j in n {
        cols.clear();
        for r in recs[start[j]..start[j + 1]] {
            for k in item[rec_start[r]..rec_start[r + 1]] {
                if mark[k] != j {
                    mark[k] = j;
                    cols.push(k);
                }
            }
        }
        sort(cols);
        add_rows(cols);
    }
    size_matrix();
}

// The rkpm method on the cloud with the physics: the records and the
// pattern. False if a shape function is singular: the support is too
// small for the cloud.
fn rkpm_setup() -> bool {
    if !build_records() { return false; }
    rkpm_pattern();
    return true;
}

// Adds v to the entry of row dof a and column dof b.
fn add(a: i64, b: i64, v: f64) { val[entry(a, b)] += v; }

// The terms of the SMOOTH record r: the stiffness at the point.
fn smooth_terms(r: i64) {
    let w = rec_w[r] * thickness;
    for i in rec_start[r]..rec_start[r + 1] {
        let j = item[i];
        let gx = item_dx[i];
        let gy = item_dy[i];
        for k in rec_start[r]..rec_start[r + 1] {
            let m = item[k];
            let hx = item_dx[k];
            let hy = item_dy[k];
            if physics == SCALAR {
                add(j, m, w * diffusivity * (gx * hx + gy * hy));
                continue;
            }
            let e = entry(2 * j, 2 * m);
            let f = entry(2 * j + 1, 2 * m);
            val[e] += w * ((lame + 2.0 * shear) * gx * hx + shear * gy * hy);
            val[e + 1] += w * (lame * gx * hy + shear * gy * hx);
            val[f] += w * (lame * gy * hx + shear * gx * hy);
            val[f + 1] += w * ((lame + 2.0 * shear) * gy * hy + shear * gx * hx);
        }
    }
}

// The terms of the DOMAIN record r: the source, and the advection.
fn domain_terms<SRC>(r: i64) {
    let w = rec_w[r];
    let f = SRC(rec_x[r], rec_y[r], time);
    let advection = ax != 0.0 || ay != 0.0;
    for i in rec_start[r]..rec_start[r + 1] {
        rhs[item[i]] += w * f * item_psi[i];
        if !advection { continue; }
        for k in rec_start[r]..rec_start[r + 1] {
            add(item[i], item[k], w * item_psi[i] * (ax * item_dx[k] + ay * item_dy[k]));
        }
    }
}

// Component c of the traction σ(Ψ_j e_p) n of the entry i of a SMOOTH
// record, or d ∂Ψ_j/∂n of the scalar physics, with the smoothed gradient.
fn traction_of(i: i64, p: i64, c: i64, n: f64[2]) -> f64 {
    let g = [item_dx[i], item_dy[i]];
    let gn = g[0] * n[0] + g[1] * n[1];
    if physics == SCALAR { return diffusivity * gn; }
    return lame * n[c] * g[p] + shear * ((if c == p { gn } else { 0.0 }) + n[p] * g[c]);
}

// The terms of the EDGE record r: the Nitsche terms of each fixed
// component, the load of each free one, and the advective flux of a free
// edge of the scalar physics.
fn edge_terms<BC>(r: i64) {
    let e = rec_edge[r];
    let n, l = edge_normal(e);
    let w = rec_w[r] * thickness;
    let beta = NITSCHE * (if physics == SCALAR { diffusivity } else { lame + 2.0 * shear }) / l;
    let o = rec_owner[r];
    for c in ncomp {
        let k = edge_fix[ncomp * e + c];
        if k == FREE {
            free_terms(r, c, w, n);
            continue;
        }
        let g = BC(k, rec_x[r], rec_y[r], time);
        rec_g[ncomp * r + c] = g;
        for i in rec_start[r]..rec_start[r + 1] {
            rhs[ncomp * item[i] + c] += w * beta * g * item_psi[i];
            for k2 in rec_start[r]..rec_start[r + 1] { add(ncomp * item[i] + c, ncomp * item[k2] + c, w * beta * item_psi[i] * item_psi[k2]); }
        }
        for i in rec_start[o]..rec_start[o + 1] {
            for p in ncomp {
                let t = w * traction_of(i, p, c, n);
                rhs[ncomp * item[i] + p] -= t * g;
                for k2 in rec_start[r]..rec_start[r + 1] {
                    add(ncomp * item[i] + p, ncomp * item[k2] + c, -t * item_psi[k2]);
                    add(ncomp * item[k2] + c, ncomp * item[i] + p, -t * item_psi[k2]);
                }
            }
        }
    }
}

// The load of the free component c at the EDGE record r: the traction or
// the flux, and (a·n) u v of the scalar physics, so that the flux
// (a u - d grad(u))·n is the load.
fn free_terms(r: i64, c: i64, w: f64, n: f64[2]) {
    let e = rec_edge[r];
    let load = edge_load[ncomp * e + c];
    let an = if physics == SCALAR { ax * n[0] + ay * n[1] } else { 0.0 };
    for i in rec_start[r]..rec_start[r + 1] {
        rhs[ncomp * item[i] + c] += w * load * item_psi[i];
        if an == 0.0 { continue; }
        for k in rec_start[r]..rec_start[r + 1] { add(item[i], item[k], -w * an * item_psi[i] * item_psi[k]); }
    }
}

// The matrix and the right side of rkpm with the source SRC and the
// fixed values BC.
fn rkpm_assemble<SRC, BC>() {
    for &v in val { v = 0.0; }
    for &v in rhs { v = 0.0; }
    for r in record_count() {
        let kind = rec_kind[r];
        if kind == SMOOTH {
            smooth_terms(r);
        } else if kind == DOMAIN {
            domain_terms(r, SRC);
        } else {
            edge_terms(r, BC);
        }
    }
}

// u at the record r: component c of Σ Ψ_j d_j, and of its gradient.
fn record_value(r: i64, c: i64) -> f64, f64, f64 {
    var v = 0.0;
    var gx = 0.0;
    var gy = 0.0;
    for i in rec_start[r]..rec_start[r + 1] {
        let d = sol[ncomp * item[i] + c];
        v += item_psi[i] * d;
        gx += item_dx[i] * d;
        gy += item_dy[i] * d;
    }
    return v, gx, gy;
}

// The reaction of the edges that `on` marks: the integral of the
// traction of each fixed component, σ(u) n - β (u - g), which the loads
// balance. Also the sum of the sizes of its parts at each Gauss point.
fn rkpm_reaction(on: u8[:]) -> f64[2], f64 {
    var out: f64[2] = [0.0, 0.0];
    var size = 0.0;
    for r in record_count() {
        if rec_kind[r] != EDGE || on[rec_edge[r]] == 0 { continue; }
        let e = rec_edge[r];
        let n, l = edge_normal(e);
        let w = rec_w[r] * thickness;
        let beta = NITSCHE * (if physics == SCALAR { diffusivity } else { lame + 2.0 * shear }) / l;
        for c in ncomp {
            if edge_fix[ncomp * e + c] == FREE { continue; }
            var t = 0.0;
            let o = rec_owner[r];
            for i in rec_start[o]..rec_start[o + 1] {
                for p in ncomp { t += traction_of(i, p, c, n) * sol[ncomp * item[i] + p]; }
            }
            let u, ux, uy = record_value(r, c);
            let part = w * (t - beta * (u - rec_g[ncomp * r + c]));
            out[c] += part;
            size += abs(part);
        }
    }
    return out, size;
}
