// The saved state of a step, and the field of another model on the nodes
// of this one: `results` and the cut of a submodel (spec §15.1, §15.2).
//
// save_state writes the mesh, the material and the solution of a step in
// a directory. load_source reads them as a second mesh, the source, beside
// the mesh of this model. source_value evaluates the source field at a
// point with the shape functions of the source cells. A grid of buckets
// over the bounding box of the source finds the cells near a point: each
// bucket lists the cells whose box of vertices meets it.
import std;
import os;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .post;

let MAX_BUCKETS = 2000000;
let INSIDE = 1.0e-6;    // the least barycentric coordinate of a point in a cell

var src_dim = 0;
var src_order = 1;
var src_per_cell = 0;
var src_ncomp = 0;
var src_lame = 0.0;
var src_shear = 0.0;
var src_coords: f64[>..] = [];
var src_conn: i64[>..] = [];
var src_sol: f64[>..] = [];
var src_lo: f64[3] = [0.0, 0.0, 0.0];
var src_step: f64[3] = [1.0, 1.0, 1.0];
var src_cells: i64[3] = [1, 1, 1];
var bucket_start: i64[>..] = [];
var bucket_items: i64[>..] = [];

// Writes the mesh, the material and the solution of step `step` in the
// directory `dir`. False if it cannot write a file.
fn save_state(dir: const u8[:], step: const u8[:]) -> bool {
    guard sol.len == dof_total() else { abort("save_state: a solution of each dof"); }
    let sizes: i64[>..] = [dim, order, per_cell, ncomp];
    let law: f64[>..] = [lame, shear];
    return save_array(dir, str(step, ".state_sizes"), sizes) && save_array(dir, str(step, ".state_law"), law) &&
           save_array(dir, str(step, ".state_coords"), coords) && save_array(dir, str(step, ".state_conn"), conn) &&
           save_array(dir, str(step, ".state_sol"), sol);
}

fn bucket_of(t: i64, x: f64) -> i64 {
    return clamp(floor((x - src_lo[t]) / src_step[t]) as! i64, 0, src_cells[t] - 1);
}

fn bucket_at(p: f64[3]) -> i64 {
    return bucket_of(0, p[0]) + src_cells[0] * (bucket_of(1, p[1]) + src_cells[1] * bucket_of(2, p[2]));
}

// Puts each source cell in the buckets that its box of vertices meets.
fn fill_buckets() {
    let cells = src_conn.len / src_per_cell;
    let buckets = src_cells[0] * src_cells[1] * src_cells[2];
    bucket_start.clear();
    for b in buckets + 1 { bucket_start.push(0); }
    var fill: i64[>..] = [];
    // Pass 0 counts the cells of each bucket, and pass 1 puts them in.
    for pass in 2 {
        if pass == 1 {
            for b in buckets { bucket_start[b + 1] += bucket_start[b]; }
            bucket_items.clear();
            for i in bucket_start[buckets] { bucket_items.push(0); }
            fill.append(bucket_start[..buckets]);
        }
        for c in cells {
            var lo: i64[3] = [src_cells[0], src_cells[1], src_cells[2]];
            var up: i64[3] = [0, 0, 0];
            for l in src_dim + 1 {
                let n = src_conn[src_per_cell * c + l];
                for t in 3 {
                    lo[t] = min(lo[t], bucket_of(t, src_coords[3 * n + t]));
                    up[t] = max(up[t], bucket_of(t, src_coords[3 * n + t]));
                }
            }
            for k in lo[2]..up[2] + 1 {
                for j in lo[1]..up[1] + 1 {
                    for i in lo[0]..up[0] + 1 {
                        let b = i + src_cells[0] * (j + src_cells[1] * k);
                        if pass == 0 {
                            bucket_start[b + 1]++;
                        } else {
                            bucket_items[fill[b]] = c;
                            fill[b]++;
                        }
                    }
                }
            }
        }
    }
}

// Reads the state that save_state wrote for step `step` in `dir` as the
// source. Each length and each value of the solution is multiplied by
// `length_scale` and `value_scale`, and each modulus by `modulus_scale`,
// into the units of this model. False if it cannot read the state.
fn load_source(dir: const u8[:], step: const u8[:], length_scale: f64, value_scale: f64, modulus_scale: f64) -> bool {
    var sizes: i64[>..] = [];
    var law: f64[>..] = [];
    guard load_array(dir, str(step, ".state_sizes"), sizes) && sizes.len == 4 &&
          load_array(dir, str(step, ".state_law"), law) && law.len == 2 &&
          load_array(dir, str(step, ".state_coords"), src_coords) && load_array(dir, str(step, ".state_conn"), src_conn) &&
          load_array(dir, str(step, ".state_sol"), src_sol) else { return false; };
    src_dim = sizes[0];
    src_order = sizes[1];
    src_per_cell = sizes[2];
    src_ncomp = sizes[3];
    guard (src_dim == 2 || src_dim == 3) && src_per_cell == node_count_of(src_dim, src_order) && src_conn.len % src_per_cell == 0 &&
          src_conn.len > 0 && src_sol.len == src_ncomp * (src_coords.len / 3) else { return false; };
    src_lame = modulus_scale * law[0];
    src_shear = modulus_scale * law[1];
    for &v in src_coords { v *= length_scale; }
    for &v in src_sol { v *= value_scale; }
    var hi: f64[3] = [0.0, 0.0, 0.0];
    for t in 3 {
        src_lo[t] = 1.0e300;
        hi[t] = -1.0e300;
        for n in src_coords.len / 3 {
            src_lo[t] = min(src_lo[t], src_coords[3 * n + t]);
            hi[t] = max(hi[t], src_coords[3 * n + t]);
        }
    }
    // About one bucket for each cell.
    var volume = 1.0;
    for t in src_dim { volume *= max(hi[t] - src_lo[t], 1.0e-300); }
    let count = min(src_conn.len / src_per_cell, MAX_BUCKETS) as f64;
    let side = pow(volume / count, 1.0 / src_dim);
    for t in 3 {
        src_cells[t] = if t < src_dim { max(ceil((hi[t] - src_lo[t]) / side) as! i64, 1) } else { 1 };
        src_step[t] = if t < src_dim { max((hi[t] - src_lo[t]) / src_cells[t], 1.0e-300) } else { 1.0 };
    }
    fill_buckets();
    return true;
}

fn source_node(c: i64, l: i64) -> i64 { return src_conn[src_per_cell * c + l]; }

// The barycentric coordinates of p in source cell c, from its vertices.
fn source_barycentric(c: i64, p: f64[3]) -> f64[4] {
    var e: f64[9] = [0.0; 9];
    var b: f64[3] = [0.0; 3];
    let v0 = source_node(c, 0);
    for r in src_dim {
        let v = source_node(c, r + 1);
        for t in src_dim { e[3 * t + r] = src_coords[3 * v + t] - src_coords[3 * v0 + t]; }
        b[r] = p[r] - src_coords[3 * v0 + r];
    }
    if src_dim == 2 { e[8] = 1.0; }
    let det = e[0] * (e[4] * e[8] - e[5] * e[7]) - e[1] * (e[3] * e[8] - e[5] * e[6]) + e[2] * (e[3] * e[7] - e[4] * e[6]);
    guard det != 0.0 else { return [-1.0, -1.0, -1.0, -1.0]; }
    var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
    // Cramer's rule: column r of e replaced by b.
    for r in src_dim {
        var m = e;
        for t in 3 { m[3 * t + r] = b[t]; }
        let xi = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
        l[r + 1] = xi;
        l[0] -= xi;
    }
    return l;
}

// The source cell that holds p, its reference point, and whether p is in
// it. A point a little outside the source mesh takes the cell of its
// bucket that is nearest in barycentric terms.
fn source_cell(p: f64[3]) -> i64, f64[3], bool {
    guard bucket_start.len > 1 else { abort("source_cell: load_source first"); }
    for t in src_dim {
        if p[t] < src_lo[t] - src_step[t] || p[t] > src_lo[t] + (src_cells[t] + 1) * src_step[t] { return -1, [0.0; 3], false; }
    }
    let b = bucket_at(p);
    var best = -1;
    var best_l = -1.0e300;
    var best_xi: f64[3] = [0.0; 3];
    for c in bucket_items[bucket_start[b]..bucket_start[b + 1]] {
        let l = source_barycentric(c, p);
        var lowest = l[0];
        for a in 1..src_dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
            best_xi = [l[1], l[2], l[3]];
        }
    }
    return best, best_xi, best >= 0 && best_l >= -INSIDE;
}

// Component `comp` of the source field at p, and whether a source cell
// holds p.
fn source_value(comp: i64, p: f64[3]) -> f64, bool {
    guard comp >= 0 && comp < src_ncomp else { abort("source_value: a component of the source"); }
    let c, xi, inside = source_cell(p);
    guard c >= 0 else { return 0.0, false; }
    var n: f64[10] = [0.0; 10];
    var dn: f64[30] = [0.0; 30];
    shape(src_dim, src_order, xi, n, dn);
    var v = 0.0;
    for l in src_per_cell { v += n[l] * src_sol[src_ncomp * source_node(c, l) + comp]; }
    return v, inside;
}

// The stress of the elastic source at p: sigma[3 i + j], and whether a
// source cell holds p.
fn source_stress(p: f64[3]) -> f64[9], bool {
    guard src_ncomp == src_dim else { abort("source_stress: an elastic source"); }
    var s: f64[9] = [0.0; 9];
    let c, xi, inside = source_cell(p);
    guard c >= 0 else { return s, false; }
    var n: f64[10] = [0.0; 10];
    var dn: f64[30] = [0.0; 30];
    shape(src_dim, src_order, xi, n, dn);
    // jac[3 r + t] = dx_t / dxi_r, and du[3 a + r] = du_a / dxi_r.
    var jac: f64[9] = [0.0; 9];
    var du: f64[9] = [0.0; 9];
    for l in src_per_cell {
        let m = source_node(c, l);
        for r in src_dim {
            for t in src_dim { jac[3 * r + t] += dn[3 * l + r] * src_coords[3 * m + t]; }
            for a in src_dim { du[3 * a + r] += dn[3 * l + r] * src_sol[src_ncomp * m + a]; }
        }
    }
    if src_dim == 2 { jac[8] = 1.0; }
    let inv = inverse_3(jac);
    // dx[3 a + t] = du_a / dx_t = Σr du_a/dxi_r dxi_r/dx_t
    var dx: f64[9] = [0.0; 9];
    for a in src_dim {
        for t in src_dim {
            for r in src_dim { dx[3 * a + t] += du[3 * a + r] * inv[3 * t + r]; }
        }
    }
    var tr = 0.0;
    for a in src_dim { tr += dx[4 * a]; }
    for i in src_dim {
        for j in src_dim {
            s[3 * i + j] = src_shear * (dx[3 * i + j] + dx[3 * j + i]);
            if i == j { s[3 * i + j] += src_lame * tr; }
        }
    }
    return s, inside;
}

// The inverse of the 3×3 matrix m, with inv[3 t + r] = (m⁻¹)[t][r] for m
// stored as m[3 r + t].
fn inverse_3(m: f64[9]) -> f64[9] {
    let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    return [(m[4] * m[8] - m[5] * m[7]) / det, (m[2] * m[7] - m[1] * m[8]) / det, (m[1] * m[5] - m[2] * m[4]) / det,
            (m[5] * m[6] - m[3] * m[8]) / det, (m[0] * m[8] - m[2] * m[6]) / det, (m[2] * m[3] - m[0] * m[5]) / det,
            (m[3] * m[7] - m[4] * m[6]) / det, (m[1] * m[6] - m[0] * m[7]) / det, (m[0] * m[4] - m[1] * m[3]) / det];
}

// Holds each component of the unknown at each node of the region at the
// source value. It gives the number of nodes outside the source mesh.
fn fix_region_source() -> i64 {
    guard physics >= 0 && src_ncomp == ncomp && dofs_sized() else { abort("fix_region_source: a physics with the components of the source"); }
    var outside = 0;
    for n in node_count() {
        if !on_region(n) { continue; }
        let p: f64[3] = [node_x(n, 0), node_x(n, 1), node_x(n, 2)];
        for a in ncomp {
            let v, inside = source_value(a, p);
            if !inside && a == 0 { outside++; }
            fixed[ncomp * n + a] = 1;
            prescribed[ncomp * n + a] = v;
        }
    }
    return outside;
}

// The resultant force ∫ sigma n dA on part k of the region, with the stress
// of this model and of the source, over the facets that have all vertices
// on the part. n is the outward normal.
fn part_resultants(k: i64) -> f64[3], f64[3] {
    guard physics == ELASTIC && src_ncomp == src_dim && k >= 0 && k < region.len else { abort("part_resultants: an elastic model and source, and a part"); }
    var own: f64[3] = [0.0; 3];
    var other: f64[3] = [0.0; 3];
    let fd = dim - 1;
    let fnodes = node_count_of(fd, order);
    var fs: f64[10] = [0.0; 10];
    var fdn: f64[30] = [0.0; 30];
    for c in cell_count() {
        for f in facet_count(dim) {
            var all = true;
            for l in fd + 1 { if !on_part(cell_node(c, facet_node(dim, f, l)), k) { all = false; } }
            if !all { continue; }
            var centre: f64[3] = [0.0; 3];
            for l in dim + 1 {
                for t in 3 { centre[t] += node_x(cell_node(c, l), t) / (dim + 1); }
            }
            for qp in quad_count(fd) {
                let fxi, w = quad_point(fd, qp);
                shape(fd, order, fxi, fs, fdn);
                var x: f64[3] = [0.0; 3];
                var tan: f64[6] = [0.0; 6];
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for t in 3 {
                        x[t] += fs[l] * node_x(n, t);
                        for r in fd { tan[3 * r + t] += fdn[3 * l + r] * node_x(n, t); }
                    }
                }
                // The area normal: its size is the area of the facet per
                // unit area of the reference facet.
                var an: f64[3] = [tan[1], -tan[0], 0.0];
                if fd == 2 { an = [tan[1] * tan[5] - tan[2] * tan[4], tan[2] * tan[3] - tan[0] * tan[5], tan[0] * tan[4] - tan[1] * tan[3]]; }
                var out = 0.0;
                for t in 3 { out += an[t] * (x[t] - centre[t]); }
                if out < 0.0 { for t in 3 { an[t] = -an[t]; } }
                let scale = w * (if fd == 1 { thickness } else { 1.0 });
                // The reference point of x in the cell: the vertices of the
                // facet are vertices of the cell.
                var xi: f64[3] = [0.0; 3];
                var lf: f64[3] = [1.0, fxi[0], fxi[1]];
                for r in fd { lf[0] -= fxi[r]; }
                for l in fd + 1 {
                    let v = facet_node(dim, f, l);
                    if v > 0 { xi[v - 1] += lf[l]; }
                }
                let sg, found = source_stress(x);
                for i in dim {
                    for j in dim {
                        own[i] += cell_field(c, xi, FIELD_STRESS, 3 * i + j) * an[j] * scale;
                        if found { other[i] += sg[3 * i + j] * an[j] * scale; }
                    }
                }
            }
        }
    }
    return own, other;
}
