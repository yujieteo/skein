import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .post;
import .edge;

let CURVE_LINE = 0;
let CURVE_CIRCLE = 1;
let CURVE_FOUND = -1.0e-6;    
let CURVE_NEAR = -0.1;        

var loc_lo: f64[3] = [0.0; 3];
var loc_step: f64[3] = [1.0; 3];
var loc_n: i64[3] = [1; 3];
var loc_start: i64[>..] = [];
var loc_items: i64[>..] = [];
var curve_missed = 0;         

fn loc_bucket(t: i64, x: f64) -> i64 {
    return clamp(floor((x - loc_lo[t]) / loc_step[t]) as! i64, 0, loc_n[t] - 1);
}

fn locate_build() {
    var hi: f64[3] = [-1.0e300; 3];
    loc_lo = [1.0e300; 3];
    for n in node_count() {
        for t in dim {
            loc_lo[t] = min(loc_lo[t], node_x(n, t));
            hi[t] = max(hi[t], node_x(n, t));
        }
    }
    for t in 3 { if t >= dim { loc_lo[t] = 0.0; hi[t] = 0.0; } }
    let per_axis: f64 = if dim == 2 { sqrt(cell_count() / 2.0) } else { pow(cell_count() / 2.0, 1.0 / 3.0) };
    let wanted = max(1, min(per_axis, 600.0) as! i64);
    for t in 3 {
        loc_n[t] = if t < dim { wanted } else { 1 };
        loc_step[t] = max((hi[t] - loc_lo[t]) / (loc_n[t]), 1.0e-300) * (1.0 + 1.0e-9);
    }
    let buckets = loc_n[0] * loc_n[1] * loc_n[2];
    loc_start.clear();
    for b in buckets + 1 { loc_start.push(0); }
    var fill: i64[>..] = [];
    for pass in 2 {
        if pass == 1 {
            for b in buckets { loc_start[b + 1] += loc_start[b]; }
            loc_items.clear();
            for i in loc_start[buckets] { loc_items.push(0); }
            fill.clear();
            fill.append(loc_start[..buckets]);
        }
        for c in cell_count() {
            if !is_live(c) { continue; }
            var lo: i64[3] = [loc_n[0], loc_n[1], loc_n[2]];
            var up: i64[3] = [0, 0, 0];
            for t in 3 { if t >= dim { lo[t] = 0; } }
            for l in dim + 1 {
                let n = cell_node(c, l);
                for t in dim {
                    lo[t] = min(lo[t], loc_bucket(t, node_x(n, t)));
                    up[t] = max(up[t], loc_bucket(t, node_x(n, t)));
                }
            }
            for k in lo[2]..up[2] + 1 {
                for j in lo[1]..up[1] + 1 {
                    for i in lo[0]..up[0] + 1 {
                        let b = i + loc_n[0] * (j + loc_n[1] * k);
                        if pass == 0 {
                            loc_start[b + 1]++;
                        } else {
                            loc_items[fill[b]] = c;
                            fill[b]++;
                        }
                    }
                }
            }
        }
    }
}

fn locate(p: f64[3]) -> i64, f64[4], bool {
    let b = loc_bucket(0, p[0]) + loc_n[0] * (loc_bucket(1, p[1]) + loc_n[1] * loc_bucket(2, p[2]));
    var best = -1;
    var best_l = -1.0e300;
    var best_b: f64[4] = [0.0; 4];
    for e in loc_start[b]..loc_start[b + 1] {
        let c = loc_items[e];
        let l = barycentric(c, p);
        var lowest = l[0];
        for a in 1..dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
            best_b = l;
        }
        if lowest >= CURVE_FOUND { break; }
    }
    return best, best_b, best >= 0 && best_l >= CURVE_NEAR;
}

fn wave_sample(x: f64, y: f64, z: f64) -> f64[8], bool {
    var r: f64[8] = [0.0; 8];
    guard sol.len == dof_total() && ncomp == 2 && loc_start.len > 0 else { return r, false; }
    let c, l, found = locate([x, y, z]);
    guard found else {
        curve_missed++;
        return r, false;
    }
    cell_geometry(c, [l[1], l[2], l[3]]);
    for k in per_cell {
        let n = cell_node(c, k);
        for a in 2 {
            r[a] += shp[k] * sol[2 * n + a];
            for t in dim { r[2 + 3 * a + t] += grad[3 * k + t] * sol[2 * n + a]; }
        }
    }
    return r, true;
}

fn curve_sample(x: f64, y: f64, z: f64) -> f64[12], bool {
    var r: f64[12] = [0.0; 12];
    if physics == EDGE {
        var found = true;
        for k in 2 {
            let v, f = edge_value_at([x, y, z], k == 1);
            for i in 6 { r[6 * k + i] = v[i]; }
            found = found && f;
        }
        if !found { curve_missed++; }
        return r, found;
    }
    let s, found = wave_sample(x, y, z);
    for i in 8 { r[i] = s[i]; }
    return r, found;
}

fn curve_integral<F>(shape_kind: i64, a: f64[3], b: f64[3], r: f64, segments: i64) -> f64[2], bool {
    guard segments > 0 && (shape_kind == CURVE_LINE || shape_kind == CURVE_CIRCLE) else { abort("curve_integral: a line or a circle and segments > 0"); }
    locate_build();
    let gauss: f64[3] = [0.5 - 0.5 * sqrt(0.6), 0.5, 0.5 + 0.5 * sqrt(0.6)];
    let weight: f64[3] = [5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0];
    var total: f64[2] = [0.0, 0.0];
    curve_missed = 0;
    var length = r * 2.0 * PI;
    var tangent: f64[3] = [0.0; 3];
    if shape_kind == CURVE_LINE {
        var sq = 0.0;
        for t in 3 {
            tangent[t] = b[t] - a[t];
            sq += tangent[t] * tangent[t];
        }
        length = sqrt(sq);
        for t in 3 { tangent[t] /= length; }
    }
    for j in segments {
        for k in 3 {
            let u = ((j) + gauss[k]) / (segments);
            var p: f64[3] = [0.0; 3];
            var n: f64[3] = [0.0; 3];
            if shape_kind == CURVE_LINE {
                for t in 3 { p[t] = a[t] + u * (b[t] - a[t]); }
                if dim == 2 {
                    n[0] = tangent[1];
                    n[1] = -tangent[0];
                }
            } else {
                let th = 2.0 * PI * u;
                n[0] = cos(th);
                n[1] = sin(th);
                for t in 3 { p[t] = a[t] + r * n[t]; }
            }
            let sample, found = curve_sample(p[0], p[1], p[2]);
            let f = F(p, n, sample);
            let ds = length * weight[k] / (segments);
            total[0] += f[0] * ds;
            total[1] += f[1] * ds;
        }
    }
    return total, curve_missed == 0;
}
