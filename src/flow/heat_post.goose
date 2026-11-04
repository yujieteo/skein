import std;
import math;
namespace flow;
import .mesh;
import .cell;
import .solve;
import .post;
import .heat;

fn fix_temperature(value: f64) {
    guard thermal && cell_count() > 0 && fixed.len == dof_count() && region.len == node_count() else { abort("fix_temperature: clear_bcs on a thermal mesh"); }
    for n in node_count() {
        if region[n] == 0 { continue; }
        fixed[tdof(n)] = 1;
        prescribed[tdof(n)] = value;
    }
}

fn add_heat_flux(flux: f64) {
    guard thermal && edge_load.len == dof_count() && region.len == node_count() else { abort("add_heat_flux: clear_bcs on a thermal mesh"); }
    let value = flux / capacity;
    for c in cell_count() {
        if int_part(c) != region_part { continue; }
        for e in 3 {
            let a = cell_node(c, EDGES[2 * e]);
            let b = cell_node(c, EDGES[2 * e + 1]);
            if edge_free[3 * c + e] == 0 || region[a] == 0 || region[b] == 0 { continue; }
            let len = sqrt((coords[2 * b] - coords[2 * a]) * (coords[2 * b] - coords[2 * a]) + (coords[2 * b + 1] - coords[2 * a + 1]) * (coords[2 * b + 1] - coords[2 * a + 1]));
            if order == 1 {
                edge_load[tdof(a)] += 0.5 * len * value;
                edge_load[tdof(b)] += 0.5 * len * value;
            } else {
                edge_load[tdof(a)] += len / 6.0 * value;
                edge_load[tdof(b)] += len / 6.0 * value;
                edge_load[tdof(cell_node(c, 3 + e))] += 4.0 * len / 6.0 * value;
            }
        }
    }
}

fn int_part(c: i64) -> i64 { return if is_solid_cell(c) { 1 } else { 0 }; }

struct Local {
    vx: f64,
    vy: f64,
    p: f64,
    t: f64,
    tx: f64,
    ty: f64,
    nx: f64,
    ny: f64,
}

fn local_at(c: i64, l: f64[3]) -> Local {
    guard c >= 0 && c < cell_count() && sol.len == dof_count() else { abort("local_at: a solve, and c < cell_count()"); }
    var f = Local { vx: 0.0, vy: 0.0, p: 0.0, t: 0.0, tx: 0.0, ty: 0.0, nx: 0.0, ny: 0.0 };
    let s = shape(geometry(c), l);
    for a in per_cell {
        let node = cell_node(c, a);
        if !is_solid_cell(c) {
            f.vx += s.n[a] * sol[vdof(node, 0)];
            f.vy += s.n[a] * sol[vdof(node, 1)];
        }
        if thermal {
            f.t += s.n[a] * sol[tdof(node)];
            f.tx += s.dn[2 * a] * sol[tdof(node)];
            f.ty += s.dn[2 * a + 1] * sol[tdof(node)];
        }
    }
    if !is_solid_cell(c) {
        for k in 3 { f.p += l[k] * sol[pdof(cell_node(c, k))]; }
    }
    return f;
}

fn boundary_value<F>(found: bool, mean: bool) -> f64 {
    guard sol.len == dof_count() && region.len == node_count() && found else { abort("boundary_value: a solve, and a region that holds a node of the mesh"); }
    let gauss: f64[3] = [0.5 - 0.5 * sqrt(0.6), 0.5, 0.5 + 0.5 * sqrt(0.6)];
    let weight: f64[3] = [5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0];
    var total = 0.0;
    var length = 0.0;
    for c in cell_count() {
        if int_part(c) != region_part { continue; }
        for e in 3 {
            let i = EDGES[2 * e];
            let j = EDGES[2 * e + 1];
            let a = cell_node(c, i);
            let b = cell_node(c, j);
            if edge_free[3 * c + e] == 0 || region[a] == 0 || region[b] == 0 { continue; }
            let tx = coords[2 * b] - coords[2 * a];
            let ty = coords[2 * b + 1] - coords[2 * a + 1];
            let len = sqrt(tx * tx + ty * ty);
            length += len;
            for q in 3 {
                var l: f64[3] = [0.0; 3];
                l[i] = 1.0 - gauss[q];
                l[j] = gauss[q];
                var f = local_at(c, l);
                f.nx = ty / len;
                f.ny = -tx / len;
                total += weight[q] * len * F(f);
            }
        }
    }
    guard length > 0.0 else { abort("boundary_value: the region holds an edge of the mesh"); }
    return if mean { total / length } else { total };
}

let LINE_SAMPLES = 8;

fn line_max<F>(x0: f64, y0: f64, x1: f64, y1: f64) -> f64 {
    guard sol.len == dof_count() else { abort("line_max: a solve"); }
    let tol = 1.0e-9;
    var best = -1.0e300;
    var found = false;
    for c in cell_count() {
        if is_solid_cell(c) { continue; }
        let la = barycentric(c, x0, y0);
        let lb = barycentric(c, x1, y1);
        var lo = 0.0;
        var hi = 1.0;
        for k in 3 {
            let d = lb[k] - la[k];
            if abs(d) < 1.0e-14 {
                if la[k] < -tol { hi = -1.0; }
            } else if d > 0.0 {
                lo = max(lo, (-tol - la[k]) / d);
            } else {
                hi = min(hi, (-tol - la[k]) / d);
            }
        }
        if lo > hi { continue; }
        for m in LINE_SAMPLES + 1 {
            let t = lo + (hi - lo) * m / LINE_SAMPLES;
            var l: f64[3] = [0.0; 3];
            for k in 3 { l[k] = max(0.0, la[k] + t * (lb[k] - la[k])); }
            best = max(best, F(local_at(c, l)));
            found = true;
        }
    }
    guard found else { abort("line_max: the line meets a cell of the fluid"); }
    return best;
}
