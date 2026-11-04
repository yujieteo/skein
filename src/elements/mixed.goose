import std;
import math;
import dictionary;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .solve;
import .post;
import .eigen;
import .edge;
import .wave;
import solvers from "../solvers";

let MIXED_TOLERANCE = 1.0e-10;
let INNER_TOLERANCE = 1.0e-13;

var facet_cells: i64[>..] = [];
var cell_facet: i64[>..] = [];
var cell_fsign: f64[>..] = [];
var mixed_k: f64[>..] = [];
var mixed_f: f64[>..] = [];
var on_pressure: u8[>..] = [];

fn facet_total() -> i64 { return facet_cells.len / 2; }

fn set_mixed(k: f64, f: f64, t: f64) {
    guard order == 1 && per_cell == dim + 1 && k > 0.0 && t > 0.0 else { abort("set_mixed: a P1 mesh, k > 0 and t > 0"); }
    physics = MIXED;
    ncomp = 1;
    thickness = if dim == 3 { 1.0 } else { t };
    clear_stages();
    build_facets();
    other_dofs = facet_total() + cell_count();
    fixed.clear();
    prescribed.clear();
    load.clear();
    sol.clear();
    for d in other_dofs {
        fixed.push(0);
        prescribed.push(0.0);
        load.push(0.0);
        sol.push(0.0);
    }
    clear_matrix();
    mixed_k.clear();
    mixed_f.clear();
    for c in cell_count() {
        mixed_k.push(k);
        mixed_f.push(f);
    }
    clear_pressure();
}

fn clear_pressure() {
    on_pressure.clear();
    for k in facet_total() { on_pressure.push(0); }
}

fn build_facets() {
    facet_cells.clear();
    cell_facet.clear();
    cell_fsign.clear();
    let nv = node_count();
    guard nv < 2000000 else { abort("build_facets: fewer than 2e6 nodes"); }
    var at = dictionary<i64, i64> {};
    for c in cell_count() {
        for f in facet_count(dim) {
            var v: i64[3] = [0, 0, 0];
            for l in dim { v[l] = cell_node(c, facet_node(dim, f, l)); }
            sort(v[..dim]);
            let key = if dim == 2 { v[0] * nv + v[1] } else { (v[0] * nv + v[1]) * nv + v[2] };
            let k = at.get_or_insert(key, facet_total());
            if k == facet_total() {
                facet_cells.push(c);
                facet_cells.push(-1);
            } else {
                facet_cells[2 * k + 1] = c;
            }
            cell_facet.push(k);
            cell_fsign.push(if facet_cells[2 * k] == c { 1.0 } else { -1.0 });
        }
    }
}

fn mixed_material(g: i64, k: f64, f: f64) -> i64 {
    guard physics == MIXED && g >= 0 && k > 0.0 else { abort("mixed_material: a mixed physics, a shape and k > 0"); }
    var count = 0;
    for c in cell_count() {
        if !cell_in(c, g) { continue; }
        mixed_k[c] = k;
        mixed_f[c] = f;
        count++;
    }
    return count;
}

fn mixed_source<F>() {
    guard physics == MIXED else { abort("mixed_source: a mixed physics"); }
    for c in cell_count() {
        var s = 0.0;
        var total = 0.0;
        for q in rich_count() {
            let xi, w = rich_point(q);
            let v: f64 = F(point_of(c, xi));
            s += w * v;
            total += w;
        }
        mixed_f[c] = s / total;
    }
}

fn opposite(f: i64) -> i64 {
    var on: i64[4] = [0, 0, 0, 0];
    for l in dim { on[facet_node(dim, f, l)] = 1; }
    for v in dim + 1 { if on[v] == 0 { return v; } }
    abort("opposite: a facet of d vertices");
}

fn cell_size(c: i64) -> f64 { return cell_geometry(c, [0.0; 3]) * (if dim == 2 { 0.5 } else { 1.0 / 6.0 }); }

fn rt_basis(c: i64, f: i64, p: f64[3], size: f64) -> f64[3] {
    let o = cell_node(c, opposite(f));
    let s = cell_fsign[(dim + 1) * c + f] / (dim * size);
    var v: f64[3] = [0.0; 3];
    for t in dim { v[t] = s * (p[t] - node_x(o, t)); }
    return v;
}

fn mixed_cell(c: i64, me: f64[16]&) {
    for &v in me { v = 0.0; }
    let size = cell_size(c);
    let nf = dim + 1;
    for q in quad_count(dim) {
        let xi, w = quad_point(dim, q);
        let p = point_of(c, xi);
        let dv = w * size * (if dim == 2 { 2.0 } else { 6.0 });
        var phi: f64[12] = [0.0; 12];
        for a in nf {
            let v = rt_basis(c, a, p, size);
            for t in 3 { phi[3 * a + t] = v[t]; }
        }
        for a in nf {
            for b in nf {
                var m = 0.0;
                for t in 3 { m += phi[3 * a + t] * phi[3 * b + t]; }
                me[nf * a + b] += m * dv / mixed_k[c];
            }
        }
    }
}

fn point_of(c: i64, xi: f64[3]) -> f64[3] {
    var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
    for r in dim {
        l[r + 1] = xi[r];
        l[0] -= xi[r];
    }
    var p: f64[3] = [0.0; 3];
    for v in dim + 1 {
        for t in 3 { p[t] += l[v] * node_x(cell_node(c, v), t); }
    }
    return p;
}

fn boundary_facet(k: i64) -> bool { return facet_cells[2 * k + 1] < 0; }

fn region_facets() -> i64[>..] {
    var out: i64[>..] = [];
    for c in cell_count() {
        for f in facet_count(dim) {
            let k = cell_facet[(dim + 1) * c + f];
            if boundary_facet(k) && facet_on_region(c, f) { out.push(k); }
        }
    }
    return out;
}

fn mixed_pressure(v: f64) -> i64 {
    guard physics == MIXED && dofs_sized() else { abort("mixed_pressure: a mixed physics and the arrays of the dofs"); }
    let ks = region_facets();
    for k in ks {
        load[k] -= v;
        on_pressure[k] = 1;
    }
    return ks.len;
}

fn mixed_normal_flux(g: f64) -> i64 {
    guard physics == MIXED && dofs_sized() else { abort("mixed_normal_flux: a mixed physics and the arrays of the dofs"); }
    var count = 0;
    for c in cell_count() {
        for f in facet_count(dim) {
            let k = cell_facet[(dim + 1) * c + f];
            if !boundary_facet(k) || !facet_on_region(c, f) { continue; }
            let n, area = facet_normal(c, f);
            fixed[k] = 1;
            prescribed[k] = g * area / thickness;
            count++;
        }
    }
    return count;
}

fn mixed_flux() -> f64 {
    guard physics == MIXED && sol.len == dof_total() else { abort("mixed_flux: a mixed physics and a solution"); }
    var s = 0.0;
    for k in region_facets() { s += sol[k]; }
    return s * thickness;
}

fn mixed_assemble(g: f64[>..]&) {
    clear_matrix();
    let nf = facet_total();
    let lf = dim + 1;
    var keys: i64[>..] = [];
    for c in cell_count() {
        for a in lf {
            for b in lf { keys.push(cell_facet[lf * c + a] * nf + cell_facet[lf * c + b]); }
        }
    }
    sort(keys);
    row_start.push(0);
    var prev = -1;
    for e in keys {
        if e == prev { continue; }
        prev = e;
        while row_start.len - 1 < e / nf { row_start.push(col.len); }
        col.push((e % nf) as i32);
        val.push(0.0);
    }
    while row_start.len < nf + 1 { row_start.push(col.len); }
    solvers::memory_check("simplex mixed cg", col.len, dof_total(), 9, 0, "");
    g.clear();
    for k in nf { g.push(if fixed[k] != 0 { prescribed[k] } else { load[k] }); }
    var me: f64[16] = [0.0; 16];
    for c in cell_count() {
        mixed_cell(c, me);
        for a in lf {
            let i = cell_facet[lf * c + a];
            for b in lf {
                let j = cell_facet[lf * c + b];
                if fixed[i] != 0 { continue; }
                if fixed[j] != 0 {
                    g[i] -= me[lf * a + b] * prescribed[j];
                } else {
                    val[entry(i, j)] += me[lf * a + b];
                }
            }
        }
    }
    for k in nf { if fixed[k] != 0 { val[entry(k, k)] = 1.0; } }
}

fn flux_multiply(x: f64[:], y: f64[:]) {
    for r in facet_total() {
        var s = 0.0;
        for e in row_start[r]..row_start[r + 1] { s += val[e] * x[col[e] as i64]; }
        y[r] = s;
    }
}

fn flux_solve(y: f64[:], x: f64[:]) -> i64 {
    let n = facet_total();
    var r: f64[>..] = [];
    var z: f64[>..] = [];
    var p: f64[>..] = [];
    var q: f64[>..] = [];
    var dinv: f64[>..] = [];
    for k in n {
        dinv.push(1.0 / val[entry(k, k)]);
        x[k] = 0.0;
        r.push(y[k]);
        z.push(y[k] * dinv[k]);
        q.push(0.0);
    }
    p.append(z);
    var rz = solvers::dot(r, z);
    let target = INNER_TOLERANCE * solvers::norm(y);
    var it = 0;
    while solvers::norm(r) > target {
        if it >= 10 * n + 100 { return -1; }
        flux_multiply(p, q);
        let alpha = rz / solvers::dot(p, q);
        for k in n {
            x[k] += alpha * p[k];
            r[k] -= alpha * q[k];
            z[k] = r[k] * dinv[k];
        }
        let rz_new = solvers::dot(r, z);
        for k in n { p[k] = z[k] + rz_new / rz * p[k]; }
        rz = rz_new;
        it++;
    }
    return it;
}

fn divergence_t(p: f64[:], y: f64[:]) {
    for &v in y { v = 0.0; }
    let lf = dim + 1;
    for c in cell_count() {
        for a in lf {
            let k = cell_facet[lf * c + a];
            if fixed[k] == 0 { y[k] -= cell_fsign[lf * c + a] * p[c]; }
        }
    }
}

fn divergence(x: f64[:], y: f64[:]) {
    let lf = dim + 1;
    for c in cell_count() {
        var s = 0.0;
        for a in lf {
            let k = cell_facet[lf * c + a];
            if fixed[k] == 0 { s -= cell_fsign[lf * c + a] * x[k]; }
        }
        y[c] = s;
    }
}

fn schur_multiply(x: f64[:], y: f64[:], t: f64[:], u: f64[:]) -> i64 {
    divergence_t(x, t);
    let its = flux_solve(t, u);
    divergence(u, y);
    return its;
}

fn mixed_solve() -> bool {
    guard physics == MIXED && dofs_sized() && mixed_k.len == cell_count() else { abort("mixed_solve: a mixed physics and the arrays of the dofs"); }
    let start = clock();
    for k in facet_total() {
        if boundary_facet(k) && fixed[k] == 0 && on_pressure[k] == 0 { fixed[k] = 1; }
    }
    var g: f64[>..] = [];
    mixed_assemble(g);
    let sizes = wave_sizes();
    let nf = facet_total();
    let nc = cell_count();
    var t = zeros(nf);
    var u = zeros(nf);
    for k in nf { t[k] = if fixed[k] != 0 { 0.0 } else { g[k] }; }
    var its = flux_solve(t, u);
    var b = zeros(nc);
    divergence(u, b);
    var dg = zeros(nc);
    let lf = dim + 1;
    for c in nc {
        var h = -mixed_f[c] * cell_size(c);
        for a in lf {
            let k = cell_facet[lf * c + a];
            if fixed[k] != 0 { h += cell_fsign[lf * c + a] * prescribed[k]; } else { dg[c] += 1.0 / val[entry(k, k)]; }
        }
        b[c] -= h;
    }
    var p = zeros(nc);
    let outer = schur_cg(b, dg, p, t, u);
    var converged = its >= 0 && outer >= 0;
    divergence_t(p, t);
    for k in nf { t[k] = if fixed[k] != 0 { 0.0 } else { g[k] - t[k] }; }
    if flux_solve(t, u) < 0 { converged = false; }
    for k in nf { sol[k] = if fixed[k] != 0 { prescribed[k] } else { u[k] }; }
    for c in nc { sol[nf + c] = p[c]; }
    guard wave_sizes() == sizes else { abort("mixed_solve: a solve grows no array of the state (rule 3)"); }
    last_iterations = max(outer, 0);
    last_seconds = clock() - start;
    return converged;
}

fn schur_cg(b: f64[:], dg: f64[:], p: f64[:], t: f64[:], u: f64[:]) -> i64 {
    let n = b.len;
    var r = zeros(n);
    var z = zeros(n);
    var d = zeros(n);
    var q = zeros(n);
    for c in n {
        p[c] = 0.0;
        r[c] = b[c];
        z[c] = r[c] / max(dg[c], 1.0e-300);
        d[c] = z[c];
    }
    var rz = solvers::dot(r, z);
    let target = MIXED_TOLERANCE * solvers::norm(b);
    var it = 0;
    while solvers::norm(r) > target {
        if it >= 10 * n + 100 || schur_multiply(d, q, t, u) < 0 { return -1; }
        let alpha = rz / solvers::dot(d, q);
        for c in n {
            p[c] += alpha * d[c];
            r[c] -= alpha * q[c];
            z[c] = r[c] / max(dg[c], 1.0e-300);
        }
        let rz_new = solvers::dot(r, z);
        for c in n { d[c] = z[c] + rz_new / rz * d[c]; }
        rz = rz_new;
        it++;
    }
    return it;
}

fn mixed_cell_value(c: i64, x: f64[3]) -> f64[4] {
    let size = cell_size(c);
    var v: f64[4] = [0.0, 0.0, 0.0, sol[facet_total() + c]];
    for f in dim + 1 {
        let phi = rt_basis(c, f, x, size);
        for t in 3 { v[t] += sol[cell_facet[(dim + 1) * c + f]] * phi[t]; }
    }
    return v;
}

fn mixed_value_at(x: f64[3]) -> f64[4], bool {
    guard physics == MIXED && sol.len == dof_total() else { abort("mixed_value_at: a mixed physics and a solution"); }
    var sum: f64[4] = [0.0; 4];
    var count = 0;
    var best = -1;
    var best_l = -1.0e300;
    for c in cell_count() {
        let l = barycentric(c, x);
        var lowest = l[0];
        for a in 1..dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
        }
        if lowest < -1.0e-9 { continue; }
        let v = mixed_cell_value(c, x);
        for t in 4 { sum[t] += v[t]; }
        count++;
    }
    if count > 0 {
        for t in 4 { sum[t] /= count; }
        return sum, true;
    }
    guard best >= 0 else { return sum, false; };
    return mixed_cell_value(best, x), best_l >= -0.1;
}

fn write_mixed_vtu(path: const u8[:], p_name: const u8[:], q_name: const u8[:]) -> bool {
    guard physics == MIXED && sol.len == dof_total() else { abort("write_mixed_vtu: a mixed physics and a solution"); }
    var v: f64[>..] = [];
    var count: f64[>..] = [];
    for n in node_count() {
        v.append([0.0, 0.0, 0.0, 0.0]);
        count.push(0.0);
    }
    for c in cell_count() {
        for a in dim + 1 {
            let n = cell_node(c, a);
            let x = mixed_cell_value(c, [node_x(n, 0), node_x(n, 1), node_x(n, 2)]);
            for t in 4 { v[4 * n + t] += x[t]; }
            count[n] += 1.0;
        }
    }
    var out: u8[>..] = [];
    vtu_mesh(out);
    format(out, "<DataArray type=\"Float64\" Name=\"", p_name, "\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for n in node_count() { format(out, v[4 * n + 3] / max(count[n], 1.0), "\n"); }
    format(out, "</DataArray>\n<DataArray type=\"Float64\" Name=\"", q_name, "\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in node_count() {
        let w = max(count[n], 1.0);
        format(out, v[4 * n] / w, " ", v[4 * n + 1] / w, " ", v[4 * n + 2] / w, "\n");
    }
    out.append("</DataArray>\n");
    out.append(VTU_END);
    return write_file(path, out);
}
