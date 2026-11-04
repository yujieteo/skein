import std;
import math;
namespace elements;
import .shape;
import .contact;
import .mesh;
import .kernels;
import .csg;
import .stages;
import .field;
import solvers from "../solvers";

let ELASTIC = 0;
let HEAT = 1;
let WAVE = 2;
let EDGE = 3;
let MIXED = 4;

let FIXED_BC = 1;
let FIXED_HELD = 2;

let SOLID_3D = 0;
let PLANE_STRESS = 1;
let PLANE_STRAIN = 2;

var physics = -1;
var ncomp = 0;
var plane = SOLID_3D;
var thickness = 1.0;
var lame = 0.0;
var lame_3d = 0.0;
var shear = 0.0;
var conductivity = 0.0;
var source = 0.0;
var volume_force: f64[3] = [0.0; 3];
var aniso_on = false;
var aniso_c: f64[36] = [0.0; 36];
var aniso_s0: f64[6] = [0.0; 6];
var aniso_t: f64[81] = [0.0; 81];

var curve_lo = 0.0;
var curve_step = 0.0;
var curve_k: f64[>..] = [];
var initial_value = 0.0;
var state_valid = false;
var keep_factor = false;
var factor_ready = false;

var fixed: u8[>..] = [];
var prescribed: f64[>..] = [];
var load: f64[>..] = [];
var sol: f64[>..] = [];

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];

var harmonic_omega = 0.0;
var other_dofs = 0;

fn dof_total() -> i64 { return if physics == EDGE || physics == MIXED { other_dofs } else { ncomp * node_count() } }

fn dofs_sized() -> bool {
    return fixed.len == dof_total() && prescribed.len == fixed.len && load.len == fixed.len && sol.len == fixed.len &&
           (row_start.len == 0 || (row_start.len == fixed.len + 1 && assembled_source.len == fixed.len));
}

fn clear_matrix() {
    row_start.clear();
    col.clear();
    val.clear();
}

fn start_physics(kind: i64, comps: i64) {
    guard per_cell > 0 && (kind == ELASTIC || kind == HEAT || kind == WAVE) && comps >= 1 && comps <= 3 else { abort("start_physics: a mesh, a physics and 1 to 3 components"); }
    physics = kind;
    ncomp = comps;
    aniso_on = false;
    harmonic_omega = 0.0;
    clear_stages();
    clear_field_terms();
    fixed.clear();
    prescribed.clear();
    load.clear();
    sol.clear();
    for d in dof_total() {
        fixed.push(0);
        prescribed.push(0.0);
        load.push(0.0);
        sol.push(0.0);
    }
    clear_matrix();
    clear_contact();
}

fn set_elastic(E: f64, nu: f64, assumption: i64, t: f64) {
    guard E > 0.0 && nu > -1.0 && nu < 0.5 && t > 0.0 && assumption >= SOLID_3D && assumption <= PLANE_STRAIN else { abort("set_elastic: E > 0, -1 < nu < 0.5, t > 0 and a plane assumption"); }
    start_physics(ELASTIC, dim);
    plane = if dim == 3 { SOLID_3D } else { assumption };
    thickness = if dim == 3 { 1.0 } else { t };
    shear = E / (2.0 * (1.0 + nu));
    lame_3d = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    lame = if plane == PLANE_STRESS { E * nu / (1.0 - nu * nu) } else { lame_3d };
    volume_force = [0.0; 3];
}

fn set_heat(k: f64, q: f64, t: f64) {
    guard k > 0.0 && t > 0.0 else { abort("set_heat: k > 0 and t > 0"); }
    start_physics(HEAT, 1);
    state_valid = false;
    curve_k.clear();
    thickness = if dim == 3 { 1.0 } else { t };
    conductivity = k;
    source = q;
}

fn set_conductivity_curve(lo: f64, hi: f64, v: f64[:]) {
    guard physics == HEAT && v.len >= 2 && hi > lo else { abort("set_conductivity_curve: a heat physics, 2 values and hi > lo"); }
    for k in v { guard k > 0.0 else { abort("set_conductivity_curve: k > 0 over the temperatures"); } }
    curve_lo = lo;
    curve_step = (hi - lo) / (v.len - 1);
    curve_k.clear();
    curve_k.append(v);
}

fn curve_value(t: f64) -> f64 {
    let n = curve_k.len;
    let x = (t - curve_lo) / curve_step;
    let i = min(max(floor(x) as! i64, 0), n - 2);
    let k = curve_k[i] + (curve_k[i + 1] - curve_k[i]) * (x - i);
    guard k > 0.0 else { abort("curve_value: k > 0 at the temperature, also beyond the ends of the curve"); }
    return k;
}

fn conductivity_at(c: i64) -> f64 {
    guard curve_k.len > 0 else { return conductivity; };
    var t = 0.0;
    for l in per_cell { t += shp[l] * sol[cell_node(c, l)]; }
    return curve_value(t);
}

var shp: f64[20];
var dref: f64[60];
var grad: f64[60];

fn cell_node(c: i64, l: i64) -> i64 {
    guard c >= 0 && c < cell_count() && l >= 0 && l < per_cell else { abort("cell_node: c < cell_count() and l < per_cell"); }
    return conn[per_cell * c + l];
}

fn cell_jacobian(c: i64, jac: f64[9]&) -> f64 {
    for l in per_cell {
        let n = cell_node(c, l);
        for r in dim {
            for t in dim { jac[3 * r + t] += dref[3 * l + r] * node_x(n, t); }
        }
    }
    if dim == 2 { return jac[0] * jac[4] - jac[1] * jac[3]; }
    return jac[0] * (jac[4] * jac[8] - jac[5] * jac[7])
         - jac[1] * (jac[3] * jac[8] - jac[5] * jac[6])
         + jac[2] * (jac[3] * jac[7] - jac[4] * jac[6]);
}

var jinv: f64[9] = [0.0; 9];

fn cell_geometry(c: i64, xi: f64[3]) -> f64 {
    guard c >= 0 && c < cell_count() && per_cell <= shp.len && 3 * per_cell <= grad.len else { abort("cell_geometry: c < cell_count() and shp and grad hold the nodes of a cell"); }
    shape(dim, order, xi, shp, dref);
    var jac: f64[9] = [0.0; 9];
    let det = cell_jacobian(c, jac);
    var inv: f64[9] = [0.0; 9];
    if dim == 2 {
        inv[0] = jac[4] / det;
        inv[1] = -jac[1] / det;
        inv[3] = -jac[3] / det;
        inv[4] = jac[0] / det;
    } else {
        inv[0] = (jac[4] * jac[8] - jac[5] * jac[7]) / det;
        inv[1] = (jac[2] * jac[7] - jac[1] * jac[8]) / det;
        inv[2] = (jac[1] * jac[5] - jac[2] * jac[4]) / det;
        inv[3] = (jac[5] * jac[6] - jac[3] * jac[8]) / det;
        inv[4] = (jac[0] * jac[8] - jac[2] * jac[6]) / det;
        inv[5] = (jac[2] * jac[3] - jac[0] * jac[5]) / det;
        inv[6] = (jac[3] * jac[7] - jac[4] * jac[6]) / det;
        inv[7] = (jac[1] * jac[6] - jac[0] * jac[7]) / det;
        inv[8] = (jac[0] * jac[4] - jac[1] * jac[3]) / det;
    }
    assert(det > 0.0);
    jinv = inv;
    for l in per_cell {
        for t in 3 {
            var g = 0.0;
            for r in dim { g += inv[3 * t + r] * dref[3 * l + r]; }
            grad[3 * l + t] = g;
        }
    }
    return det;
}

var ke: f64[3600];
var fe: f64[60];

fn build_pattern() {
    guard per_cell > 0 && ncomp > 0 else { abort("build_pattern: a mesh and a physics"); }
    clear_matrix();
    let nn = node_count();
    var keys: i64[>..] = [];
    for c in cell_count() {
        for a in per_cell {
            for b in per_cell { keys.push(cell_node(c, a) * nn + cell_node(c, b)); }
        }
    }
    if physics == ELASTIC { tie_keys(keys); }
    sort(keys);
    row_start.push(0);
    var first = 0;
    for a in nn {
        var last = first;
        while last < keys.len && keys[last] / nn == a { last++; }
        for p in ncomp {
            var prev = -1;
            for e in first..last {
                if keys[e] == prev { continue; }
                prev = keys[e];
                let b = keys[e] % nn;
                for q in ncomp {
                    col.push((ncomp * b + q) as i32);
                    val.push(0.0);
                }
            }
            row_start.push(col.len);
        }
        first = last;
    }
}

fn entry(r: i64, c: i64) -> i64 {
    guard r >= 0 && r + 1 < row_start.len && c >= 0 && c + 1 < row_start.len else { abort("entry: r and c are rows of the matrix"); }
    var lo = row_start[r];
    var hi = row_start[r + 1] - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if col[mid] < c { lo = mid + 1; } else { hi = mid; }
    }
    assert(col[lo] == c);
    return lo;
}

var assembled_source: f64[>..] = [];

fn assemble(path: u8[:], vectors: i64) {
    guard physics >= 0 && dofs_sized() else { abort("assemble: a physics and the arrays of the dofs"); }
    build_pattern();
    internal.clear();
    if !plan_direct(path, vectors) { solvers::memory_check(path, col.len, dof_total(), vectors, worker_mesh_bytes(), ""); }
    assembled_source.clear();
    for d in dof_total() { assembled_source.push(0.0); }
    fill_cells();
    add_films();
}

fn reassemble() {
    guard physics >= 0 && dofs_sized() && row_start.len > 0 else { abort("reassemble: a physics, the arrays of the dofs and a matrix"); }
    for e in val.len { val[e] = 0.0; }
    for d in assembled_source.len { assembled_source[d] = 0.0; }
    fill_cells();
    add_films();
}

fn fill_cells() {
    let m = per_cell * ncomp;
    let workers = min(solvers::case_threads(), (cell_count() + solvers::CELL_BLOCK - 1) / solvers::CELL_BLOCK);
    if workers < 2 {
        for c in cell_count() {
            if !is_live(c) { continue; }
            cell_terms(c);
            add_cell(c, ke[..m * m], fe[..m]);
        }
        use_law_0();
        return;
    }
    var ids: i64[>..] = [];
    for w in workers { ids.push(thread_spawn(cell_worker)); }
    solvers::gather_cells(cell_count(), workers, m * m + m, add_cells);
    for id in ids { thread_wait(id); }
}

var direct = false;
var direct_path: u8[>..] = [];
var direct_vectors = 0;

fn plan_direct(path: u8[:], vectors: i64) -> bool {
    var place: i64[>..] = [];
    for n in node_count() { place.push(n); }
    direct_path.clear();
    direct_path.append(path);
    direct_vectors = vectors;
    direct = (path == "simplex elastic cg amg" || path == "simplex cg amg") && on_lattice &&
             solvers::nested_plan(lat, place, ncomp, row_start, col) && solvers::nested_pays(vectors);
    if direct { solvers::memory_check("lattice ldlt", solvers::nested_entries(), dof_total(), vectors, worker_mesh_bytes(), ""); }
    return direct;
}

fn worker_mesh_bytes() -> i64 {
    return 8 * (coords.len + conn.len + cell_law.len + birth_u.len) + dead.len + born.len;
}

var incompatible_modes = false;

let MODES = 9;
let NODE_DOFS = 24;

fn pair_block(a: f64[3], b: f64[3], p: i64, q: i64) -> f64 {
    if aniso_on {
        var w = 0.0;
        for r in 3 {
            for t in 3 { w += a[r] * b[t] * aniso_t[27 * p + 9 * r + 3 * q + t]; }
        }
        return w;
    }
    var v = lame * a[p] * b[q] + shear * a[q] * b[p];
    if p == q { v += shear * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]); }
    return v;
}

fn condense_incompatible(c: i64) {
    guard physics == ELASTIC && dim == 3 && per_cell == 8 && ncomp == 3 else { abort("condense_incompatible: an elastic hexahedron"); }
    let det0 = cell_geometry(c, [0.5, 0.5, 0.5]);
    let inv0 = jinv;
    var kaa: f64[81] = [0.0; 81];
    var kau: f64[216] = [0.0; 216];
    for qp in quad_count(3) {
        let xi, w = quad_point(3, qp);
        let det = cell_geometry(c, xi);
        let dv = w * det * thickness;
        var gm: f64[9] = [0.0; 9];
        for m in 3 {
            let dp = -4.0 * (2.0 * xi[m] - 1.0) * det0 / det;
            for t in 3 { gm[3 * m + t] = inv0[3 * t + m] * dp; }
        }
        for m in 3 {
            let am: f64[3] = [gm[3 * m], gm[3 * m + 1], gm[3 * m + 2]];
            for n in 3 {
                let an: f64[3] = [gm[3 * n], gm[3 * n + 1], gm[3 * n + 2]];
                for q in 3 {
                    for r in 3 { kaa[9 * (3 * m + q) + 3 * n + r] += pair_block(am, an, q, r) * dv; }
                }
            }
            for i in 8 {
                let gi: f64[3] = [grad[3 * i], grad[3 * i + 1], grad[3 * i + 2]];
                for q in 3 {
                    for p in 3 { kau[NODE_DOFS * (3 * m + q) + 3 * i + p] += pair_block(am, gi, q, p) * dv; }
                }
            }
        }
    }
    var x: f64[216] = kau;
    for k in MODES {
        let piv = kaa[10 * k];
        assert(piv > 0.0);
        for r in k + 1..MODES {
            let f = kaa[9 * r + k] / piv;
            guard f != 0.0 else { continue; }
            for j in k..MODES { kaa[9 * r + j] -= f * kaa[9 * k + j]; }
            for j in NODE_DOFS { x[NODE_DOFS * r + j] -= f * x[NODE_DOFS * k + j]; }
        }
    }
    for kk in MODES {
        let k = MODES - 1 - kk;
        for j in NODE_DOFS {
            var s = x[NODE_DOFS * k + j];
            for t in k + 1..MODES { s -= kaa[9 * k + t] * x[NODE_DOFS * t + j]; }
            x[NODE_DOFS * k + j] = s / kaa[10 * k];
        }
    }
    for a in NODE_DOFS {
        for b in NODE_DOFS {
            var s = 0.0;
            for r in MODES { s += kau[NODE_DOFS * r + a] * x[NODE_DOFS * r + b]; }
            ke[NODE_DOFS * a + b] -= s;
        }
    }
}

fn aniso_cell(c: i64) {
    let m = per_cell * ncomp;
    for &v in ke { v = 0.0; }
    for &v in fe { v = 0.0; }
    for qp in quad_count(dim) {
        let xi, w = quad_point(dim, qp);
        let dv = w * cell_geometry(c, xi) * thickness;
        for i in per_cell {
            let gi: f64[3] = [grad[3 * i], grad[3 * i + 1], grad[3 * i + 2]];
            for p in ncomp {
                fe[ncomp * i + p] += volume_force[p] * shp[i] * dv;
                for a in 3 { fe[ncomp * i + p] -= aniso_s0_tensor(p, a) * gi[a] * dv; }
            }
            for j in per_cell {
                let gj: f64[3] = [grad[3 * j], grad[3 * j + 1], grad[3 * j + 2]];
                for p in ncomp {
                    for q in ncomp { ke[m * (ncomp * i + p) + ncomp * j + q] += pair_block(gi, gj, p, q) * dv; }
                }
            }
        }
    }
}

fn aniso_s0_tensor(p: i64, a: i64) -> f64 { return aniso_s0[voigt_of(p, a)]; }

fn voigt_of(i: i64, j: i64) -> i64 { return if i == j { i } else if i + j == 1 { 3 } else if i + j == 3 { 4 } else { 5 }; }

fn set_aniso(c: f64[36], s0: f64[6]) {
    aniso_on = true;
    aniso_c = c;
    aniso_s0 = s0;
    for p in 3 {
        for a in 3 {
            for q in 3 {
                for b in 3 {
                    let v = voigt_of(p, a);
                    let w = voigt_of(q, b);
                    aniso_t[27 * p + 9 * a + 3 * q + b] = c[6 * v + w] * (if q == b { 1.0 } else { 0.5 });
                }
            }
        }
    }
}

fn aniso_stress(du: f64[9], i: i64, j: i64) -> f64 {
    let v = voigt_of(i, j);
    let e: f64[6] = [du[0], du[4], du[8], 0.5 * (du[1] + du[3]), 0.5 * (du[5] + du[7]), 0.5 * (du[2] + du[6])];
    var s = aniso_s0[v];
    for w in 6 { s += aniso_c[6 * v + w] * e[w]; }
    return s;
}

fn cell_terms(c: i64) {
    let m = per_cell * ncomp;
    use_law(c);
    if aniso_on { aniso_cell(c); } else { cell_matrix(c); }
    if incompatible_modes && physics == ELASTIC && hex_cells && dim == 3 { condense_incompatible(c); }
    if ft_kind.len > 0 { add_field_terms(c); }
    if born.len > 0 && born[c] != 0 {
        for r in m {
            var s = 0.0;
            for k in m { s += ke[m * r + k] * birth_u[m * c + k]; }
            fe[r] += s;
        }
    }
}

thread_fn cell_worker() {
    let m = per_cell * ncomp;
    var data: f64[>..] = [];
    loop {
        let job = qget<solvers::CellJob>();
        guard job.lo >= 0 else { break; };
        data.clear();
        for c in job.lo..job.hi {
            if is_live(c) { cell_terms(c); }
            data.append(ke[..m * m]);
            data.append(fe[..m]);
        }
        qput(solvers::CellData { lo: job.lo, bad: -1, data: data[..] });
    }
}

fn add_cells(lo: i64, hi: i64, d: f64[:]) {
    let m = per_cell * ncomp;
    for c in lo..hi {
        if !is_live(c) { continue; }
        let at = (c - lo) * (m * m + m);
        add_cell(c, d[at..at + m * m], d[at + m * m..at + m * m + m]);
    }
}

fn add_cell(c: i64, k: f64[:], f: f64[:]) {
    let m = per_cell * ncomp;
    guard k.len == m * m && f.len == m else { abort("add_cell: the matrix and the source of the nodes of a cell"); }
    for i in per_cell {
        let ni = cell_node(c, i);
        for p in ncomp { assembled_source[ncomp * ni + p] += f[ncomp * i + p]; }
        for j in per_cell {
            let nj = cell_node(c, j);
            for p in ncomp {
                let e = entry(ncomp * ni + p, ncomp * nj);
                for q in ncomp { val[e + q] += k[m * (ncomp * i + p) + ncomp * j + q]; }
            }
        }
    }
}

struct Part {
    hole: bool,
    axis: i64,
    value: f64,
    c: f64[2],
    r: f64,
    shape: i64,
    node: i64,
    rel: i64,
    la: i64,
    lv: f64,
    lr: i64,
}

var region: Part[..16] = [];
var region_out: Part[..8] = [];
var region_body = -1;

fn on_part(n: i64, k: i64) -> bool {
    guard n >= 0 && n < node_count() && k >= 0 && k < region.len else { abort("on_part: n < node_count() and k < region.len"); }
    return part_holds(region[k], n);
}

fn part_holds(p: Part, n: i64) -> bool {
    if p.node >= 0 { return n == p.node; }
    let tol = geometry_tolerance();
    if p.rel != 0 {
        let dx = node_x(n, 0) - p.c[0];
        let dy = node_x(n, 1) - p.c[1];
        let v = if p.hole { sqrt(dx * dx + dy * dy) } else { node_x(n, p.axis) };
        let at = if p.hole { p.r } else { p.value };
        if p.rel == 2 { return v >= at - tol; }
        if p.rel == 1 { return v > at + tol; }
        if p.rel == -2 { return v <= at + tol; }
        return v < at - tol;
    }
    if p.shape >= 0 { return shape_gap(p.shape, n) <= tol; }
    let g = if region_body >= 0 { body_shapes[region_body] } else { body_shape };
    if !p.hole && p.la >= 0 {
        let v = node_x(n, p.la);
        let cut = if p.lr == 2 { v >= p.lv - tol } else if p.lr == 1 { v > p.lv + tol } else if p.lr == -2 { v <= p.lv + tol } else { v < p.lv - tol };
        if cut { return false; }
    }
    if !p.hole { return abs(node_x(n, p.axis) - p.value) <= tol && (g < 0 || shape_gap(g, n) <= tol); }
    let dx = node_x(n, 0) - p.c[0];
    let dy = node_x(n, 1) - p.c[1];
    return abs(sqrt(dx * dx + dy * dy) - p.r) <= tol;
}

fn on_region(n: i64) -> bool {
    guard n >= 0 && n < node_count() else { abort("on_region: n < node_count()"); }
    if region_body >= 0 && node_body[n] != region_body { return false; }
    for p in region_out { if part_holds(p, n) { return false; } }
    for k in region.len { if on_part(n, k) { return true; } }
    return false;
}

fn region_node_count() -> i64 {
    var count = 0;
    for n in node_count() { if on_region(n) { count++; } }
    return count;
}

fn fix_region(comp: i64, v: f64) {
    guard physics >= 0 && comp >= -1 && comp < ncomp && dofs_sized() else { abort("fix_region: a physics, -1 <= comp < ncomp and the arrays of the dofs"); }
    for n in node_count() {
        if !on_region(n) { continue; }
        for p in ncomp {
            if comp >= 0 && p != comp { continue; }
            fixed[ncomp * n + p] = FIXED_BC;
            prescribed[ncomp * n + p] = v;
        }
    }
}

fn fix_region_with<VALUE>(comp: i64) {
    guard physics >= 0 && comp >= 0 && comp < ncomp && dofs_sized() else { abort("fix_region_with: a physics, 0 <= comp < ncomp and the arrays of the dofs"); }
    for n in node_count() {
        if !on_region(n) { continue; }
        fixed[ncomp * n + comp] = FIXED_BC;
        prescribed[ncomp * n + comp] = VALUE([node_x(n, 0), node_x(n, 1), node_x(n, 2)]);
    }
}

fn add_region_load(gv: f64[3]) { add_facet_load(gv, 0.0); }

fn add_region_pressure(p: f64) { add_facet_load([0.0; 3], p); }

fn facet_normal(c: i64, x: f64[3], tan: f64[6]) -> f64[3] {
    var nv: f64[3] = if dim == 2 { [tan[1], -tan[0], 0.0] } else {
        [tan[1] * tan[5] - tan[2] * tan[4], tan[2] * tan[3] - tan[0] * tan[5], tan[0] * tan[4] - tan[1] * tan[3]] };
    let count = node_count_of(dim, order);
    var side = 0.0;
    for l in count {
        let n = cell_node(c, l);
        for t in 3 { side += nv[t] * (x[t] - node_x(n, t)) / count; }
    }
    let size = sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
    let s = if side < 0.0 { -1.0 / size } else { 1.0 / size };
    return [s * nv[0], s * nv[1], s * nv[2]];
}

fn add_facet_load(gv: f64[3], p: f64) {
    guard physics >= 0 && ncomp <= gv.len && dofs_sized() else { abort("add_region_load: a physics and the arrays of the dofs"); }
    let fd = dim - 1;
    let fnodes = node_count_of(fd, order);
    var fn_shape: f64[20] = [0.0; 20];
    var fn_dref: f64[60] = [0.0; 60];
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            for qp in quad_count(fd) {
                let xi, w = quad_point(fd, qp);
                shape(fd, order, xi, fn_shape, fn_dref);
                var tan: f64[6] = [0.0; 6];
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for r in fd {
                        for t in 3 { tan[3 * r + t] += fn_dref[3 * l + r] * node_x(n, t); }
                    }
                }
                var da = 0.0;
                if fd == 1 {
                    da = sqrt(tan[0] * tan[0] + tan[1] * tan[1]) * thickness;
                } else {
                    let nx = tan[1] * tan[5] - tan[2] * tan[4];
                    let ny = tan[2] * tan[3] - tan[0] * tan[5];
                    let nz = tan[0] * tan[4] - tan[1] * tan[3];
                    da = sqrt(nx * nx + ny * ny + nz * nz);
                }
                var tv = gv;
                if p != 0.0 {
                    var x: f64[3] = [0.0; 3];
                    for l in fnodes {
                        for t in 3 { x[t] += fn_shape[l] * node_x(cell_node(c, facet_node(dim, f, l)), t); }
                    }
                    tv = gv - p * facet_normal(c, x, tan);
                }
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for k in ncomp { load[ncomp * n + k] += tv[k] * fn_shape[l] * w * da; }
                }
            }
        }
    }
}

fn add_region_traction_with<G>() {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("add_region_traction_with: a physics and the arrays of the dofs"); }
    let fd = dim - 1;
    let fnodes = node_count_of(fd, order);
    var fn_shape: f64[20] = [0.0; 20];
    var fn_dref: f64[60] = [0.0; 60];
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            for qp in quad_count(fd) {
                let xi, w = quad_point(fd, qp);
                shape(fd, order, xi, fn_shape, fn_dref);
                var tan: f64[6] = [0.0; 6];
                var x: f64[3] = [0.0; 3];
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for t in 3 { x[t] += fn_shape[l] * node_x(n, t); }
                    for r in fd {
                        for t in 3 { tan[3 * r + t] += fn_dref[3 * l + r] * node_x(n, t); }
                    }
                }
                var da = 0.0;
                if fd == 1 {
                    da = sqrt(tan[0] * tan[0] + tan[1] * tan[1]) * thickness;
                } else {
                    let nx = tan[1] * tan[5] - tan[2] * tan[4];
                    let ny = tan[2] * tan[3] - tan[0] * tan[5];
                    let nz = tan[0] * tan[4] - tan[1] * tan[3];
                    da = sqrt(nx * nx + ny * ny + nz * nz);
                }
                let tv: f64[3] = G(x);
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for k in ncomp { load[ncomp * n + k] += tv[k] * fn_shape[l] * w * da; }
                }
            }
        }
    }
}

fn region_facet(c: i64, f: i64) -> bool {
    if region_body >= 0 && cell_body[c] != region_body { return false; }
    for k in region.len {
        var all = true;
        for l in facet_vertices(dim) { if !on_part(cell_node(c, facet_node(dim, f, l)), k) { all = false; } }
        if all { return true; }
    }
    return false;
}

fn facet_point(c: i64, f: i64, qp: i64, fn_shape: f64[20]&) -> f64 {
    let fd = dim - 1;
    var fn_dref: f64[60] = [0.0; 60];
    let xi, w = quad_point(fd, qp);
    shape(fd, order, xi, fn_shape, fn_dref);
    var tan: f64[6] = [0.0; 6];
    for l in node_count_of(fd, order) {
        let n = cell_node(c, facet_node(dim, f, l));
        for r in fd {
            for t in 3 { tan[3 * r + t] += fn_dref[3 * l + r] * node_x(n, t); }
        }
    }
    if fd == 1 { return w * sqrt(tan[0] * tan[0] + tan[1] * tan[1]) * thickness; }
    let nx = tan[1] * tan[5] - tan[2] * tan[4];
    let ny = tan[2] * tan[3] - tan[0] * tan[5];
    let nz = tan[0] * tan[4] - tan[1] * tan[3];
    return w * sqrt(nx * nx + ny * ny + nz * nz);
}

fn add_node_force(f: f64[3]) {
    guard physics >= 0 && dofs_sized() else { abort("add_node_force: a physics and the arrays of the dofs"); }
    for n in node_count() {
        if !on_region(n) { continue; }
        for p in ncomp { load[ncomp * n + p] += f[p]; }
    }
}

fn body_mean(what: i64, comp: i64) -> f64 {
    var sum = 0.0;
    var volume = 0.0;
    for c in cell_count() {
        if !is_live(c) { continue; }
        for qp in quad_count(dim) {
            let xi, w = quad_point(dim, qp);
            let dv = cell_geometry(c, xi) * w;
            sum += cell_field(c, xi, what, comp) * dv;
            volume += dv;
        }
    }
    return if volume > 0.0 { sum / volume } else { 0.0 };
}

fn box_mean(what: i64, comp: i64, lo: f64[3], hi: f64[3]) -> f64 {
    var sum = 0.0;
    var volume = 0.0;
    let tol = geometry_tolerance();
    for c in cell_count() {
        if !is_live(c) { continue; }
        for k in 4 {
            for j in 4 {
                for i in 4 {
                    let xi: f64[3] = [(i + 0.5) / 4.0, (j + 0.5) / 4.0, (k + 0.5) / 4.0];
                    let dv = cell_geometry(c, xi);
                    var x: f64[3] = [0.0; 3];
                    for l in per_cell { for t in 3 { x[t] += shp[l] * node_x(cell_node(c, l), t); } }
                    var inside = true;
                    for t in 3 { if x[t] < lo[t] - tol || x[t] > hi[t] + tol { inside = false; } }
                    if !inside { continue; }
                    sum += cell_field(c, xi, what, comp) * dv;
                    volume += dv;
                }
            }
        }
    }
    return if volume > 0.0 { sum / volume } else { 0.0 / 0.0 };
}

var film_cell: i64[>..] = [];
var film_facet: i64[>..] = [];
var film_h: f64[>..] = [];

fn add_film(h: f64, t_inf: f64) {
    guard physics == HEAT && h > 0.0 && dofs_sized() else { abort("add_film: a heat physics and h > 0"); }
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            film_cell.push(c);
            film_facet.push(f);
            film_h.push(h);
        }
    }
    add_region_load([h * t_inf, 0.0, 0.0]);
}

fn clear_films() {
    film_cell.clear();
    film_facet.clear();
    film_h.clear();
}

var amb_cell: i64[>..] = [];
var amb_facet: i64[>..] = [];
var amb_c: f64[>..] = [];
var amb_t: f64[>..] = [];

fn add_ambient(c: f64, t_inf: f64) {
    guard physics == HEAT && c > 0.0 && t_inf >= 0.0 && dofs_sized() else { abort("add_ambient: a heat physics, c > 0 and t_inf >= 0"); }
    for cell in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(cell, f) { continue; }
            amb_cell.push(cell);
            amb_facet.push(f);
            amb_c.push(c);
            amb_t.push(t_inf);
        }
    }
}

fn clear_ambient() {
    amb_cell.clear();
    amb_facet.clear();
    amb_c.clear();
    amb_t.clear();
}

fn facet_temperature(c: i64, f: i64, fn_shape: f64[20]&) -> f64 {
    var t = 0.0;
    for a in node_count_of(dim - 1, order) { t += fn_shape[a] * sol[cell_node(c, facet_node(dim, f, a))]; }
    return t;
}

fn ambient_residual(r: f64[:]) {
    var fn_shape: f64[20] = [0.0; 20];
    let fnodes = node_count_of(dim - 1, order);
    for c, i in amb_cell {
        let f = amb_facet[i];
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape);
            let t = facet_temperature(c, f, fn_shape);
            let g = amb_c[i] * (t * t * t * t - amb_t[i] * amb_t[i] * amb_t[i] * amb_t[i]);
            for a in fnodes { r[cell_node(c, facet_node(dim, f, a))] += fn_shape[a] * g * wda; }
        }
    }
}

fn ambient_tangent(m: f64[:]) {
    var fn_shape: f64[20] = [0.0; 20];
    let fnodes = node_count_of(dim - 1, order);
    for c, i in amb_cell {
        let f = amb_facet[i];
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape);
            let t = facet_temperature(c, f, fn_shape);
            let h = 4.0 * amb_c[i] * t * t * t * wda;
            for a in fnodes {
                let na = cell_node(c, facet_node(dim, f, a));
                for b in fnodes { m[entry(na, cell_node(c, facet_node(dim, f, b)))] += fn_shape[a] * fn_shape[b] * h; }
            }
        }
    }
}

fn region_ambient_heat() -> f64 {
    var fn_shape: f64[20] = [0.0; 20];
    var heat = 0.0;
    for c, i in amb_cell {
        let f = amb_facet[i];
        if !region_facet(c, f) { continue; }
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape);
            let t = facet_temperature(c, f, fn_shape);
            heat += amb_c[i] * (t * t * t * t - amb_t[i] * amb_t[i] * amb_t[i] * amb_t[i]) * wda;
        }
    }
    return heat;
}

fn film_heat() -> f64 {
    var fn_shape: f64[20] = [0.0; 20];
    var heat = 0.0;
    for c, i in film_cell {
        let f = film_facet[i];
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape) * film_h[i];
            for a in node_count_of(dim - 1, order) { heat += fn_shape[a] * sol[cell_node(c, facet_node(dim, f, a))] * wda; }
        }
    }
    return heat;
}

fn add_films() {
    var fn_shape: f64[20] = [0.0; 20];
    let fnodes = node_count_of(dim - 1, order);
    for c, i in film_cell {
        let f = film_facet[i];
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape) * film_h[i];
            for a in fnodes {
                let na = cell_node(c, facet_node(dim, f, a));
                for b in fnodes { val[entry(na, cell_node(c, facet_node(dim, f, b)))] += fn_shape[a] * fn_shape[b] * wda; }
            }
        }
    }
}

fn clear_conditions() {
    held.clear();
    held.append(sol);
    for &v in fixed { v = 0; }
    for &v in prescribed { v = 0.0; }
    for &v in load { v = 0.0; }
    for &v in sol { v = 0.0; }
}

var internal: f64[>..] = [];

fn row_force(d: i64) -> f64, f64 {
    if internal.len == sol.len { return internal[d], abs(internal[d]); }
    var ku = 0.0;
    var terms = 0.0;
    for e in row_start[d]..row_start[d + 1] {
        ku += val[e] * sol[col[e]];
        terms += abs(val[e] * sol[col[e]]);
    }
    return ku, terms;
}

fn region_reaction() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("region_reaction: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return sum; }
    for n in node_count() {
        if !on_region(n) { continue; }
        for p in ncomp {
            let d = ncomp * n + p;
            if fixed[d] != FIXED_BC { continue; }
            let ku, terms = row_force(d);
            sum[p] += ku - load[d] - assembled_source[d];
        }
    }
    if physics == HEAT { sum[0] -= region_ambient_heat(); }
    return sum;
}

fn region_mean(comp: i64) -> f64 {
    guard physics >= 0 && comp >= 0 && comp < ncomp && dofs_sized() else { abort("region_mean: a physics, 0 <= comp < ncomp and the arrays of the dofs"); }
    var kept: f64[>..] = [];
    kept.append(load);
    for &v in load { v = 0.0; }
    add_region_load([1.0; 3]);
    var sum = 0.0;
    var area = 0.0;
    for n in node_count() {
        let w = load[ncomp * n + comp];
        sum += w * sol[ncomp * n + comp];
        area += w;
    }
    for d in load.len { load[d] = kept[d]; }
    return if area > 0.0 { sum / area } else { 0.0 / 0.0 };
}

fn all_reaction() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("all_reaction: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return sum; }
    for d in sol.len {
        if fixed[d] != FIXED_BC { continue; }
        let ku, terms = row_force(d);
        sum[d % ncomp] += ku - load[d] - assembled_source[d];
    }
    return sum;
}

fn load_magnitude() -> f64 {
    guard dofs_sized() else { abort("load_magnitude: the arrays of the dofs"); }
    var sum = 0.0;
    for v in load { sum += abs(v); }
    for v in assembled_source { sum += abs(v); }
    guard row_start.len > 0 else { return sum; }
    var terms = 0.0;
    for d in sol.len {
        if fixed[d] != FIXED_BC { continue; }
        let ku, t = row_force(d);
        terms += t;
        sum += abs(ku - load[d] - assembled_source[d]);
    }
    return max(sum, terms);
}

fn total_load() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("total_load: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    for v, d in load { sum[d % ncomp] += v; }
    for v, d in assembled_source { sum[d % ncomp] += v; }
    sum[0] -= film_heat();
    return sum;
}
