import std;
import math;
import dictionary;
import os;
namespace flow;
import solvers from "../solvers";
import .mesh;
import .cell;
import .solve;

let FSI_MAX = 50;
let FSI_RTOL = 1.0e-6;
let OMEGA_START = 0.5;
let MOTION_RTOL = 1.0e-10;

var coupled = false;
var young = 1.0;
var poisson = 0.3;
var fsi_iterations = 0;

var s_vertices = 0;
var s_coords: f64[>..] = [];
var s_conn: i64[>..] = [];
var s_boundary: u8[>..] = [];
var s_region: u8[>..] = [];
var s_fixed: u8[>..] = [];
var s_u: f64[>..] = [];
var s_position: i64[>..] = [];
var s_width = 0;
var s_band: f64[>..] = [];
var link: i64[>..] = [];

var rest: f64[>..] = [];
var m_rs: i64[>..] = [];
var m_col: i32[>..] = [];
var m_full: f64[>..] = [];
var m_val: f64[>..] = [];
var m_lu: f64[>..] = [];
var m_diag: i64[>..] = [];
var m_fixed: u8[>..] = [];

var finite = false;
var s_state: f64[>..] = [];

fn solid_finite() { finite = true; }

fn solid(e: f64, nu: f64) {
    guard e > 0.0 && nu > -1.0 && nu < 0.5 else { abort("solid: E > 0 and -1 < nu < 0.5"); }
    young = e;
    poisson = nu;
}

fn solid_nodes() -> i64 { return s_coords.len / 2; }
fn solid_cells() -> i64 { return s_conn.len / 6; }

fn read_solid(dir: const u8[:]) -> bool {
    var a: u8[>..] = [];
    var b: u8[>..] = [];
    guard read_file(str(dir, "/solid.xy"), a) && read_file(str(dir, "/solid.tri"), b) else { return false; }
    let xy, xy_ok = from_bytes<f64[>..]>(a);
    let tri, tri_ok = from_bytes<i64[>..]>(b);
    guard xy_ok && tri_ok && tri.len > 0 else { return false; }
    solid_from_arrays(xy, tri);
    return true;
}

fn solid_from_arrays(xy: f64[:], tri: i64[:]) {
    guard xy.len % 2 == 0 && tri.len % 3 == 0 && tri.len > 0 else { abort("solid_from_arrays: vertices and triangles"); }
    s_vertices = xy.len / 2;
    s_coords.clear();
    s_coords.append(xy);
    s_conn.clear();
    var mid = dictionary<i64, i64> {};
    var count = dictionary<i64, i64> {};
    for c in tri.len / 3 {
        var v: i64[3] = [tri[3 * c], tri[3 * c + 1], tri[3 * c + 2]];
        guard v[0] >= 0 && v[1] >= 0 && v[2] >= 0 && v[0] < s_vertices && v[1] < s_vertices && v[2] < s_vertices else {
            abort("solid_from_arrays: each vertex of a triangle < xy.len / 2");
        }
        let area = (s_coords[2 * v[1]] - s_coords[2 * v[0]]) * (s_coords[2 * v[2] + 1] - s_coords[2 * v[0] + 1]) -
                   (s_coords[2 * v[1] + 1] - s_coords[2 * v[0] + 1]) * (s_coords[2 * v[2]] - s_coords[2 * v[0]]);
        if area < 0.0 {
            let t = v[1];
            v[1] = v[2];
            v[2] = t;
        }
        s_conn.append(v[..]);
        for e in 3 {
            let a = v[EDGES[2 * e]];
            let b = v[EDGES[2 * e + 1]];
            let key = min(a, b) * s_vertices + max(a, b);
            count.update(key, 0) { it += 1; };
            let m = mid.get_or_insert(key, solid_nodes());
            if m == solid_nodes() {
                s_coords.push(0.5 * (s_coords[2 * a] + s_coords[2 * b]));
                s_coords.push(0.5 * (s_coords[2 * a + 1] + s_coords[2 * b + 1]));
            }
            s_conn.push(m);
        }
    }
    s_boundary.clear();
    s_region.clear();
    for n in solid_nodes() {
        s_boundary.push(0);
        s_region.push(0);
    }
    count.each() { key, k =>
        if k == 1 {
            s_boundary[key / s_vertices] = 1;
            s_boundary[key % s_vertices] = 1;
            s_boundary[mid.get_or_insert(key, 0)] = 1;
        }
    };
    solid_pattern();
    coupled = true;
}

fn node_csr(pairs: i64[:], n: i64, start: i64[>..]&, next: i64[>..]&) {
    start.clear();
    next.clear();
    start.push(0);
    var last = -1;
    for e in pairs {
        if e == last { continue; }
        last = e;
        while start.len <= e / n { start.push(next.len); }
        next.push(e % n);
    }
    while start.len <= n { start.push(next.len); }
}

fn solid_pattern() {
    let n = solid_nodes();
    var pairs: i64[>..] = [];
    for c in solid_cells() {
        for a in 6 {
            for b in 6 { pairs.push(s_conn[6 * c + a] * n + s_conn[6 * c + b]); }
        }
    }
    sort(pairs);
    var start: i64[>..] = [];
    var next: i64[>..] = [];
    node_csr(pairs, n, start, next);
    var seq: i64[>..] = [];
    breadth_first(0, start, next, seq);
    breadth_first(seq[n - 1], start, next, seq);
    s_position.clear();
    for d in 2 * n { s_position.push(0); }
    for i in n {
        s_position[2 * seq[i]] = 2 * i;
        s_position[2 * seq[i] + 1] = 2 * i + 1;
    }
    var width = 0;
    for e in pairs { width = max(width, abs(s_position[2 * (e / n)] - s_position[2 * (e % n)]) + 2); }
    s_width = width;
    s_fixed.clear();
    s_u.clear();
    for d in 2 * n {
        s_fixed.push(0);
        s_u.push(0.0);
    }
    solvers::zeros(s_band, 2 * n * (s_width + 1));
}

fn breadth_first(from: i64, start: i64[:], next: i64[:], seq: i64[>..]&) {
    let n = start.len - 1;
    var seen: u8[>..] = [];
    for i in n { seen.push(0); }
    seq.clear();
    seq.push(from);
    seen[from] = 1;
    var head = 0;
    var k = 0;
    while seq.len < n {
        if head == seq.len {
            while seen[k] != 0 { k++; }
            seq.push(k);
            seen[k] = 1;
        }
        let v = seq[head];
        head++;
        for e in start[v]..start[v + 1] {
            if seen[next[e]] == 0 {
                seen[next[e]] = 1;
                seq.push(next[e]);
            }
        }
    }
}

fn select_solid_boundary<P>() -> bool {
    guard solid_cells() > 0 else { abort("select_solid_boundary: a solid"); }
    var found = false;
    for n in solid_nodes() {
        s_region[n] = if s_boundary[n] != 0 && P(s_coords[2 * n], s_coords[2 * n + 1]) { 1 } else { 0 };
    }
    for c in solid_cells() {
        for e in 3 {
            let m = s_conn[6 * c + 3 + e];
            if s_boundary[m] != 0 && s_region[s_conn[6 * c + EDGES[2 * e]]] != 0 && s_region[s_conn[6 * c + EDGES[2 * e + 1]]] != 0 { s_region[m] = 1; }
        }
    }
    for n in solid_nodes() { if s_region[n] != 0 { found = true; } }
    return found;
}

fn clear_solid_bcs() {
    for d in s_fixed.len { s_fixed[d] = 0; }
}

fn fix_solid() {
    guard s_region.len == solid_nodes() else { abort("fix_solid: a solid"); }
    for n in solid_nodes() {
        if s_region[n] != 0 {
            s_fixed[2 * n] = 1;
            s_fixed[2 * n + 1] = 1;
        }
    }
}

fn solid_gradients(c: i64) -> f64[6], f64 {
    let a = s_conn[6 * c];
    let b = s_conn[6 * c + 1];
    let d = s_conn[6 * c + 2];
    let x: f64[3] = [s_coords[2 * a], s_coords[2 * b], s_coords[2 * d]];
    let y: f64[3] = [s_coords[2 * a + 1], s_coords[2 * b + 1], s_coords[2 * d + 1]];
    let twice = (x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]);
    var g: f64[6] = [0.0; 6];
    for k in 3 {
        let i = (k + 1) % 3;
        let j = (k + 2) % 3;
        g[2 * k] = (y[i] - y[j]) / twice;
        g[2 * k + 1] = (x[j] - x[i]) / twice;
    }
    return g, 0.5 * twice;
}

fn p2_gradients(g: f64[6], l: f64[3]) -> f64[12] {
    var d: f64[12] = [0.0; 12];
    for k in 3 {
        for t in 2 { d[2 * k + t] = (4.0 * l[k] - 1.0) * g[2 * k + t]; }
    }
    for e in 3 {
        let a = EDGES[2 * e];
        let b = EDGES[2 * e + 1];
        for t in 2 { d[2 * (3 + e) + t] = 4.0 * (l[a] * g[2 * b + t] + l[b] * g[2 * a + t]); }
    }
    return d;
}

fn solid_cell(c: i64, k: f64[144]&) {
    if finite {
        svk_tangent(c, k);
        return;
    }
    let g, area = solid_gradients(c);
    let lam = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let shear = young / (2.0 * (1.0 + poisson));
    for i in 144 { k[i] = 0.0; }
    for q in 3 {
        var l: f64[3] = [1.0 / 6.0; 3];
        l[q] = 2.0 / 3.0;
        let d = p2_gradients(g, l);
        let w = area / 3.0;
        for a in 6 {
            let ax = d[2 * a];
            let ay = d[2 * a + 1];
            for b in 6 {
                let bx = d[2 * b];
                let by = d[2 * b + 1];
                k[24 * a + 2 * b] += w * ((lam + 2.0 * shear) * ax * bx + shear * ay * by);
                k[24 * a + 2 * b + 1] += w * (lam * ax * by + shear * ay * bx);
                k[24 * a + 12 + 2 * b] += w * (lam * ay * bx + shear * ax * by);
                k[24 * a + 12 + 2 * b + 1] += w * ((lam + 2.0 * shear) * ay * by + shear * ax * bx);
            }
        }
    }
}

fn svk_force(c: i64, ue: f64[12], f: f64[12]&) {
    let g, area = solid_gradients(c);
    let lam = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let shear = young / (2.0 * (1.0 + poisson));
    for i in 12 { f[i] = 0.0; }
    for q in 3 {
        var l: f64[3] = [1.0 / 6.0; 3];
        l[q] = 2.0 / 3.0;
        let d = p2_gradients(g, l);
        var h: f64[4] = [0.0; 4];
        for a in 6 {
            for i in 2 { for j in 2 { h[2 * i + j] += ue[2 * a + i] * d[2 * a + j]; } }
        }
        let fg: f64[4] = [1.0 + h[0], h[1], h[2], 1.0 + h[3]];
        var e: f64[4] = [0.0; 4];
        for i in 2 { for j in 2 {
            var t = 0.0;
            for m in 2 { t += fg[2 * m + i] * fg[2 * m + j]; }
            let eye = if i == j { 1.0 } else { 0.0 };
            e[2 * i + j] = 0.5 * (t - eye);
        } }
        let tr = e[0] + e[3];
        var sg: f64[4] = [0.0; 4];
        for i in 4 { sg[i] = 2.0 * shear * e[i]; }
        sg[0] += lam * tr;
        sg[3] += lam * tr;
        for a in 6 {
            for i in 2 { for j in 2 { for m in 2 {
                f[2 * a + i] += area / 3.0 * fg[2 * i + m] * sg[2 * m + j] * d[2 * a + j];
            } } }
        }
    }
}

fn svk_tangent(c: i64, k: f64[144]&) {
    var ue: f64[12] = [0.0; 12];
    for a in 12 { ue[a] = s_state[2 * s_conn[6 * c + a / 2] + a % 2]; }
    let step = 1.0e-7;
    var fp: f64[12] = [0.0; 12];
    var fm: f64[12] = [0.0; 12];
    for j in 12 {
        let keep = ue[j];
        ue[j] = keep + step;
        svk_force(c, ue, fp);
        ue[j] = keep - step;
        svk_force(c, ue, fm);
        ue[j] = keep;
        for i in 12 { k[12 * i + j] = (fp[i] - fm[i]) / (2.0 * step); }
    }
}

fn svk_internal(u: f64[:], g: f64[:]) {
    for d in g.len { g[d] = 0.0; }
    var ue: f64[12] = [0.0; 12];
    var f: f64[12] = [0.0; 12];
    for c in solid_cells() {
        for a in 12 { ue[a] = u[2 * s_conn[6 * c + a / 2] + a % 2]; }
        svk_force(c, ue, f);
        for a in 12 { g[2 * s_conn[6 * c + a / 2] + a % 2] += f[a]; }
    }
}

fn band(i: i64, j: i64) -> i64 { return (s_width + 1) * i + i - j; }

fn factor_solid() -> bool {
    for e in s_band.len { s_band[e] = 0.0; }
    var k: f64[144] = [0.0; 144];
    for c in solid_cells() {
        solid_cell(c, k);
        for a in 12 {
            let row = 2 * s_conn[6 * c + a / 2] + a % 2;
            for b in 12 {
                let column = 2 * s_conn[6 * c + b / 2] + b % 2;
                let i = s_position[row];
                let j = s_position[column];
                if j > i { continue; }
                let held = s_fixed[row] != 0 || s_fixed[column] != 0;
                s_band[band(i, j)] = if held { if i == j { 1.0 } else { 0.0 } } else { s_band[band(i, j)] + k[12 * a + b] };
            }
        }
    }
    let n = s_fixed.len;
    for i in n {
        for j in max(0, i - s_width + 1)..i + 1 {
            var t = s_band[band(i, j)];
            for m in max(0, i - s_width + 1)..j { t -= s_band[band(i, m)] * s_band[band(j, m)]; }
            if i == j {
                if t <= 0.0 { return false; }
                s_band[band(i, i)] = sqrt(t);
            } else {
                s_band[band(i, j)] = t / s_band[band(j, j)];
            }
        }
    }
    return true;
}

fn solve_solid(load: f64[:], u: f64[:]) {
    if !finite {
        solve_linear_solid(load, u);
        return;
    }
    let n = load.len;
    var g: f64[>..] = [];
    var r: f64[>..] = [];
    var du: f64[>..] = [];
    solvers::zeros(g, n);
    solvers::zeros(r, n);
    solvers::zeros(du, n);
    s_state.clear();
    s_state.append(u[..]);
    for it in 20 {
        for d in n { s_state[d] = u[d]; }
        guard factor_solid() else { abort("solve_solid: the tangent is not positive definite"); }
        svk_internal(u, g);
        for d in n { r[d] = load[d] - g[d]; }
        solve_linear_solid(r, du);
        for d in n { u[d] += du[d]; }
        if solvers::norm(du) <= 1.0e-10 * solvers::norm(u) + 1.0e-16 { return; }
    }
    abort("solve_solid: the Newton iteration did not converge in 20 steps");
}

fn solve_linear_solid(load: f64[:], u: f64[:]) {
    guard load.len == s_fixed.len && u.len == s_fixed.len else { abort("solve_solid: a load and a u for each dof"); }
    let n = load.len;
    var y: f64[>..] = [];
    solvers::zeros(y, n);
    for d in n { y[s_position[d]] = if s_fixed[d] != 0 { 0.0 } else { load[d] }; }
    for i in n {
        for m in max(0, i - s_width + 1)..i { y[i] -= s_band[band(i, m)] * y[m]; }
        y[i] /= s_band[band(i, i)];
    }
    for r in n {
        let i = n - 1 - r;
        y[i] /= s_band[band(i, i)];
        for m in max(0, i - s_width + 1)..i { y[m] -= s_band[band(i, m)] * y[i]; }
    }
    for d in n { u[d] = y[s_position[d]]; }
}

fn link_interface() -> bool {
    let tol = tolerance();
    link.clear();
    var found = false;
    for s in solid_nodes() {
        var f = -1;
        if s_boundary[s] != 0 {
            let x = s_coords[2 * s];
            let y = s_coords[2 * s + 1];
            for n in node_count() {
                if on_boundary[n] != 0 && abs(coords[2 * n] - x) <= tol && abs(coords[2 * n + 1] - y) <= tol { f = n; }
            }
        }
        link.push(f);
        if f < 0 { continue; }
        found = true;
        for t in 2 {
            fixed[vdof(f, t)] = 1;
            prescribed[vdof(f, t)] = 0.0;
        }
    }
    return found;
}

fn motion_setup() -> bool {
    let n = vertex_count;
    var pairs: i64[>..] = [];
    for c in cell_count() {
        for a in 3 {
            for b in 3 { pairs.push(cell_node(c, a) * n + cell_node(c, b)); }
        }
    }
    sort(pairs);
    var start: i64[>..] = [];
    var next: i64[>..] = [];
    node_csr(pairs, n, start, next);
    m_rs.clear();
    m_rs.append(start[..]);
    m_col.clear();
    for e in next { m_col.push(e as! i32); }
    solvers::zeros(m_full, m_col.len);
    for c in cell_count() {
        let g = geometry(c);
        for a in 3 {
            for b in 3 {
                let e = solvers::entry(m_rs, m_col, cell_node(c, a), cell_node(c, b));
                m_full[e] += g.area * (g.gl[2 * a] * g.gl[2 * b] + g.gl[2 * a + 1] * g.gl[2 * b + 1]);
            }
        }
    }
    m_fixed.clear();
    m_diag.clear();
    for v in n {
        m_fixed.push(on_boundary[v]);
        m_diag.push(0);
    }
    m_val.clear();
    m_val.append(m_full[..]);
    solvers::zeros(m_lu, m_col.len);
    var b: f64[>..] = [];
    solvers::zeros(b, n);
    solvers::constrain(m_rs, m_col, m_val, m_fixed, b);
    return solvers::ilu0(m_rs, m_col, m_val, m_lu, m_diag);
}

fn move_mesh() -> bool {
    let n = vertex_count;
    var g: f64[>..] = [];
    var b: f64[>..] = [];
    var x: f64[>..] = [];
    for t in 2 {
        solvers::zeros(g, n);
        solvers::zeros(b, n);
        solvers::zeros(x, n);
        for s in solid_nodes() {
            if link[s] >= 0 && link[s] < n { g[link[s]] = s_u[2 * s + t]; }
        }
        solvers::multiply(m_rs, m_col, m_full, g, b);
        for v in n { b[v] = if m_fixed[v] != 0 { 0.0 } else { -b[v] }; }
        let its, solved = solvers::gmres(m_rs, m_col, m_val, b, x, MOTION_RTOL, LINEAR_RESTART, LINEAR_MAX) { r, z =>
            solvers::ilu0_apply(m_rs, m_col, m_lu, m_diag, r, z);
        };
        if !solved { return false; }
        for v in n { coords[2 * v + t] = rest[2 * v + t] + x[v] + g[v]; }
    }
    for c in cell_count() {
        for e in per_cell - 3 {
            let m = cell_node(c, 3 + e);
            let a = cell_node(c, EDGES[2 * e]);
            let d = cell_node(c, EDGES[2 * e + 1]);
            for t in 2 { coords[2 * m + t] = 0.5 * (coords[2 * a + t] + coords[2 * d + t]); }
        }
    }
    return true;
}

fn solve_fsi(rtol: f64) -> bool {
    guard coupled && !turbulent && fixed.len == dof_count() && rtol > 0.0 else { abort("solve_fsi: a solid, a laminar flow, clear_bcs, and rtol > 0"); }
    let start = clock();
    guard link_interface() else { abort("solve_fsi: no node of the solid is on the fluid mesh"); }
    solvers::zeros(s_state, s_fixed.len);
    if !factor_solid() || !motion_setup() { return false; }
    rest.clear();
    rest.append(coords[..]);
    let ns = s_u.len;
    var load: f64[>..] = [];
    var trial: f64[>..] = [];
    var r: f64[>..] = [];
    var last: f64[>..] = [];
    solvers::zeros(load, ns);
    solvers::zeros(trial, ns);
    solvers::zeros(r, ns);
    solvers::zeros(last, ns);
    for d in ns { s_u[d] = 0.0; }
    var omega = OMEGA_START;
    var steps = 0;
    var its = 0;
    var converged = false;
    warm = false;
    fsi_iterations = 0;
    while fsi_iterations < FSI_MAX && move_mesh() {
        let flow_ok = solve_steady(rtol);
        steps += newton_steps;
        its += linear_iterations;
        if !flow_ok { break; }
        warm = true;
        for s in solid_nodes() {
            for t in 2 { load[2 * s + t] = if link[s] >= 0 { -raw[vdof(link[s], t)] } else { 0.0 }; }
        }
        for d in ns { trial[d] = s_u[d]; }
        solve_solid(load, trial);
        for d in ns { r[d] = trial[d] - s_u[d]; }
        fsi_iterations++;
        if solvers::norm(r) <= FSI_RTOL * solvers::norm(trial) {
            for d in ns { s_u[d] = trial[d]; }
            converged = true;
            break;
        }
        if fsi_iterations > 1 { omega = aitken(omega, r, last); }
        for d in ns {
            s_u[d] += omega * r[d];
            last[d] = r[d];
        }
    }
    warm = false;
    for i in coords.len { coords[i] = rest[i]; }
    newton_steps = steps;
    linear_iterations = its;
    last_seconds = clock() - start;
    return converged;
}

fn aitken(omega: f64, r: f64[:], last: f64[:]) -> f64 {
    var top = 0.0;
    var bottom = 0.0;
    for d in r.len {
        top += last[d] * (r[d] - last[d]);
        bottom += (r[d] - last[d]) * (r[d] - last[d]);
    }
    return if bottom > 0.0 { -omega * top / bottom } else { omega };
}

fn solid_value_at(comp: i64, x: f64, y: f64) -> f64, bool {
    guard (comp == 0 || comp == 1) && solid_cells() > 0 else { abort("solid_value_at: comp 0 or 1, and a solid"); }
    let tol = 1.0e-9;
    for c in solid_cells() {
        let g, area = solid_gradients(c);
        var l: f64[3] = [0.0; 3];
        for k in 3 {
            let i = s_conn[6 * c + (k + 1) % 3];
            l[k] = g[2 * k] * (x - s_coords[2 * i]) + g[2 * k + 1] * (y - s_coords[2 * i + 1]);
        }
        if l[0] < -tol || l[1] < -tol || l[2] < -tol { continue; }
        var v = 0.0;
        for k in 3 { v += l[k] * (2.0 * l[k] - 1.0) * s_u[2 * s_conn[6 * c + k] + comp]; }
        for e in 3 { v += 4.0 * l[EDGES[2 * e]] * l[EDGES[2 * e + 1]] * s_u[2 * s_conn[6 * c + 3 + e] + comp]; }
        return v, true;
    }
    return 0.0, false;
}

fn write_solid_vtu(path: const u8[:], x_ref: f64) -> bool {
    let nodes = solid_nodes();
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", solid_cells(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, x_ref * s_coords[2 * n], " ", x_ref * s_coords[2 * n + 1], " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in solid_cells() {
        for l in 6 { format(out, s_conn[6 * c + l], if l < 5 { " " } else { "\n" }); }
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in solid_cells() { format(out, 6 * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in solid_cells() { out.append("22\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    out.append("<DataArray type=\"Float64\" Name=\"u\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, x_ref * s_u[2 * n], " ", x_ref * s_u[2 * n + 1], " 0\n"); }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
