import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .solve;
import .stages;
import solvers from "../solvers";

let MASS = 0;
let GEOMETRIC = 1;
let EIGEN_TOLERANCE = 1.0e-8;
let EIGEN_ITERATIONS = 500;
let FREE_BODY_SHIFT = 1.0e-6;
let LOBPCG_BLOCKS = 14;

var density = 0.0;
var val_b: f64[>..] = [];
var eigen_kind = MASS;
var eigenvalues: f64[>..] = [];
var modes: f64[>..] = [];

let GAUSS_X: f64[4] = [0.0694318442029737, 0.33000947820757187, 0.6699905217924281, 0.9305681557970263];
let GAUSS_W: f64[4] = [0.17392742256872692, 0.3260725774312731, 0.3260725774312731, 0.17392742256872692];

fn rich_count() -> i64 { return if hex_cells { quad_count(dim) } else if dim == 2 { 16 } else { 64 }; }

fn rich_point(q: i64) -> f64[3], f64 {
    if hex_cells { return quad_point(dim, q); }
    let u = GAUSS_X[q % 4];
    let v = GAUSS_X[q / 4 % 4];
    let w = GAUSS_X[q / 16];
    let weight = GAUSS_W[q % 4] * GAUSS_W[q / 4 % 4] * (if dim == 2 { 1.0 } else { GAUSS_W[q / 16] });
    if dim == 2 { return [u, v * (1.0 - u), 0.0], weight * (1.0 - u); }
    return [u, v * (1.0 - u), w * (1.0 - u) * (1.0 - v)], weight * (1.0 - u) * (1.0 - u) * (1.0 - v);
}

fn stress_here(c: i64) -> f64[9] {
    let b = born.len > 0 && born[c] != 0;
    var du: f64[9] = [0.0; 9];
    for l in per_cell {
        let n = cell_node(c, l);
        for a in ncomp {
            let v = held[ncomp * n + a] - (if b { birth_u[(per_cell * c + l) * ncomp + a] } else { 0.0 });
            for t in dim { du[3 * a + t] += v * grad[3 * l + t]; }
        }
    }
    var tr = 0.0;
    for a in dim { tr += du[4 * a]; }
    var s: f64[9] = [0.0; 9];
    for i in dim {
        for j in dim { s[3 * i + j] = shear * (du[3 * i + j] + du[3 * j + i]) + (if i == j { lame * tr } else { 0.0 }); }
    }
    return s;
}

fn assemble_second(kind: i64) {
    guard physics == ELASTIC && val.len > 0 else { abort("assemble_second: an elastic physics and K"); }
    val_b.clear();
    for e in val.len { val_b.push(0.0); }
    var me: f64[100] = [0.0; 100];
    for c in cell_count() {
        if !is_live(c) { continue; }
        use_law(c);
        for &v in me { v = 0.0; }
        for q in rich_count() {
            let xi, w = rich_point(q);
            let dv = w * cell_geometry(c, xi) * thickness;
            let s = if kind == GEOMETRIC { stress_here(c) } else { [0.0; 9] };
            for i in per_cell {
                for j in per_cell {
                    var g = 0.0;
                    if kind == MASS {
                        g = density * shp[i] * shp[j];
                    } else {
                        for k in dim {
                            for l in dim { g += grad[3 * i + k] * s[3 * k + l] * grad[3 * j + l]; }
                        }
                    }
                    me[per_cell * i + j] += g * dv;
                }
            }
        }
        for i in per_cell {
            for j in per_cell {
                for p in ncomp { val_b[entry(ncomp * cell_node(c, i) + p, ncomp * cell_node(c, j)) + p] += me[per_cell * i + j]; }
            }
        }
    }
    use_law_0();
}

fn apply(second: bool, x: f64[:], y: f64[:]) {
    let a = eigen_kind == MASS;
    let use_b = if second { a } else { !a };
    for r in y.len {
        var s = 0.0;
        if fixed[r] == 0 {
            if use_b {
                for e in row_start[r]..row_start[r + 1] { s += val_b[e] * x[col[e]]; }
            } else {
                for e in row_start[r]..row_start[r + 1] { s += val[e] * x[col[e]]; }
            }
        }
        y[r] = s;
    }
}

fn symmetric_eigen(h: f64[:], k: i64) -> f64[>..], f64[>..] {
    var a: f64[>..] = [];
    a.append(h);
    var v: f64[>..] = [];
    for i in k * k { v.push(if i / k == i % k { 1.0 } else { 0.0 }); }
    for sweep in 100 {
        var off = 0.0;
        var all = 0.0;
        for i in k * k {
            all += a[i] * a[i];
            if i / k != i % k { off += a[i] * a[i]; }
        }
        if off <= 1.0e-30 * all { break; }
        for p in k {
            for q in p + 1..k {
                if abs(a[k * p + q]) <= 1.0e-300 { continue; }
                let theta = (a[k * q + q] - a[k * p + p]) / (2.0 * a[k * p + q]);
                let t = (if theta >= 0.0 { 1.0 } else { -1.0 }) / (abs(theta) + sqrt(theta * theta + 1.0));
                let cs = 1.0 / sqrt(t * t + 1.0);
                let sn = t * cs;
                rotate(a, v, k, p, q, cs, sn);
            }
        }
    }
    return sorted_pairs(a, v, k);
}

fn rotate(a: f64[>..]&, v: f64[>..]&, k: i64, p: i64, q: i64, cs: f64, sn: f64) {
    for r in k {
        let x = a[k * r + p];
        let y = a[k * r + q];
        a[k * r + p] = cs * x - sn * y;
        a[k * r + q] = sn * x + cs * y;
    }
    for r in k {
        let x = a[k * p + r];
        let y = a[k * q + r];
        a[k * p + r] = cs * x - sn * y;
        a[k * q + r] = sn * x + cs * y;
    }
    for r in k {
        let x = v[k * r + p];
        let y = v[k * r + q];
        v[k * r + p] = cs * x - sn * y;
        v[k * r + q] = sn * x + cs * y;
    }
}

fn sorted_pairs(a: f64[:], v: f64[:], k: i64) -> f64[>..], f64[>..] {
    var order: i64[>..] = [];
    for i in k { order.push(i); }
    for i in k {
        for j in i + 1..k {
            if a[k * order[j] + order[j]] < a[k * order[i] + order[i]] {
                let t = order[i];
                order[i] = order[j];
                order[j] = t;
            }
        }
    }
    var values: f64[>..] = [];
    var vectors: f64[>..] = [];
    for i in k { values.push(a[k * order[i] + order[i]]); }
    for r in k {
        for j in k { vectors.push(v[k * r + order[j]]); }
    }
    return values, vectors;
}

var basis: f64[>..] = [];
var basis_b: f64[>..] = [];
var basis_a: f64[>..] = [];

fn dot_of(a: f64[:], i: i64, b: f64[:], j: i64, n: i64) -> f64 {
    var s = 0.0;
    for d in n { s += a[n * i + d] * b[n * j + d]; }
    return s;
}

fn add_to_basis(v: f64[:]) -> bool {
    let n = v.len;
    var u: f64[>..] = [];
    u.append(v);
    var bu: f64[>..] = [];
    for d in n { bu.push(0.0); }
    apply(true, u, bu);
    var first = 0.0;
    for d in n { first += u[d] * bu[d]; }
    guard first > 0.0 else { return false; }
    let k = basis.len / n;
    for pass in 2 {
        for i in k {
            var c = 0.0;
            for d in n { c += u[d] * basis_b[n * i + d]; }
            for d in n {
                u[d] -= c * basis[n * i + d];
                bu[d] -= c * basis_b[n * i + d];
            }
        }
    }
    var norm = 0.0;
    for d in n { norm += u[d] * bu[d]; }
    guard norm > 1.0e-16 * first else { return false; }
    let f = 1.0 / sqrt(norm);
    for &x in u { x *= f; }
    apply(true, u, bu);
    basis.append(u);
    basis_b.append(bu);
    apply(false, u, bu);
    basis_a.append(bu);
    return true;
}

fn ritz_vectors(src: f64[:], y: f64[:], k: i64, m: i64, from: i64, n: i64) -> f64[>..] {
    var out: f64[>..] = [];
    for j in m {
        for d in n { out.push(0.0); }
        for i in from..k {
            let f = y[k * i + j];
            for d in n { out[n * j + d] += f * src[n * i + d]; }
        }
    }
    return out;
}

fn start_block(m: i64, n: i64) -> f64[>..] {
    var x: f64[>..] = [];
    for i in m * n {
        let t = sin(12.9898 * (i % n) + 78.233 * (i / n) + 1.0) * 43758.5453;
        x.push(if fixed[i % n] == 0 { t - floor(t) - 0.5 } else { 0.0 });
    }
    return x;
}

var eigen_its = 0;

fn lobpcg(count: i64) -> bool {
    let n = dof_total();
    var free = 0;
    for d in n { if fixed[d] == 0 { free++; } }
    let m = min(count + max(count / 2, 2), free);
    guard count >= 1 && count <= m else { return false; }
    var x = start_block(m, n);
    var w: f64[>..] = [];
    var p: f64[>..] = [];
    var converged = false;
    eigen_its = 0;
    while !converged && eigen_its < EIGEN_ITERATIONS {
        basis.clear();
        basis_b.clear();
        basis_a.clear();
        for j in m { add_to_basis(x[n * j..n * (j + 1)]); }
        let nx = basis.len / n;
        for j in w.len / n { add_to_basis(w[n * j..n * (j + 1)]); }
        for j in p.len / n { add_to_basis(p[n * j..n * (j + 1)]); }
        let k = basis.len / n;
        guard nx == m else { return false; }
        let theta, y = ritz(k, n);
        x = ritz_vectors(basis, y, k, m, 0, n);
        let ax = ritz_vectors(basis_a, y, k, m, 0, n);
        let bx = ritz_vectors(basis_b, y, k, m, 0, n);
        p = if k > nx { ritz_vectors(basis, y, k, m, nx, n) } else { [] };
        converged, w = residuals(theta, ax, bx, count, m, n);
        eigen_its++;
        eigenvalues.clear();
        eigenvalues.append(theta[..count]);
    }
    modes.clear();
    modes.append(x[..count * n]);
    return converged;
}

fn ritz(k: i64, n: i64) -> f64[>..], f64[>..] {
    var h: f64[>..] = [];
    var l: f64[>..] = [];
    for i in k {
        for j in k {
            h.push(0.5 * (dot_of(basis, i, basis_a, j, n) + dot_of(basis, j, basis_a, i, n)));
            l.push(0.5 * (dot_of(basis, i, basis_b, j, n) + dot_of(basis, j, basis_b, i, n)));
        }
    }
    for j in k {
        for i in j..k {
            var s = l[k * i + j];
            for t in j { s -= l[k * i + t] * l[k * j + t]; }
            l[k * i + j] = if i == j { sqrt(max(s, 1.0e-300)) } else { s / l[k * j + j] };
        }
    }
    for c in k {
        for i in k {
            var s = h[k * i + c];
            for t in i { s -= l[k * i + t] * h[k * t + c]; }
            h[k * i + c] = s / l[k * i + i];
        }
    }
    for r in k {
        for i in k {
            var s = h[k * r + i];
            for t in i { s -= l[k * i + t] * h[k * r + t]; }
            h[k * r + i] = s / l[k * i + i];
        }
    }
    let theta, z = symmetric_eigen(h, k);
    var y: f64[>..] = [];
    y.append(z);
    for c in k {
        for ii in k {
            let i = k - 1 - ii;
            var s = y[k * i + c];
            for t in i + 1..k { s -= l[k * t + i] * y[k * t + c]; }
            y[k * i + c] = s / l[k * i + i];
        }
    }
    return theta, y;
}

fn residuals(theta: f64[:], ax: f64[:], bx: f64[:], count: i64, m: i64, n: i64) -> bool, f64[>..] {
    var w: f64[>..] = [];
    var r: f64[>..] = [];
    var z: f64[>..] = [];
    for d in n {
        r.push(0.0);
        z.push(0.0);
    }
    var all = true;
    for j in m {
        var r2 = 0.0;
        var a2 = 0.0;
        var b2 = 0.0;
        for d in n {
            r[d] = ax[n * j + d] - theta[j] * bx[n * j + d];
            r2 += r[d] * r[d];
            a2 += ax[n * j + d] * ax[n * j + d];
            b2 += bx[n * j + d] * bx[n * j + d];
        }
        let done = sqrt(r2) <= EIGEN_TOLERANCE * (sqrt(a2) + abs(theta[j]) * sqrt(b2));
        if j < count && !done { all = false; }
        if done && j < count { continue; }
        precondition(r, z);
        w.append(z);
    }
    return all, w;
}

fn free_body_shift(kind: i64) -> f64 {
    guard kind == MASS else { return 0.0; }
    for d in dof_total() { if fixed[d] != 0 { return 0.0; } }
    var dk = 0.0;
    var dm = 0.0;
    for r in dof_total() {
        dk += val[entry(r, r)];
        dm += val_b[entry(r, r)];
    }
    let s = FREE_BODY_SHIFT * dk / dm;
    for e in val.len { val[e] += s * val_b[e]; }
    return s;
}

fn eigen_solve(kind: i64, count: i64) -> bool {
    guard physics == ELASTIC && dofs_sized() && count >= 1 else { abort("eigen_solve: an elastic physics, the arrays of the dofs and count >= 1"); }
    guard kind == GEOMETRIC || density > 0.0 else { abort("eigen_solve: a density > 0 for a modal step"); }
    guard kind == MASS || held.len == sol.len else { abort("eigen_solve: the solution of a step before a buckling step"); }
    let start = clock();
    let m = count + max(count / 2, 2);
    assemble(if kind == MASS { "simplex modal lobpcg amg" } else { "simplex buckling lobpcg amg" }, LOBPCG_BLOCKS * m + 2 * mode_count());
    hold_dead_nodes();
    eigen_kind = kind;
    assemble_second(kind);
    let shift = free_body_shift(kind);
    let factored = setup_preconditioner(0);
    let done = factored && lobpcg(count);
    for t in eigenvalues.len { eigenvalues[t] -= shift; }
    last_iterations = eigen_its;
    last_seconds = clock() - start;
    return done;
}

fn mode_value(k: i64) -> f64 {
    guard k >= 1 && k <= eigenvalues.len else { abort("mode_value: 1 <= k <= the count of the last eigen step"); }
    let t = eigenvalues[k - 1];
    return if eigen_kind == MASS { sqrt(max(t, 0.0)) / TAU } else { -1.0 / t };
}
