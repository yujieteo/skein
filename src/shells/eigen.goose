import std;
import math;
import os;
namespace shells;
import .model;
import .beam;
import .solve;
import solvers from "../solvers";

let MODAL = 0;
let BUCKLING = 1;
let EIGEN_TOLERANCE = 1.0e-6;
let EIGEN_ITERATIONS = 100;

var eigen_kind = MODAL;
var eigenvalues: f64[>..] = [];
var mode_shapes: f64[>..] = [];
var eigen_iterations = 0;
var eigen_error: u8[>..] = [];
var too_few = false;
var second: f64[>..] = [];

fn modes(count: i64) -> bool {
    guard count >= 1 && beams.len > 0 else { abort("modes: count >= 1, and a beam"); }
    for b in beams {
        if !(beam_sections[b.section].rho > 0.0) {
            return eigen_fail(str("beam section ", b.section, " has no density (beam_density)"));
        }
    }
    return eigen_solve(MODAL, count);
}

fn buckling(count: i64) -> bool {
    guard count >= 1 && beams.len > 0 && disp.len == dof_count() && !disp_is_mode else { abort("buckling: count >= 1, a beam, and a static solution for each dof (show_mode replaced it)"); }
    return eigen_solve(BUCKLING, count);
}

fn eigenvalue(k: i64) -> f64 {
    guard k >= 1 && k <= eigenvalues.len else { abort("eigenvalue: 1 <= k <= the count of the last eigen solve"); }
    let l = eigenvalues[k - 1];
    return if eigen_kind == MODAL { sqrt(max(l, 0.0)) / (2.0 * PI) } else { l };
}

fn show_mode(k: i64) {
    guard k >= 1 && k * dof_count() <= mode_shapes.len && disp.len == dof_count() else { abort("show_mode: 1 <= k <= the count of the last eigen solve"); }
    let n = dof_count();
    disp_is_mode = true;
    var big: f64[2] = [0.0, 0.0];
    for d in n {
        let v = mode_shapes[n * (k - 1) + d];
        disp[d] = v;
        if abs(v) > abs(big[d % 6 / 3]) { big[d % 6 / 3] = v; }
    }
    let f = if big[0] != 0.0 { big[0] } else { big[1] };
    if f != 0.0 { for &v in disp { v /= f; } }
}

fn eigen_fail(why: u8[:]) -> bool {
    eigen_error.clear();
    eigen_error.append(why);
    return false;
}

fn too_few_fail(count: i64) -> bool {
    return eigen_fail(str("fewer than ", count, " positive eigenvalues: the load puts no beam in compression?"));
}

fn eigen_solve(kind: i64, count: i64) -> bool {
    let t0 = clock();
    eigen_kind = kind;
    eigen_error.clear();
    eigenvalues.clear();
    mode_shapes.clear();
    eigen_iterations = 0;
    assemble();
    guard factor() else { return eigen_fail(str("K is singular at dof ", singular_dof, ": a mechanism")); };
    let n = dof_count();
    var free = 0;
    for f in fixed { if f == 0 { free++; } }
    guard count <= free else { return eigen_fail(str(count, " modes, but only ", free, " free dofs")); };
    var q = min(min(2 * count, count + 8), free);
    make_second(kind);
    let sizes = state_sizes();
    var ok = false;
    let first_q = q;
    while true {
        solvers::memory_check("shell skyline ldlt", 0, n, 3 * q + (144 * beams.len + n - 1) / n, 0, ", or fewer modes");
        eigen_iterations = 0;
        ok = subspace(count, q, n);
        if !ok && !too_few && q > first_q { ok = too_few_fail(count); break; }
        if ok || eigen_kind != BUCKLING || !too_few || q == free { break; }
        q = min(2 * q, free);
    }
    guard state_sizes() == sizes else { abort("eigen_solve: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - t0;
    return ok;
}

fn make_second(kind: i64) {
    guard beams.len > 0 && disp.len == dof_count() else { abort("make_second: a beam, and a solution for each dof"); }
    second.clear();
    var m: f64[144] = [0.0; 144];
    for b in beams.len {
        if kind == MODAL {
            beam_mass(b, m);
        } else {
            let f = beam_forces(b);
            beam_geometric(b, 0.5 * (f[6] - f[0]), m);
            for &v in m { v = -v; }
        }
        for v in m { second.push(v); }
    }
}

fn apply_second(x: f64[:], y: f64[:]) {
    guard x.len == dof_count() && y.len == x.len && second.len == 144 * beams.len else { abort("apply_second: x and y have each dof, and B has each beam"); }
    for &v in y { v = 0.0; }
    var d: i64[12] = [0; 12];
    for b in beams.len {
        for l in 2 { for k in 6 { d[6 * l + k] = 6 * beams[b].n[l] + k; } }
        for a in 12 {
            if fixed[d[a]] != 0 { continue; }
            var s = 0.0;
            for c in 12 { if fixed[d[c]] == 0 { s += second[144 * b + 12 * a + c] * x[d[c]]; } }
            y[d[a]] += s;
        }
    }
}

fn back_solve(v: f64[:]) {
    guard v.len == dof_count() && top.len == v.len else { abort("back_solve: v and the profile have each dof"); }
    let n = v.len;
    for j in n {
        var s = v[j];
        for i in top[j]..j { s -= sky[at(i, j)] * v[i]; }
        v[j] = s;
    }
    for j in n { v[j] /= sky[at(j, j)]; }
    for jj in n {
        let j = n - 1 - jj;
        let x = v[j];
        for i in top[j]..j { v[i] -= sky[at(i, j)] * x; }
    }
}

fn start_vectors(q: i64, n: i64) -> f64[>..] {
    var x: f64[>..] = [];
    for i in q * n {
        let t = sin(12.9898 * (i % n) + 78.233 * (i / n) + 1.0) * 43758.5453;
        x.push(if fixed[i % n] == 0 { t - floor(t) - 0.5 } else { 0.0 });
    }
    return x;
}

fn project(a: f64[:], b: f64[:], h: f64[:], q: i64, n: i64) {
    guard a.len == q * n && b.len == a.len && h.len == q * q else { abort("project: q vectors of n in a and b, and q × q in h"); }
    for i in q {
        for j in q {
            var s = 0.0;
            for d in n { s += a[n * i + d] * b[n * j + d]; }
            h[q * i + j] = s;
        }
    }
    for i in q {
        for j in i + 1..q {
            let s = 0.5 * (h[q * i + j] + h[q * j + i]);
            h[q * i + j] = s;
            h[q * j + i] = s;
        }
    }
}

fn subspace(count: i64, q: i64, n: i64) -> bool {
    guard count >= 1 && q >= count && n == dof_count() else { abort("subspace: 1 <= count <= q, and n is the dof count"); }
    too_few = false;
    var x = start_vectors(q, n);
    var xb: f64[>..] = [];
    var y: f64[>..] = [];
    xb.append(x);
    y.append(x);
    var kr: f64[>..] = [];
    var br: f64[>..] = [];
    for i in q * q {
        kr.push(0.0);
        br.push(0.0);
    }
    var mu: f64[>..] = [];
    var worst = 0.0;
    var done = false;
    while !done && eigen_iterations < EIGEN_ITERATIONS {
        for j in q {
            apply_second(x[n * j..n * (j + 1)], y[n * j..n * (j + 1)]);
            for d in n { xb[n * j + d] = y[n * j + d]; }
            back_solve(xb[n * j..n * (j + 1)]);
        }
        project(xb, y, kr, q, n);
        if mu.len == q {
            worst = 0.0;
            for i in count { worst = max(worst, sqrt(max(kr[q * i + i] - mu[i] * mu[i], 0.0)) / abs(mu[i])); }
            done = worst <= EIGEN_TOLERANCE;
            if done { break; }
        }
        for j in q { apply_second(xb[n * j..n * (j + 1)], y[n * j..n * (j + 1)]); }
        project(xb, y, br, q, n);
        let z, good = ritz(kr, br, q, mu);
        guard good else { return eigen_fail("the vectors of the subspace are not independent (B x = 0 on each?)"); };
        for j in q {
            for d in n {
                var s = 0.0;
                for i in q { s += xb[n * i + d] * z[q * i + j]; }
                x[n * j + d] = s;
            }
        }
        eigen_iterations++;
    }
    too_few = mu.len == q && !(mu[count - 1] > 0.0);
    guard done else { return eigen_fail(str("no convergence in ", EIGEN_ITERATIONS, " iterations: residual ", worst)); };
    guard !too_few else { return too_few_fail(count); };
    for i in count { eigenvalues.push(1.0 / mu[i]); }
    mode_shapes.append(x[..count * n]);
    return true;
}

fn ritz(kr: f64[:], br: f64[:], q: i64, mu: f64[>..]&) -> f64[>..], bool {
    guard kr.len == q * q && br.len == kr.len else { abort("ritz: kr and br are q × q"); }
    var l: f64[>..] = [];
    l.append(kr);
    for j in q {
        for i in j..q {
            var s = l[q * i + j];
            for t in j { s -= l[q * i + t] * l[q * j + t]; }
            if i == j && !(s > 0.0) { return [], false; }
            l[q * i + j] = if i == j { sqrt(s) } else { s / l[q * j + j] };
        }
    }
    var c: f64[>..] = [];
    c.append(br);
    for col in q {
        for i in q {
            var s = c[q * i + col];
            for t in i { s -= l[q * i + t] * c[q * t + col]; }
            c[q * i + col] = s / l[q * i + i];
        }
    }
    for r in q {
        for i in q {
            var s = c[q * r + i];
            for t in i { s -= l[q * i + t] * c[q * r + t]; }
            c[q * r + i] = s / l[q * i + i];
        }
    }
    var v = jacobi_eigen(c, q, mu);
    for col in q {
        for ii in q {
            let i = q - 1 - ii;
            var s = v[q * i + col];
            for t in i + 1..q { s -= l[q * t + i] * v[q * t + col]; }
            v[q * i + col] = s / l[q * i + i];
        }
    }
    return v, true;
}

fn jacobi_eigen(c: f64[>..]&, q: i64, mu: f64[>..]&) -> f64[>..] {
    var v: f64[>..] = [];
    for i in q * q { v.push(if i / q == i % q { 1.0 } else { 0.0 }); }
    for sweep in 100 {
        var off = 0.0;
        var all = 0.0;
        for i in q * q {
            all += c[i] * c[i];
            if i / q != i % q { off += c[i] * c[i]; }
        }
        if off <= 1.0e-30 * all { break; }
        for p in q {
            for r in p + 1..q {
                if abs(c[q * p + r]) <= 1.0e-300 { continue; }
                let theta = (c[q * r + r] - c[q * p + p]) / (2.0 * c[q * p + r]);
                let t = (if theta >= 0.0 { 1.0 } else { -1.0 }) / (abs(theta) + sqrt(theta * theta + 1.0));
                let cs = 1.0 / sqrt(t * t + 1.0);
                rotate_pair(c, v, q, p, r, cs, t * cs);
            }
        }
    }
    var order: i64[>..] = [];
    for i in q { order.push(i); }
    for i in q {
        for j in i + 1..q {
            if c[q * order[j] + order[j]] > c[q * order[i] + order[i]] {
                let t = order[i];
                order[i] = order[j];
                order[j] = t;
            }
        }
    }
    mu.clear();
    for i in q { mu.push(c[q * order[i] + order[i]]); }
    var z: f64[>..] = [];
    for r in q { for j in q { z.push(v[q * r + order[j]]); } }
    return z;
}

fn rotate_pair(c: f64[>..]&, v: f64[>..]&, q: i64, p: i64, r: i64, cs: f64, sn: f64) {
    for k in q {
        let a = c[q * k + p];
        let b = c[q * k + r];
        c[q * k + p] = cs * a - sn * b;
        c[q * k + r] = sn * a + cs * b;
    }
    for k in q {
        let a = c[q * p + k];
        let b = c[q * r + k];
        c[q * p + k] = cs * a - sn * b;
        c[q * r + k] = sn * a + cs * b;
    }
    for k in q {
        let a = v[q * k + p];
        let b = v[q * k + r];
        v[q * k + p] = cs * a - sn * b;
        v[q * k + r] = sn * a + cs * b;
    }
}
