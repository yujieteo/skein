// The `modal` and `buckling` steps of the beams (spec §14.2): the lowest
// eigenpairs of K x = λ B x by subspace iteration (K. J. Bathe, Finite
// Element Procedures, 1996, §11.6).
//
//   modal      B = M, the consistent mass: λ = ω²
//   buckling   B = -K_G, the geometric stiffness of the axial forces of the
//              last static solve: λ is the load factor of its load
//
// K is the skyline LDLᵀ of solve.goose, with the supports of the model. Each
// iteration is X̄ = K⁻¹ B X, then the Ritz pairs of the q vectors of X̄:
// (X̄ᵀ B X̄) z = μ (X̄ᵀ K X̄) z, with μ = 1/λ. The largest μ give the lowest
// λ, and X = X̄ Z. B can be indefinite (a beam in tension), but K is
// positive definite, so the Ritz problem is a symmetric one.
//
// Only the beams have B. A shell cell can add its mass and its K_G to
// `second` with the same layout later, 576 values for each cell.
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
let EIGEN_TOLERANCE = 1.0e-6;   // the residual of K⁻¹ B x - μ x in the K norm, relative to μ
let EIGEN_ITERATIONS = 100;

var eigen_kind = MODAL;
var eigenvalues: f64[>..] = [];  // λ of each mode of the last eigen solve: ω² or the load factor
var mode_shapes: f64[>..] = [];  // mode k at dof_count() k, with xᵀ K x = 1
var eigen_iterations = 0;        // the iterations of the last eigen solve
var eigen_error: u8[>..] = [];   // why the last eigen solve failed, or ""
var second: f64[>..] = [];       // B of each beam in global dofs: 144 values for each beam

// The `count` lowest natural modes of the beams. Each beam section needs a
// density (beam_density). False, with the reason in eigen_error, if K is
// singular or the iteration does not converge.
fn modes(count: i64) -> bool {
    guard count >= 1 && beams.len > 0 else { abort("modes: count >= 1, and a beam"); }
    for b in beams {
        if !(beam_sections[b.section].rho > 0.0) {
            return eigen_fail(str("beam section ", b.section, " has no density (beam_density)"));
        }
    }
    return eigen_solve(MODAL, count);
}

// The `count` lowest buckling modes of the beams under the load of the last
// solve. Call it after solve(), with the same supports.
fn buckling(count: i64) -> bool {
    guard count >= 1 && beams.len > 0 && disp.len == dof_count() else { abort("buckling: count >= 1, a beam, and a solution for each dof"); }
    return eigen_solve(BUCKLING, count);
}

// Eigenvalue k (1 ..) of the last eigen solve, as Skein's eigenvalue(k):
// the frequency ω/(2π) of a modal solve, the load factor of a buckling solve.
fn eigenvalue(k: i64) -> f64 {
    guard k >= 1 && k <= eigenvalues.len else { abort("eigenvalue: 1 <= k <= the count of the last eigen solve"); }
    let l = eigenvalues[k - 1];
    return if eigen_kind == MODAL { sqrt(max(l, 0.0)) / (2.0 * PI) } else { l };
}

// Copies mode k (1 ..) into disp, for write_vtu, with its largest
// displacement (or rotation, if it has no displacement) 1. It replaces the
// static solution, so call buckling first.
fn show_mode(k: i64) {
    guard k >= 1 && k * dof_count() <= mode_shapes.len && disp.len == dof_count() else { abort("show_mode: 1 <= k <= the count of the last eigen solve"); }
    let n = dof_count();
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

// Assembles and factors K, makes B, and runs the subspace iteration.
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
    // q vectors as Bathe: min(2 p, p + 8), at most the free dofs.
    let q = min(min(2 * count, count + 8), free);
    // X, X̄ and B X: 3 q vectors, and B of each beam.
    solvers::memory_check("shell skyline ldlt", 0, n, 3 * q + (144 * beams.len + n - 1) / n, 0, ", or fewer modes");
    make_second(kind);
    let sizes = state_sizes();
    let ok = subspace(count, q, n);
    guard state_sizes() == sizes else { abort("eigen_solve: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - t0;
    return ok;
}

// B of each beam into second: M, or -K_G of the axial force of disp.
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

// y = B x on the free dofs, 0 on the fixed dofs.
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

// v = K⁻¹ v with the factors of factor() in sky.
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

// The start vectors: a fixed pseudo-random value on each free dof.
fn start_vectors(q: i64, n: i64) -> f64[>..] {
    var x: f64[>..] = [];
    for i in q * n {
        let t = sin(12.9898 * (i % n) + 78.233 * (i / n) + 1.0) * 43758.5453;
        x.push(if fixed[i % n] == 0 { t - floor(t) - 0.5 } else { 0.0 });
    }
    return x;
}

// h[q i + j] = a_i · b_j of the q vectors of a and b, made symmetric.
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

// The subspace iteration for the `count` lowest λ with q vectors. With x_i
// a Ritz vector (x_iᵀ K x_i = 1, x_iᵀ B x_i = μ_i) and x̄_i = K⁻¹ B x_i, the
// residual is |x̄_i - μ_i x_i|²_K = x̄_iᵀ B x_i - μ_i², which is diagonal i of
// the next X̄ᵀ K X̄ = X̄ᵀ B X.
fn subspace(count: i64, q: i64, n: i64) -> bool {
    guard count >= 1 && q >= count && n == dof_count() else { abort("subspace: 1 <= count <= q, and n is the dof count"); }
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
    guard done else { return eigen_fail(str("no convergence in ", EIGEN_ITERATIONS, " iterations: residual ", worst)); };
    guard mu[count - 1] > 0.0 else { return eigen_fail(str("fewer than ", count, " positive eigenvalues: the load puts no beam in compression?")); };
    for i in count { eigenvalues.push(1.0 / mu[i]); }
    mode_shapes.append(x[..count * n]);
    return true;
}

// The Ritz pairs of br z = μ kr z: μ in decreasing order into mu, and the
// vectors z (column j, zᵀ kr z = I). With kr = L Lᵀ (Cholesky), it is
// C v = μ v for C = L⁻¹ br L⁻ᵀ, and z = L⁻ᵀ v. False if kr is not
// positive definite.
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
    // c = L⁻¹ c L⁻ᵀ: forward substitution on the columns, then on the rows.
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
    // z = L⁻ᵀ v: back substitution on each column.
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

// The eigenvalues of the symmetric q × q matrix c into mu, in decreasing
// order, and the eigenvectors (column j): cyclic Jacobi rotations.
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
    // Sorts the pairs: the largest μ first.
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

// c = Jᵀ c J and v = v J, with J the rotation of (p, r).
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
