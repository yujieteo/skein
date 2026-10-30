// The Krylov solvers: CG for symmetric positive definite systems and
// restarted GMRES(m) for the others. Each takes the preconditioner as a
// function value P(r, z) that sets z = M⁻¹ r, for example
//
//     solvers::gmres(row_start, col, val, b, x, 1.0e-10, 50, 1000) { r, z =>
//         solvers::ilu0_apply(row_start, col, lu, diag, r, z);
//     };
//
// Both start from the x that they get, and stop when the norm of the
// residual b - A x is at most rtol times the norm of b. Each returns the
// number of iterations (matrix-vector products in the loop) and true if
// it converged.
import std;
import math;
namespace solvers;
import .sparse;

// Preconditioned CG (Saad, Algorithm 9.1).
fn cg<P>(row_start: i64[:], col: i32[:], val: f64[:], b: f64[:], x: f64[:],
         rtol: f64, max_iterations: i64) -> i64, bool {
    let n = b.len;
    var r: f64[>..] = [];
    var z: f64[>..] = [];
    var p: f64[>..] = [];
    var q: f64[>..] = [];
    zeros(r, n);
    zeros(z, n);
    zeros(p, n);
    zeros(q, n);
    multiply(row_start, col, val, x, q);
    for d in n { r[d] = b[d] - q[d]; }
    let target = rtol * norm(b);
    if norm(r) <= target { return 0, true; }
    P(r[..], z[..]);
    for d in n { p[d] = z[d]; }
    var rz = dot(r, z);
    var it = 0;
    while it < max_iterations {
        multiply(row_start, col, val, p, q);
        let alpha = rz / dot(p, q);
        for d in n {
            x[d] += alpha * p[d];
            r[d] -= alpha * q[d];
        }
        it++;
        if norm(r) <= target { return it, true; }
        P(r[..], z[..]);
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    return it, false;
}

// Right-preconditioned GMRES with a restart after m iterations (Saad,
// Algorithm 9.5), with Givens rotations for the least-squares problem and
// modified Gram-Schmidt for the basis. With the right preconditioner, the
// residual that it measures is the true residual b - A x.
fn gmres<P>(row_start: i64[:], col: i32[:], val: f64[:], b: f64[:], x: f64[:],
            rtol: f64, m: i64, max_iterations: i64) -> i64, bool {
    assert(m > 0);
    let n = b.len;
    var basis: f64[>..] = [];       // m + 1 vectors of n, one after the other
    zeros(basis, (m + 1) * n);
    var h: f64[>..] = [];           // (m + 1) × m, h[(m + 1) j + i] for row i and column j
    zeros(h, (m + 1) * m);
    var cs: f64[>..] = [];
    var sn: f64[>..] = [];
    var g: f64[>..] = [];
    var y: f64[>..] = [];
    zeros(cs, m);
    zeros(sn, m);
    zeros(g, m + 1);
    zeros(y, m);
    var r: f64[>..] = [];
    var z: f64[>..] = [];
    zeros(r, n);
    zeros(z, n);
    let target = rtol * norm(b);
    var it = 0;
    while true {
        multiply(row_start, col, val, x, r);
        for d in n { r[d] = b[d] - r[d]; }
        let beta = norm(r);
        if beta <= target { return it, true; }
        if it >= max_iterations { return it, false; }
        for d in n { basis[d] = r[d] / beta; }
        for i in m + 1 { g[i] = 0.0; }
        g[0] = beta;
        var k = 0;          // the number of columns of this cycle
        while k < m && it < max_iterations {
            // w = A M⁻¹ v_k, in the slot of v_(k+1)
            let vk = basis[k * n..(k + 1) * n];
            let w = basis[(k + 1) * n..(k + 2) * n];
            P(vk, z[..]);
            multiply(row_start, col, val, z, w);
            let hk = (m + 1) * k;
            for i in k + 1 {
                let vi = basis[i * n..(i + 1) * n];
                let hik = dot(w, vi);
                h[hk + i] = hik;
                for d in n { w[d] -= hik * vi[d]; }
            }
            let wn = norm(w);
            h[hk + k + 1] = wn;
            if wn > 0.0 { for d in n { w[d] /= wn; } }
            for i in k {
                let a = h[hk + i];
                let c = h[hk + i + 1];
                h[hk + i] = cs[i] * a + sn[i] * c;
                h[hk + i + 1] = -sn[i] * a + cs[i] * c;
            }
            let a = h[hk + k];
            let c = h[hk + k + 1];
            let rho = sqrt(a * a + c * c);
            cs[k] = if rho == 0.0 { 1.0 } else { a / rho };
            sn[k] = if rho == 0.0 { 0.0 } else { c / rho };
            h[hk + k] = rho;
            h[hk + k + 1] = 0.0;
            g[k + 1] = -sn[k] * g[k];
            g[k] = cs[k] * g[k];
            k++;
            it++;
            if abs(g[k]) <= target || wn == 0.0 { break; }
        }
        // x += M⁻¹ V y, with H y = g by back substitution
        for t in k {
            let i = k - 1 - t;
            var s = g[i];
            for j in i + 1..k { s -= h[(m + 1) * j + i] * y[j]; }
            y[i] = s / h[(m + 1) * i + i];
        }
        for d in n { r[d] = 0.0; }
        for j in k {
            let vj = basis[j * n..(j + 1) * n];
            for d in n { r[d] += y[j] * vj[d]; }
        }
        P(r[..], z[..]);
        for d in n { x[d] += z[d]; }
    }
    return it, false;
}
