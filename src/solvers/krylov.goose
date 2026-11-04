import std;
import math;
namespace solvers;
import .sparse;
import .pool;

fn rotate(hk: f64[:], cs: f64[:], sn: f64[:], g: f64[:], k: i64) {
    guard k >= 0 && hk.len == k + 2 && cs.len > k && sn.len > k && g.len > k + 1 else { abort("rotate: hk has k + 2 entries, cs and sn more than k, g more than k + 1"); }
    for i in k {
        let a = hk[i];
        let c = hk[i + 1];
        hk[i] = cs[i] * a + sn[i] * c;
        hk[i + 1] = -sn[i] * a + cs[i] * c;
    }
    let a = hk[k];
    let c = hk[k + 1];
    let rho = sqrt(a * a + c * c);
    cs[k] = if rho == 0.0 { 1.0 } else { a / rho };
    sn[k] = if rho == 0.0 { 0.0 } else { c / rho };
    hk[k] = rho;
    hk[k + 1] = 0.0;
    g[k + 1] = -sn[k] * g[k];
    g[k] = cs[k] * g[k];
}

fn combine(h: f64[:], g: f64[:], y: f64[:], basis: f64[:], r: f64[:], m: i64, k: i64) {
    guard k >= 0 && k <= m && h.len == (m + 1) * m && g.len >= k && y.len >= k && basis.len >= k * r.len else { abort("combine: 0 <= k <= m, h is (m + 1) x m, g and y have k entries, basis k vectors"); }
    let n = r.len;
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
}

fn gmres<P>(row_start: i64[:], col: i32[:], val: f64[:], b: f64[:], x: f64[:],
            rtol: f64, m: i64, max_iterations: i64) -> i64, bool {
    guard m > 0 && x.len == b.len && row_start.len == b.len + 1 && col.len == val.len else { abort("gmres: m > 0, x and b one size, row_start b.len + 1 entries, col and val one size"); }
    let n = b.len;
    var basis: f64[>..] = [];
    zeros(basis, (m + 1) * n);
    var h: f64[>..] = [];
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
        rows_multiply(row_start, col, val, x, r);
        for d in n { r[d] = b[d] - r[d]; }
        let beta = norm(r);
        if beta <= target { return it, true; }
        if it >= max_iterations { return it, false; }
        for d in n { basis[d] = r[d] / beta; }
        for i in m + 1 { g[i] = 0.0; }
        g[0] = beta;
        var k = 0;
        while k < m && it < max_iterations {
            let vk = basis[k * n..(k + 1) * n];
            let w = basis[(k + 1) * n..(k + 2) * n];
            P(vk, z[..]);
            rows_multiply(row_start, col, val, z, w);
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
            rotate(h[hk..hk + k + 2], cs, sn, g, k);
            k++;
            it++;
            if abs(g[k]) <= target || wn == 0.0 { break; }
        }
        combine(h, g, y, basis, r, m, k);
        P(r[..], z[..]);
        for d in n { x[d] += z[d]; }
    }
    return it, false;
}
