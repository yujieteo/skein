import std;
import math;
namespace solvers;

fn band_widths(row_start: i64[:], col: i32[:]) -> i64, i64 {
    guard row_start.len > 0 && row_start[row_start.len - 1] == col.len else { abort("band_widths: row_start ends at col.len"); }
    var kl = 0;
    var ku = 0;
    for r in row_start.len - 1 {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            kl = max(kl, r - c);
            ku = max(ku, c - r);
        }
    }
    return kl, ku;
}

fn band_size(n: i64, kl: i64, ku: i64) -> i64 { return n * (2 * kl + ku + 1); }

fn band_solve(row_start: i64[:], col: i32[:], val: f64[:], b: f64[:], x: f64[:], kl: i64, ku: i64, band: f64[:]) -> bool {
    let n = row_start.len - 1;
    let w = 2 * kl + ku + 1;
    guard kl >= 0 && ku >= 0 && band.len == n * w && b.len == n && x.len == n && col.len == val.len else {
        abort("band_solve: band has band_size entries, and b and x have one entry for each row");
    }
    for i in band.len { band[i] = 0.0; }
    for r in n {
        for e in row_start[r]..row_start[r + 1] { band[r * w + (col[e] - r + kl)] = val[e]; }
        x[r] = b[r];
    }
    for k in n {
        let last = min(n - 1, k + kl);
        var p = k;
        for i in k + 1..last + 1 {
            if abs(band[i * w + (k - i + kl)]) > abs(band[p * w + (k - p + kl)]) { p = i; }
        }
        if band[p * w + (k - p + kl)] == 0.0 { return false; }
        if p != k { swap_rows(band, w, k, p, kl, ku, x); }
        let pivot = band[k * w + kl];
        let width = min(n - 1, k + ku + kl);
        for i in k + 1..last + 1 {
            let f = band[i * w + (k - i + kl)] / pivot;
            if f == 0.0 { continue; }
            for c in k..width + 1 { band[i * w + (c - i + kl)] -= f * band[k * w + (c - k + kl)]; }
            x[i] -= f * x[k];
        }
    }
    for kk in n {
        let k = n - 1 - kk;
        var s = x[k];
        for c in k + 1..min(n - 1, k + ku + kl) + 1 { s -= band[k * w + (c - k + kl)] * x[c]; }
        x[k] = s / band[k * w + kl];
    }
    return true;
}

fn swap_rows(band: f64[:], w: i64, k: i64, p: i64, kl: i64, ku: i64, x: f64[:]) {
    let n = x.len;
    for c in k..min(n - 1, k + ku + kl) + 1 {
        let a = k * w + (c - k + kl);
        let b = p * w + (c - p + kl);
        let t = band[a];
        band[a] = band[b];
        band[b] = t;
    }
    let t = x[k];
    x[k] = x[p];
    x[p] = t;
}
