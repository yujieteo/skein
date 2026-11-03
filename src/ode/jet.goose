// Truncated Taylor series in time ("jets") of the package ode. A jet is
// the coefficients c_k = f^(k)(t) / k! of a function at one time, k = 0 to
// JN - 1. The arithmetic of two jets is the arithmetic of the series, so
// `ddt` of any expression of the unknowns is exact. A jet loses its last
// coefficient at each ddt: a model uses at most JN - 1 nested ddt.
import std;
import math;
namespace ode;

let JN = 5;

fn jconst(c: f64) -> f64[5] { return [c, 0.0, 0.0, 0.0, 0.0]; }

fn jadd(a: f64[5], b: f64[5]) -> f64[5] { return [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4]]; }

fn jsub(a: f64[5], b: f64[5]) -> f64[5] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3], a[4] - b[4]]; }

fn jneg(a: f64[5]) -> f64[5] { return [-a[0], -a[1], -a[2], -a[3], -a[4]]; }

fn jmul(a: f64[5], b: f64[5]) -> f64[5] {
    var r: f64[5] = [0.0; 5];
    for k in JN {
        for i in k + 1 { r[k] += a[i] * b[k - i]; }
    }
    return r;
}

// a / b: the series q with q b = a.
fn jdiv(a: f64[5], b: f64[5]) -> f64[5] {
    var q: f64[5] = [0.0; 5];
    for k in JN {
        var s = a[k];
        for i in 1..k + 1 { s -= b[i] * q[k - i]; }
        q[k] = s / b[0];
    }
    return q;
}

// d/dt: c'_k = (k + 1) c_(k+1). The top coefficient is lost.
fn jddt(a: f64[5]) -> f64[5] { return [a[1], 2.0 * a[2], 3.0 * a[3], 4.0 * a[4], 0.0]; }

// sin and cos of a: s' = a' c and c' = -a' s.
fn jsin(a: f64[5]) -> f64[5] {
    var s: f64[5] = [sin(a[0]), 0.0, 0.0, 0.0, 0.0];
    var c: f64[5] = [cos(a[0]), 0.0, 0.0, 0.0, 0.0];
    for k in 1..JN {
        var ds = 0.0;
        var dc = 0.0;
        for i in 1..k + 1 {
            ds += i * a[i] * c[k - i];
            dc -= i * a[i] * s[k - i];
        }
        s[k] = ds / k;
        c[k] = dc / k;
    }
    return s;
}

fn jcos(a: f64[5]) -> f64[5] {
    var s: f64[5] = [sin(a[0]), 0.0, 0.0, 0.0, 0.0];
    var c: f64[5] = [cos(a[0]), 0.0, 0.0, 0.0, 0.0];
    for k in 1..JN {
        var ds = 0.0;
        var dc = 0.0;
        for i in 1..k + 1 {
            ds += i * a[i] * c[k - i];
            dc -= i * a[i] * s[k - i];
        }
        s[k] = ds / k;
        c[k] = dc / k;
    }
    return c;
}

// exp(a): e' = a' e.
fn jexp(a: f64[5]) -> f64[5] {
    var e: f64[5] = [exp(a[0]), 0.0, 0.0, 0.0, 0.0];
    for k in 1..JN {
        var d = 0.0;
        for i in 1..k + 1 { d += i * a[i] * e[k - i]; }
        e[k] = d / k;
    }
    return e;
}

// ln(a): a l' = a'.
fn jlog(a: f64[5]) -> f64[5] {
    var l: f64[5] = [log(a[0]), 0.0, 0.0, 0.0, 0.0];
    for k in 1..JN {
        var d = k * a[k];
        for i in 1..k { d -= i * l[i] * a[k - i]; }
        l[k] = d / (k * a[0]);
    }
    return l;
}

// sqrt(a): s s = a.
fn jsqrt(a: f64[5]) -> f64[5] {
    var s: f64[5] = [sqrt(a[0]), 0.0, 0.0, 0.0, 0.0];
    for k in 1..JN {
        var d = a[k];
        for i in 1..k { d -= s[i] * s[k - i]; }
        s[k] = d / (2.0 * s[0]);
    }
    return s;
}

// a to the whole power n.
fn jpow(a: f64[5], n: i64) -> f64[5] {
    var r = jconst(1.0);
    for i in abs(n) { r = jmul(r, a); }
    return if n < 0 { jdiv(jconst(1.0), r) } else { r };
}
