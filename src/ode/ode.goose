import std;
import math;
namespace ode;
import .jet;
import solvers from "../solvers";

let MAX_UNKNOWNS = 16;
let MAX_POINTS = 6;
let NEWTON_MAX = 50;
let NEWTON_TOLERANCE = 1.0e-9;
let DIFFERENCE = 1.0e-7;

var count = 0;
var step_size = 0.0;
var step_points = MAX_POINTS;
var step_count = 0;
var now = 0;
var history: f64[>..] = [];
var uj: f64[5][16] = [[0.0; 5]; 16];
var tj: f64[5] = [0.0; 5];
var resid: f64[16] = [0.0; 16];
var weights: f64[5][6] = [[0.0; 5]; 6];
var weights_points = 0;
var solves = 0;
var newton_steps = 0;
var last_seconds = 0.0;

fn begin(unknowns: i64, dt: f64, points: i64, steps: i64) {
    guard unknowns >= 1 && unknowns <= MAX_UNKNOWNS && dt > 0.0 && points >= 2 && points <= MAX_POINTS && steps >= 1 else {
        abort("begin: 1 to 16 unknowns, dt > 0, 2 to 6 points and steps >= 1");
    }
    solvers::memory_check("ode backward differences", 0, unknowns, steps + 8, 0, "");
    count = unknowns;
    step_size = dt;
    step_points = points;
    step_count = steps;
    now = 0;
    solves = 0;
    newton_steps = 0;
    weights_points = 0;
    history.clear();
    for i in unknowns * (steps + 1) { history.push(0.0); }
}

fn initial(i: i64, v: f64) {
    guard i >= 0 && i < count && now == 0 else { abort("initial: an unknown, before the first step"); }
    history[i] = v;
}

fn value(i: i64) -> f64 {
    guard i >= 0 && i < count else { abort("value: an unknown"); }
    return history[now * count + i];
}

fn time() -> f64 { return now * step_size; }

fn stencil(points: i64) {
    if weights_points == points { return; }
    for j in points {
        var poly: f64[6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        var degree = 0;
        var denominator = 1.0;
        for m in points {
            if m == j { continue; }
            var k = degree + 1;
            while k >= 0 {
                poly[k] = (if k > 0 { poly[k - 1] } else { 0.0 }) + m * poly[k];
                k--;
            }
            degree++;
            denominator *= m - j;
        }
        var scale = 1.0;
        for k in JN {
            weights[j][k] = poly[k] / denominator / scale;
            scale *= step_size;
        }
    }
    weights_points = points;
}

fn load_jets(n: i64, z: f64[:]) {
    let points = min(n, step_points - 1) + 1;
    stencil(points);
    for i in count {
        for k in JN {
            var c = weights[0][k] * z[i];
            for j in 1..points { c += weights[j][k] * history[(n - j) * count + i]; }
            uj[i][k] = c;
        }
    }
    tj = [n * step_size, 1.0, 0.0, 0.0, 0.0];
}

fn dense_solve(a: f64[:], b: f64[:]) -> bool {
    let n = count;
    for c in n {
        var best = c;
        for r in c + 1..n { if abs(a[r * n + c]) > abs(a[best * n + c]) { best = r; } }
        if abs(a[best * n + c]) < 1.0e-300 { return false; }
        if best != c {
            for k in n {
                let t = a[c * n + k];
                a[c * n + k] = a[best * n + k];
                a[best * n + k] = t;
            }
            let t = b[c];
            b[c] = b[best];
            b[best] = t;
        }
        for r in c + 1..n {
            let f = a[r * n + c] / a[c * n + c];
            for k in c..n { a[r * n + k] -= f * a[c * n + k]; }
            b[r] -= f * b[c];
        }
    }
    var c = n - 1;
    while c >= 0 {
        var s = b[c];
        for k in c + 1..n { s -= a[c * n + k] * b[k]; }
        b[c] = s / a[c * n + c];
        c--;
    }
    return true;
}

fn residual_at<RES>(n: i64, z: f64[:], r: f64[:]) {
    load_jets(n, z);
    RES(n * step_size);
    for k in count { r[k] = resid[k]; }
}

fn jacobian<RES>(n: i64, z: f64[:], r: f64[:], a: f64[:]) {
    var zj: f64[>..] = [];
    var rj: f64[>..] = [];
    for i in count {
        zj.push(z[i]);
        rj.push(0.0);
    }
    for j in count {
        let d = DIFFERENCE * (1.0 + abs(z[j]));
        zj[j] = z[j] + d;
        residual_at(n, zj, rj) { t => RES(t); };
        zj[j] = z[j];
        for i in count { a[i * count + j] = (rj[i] - r[i]) / d; }
    }
}

fn advance<RES>(n: i64) -> bool {
    guard n == now + 1 && n <= step_count else { abort("advance: the next step of the run"); }
    let start = clock();
    var z: f64[>..] = [];
    var r: f64[>..] = [];
    var a: f64[>..] = [];
    var b: f64[>..] = [];
    for i in count {
        z.push(history[(n - 1) * count + i]);
        r.push(0.0);
        b.push(0.0);
    }
    for i in count * count { a.push(0.0); }
    var converged = false;
    var it = 0;
    while it < NEWTON_MAX && !converged {
        residual_at(n, z, r) { t => RES(t); };
        jacobian(n, z, r, a) { t => RES(t); };
        for i in count { b[i] = -r[i]; }
        if !dense_solve(a, b) { return false; }
        converged = true;
        for i in count {
            z[i] += b[i];
            if abs(b[i]) > NEWTON_TOLERANCE * (1.0 + abs(z[i])) { converged = false; }
        }
        it++;
        newton_steps++;
    }
    for i in count { history[n * count + i] = z[i]; }
    now = n;
    solves++;
    last_seconds += clock() - start;
    return converged;
}
