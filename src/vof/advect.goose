import std;
import math;
namespace vof;
import .mesh;

let HRIC_LOW = 0.3;
let HRIC_HIGH = 0.7;

var af: f64[>..] = [];
var cour: f64[>..] = [];
var alpha_new: f64[>..] = [];
var alpha_s: f64[>..] = [];
var normal: f64[>..] = [];
var kappa: f64[>..] = [];

fn face_cells(f: i64) -> i64[4] {
    if f < x_face_count() {
        let i = f % (nx + 1);
        let j = f / (nx + 1);
        return [if i > 1 { i - 2 + nx * j } else { -1 }, if i > 0 { i - 1 + nx * j } else { -1 },
                if i < nx { i + nx * j } else { -1 }, if i + 1 < nx { i + 1 + nx * j } else { -1 }];
    }
    let k = f - x_face_count();
    let i = k % nx;
    let j = k / nx;
    return [if j > 1 { i + nx * (j - 2) } else { -1 }, if j > 0 { i + nx * (j - 1) } else { -1 },
            if j < ny { i + nx * j } else { -1 }, if j + 1 < ny { i + nx * (j + 1) } else { -1 }];
}

fn face_axis(f: i64) -> i64 { return if f < x_face_count() { 0 } else { 1 }; }

fn face_area(f: i64) -> f64 { return if f < x_face_count() { hy } else { hx }; }

fn cell_volume() -> f64 { return hx * hy; }

fn beside(q: f64[:], c: i64, s: i64) -> f64 {
    let n = neighbour(c, s);
    return if n >= 0 { q[n] } else { q[c] };
}

fn cell_gradient(q: f64[:], c: i64) -> f64[2] {
    let gx = (beside(q, c, 1) - beside(q, c, 0)) / (if neighbour(c, 0) >= 0 && neighbour(c, 1) >= 0 { 2.0 * hx } else { hx });
    let gy = (beside(q, c, 3) - beside(q, c, 2)) / (if neighbour(c, 2) >= 0 && neighbour(c, 3) >= 0 { 2.0 * hy } else { hy });
    return [gx, gy];
}

fn courant(dt: f64) -> f64 {
    var most = 0.0;
    for c in cell_count() {
        var out = 0.0;
        for s in 4 {
            let f = face_of(c, s);
            out += max(0.0, outward(s) * uf[f]) * face_area(f);
        }
        cour[c] = dt * out / cell_volume();
        most = max(most, cour[c]);
    }
    return most;
}

fn face_alpha(f: i64) -> f64 {
    let fc = face_cells(f);
    let up = uf[f] >= 0.0;
    let d = if up { fc[1] } else { fc[2] };
    let a = if up { fc[2] } else { fc[1] };
    let far = if up { fc[0] } else { fc[3] };
    if d < 0 { return alpha[a]; }
    if a < 0 || far < 0 { return alpha[d]; }
    let den = alpha[a] - alpha[far];
    if abs(den) < 1.0e-12 { return alpha[d]; }
    let c = (alpha[d] - alpha[far]) / den;
    if c <= 0.0 || c >= 1.0 { return alpha[d]; }
    let g = cell_gradient(alpha, d);
    let size = sqrt(g[0] * g[0] + g[1] * g[1]);
    let cos_theta = if size > 1.0e-12 { abs(g[face_axis(f)]) / size } else { 1.0 };
    let w = sqrt(cos_theta);
    var cf = (if c < 0.5 { 2.0 * c } else { 1.0 }) * w + c * (1.0 - w);
    if cour[d] >= HRIC_HIGH { cf = c; }
    else if cour[d] > HRIC_LOW { cf = c + (cf - c) * (HRIC_HIGH - cour[d]) / (HRIC_HIGH - HRIC_LOW); }
    return alpha[far] + cf * den;
}

fn face_alphas() {
    for f in face_count() {
        let fc = face_cells(f);
        af[f] = if fc[1] < 0 || fc[2] < 0 { alpha[if fc[1] < 0 { fc[2] } else { fc[1] }] } else { face_alpha(f) };
    }
}

fn advect(dt: f64) {
    face_alphas();
    for c in cell_count() {
        var sum = 0.0;
        var total = 0.0;
        for s in 4 {
            let f = face_of(c, s);
            let flux = outward(s) * uf[f] * face_area(f);
            sum += flux * af[f];
            total += flux;
        }
        alpha_new[c] = clamp(alpha[c] - dt * (sum - alpha[c] * total) / cell_volume(), 0.0, 1.0);
    }
    for c in cell_count() { alpha[c] = alpha_new[c]; }
}

fn smooth_alpha() {
    for c in cell_count() {
        var sum = 0.0;
        for dj in 3 {
            for di in 3 {
                let i = clamp(c % nx + di - 1, 0, nx - 1);
                let j = clamp(c / nx + dj - 1, 0, ny - 1);
                let w = (if di == 1 { 2.0 } else { 1.0 }) * (if dj == 1 { 2.0 } else { 1.0 });
                sum += w * alpha[i + nx * j];
            }
        }
        alpha_s[c] = sum / 16.0;
    }
}

fn curvature() {
    smooth_alpha();
    for c in cell_count() {
        let g = cell_gradient(alpha_s, c);
        let size = sqrt(g[0] * g[0] + g[1] * g[1]);
        let on = size > 1.0e-6 / min(hx, hy);
        normal[2 * c] = if on { g[0] / size } else { 0.0 };
        normal[2 * c + 1] = if on { g[1] / size } else { 0.0 };
    }
    for c in cell_count() {
        let ex = beside_component(c, 1, 0) - beside_component(c, 0, 0);
        let ey = beside_component(c, 3, 1) - beside_component(c, 2, 1);
        let dx = if neighbour(c, 0) >= 0 && neighbour(c, 1) >= 0 { 2.0 * hx } else { hx };
        let dy = if neighbour(c, 2) >= 0 && neighbour(c, 3) >= 0 { 2.0 * hy } else { hy };
        kappa[c] = -(ex / dx + ey / dy);
    }
}

fn beside_component(c: i64, s: i64, k: i64) -> f64 {
    let n = neighbour(c, s);
    return normal[2 * (if n >= 0 { n } else { c }) + k];
}
