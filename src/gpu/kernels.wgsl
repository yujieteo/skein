// The compute kernels of the GPU path of the page of `skein build --target
// wasm`. src/gpu/gpu.goose makes the buffers and the bind groups, and
// records the dispatches; each kernel does one part of a step, in f32.
//
// A dispatch of n threads has ceil(n / 64) workgroups, in rows of at most
// 65535 (src/gpu/gpu.js), so thread i is gid.x + 64 nwg.x gid.y.
//
//   central_difference  cd_drift, cd_kick (src/elements/explicit.goose)
//   ssprk3, dg          rk_stage (src/grid/solve.goose)
//   ssprk3, fv          fv_grad, fv_minmax, fv_venkat, fv_stage (src/grid/fv.goose)
//   cg, Jacobi          cg_start, cg_matvec, cg_update, cg_direction, cg_shift, reduce
//
// Each kernel binds at most 8 storage buffers, the default limit of
// WebGPU (maxStorageBuffersPerShaderStage).

const WG: u32 = 64u;

// The parameters of the kernel: counts and the codes of the settings as
// f32, exact below 2^24.
@group(0) @binding(0) var<storage, read> P: array<f32>;
// A sparse matrix in CSR.
@group(0) @binding(1) var<storage, read> rows: array<u32>;
@group(0) @binding(2) var<storage, read> cols: array<u32>;
@group(0) @binding(3) var<storage, read> vals: array<f32>;
// central_difference: the state and f, with the rows of vals and f scaled
// by 1/m. cg: the Jacobi preconditioner w.
@group(0) @binding(4) var<storage, read_write> u: array<f32>;
@group(0) @binding(5) var<storage, read_write> v: array<f32>;
@group(0) @binding(6) var<storage, read_write> acc: array<f32>;
@group(0) @binding(7) var<storage, read> w: array<f32>;
@group(0) @binding(8) var<storage, read> f: array<f32>;
// A stage of ssprk3: xo = c0 x0 + c1 (xv - dt R(xv)).
@group(0) @binding(10) var<storage, read> x0: array<f32>;
@group(0) @binding(11) var<storage, read> xv: array<f32>;
@group(0) @binding(12) var<storage, read_write> xo: array<f32>;
// fv: the value of each boundary face (4 n) and then the source of each
// cell (n), the face_bc of src/grid/mesh.goose,
// the gradient and the limiter of each cell, and min and max of xv.
@group(0) @binding(13) var<storage, read> gface: array<f32>;
@group(0) @binding(14) var<storage, read> fbc: array<i32>;
@group(0) @binding(15) var<storage, read_write> grad: array<f32>;
@group(0) @binding(16) var<storage, read_write> phi: array<f32>;
@group(0) @binding(17) var<storage, read_write> stats: array<f32>;
// cg: the vectors, the sums of each workgroup and the scalars
// (0 r·z, 1 r·r, 2 p·q, 3 the new r·z), with z = w r.
@group(0) @binding(20) var<storage, read_write> cx: array<f32>;
@group(0) @binding(21) var<storage, read_write> cr: array<f32>;
@group(0) @binding(22) var<storage, read_write> cp: array<f32>;
@group(0) @binding(23) var<storage, read_write> cq: array<f32>;
@group(0) @binding(25) var<storage, read_write> part: array<f32>;
@group(0) @binding(26) var<storage, read_write> sc: array<f32>;

var<workgroup> sh0: array<f32, WG>;
var<workgroup> sh1: array<f32, WG>;

fn thread(g: vec3u, n: vec3u) -> u32 { return g.x + g.y * n.x * WG; }

fn count() -> u32 { return u32(P[0]); }

// ---- central_difference: P = (n, dt) ----

@compute @workgroup_size(WG)
fn cd_drift(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let i = thread(g, nw);
    if (i >= count()) { return; }
    let h = P[1];
    v[i] += 0.5 * h * acc[i];
    u[i] += h * v[i];
}

@compute @workgroup_size(WG)
fn cd_kick(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let i = thread(g, nw);
    if (i >= count()) { return; }
    var s = 0.0;
    for (var e = rows[i]; e < rows[i + 1u]; e++) { s += vals[e] * u[cols[e]]; }
    acc[i] = f[i] - s;
    v[i] += 0.5 * P[1] * acc[i];
}

// ---- a stage of ssprk3 of dg: R(xv) = A xv - b, P = (n, c0, c1, dt) ----
// The rows of vals and f are scaled by the stage weight of the dof.

@compute @workgroup_size(WG)
fn rk_stage(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let i = thread(g, nw);
    if (i >= count()) { return; }
    var s = 0.0;
    for (var e = rows[i]; e < rows[i + 1u]; e++) { s += vals[e] * xv[cols[e]]; }
    xo[i] = P[1] * x0[i] + P[2] * (xv[i] - P[3] * (s - f[i]));
}

// ---- fv on the grid of nx × ny rectangles ----
// P = (n, nx, ny, hx, hy, ax, ay, d, flux, gradient, limiter, kh³, c0, c1, dt, 1/(ct V))
// The codes: flux 0 upwind, 1 linear_upwind, 2 central; gradient 0
// green_gauss, 1 least_squares; limiter 0 none, 1 minmod, 2 van_leer,
// 3 venkatakrishnan.

fn normal(s: u32) -> vec2f {
    if (s == 0u) { return vec2f(-1.0, 0.0); }
    if (s == 1u) { return vec2f(1.0, 0.0); }
    if (s == 2u) { return vec2f(0.0, -1.0); }
    return vec2f(0.0, 1.0);
}

fn side_length(s: u32) -> f32 { return select(P[3], P[4], s < 2u); }
fn side_distance(s: u32) -> f32 { return select(P[4], P[3], s < 2u); }

// The cell beside side s of cell c, or -1 on the boundary.
fn neighbour(c: u32, s: u32) -> i32 {
    let nx = u32(P[1]);
    let ny = u32(P[2]);
    let i = c % nx;
    let j = c / nx;
    if (s == 0u) { return select(-1, i32(c) - 1, i > 0u); }
    if (s == 1u) { return select(-1, i32(c) + 1, i + 1u < nx); }
    if (s == 2u) { return select(-1, i32(c) - i32(nx), j > 0u); }
    return select(-1, i32(c + nx), j + 1u < ny);
}

fn normal_velocity(s: u32) -> f32 { return dot(vec2f(P[5], P[6]), normal(s)); }

@compute @workgroup_size(WG)
fn fv_grad(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let c = thread(g, nw);
    if (c >= count()) { return; }
    let gg = P[9] == 0.0;
    var a = vec3f(0.0);
    var gr = vec2f(0.0);
    for (var s = 0u; s < 4u; s++) {
        let nb = neighbour(c, s);
        let b = fbc[4u * c + s];
        let ns = normal(s);
        if (nb < 0 && b < 0) {
            if (gg) { gr += xv[c] * ns * side_length(s); }
            continue;
        }
        let uf = select(gface[4u * c + s], xv[u32(max(nb, 0))], nb >= 0);
        let h = select(0.5 * side_distance(s), side_distance(s), nb >= 0);
        if (gg) {
            let um = select(uf, 0.5 * (xv[c] + uf), nb >= 0);
            gr += um * ns * side_length(s);
            continue;
        }
        let dx = h * ns.x;
        let dy = h * ns.y;
        let wt = 1.0 / (h * h);
        a += wt * vec3f(dx * dx, dx * dy, dy * dy);
        gr += wt * vec2f(dx, dy) * (uf - xv[c]);
    }
    var res = vec2f(0.0);
    if (gg) {
        res = gr / (P[3] * P[4]);
    } else {
        let det = a.x * a.z - a.y * a.y;
        if (det > 1.0e-12 * (a.x * a.x + a.z * a.z)) {
            res = vec2f(a.z * gr.x - a.y * gr.y, a.x * gr.y - a.y * gr.x) / det;
        } else {
            if (a.x > 0.0) { res.x = gr.x / a.x; }
            if (a.z > 0.0) { res.y = gr.y / a.z; }
        }
    }
    grad[2u * c] = res.x;
    grad[2u * c + 1u] = res.y;
}

// min and max of xv in stats, by one workgroup.
@compute @workgroup_size(WG)
fn fv_minmax(@builtin(local_invocation_index) li: u32) {
    let n = count();
    var lo = xv[0];
    var hi = xv[0];
    for (var c = li; c < n; c += WG) {
        lo = min(lo, xv[c]);
        hi = max(hi, xv[c]);
    }
    sh0[li] = lo;
    sh1[li] = hi;
    workgroupBarrier();
    for (var k = WG / 2u; k > 0u; k /= 2u) {
        if (li < k) {
            sh0[li] = min(sh0[li], sh0[li + k]);
            sh1[li] = max(sh1[li], sh1[li + k]);
        }
        workgroupBarrier();
    }
    if (li == 0u) {
        stats[0] = sh0[0];
        stats[1] = sh1[0];
    }
}

// The limiter of Venkatakrishnan of each cell.
@compute @workgroup_size(WG)
fn fv_venkat(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let c = thread(g, nw);
    if (c >= count()) { return; }
    let e2 = max(P[11] * (stats[1] - stats[0]) * (stats[1] - stats[0]), 1.0e-30);
    var umax = xv[c];
    var umin = xv[c];
    for (var s = 0u; s < 4u; s++) {
        let nb = neighbour(c, s);
        if (nb < 0 && fbc[4u * c + s] < 0) { continue; }
        let val = select(2.0 * gface[4u * c + s] - xv[c], xv[u32(max(nb, 0))], nb >= 0);
        umax = max(umax, val);
        umin = min(umin, val);
    }
    var fl = 1.0;
    let gc = vec2f(grad[2u * c], grad[2u * c + 1u]);
    for (var s = 0u; s < 4u; s++) {
        let d2 = 0.5 * side_distance(s) * dot(gc, normal(s));
        if (d2 == 0.0) { continue; }
        let d1 = select(umin - xv[c], umax - xv[c], d2 > 0.0);
        fl = min(fl, (d1 * d1 + e2 + 2.0 * d1 * d2) / (d1 * d1 + 2.0 * d2 * d2 + d1 * d2 + e2));
    }
    phi[c] = fl;
}

fn psi(r: f32) -> f32 {
    let lim = P[10];
    if (lim == 1.0) { return max(0.0, min(1.0, r)); }
    if (lim == 2.0) { return (r + abs(r)) / (1.0 + abs(r)); }
    return 0.5 * (1.0 + r);
}

// The advected value on the face of the upwind cell c, with the value ud
// downwind at the distance h along normal s.
fn upwind_value(c: u32, s: u32, ud: f32, h: f32) -> f32 {
    let fl = P[8];
    if (fl == 0.0) { return xv[c]; }
    if (fl == 2.0) { return 0.5 * (xv[c] + ud); }
    let gs = h * dot(vec2f(grad[2u * c], grad[2u * c + 1u]), normal(s));
    if (P[10] == 3.0) { return xv[c] + 0.5 * phi[c] * gs; }
    if (P[10] == 0.0) { return xv[c] + 0.5 * gs; }
    let du = ud - xv[c];
    if (abs(du) <= 1.0e-14 * (abs(xv[c]) + abs(ud))) { return xv[c]; }
    return xv[c] + 0.5 * psi(2.0 * gs / du - 1.0) * du;
}

// The residual of cell c (the outward flux less the source) and the stage.
@compute @workgroup_size(WG)
fn fv_stage(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let c = thread(g, nw);
    if (c >= count()) { return; }
    let d = P[7];
    var r = -gface[4u * count() + c];
    for (var s = 0u; s < 4u; s++) {
        let nb = neighbour(c, s);
        let an = normal_velocity(s);
        let h = side_distance(s);
        if (nb >= 0) {
            let m = u32(nb);
            let uf = select(upwind_value(m, s ^ 1u, xv[c], h), upwind_value(c, s, xv[m], h), an >= 0.0);
            r += (an * uf - d * (xv[m] - xv[c]) / h) * side_length(s);
        } else if (fbc[4u * c + s] >= 0) {
            let gv = gface[4u * c + s];
            let uf = select(gv, upwind_value(c, s, 2.0 * gv - xv[c], h), an >= 0.0);
            r += (an * uf - d * (gv - xv[c]) / (0.5 * h)) * side_length(s);
        }
    }
    xo[c] = P[12] * x0[c] + P[13] * (xv[c] - P[14] * r * P[15]);
}

// ---- cg with the Jacobi preconditioner w (0 on the fixed dofs): P = (n) ----

// The sums of a and b over the workgroup, in part.
fn sums(li: u32, wid: u32, a: f32, b: f32) {
    sh0[li] = a;
    sh1[li] = b;
    workgroupBarrier();
    for (var k = WG / 2u; k > 0u; k /= 2u) {
        if (li < k) {
            sh0[li] += sh0[li + k];
            sh1[li] += sh1[li + k];
        }
        workgroupBarrier();
    }
    if (li == 0u) {
        part[2u * wid] = sh0[0];
        part[2u * wid + 1u] = sh1[0];
    }
}

@compute @workgroup_size(WG)
fn cg_start(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u,
            @builtin(local_invocation_index) li: u32, @builtin(workgroup_id) wg: vec3u) {
    let i = thread(g, nw);
    var a = 0.0;
    var b = 0.0;
    if (i < count()) {
        cp[i] = w[i] * cr[i];
        a = cr[i] * cp[i];
        b = cr[i] * cr[i];
    }
    sums(li, wg.x + wg.y * nw.x, a, b);
}

@compute @workgroup_size(WG)
fn cg_matvec(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u,
             @builtin(local_invocation_index) li: u32, @builtin(workgroup_id) wg: vec3u) {
    let i = thread(g, nw);
    var a = 0.0;
    if (i < count()) {
        var s = 0.0;
        if (w[i] > 0.0) {
            for (var e = rows[i]; e < rows[i + 1u]; e++) { s += vals[e] * cp[cols[e]]; }
        }
        cq[i] = s;
        a = cp[i] * s;
    }
    sums(li, wg.x + wg.y * nw.x, a, 0.0);
}

@compute @workgroup_size(WG)
fn cg_update(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u,
             @builtin(local_invocation_index) li: u32, @builtin(workgroup_id) wg: vec3u) {
    let i = thread(g, nw);
    let alpha = select(0.0, sc[0] / sc[2], sc[2] != 0.0);
    var a = 0.0;
    var b = 0.0;
    if (i < count()) {
        cx[i] += alpha * cp[i];
        cr[i] -= alpha * cq[i];
        a = w[i] * cr[i] * cr[i];
        b = cr[i] * cr[i];
    }
    sums(li, wg.x + wg.y * nw.x, a, b);
}

@compute @workgroup_size(WG)
fn cg_direction(@builtin(global_invocation_id) g: vec3u, @builtin(num_workgroups) nw: vec3u) {
    let i = thread(g, nw);
    if (i >= count()) { return; }
    let beta = select(0.0, sc[3] / sc[0], sc[0] != 0.0);
    cp[i] = w[i] * cr[i] + beta * cp[i];
}

@compute @workgroup_size(1)
fn cg_shift() { sc[0] = sc[3]; }

// The sums of the P[0] workgroups in part, by one workgroup, in sc[P[1]]
// and sc[P[2]].
@compute @workgroup_size(WG)
fn reduce(@builtin(local_invocation_index) li: u32) {
    let groups = count();
    var a = 0.0;
    var b = 0.0;
    for (var k = li; k < groups; k += WG) {
        a += part[2u * k];
        b += part[2u * k + 1u];
    }
    sh0[li] = a;
    sh1[li] = b;
    workgroupBarrier();
    for (var k = WG / 2u; k > 0u; k /= 2u) {
        if (li < k) {
            sh0[li] += sh0[li + k];
            sh1[li] += sh1[li + k];
        }
        workgroupBarrier();
    }
    if (li == 0u) {
        sc[u32(P[1])] = sh0[0];
        sc[u32(P[2])] = sh1[0];
    }
}
