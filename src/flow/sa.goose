import std;
import math;
namespace flow;
import .mesh;
import .cell;
import .rans;

var sa_sigma = 0.66666;
var sa_kappa = 0.41;
var sa_cb1 = 0.1355;
var sa_cb2 = 0.622;
var sa_cw2 = 0.3;
var sa_cw3 = 2.0;
var sa_cv1 = 7.1;
var sa_cs = 0.3;
var sa_cw1 = 0.1355 / (0.41 * 0.41) + 1.622 / 0.66666;

var sa_start_v: f64[2] = [0.0, 0.0];
var sa_start_nt = 0.0;

fn spalart_start(vx: f64, vy: f64, nt: f64) {
    sa_start_v = [vx, vy];
    sa_start_nt = nt;
}

let N0 = 9;

fn spalart_constants(sigma: f64, kappa: f64, cb1: f64, cb2: f64, cw1: f64, cw2: f64, cw3: f64, cv1: f64, cs: f64) {
    guard sigma > 0.0 && kappa > 0.0 && cb1 > 0.0 && cw3 > 0.0 && cv1 > 0.0 else { abort("spalart_constants: sigma, kappa, cb1, cw3 and cv1 > 0"); }
    sa_sigma = sigma;
    sa_kappa = kappa;
    sa_cb1 = cb1;
    sa_cb2 = cb2;
    sa_cw2 = cw2;
    sa_cw3 = cw3;
    sa_cv1 = cv1;
    sa_cs = cs;
    sa_cw1 = cw1;
}

struct Spalart {
    nut: f64,
    s_t: f64,
    fw: f64,
}

fn spalart(nt: f64, d: f64, omega: f64) -> Spalart {
    let nu = mu / rho;
    let n = max(nt, 0.0);
    let dp = max(d, 1.0e-12 * diameter);
    let chi = n / nu;
    let c3 = chi * chi * chi;
    let fv1 = c3 / (c3 + sa_cv1 * sa_cv1 * sa_cv1);
    let fv2 = 1.0 - chi / (1.0 + chi * fv1);
    let k2d2 = sa_kappa * sa_kappa * dp * dp;
    let s_t = max(omega + fv2 * n / k2d2, sa_cs * omega);
    let r = min(n / (max(s_t, 1.0e-30) * k2d2), 10.0);
    let r3 = r * r * r;
    let g = r + sa_cw2 * (r3 * r3 - r);
    let g6 = g * g * g * g * g * g;
    let c6 = sa_cw3 * sa_cw3 * sa_cw3 * sa_cw3 * sa_cw3 * sa_cw3;
    let fw = g * pow((1.0 + c6) / (g6 + c6), 1.0 / 6.0);
    return Spalart { nut: n * fv1, s_t: s_t, fw: fw };
}

fn sa_residual(c: i64, u: f64[15], r: f64[15]&) { sa_rows(c, u, u, r); }

fn sa_rows(c: i64, u: f64[15], ub: f64[15], r: f64[15]&) {
    guard turbulent && spalart_on && order == 1 && c >= 0 && c < cell_count() && wall_d.len == vertex_count else { abort("sa_residual: a Spalart-Allmaras P1 flow after clear_bcs"); }
    let g = geometry(c);
    for d in 15 { r[d] = 0.0; }
    var gv: f64[4] = [0.0; 4];
    var gn: f64[2] = [0.0; 2];
    var gp: f64[2] = [0.0; 2];
    for a in 3 {
        for j in 2 {
            gv[j] += g.gl[2 * a + j] * u[2 * a];
            gv[2 + j] += g.gl[2 * a + j] * u[2 * a + 1];
            gp[j] += g.gl[2 * a + j] * u[6 + a];
            gn[j] += g.gl[2 * a + j] * u[N0 + a];
        }
    }
    let omega = abs(gv[1] - gv[2]);
    let div = gv[0] + gv[3];
    for q in 7 {
        let l, wq = quad_point(q);
        sa_point(c, g, l, wq * g.area, u, ub, gv, gp, gn, omega, div, r);
    }
}

fn sa_point(c: i64, g: Geometry, l: f64[3], w: f64, u: f64[15], ub: f64[15], gv: f64[4], gp: f64[2], gn: f64[2],
            omega: f64, div: f64, r: f64[15]&) {
    var v: f64[2] = [0.0; 2];
    var p = 0.0;
    var nt = 0.0;
    var d = 0.0;
    for a in 3 {
        v[0] += l[a] * u[2 * a];
        v[1] += l[a] * u[2 * a + 1];
        p += l[a] * u[6 + a];
        nt += l[a] * u[N0 + a];
        d += l[a] * wall_d[cell_node(c, a)];
    }
    let t = spalart(nt, d, omega);
    let nu = mu / rho;
    let nu_eff = nu + t.nut;
    let speed2 = v[0] * v[0] + v[1] * v[1];
    let tm = supg_tau(g.h, speed2, nu_eff, 0.0);
    var res: f64[2] = [0.0; 2];
    for i in 2 { res[i] = rho * (gv[2 * i] * v[0] + gv[2 * i + 1] * v[1]) + gp[i] - force[i]; }
    let n = max(nt, 0.0);
    let dp = max(d, 1.0e-12 * diameter);
    let gn2 = gn[0] * gn[0] + gn[1] * gn[1];
    let pb = production(c, g, l, ub);
    let source = pb - sa_cw1 * t.fw * n * n / (dp * dp);
    let rn = v[0] * gn[0] + v[1] * gn[1] - sa_cb2 / sa_sigma * gn2 - source;
    let visc = (nu + n) / sa_sigma;
    let tn = supg_tau(g.h, speed2, visc, 0.0);
    let tc = 1.0 / (tm * g.trace_g);
    for a in 3 {
        let gx = g.gl[2 * a];
        let gy = g.gl[2 * a + 1];
        let streamline = v[0] * gx + v[1] * gy;
        for i in 2 {
            var row = l[a] * (res[i] - gp[i]) - p * g.gl[2 * a + i];
            row += rho * nu_eff * (gx * (gv[2 * i] + gv[i]) + gy * (gv[2 * i + 1] + gv[2 + i]));
            row += tm * streamline * res[i] + tc * rho * g.gl[2 * a + i] * div;
            r[2 * a + i] += w * row;
        }
        r[6 + a] += w * (l[a] * div + tm * (gx * res[0] + gy * res[1]) / rho);
        r[N0 + a] += w * (l[a] * (v[0] * gn[0] + v[1] * gn[1] - sa_cb2 / sa_sigma * gn2 - source)
                          + visc * (gx * gn[0] + gy * gn[1]) + tn * streamline * rn);
    }
}

fn production(c: i64, g: Geometry, l: f64[3], u: f64[15]) -> f64 {
    var gv: f64[4] = [0.0; 4];
    var nt = 0.0;
    var d = 0.0;
    for a in 3 {
        for j in 2 {
            gv[j] += g.gl[2 * a + j] * u[2 * a];
            gv[2 + j] += g.gl[2 * a + j] * u[2 * a + 1];
        }
        nt += l[a] * u[N0 + a];
        d += l[a] * wall_d[cell_node(c, a)];
    }
    let t = spalart(nt, d, abs(gv[1] - gv[2]));
    return sa_cb1 * t.s_t * max(nt, 0.0);
}

fn sa_jacobian(c: i64, u: f64[15], r: f64[15]&, k: f64[225]&) {
    sa_residual(c, u, r);
    var shifted = u;
    var rs: f64[15] = [0.0; 15];
    for b in 12 {
        let step = 1.0e-7 * max(1.0, abs(u[b]));
        shifted[b] = u[b] + step;
        sa_rows(c, shifted, u, rs);
        shifted[b] = u[b];
        for a in 12 { k[15 * a + b] = (rs[a] - r[a]) / step; }
    }
}
