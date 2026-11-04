import std;
import math;
namespace flow;
import .mesh;
import .cell;
import .rans;

var conductivity = 1.0;
var capacity = 1.0;
var source = 0.0;
var solid_conductivity = 1.0;
var solid_source = 0.0;
var t_start = 0.0;
var edge_load: f64[>..] = [];

fn thermal_fluid(k: f64, rho_c: f64, q: f64) {
    guard k > 0.0 && rho_c > 0.0 else { abort("thermal_fluid: k > 0 and rho c > 0"); }
    conductivity = k;
    capacity = rho_c;
    source = q;
}

fn thermal_solid(k: f64, q: f64) {
    guard k > 0.0 else { abort("thermal_solid: k > 0"); }
    solid_conductivity = k;
    solid_source = q;
}

fn buoyant_force(fx: f64, fy: f64) {
    buoyancy[0] = fx;
    buoyancy[1] = fy;
}

fn temperature_start(t: f64) { t_start = t; }

fn thermal_residual(c: i64, u: f64[21], r: f64[21]&) {
    guard thermal && c >= 0 && c < cell_count() else { abort("thermal_residual: a thermal flow and c < cell_count()"); }
    for d in MAX_CELL_DOFS { r[d] = 0.0; }
    let n = 2 * per_cell + 3;
    var tl: f64[6] = [0.0; 6];
    for a in per_cell { tl[a] = u[n + a]; }
    let g = geometry(c);
    if is_solid_cell(c) {
        solid_rows(g, tl, n, r);
        return;
    }
    var u15: f64[15] = [0.0; 15];
    var r15: f64[15] = [0.0; 15];
    for d in n { u15[d] = u[d]; }
    cell_residual_t(c, u15, tl, r15);
    for d in n { r[d] = r15[d]; }
    fluid_rows(g, u, tl, n, r);
}

fn solid_rows(g: Geometry, tl: f64[6], n: i64, r: f64[21]&) {
    for q in 7 {
        let l, wq = quad_point(q);
        let w = wq * g.area;
        let s = shape(g, l);
        var gt: f64[2] = [0.0; 2];
        for a in per_cell {
            for j in 2 { gt[j] += s.dn[2 * a + j] * tl[a]; }
        }
        for a in per_cell {
            r[n + a] += w / capacity * (solid_conductivity * (s.dn[2 * a] * gt[0] + s.dn[2 * a + 1] * gt[1]) - s.n[a] * solid_source);
        }
    }
}

fn fluid_rows(g: Geometry, u: f64[21], tl: f64[6], n: i64, r: f64[21]&) {
    let alpha = conductivity / capacity;
    for q in 7 {
        let l, wq = quad_point(q);
        let w = wq * g.area;
        let s = shape(g, l);
        var v: f64[2] = [0.0; 2];
        var gt: f64[2] = [0.0; 2];
        var lap = 0.0;
        for a in per_cell {
            for i in 2 {
                v[i] += s.n[a] * u[2 * a + i];
                gt[i] += s.dn[2 * a + i] * tl[a];
            }
            lap += (s.hs[3 * a] + s.hs[3 * a + 2]) * tl[a];
        }
        let advect = v[0] * gt[0] + v[1] * gt[1];
        let tau = supg_tau(g.h, v[0] * v[0] + v[1] * v[1], alpha, 0.0);
        let res = advect - alpha * lap - source / capacity;
        for a in per_cell {
            let streamline = v[0] * s.dn[2 * a] + v[1] * s.dn[2 * a + 1];
            let diffuse = alpha * (s.dn[2 * a] * gt[0] + s.dn[2 * a + 1] * gt[1]);
            r[n + a] += w * (s.n[a] * (advect - source / capacity) + diffuse + tau * streamline * res);
        }
    }
}

fn thermal_jacobian(c: i64, u: f64[21], r: f64[21]&, k: f64[441]&) {
    thermal_residual(c, u, r);
    let m = cell_dof_count();
    var shifted = u;
    var rs: f64[21] = [0.0; 21];
    for b in m {
        let step = 1.0e-7 * max(1.0, abs(u[b]));
        shifted[b] = u[b] + step;
        thermal_residual(c, shifted, rs);
        shifted[b] = u[b];
        for a in m { k[21 * a + b] = (rs[a] - r[a]) / step; }
    }
}
