import std;
import math;
import os;
import solvers from "../solvers";
namespace vof;
import .mesh;
import .advect;
import .step;

fn boundary_cell(f: i64) -> i64, i64 {
    let fc = face_cells(f);
    let x = face_axis(f) == 0;
    if fc[2] < 0 { return fc[1], if x { 1 } else { 3 }; }
    return fc[2], if x { 0 } else { 2 };
}

fn face_fluxes(dt: f64) {
    for f in face_count() {
        let fc = face_cells(f);
        let ax = face_axis(f);
        let h = if ax == 0 { hx } else { hy };
        if fc[1] >= 0 && fc[2] >= 0 {
            let rf = 0.5 * (rho_c[fc[1]] + rho_c[fc[2]]);
            let kf = 0.5 * (kappa[fc[1]] + kappa[fc[2]]);
            let force = gravity[ax] + tension * kf * (alpha[fc[2]] - alpha[fc[1]]) / (h * rf);
            ui[f] = 0.5 * (vstar[2 * fc[1] + ax] + vstar[2 * fc[2] + ax]);
            ustar[f] = ui[f] + dt * (force - (pres[fc[2]] - pres[fc[1]]) / (h * rf));
            coef[f] = face_area(f) * dt / (rf * h);
            continue;
        }
        let e, s = boundary_cell(f);
        guard side_kind[4 * e + s] == OPEN else { continue; };
        ui[f] = vstar[2 * e + ax];
        ustar[f] = ui[f] + dt * (gravity[ax] + outward(s) * pres[e] / (0.5 * h * rho_c[e]));
        coef[f] = face_area(f) * dt / (rho_c[e] * 0.5 * h);
    }
}

fn flux_out(c: i64, s: i64) -> f64 {
    let f = face_of(c, s);
    return outward(s) * (if side_kind[4 * c + s] >= 0 { uf[f] } else { ustar[f] }) * face_area(f);
}

fn pressure_matrix() {
    var open = false;
    for c in cell_count() {
        var e = row_start[c];
        var at = -1;
        var sum = 0.0;
        var d = 0.0;
        for s in [2, 0, 4, 1, 3] {
            if s == 4 {
                at = e;
                e++;
                continue;
            }
            let f = face_of(c, s);
            sum -= flux_out(c, s);
            if neighbour(c, s) >= 0 {
                val[e] = -coef[f];
                d += coef[f];
                e++;
            } else if side_kind[4 * c + s] == OPEN {
                d += coef[f];
                open = true;
            }
        }
        val[at] = d;
        rhs[c] = sum;
    }
    if !open { val[row_start[0] + neighbour_count(0)] *= 1.001; }
}

fn neighbour_count(c: i64) -> i64 {
    var n = 0;
    for s in [2, 0] { if neighbour(c, s) >= 0 { n++; } }
    return n;
}

fn pressure_solve() -> bool {
    for c in cell_count() { dp[c] = 0.0; }
    let target = CG_RTOL * solvers::norm(rhs);
    guard target > 0.0 else { return true; };
    guard solvers::amg_setup(row_start, col, val, 1) else { return false; };
    for c in cell_count() { cg_r[c] = rhs[c]; }
    solvers::amg_apply(cg_r, cg_z);
    for c in cell_count() { cg_p[c] = cg_z[c]; }
    var rz = solvers::dot(cg_r, cg_z);
    for it in CG_MAX {
        solvers::multiply(row_start, col, val, cg_p, cg_q);
        let step = rz / solvers::dot(cg_p, cg_q);
        for c in cell_count() {
            dp[c] += step * cg_p[c];
            cg_r[c] -= step * cg_q[c];
        }
        linear_iterations++;
        if solvers::norm(cg_r) <= target { return true; }
        solvers::amg_apply(cg_r, cg_z);
        let next = solvers::dot(cg_r, cg_z);
        for c in cell_count() { cg_p[c] = cg_z[c] + (next / rz) * cg_p[c]; }
        rz = next;
    }
    return false;
}

fn correct() {
    for f in face_count() {
        let fc = face_cells(f);
        if fc[1] >= 0 && fc[2] >= 0 {
            uf[f] = ustar[f] - (coef[f] / face_area(f)) * (dp[fc[2]] - dp[fc[1]]);
            continue;
        }
        let e, s = boundary_cell(f);
        if side_kind[4 * e + s] == OPEN { uf[f] = ustar[f] + outward(s) * (coef[f] / face_area(f)) * dp[e]; }
    }
    for c in cell_count() {
        pres[c] += dp[c];
        for i in 2 {
            let lo = face_of(c, 2 * i);
            let hi = face_of(c, 2 * i + 1);
            vel[2 * c + i] = vstar[2 * c + i] + 0.5 * (uf[lo] - ui[lo] + uf[hi] - ui[hi]);
        }
    }
}

fn boundary_fluxes() {
    for c in cell_count() {
        for s in 4 {
            if side_kind[4 * c + s] >= 0 { uf[face_of(c, s)] = side_normal_velocity(c, s); }
        }
    }
}
