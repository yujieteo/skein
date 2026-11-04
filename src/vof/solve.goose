import std;
import math;
import os;
import solvers from "../solvers";
namespace vof;
import .mesh;
import .advect;
import .step;
import .project;

let STEP_COURANT = 1;
let STEP_FINITE = 2;
let STEP_LINEAR = 3;

fn state_sizes() -> i64[24] {
    return [side_kind.len, selected.len, side_vx.len, side_vy.len, alpha.len, vel.len, pres.len, uf.len, af.len, cour.len,
            alpha_new.len, alpha_s.len, normal.len, kappa.len, rho_c.len, mu_c.len, rho_adv.len, vstar.len, gradv.len,
            ustar.len, ui.len, coef.len, row_start.len + col.len + val.len, rhs.len + dp.len + cg_r.len + cg_z.len + cg_p.len + cg_q.len];
}

fn advance<PROPS, VALUE>(dt: f64) -> i64 {
    set_boundary_values(time + dt, VALUE);
    boundary_fluxes();
    let co = courant(dt);
    max_courant = max(max_courant, co);
    if !(co <= 1.0) { return STEP_COURANT; }
    advect(dt);
    properties(PROPS);
    if tension > 0.0 { curvature(); }
    predict(dt);
    face_fluxes(dt);
    pressure_matrix();
    if !pressure_solve() { return STEP_LINEAR; }
    correct();
    time += dt;
    var finite = true;
    for v in vel { if !(abs(v) < 1.0e30) { finite = false; } }
    return if finite { 0 } else { STEP_FINITE };
}

fn step_fault(code: i64, dt: f64) {
    if code == STEP_COURANT {
        write_stderr(str("skein run: the Courant number of the flow reaches ", max_courant, " at t = ", time, " s, above 1: use a step dt below ",
                         dt / max_courant, " s (the step of ", dt, " s is too long for the velocity and the cells)\n"));
    } else if code == STEP_LINEAR {
        write_stderr(str("skein run: the pressure solve does not converge at t = ", time, " s after ", CG_MAX, " iterations: use a shorter step or a coarser mesh\n"));
    } else {
        write_stderr(str("skein run: the velocity is not finite at t = ", time, " s: use a shorter step\n"));
    }
}

fn solve_transient<PROPS, VALUE, AFTER>(end: f64, dt: f64) -> bool {
    guard alpha.len == cell_count() && rhs.len == cell_count() && end > time && dt > 0.0 else { abort("solve_transient: a mesh, a start, end after the time of the state, and dt > 0"); }
    let start = clock();
    let sizes = state_sizes();
    solvers::memory_check("vof cg amg", col.len, cell_count(), 12, 0, "");
    let steps = max(1, ceil((end - time) / dt - 1.0e-9) as! i64);
    let k = (end - time) / steps;
    linear_iterations = 0;
    max_courant = 0.0;
    for i in steps {
        let code = advance(k, PROPS, VALUE);
        guard code == 0 else {
            step_fault(code, k);
            return false;
        };
        AFTER(i + 1, time);
    }
    steps_done = steps;
    guard state_sizes() == sizes else { abort("solve_transient: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - start;
    return true;
}
