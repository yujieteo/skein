// The explicit step `central_difference` of linear elasticity (spec §14.2):
// M a = f - K u with the lumped mass M, from rest at u = 0 on the free dofs
// and the prescribed values on the fixed dofs. Each step is the central
// difference in the form of velocity Verlet, with one product K u:
//
//     v += dt/2 a;  u += dt v;  a = M⁻¹ (f - K u);  v += dt/2 a
//
// M is the consistent mass of eigen.goose with each diagonal scaled so
// that each component keeps the total mass (Hinton, Rock and Zienkiewicz
// 1976). All simplices of one order have the same ratio of the diagonal
// to the total, so the global scale is the scale of each cell.
//
// The stable step is 2/ω_max. Gershgorin's bound on M⁻¹ K gives
// ω_max² <= max_d Σ_e |K_de| / m_d, so 0.9 times 2/ω of that bound is
// stable. The page of `skein build --target wasm` runs the same steps on
// the GPU (src/gpu/): cd_prepare sets up both.
import std;
import math;
import os;
namespace elements;
import .mesh;
import .assemble;
import .solve;
import .eigen;
import .stages;
import .contact;
import solvers from "../solvers";

let CD_SAFETY = 0.9;    // the fraction of Gershgorin's stable step that `dt: stable` takes

var minv: f64[>..] = [];        // 1/m of each dof, 0 on the fixed dofs
var velocity: f64[>..] = [];
var accel: f64[>..] = [];
var force: f64[>..] = [];       // f of each free dof
var cd_dt = 0.0;
var cd_steps = 0;

// Assembles K, f and the lumped mass, puts the state at rest, and sets the
// step: dt, or the stable step for dt <= 0, shortened so that a whole
// number of steps ends at `end`. False if a free dof has no mass.
fn cd_prepare(end: f64, dt: f64) -> bool {
    guard physics == ELASTIC && dofs_sized() && density > 0.0 && end > 0.0 else { abort("cd_prepare: an elastic physics, a density > 0 and end > 0"); }
    assemble("simplex central_difference", 4);
    hold_dead_nodes();
    assemble_second(MASS);
    let n = sol.len;
    solvers::zeros(minv, n);
    solvers::zeros(velocity, n);
    solvers::zeros(accel, n);
    solvers::zeros(force, n);
    var total: f64[3] = [0.0; 3];
    var diagonal: f64[3] = [0.0; 3];
    for d in n {
        for e in row_start[d]..row_start[d + 1] { total[d % ncomp] += val_b[e]; }
        diagonal[d % ncomp] += val_b[entry(d, d)];
    }
    // the rows of the contact terms
    var tie: f64[>..] = [];
    solvers::zeros(tie, n);
    contact_rows(tie);
    var rate = 0.0;
    for d in n {
        sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 };
        if fixed[d] != 0 { continue; }
        let m = val_b[entry(d, d)] * total[d % ncomp] / diagonal[d % ncomp];
        guard m > 0.0 else { return false; };
        minv[d] = 1.0 / m;
        force[d] = load[d] + assembled_source[d];
        var s = tie[d];
        for e in row_start[d]..row_start[d + 1] { s += abs(val[e]); }
        rate = max(rate, s * minv[d]);
    }
    let h = if dt > 0.0 { dt } else if rate > 0.0 { CD_SAFETY * 2.0 / sqrt(rate) } else { end };
    cd_steps = max(1, ceil(end / h - 1.0e-9) as! i64);
    cd_dt = end / cd_steps;
    return true;
}

// K u and the contact forces at u into out, with the state of the
// contact updated.
fn cd_force(u: f64[:], out: f64[:]) {
    multiply(u, out);
    contact_terms(u, out, false);
    commit_contact();
}

// a = M⁻¹ (f - K u), from the product K u in accel.
fn cd_accel() {
    for d in sol.len { accel[d] = minv[d] * (force[d] - accel[d]); }
}

// The steps of cd_prepare on the CPU. False if the state is not finite.
fn cd_run() -> bool {
    let n = sol.len;
    let h = cd_dt;
    solvers::rows_start(row_start, col, val);
    cd_force(sol, accel);
    cd_accel();
    for i in cd_steps {
        for d in n {
            velocity[d] += 0.5 * h * accel[d];
            sol[d] += h * velocity[d];
        }
        cd_force(sol, accel);
        cd_accel();
        for d in n { velocity[d] += 0.5 * h * accel[d]; }
    }
    solvers::rows_stop();
    var finite = true;
    for v in sol { if !(abs(v) < 1.0e300) { finite = false; } }
    return finite;
}

// The step to `end` with dt (dt <= 0 for the stable step) on the CPU.
fn central_difference_cpu(end: f64, dt: f64) -> bool {
    let start = clock();
    var done = cd_prepare(end, dt);
    let sizes = state_sizes();
    done = done && cd_run();
    guard state_sizes() == sizes else { abort("central_difference: a solve grows no array of the state (rule 3)"); }
    last_iterations = cd_steps;
    last_seconds = clock() - start;
    return done;
}
