// The nonlinear solid of spec §10.6 on the simplex mesh: a law with state
// at each integration point (material.goose), in the total Lagrangian
// form on the reference mesh. A cell gives the internal force
//
//     f_int(a i) = ∫ F_iJ S_JK ∂N_a/∂X_K dV
//
// and, for a static step, the tangent ∫ Bᵀ D B + ∇N_a · S ∇N_b δ_ij dV,
// with F = I and no second term at small strain.
//
// A static step takes `ramp` increments. Increment k moves the fixed
// values from those of the end of the step before to those of the step,
// and the loads likewise, by k / ramp. Each increment is Newton's method
// with the skyline of direct.goose on the tangent, or CG and AMG
// (solve.goose) when the skyline is too large: the first solve also
// moves the fixed dofs, so the tangent predicts the change of the free
// dofs. After the first increment, an increment starts from u plus the
// change of u in the increment before, and its solves need less accuracy.
// The state of each point changes only at the end of an increment.
//
// The explicit step is central_difference (explicit.goose) with the
// internal force in place of K u, from the velocity of the unknown v of
// the first-order form (spec §14.1). Each step updates the state.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .solve;
import .stages;
import .material;
import .explicit;
import .direct;
import solvers from "../solvers";

let SOLID_NONE = 0;
let SOLID_J2 = 1;           // J2 at small strain
let SOLID_J2_FINITE = 2;    // J2 at finite strain, the multiplicative split
let SOLID_ENERGY = 3;       // an energy of the program

let NEWTON_MAX = 30;        // the Newton iterations of an increment
let NEWTON_TOL = 1.0e-9;    // the residual of the free dofs, relative to the forces

var solid_law = SOLID_NONE;
var j2 = J2 { lame: 0.0, shear: 0.0, yield0: 0.0, hard: 0.0 };
var pt_state: f64[>..] = [];    // 6 per point: εp or Cp⁻¹ (Voigt), at the end of the last increment or step
var pt_alpha: f64[>..] = [];
var trial_state: f64[>..] = []; // the same at the iterate of Newton's method
var trial_alpha: f64[>..] = [];
var pt_stress: f64[>..] = [];   // 6 per point: σ, τ or S (Voigt) of the last evaluation
var nl_u: f64[>..] = [];        // the displacement of the iterate
var nl_start: f64[>..] = [];    // u at the start of the step
var nl_target: f64[>..] = [];   // the fixed values of the step
var nl_load: f64[>..] = [];     // the loads of the step
var nl_last_load: f64[>..] = [];    // the loads of the end of the step before
var nl_base: f64[>..] = [];     // u at the start of the increment
var nl_step: f64[>..] = [];     // the change of u in the increment before
var newton_steps = 0;           // the Newton iterations of the last static step
var newton_direct = false;      // the step solves with the skyline
var pt_s: f64[6] = [0.0; 6];    // the stress and the tangent of the last point
var pt_d: f64[36] = [0.0; 36];

fn point_total() -> i64 { return cell_count() * quad_count(dim); }

// The law of the solid. J2 needs the elastic law first (set_elastic).
fn set_solid(law: i64, yield0: f64, hard: f64) {
    guard physics == ELASTIC && law > SOLID_NONE && law <= SOLID_ENERGY && plane != PLANE_STRESS && cell_law.len == 0 &&
          (law == SOLID_ENERGY || yield0 > 0.0 && hard >= 0.0) else { abort("set_solid: an elastic physics, not plane stress, one law, and a yield stress > 0 and a hardening >= 0 for J2"); }
    solid_law = law;
    j2 = J2 { lame: lame_3d, shear: shear, yield0: yield0, hard: hard };
    pt_state.clear();
    pt_alpha.clear();
    trial_state.clear();
    trial_alpha.clear();
    pt_stress.clear();
    let start: f64[6] = if law == SOLID_J2_FINITE { [1.0, 1.0, 1.0, 0.0, 0.0, 0.0] } else { [0.0; 6] };
    for q in point_total() {
        pt_state.append(start);
        trial_state.append(start);
        pt_alpha.push(0.0);
        trial_alpha.push(0.0);
        pt_stress.append([0.0; 6]);
    }
    nl_u.clear();
    nl_start.clear();
    nl_target.clear();
    nl_load.clear();
    nl_last_load.clear();
    nl_base.clear();
    nl_step.clear();
    for d in dof_total() {
        nl_base.push(0.0);
        nl_step.push(0.0);
        nl_u.push(0.0);
        nl_start.push(0.0);
        nl_target.push(0.0);
        nl_load.push(0.0);
        nl_last_load.push(0.0);
    }
}

// The deformation gradient F = I + Σ u_l ⊗ ∇N_l of cell c, from the
// gradients of the last cell_geometry.
fn cell_defgrad(c: i64, u: f64[:]) -> f64[9] {
    var f = EYE;
    for l in per_cell {
        let n = cell_node(c, l);
        for i in ncomp {
            for j in dim { f[3 * i + j] += u[ncomp * n + i] * grad[3 * l + j]; }
        }
    }
    return f;
}

// The stress and the tangent of point q at F into pt_s and pt_d, from the
// state of the end of the last increment, which the trial state gets.
fn point_response<ENERGY>(q: i64, f: f64[9]) {
    var st: f64[6] = [0.0; 6];
    for k in 6 { st[k] = pt_state[6 * q + k]; }
    var alpha = pt_alpha[q];
    if solid_law == SOLID_J2 {
        var eps = (f + mat_t(f)) * 0.5;
        for i in 3 { eps[4 * i] -= 1.0; }
        j2_small(j2, eps, st, alpha, pt_s, pt_d);
    } else {
        let c, e = principal_c(f);
        var w1: f64[3] = [0.0; 3];
        var w2: f64[9] = [0.0; 9];
        energy_derivatives(c, dim, w1, w2, ENERGY);
        energy_stress(c, e, dim, w1, w2, pt_s, pt_d);
    }
    for k in 6 {
        trial_state[6 * q + k] = st[k];
        pt_stress[6 * q + k] = pt_s[k];
    }
    trial_alpha[q] = alpha;
}

// The principal values and directions of C = Fᵀ F. In 2d the direction
// out of the plane comes last, so the first two are those of the plane.
fn principal_c(f: f64[9]) -> f64[3], f64[9] {
    var c, e = sym_eigen(mat_mul(mat_t(f), f));
    guard dim == 2 else { return c, e; };
    var k = 2;
    for j in 2 { if abs(e[6 + j]) > abs(e[6 + k]) { k = j; } }
    guard k != 2 else { return c, e; };
    let ck = c[k];
    c[k] = c[2];
    c[2] = ck;
    for i in 3 {
        let x = e[3 * i + k];
        e[3 * i + k] = e[3 * i + 2];
        e[3 * i + 2] = x;
    }
    return c, e;
}

// Adds ∫ P ∇N_a dV of the first Piola-Kirchhoff stress P to fe.
fn add_point_force(p: f64[9], dv: f64) {
    for l in per_cell {
        for i in ncomp {
            var s = 0.0;
            for k in dim { s += p[3 * i + k] * grad[3 * l + k]; }
            fe[ncomp * l + i] += s * dv;
        }
    }
}

// Adds ∫ Bᵀ D B dV, and with `finite` the term of the stress, at the last
// point to ke.
fn add_point_tangent(f: f64[9], dv: f64, finite: bool) {
    let m = per_cell * ncomp;
    var b: f64[180] = [0.0; 180];  // b[m v + r] for the Voigt row v and the dof r
    var db: f64[180] = [0.0; 180];
    for l in per_cell {
        for i in ncomp {
            let r = ncomp * l + i;
            for v in 6 {
                let x = VOIGT_I[v];
                let y = VOIGT_J[v];
                b[m * v + r] = f[3 * i + x] * grad[3 * l + y] + (if x == y { 0.0 } else { f[3 * i + y] * grad[3 * l + x] });
            }
        }
    }
    for v in 6 {
        for r in m {
            var s = 0.0;
            for w in 6 { s += pt_d[6 * v + w] * b[m * w + r]; }
            db[m * v + r] = s;
        }
    }
    for r in m {
        for k in m {
            var s = 0.0;
            for v in 6 { s += b[m * v + r] * db[m * v + k]; }
            ke[m * r + k] += s * dv;
        }
    }
    guard finite else { return; };
    let sm = from_voigt(pt_s);
    for a in per_cell {
        for c in per_cell {
            var g = 0.0;
            for x in dim {
                for y in dim { g += grad[3 * a + x] * sm[3 * x + y] * grad[3 * c + y]; }
            }
            for i in ncomp { ke[m * (ncomp * a + i) + ncomp * c + i] += g * dv; }
        }
    }
}

// The internal force fe and, with `tangent`, the tangent ke of cell c at
// the displacement u.
fn solid_cell<ENERGY>(c: i64, u: f64[:], tangent: bool) {
    let m = per_cell * ncomp;
    for k in m * m { ke[k] = 0.0; }
    for k in m { fe[k] = 0.0; }
    let finite = solid_law != SOLID_J2;
    for qp in quad_count(dim) {
        let xi, w = quad_point(dim, qp);
        let dv = w * cell_geometry(c, xi) * thickness;
        let f = cell_defgrad(c, u);
        point_response(c * quad_count(dim) + qp, f, ENERGY);
        let fk = if finite { f } else { EYE };
        add_point_force(mat_mul(fk, from_voigt(pt_s)), dv);
        if tangent { add_point_tangent(fk, dv, finite); }
    }
}

// The tangent in val and the internal force in `internal`, at nl_u.
fn solid_assemble<ENERGY>() {
    let m = per_cell * ncomp;
    for &v in val { v = 0.0; }
    for &v in assembled_source { v = 0.0; }
    for c in cell_count() {
        solid_cell(c, nl_u, true, ENERGY);
        add_cell(c, ke[..m * m], fe[..m]);
    }
    for d in internal.len {
        internal[d] = assembled_source[d];
        assembled_source[d] = 0.0;
    }
}

// The energy of a law with no energy.
fn no_energy(c: Hd[3]) -> Hd { return hd(0.0); }

// A static step of the nonlinear solid in `ramp` increments, with the
// energy ENERGY of the program for SOLID_ENERGY. False if an increment
// does not converge.
fn solve_solid<ENERGY>(ramp: i64) -> bool {
    guard (solid_law == SOLID_J2 || solid_law == SOLID_ENERGY) && ramp >= 1 && dofs_sized() && nl_u.len == sol.len else { abort("solve_solid: set_solid of J2 at small strain or of an energy, ramp >= 1 and the arrays of the dofs"); }
    let start = clock();
    build_pattern();
    newton_direct = direct_setup();
    if !newton_direct { solvers::memory_check("simplex newton cg amg", col.len, dof_total(), 6, worker_mesh_bytes(), ""); }
    assembled_source.clear();
    internal.clear();
    for d in dof_total() {
        assembled_source.push(0.0);
        internal.push(0.0);
    }
    for d in sol.len {
        nl_start[d] = if held.len == sol.len { held[d] } else { 0.0 };
        nl_u[d] = nl_start[d];
        nl_target[d] = prescribed[d];
        nl_load[d] = load[d];
        nl_step[d] = 0.0;
    }
    let sizes = state_sizes();
    var its = 0;
    newton_steps = 0;
    var done = true;
    for k in 1..ramp + 1 {
        if !done { break; }
        let lam = k * 1.0 / ramp;
        for d in sol.len {
            nl_base[d] = nl_u[d];
            if k > 1 { nl_u[d] = if fixed[d] != 0 { nl_start[d] + lam * (nl_target[d] - nl_start[d]) } else { nl_u[d] + nl_step[d] }; }
        }
        done = increment(lam, k > 1, start, its, ENERGY);
        if done { commit_points(); }
        for d in sol.len { nl_step[d] = nl_u[d] - nl_base[d]; }
    }
    guard state_sizes() == sizes else { abort("solve_solid: a solve grows no array of the state (rule 3)"); }
    for d in sol.len {
        sol[d] = nl_u[d];
        load[d] = nl_load[d];
        nl_last_load[d] = nl_load[d];
    }
    last_iterations = its;
    last_seconds = clock() - start;
    return done;
}

// The state of the end of the increment.
fn commit_points() {
    for k in pt_state.len { pt_state[k] = trial_state[k]; }
    for q in pt_alpha.len { pt_alpha[q] = trial_alpha[q]; }
}

// Newton's method for the increment to the fraction lam of the step. its
// gets the CG iterations. A solve takes the residual to 1e-4 of its start,
// but the first solve of an increment with no prediction takes it to
// CG_TOLERANCE: it moves the fixed dofs alone, and a coarse solve of that
// can turn cells inside out.
fn increment<ENERGY>(lam: f64, predicted: bool, start: f64, its: i64&) -> bool {
    let n = sol.len;
    for it in NEWTON_MAX {
        solid_assemble(ENERGY);
        var r2 = 0.0;
        var scale = 0.0;
        var moved = 0.0;
        for d in n {
            let f = nl_last_load[d] + lam * (nl_load[d] - nl_last_load[d]);
            scale += f * f + internal[d] * internal[d];
            let target = nl_start[d] + lam * (nl_target[d] - nl_start[d]);
            sol[d] = if fixed[d] != 0 { target - nl_u[d] } else { 0.0 };
            load[d] = if fixed[d] != 0 { 0.0 } else { f - internal[d] };
            moved = max(moved, abs(sol[d]));
            if fixed[d] == 0 { r2 += load[d] * load[d]; }
        }
        if it > 0 && moved == 0.0 && sqrt(r2) <= NEWTON_TOL * sqrt(scale) { return true; }
        cg_tolerance = if it == 0 && !predicted { CG_TOLERANCE } else { 1.0e-2 };
        let solved = if newton_direct { direct_solve(start) } else { cg_solve(start, 0, false) };
        cg_tolerance = CG_TOLERANCE;
        guard solved else { return false; };
        its += last_iterations;
        newton_steps++;
        for d in n { nl_u[d] += sol[d]; }
        if fixed.len > 0 {
            for d in n { if fixed[d] != 0 { nl_u[d] = nl_start[d] + lam * (nl_target[d] - nl_start[d]); } }
        }
    }
    return false;
}

// --- the explicit step --------------------------------------------------------

// The internal force of J2 at finite strain at the displacement u into
// out, with the state of each point updated.
fn explicit_force(u: f64[:], out: f64[:]) {
    for &v in out { v = 0.0; }
    let m = per_cell * ncomp;
    for c in cell_count() {
        for k in m { fe[k] = 0.0; }
        for qp in quad_count(dim) {
            let xi, w = quad_point(dim, qp);
            let dv = w * cell_geometry(c, xi) * thickness;
            let q = c * quad_count(dim) + qp;
            var cpi: f64[6] = [0.0; 6];
            for k in 6 { cpi[k] = pt_state[6 * q + k]; }
            var tau: f64[6] = [0.0; 6];
            var alpha = pt_alpha[q];
            let p = j2_finite(j2, cell_defgrad(c, u), cpi, alpha, tau);
            for k in 6 {
                pt_state[6 * q + k] = cpi[k];
                pt_stress[6 * q + k] = tau[k];
            }
            pt_alpha[q] = alpha;
            add_point_force(p, dv);
        }
        for l in per_cell {
            for i in ncomp { out[ncomp * cell_node(c, l) + i] += fe[ncomp * l + i]; }
        }
    }
}

// The start of an explicit step of J2 at finite strain to the time `end`
// with the step dt (dt <= 0 for the stable step of the elastic law): the
// velocity v0 on the free dofs, and the acceleration. The GPU path
// (src/gpu/) starts with it too.
fn solid_explicit_start(end: f64, dt: f64, v0: f64[3]) -> bool {
    guard solid_law == SOLID_J2_FINITE else { abort("solid_explicit: set_solid of J2 at finite strain"); }
    guard cd_prepare(end, dt) else { return false; };
    for d in sol.len { velocity[d] = if fixed[d] != 0 { 0.0 } else { v0[d % ncomp] }; }
    explicit_force(sol, accel);
    cd_accel();
    return true;
}

// The explicit step of J2 at finite strain from solid_explicit_start.
// False if the state is not finite.
fn solid_explicit(end: f64, dt: f64, v0: f64[3]) -> bool {
    let start = clock();
    guard solid_explicit_start(end, dt, v0) else { return false; };
    let n = sol.len;
    let sizes = state_sizes();
    let h = cd_dt;
    var steps = 0;
    var finite = true;
    while steps < cd_steps && finite {
        for d in n {
            velocity[d] += 0.5 * h * accel[d];
            sol[d] += h * velocity[d];
        }
        explicit_force(sol, accel);
        cd_accel();
        for d in n { velocity[d] += 0.5 * h * accel[d]; }
        steps++;
        // Reason: a step above the stable step gives NaN, and the steps of
        // NaN are slow (the Taylor bar at 0.8 mm took 13 min to 40 us).
        if steps % 100 == 0 || steps == cd_steps { finite = all_finite(sol); }
    }
    guard state_sizes() == sizes else { abort("solid_explicit: a solve grows no array of the state (rule 3)"); }
    last_iterations = steps;
    last_seconds = clock() - start;
    return finite;
}

fn all_finite(v: f64[:]) -> bool {
    for x in v { if !(abs(x) < 1.0e300) { return false; } }
    return true;
}
