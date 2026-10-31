// The physics, the method and the conditions of the cloud.
//
// The physics is one scalar unknown u with
//
//     div(a u) - d div(grad(u)) = f(x, y)
//
// and a, d constant, or linear elasticity in plane stress or plane strain
// with the unknown (u_x, u_y). A condition k fixes a component on the
// selected boundary edges at the value BC(k, x, y, t) of the solve. A
// boundary edge that no condition holds has no flux: (a u - d grad(u))·n
// = 0, or the traction or the outward flux of load_selected (0 if none).
import std;
import math;
namespace meshless;
import .cloud;
import .rk;
import .rkpm;
import .gfd;

let SCALAR = 0;
let ELASTIC = 1;
let RKPM = 0;
let GFD = 1;

var physics = -1;
var ncomp = 1;
var ax = 0.0;
var ay = 0.0;
var diffusivity = 0.0;
var lame = 0.0;          // lambda, or its plane stress value 2 lambda mu / (lambda + 2 mu)
var lame_3d = 0.0;       // lambda of the 3D law, for sigma_zz in plane strain
var shear = 0.0;
var plane_stress = true;
var thickness = 1.0;     // multiplies each integral

var method = -1;
var stencil = 0;         // gfd: the points of each stencil

var sol: f64[>..] = [];  // rkpm: the coefficient of each point; gfd: the value at each point. ncomp each.
var time = 0.0;

// The physics is unset until equation or elastic.
fn start_physics() {
    physics = -1;
    method = -1;
    sol.clear();
}

fn dof_count() -> i64 { return ncomp * point_count(); }

// Each point is a cell of the solve records.
fn cell_count() -> i64 { return point_count(); }

// The scalar equation div(a u) - d div(grad(u)) = f.
fn equation(a_x: f64, a_y: f64, d: f64) {
    guard point_count() > 0 && d > 0.0 else { abort("equation: a cloud and d > 0"); }
    physics = SCALAR;
    ncomp = 1;
    ax = a_x;
    ay = a_y;
    diffusivity = d;
    thickness = 1.0;
    size_state();
}

// Linear elasticity with Young's modulus E and Poisson's ratio nu, in
// plane stress (`stress`) or plane strain, with the thickness t.
fn elastic(E: f64, nu: f64, stress: bool, t: f64) {
    guard point_count() > 0 && E > 0.0 && nu > -1.0 && nu < 0.5 && t > 0.0 else { abort("elastic: a cloud, E > 0, -1 < nu < 0.5 and t > 0"); }
    physics = ELASTIC;
    ncomp = 2;
    shear = E / (2.0 * (1.0 + nu));
    lame_3d = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    lame = if stress { 2.0 * lame_3d * shear / (lame_3d + 2.0 * shear) } else { lame_3d };
    plane_stress = stress;
    thickness = t;
    size_state();
}

fn size_state() {
    method = -1;
    sol.clear();
    for d in dof_count() { sol.push(0.0); }
    edge_fix.clear();
    edge_load.clear();
    for e in ncomp * edge_count() {
        edge_fix.push(FREE);
        edge_load.push(0.0);
    }
    time = 0.0;
}

// rkpm with the linear basis and the support `factor` times the spacing
// of each point. False if a shape function is singular: the support is
// too small for the cloud.
fn rkpm(factor: f64) -> bool {
    guard physics >= 0 && factor > 1.0 else { abort("rkpm: a physics and a support factor > 1"); }
    method = RKPM;
    set_support(factor);
    return rkpm_setup();
}

// gfd with `points` points in each stencil. The scalar physics only.
// False if a stencil cannot fit the derivatives.
fn gfd(points: i64) -> bool {
    guard physics == SCALAR && points >= 6 && points < point_count() else { abort("gfd: the scalar physics, and 6 <= points < the points of the cloud"); }
    method = GFD;
    stencil = points;
    // The buckets of the stencil search: the spacing of each point.
    set_support(1.0);
    return gfd_setup();
}

// Frees each boundary edge, and removes the loads.
fn clear_bcs() {
    guard method >= 0 && edge_fix.len == ncomp * edge_count() else { abort("clear_bcs: a physics and a method"); }
    for &k in edge_fix { k = FREE; }
    for &v in edge_load { v = 0.0; }
}

// Selects the boundary edges whose midpoint (x, y) has P(x, y). False if
// it selects none.
fn select_faces<P>() -> bool {
    guard selected.len == edge_count() else { abort("select_faces: a cloud"); }
    var found = false;
    for e in edge_count() {
        let m = edge_point(e, 0.5);
        let on = P(m[0], m[1]);
        selected[e] = if on { 1 } else { 0 };
        if on { found = true; }
    }
    return found;
}

// Selects the boundary edges with both points in `mask` (1 for each point
// of the region). False if it selects none.
fn select_points(mask: u8[:]) -> bool {
    guard mask.len == point_count() && selected.len == edge_count() else { abort("select_points: a mask of the points of the cloud"); }
    var found = false;
    for e in edge_count() {
        let on = mask[edges[2 * e]] != 0 && mask[edges[2 * e + 1]] != 0;
        selected[e] = if on { 1 } else { 0 };
        if on { found = true; }
    }
    return found;
}

// Puts the condition k on component `comp` of each selected edge.
fn fix_selected(k: i64) { fix_component(0, k); }

fn fix_component(comp: i64, k: i64) {
    guard k >= 0 && comp >= 0 && comp < ncomp && method >= 0 && edge_fix.len == ncomp * edge_count() else { abort("fix_component: k >= 0, a component of the unknown, and a method"); }
    for e in edge_count() { if selected[e] != 0 { edge_fix[ncomp * e + comp] = k; } }
}

// Adds the value v to the load of component `comp` on each selected edge:
// a traction, or the outward flux of the scalar physics.
fn load_selected(comp: i64, v: f64) {
    guard comp >= 0 && comp < ncomp && method >= 0 && edge_load.len == ncomp * edge_count() else { abort("load_selected: a component of the unknown and a method"); }
    for e in edge_count() { if selected[e] != 0 { edge_load[ncomp * e + comp] += v; } }
}
