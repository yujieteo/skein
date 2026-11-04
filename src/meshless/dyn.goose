import std;
import math;
import os;
import solvers from "../solvers";
namespace meshless;
import .cloud;
import .rk;
import .model;
import .rkpm;

let NITSCHE_EXPLICIT = 4.0;

var mass: f64[>..] = [];
var vhalf: f64[>..] = [];
var force: f64[>..] = [];
var rec_p: f64[>..] = [];
var node_start: i64[>..] = [];
var node_item: i64[>..] = [];
var node_psi: f64[>..] = [];
var modulus: f64[>..] = [];
var dyn_steps = 0;
var dyn_seconds = 0.0;
var dyn_courant = 0.0;
var started = false;
var divergence = 0.0;

fn explicit_setup(factor: f64, rho: f64) -> bool {
    guard point_count() > 0 && factor > 1.0 && rho > 0.0 else { abort("explicit_setup: a cloud, a support factor > 1 and rho > 0"); }
    physics = ELASTIC;
    ncomp = 2;
    thickness = 1.0;
    size_state();
    set_support(factor);
    if !build_records() { return false; }
    method = RKPM;
    mass.clear();
    for p in point_count() { mass.push(0.0); }
    for r in record_count() {
        if rec_kind[r] == SMOOTH { mass[tri[r]] += rho * rec_w[r]; }
    }
    vhalf.clear();
    force.clear();
    rec_p.clear();
    for d in dof_count() {
        vhalf.push(0.0);
        force.push(0.0);
    }
    for r in 4 * record_count() { rec_p.push(0.0); }
    node_start.clear();
    node_item.clear();
    node_psi.clear();
    node_start.push(0);
    for p in point_count() {
        if !shape_at(px(p), py(p)) { return false; }
        node_item.append(near);
        node_psi.append(psi);
        node_start.push(node_item.len);
    }
    started = false;
    divergence = 1000.0 * max(hi[0] - lo[0], hi[1] - lo[1]);
    return true;
}

fn start_displacement<U0>() {
    guard sol.len == dof_count() && mass.len == point_count() else { abort("start_displacement: explicit_setup"); }
    for p in point_count() {
        let u = U0(px(p), py(p));
        sol[2 * p] = u[0];
        sol[2 * p + 1] = u[1];
    }
    for &v in vhalf { v = 0.0; }
    started = false;
    time = 0.0;
}

fn node_u(c: i64, p: i64) -> f64 {
    var v = 0.0;
    for i in node_start[p]..node_start[p + 1] { v += node_psi[i] * sol[2 * node_item[i] + c]; }
    return v;
}

fn node_gradient(p: i64) -> f64[4] {
    var f: f64[4] = [1.0, 0.0, 0.0, 1.0];
    let r = first_record_of(p);
    for i in rec_start[r]..rec_start[r + 1] {
        let j = item[i];
        f[0] += item_dx[i] * sol[2 * j];
        f[1] += item_dy[i] * sol[2 * j];
        f[2] += item_dx[i] * sol[2 * j + 1];
        f[3] += item_dy[i] * sol[2 * j + 1];
    }
    return f;
}

fn first_record_of(p: i64) -> i64 {
    let t = point_tris[point_start[p]];
    var k = 0;
    while tri[3 * t + k] != p { k++; }
    return 3 * t + k;
}

fn internal_force<STRESS>() {
    for &v in force { v = 0.0; }
    for r in record_count() {
        if rec_kind[r] != SMOOTH { continue; }
        var f: f64[4] = [1.0, 0.0, 0.0, 1.0];
        for i in rec_start[r]..rec_start[r + 1] {
            let j = item[i];
            f[0] += item_dx[i] * sol[2 * j];
            f[1] += item_dy[i] * sol[2 * j];
            f[2] += item_dx[i] * sol[2 * j + 1];
            f[3] += item_dy[i] * sol[2 * j + 1];
        }
        let s = STRESS(f);
        for k in 4 { rec_p[4 * r + k] = s[k]; }
        let w = rec_w[r];
        for i in rec_start[r]..rec_start[r + 1] {
            let j = item[i];
            force[2 * j] -= w * (s[0] * item_dx[i] + s[1] * item_dy[i]);
            force[2 * j + 1] -= w * (s[2] * item_dx[i] + s[3] * item_dy[i]);
        }
    }
}

fn wall_modulus<STRESS>(n: f64[2], c: i64) -> f64 {
    let h = 1.0e-6;
    var f: f64[4] = [1.0, 0.0, 0.0, 1.0];
    let s0 = STRESS(f);
    f[2 * c] += h * n[0];
    f[2 * c + 1] += h * n[1];
    let s1 = STRESS(f);
    var d = 0.0;
    for k in 2 { d += (s1[2 * c + k] - s0[2 * c + k]) * n[k]; }
    return abs(d) / h;
}

fn wall_force<STRESS, BC>() {
    for r in record_count() {
        if rec_kind[r] != EDGE { continue; }
        let e = rec_edge[r];
        let n, l = edge_normal(e);
        let o = rec_owner[r];
        for c in 2 {
            let k = edge_fix[2 * e + c];
            if k == FREE { continue; }
            let g = BC(k, rec_x[r], rec_y[r], time);
            var u = 0.0;
            for i in rec_start[r]..rec_start[r + 1] { u += item_psi[i] * sol[2 * item[i] + c]; }
            let beta = NITSCHE_EXPLICIT * wall_modulus(n, c, STRESS) / l;
            let t = rec_p[4 * o + 2 * c] * n[0] + rec_p[4 * o + 2 * c + 1] * n[1];
            let q = rec_w[r] * (t - beta * (u - g));
            for i in rec_start[r]..rec_start[r + 1] { force[2 * item[i] + c] += q * item_psi[i]; }
        }
    }
}

fn dyn_state_sizes() -> i64[10] {
    return [sol.len, vhalf.len, force.len, rec_p.len, mass.len, item.len, rec_kind.len, node_item.len, edge_fix.len, rec_g.len];
}

fn advance<STRESS, BC>(dt: f64, bx: f64, by: f64) -> bool {
    internal_force(STRESS);
    wall_force(STRESS, BC);
    var finite = true;
    for p in point_count() {
        let ax = force[2 * p] / mass[p] + bx;
        let ay = force[2 * p + 1] / mass[p] + by;
        let half = if started { 1.0 } else { 0.5 };
        vhalf[2 * p] += half * dt * ax;
        vhalf[2 * p + 1] += half * dt * ay;
        sol[2 * p] += dt * vhalf[2 * p];
        sol[2 * p + 1] += dt * vhalf[2 * p + 1];
        if !(abs(sol[2 * p]) < divergence && abs(sol[2 * p + 1]) < divergence) { finite = false; }
    }
    started = true;
    time += dt;
    return finite;
}

fn solve_dynamic<STRESS, BC, AFTER>(end: f64, dt: f64, bx: f64, by: f64) -> bool {
    guard method == RKPM && mass.len == point_count() && end > time && dt > 0.0 else { abort("solve_dynamic: explicit_setup, an end after the time of the state, and dt > 0"); }
    let start = clock();
    let sizes = dyn_state_sizes();
    solvers::memory_check("meshless explicit", item.len, dof_count(), 4, 0, "");
    let steps = max(1, ceil((end - time) / dt - 1.0e-9) as! i64);
    let k = (end - time) / steps;
    for i in steps {
        if !advance(k, bx, by, STRESS, BC) {
            write_stderr(str("skein run: the displacement is not finite at t = ", time, " s: use a shorter step\n"));
            return false;
        }
        AFTER(i + 1, time);
    }
    dyn_steps = steps;
    guard dyn_state_sizes() == sizes else { abort("solve_dynamic: a solve grows no array of the state (rule 3)"); }
    dyn_seconds = clock() - start;
    return true;
}

fn write_dynamic_vtu(path: const u8[:], show_u: bool, names: const u8[:], values: f64[:]) -> bool {
    guard method == RKPM && sol.len == dof_count() else { abort("write_dynamic_vtu: explicit_setup"); }
    let n = point_count();
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", n, "\" NumberOfCells=\"", triangle_count(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for p in n { format(out, px(p) + node_u(0, p), " ", py(p) + node_u(1, p), " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for t in triangle_count() { format(out, tri[3 * t], " ", tri[3 * t + 1], " ", tri[3 * t + 2], "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for t in triangle_count() { format(out, 3 * (t + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for t in triangle_count() { out.append("5\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    if show_u {
        out.append("<DataArray type=\"Float64\" Name=\"u\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for p in n { format(out, node_u(0, p), " ", node_u(1, p), " 0\n"); }
        out.append("</DataArray>\n");
    }
    var field = 0;
    var from = 0;
    for at in names.len + 1 {
        if at < names.len && names[at] != ',' { continue; }
        if at == from { continue; }
        format(out, "<DataArray type=\"Float64\" Name=\"", names[from..at], "\" format=\"ascii\">\n");
        for p in n { format(out, values[field * n + p], "\n"); }
        out.append("</DataArray>\n");
        field++;
        from = at + 1;
    }
    out.append("</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
