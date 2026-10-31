// The dofs, the loads on face sets, the sparse matrix (CSR in three
// arrays), the assembly over the cells of the mesh, the linear solver
// interface with its first entry (CG with a Jacobi preconditioner), the
// Newton interface, the reactions and the values at a point.
//
// The matrix holds all dofs, with the fixed dofs too, so that the
// reactions are K u - f on the fixed dofs. The linear operator of a solve
// is K on the free dofs: it holds its result at 0 on the fixed dofs.
import std;
import math;
import os;
namespace runtime;
import .fem;
import .solid;
import .hex8i;
import .lane_solvers;
import solvers from "../solvers";

var fixed: u8[>..] = [];            // 1 for each dof that is held at 0
var load: f64[>..] = [];            // the external force on each dof
var disp: f64[>..] = [];            // the value of each dof: the displacement of a solid

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];

var last_iterations = 0;
var last_seconds = 0.0;

let CG_TOLERANCE = 1.0e-10;   // the residual relative to the load

// The linear solver and the preconditioner of the next solve, by name.
// "cg" with "jacobi" is the built-in entry. The solvers lane gives the
// others (lane_solvers.goose).
var solver_name: u8[..64] = "cg";
var preconditioner_name: u8[..64] = "jacobi";

fn dof_total() -> i64 { return dofs_per_node * node_count(); }

// Sets all dofs free, with no load and no solution, and clears the matrix.
fn clear_dofs() {
    fixed.clear();
    load.clear();
    disp.clear();
    for d in dof_total() {
        fixed.push(0);
        load.push(0.0);
        disp.push(0.0);
    }
    clear_matrix();
}

fn clear_matrix() {
    row_start.clear();
    col.clear();
    val.clear();
}

fn clear_solution() {
    for &v in disp { v = 0.0; }
}

fn clear_load_vector() {
    for &v in load { v = 0.0; }
    clear_solution();
}

// The sum of the loads on the dofs 0 to 2 of each node.
fn load_sum() -> f64, f64, f64 {
    guard load.len == dof_total() else { abort("load_sum: the loads of the dofs"); }
    var s: f64[3] = [0.0, 0.0, 0.0];
    for n in node_count() {
        for k in min(dofs_per_node, 3) { s[k] += load[dofs_per_node * n + k]; }
    }
    return s[0], s[1], s[2];
}

// Holds all dofs of the nodes of node set s at 0.
fn fix_set(s: i64, first: i64, last: i64) {
    guard s >= 0 && s < set_is_face.len && !set_is_face[s] && first >= 0 && first < last && last <= dofs_per_node && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fix_set: a node set, 0 <= first < last <= dofs_per_node and the dofs"); }
    for n in set_members(s) {
        for p in first..last { fixed[dofs_per_node * n + p] = 1; }
    }
    clear_solution();
}

// --- loads -------------------------------------------------------------------------

var face_w: f64[>..] = [];

// Adds a uniform traction (force per area) (tx, ty, tz) + sn n on the
// faces of face set s, with n the unit normal out of the body: the
// consistent nodal forces of cell_face_weights. False, with a message, if
// a cell type has no face routine.
fn add_set_traction(s: i64, tx: f64, ty: f64, tz: f64, sn: f64) -> bool {
    guard s >= 0 && s < set_is_face.len && set_is_face[s] && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("add_set_traction: a face set and the dofs"); }
    let t: f64[3] = [tx, ty, tz];
    for item in set_members(s) {
        let c = item / FACE_STRIDE;
        if !cell_face_weights(c, item % FACE_STRIDE, face_w) {
            write_stderr(str("runtime: the cell type ", CELL_NAMES[cell_type[c]], " has no face routine\n"));
            return false;
        }
        for a in face_w.len / 4 {
            let d = dofs_per_node * cell_node(c, a);
            for k in min(dofs_per_node, 3) { load[d + k] += face_w[4 * a] * t[k] + face_w[4 * a + 1 + k] * sn; }
        }
    }
    clear_solution();
    return true;
}

// The area of the faces of face set s that have a face routine.
fn set_area(s: i64) -> f64 {
    guard s >= 0 && s < set_is_face.len && set_is_face[s] else { abort("set_area: a face set"); }
    var area = 0.0;
    for item in set_members(s) {
        if cell_face_weights(item / FACE_STRIDE, item % FACE_STRIDE, face_w) {
            for a in face_w.len / 4 { area += face_w[4 * a]; }
        }
    }
    return area;
}

// Adds the force (fx, fy, fz) at the node nearest to (x, y, z).
fn add_point_force(x: f64, y: f64, z: f64, fx: f64, fy: f64, fz: f64) {
    guard node_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("add_point_force: a mesh and its dofs"); }
    var best = 0;
    var best_d = -1.0;
    for n in node_count() {
        let dx = coords[3 * n] - x;
        let dy = coords[3 * n + 1] - y;
        let dz = coords[3 * n + 2] - z;
        let d = dx * dx + dy * dy + dz * dz;
        if best_d < 0.0 || d < best_d { best = n; best_d = d; }
    }
    let f: f64[3] = [fx, fy, fz];
    for k in min(dofs_per_node, 3) { load[dofs_per_node * best + k] += f[k]; }
    clear_solution();
}

// --- the sparsity pattern ------------------------------------------------------

// Two nodes are coupled if they share a cell. All dofs of a node have the
// same columns, in increasing order.
fn build_pattern() {
    guard dofs_per_node > 0 && cell_start.len == cell_count() + 1 && cell_nodes.len == cell_start[cell_count()] else { abort("build_pattern: the cell arrays agree"); }
    clear_matrix();
    let nodes = node_count();
    // The cells of each node, in CSR.
    var first: i64[>..] = [];
    for i in nodes + 1 { first.push(0); }
    for n in cell_nodes { first[n + 1]++; }
    for i in nodes { first[i + 1] += first[i]; }
    var fill: i64[>..] = [];
    for i in nodes { fill.push(first[i]); }
    var cells_of: i64[>..] = [];
    for i in cell_nodes.len { cells_of.push(0); }
    for c in cell_count() {
        for e in cell_start[c]..cell_start[c + 1] {
            let n = cell_nodes[e];
            cells_of[fill[n]] = c;
            fill[n]++;
        }
    }
    let dpn = dofs_per_node;
    var near: i64[>..] = [];
    row_start.push(0);
    for i in nodes {
        near.clear();
        for c in cells_of[first[i]..first[i + 1]] {
            for e in cell_start[c]..cell_start[c + 1] { near.push(cell_nodes[e]); }
        }
        sort(near);
        for p in dpn {
            var last = -1;
            for m in near {
                if m == last { continue; }
                last = m;
                for q in dpn {
                    col.push((dpn * m + q) as i32);
                    val.push(0.0);
                }
            }
            row_start.push(col.len);
        }
    }
}

// The position of column c in row r.
fn entry(r: i64, c: i64) -> i64 {
    guard r >= 0 && c >= 0 && r + 1 < row_start.len && c + 1 < row_start.len else { abort("entry: a row and a column of the matrix"); }
    var lo = row_start[r];
    var hi = row_start[r + 1] - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if col[mid] < c { lo = mid + 1; } else { hi = mid; }
    }
    assert(col[lo] == c);
    return lo;
}

var ke: f64[>..] = [];

// Adds the element matrix of each cell into the matrix, on the pool of the
// threads of the case (src/solvers/pool.goose) if it has more than one
// thread. False, with a message, if a cell type has no element routine in
// the registry.
fn assemble() -> bool {
    guard cell_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("assemble: a mesh and its dofs"); }
    build_pattern();
    let workers = min(solvers::case_threads(), (cell_count() + solvers::CELL_BLOCK - 1) / solvers::CELL_BLOCK);
    var bad = -1;
    if workers < 2 {
        for c in cell_count() {
            guard cell_matrix(c, ke) else { bad = c; break; };
            add_cell(c, ke);
        }
    } else {
        var largest = 0;
        for c in cell_count() { largest = max(largest, cell_start[c + 1] - cell_start[c]); }
        var ids: i64[>..] = [];
        for w in workers { ids.push(thread_spawn(cell_worker)); }
        bad = solvers::gather_cells(cell_count(), workers, (dofs_per_node * largest) * (dofs_per_node * largest), add_cells);
        for id in ids { thread_wait(id); }
    }
    if bad >= 0 {
        write_stderr(str("runtime: the cell type ", CELL_NAMES[cell_type[bad]], " has no element routine\n"));
        return false;
    }
    return true;
}

// A worker of the assembly: the element matrices of the cells of each
// task, one after the other. Its copy of the mesh is the copy of the
// globals at its spawn.
thread_fn cell_worker() {
    var data: f64[>..] = [];
    loop {
        let job = qget<solvers::CellJob>();
        guard job.lo >= 0 else { break; };
        data.clear();
        var bad = -1;
        for c in job.lo..job.hi {
            if bad < 0 && !cell_matrix(c, ke) { bad = c; }
            data.append(ke);
        }
        qput(solvers::CellData { lo: job.lo, bad: bad, data: data[..] });
    }
}

// Adds the element matrices of the cells lo <= c < hi, one after the
// other in m, into the matrix.
fn add_cells(lo: i64, hi: i64, m: f64[:]) {
    var at = 0;
    for c in lo..hi {
        let size = dofs_per_node * (cell_start[c + 1] - cell_start[c]);
        add_cell(c, m[at..at + size * size]);
        at += size * size;
    }
}

// Adds the element matrix k of cell c into the matrix.
fn add_cell(c: i64, k: f64[:]) {
    let dpn = dofs_per_node;
    let nn = cell_start[c + 1] - cell_start[c];
    let size = dpn * nn;
    guard k.len == size * size else { abort("add_cell: the matrix of the nodes of the cell"); }
    for i in nn {
        let ni = cell_node(c, i);
        for j in nn {
            let nj = cell_node(c, j);
            for p in dpn {
                let e = entry(dpn * ni + p, dpn * nj);
                for q in dpn { val[e + q] += k[size * (dpn * i + p) + dpn * j + q]; }
            }
        }
    }
}

// --- the linear solver interface -------------------------------------------------

// A linear operator is a function op(x, y) that sets y = A x. A
// preconditioner is a function pc(r, z) that sets z = M⁻¹ r. A Krylov
// solver takes both as static function values, with this signature:
//
//     fn name<OP, PC>(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, i64
//
// and a call passes the functions after the arguments, as in
// `cg(b, x, tol, max_it, free_operator, jacobi)`.
// It solves A x = b from the start value in x, until |b - A x| <= tol |b|,
// and gives whether it converged and the number of iterations.

// y = K x, on the workers of solvers::rows_start during a solve.
fn multiply(x: f64[:], y: f64[:]) {
    guard x.len == y.len && row_start.len == y.len + 1 && col.len == val.len else { abort("multiply: x and y of the size of the matrix"); }
    solvers::rows_multiply(row_start, col, val, x, y);
}

// y = K x on the free dofs, and 0 on the fixed dofs.
fn free_operator(x: f64[:], y: f64[:]) {
    guard fixed.len == y.len else { abort("free_operator: y of the size of the dofs"); }
    multiply(x, y);
    for d in y.len { if fixed[d] != 0 { y[d] = 0.0; } }
}

var inv_diag: f64[>..] = [];

// z = D⁻¹ r on the free dofs, with D the diagonal of K, and 0 on the
// fixed dofs. solve_free sets inv_diag.
fn jacobi(r: f64[:], z: f64[:]) {
    guard inv_diag.len == r.len && z.len == r.len else { abort("jacobi: r, z and inv_diag of one size"); }
    for d in r.len { z[d] = inv_diag[d] * r[d]; }
}

fn dot(a: f64[:], b: f64[:]) -> f64 {
    guard a.len == b.len else { abort("dot: a and b of one size"); }
    var s = 0.0;
    for v, i in a { s += v * b[i]; }
    return s;
}

fn zeros(n: i64) -> f64[>..] {
    guard n >= 0 else { abort("zeros: n >= 0"); }
    var v: f64[>..] = [];
    for i in n { v.push(0.0); }
    return v;
}

// The preconditioned conjugate gradient method, for a symmetric positive
// definite operator.
fn cg<OP, PC>(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, i64 {
    guard x.len == b.len && max_it >= 0 else { abort("cg: x of the size of b and max_it >= 0"); }
    let n = b.len;
    var r = zeros(n);
    var z = zeros(n);
    var p = zeros(n);
    var q = zeros(n);
    OP(x, q);
    for d in n { r[d] = b[d] - q[d]; }
    let target = tol * sqrt(dot(b, b));
    PC(r, z);
    for d in n { p[d] = z[d]; }
    var rz = dot(r, z);
    var converged = sqrt(dot(r, r)) <= target;
    var it = 0;
    while !converged && it < max_it {
        OP(p, q);
        let alpha = rz / dot(p, q);
        for d in n {
            x[d] += alpha * p[d];
            r[d] -= alpha * q[d];
        }
        it++;
        if sqrt(dot(r, r)) <= target { converged = true; break; }
        PC(r, z);
        let rz_new = dot(r, z);
        let beta = rz_new / rz;
        rz = rz_new;
        for d in n { p[d] = z[d] + beta * p[d]; }
    }
    return converged, it;
}

// Solves K x = b on the free dofs with the solver and the preconditioner
// of solver_name and preconditioner_name. b must be 0 on the fixed dofs.
// It sets last_iterations. It grows no array of the state (state_sizes).
fn solve_free(b: f64[:], x: f64[>..]&) -> bool {
    guard b.len == x.len && b.len == fixed.len && row_start.len == b.len + 1 else { abort("solve_free: b and x of the size of the matrix"); }
    inv_diag.clear();
    for d in b.len { inv_diag.push(if fixed[d] == 0 { 1.0 / val[entry(d, d)] } else { 0.0 }); }
    let sizes = state_sizes();
    let max_it = 10 * b.len;
    var converged = false;
    var found = true;
    if solver_name == "cg" && preconditioner_name == "jacobi" {
        solvers::rows_start(row_start, col, val);
        converged, last_iterations = cg(b, x, CG_TOLERANCE, max_it, free_operator, jacobi);
        solvers::rows_stop();
    } else {
        found, converged, last_iterations = solvers_solve(solver_name, preconditioner_name, b, x, CG_TOLERANCE, max_it);
    }
    guard state_sizes() == sizes else { abort("solve_free: a solve grows no array of the state (rule 3)"); }
    if !found {
        write_stderr(str("runtime: no linear solver ", solver_name, " with the preconditioner ",
                         preconditioner_name, "\n"));
    }
    return found && converged;
}

// The length of each array of the state of the runtime. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[30] {
    return [coords.len, cell_type.len, cell_start.len, cell_nodes.len, point_fields.len, set_names.len,
            set_is_face.len, set_start.len, set_items.len, fixed.len, load.len, disp.len, row_start.len,
            col.len, val.len, face_w.len, ke.len, inv_diag.len, point_n.len, point_dn.len, quad_xi.len,
            quad_w.len, shape_n.len, shape_dn.len, grad.len, centre_n.len, centre_dn.len, full.len,
            condensed.len, mode_grad.len];
}

// Assembles the stiffness matrix and solves K u = f.
fn solve_linear_static() -> bool {
    guard young > 0.0 && cell_count() > 0 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("solve_linear_static: a material, a mesh and its dofs"); }
    let start = clock();
    var converged = assemble();
    if converged {
        var b = zeros(disp.len);
        for d in disp.len {
            disp[d] = 0.0;
            if fixed[d] == 0 { b[d] = load[d]; }
        }
        converged = solve_free(b, disp);
    }
    last_seconds = clock() - start;
    return converged;
}

// --- the Newton interface ----------------------------------------------------------

// A nonlinear problem gives two functions as static function values:
//
//     res(u: f64[:], r: f64[:])   sets r = R(u), the residual, on all dofs
//     jac(u: f64[:])              assembles K = ∂R/∂u into the matrix
//
// and a call passes them after the arguments: `newton(u, tol, max_it, res, jac)`.
// newton solves R(u) = 0 on the free dofs from the start value in u, until
// |R(u)| <= tol |R(u0)|, with the linear solver of solve_free. It gives
// whether it converged and the number of Newton iterations.
fn newton<RES, JAC>(u: f64[>..]&, tol: f64, max_it: i64) -> bool, i64 {
    guard u.len == fixed.len && max_it >= 0 else { abort("newton: u of the size of the dofs and max_it >= 0"); }
    let n = u.len;
    var r = zeros(n);
    var du = zeros(n);
    var it = 0;
    var r0 = -1.0;
    while it < max_it {
        RES(u, r);
        for d in n { if fixed[d] != 0 { r[d] = 0.0; } }
        let size = sqrt(dot(r, r));
        if r0 < 0.0 { r0 = size; }
        if size <= tol * r0 { return true, it; }
        JAC(u);
        for d in n {
            r[d] = -r[d];
            du[d] = 0.0;
        }
        if !solve_free(r, du) { return false, it; }
        for d in n { u[d] += du[d]; }
        it++;
    }
    return false, it;
}

// --- reactions -------------------------------------------------------------------

// The force of the supports on the body: K u - f on the fixed dofs of the
// nodes of node set s. With equilibrium, the reactions and the loads sum
// to 0.
fn set_reaction(s: i64) -> f64, f64, f64 {
    guard s >= 0 && s < set_is_face.len && !set_is_face[s] && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("set_reaction: a node set and the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return 0.0, 0.0, 0.0; }
    for node in set_members(s) { add_node_reaction(node, sum); }
    return sum[0], sum[1], sum[2];
}

// The reaction on all fixed nodes. A node in two fixed sets counts once.
fn fixed_reaction() -> f64, f64, f64 {
    guard fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("fixed_reaction: the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return 0.0, 0.0, 0.0; }
    for node in node_count() { add_node_reaction(node, sum); }
    return sum[0], sum[1], sum[2];
}

// Adds K u - f on the fixed dofs 0 to 2 of `node` to sum.
fn add_node_reaction(node: i64, sum: f64[3]&) {
    guard node >= 0 && node < node_count() && row_start.len == fixed.len + 1 && fixed.len == dof_total() && load.len == fixed.len && disp.len == fixed.len else { abort("add_node_reaction: a node, the dofs and the matrix"); }
    for c in min(dofs_per_node, 3) {
        let d = dofs_per_node * node + c;
        if fixed[d] == 0 { continue; }
        var ku = 0.0;
        for e in row_start[d]..row_start[d + 1] { ku += val[e] * disp[col[e]]; }
        sum[c] += ku - load[d];
    }
}

// --- points ----------------------------------------------------------------------

var point_n: f64[>..] = [];
var point_dn: f64[>..] = [];
var point_xi: f64[3] = [0.0, 0.0, 0.0];

// Whether p is in the box that holds the nodes of cell c.
fn in_cell_bounds(c: i64, p: f64[3]) -> bool {
    guard c >= 0 && c < cell_count() else { abort("in_cell_bounds: a cell"); }
    let first = cell_node(c, 0);
    var lo: f64[3] = [coords[3 * first], coords[3 * first + 1], coords[3 * first + 2]];
    var hi: f64[3] = lo;
    for e in cell_start[c]..cell_start[c + 1] {
        let n = cell_nodes[e];
        for s in 3 {
            lo[s] = min(lo[s], coords[3 * n + s]);
            hi[s] = max(hi[s], coords[3 * n + s]);
        }
    }
    let tol = 1.0e-9 * max(hi[0] - lo[0], max(hi[1] - lo[1], hi[2] - lo[2]));
    for s in 3 {
        if p[s] < lo[s] - tol || p[s] > hi[s] + tol { return false; }
    }
    return true;
}

// The reference point xi of p in the 3D cell c, by Newton's method from
// xi = 0. False if it does not converge.
fn reference_point(c: i64, p: f64[3], xi: f64[3]&) -> bool {
    guard c >= 0 && c < cell_count() else { abort("reference_point: a cell"); }
    let t = cell_type[c] as i64;
    xi = [0.0, 0.0, 0.0];
    for it in 20 {
        guard cell_shape(t, xi, point_n, point_dn) else { return false; }
        var x: f64[3] = [0.0, 0.0, 0.0];
        for w, a in point_n {
            let n = cell_node(c, a);
            for s in 3 { x[s] += w * coords[3 * n + s]; }
        }
        var jac: f64[9] = [0.0; 9];
        cell_jacobian(c, point_dn, jac);
        var inv: f64[9] = [0.0; 9];
        guard inverse3(jac, inv) != 0.0 else { return false; }
        var step = 0.0;
        for r in 3 {
            let d = inv[r] * (p[0] - x[0]) + inv[3 + r] * (p[1] - x[1]) + inv[6 + r] * (p[2] - x[2]);
            xi[r] += d;
            step = max(step, abs(d));
        }
        if step <= 1.0e-12 { return true; }
    }
    return false;
}

// The last 3D cell that holds (x, y, z), with the reference point in
// point_xi and the shape functions there in point_n and point_dn, or -1. On a face between two cells of a lattice
// mesh, the last cell is the cell on the high side of each axis.
fn find_cell(x: f64, y: f64, z: f64) -> i64 {
    guard cell_start.len == cell_count() + 1 else { abort("find_cell: the cell arrays agree"); }
    let p: f64[3] = [x, y, z];
    var xi: f64[3] = [0.0, 0.0, 0.0];
    var c = cell_count() - 1;
    while c >= 0 {
        let t = cell_type[c] as i64;
        if CELL_DIMS[t] == 3 && in_cell_bounds(c, p) && reference_point(c, p, xi) && cell_inside(t, xi) {
            cell_shape(t, xi, point_n, point_dn);
            point_xi = xi;
            return c;
        }
        c--;
    }
    return -1;
}

// The dofs 0 to 2 at (x, y, z), interpolated in the cell of find_cell. A
// point outside the 3D cells gives 0, 0, 0.
fn interpolate_u(x: f64, y: f64, z: f64) -> f64, f64, f64 {
    guard disp.len == dof_total() else { abort("interpolate_u: the solution of the dofs"); }
    guard disp.len > 0 else { return 0.0, 0.0, 0.0; }
    let c = find_cell(x, y, z);
    guard c >= 0 else { return 0.0, 0.0, 0.0; }
    var u: f64[3] = [0.0, 0.0, 0.0];
    for w, a in point_n {
        let d = dofs_per_node * cell_node(c, a);
        for k in min(dofs_per_node, 3) { u[k] += w * disp[d + k]; }
    }
    return u[0], u[1], u[2];
}

// The linear elastic stress at (x, y, z) from the gradient of the dofs 0
// to 2 in the cell of find_cell: xx, yy, zz, yz, xz and xy. A point
// outside the 3D cells gives 0.
fn stress_at_point(x: f64, y: f64, z: f64) -> f64, f64, f64, f64, f64, f64 {
    guard disp.len == dof_total() else { abort("stress_at_point: the solution of the dofs"); }
    guard disp.len > 0 && dofs_per_node >= 3 else { return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; }
    let c = find_cell(x, y, z);
    guard c >= 0 else { return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; }
    var jac: f64[9] = [0.0; 9];
    cell_jacobian(c, point_dn, jac);
    var inv: f64[9] = [0.0; 9];
    let det = inverse3(jac, inv);
    guard det != 0.0 else { return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; }
    // g[3 k + s] = ∂u_k / ∂x_s
    var g: f64[9] = [0.0; 9];
    for a in point_n.len {
        let d = dofs_per_node * cell_node(c, a);
        for s in 3 {
            let ds = inv[3 * s] * point_dn[3 * a] + inv[3 * s + 1] * point_dn[3 * a + 1]
                + inv[3 * s + 2] * point_dn[3 * a + 2];
            for k in 3 { g[3 * k + s] += disp[d + k] * ds; }
        }
    }
    if cell_type[c] as i64 == HEX8I { hex8i_add_modes(c, point_xi, det, g); }
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    let tr = lambda * (g[0] + g[4] + g[8]);
    return tr + 2.0 * mu * g[0], tr + 2.0 * mu * g[4], tr + 2.0 * mu * g[8],
        mu * (g[5] + g[7]), mu * (g[2] + g[6]), mu * (g[1] + g[3]);
}
