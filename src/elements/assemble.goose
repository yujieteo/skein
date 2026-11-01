// The element routines of linear elasticity and heat conduction on the
// simplex mesh, the assembly into a sparse matrix (CSR in three arrays),
// the boundary conditions and the reactions.
//
// Each node has `ncomp` dofs: dim for elasticity, 1 for heat. Dof d is
// component d % ncomp of node d / ncomp. The matrix holds all dofs, with the
// fixed dofs too, so that the reactions are K u - f on the fixed dofs.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .kernels;   // cell_matrix
import .csg;
import .stages;
import solvers from "../solvers";

let ELASTIC = 0;
let HEAT = 1;
let WAVE = 2;      // the complex scalar wave of wave.goose
let EDGE = 3;      // the edge elements of edge.goose
let MIXED = 4;     // the mixed Darcy flow of mixed.goose

// The plane assumption of an elastic 2D model.
let SOLID_3D = 0;
let PLANE_STRESS = 1;
let PLANE_STRAIN = 2;

var physics = -1;
var ncomp = 0;
var plane = SOLID_3D;
var thickness = 1.0;    // multiplies the integrals of a 2D model
var lame = 0.0;         // lambda, or its plane stress value 2 lambda mu / (lambda + 2 mu)
var lame_3d = 0.0;      // lambda of the 3D law, for sigma_zz in plane strain
var shear = 0.0;        // mu
var conductivity = 0.0;
var source = 0.0;       // the heat source per volume
var volume_force: f64[3] = [0.0; 3];   // the elastic body force per volume

var fixed: u8[>..] = [];
var prescribed: f64[>..] = [];
var load: f64[>..] = [];
var sol: f64[>..] = [];

var row_start: i64[>..] = [];
var col: i32[>..] = [];
var val: f64[>..] = [];

var other_dofs = 0;     // the dofs of an edge or a mixed physics, which are not on the nodes

fn dof_total() -> i64 { return if physics == EDGE || physics == MIXED { other_dofs } else { ncomp * node_count() } }

// Whether each array of the dofs has one entry for each dof, and the
// matrix, if there is one, one row for each dof.
fn dofs_sized() -> bool {
    return fixed.len == dof_total() && prescribed.len == fixed.len && load.len == fixed.len && sol.len == fixed.len &&
           (row_start.len == 0 || (row_start.len == fixed.len + 1 && assembled_source.len == fixed.len));
}

fn clear_matrix() {
    row_start.clear();
    col.clear();
    val.clear();
}

fn start_physics(kind: i64, comps: i64) {
    guard per_cell > 0 && (kind == ELASTIC || kind == HEAT || kind == WAVE) && comps >= 1 && comps <= 3 else { abort("start_physics: a mesh, a physics and 1 to 3 components"); }
    physics = kind;
    ncomp = comps;
    clear_stages();
    fixed.clear();
    prescribed.clear();
    load.clear();
    sol.clear();
    for d in dof_total() {
        fixed.push(0);
        prescribed.push(0.0);
        load.push(0.0);
        sol.push(0.0);
    }
    clear_matrix();
}

fn set_elastic(E: f64, nu: f64, assumption: i64, t: f64) {
    guard E > 0.0 && nu > -1.0 && nu < 0.5 && t > 0.0 && assumption >= SOLID_3D && assumption <= PLANE_STRAIN else { abort("set_elastic: E > 0, -1 < nu < 0.5, t > 0 and a plane assumption"); }
    start_physics(ELASTIC, dim);
    plane = if dim == 3 { SOLID_3D } else { assumption };
    thickness = if dim == 3 { 1.0 } else { t };
    shear = E / (2.0 * (1.0 + nu));
    lame_3d = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    lame = if plane == PLANE_STRESS { E * nu / (1.0 - nu * nu) } else { lame_3d };
    volume_force = [0.0; 3];
}

fn set_heat(k: f64, q: f64, t: f64) {
    guard k > 0.0 && t > 0.0 else { abort("set_heat: k > 0 and t > 0"); }
    start_physics(HEAT, 1);
    thickness = if dim == 3 { 1.0 } else { t };
    conductivity = k;
    source = q;
}

// --- the geometry of a cell at a point ---------------------------------------

var shp: f64[10];
var dref: f64[30];
var grad: f64[30];      // dN_l / dx_t at grad[3 l + t]

fn cell_node(c: i64, l: i64) -> i64 {
    guard c >= 0 && c < cell_count() && l >= 0 && l < per_cell else { abort("cell_node: c < cell_count() and l < per_cell"); }
    return conn[per_cell * c + l];
}

// Fills shp and grad at reference point xi of cell c, and gives the
// determinant of the Jacobian.
fn cell_geometry(c: i64, xi: f64[3]) -> f64 {
    guard c >= 0 && c < cell_count() && per_cell <= shp.len && 3 * per_cell <= grad.len else { abort("cell_geometry: c < cell_count() and shp and grad hold the nodes of a cell"); }
    shape(dim, order, xi, shp, dref);
    // jac[3 r + t] = dx_t / dxi_r
    var jac: f64[9] = [0.0; 9];
    for l in per_cell {
        let n = cell_node(c, l);
        for r in dim {
            for t in dim { jac[3 * r + t] += dref[3 * l + r] * node_x(n, t); }
        }
    }
    var inv: f64[9] = [0.0; 9];
    var det = 0.0;
    if dim == 2 {
        det = jac[0] * jac[4] - jac[1] * jac[3];
        inv[0] = jac[4] / det;
        inv[1] = -jac[1] / det;
        inv[3] = -jac[3] / det;
        inv[4] = jac[0] / det;
    } else {
        det = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7])
            - jac[1] * (jac[3] * jac[8] - jac[5] * jac[6])
            + jac[2] * (jac[3] * jac[7] - jac[4] * jac[6]);
        inv[0] = (jac[4] * jac[8] - jac[5] * jac[7]) / det;
        inv[1] = (jac[2] * jac[7] - jac[1] * jac[8]) / det;
        inv[2] = (jac[1] * jac[5] - jac[2] * jac[4]) / det;
        inv[3] = (jac[5] * jac[6] - jac[3] * jac[8]) / det;
        inv[4] = (jac[0] * jac[8] - jac[2] * jac[6]) / det;
        inv[5] = (jac[2] * jac[3] - jac[0] * jac[5]) / det;
        inv[6] = (jac[3] * jac[7] - jac[4] * jac[6]) / det;
        inv[7] = (jac[1] * jac[6] - jac[0] * jac[7]) / det;
        inv[8] = (jac[0] * jac[4] - jac[1] * jac[3]) / det;
    }
    assert(det > 0.0);
    // inv = jac⁻¹, so dN/dx_t = Σr inv[3 t + r] dN/dxi_r
    for l in per_cell {
        for t in 3 {
            var g = 0.0;
            for r in dim { g += inv[3 * t + r] * dref[3 * l + r]; }
            grad[3 * l + t] = g;
        }
    }
    return det;
}

// --- the cell matrix ---------------------------------------------------------

// cell_matrix(c) of kernels.goose fills ke and fe for cell c.
var ke: f64[900];       // (per_cell ncomp)², row ncomp l + p for component p of node l
var fe: f64[30];

// --- the sparsity pattern and the assembly -----------------------------------

// Two nodes are neighbours if they share a cell. The ncomp rows of a node
// have the same columns, in increasing order.
fn build_pattern() {
    guard per_cell > 0 && ncomp > 0 else { abort("build_pattern: a mesh and a physics"); }
    clear_matrix();
    let nn = node_count();
    var keys: i64[>..] = [];
    for c in cell_count() {
        for a in per_cell {
            for b in per_cell { keys.push(cell_node(c, a) * nn + cell_node(c, b)); }
        }
    }
    sort(keys);
    row_start.push(0);
    var first = 0;
    for a in nn {
        var last = first;
        while last < keys.len && keys[last] / nn == a { last++; }
        for p in ncomp {
            var prev = -1;
            for e in first..last {
                if keys[e] == prev { continue; }
                prev = keys[e];
                let b = keys[e] % nn;
                for q in ncomp {
                    col.push((ncomp * b + q) as i32);
                    val.push(0.0);
                }
            }
            row_start.push(col.len);
        }
        first = last;
    }
}

// The position of column c in row r.
fn entry(r: i64, c: i64) -> i64 {
    guard r >= 0 && r + 1 < row_start.len && c >= 0 && c + 1 < row_start.len else { abort("entry: r and c are rows of the matrix"); }
    var lo = row_start[r];
    var hi = row_start[r + 1] - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if col[mid] < c { lo = mid + 1; } else { hi = mid; }
    }
    assert(col[lo] == c);
    return lo;
}

// The matrix, and the source vector added to `load`.
var assembled_source: f64[>..] = [];

// The matrix and the source vector, on the pool of the threads of the
// case (src/solvers/pool.goose) if it has more than one thread. `path` and
// `vectors` are the solution path and the vectors of its solver, for the
// memory check after the pattern (solvers::memory_check).
fn assemble(path: u8[:], vectors: i64) {
    guard physics >= 0 && dofs_sized() else { abort("assemble: a physics and the arrays of the dofs"); }
    build_pattern();
    internal.clear();
    solvers::memory_check(path, col.len, dof_total(), vectors, worker_mesh_bytes(), "");
    assembled_source.clear();
    for d in dof_total() { assembled_source.push(0.0); }
    let m = per_cell * ncomp;
    let workers = min(solvers::case_threads(), (cell_count() + solvers::CELL_BLOCK - 1) / solvers::CELL_BLOCK);
    if workers < 2 {
        for c in cell_count() {
            if !is_live(c) { continue; }
            cell_terms(c);
            add_cell(c, ke[..m * m], fe[..m]);
        }
        use_law_0();
        return;
    }
    var ids: i64[>..] = [];
    for w in workers { ids.push(thread_spawn(cell_worker)); }
    solvers::gather_cells(cell_count(), workers, m * m + m, add_cells);
    for id in ids { thread_wait(id); }
}

// The bytes of the globals that each cell_worker copies at its spawn.
fn worker_mesh_bytes() -> i64 {
    return 8 * (coords.len + conn.len + cell_law.len + birth_u.len) + dead.len + born.len;
}

// The matrix ke and the source fe of cell c, with the law of its region.
// A born cell adds K_e u_b on its nodes (stages.goose).
fn cell_terms(c: i64) {
    let m = per_cell * ncomp;
    use_law(c);
    cell_matrix(c);
    if born.len > 0 && born[c] != 0 {
        for r in m {
            var s = 0.0;
            for k in m { s += ke[m * r + k] * birth_u[m * c + k]; }
            fe[r] += s;
        }
    }
}

// A worker of the assembly: the matrix and then the source of each cell of
// each task, one cell after the other. Its copy of the mesh is the copy of
// the globals at its spawn. A dead cell sends the values of the cell before
// it, and add_cells skips it.
thread_fn cell_worker() {
    let m = per_cell * ncomp;
    var data: f64[>..] = [];
    loop {
        let job = qget<solvers::CellJob>();
        guard job.lo >= 0 else { break; };
        data.clear();
        for c in job.lo..job.hi {
            if is_live(c) { cell_terms(c); }
            data.append(ke[..m * m]);
            data.append(fe[..m]);
        }
        qput(solvers::CellData { lo: job.lo, bad: -1, data: data[..] });
    }
}

// Adds the matrix and the source of the live cells lo <= c < hi, one cell
// after the other in d.
fn add_cells(lo: i64, hi: i64, d: f64[:]) {
    let m = per_cell * ncomp;
    for c in lo..hi {
        if !is_live(c) { continue; }
        let at = (c - lo) * (m * m + m);
        add_cell(c, d[at..at + m * m], d[at + m * m..at + m * m + m]);
    }
}

// Adds the matrix k and the source f of cell c.
fn add_cell(c: i64, k: f64[:], f: f64[:]) {
    let m = per_cell * ncomp;
    guard k.len == m * m && f.len == m else { abort("add_cell: the matrix and the source of the nodes of a cell"); }
    for i in per_cell {
        let ni = cell_node(c, i);
        for p in ncomp { assembled_source[ncomp * ni + p] += f[ncomp * i + p]; }
        for j in per_cell {
            let nj = cell_node(c, j);
            for p in ncomp {
                let e = entry(ncomp * ni + p, ncomp * nj);
                for q in ncomp { val[e + q] += k[m * (ncomp * i + p) + ncomp * j + q]; }
            }
        }
    }
}

// --- boundary regions --------------------------------------------------------

// A region of the boundary: a union of parts. A part is the plane where
// coordinate `axis` is `value`, the surface of the hole at distance
// `radius` from the z axis through (cx, cy), or the surface of the shape
// `shape` of src/mesher/. On the mesh of a shape (csg.goose), a plane
// holds only the nodes on the surface of that shape.
struct Part {
    hole: bool,
    axis: i64,
    value: f64,
    c: f64[2],
    r: f64,
    shape: i64,
}

var region: Part[..8] = [];

fn on_part(n: i64, k: i64) -> bool {
    guard n >= 0 && n < node_count() && k >= 0 && k < region.len else { abort("on_part: n < node_count() and k < region.len"); }
    let tol = geometry_tolerance();
    let p = region[k];
    if p.shape >= 0 { return shape_gap(p.shape, n) <= tol; }
    if !p.hole { return abs(node_x(n, p.axis) - p.value) <= tol && (body_shape < 0 || shape_gap(body_shape, n) <= tol); }
    let dx = node_x(n, 0) - p.c[0];
    let dy = node_x(n, 1) - p.c[1];
    return abs(sqrt(dx * dx + dy * dy) - p.r) <= tol;
}

fn on_region(n: i64) -> bool {
    guard n >= 0 && n < node_count() else { abort("on_region: n < node_count()"); }
    for k in region.len { if on_part(n, k) { return true; } }
    return false;
}

fn region_node_count() -> i64 {
    var count = 0;
    for n in node_count() { if on_region(n) { count++; } }
    return count;
}

// Holds component `comp` (or all components if comp < 0) of each node of
// the region at `v`.
fn fix_region(comp: i64, v: f64) {
    guard physics >= 0 && comp >= -1 && comp < ncomp && dofs_sized() else { abort("fix_region: a physics, -1 <= comp < ncomp and the arrays of the dofs"); }
    for n in node_count() {
        if !on_region(n) { continue; }
        for p in ncomp {
            if comp >= 0 && p != comp { continue; }
            fixed[ncomp * n + p] = 1;
            prescribed[ncomp * n + p] = v;
        }
    }
}

// Holds component `comp` of each node of the region at VALUE(x, y, z).
fn fix_region_with<VALUE>(comp: i64) {
    guard physics >= 0 && comp >= 0 && comp < ncomp && dofs_sized() else { abort("fix_region_with: a physics, 0 <= comp < ncomp and the arrays of the dofs"); }
    for n in node_count() {
        if !on_region(n) { continue; }
        fixed[ncomp * n + comp] = 1;
        prescribed[ncomp * n + comp] = VALUE(node_x(n, 0), node_x(n, 1), node_x(n, 2));
    }
}

// Adds ∫ g Ni dA over the facets of the region to component p of node i,
// for g = gv[p].
fn add_region_load(gv: f64[3]) {
    guard physics >= 0 && ncomp <= gv.len && dofs_sized() else { abort("add_region_load: a physics and the arrays of the dofs"); }
    let fd = dim - 1;
    let fnodes = node_count_of(fd, order);
    var fn_shape: f64[10] = [0.0; 10];
    var fn_dref: f64[30] = [0.0; 30];
    for c in cell_count() {
        for f in facet_count(dim) {
            // A facet with its vertices on two parts can cut a corner of
            // the body, so all vertices must be on one part.
            var on = false;
            for k in region.len {
                var all = true;
                for l in fd + 1 { if !on_part(cell_node(c, facet_node(dim, f, l)), k) { all = false; } }
                if all { on = true; }
            }
            if !on { continue; }
            for qp in quad_count(fd) {
                let xi, w = quad_point(fd, qp);
                shape(fd, order, xi, fn_shape, fn_dref);
                // The tangents dx/dxi_r of the facet.
                var tan: f64[6] = [0.0; 6];
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for r in fd {
                        for t in 3 { tan[3 * r + t] += fn_dref[3 * l + r] * node_x(n, t); }
                    }
                }
                var da = 0.0;
                if fd == 1 {
                    da = sqrt(tan[0] * tan[0] + tan[1] * tan[1]) * thickness;
                } else {
                    let nx = tan[1] * tan[5] - tan[2] * tan[4];
                    let ny = tan[2] * tan[3] - tan[0] * tan[5];
                    let nz = tan[0] * tan[4] - tan[1] * tan[3];
                    da = sqrt(nx * nx + ny * ny + nz * nz);
                }
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for p in ncomp { load[ncomp * n + p] += gv[p] * fn_shape[l] * w * da; }
                }
            }
        }
    }
}

fn clear_conditions() {
    held.clear();
    held.append(sol);
    for &v in fixed { v = 0; }
    for &v in prescribed { v = 0.0; }
    for &v in load { v = 0.0; }
    for &v in sol { v = 0.0; }
}

// --- reactions -----------------------------------------------------------------

// The internal force of each dof after a static step of the nonlinear
// solid (solid.goose), in place of K u; empty after a linear solve.
var internal: f64[>..] = [];

// The force K u of dof d, or its internal force, and the sum of the sizes
// of the terms of the sum K u (the size of its round-off).
fn row_force(d: i64) -> f64, f64 {
    if internal.len == sol.len { return internal[d], abs(internal[d]); }
    var ku = 0.0;
    var terms = 0.0;
    for e in row_start[d]..row_start[d + 1] {
        ku += val[e] * sol[col[e]];
        terms += abs(val[e] * sol[col[e]]);
    }
    return ku, terms;
}

// The force of the supports on the body (or the heat that flows in through
// them): K u - f on the fixed dofs of the region. With equilibrium, the
// reactions and the applied loads sum to 0.
fn region_reaction() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("region_reaction: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return sum; }
    for n in node_count() {
        if !on_region(n) { continue; }
        for p in ncomp {
            let d = ncomp * n + p;
            if fixed[d] != 1 { continue; }
            let ku, terms = row_force(d);
            sum[p] += ku - load[d] - assembled_source[d];
        }
    }
    return sum;
}

// The sum of the reactions on all fixed dofs.
fn all_reaction() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("all_reaction: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    guard row_start.len > 0 else { return sum; }
    for d in sol.len {
        if fixed[d] != 1 { continue; }
        let ku, terms = row_force(d);
        sum[d % ncomp] += ku - load[d] - assembled_source[d];
    }
    return sum;
}

// The scale of the check `equilibrium`: the larger of the sum of the sizes
// of the nodal loads, sources and reactions, and the sum of the sizes of
// the terms K_ij u_j on the fixed rows. Loads that cancel do not make the
// first 0. The reactions keep it above 0 when only fixed values drive the
// model. The second is the size of the round-off of K u, which is not 0
// when one nonzero value at all fixed dofs drives nothing.
fn load_magnitude() -> f64 {
    guard dofs_sized() else { abort("load_magnitude: the arrays of the dofs"); }
    var sum = 0.0;
    for v in load { sum += abs(v); }
    for v in assembled_source { sum += abs(v); }
    guard row_start.len > 0 else { return sum; }
    var terms = 0.0;
    for d in sol.len {
        if fixed[d] != 1 { continue; }
        let ku, t = row_force(d);
        terms += t;
        sum += abs(ku - load[d] - assembled_source[d]);
    }
    return max(sum, terms);
}

// The sum of the applied loads and sources, by component.
fn total_load() -> f64[3] {
    guard physics >= 0 && ncomp <= 3 && dofs_sized() else { abort("total_load: a physics and the arrays of the dofs"); }
    var sum: f64[3] = [0.0, 0.0, 0.0];
    for v, d in load { sum[d % ncomp] += v; }
    for v, d in assembled_source { sum[d % ncomp] += v; }
    return sum;
}
