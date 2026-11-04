// The complex scalar wave physics of a `harmonic` step (spec §14.2):
// div(grad(u)) + c u = s, for a `Complex` unknown u at one frequency, with
// the structural damping of the step, a PML (spec §13) and the ties of
// `bloch` groups (spec §12.6).
//
// A complex value a + ib is two dofs of the node: component 0 is a and
// component 1 is b. A complex matrix entry a + ib is the 2 × 2 block
// [a, -b; b, a] of the real matrix, so the real system is the complex one.
//
// The weak form is ∫ (1 + 2iζ) Λ∇u·∇v - c det(S) u v = -∫ det(S) s v +
// ∫ g v dA, with g the flux dot(grad(u), n) of a natural condition
// (add_region_load). In the PML, S = diag(s_j) stretches each coordinate
// and Λ = det(S) S⁻². Outside it, S = I.
//
// A tied node of a `bloch` group takes the value of its master times the
// phase exp(iθ). The assembly adds each term of a tied node to its master,
// times the phase (the test function too, with the conjugate phase), so
// the row of the tied node is the identity. A direct sparse LU
// (src/solvers/slu.goose) solves the real system, then the tied nodes take their values.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .eigen;
import .solve;
import solvers from "../solvers";

let WAVE_TOLERANCE = 1.0e-10;  // the GMRES residual relative to the right side
let PML_REFLECTION = 1.0e-6;    // the reflection of the layer at normal incidence

var wave_c: f64[2] = [0.0, 0.0];   // c, the coefficient of u
var wave_s: f64[2] = [0.0, 0.0];   // s, the source
var wave_damping = 0.0;            // ζ
var wave_cell_terms: f64[>..] = [];  // for each cell: c (re, im), s (re, im) and the diffusion ν
var cell_mark: i64[>..] = [];
var pml_shape = -1;                // the shape of the layer, or -1
var pml_flag: u8[>..] = [];
var pml_lo: f64[3] = [0.0; 3];     // the rest of the body: the box from pml_lo to pml_hi
var pml_hi: f64[3] = [0.0; 3];
var pml_h = 0.0;                   // the thickness of the layer
var tie: i64[>..] = [];            // the master of each node, or -1
var tie_phase: f64[>..] = [];      // cos θ and sin θ of each node
var low_nodes: i64[>..] = [];      // the nodes of the low region of the next bloch
var robin_cell: i64[>..] = [];     // the facets of the conditions dot(grad(u), n) = z u: cell, facet and z
var robin_facet: i64[>..] = [];
var robin_re: f64[>..] = [];
var robin_im: f64[>..] = [];

fn cmul(a: f64[2], b: f64[2]) -> f64[2] { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }

fn cdiv(a: f64[2], b: f64[2]) -> f64[2] {
    let d = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
}

fn conj(a: f64[2]) -> f64[2] { return [a[0], -a[1]]; }

// Adds the condition dot(grad(u), n) = z u, with z complex and constant, on
// the facets of the region (the exact outgoing wave of a sphere, an
// impedance). It adds -∫ z u v dA to the left side.
fn add_robin(z: f64[2]) {
    guard physics == WAVE && dofs_sized() else { abort("add_robin: a wave physics"); }
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            robin_cell.push(c);
            robin_facet.push(f);
            robin_re.push(z[0]);
            robin_im.push(z[1]);
        }
    }
}

fn clear_robin() {
    robin_cell.clear();
    robin_facet.clear();
    robin_re.clear();
    robin_im.clear();
}

// The matrix of the Robin facets, on the masters.
fn robin_assemble() {
    var fn_shape: f64[20] = [0.0; 20];
    let fnodes = node_count_of(dim - 1, order);
    for c, i in robin_cell {
        let f = robin_facet[i];
        for qp in quad_count(dim - 1) {
            let wda = facet_point(c, f, qp, fn_shape);
            for a in fnodes {
                let mi, phi = master(cell_node(c, facet_node(dim, f, a)));
                for b in fnodes {
                    let mj, pj = master(cell_node(c, facet_node(dim, f, b)));
                    let z = cmul(cmul(conj(phi), [-robin_re[i] * fn_shape[a] * fn_shape[b] * wda, -robin_im[i] * fn_shape[a] * fn_shape[b] * wda]), pj);
                    add_block(mi, mj, z);
                }
            }
        }
    }
}

// Starts the wave physics on the mesh. In 2D, t is the thickness.
fn set_wave(t: f64) {
    guard t > 0.0 else { abort("set_wave: t > 0"); }
    start_physics(WAVE, 2);
    thickness = if dim == 3 { 1.0 } else { t };
    wave_c = [0.0, 0.0];
    wave_s = [0.0, 0.0];
    wave_damping = 0.0;
    pml_shape = -1;
    tie.clear();
    tie_phase.clear();
    for n in node_count() {
        tie.push(-1);
        tie_phase.append([1.0, 0.0]);
    }
}

// The terms of one frequency: c, s and ζ.
fn set_wave_terms(c: f64[2], s: f64[2], damping: f64) {
    guard physics == WAVE && damping >= 0.0 else { abort("set_wave_terms: a wave physics and damping >= 0"); }
    wave_c = c;
    wave_s = s;
    wave_damping = damping;
    wave_cell_terms.clear();
    cell_mark.clear();
    for e in cell_count() {
        wave_cell_terms.append([c[0], c[1], s[0], s[1], 1.0]);
        cell_mark.push(-1);
    }
}

// Gives all cells the diffusion nu.
fn wave_material_all(nu: f64) {
    guard physics == WAVE && nu > 0.0 && wave_cell_terms.len == 5 * cell_count() else { abort("wave_material_all: a wave physics with its terms set, and nu > 0"); }
    for e in cell_count() { wave_cell_terms[5 * e + 4] = nu; }
}

// The integral of u over the cells in the shape g: its real and imaginary part, then the area of the cells.
fn wave_cells_integral(g: i64) -> f64[3] {
    guard physics == WAVE && g >= 0 else { abort("wave_cells_integral: a wave physics and a shape"); }
    var total: f64[3] = [0.0; 3];
    for c in cell_count() {
        if !cell_in(c, g) { continue; }
        for q in rich_count() {
            let xi, w = rich_point(q);
            let dv = w * cell_geometry(c, xi) * thickness;
            var u: f64[2] = [0.0, 0.0];
            for l in per_cell {
                u[0] += shp[l] * sol[2 * cell_node(c, l)];
                u[1] += shp[l] * sol[2 * cell_node(c, l) + 1];
            }
            total[0] += u[0] * dv;
            total[1] += u[1] * dv;
            total[2] += dv;
        }
    }
    return total;
}

// Gives the cells in the shape g the diffusion ν, the coefficient c and the
// source s, for div(ν grad(u)) + c u = s. It gives the number of these cells.
fn wave_material(g: i64, nu: f64, c: f64[2], s: f64[2]) -> i64 {
    guard physics == WAVE && g >= 0 && nu > 0.0 && wave_cell_terms.len == 5 * cell_count() else {
        abort("wave_material: a wave physics with its terms set, a shape and nu > 0");
    }
    var count = 0;
    for e in cell_count() {
        if !cell_in(e, g) { continue; }
        wave_cell_terms[5 * e] = c[0];
        wave_cell_terms[5 * e + 1] = c[1];
        wave_cell_terms[5 * e + 2] = s[0];
        wave_cell_terms[5 * e + 3] = s[1];
        wave_cell_terms[5 * e + 4] = nu;
        cell_mark[e] = g;
        count++;
    }
    return count;
}

fn wave_source<S>(mark: i64) {
    guard physics == WAVE && dofs_sized() && cell_mark.len == cell_count() else { abort("wave_source: a wave physics with its terms set"); }
    for c in cell_count() {
        if !source_cell(c, mark) { continue; }
        for q in rich_count() {
            let xi, w = rich_point(q);
            let dv = w * cell_geometry(c, xi) * thickness;
            var b: f64[4] = [1.0, 0.0, 0.0, 0.0];
            for r in dim {
                b[r + 1] = xi[r];
                b[0] -= xi[r];
            }
            let p = source_point(c, b);
            let st = stretch(c, p, pml_k(c));
            var det: f64[2] = [1.0, 0.0];
            for j in dim { det = cmul(det, [st[2 * j], st[2 * j + 1]]); }
            let f = cmul(S(p), det);
            for l in per_cell {
                load[2 * cell_node(c, l)] -= f[0] * shp[l] * dv;
                load[2 * cell_node(c, l) + 1] -= f[1] * shp[l] * dv;
            }
        }
    }
}

// The layer of shape g, past the box from lo to hi, of thickness h.
fn set_pml(g: i64, lo: f64[3], hi: f64[3], h: f64) {
    guard (physics == WAVE || physics == EDGE) && g >= 0 && h > 0.0 else { abort("set_pml: a wave or edge physics, a shape and h > 0"); }
    pml_shape = g;
    pml_flag.clear();
    for c in cell_count() { pml_flag.push(if cell_in(c, g) { 1 } else { 0 }); }
    pml_lo = lo;
    pml_hi = hi;
    pml_h = h;
}

// Keeps the nodes of the region, as the low region of the next tie.
fn keep_low() -> bool {
    low_nodes.clear();
    for n in node_count() { if on_region(n) { low_nodes.push(n); } }
    return low_nodes.len > 0;
}

// Ties each node of the region to the node of the low region at its point
// less `period`, with the phase θ. False if a node has no partner.
fn tie_high(period: f64[3], theta: f64) -> bool {
    guard physics == WAVE && tie.len == node_count() else { abort("tie_high: a wave physics"); }
    let tol = geometry_tolerance();
    var all = true;
    for n in node_count() {
        if !on_region(n) { continue; }
        var best = -1;
        var best_d = tol;
        for m in low_nodes {
            var d = 0.0;
            for t in 3 { d = max(d, abs(node_x(n, t) - period[t] - node_x(m, t))); }
            if d <= best_d && m != n {
                best = m;
                best_d = d;
            }
        }
        if best < 0 {
            all = false;
            continue;
        }
        tie[n] = best;
        tie_phase[2 * n] = cos(theta);
        tie_phase[2 * n + 1] = sin(theta);
    }
    return all;
}

// The master of node n through the ties, and the product of the phases.
fn master(n: i64) -> i64, f64[2] {
    var m = n;
    var p: f64[2] = [1.0, 0.0];
    var steps = 0;
    while tie[m] >= 0 && steps < 4 {
        p = cmul(p, [tie_phase[2 * m], tie_phase[2 * m + 1]]);
        m = tie[m];
        steps++;
    }
    assert(tie[m] < 0);
    return m, p;
}

// --- the assembly --------------------------------------------------------------

// The pattern of the masters: two masters are neighbours if their nodes
// share a cell. The row of a tied node holds only its diagonal.
fn wave_pattern() {
    clear_matrix();
    let nn = node_count();
    var keys: i64[>..] = [];
    for c in cell_count() {
        for a in per_cell {
            let ma, pa = master(cell_node(c, a));
            for b in per_cell {
                let mb, pb = master(cell_node(c, b));
                keys.push(ma * nn + mb);
            }
        }
    }
    for n in nn { if tie[n] >= 0 { keys.push(n * nn + n); } }
    sort(keys);
    row_start.push(0);
    var first = 0;
    for a in nn {
        var last = first;
        while last < keys.len && keys[last] / nn == a { last++; }
        for p in 2 {
            var prev = -1;
            for e in first..last {
                if keys[e] == prev { continue; }
                prev = keys[e];
                let b = keys[e] % nn;
                for q in 2 {
                    col.push((2 * b + q) as i32);
                    val.push(0.0);
                }
            }
            row_start.push(col.len);
        }
        first = last;
    }
}

// Adds the complex value z at row node a and column node b.
fn add_block(a: i64, b: i64, z: f64[2]) {
    let e0 = entry(2 * a, 2 * b);
    let e1 = entry(2 * a + 1, 2 * b);
    val[e0] += z[0];
    val[e0 + 1] -= z[1];
    val[e1] += z[1];
    val[e1 + 1] += z[0];
}

fn pml_k(c: i64) -> f64 {
    if physics == EDGE { return if edge_c[2 * c] < 0.0 { sqrt(-edge_c[2 * c] / edge_nu[c]) } else { 0.0 }; }
    return if wave_cell_terms[5 * c] > 0.0 { sqrt(wave_cell_terms[5 * c] / wave_cell_terms[5 * c + 4]) } else { 0.0 };
}

fn source_cell(c: i64, mark: i64) -> bool { return mark == -2 || cell_mark[c] == mark; }

fn source_point(c: i64, l: f64[4]) -> f64[3] {
    var p: f64[3] = [0.0; 3];
    for v in dim + 1 {
        for t in 3 { p[t] += l[v] * node_x(cell_node(c, v), t); }
    }
    return p;
}

// The stretch s_j of each axis at point p of cell c, in the PML.
fn stretch(c: i64, p: f64[3], k: f64) -> f64[6] {
    var s: f64[6] = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
    guard pml_shape >= 0 && pml_flag[c] != 0 else { return s; }
    guard k > 0.0 else { abort("stretch: a cell of the PML layer has a wave number k > 0"); }
    let top = 3.0 * log(1.0 / PML_REFLECTION) / (2.0 * k * pml_h);
    for j in dim {
        let d = max(0.0, max(pml_lo[j] - p[j], p[j] - pml_hi[j]));
        s[2 * j + 1] = -top * (d / pml_h) * (d / pml_h);
    }
    return s;
}

// The complex matrix me and source fe of cell c, node by node.
fn wave_cell(c: i64, me: f64[:], fe: f64[:]) {
    for &v in me { v = 0.0; }
    for &v in fe { v = 0.0; }
    let damp: f64[2] = [1.0, 2.0 * wave_damping];
    let nu = wave_cell_terms[5 * c + 4];
    for q in rich_count() {
        let xi, w = rich_point(q);
        let dv = w * cell_geometry(c, xi) * thickness;
        var p: f64[3] = [0.0; 3];
        for l in per_cell {
            for t in 3 { p[t] += shp[l] * node_x(cell_node(c, l), t); }
        }
        let s = stretch(c, p, pml_k(c));
        var det: f64[2] = [1.0, 0.0];
        for j in dim { det = cmul(det, [s[2 * j], s[2 * j + 1]]); }
        var lam: f64[6] = [0.0; 6];
        for j in dim {
            let sj: f64[2] = [s[2 * j], s[2 * j + 1]];
            let l = cmul(damp, cdiv(det, cmul(sj, sj)));
            lam[2 * j] = nu * l[0];
            lam[2 * j + 1] = nu * l[1];
        }
        let cc: f64[2] = [wave_cell_terms[5 * c], wave_cell_terms[5 * c + 1]];
        let ss: f64[2] = [wave_cell_terms[5 * c + 2], wave_cell_terms[5 * c + 3]];
        let m = cmul(cc, det);
        let f = cmul(ss, det);
        for i in per_cell {
            fe[2 * i] -= f[0] * shp[i] * dv;
            fe[2 * i + 1] -= f[1] * shp[i] * dv;
            for j in per_cell {
                var a: f64[2] = [-m[0] * shp[i] * shp[j], -m[1] * shp[i] * shp[j]];
                for t in dim {
                    let g = grad[3 * i + t] * grad[3 * j + t];
                    a[0] += lam[2 * t] * g;
                    a[1] += lam[2 * t + 1] * g;
                }
                me[2 * (per_cell * i + j)] += a[0] * dv;
                me[2 * (per_cell * i + j) + 1] += a[1] * dv;
            }
        }
    }
}

// The matrix and the right side of the masters, in val and assembled_source.
// The natural loads of the tied nodes move to their masters.
fn wave_assemble() {
    wave_pattern();
    // The restart vectors of GMRES, and its (m + 1) × m Hessenberg matrix
    // (8 MB for m = 1000) as vectors of the dofs.
    wave_choose_path();
    assembled_source.clear();
    for d in dof_total() { assembled_source.push(0.0); }
    var me: f64[200] = [0.0; 200];
    var fe: f64[20] = [0.0; 20];
    for c in cell_count() {
        wave_cell(c, me[..2 * per_cell * per_cell], fe[..2 * per_cell]);
        for i in per_cell {
            let mi, phi = master(cell_node(c, i));
            let f = cmul(conj(phi), [fe[2 * i], fe[2 * i + 1]]);
            assembled_source[2 * mi] += f[0];
            assembled_source[2 * mi + 1] += f[1];
            for j in per_cell {
                let mj, pj = master(cell_node(c, j));
                let k = per_cell * i + j;
                add_block(mi, mj, cmul(cmul(conj(phi), [me[2 * k], me[2 * k + 1]]), pj));
            }
        }
    }
    robin_assemble();
    for n in node_count() {
        if tie[n] < 0 { continue; }
        let m, p = master(n);
        let f = cmul(conj(p), [load[2 * n], load[2 * n + 1]]);
        load[2 * m] += f[0];
        load[2 * m + 1] += f[1];
        load[2 * n] = 0.0;
        load[2 * n + 1] = 0.0;
        add_block(n, n, [1.0, 0.0]);
    }
}

// --- the solve -------------------------------------------------------------------

// Whether dof d is held: fixed on a node with no tie, or on a tied node.
fn held_dof(d: i64) -> bool { return tie[d / 2] >= 0 || fixed[d] != 0; }

// The right side, with the fixed values in it and the held rows the
// identity.
fn wave_rhs() -> f64[>..] {
    let n = sol.len;
    var b: f64[>..] = [];
    for d in n {
        let v = if tie[d / 2] >= 0 { 0.0 } else if fixed[d] != 0 { prescribed[d] } else { 0.0 };
        sol[d] = v;
        b.push(if held_dof(d) { v } else { load[d] + assembled_source[d] });
    }
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            if held_dof(r) {
                val[e] = if c == r { 1.0 } else { 0.0 };
            } else if held_dof(c) {
                b[r] -= val[e] * sol[c];
                val[e] = 0.0;
            }
        }
    }
    return b;
}

let WAVE_SECONDS = 900.0;       // the longest estimated time of a direct solve (15 minutes). Reason: a run that cannot finish stops with a message and never looks like a hang

// The solution path of a system with its pattern in row_start and col: the
// direct solve of src/solvers/slu.goose, whose factor is also the
// preconditioner of GMRES. Prints the time estimate if it is over 20 s.
// Stops before the solve, with the estimate and the fix, if the estimate is
// over WAVE_SECONDS. Checks the memory of the path. Reason: GMRES with
// ILU(0) does not converge on the indefinite systems of a fine mesh
// (waveguide_step at 1.4 mm ran for more than 5 minutes).
fn wave_choose_path() {
    guard solvers::slu_plan(row_start, col, 2) else {
        write_stderr(str("skein run: the direct solve cannot plan the system of ", dof_total(), " dofs. Use a coarser mesh\n"));
        exit(1);
    };
    let seconds = solvers::slu_work() / solvers::SLU_RATE;
    if seconds > WAVE_SECONDS {
        write_stderr(str("skein run: the direct solve of ", dof_total(), " dofs takes about ", round(seconds) as! i64,
                         " seconds, and the limit is ", round(WAVE_SECONDS) as! i64, " seconds. Use a coarser mesh\n"));
        exit(1);
    }
    solvers::memory_check("simplex wave lu", solvers::slu_entry_count(), dof_total(), 40, 0, "");
    if seconds > 20.0 {
        write_stderr(str("skein run: the direct solve of ", dof_total(), " dofs takes about ", round(seconds) as! i64, " seconds\n"));
    }
    solvers::slu_build();
}

// The solve of b with the factor of the direct solve: x in sol. False if
// GMRES does not converge.
fn wave_krylov(b: f64[:]) -> i64, bool {
    var its = 0;
    var converged = solvers::slu_factor(row_start, col, val);
    if converged {
        solvers::rows_start(row_start, col, val);
        its, converged = solvers::gmres(row_start, col, val, b, sol, WAVE_TOLERANCE, 30, 200) { r, z => solvers::slu_apply(r, z); };
        solvers::rows_stop();
    }
    solvers::slu_free();
    return its, converged;
}

// Assembles and solves the system of one frequency. False if GMRES does
// not converge.
fn wave_solve() -> bool {
    guard physics == WAVE && dofs_sized() && tie.len == node_count() else { abort("wave_solve: a wave physics and the arrays of the dofs"); }
    let start = clock();
    wave_assemble();
    let sizes = wave_sizes();
    let b = wave_rhs();
    let its, converged = wave_krylov(b);
    for n in node_count() {
        if tie[n] < 0 { continue; }
        let m, p = master(n);
        let v = cmul(p, [sol[2 * m], sol[2 * m + 1]]);
        sol[2 * n] = v[0];
        sol[2 * n + 1] = v[1];
    }
    guard wave_sizes() == sizes else { abort("wave_solve: a solve grows no array of the state (rule 3)"); }
    last_iterations = its;
    last_seconds = clock() - start;
    return converged;
}

// The length of each array of the state of the solve (rule 3).
fn wave_sizes() -> i64[8] {
    return [fixed.len, prescribed.len, load.len, sol.len, row_start.len, col.len, val.len, assembled_source.len];
}
