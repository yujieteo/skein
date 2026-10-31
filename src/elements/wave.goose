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
// the row of the tied node is the identity. GMRES with ILU(0) solves the
// real system, then the tied nodes take their values.
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
let WAVE_RESTART = 1000;        // the largest restart: GMRES(200) stalls on a PML duct of 8 rows of cells
let WAVE_KRYLOV = 50000000;     // the largest restart times the dofs: 400 MB of Krylov vectors
let WAVE_ITERATIONS = 20000;
let PML_REFLECTION = 1.0e-6;    // the reflection of the layer at normal incidence

var wave_c: f64[2] = [0.0, 0.0];   // c, the coefficient of u
var wave_s: f64[2] = [0.0, 0.0];   // s, the source
var wave_damping = 0.0;            // ζ
var pml_shape = -1;                // the shape of the layer, or -1
var pml_lo: f64[3] = [0.0; 3];     // the rest of the body: the box from pml_lo to pml_hi
var pml_hi: f64[3] = [0.0; 3];
var pml_h = 0.0;                   // the thickness of the layer
var tie: i64[>..] = [];            // the master of each node, or -1
var tie_phase: f64[>..] = [];      // cos θ and sin θ of each node
var low_nodes: i64[>..] = [];      // the nodes of the low region of the next bloch

fn cmul(a: f64[2], b: f64[2]) -> f64[2] { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }

fn cdiv(a: f64[2], b: f64[2]) -> f64[2] {
    let d = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
}

fn conj(a: f64[2]) -> f64[2] { return [a[0], -a[1]]; }

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
}

// The layer of shape g, past the box from lo to hi, of thickness h.
fn set_pml(g: i64, lo: f64[3], hi: f64[3], h: f64) {
    guard physics == WAVE && g >= 0 && h > 0.0 && wave_c[0] > 0.0 else { abort("set_pml: a wave physics with c > 0, a shape and h > 0"); }
    pml_shape = g;
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

// The stretch s_j of each axis at point p of cell c, in the PML.
fn stretch(c: i64, p: f64[3]) -> f64[6] {
    var s: f64[6] = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
    guard pml_shape >= 0 && cell_in(c, pml_shape) else { return s; }
    let k = sqrt(wave_c[0]);
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
    for q in rich_count() {
        let xi, w = rich_point(q);
        let dv = w * cell_geometry(c, xi) * thickness;
        var p: f64[3] = [0.0; 3];
        for l in per_cell {
            for t in 3 { p[t] += shp[l] * node_x(cell_node(c, l), t); }
        }
        let s = stretch(c, p);
        var det: f64[2] = [1.0, 0.0];
        for j in dim { det = cmul(det, [s[2 * j], s[2 * j + 1]]); }
        var lam: f64[6] = [0.0; 6];
        for j in dim {
            let sj: f64[2] = [s[2 * j], s[2 * j + 1]];
            let l = cmul(damp, cdiv(det, cmul(sj, sj)));
            lam[2 * j] = l[0];
            lam[2 * j + 1] = l[1];
        }
        let m = cmul(wave_c, det);
        let f = cmul(wave_s, det);
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
    solvers::memory_check("simplex wave gmres ilu0", col.len, dof_total(), wave_restart(dof_total()) + 8, "");
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

// The GMRES restart of a system of n dofs.
fn wave_restart(n: i64) -> i64 { return min(WAVE_RESTART, max(50, WAVE_KRYLOV / max(n, 1))); }

// Assembles and solves the system of one frequency. False if GMRES does
// not converge.
fn wave_solve() -> bool {
    guard physics == WAVE && dofs_sized() && tie.len == node_count() else { abort("wave_solve: a wave physics and the arrays of the dofs"); }
    let start = clock();
    wave_assemble();
    let sizes = wave_sizes();
    let b = wave_rhs();
    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    var converged = solvers::ilu0(row_start, col, val, lu, diag);
    var its = 0;
    if converged {
        let restart = wave_restart(b.len);
        solvers::rows_start(row_start, col, val);
        its, converged = solvers::gmres(row_start, col, val, b, sol, WAVE_TOLERANCE, restart, WAVE_ITERATIONS) { r, z =>
            solvers::ilu0_apply(row_start, col, lu, diag, r, z);
        };
        solvers::rows_stop();
    }
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
