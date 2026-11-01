// The assembly and a direct solve: a skyline (profile) LDLᵀ factorization,
// as COLSOL in K. J. Bathe, Finite Element Procedures, 1996, §8.2.3.
//
// A shell model with a small thickness has a stiffness matrix with a large
// condition number (about (L/h)²), so CG with Jacobi converges slowly. A
// direct solve does not depend on it.
//
// Column j holds the rows top[j]..j. A fixed dof keeps a 1 on the diagonal
// and no other entry, with a load of 0.
import std;
import math;
import os;
namespace shells;
import .laminate;
import .model;
import .mitc4;
import .mitc9;
import .beam;
import .mesh;
import .dissect;
import solvers from "../solvers";

var top: i64[>..] = [];
var start: i64[>..] = [];      // the index in sky of the entry (top[j], j)
var sky: f64[>..] = [];

var last_seconds = 0.0;
var singular_dof = -1;          // the first dof with no stiffness, or -1

// The direct solves: the skyline (this file) or the nested dissection of
// one grid (dissect.goose). AUTO takes the dissection for a grid of order
// 2, which has a wide skyline, and the skyline for all else.
let SOLVE_AUTO = 0;
let SOLVE_SKYLINE = 1;
let SOLVE_DISSECTION = 2;
var solve_method = SOLVE_AUTO;
var last_method = SOLVE_SKYLINE;    // the solve that solve() used last

fn cell_dofs(c: i64, beam_cell: bool) -> i64[54], i64 {
    guard c >= 0 && c < (if beam_cell { beams.len } else { shells.len }) else { abort("cell_dofs: 0 <= c < the count of beams or of shells"); }
    var d: i64[54] = [0; 54];
    let nodes = if beam_cell { 2 } else { shells[c].count };
    for l in nodes {
        let n = if beam_cell { beams[c].n[l] } else { shells[c].n[l] };
        for k in 6 { d[6 * l + k] = 6 * n + k; }
    }
    return d, 6 * nodes;
}

fn profile() {
    guard shells.len + beams.len > 0 && fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("profile: a cell, and a fixed, load and solution value for each dof"); }
    let n = dof_count();
    top.clear();
    for j in n { top.push(j); }
    for pass in 2 {
        let count = if pass == 0 { shells.len } else { beams.len };
        for c in count {
            let d, m = cell_dofs(c, pass == 1);
            var lo = d[0];
            for i in m { lo = min(lo, d[i]); }
            for i in m { top[d[i]] = min(top[d[i]], lo); }
        }
    }
    start.clear();
    var s = 0;
    for j in n {
        start.push(s);
        s += j - top[j] + 1;
    }
    solvers::memory_check("shell skyline ldlt", s, n, 1, 0, "");
    sky.clear();
    for i in s { sky.push(0.0); }
}

fn at(i: i64, j: i64) -> i64 {
    guard j >= 0 && j < top.len && i >= top[j] && i <= j else { abort("at: (i, j) is in the profile"); }
    return start[j] + i - top[j];
}

fn scatter(d: i64[54], m: i64, k: f64[:]) {
    guard (m == 12 || m == 24 || m == 54) && k.len >= m * m && top.len == fixed.len else { abort("scatter: m is 12, 24 or 54, k has m x m entries, and the profile is made"); }
    for a in m {
        let i = d[a];
        if fixed[i] != 0 { continue; }
        for b in m {
            let j = d[b];
            if fixed[j] != 0 || i > j { continue; }
            sky[at(i, j)] += k[m * a + b];
        }
    }
}

// The stiffness of shell cell c: ke of MITC4 or ke9 of MITC9.
fn cell_stiffness(c: i64) {
    guard c >= 0 && c < shells.len else { abort("cell_stiffness: 0 <= c < shells.len"); }
    if shells[c].count == 9 { mitc9_stiffness(c); } else { shell_stiffness(c); }
}

// Holds each dof of a node that no cell uses: a node in a cutout.
fn fix_free_nodes() {
    var used: u8[>..] = [];
    for n in node_count() { used.push(0); }
    for c in shells { for l in c.count { used[c.n[l]] = 1; } }
    for b in beams { for n in b.n { used[n] = 1; } }
    for n in node_count() { if used[n] == 0 { for d in 6 { fix(n, d); } } }
}

fn assemble() {
    guard directors.len == coords.len && fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("assemble: a director, fixed, load and solution value for each node and dof"); }
    mean_directors();
    fix_free_nodes();
    profile();
    var kbeam: f64[144] = [0.0; 144];
    for c in shells.len {
        cell_stiffness(c);
        let d, m = cell_dofs(c, false);
        scatter(d, m, if m == 54 { ke9[..] } else { ke[..] });
    }
    for b in beams.len {
        beam_stiffness(b, kbeam);
        let d, m = cell_dofs(b, true);
        scatter(d, m, kbeam);
    }
    for j in dof_count() { if fixed[j] != 0 { sky[at(j, j)] = 1.0; } }
}

// sky = L D Lᵀ in place: the diagonal holds D, column j above it holds Lᵀ.
fn factor() -> bool {
    guard top.len == dof_count() && start.len == top.len else { abort("factor: the profile has each dof"); }
    let n = dof_count();
    singular_dof = -1;
    for j in n {
        let tj = top[j];
        let bj = start[j] - tj;
        for i in tj + 1..j {
            let bi = start[i] - top[i];
            var s = sky[bj + i];
            for r in max(top[i], tj)..i { s -= sky[bi + r] * sky[bj + r]; }
            sky[bj + i] = s;
        }
        var d = sky[at(j, j)];
        let scale_j = abs(d);
        for r in tj..j {
            let g = sky[at(r, j)];
            let l = g / sky[at(r, r)];
            sky[at(r, j)] = l;
            d -= l * g;
        }
        if !(d > 1.0e-12 * scale_j) {
            singular_dof = j;
            return false;
        }
        sky[at(j, j)] = d;
    }
    return true;
}

fn substitute() {
    guard top.len == dof_count() && start.len == top.len && fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("substitute: the profile, the loads and the solution have each dof"); }
    let n = dof_count();
    for j in n {
        disp[j] = if fixed[j] != 0 { 0.0 } else { load[j] };
    }
    for j in n {
        var s = disp[j];
        for i in top[j]..j { s -= sky[at(i, j)] * disp[i]; }
        disp[j] = s;
    }
    for j in n { disp[j] /= sky[at(j, j)]; }
    for jj in n {
        let j = n - 1 - jj;
        let x = disp[j];
        for i in top[j]..j { disp[i] -= sky[at(i, j)] * x; }
    }
}

// Assembles and solves K u = f. Returns false if the matrix is singular:
// a mechanism or a dof with no stiffness (singular_dof).
fn solve() -> bool {
    let t0 = clock();
    let g = dissect_grid();
    guard solve_method != SOLVE_DISSECTION || g >= 0 else { abort("solve: the dissection solves only a model of one grid, with each beam in a cell"); }
    let dissect = g >= 0 && (solve_method == SOLVE_DISSECTION || (solve_method == SOLVE_AUTO && grids[g].order == 2));
    last_method = if dissect { SOLVE_DISSECTION } else { SOLVE_SKYLINE };
    if dissect {
        guard directors.len == coords.len && fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("solve: a director, fixed, load and solution value for each node and dof"); }
        mean_directors();
        fix_free_nodes();
        dissect_prepare(g);
    } else {
        assemble();
    }
    let sizes = state_sizes();
    let ok = if dissect { dissect_factor(g) } else { factor() };
    if ok && dissect { dissect_substitute(); }
    if ok && !dissect { substitute(); }
    guard state_sizes() == sizes else { abort("solve: a solve grows no array of the state (rule 3)"); }
    last_seconds = clock() - t0;
    return ok;
}

// K u - f on each dof, from each cell again.
fn out_of_balance() -> f64[>..] {
    guard fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("out_of_balance: a fixed, load and solution value for each dof"); }
    var r: f64[>..] = [];
    for d in dof_count() { r.push(-load[d]); }
    var kbeam: f64[144] = [0.0; 144];
    for pass in 2 {
        let count = if pass == 0 { shells.len } else { beams.len };
        for c in count {
            if pass == 0 { cell_stiffness(c); } else { beam_stiffness(c, kbeam); }
            let d, m = cell_dofs(c, pass == 1);
            for a in m {
                var s = 0.0;
                for b in m {
                    let k = if pass == 1 { kbeam[m * a + b] } else if m == 54 { ke9[m * a + b] } else { ke[m * a + b] };
                    s += k * disp[d[b]];
                }
                r[d[a]] += s;
            }
        }
    }
    return r;
}

// The reaction forces: K u - f on the fixed dofs. Returns the reaction of
// each dof (0 on a free dof).
fn reactions() -> f64[>..] {
    var r = out_of_balance();
    for d in dof_count() { if fixed[d] == 0 { r[d] = 0.0; } }
    return r;
}

// |K u - f| / |f| on the free dofs: the check of a solve.
fn relative_residual() -> f64 {
    let r = out_of_balance();
    var rr = 0.0;
    var ff = 0.0;
    for d in dof_count() {
        if fixed[d] != 0 { continue; }
        rr += r[d] * r[d];
        ff += load[d] * load[d];
    }
    return sqrt(rr / max(ff, 1.0e-300));
}

// The sum of the reaction forces r (from reactions) on the nodes.
fn reaction_sum(r: f64[:], nodes: i64[:]) -> f64, f64, f64 {
    guard r.len == dof_count() else { abort("reaction_sum: r has each dof"); }
    var s: f64[3] = [0.0, 0.0, 0.0];
    for n in nodes { for k in 3 { s[k] += r[6 * n + k]; } }
    return s[0], s[1], s[2];
}

// The sizes of the loads and of the reactions, for the limit of the
// equilibrium check.
fn load_scale(r: f64[:]) -> f64 {
    guard r.len == load.len else { abort("load_scale: r and load have one size"); }
    var s = 0.0;
    for d in load.len { if d % 6 < 3 { s += abs(load[d]) + abs(r[d]); } }
    return s;
}

fn profile_size() -> i64 { return sky.len; }

// The length of each array of the state of the package. A solve grows none
// of them (rule 3 of the Power of Ten, docs/tradeoffs.md §6).
fn state_sizes() -> i64[27] {
    return [coords.len, directors.len, shells.len, beams.len, beam_sections.len, fixed.len, load.len,
            disp.len, grids.len, sections.len, plies.len, top.len, start.len, sky.len, grid_lines.len, grid_cells.len,
            fronts.len, front_nodes.len, lfactor.len, front_buf.len, front_panel.len, front_stack.len, front_map.len,
            front_pos.len, kdiag.len, leaf_beam_start.len, leaf_beams.len];
}
