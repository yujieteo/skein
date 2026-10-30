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
import .beam;

var top: i64[>..] = [];
var start: i64[>..] = [];      // the index in sky of the entry (top[j], j)
var sky: f64[>..] = [];

var last_seconds = 0.0;
var singular_dof = -1;          // the first dof with no stiffness, or -1

fn cell_dofs(c: i64, beam_cell: bool) -> i64[24], i64 {
    var d: i64[24] = [0; 24];
    let nodes = if beam_cell { 2 } else { 4 };
    for l in nodes {
        let n = if beam_cell { beams[c].n[l] } else { shells[c].n[l] };
        for k in 6 { d[6 * l + k] = 6 * n + k; }
    }
    return d, 6 * nodes;
}

fn profile() {
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
    sky.clear();
    for i in s { sky.push(0.0); }
}

fn at(i: i64, j: i64) -> i64 { return start[j] + i - top[j]; }

fn scatter(d: i64[24], m: i64, k: f64[:]) {
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

// Holds each dof of a node that no cell uses: a node in a cutout.
fn fix_free_nodes() {
    var used: u8[>..] = [];
    for n in node_count() { used.push(0); }
    for c in shells { for n in c.n { used[n] = 1; } }
    for b in beams { for n in b.n { used[n] = 1; } }
    for n in node_count() { if used[n] == 0 { for d in 6 { fix(n, d); } } }
}

fn assemble() {
    mean_directors();
    fix_free_nodes();
    profile();
    var kbeam: f64[144] = [0.0; 144];
    for c in shells.len {
        shell_stiffness(c);
        let d, m = cell_dofs(c, false);
        scatter(d, m, ke);
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
    assemble();
    let ok = factor();
    if ok { substitute(); }
    last_seconds = clock() - t0;
    return ok;
}

// The reaction forces: K u - f on the fixed dofs, from each cell again.
// Returns the reaction of each dof (0 on a free dof).
fn reactions() -> f64[>..] {
    var r: f64[>..] = [];
    for d in dof_count() { r.push(0.0); }
    var kbeam: f64[144] = [0.0; 144];
    for pass in 2 {
        let count = if pass == 0 { shells.len } else { beams.len };
        for c in count {
            if pass == 0 { shell_stiffness(c); } else { beam_stiffness(c, kbeam); }
            let d, m = cell_dofs(c, pass == 1);
            for a in m {
                if fixed[d[a]] == 0 { continue; }
                var s = 0.0;
                for b in m {
                    let k = if pass == 0 { ke[m * a + b] } else { kbeam[m * a + b] };
                    s += k * disp[d[b]];
                }
                r[d[a]] += s;
            }
        }
    }
    for d in dof_count() { if fixed[d] != 0 { r[d] -= load[d]; } }
    return r;
}

// The sum of the reaction forces r (from reactions) on the nodes.
fn reaction_sum(r: f64[:], nodes: i64[:]) -> f64, f64, f64 {
    var s: f64[3] = [0.0, 0.0, 0.0];
    for n in nodes { for k in 3 { s[k] += r[6 * n + k]; } }
    return s[0], s[1], s[2];
}

// The sizes of the loads and of the reactions, for the limit of the
// equilibrium check.
fn load_scale(r: f64[:]) -> f64 {
    var s = 0.0;
    for d in load.len { if d % 6 < 3 { s += abs(load[d]) + abs(r[d]); } }
    return s;
}

fn profile_size() -> i64 { return sky.len; }
