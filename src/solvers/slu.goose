// A direct solve of a sparse system with a symmetric pattern and no
// lattice: the LU factors of A without pivoting, on a nested dissection
// of the graph of A. It is the solve of the indefinite systems of the
// wave and edge physics (a complex dof is 2 dofs, `block` = 2), where
// GMRES with ILU(0) does not converge on a fine mesh. References: A.
// George (1973), "Nested dissection of a regular finite element mesh",
// SIAM J. Numer. Anal. 10, 345-363 (here with the level sets of a
// breadth-first search as the separators), and T. A. Davis (2006), Direct
// Methods for Sparse Linear Systems, chapters 4 and 6 (the elimination
// tree, the pattern of each row, and the up-looking factorization).
//
// slu_plan orders the blocks of dofs and counts the entries of L and the
// multiplies of the factorization, so that the caller checks the memory
// and the time before slu_build sizes the factors. slu_factor puts the
// values in. A pivot that is almost 0 takes the size 1e-14 of the largest
// entry of its row, so the factor is a preconditioner of GMRES even for
// a matrix that has no LU without pivoting.
import std;
import math;
namespace solvers;

let SLU_LEAF = 96;          // the blocks of a part that the dissection does not cut
let SLU_PIVOT = 1.0e-14;    // the smallest pivot, relative to the largest entry of its row
let SLU_RATE = 4.0e8;       // Reason: the multiplies and adds for each second of slu_factor on the desktop (waveguide_step, 1.4 mm)

var slu_n = 0;                   // the dofs
var slu_perm: i64[>..] = [];     // the dof of each place of the order
var slu_inv: i64[>..] = [];      // the place of each dof
var slu_lo_start: i64[>..] = []; // the entries below the diagonal of the permuted pattern, by row
var slu_lo_col: i32[>..] = [];
var slu_parent: i64[>..] = [];   // the elimination tree
var slu_count: i64[>..] = [];    // the entries of each column of L, then (for the rows) of each row
var slu_rowcount: i64[>..] = [];
var slu_lc_start: i64[>..] = []; // L by columns: the rows of column k, and the columns of U by rows
var slu_lc_row: i32[>..] = [];
var slu_lr_start: i64[>..] = []; // L by rows: the columns of row i
var slu_lr_col: i32[>..] = [];
var slu_lval: f64[>..] = [];     // L, in the order of slu_lr_col
var slu_uval: f64[>..] = [];     // U without the diagonal, in the order of slu_lc_row
var slu_diag: f64[>..] = [];
var slu_w: f64[>..] = [];
var slu_x: f64[>..] = [];
var slu_entries = 0;             // the entries of L (and of U)
var slu_flops = 0.0;             // the multiplies and adds of slu_factor
var slu_small = 0;               // the pivots that slu_factor replaced

fn slu_entry_count() -> i64 { return slu_entries; }
fn slu_work() -> f64 { return slu_flops; }
fn slu_replaced() -> i64 { return slu_small; }

fn slu_ints(v: i64[>..]&, n: i64, x: i64) {
    v.clear();
    for i in n { v.push(x); }
}

// --- the order -----------------------------------------------------------------------------

var dis_start: i64[>..] = [];    // the neighbours of each width
var dis_adj: i64[>..] = [];
var dis_part: i64[>..] = [];     // the part that holds each width
var dis_level: i64[>..] = [];
var dis_queue: i64[>..] = [];
var dis_vs: i64[>..] = [];       // the blocks of the parts, one part after the other
var dis_order: i64[>..] = [];    // the width at each place
var dis_stack: i64[>..] = [];    // parts to cut: first in vs, end in vs, first place

fn slu_block_graph(row_start: i64[:], col: i32[:], width: i64) -> i64 {
    let nb = (row_start.len - 1) / width;
    var mark: i64[>..] = [];
    slu_ints(mark, nb, -1);
    dis_start.clear();
    dis_adj.clear();
    for b in nb {
        dis_start.push(dis_adj.len);
        for d in width {
            let r = b * width + d;
            for e in row_start[r]..row_start[r + 1] {
                let c = col[e] / width;
                if c != b && mark[c] != b {
                    mark[c] = b;
                    dis_adj.push(c);
                }
            }
        }
    }
    dis_start.push(dis_adj.len);
    return nb;
}

// The breadth-first search from `root` in the part `id` (the blocks
// vs[s..e]). The queue holds the blocks that it reached, level by level.
// Returns the number of blocks and the number of levels.
fn dis_bfs(root: i64, id: i64, s: i64, e: i64) -> i64, i64 {
    for t in s..e { dis_level[dis_vs[t]] = -1; }
    dis_queue.clear();
    dis_queue.push(root);
    dis_level[root] = 0;
    var head = 0;
    var levels = 1;
    while head < dis_queue.len {
        let v = dis_queue[head];
        head++;
        for k in dis_start[v]..dis_start[v + 1] {
            let w = dis_adj[k];
            if dis_part[w] == id && dis_level[w] < 0 {
                dis_level[w] = dis_level[v] + 1;
                levels = max(levels, dis_level[w] + 1);
                dis_queue.push(w);
            }
        }
    }
    return dis_queue.len, levels;
}

// The level that cuts the part in two: of the levels that leave 1/4 to
// 3/4 of the blocks before them, the smallest.
fn dis_separator(levels: i64, size: i64) -> i64 {
    var count: i64[>..] = [];
    slu_ints(count, levels, 0);
    for v in dis_queue { count[dis_level[v]]++; }
    var before = count[0];
    var best = -1;
    var near = 1;
    var gap = size;
    for l in 1..levels - 1 {
        if 4 * before >= size && 4 * before <= 3 * size && (best < 0 || count[l] < count[best]) { best = l; }
        if abs(2 * before - size) < gap {
            gap = abs(2 * before - size);
            near = l;
        }
        before += count[l];
    }
    return if best >= 0 { best } else { near };
}

// Puts the blocks of the queue at the places lo.. in the order of the queue.
fn dis_leaf(s: i64, lo: i64) {
    for t in dis_queue.len {
        dis_vs[s + t] = dis_queue[t];
        dis_order[lo + t] = dis_queue[t];
    }
}

// Cuts the part at the level `ls`: the blocks before it and after it go on the stack,
// and the blocks of the level take the last places of the part.
fn dis_cut(s: i64, lo: i64, ls: i64) {
    let size = dis_queue.len;
    var a = 0;
    var sep = 0;
    for v in dis_queue {
        if dis_level[v] < ls { a++; }
        if dis_level[v] == ls { sep++; }
    }
    let b = size - a - sep;
    var ia = 0;
    var ib = 0;
    var is = 0;
    for v in dis_queue {
        if dis_level[v] < ls {
            dis_vs[s + ia] = v;
            ia++;
        } else if dis_level[v] == ls {
            dis_order[lo + a + b + is] = v;
            is++;
        } else {
            dis_vs[s + a + ib] = v;
            ib++;
        }
    }
    dis_stack.push(s + a);
    dis_stack.push(s + a + b);
    dis_stack.push(lo + a);
    dis_stack.push(s);
    dis_stack.push(s + a);
    dis_stack.push(lo);
}

// Cuts one part (the blocks vs[s..e], from place lo) or orders it.
fn dis_part_step(s: i64, e: i64, lo: i64, id: i64) {
    let size = e - s;
    for t in s..e { dis_part[dis_vs[t]] = id; }
    var reached, levels = dis_bfs(dis_vs[s], id, s, e);
    let far = dis_queue[reached - 1];
    reached, levels = dis_bfs(far, id, s, e);
    if reached < size {
        // A part with 2 pieces: the reached blocks, and the others.
        for t in s..e {
            if dis_level[dis_vs[t]] < 0 { dis_queue.push(dis_vs[t]); }
        }
        for t in size { dis_vs[s + t] = dis_queue[t]; }
        dis_stack.push(s);
        dis_stack.push(s + reached);
        dis_stack.push(lo);
        dis_stack.push(s + reached);
        dis_stack.push(e);
        dis_stack.push(lo + reached);
        return;
    }
    if size <= SLU_LEAF || levels < 3 {
        dis_leaf(s, lo);
        return;
    }
    dis_cut(s, lo, dis_separator(levels, size));
}

// The order of the nb blocks.
fn slu_dissect(nb: i64) {
    slu_ints(dis_part, nb, -1);
    slu_ints(dis_level, nb, -1);
    slu_ints(dis_order, nb, 0);
    dis_vs.clear();
    for b in nb { dis_vs.push(b); }
    dis_stack.clear();
    dis_stack.push(0);
    dis_stack.push(nb);
    dis_stack.push(0);
    var id = 0;
    while dis_stack.len > 0 {
        let lo = dis_stack[dis_stack.len - 1];
        let e = dis_stack[dis_stack.len - 2];
        let s = dis_stack[dis_stack.len - 3];
        dis_stack.pop();
        dis_stack.pop();
        dis_stack.pop();
        if e > s {
            dis_part_step(s, e, lo, id);
            id++;
        }
    }
}

// --- the pattern of L ------------------------------------------------------------------------

// The entries below the diagonal of the permuted pattern, from A and from
// its transpose, by row (some twice).
fn slu_lower(row_start: i64[:], col: i32[:]) {
    let n = slu_n;
    var fill: i64[>..] = [];
    slu_ints(fill, n, 0);
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let pr = slu_inv[r];
            let pc = slu_inv[col[e]];
            if pc < pr { fill[pr]++; } else if pc > pr { fill[pc]++; }
        }
    }
    slu_lo_start.clear();
    slu_lo_start.push(0);
    for i in n { slu_lo_start.push(slu_lo_start[i] + fill[i]); }
    slu_lo_col.clear();
    for t in slu_lo_start[n] { slu_lo_col.push(0); }
    for i in n { fill[i] = slu_lo_start[i]; }
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let pr = slu_inv[r];
            let pc = slu_inv[col[e]];
            if pc < pr {
                slu_lo_col[fill[pr]] = pc as i32;
                fill[pr]++;
            } else if pc > pr {
                slu_lo_col[fill[pc]] = pr as i32;
                fill[pc]++;
            }
        }
    }
}

// Visits the columns of the pattern of row i of L: with `into` < 0 it counts
// them, else it puts i in the column. Returns the entries of the row.
fn slu_reach(i: i64, mark: i64[:], into: bool, at: i64[:]) -> i64 {
    var found = 0;
    mark[i] = i;
    for t in slu_lo_start[i]..slu_lo_start[i + 1] {
        var k = slu_lo_col[t] as i64;
        while mark[k] != i {
            if into {
                slu_lc_row[at[k]] = i as i32;
                at[k]++;
            } else {
                slu_count[k]++;
            }
            found++;
            mark[k] = i;
            k = slu_parent[k];
        }
    }
    return found;
}

// The elimination tree, and the number of entries in each column and row of L.
fn slu_tree() {
    let n = slu_n;
    var ancestor: i64[>..] = [];
    var mark: i64[>..] = [];
    var nothing: i64[>..] = [];
    slu_ints(ancestor, n, -1);
    slu_ints(slu_parent, n, -1);
    slu_ints(mark, n, -1);
    slu_ints(slu_count, n, 0);
    slu_rowcount.clear();
    for i in n {
        for t in slu_lo_start[i]..slu_lo_start[i + 1] {
            var a = slu_lo_col[t] as i64;
            while a != -1 && a < i {
                let next = ancestor[a];
                ancestor[a] = i;
                if next == -1 { slu_parent[a] = i; }
                a = next;
            }
        }
        slu_rowcount.push(slu_reach(i, mark, false, nothing));
    }
    slu_entries = 0;
    slu_flops = 0.0;
    for k in n {
        slu_entries += slu_count[k];
        slu_flops += 2.0 * slu_count[k] * slu_count[k];
    }
}

// Plans the factorization of A, with `width` dofs for each node of the
// graph: the order, the entries of L and the work. False if A has no
// rows.
fn slu_plan(row_start: i64[:], col: i32[:], width: i64) -> bool {
    guard row_start.len > 1 && width > 0 && (row_start.len - 1) % width == 0 && row_start[row_start.len - 1] == col.len else {
        abort("slu_plan: row_start ends at col.len, and the rows are a whole number of blocks");
    }
    slu_n = row_start.len - 1;
    let nb = slu_block_graph(row_start, col, width);
    slu_dissect(nb);
    slu_perm.clear();
    for p in nb {
        for d in width { slu_perm.push(dis_order[p] * width + d); }
    }
    slu_ints(slu_inv, slu_n, 0);
    for i in slu_n { slu_inv[slu_perm[i]] = i; }
    slu_lower(row_start, col);
    slu_tree();
    return true;
}

// Sizes the factors: L by columns and by rows.
fn slu_build() {
    let n = slu_n;
    slu_lc_start.clear();
    slu_lc_start.push(0);
    for k in n { slu_lc_start.push(slu_lc_start[k] + slu_count[k]); }
    slu_lr_start.clear();
    slu_lr_start.push(0);
    for i in n { slu_lr_start.push(slu_lr_start[i] + slu_rowcount[i]); }
    slu_lc_row.clear();
    slu_lr_col.clear();
    for e in slu_entries {
        slu_lc_row.push(0);
        slu_lr_col.push(0);
    }
    var at: i64[>..] = [];
    var mark: i64[>..] = [];
    slu_ints(mark, n, -1);
    for k in n { at.push(slu_lc_start[k]); }
    for i in n { slu_reach(i, mark, true, at); }
    var into: i64[>..] = [];
    for i in n { into.push(slu_lr_start[i]); }
    for k in n {
        for u in slu_lc_start[k]..slu_lc_start[k + 1] {
            let i = slu_lc_row[u] as i64;
            slu_lr_col[into[i]] = k as i32;
            into[i]++;
        }
    }
    zeros(slu_lval, slu_entries);
    zeros(slu_uval, slu_entries);
    zeros(slu_diag, n);
    zeros(slu_w, n);
    zeros(slu_x, n);
}

// --- the values -----------------------------------------------------------------------------

// Row i of L and U from row i of A (the up-looking factorization). False if
// a row of A is 0.
fn slu_row(i: i64, row_start: i64[:], col: i32[:], val: f64[:]) -> bool {
    let r = slu_perm[i];
    var scale = 0.0;
    for e in row_start[r]..row_start[r + 1] {
        slu_w[slu_inv[col[e]]] += val[e];
        scale = max(scale, abs(val[e]));
    }
    guard scale > 0.0 else { return false; }
    for t in slu_lr_start[i]..slu_lr_start[i + 1] {
        let k = slu_lr_col[t] as i64;
        let l = slu_w[k] / slu_diag[k];
        slu_w[k] = 0.0;
        slu_lval[t] = l;
        if l == 0.0 { continue; }
        for u in slu_lc_start[k]..slu_lc_start[k + 1] { slu_w[slu_lc_row[u]] -= l * slu_uval[u]; }
    }
    var p = slu_w[i];
    slu_w[i] = 0.0;
    if abs(p) < SLU_PIVOT * scale {
        p = if p < 0.0 { -SLU_PIVOT * scale } else { SLU_PIVOT * scale };
        slu_small++;
    }
    slu_diag[i] = p;
    for u in slu_lc_start[i]..slu_lc_start[i + 1] {
        slu_uval[u] = slu_w[slu_lc_row[u]];
        slu_w[slu_lc_row[u]] = 0.0;
    }
    return true;
}

// Factors A. False if a row of A is 0.
fn slu_factor(row_start: i64[:], col: i32[:], val: f64[:]) -> bool {
    guard slu_n == row_start.len - 1 && slu_diag.len == slu_n else { abort("slu_factor: slu_plan and slu_build come first, for this matrix"); }
    slu_small = 0;
    for i in slu_n {
        if !slu_row(i, row_start, col, val) { return false; }
    }
    return true;
}

// z = A⁻¹ r with the factors.
fn slu_apply(r: f64[:], z: f64[:]) {
    guard r.len == slu_n && z.len == slu_n else { abort("slu_apply: r and z have one entry for each dof"); }
    let n = slu_n;
    for i in n {
        var s = r[slu_perm[i]];
        for t in slu_lr_start[i]..slu_lr_start[i + 1] { s -= slu_lval[t] * slu_x[slu_lr_col[t]]; }
        slu_x[i] = s;
    }
    for q in n {
        let i = n - 1 - q;
        var s = slu_x[i];
        for u in slu_lc_start[i]..slu_lc_start[i + 1] { s -= slu_uval[u] * slu_x[slu_lc_row[u]]; }
        slu_x[i] = s / slu_diag[i];
    }
    for i in n { z[slu_perm[i]] = slu_x[i]; }
}

// Frees the arrays of the factors and of the plan.
fn slu_free() {
    slu_n = 0;
    slu_entries = 0;
    slu_lval.clear();
    slu_uval.clear();
    slu_diag.clear();
    slu_w.clear();
    slu_x.clear();
    slu_lc_row.clear();
    slu_lr_col.clear();
    slu_lo_col.clear();
}
