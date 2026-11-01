// Smoothed aggregation AMG (Vaněk, Mandel and Brezina 1996), the
// preconditioner of CG for a symmetric positive definite matrix (spec
// §18.4, docs/implementation.md §6).
//
// amg_setup builds the levels. On each level:
//
// 1. The dofs are in blocks of `bs` dofs, one block for each node
//    (bs = 3 for a 3D displacement). Node j is strong for node i if the
//    sum of the squares of the block A_ij is more than θ² times
//    sqrt(|A_ii|² |A_jj|²), with the same sums, or at least RELATIVE²
//    times the largest sum of a block of row i. θ is THETA on the fine
//    level and halves on each coarser level. The relative test keeps the
//    coarse levels of P2 cells coarsening: their rows couple many nodes,
//    each by a small part of the diagonal.
// 2. The aggregates, in three passes: a node whose strong neighbours are
//    all free makes an aggregate with them; a free node joins the
//    aggregate of a strong neighbour; a node still free makes an aggregate
//    with its free strong neighbours.
// 3. The tentative prolongator P0 holds the near null space B on each
//    aggregate: B of the rows of the aggregate is Q R (Gram-Schmidt), the
//    columns of Q are the columns of P0 for the aggregate, and R is B of
//    the coarse level. Each coarse node has a block of nb dofs, one for
//    each column of B. For elasticity B is the rigid-body modes (nb = 6 in
//    3D, 3 in 2D); amg_setup takes the constant of each component
//    (nb = bs). A column that Gram-Schmidt finds dependent on the others
//    is zero: its coarse dof has an empty row, which the last level holds
//    at 0.
// 4. P = (I - ω D⁻¹ A) P0, with ω = 4 / (3 ρ) and ρ the largest eigenvalue
//    of D⁻¹ A, from power steps. The next level is Pᵀ A P.
//
// The levels stop at COARSE dofs or fewer, or when the aggregates do not
// make the level smaller. A dense Cholesky factor solves the last level.
// amg_apply is one V-cycle: a forward Gauss-Seidel sweep before the coarse
// correction and a backward sweep after it, so the preconditioner is
// symmetric. The levels are loops with a fixed maximum, not recursion
// (docs/tradeoffs.md §6).
//
// The levels are the state of the package between amg_setup and the
// calls of amg_apply. Each level l is in the flat arrays below, from the
// offsets of level l.
import std;
import math;
namespace solvers;
import .sparse;

let AMG_MAX_LEVELS = 12;
let AMG_COARSE = 600;
let THETA = 0.08;
let RELATIVE = 0.5;
let POWER_STEPS = 15;

var amg_levels = 0;
var amg_n: i64[>..] = [];       // the dofs of each level
var amg_rs: i64[>..] = [];      // the row starts of level l from amg_rs_at[l]
var amg_rs_at: i64[>..] = [];
var amg_col: i32[>..] = [];     // the entries of level l from amg_e_at[l] to amg_e_at[l + 1]
var amg_val: f64[>..] = [];
var amg_e_at: i64[>..] = [];
var amg_prs: i64[>..] = [];     // the prolongator from level l + 1 to level l, the same way
var amg_prs_at: i64[>..] = [];
var amg_pcol: i32[>..] = [];
var amg_pval: f64[>..] = [];
var amg_pe_at: i64[>..] = [];
var amg_x: f64[>..] = [];       // the vectors of level l from amg_v_at[l]
var amg_b: f64[>..] = [];
var amg_res: f64[>..] = [];
var amg_diag: f64[>..] = [];
var amg_v_at: i64[>..] = [];
var amg_chol: f64[>..] = [];    // the dense Cholesky factor of the last level, by rows

fn level_rs(l: i64) -> const i64[:] { return amg_rs[amg_rs_at[l]..amg_rs_at[l] + amg_n[l] + 1]; }
fn level_col(l: i64) -> const i32[:] { return amg_col[amg_e_at[l]..amg_e_at[l + 1]]; }
fn level_val(l: i64) -> const f64[:] { return amg_val[amg_e_at[l]..amg_e_at[l + 1]]; }

// Appends a level of n dofs with the matrix (rs, col, val).
fn push_level(n: i64, rs: i64[:], col: i32[:], val: f64[:]) {
    guard rs.len == n + 1 && col.len == val.len && rs[n] == col.len else { abort("push_level: n + 1 row starts that end at col.len"); }
    amg_n.push(n);
    amg_rs_at.push(amg_rs.len);
    amg_rs.append(rs);
    if amg_e_at.len == 0 { amg_e_at.push(0); }
    amg_col.append(col);
    amg_val.append(val);
    amg_e_at.push(amg_col.len);
    amg_v_at.push(amg_x.len);
    for r in n {
        amg_x.push(0.0);
        amg_b.push(0.0);
        amg_res.push(0.0);
        var d = 0.0;
        for e in rs[r]..rs[r + 1] { if col[e] == r { d = val[e]; } }
        amg_diag.push(if d > 0.0 { d } else { 1.0 });
    }
    amg_levels = amg_n.len;
}

// The strong neighbours of each node of level l: adj[adj_start[i]..adj_start[i + 1]].
fn strong_graph(l: i64, bs: i64) -> i64[>..], i64[>..] {
    let rs = level_rs(l);
    let col = level_col(l);
    let val = level_val(l);
    let nodes = amg_n[l] / bs;
    // The sum of squares of each block A_ij of node row i, in `sum`.
    var sum: f64[>..] = [];
    var seen: i64[>..] = [];
    for j in nodes {
        sum.push(0.0);
        seen.push(-1);
    }
    var self_sum: f64[>..] = [];
    var strong_start: i64[>..] = [0];
    var strong: i64[>..] = [];
    var strong_sum: f64[>..] = [];
    var touched: i64[>..] = [];
    for i in nodes {
        touched.clear();
        for r in i * bs..(i + 1) * bs {
            for e in rs[r]..rs[r + 1] {
                let j = col[e] / bs;
                if seen[j] != i {
                    seen[j] = i;
                    sum[j] = 0.0;
                    touched.push(j);
                }
                sum[j] += val[e] * val[e];
            }
        }
        self_sum.push(if seen[i] == i { sum[i] } else { 0.0 });
        for j in touched {
            if j != i {
                strong.push(j);
                strong_sum.push(sum[j]);
            }
        }
        strong_start.push(strong.len);
    }
    // Keep the strong pairs. The limit halves at each level, as in Vaněk,
    // Mandel and Brezina: the coarse levels couple more nodes, each less.
    let theta = THETA * pow(0.5, l as f64);
    var adj_start: i64[>..] = [0];
    var adj: i64[>..] = [];
    for i in nodes {
        var most = 0.0;
        for k in strong_start[i]..strong_start[i + 1] { most = max(most, strong_sum[k]); }
        for k in strong_start[i]..strong_start[i + 1] {
            let j = strong[k];
            if strong_sum[k] > theta * theta * sqrt(self_sum[i] * self_sum[j]) || strong_sum[k] >= RELATIVE * RELATIVE * most { adj.push(j); }
        }
        adj_start.push(adj.len);
    }
    return adj_start, adj;
}

// The aggregate of each node of level l, and the count of aggregates.
fn aggregates(l: i64, bs: i64) -> i64[>..], i64 {
    let nodes = amg_n[l] / bs;
    let adj_start, adj = strong_graph(l, bs);
    var agg: i64[>..] = [];
    for i in nodes { agg.push(-1); }
    var count = 0;
    for i in nodes {
        guard agg[i] < 0 && adj_start[i + 1] > adj_start[i] else { continue; }
        var free = true;
        for j in adj[adj_start[i]..adj_start[i + 1]] { if agg[j] >= 0 { free = false; } }
        guard free else { continue; }
        agg[i] = count;
        for j in adj[adj_start[i]..adj_start[i + 1]] { agg[j] = count; }
        count++;
    }
    var joined: i64[>..] = [];
    joined.append(agg);
    for i in nodes {
        guard agg[i] < 0 else { continue; }
        for j in adj[adj_start[i]..adj_start[i + 1]] {
            if agg[j] >= 0 {
                joined[i] = agg[j];
                break;
            }
        }
    }
    agg.clear();
    agg.append(joined);
    for i in nodes {
        guard agg[i] < 0 else { continue; }
        agg[i] = count;
        for j in adj[adj_start[i]..adj_start[i + 1]] { if agg[j] < 0 { agg[j] = count; } }
        count++;
    }
    return agg, count;
}

// The largest eigenvalue of D⁻¹ A of level l, from power steps.
fn spectral_radius(l: i64) -> f64 {
    let n = amg_n[l];
    let d = amg_diag[amg_v_at[l]..amg_v_at[l] + n];
    var x: f64[>..] = [];
    var y: f64[>..] = [];
    for i in n {
        x.push(1.0 + ((i * 7919) % 13) / 13.0);
        y.push(0.0);
    }
    var rho = 1.0;
    for k in POWER_STEPS {
        let s = norm(x);
        guard s > 0.0 else { return rho; }
        for &v in x { v /= s; }
        multiply(level_rs(l), level_col(l), level_val(l), x, y);
        for i in n { y[i] /= d[i]; }
        rho = dot(x, y);
        for i in n { x[i] = y[i]; }
    }
    return max(rho, norm(x));
}

// C = A B for CSR A (rows of A) and B (nb columns).
fn spmm(ars: i64[:], acol: i32[:], aval: f64[:], brs: i64[:], bcol: i32[:], bval: f64[:], nb: i64,
        crs: i64[>..]&, ccol: i32[>..]&, cval: f64[>..]&) {
    guard ars.len >= 1 && brs.len >= 1 && nb >= 0 else { abort("spmm: two CSR matrices"); }
    crs.clear();
    ccol.clear();
    cval.clear();
    crs.push(0);
    var at: i64[>..] = [];
    for c in nb { at.push(-1); }
    for r in ars.len - 1 {
        let row = ccol.len;
        for e in ars[r]..ars[r + 1] {
            let k = acol[e] as i64;
            for f in brs[k]..brs[k + 1] {
                let c = bcol[f] as i64;
                if at[c] < row {
                    at[c] = ccol.len;
                    ccol.push(c as i32);
                    cval.push(aval[e] * bval[f]);
                } else {
                    cval[at[c]] += aval[e] * bval[f];
                }
            }
        }
        crs.push(ccol.len);
    }
}

// The transpose of the CSR matrix (rs, col, val) with nc columns.
fn transpose(rs: i64[:], col: i32[:], val: f64[:], nc: i64, trs: i64[>..]&, tcol: i32[>..]&, tval: f64[>..]&) {
    guard rs.len >= 1 && col.len == val.len else { abort("transpose: a CSR matrix"); }
    trs.clear();
    for c in nc + 1 { trs.push(0); }
    for c in col { trs[c + 1]++; }
    for c in nc { trs[c + 1] += trs[c]; }
    tcol.clear();
    tval.clear();
    for e in col.len {
        tcol.push(0);
        tval.push(0.0);
    }
    var fill: i64[>..] = [];
    fill.append(trs[..nc]);
    for r in rs.len - 1 {
        for e in rs[r]..rs[r + 1] {
            let c = col[e] as i64;
            tcol[fill[c]] = r as i32;
            tval[fill[c]] = val[e];
            fill[c]++;
        }
    }
}

// The dense Cholesky factor of the last level into amg_chol. False if the
// level is not positive definite.
fn factor_last() -> bool {
    let l = amg_levels - 1;
    let n = amg_n[l];
    let rs = level_rs(l);
    let col = level_col(l);
    let val = level_val(l);
    amg_chol.clear();
    for i in n * n { amg_chol.push(0.0); }
    for r in n {
        for e in rs[r]..rs[r + 1] { amg_chol[n * r + col[e]] += val[e]; }
        // The empty row of a dependent column of P0 (tentative) holds its dof at 0.
        if rs[r] == rs[r + 1] { amg_chol[n * r + r] = 1.0; }
    }
    for j in n {
        var d = amg_chol[n * j + j];
        for k in j { d -= amg_chol[n * j + k] * amg_chol[n * j + k]; }
        guard d > 0.0 else { return false; }
        let s = sqrt(d);
        amg_chol[n * j + j] = s;
        for i in j + 1..n {
            var v = amg_chol[n * i + j];
            for k in j { v -= amg_chol[n * i + k] * amg_chol[n * j + k]; }
            amg_chol[n * i + j] = v / s;
        }
    }
    return true;
}

fn clear_levels() {
    amg_levels = 0;
    amg_n.clear();
    amg_rs.clear();
    amg_rs_at.clear();
    amg_col.clear();
    amg_val.clear();
    amg_e_at.clear();
    amg_prs.clear();
    amg_prs_at.clear();
    amg_pcol.clear();
    amg_pval.clear();
    amg_pe_at.clear();
    amg_x.clear();
    amg_b.clear();
    amg_res.clear();
    amg_diag.clear();
    amg_v_at.clear();
}

// Builds the levels for the matrix (row_start, col, val) with dofs in
// blocks of `bs`, and the constant of each component as the near null
// space. False if the last level is not positive definite.
fn amg_setup(row_start: i64[:], col: i32[:], val: f64[:], bs: i64) -> bool {
    var modes: f64[>..] = [];
    for r in row_start.len - 1 {
        for c in bs { modes.push(if r % bs == c { 1.0 } else { 0.0 }); }
    }
    return amg_setup_modes(row_start, col, val, bs, modes, bs);
}

// The same with the near null space `modes`: nb values for each dof, by
// rows.
fn amg_setup_modes(row_start: i64[:], col: i32[:], val: f64[:], bs: i64, modes: f64[:], nb: i64) -> bool {
    guard bs >= 1 && nb >= 1 && row_start.len >= 2 && (row_start.len - 1) % bs == 0 && col.len == val.len &&
          modes.len == nb * (row_start.len - 1) else { abort("amg_setup_modes: a CSR matrix with the dofs in blocks, and nb modes of each dof"); }
    clear_levels();
    push_level(row_start.len - 1, row_start, col, val);
    var ns: f64[>..] = [];
    ns.append(modes);
    var level_bs = bs;
    var p0rs: i64[>..] = [];
    var p0col: i32[>..] = [];
    var p0val: f64[>..] = [];
    var prs: i64[>..] = [];
    var pcol: i32[>..] = [];
    var pval: f64[>..] = [];
    var trs: i64[>..] = [];
    var tcol: i32[>..] = [];
    var tval: f64[>..] = [];
    var crs: i64[>..] = [];
    var ccol: i32[>..] = [];
    var cval: f64[>..] = [];
    for l in AMG_MAX_LEVELS - 1 {
        let n = amg_n[l];
        guard n > AMG_COARSE else { break; }
        let agg, count = aggregates(l, level_bs);
        let nc = count * nb;
        guard nc < n else { break; }
        let coarse = tentative(agg, count, level_bs, nb, ns, p0rs, p0col, p0val);
        ns.clear();
        ns.append(coarse);
        smooth(l, nc, p0rs, p0col, p0val, prs, pcol, pval);
        // The next level: Pᵀ (A P).
        spmm(level_rs(l), level_col(l), level_val(l), prs, pcol, pval, nc, p0rs, p0col, p0val);
        transpose(prs, pcol, pval, nc, trs, tcol, tval);
        spmm(trs, tcol, tval, p0rs, p0col, p0val, nc, crs, ccol, cval);
        if amg_pe_at.len == 0 { amg_pe_at.push(0); }
        amg_prs_at.push(amg_prs.len);
        amg_prs.append(prs);
        amg_pcol.append(pcol);
        amg_pval.append(pval);
        amg_pe_at.push(amg_pcol.len);
        push_level(nc, crs, ccol, cval);
        level_bs = nb;
    }
    return amg_n[amg_levels - 1] > 3000 || factor_last();
}

// P0 of the aggregates `agg` of the nodes of a level with blocks of bs
// dofs and the near null space ns (nb values for each dof), into p0. It
// gives the near null space of the coarse level.
fn tentative(agg: i64[:], count: i64, bs: i64, nb: i64, ns: f64[:], p0rs: i64[>..]&, p0col: i32[>..]&, p0val: f64[>..]&) -> f64[>..] {
    // The nodes of each aggregate: a counting sort.
    var start: i64[>..] = [];
    for a in count + 1 { start.push(0); }
    for a in agg { start[a + 1]++; }
    for a in count { start[a + 1] += start[a]; }
    var fill: i64[>..] = [];
    fill.append(start[..count]);
    var nodes: i64[>..] = [];
    for i in agg.len { nodes.push(0); }
    for at, i in agg {
        nodes[fill[at]] = i;
        fill[at]++;
    }
    // q[nb r + k]: column k of Q at dof r.
    var q: f64[>..] = [];
    q.append(ns);
    var coarse: f64[>..] = [];
    for i in count * nb * nb { coarse.push(0.0); }
    for a in count {
        for k in nb {
            var first = 0.0;
            for i in nodes[start[a]..start[a + 1]] {
                for c in bs { first += q[nb * (bs * i + c) + k] * q[nb * (bs * i + c) + k]; }
            }
            for j in k {
                var r = 0.0;
                for i in nodes[start[a]..start[a + 1]] {
                    for c in bs { r += q[nb * (bs * i + c) + j] * q[nb * (bs * i + c) + k]; }
                }
                coarse[nb * (nb * a + j) + k] = r;
                for i in nodes[start[a]..start[a + 1]] {
                    for c in bs { q[nb * (bs * i + c) + k] -= r * q[nb * (bs * i + c) + j]; }
                }
            }
            var norm = 0.0;
            for i in nodes[start[a]..start[a + 1]] {
                for c in bs { norm += q[nb * (bs * i + c) + k] * q[nb * (bs * i + c) + k]; }
            }
            let keep = norm > 1.0e-20 * first && norm > 0.0;
            coarse[nb * (nb * a + k) + k] = if keep { sqrt(norm) } else { 0.0 };
            for i in nodes[start[a]..start[a + 1]] {
                for c in bs { q[nb * (bs * i + c) + k] = if keep { q[nb * (bs * i + c) + k] / sqrt(norm) } else { 0.0 }; }
            }
        }
    }
    p0rs.clear();
    p0col.clear();
    p0val.clear();
    p0rs.push(0);
    for r in agg.len * bs {
        for k in nb {
            if q[nb * r + k] == 0.0 { continue; }
            p0col.push((nb * agg[r / bs] + k) as i32);
            p0val.push(q[nb * r + k]);
        }
        p0rs.push(p0col.len);
    }
    return coarse;
}

// P = P0 - ω D⁻¹ A P0 of level l, with ω = 4 / (3 ρ(D⁻¹ A)).
fn smooth(l: i64, nc: i64, p0rs: i64[:], p0col: i32[:], p0val: f64[:], prs: i64[>..]&, pcol: i32[>..]&, pval: f64[>..]&) {
    let n = amg_n[l];
    let omega = 4.0 / (3.0 * spectral_radius(l));
    spmm(level_rs(l), level_col(l), level_val(l), p0rs, p0col, p0val, nc, prs, pcol, pval);
    let d = amg_diag[amg_v_at[l]..amg_v_at[l] + n];
    for r in n {
        for e in prs[r]..prs[r + 1] { pval[e] *= -omega / d[r]; }
        // The columns of P0 are in the row of A P0: A has its diagonal.
        for f in p0rs[r]..p0rs[r + 1] {
            for e in prs[r]..prs[r + 1] {
                if pcol[e] == p0col[f] { pval[e] += p0val[f]; }
            }
        }
    }
}

// One Gauss-Seidel sweep on level l, forward or backward.
fn sweep(l: i64, forward: bool) {
    let n = amg_n[l];
    let rs = level_rs(l);
    let col = level_col(l);
    let val = level_val(l);
    let at = amg_v_at[l];
    for k in n {
        let r = if forward { k } else { n - 1 - k };
        var s = amg_b[at + r];
        var d = 0.0;
        for e in rs[r]..rs[r + 1] {
            let c = col[e] as i64;
            if c == r { d = val[e]; } else { s -= val[e] * amg_x[at + c]; }
        }
        amg_x[at + r] = if d > 0.0 { s / d } else { 0.0 };
    }
}

// The solve of the last level: the Cholesky factor, or sweeps if the level
// is too large to factor.
fn solve_last() {
    let l = amg_levels - 1;
    let n = amg_n[l];
    let at = amg_v_at[l];
    guard amg_chol.len == n * n else {
        for k in 10 {
            sweep(l, true);
            sweep(l, false);
        }
        return;
    }
    for i in n {
        var s = amg_b[at + i];
        for k in i { s -= amg_chol[n * i + k] * amg_x[at + k]; }
        amg_x[at + i] = s / amg_chol[n * i + i];
    }
    for k in n {
        let i = n - 1 - k;
        var s = amg_x[at + i];
        for j in i + 1..n { s -= amg_chol[n * j + i] * amg_x[at + j]; }
        amg_x[at + i] = s / amg_chol[n * i + i];
    }
}

// z = M⁻¹ r: one V-cycle from z = 0.
fn amg_apply(r: f64[:], z: f64[:]) {
    guard amg_levels > 0 && r.len == amg_n[0] && z.len == r.len else { abort("amg_apply: amg_setup, and r and z of the size of its matrix"); }
    for i in r.len {
        amg_b[i] = r[i];
        amg_x[i] = 0.0;
    }
    for l in amg_levels - 1 {
        let n = amg_n[l];
        let at = amg_v_at[l];
        let next = amg_v_at[l + 1];
        sweep(l, true);
        multiply(level_rs(l), level_col(l), level_val(l), amg_x[at..at + n], amg_res[at..at + n]);
        for i in n { amg_res[at + i] = amg_b[at + i] - amg_res[at + i]; }
        // b of level l + 1 = Pᵀ res.
        for i in amg_n[l + 1] {
            amg_b[next + i] = 0.0;
            amg_x[next + i] = 0.0;
        }
        let prs = amg_prs[amg_prs_at[l]..amg_prs_at[l] + n + 1];
        let pcol = amg_pcol[amg_pe_at[l]..amg_pe_at[l + 1]];
        let pval = amg_pval[amg_pe_at[l]..amg_pe_at[l + 1]];
        for i in n {
            for e in prs[i]..prs[i + 1] { amg_b[next + pcol[e]] += pval[e] * amg_res[at + i]; }
        }
    }
    solve_last();
    for k in amg_levels - 1 {
        let l = amg_levels - 2 - k;
        let n = amg_n[l];
        let at = amg_v_at[l];
        let next = amg_v_at[l + 1];
        let prs = amg_prs[amg_prs_at[l]..amg_prs_at[l] + n + 1];
        let pcol = amg_pcol[amg_pe_at[l]..amg_pe_at[l + 1]];
        let pval = amg_pval[amg_pe_at[l]..amg_pe_at[l + 1]];
        for i in n {
            var s = 0.0;
            for e in prs[i]..prs[i + 1] { s += pval[e] * amg_x[next + pcol[e]]; }
            amg_x[at + i] += s;
        }
        sweep(l, false);
    }
    for i in r.len { z[i] = amg_x[i]; }
}

// The levels of the last amg_setup without the fine level: the dofs of
// each level in n, the matrices of the coarse levels (rs, col, val), the
// prolongators (prs, pcol, pval) and the factor of the last level. Each
// row start counts from its own level. warm.goose keeps them.
fn amg_coarse(n: i64[>..]&, rs: i64[>..]&, col: i32[>..]&, val: f64[>..]&, prs: i64[>..]&, pcol: i32[>..]&,
              pval: f64[>..]&, chol: f64[>..]&) {
    guard amg_levels > 0 else { abort("amg_coarse: amg_setup"); }
    n.clear();
    n.append(amg_n);
    rs.clear();
    col.clear();
    val.clear();
    for l in 1..amg_levels {
        rs.append(level_rs(l));
        col.append(level_col(l));
        val.append(level_val(l));
    }
    prs.clear();
    prs.append(amg_prs);
    pcol.clear();
    pcol.append(amg_pcol);
    pval.clear();
    pval.append(amg_pval);
    chol.clear();
    chol.append(amg_chol);
}

// The levels of amg_coarse again, with the matrix (row_start, col, val) as
// the fine level. False if the arrays do not make the levels.
fn amg_restore(row_start: i64[:], col: i32[:], val: f64[:], n: i64[:], rs: i64[:], cl: i32[:], vl: f64[:],
               prs: i64[:], pcol: i32[:], pval: f64[:], chol: f64[:]) -> bool {
    guard n.len >= 1 && n.len <= AMG_MAX_LEVELS && n[0] == row_start.len - 1 && col.len == val.len else { return false; }
    // The sizes first, so that a file of another program fails here.
    var rows = 0;
    var entries = 0;
    var prows = 0;
    var pentries = 0;
    for l in n.len {
        guard n[l] >= 1 else { return false; }
        if l >= 1 {
            guard rows + n[l] + 1 <= rs.len && rs[rows] == 0 && rs[rows + n[l]] >= 0 else { return false; }
            entries += rs[rows + n[l]];
            rows += n[l] + 1;
        }
        if l + 1 < n.len {
            guard prows + n[l] + 1 <= prs.len && prs[prows] == 0 && prs[prows + n[l]] >= 0 else { return false; }
            pentries += prs[prows + n[l]];
            prows += n[l] + 1;
        }
    }
    let last = n[n.len - 1];
    guard rows == rs.len && entries == cl.len && cl.len == vl.len && prows == prs.len && pentries == pcol.len &&
          pcol.len == pval.len && (chol.len == 0 || chol.len == last * last) else { return false; }
    for l in n.len - 1 {
        guard n[l + 1] < n[l] else { return false; }
    }
    clear_levels();
    push_level(n[0], row_start, col, val);
    rows = 0;
    entries = 0;
    for l in 1..n.len {
        let r = rs[rows..rows + n[l] + 1];
        let count = r[n[l]];
        for e in entries..entries + count { guard cl[e] >= 0 && cl[e] < n[l] else { return false; } }
        push_level(n[l], r, cl[entries..entries + count], vl[entries..entries + count]);
        rows += n[l] + 1;
        entries += count;
    }
    amg_pe_at.push(0);
    prows = 0;
    pentries = 0;
    for l in n.len - 1 {
        let r = prs[prows..prows + n[l] + 1];
        let count = r[n[l]];
        for e in pentries..pentries + count { guard pcol[e] >= 0 && pcol[e] < n[l + 1] else { return false; } }
        amg_prs_at.push(amg_prs.len);
        amg_prs.append(r);
        amg_pcol.append(pcol[pentries..pentries + count]);
        amg_pval.append(pval[pentries..pentries + count]);
        amg_pe_at.push(amg_pcol.len);
        prows += n[l] + 1;
        pentries += count;
    }
    amg_chol.clear();
    amg_chol.append(chol);
    return true;
}
