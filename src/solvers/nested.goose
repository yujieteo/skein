// A direct solve of K x = b for a symmetric positive definite K on the
// nodes of a lattice: a multifrontal LDLᵀ factorization on a nested
// dissection of the lattice, as src/shells/dissect.goose does for the 2D
// grid of a shell. References: A. George (1973), "Nested dissection of a
// regular finite element mesh", SIAM J. Numer. Anal. 10, 345-363, and
// I. S. Duff and J. K. Reid (1983), "The multifrontal solution of
// indefinite sparse symmetric linear equations", ACM TOMS 9, 302-325.
//
// The caller gives the lattice: its points along each axis, and the point
// i + lat[0] (j + lat[1] k) of each node. A point need not have a node, so
// the grid of a union of boxes fits. K is CSR (solvers.goose) with `dpn`
// dofs for each node: node n has the dofs dpn n to dpn n + dpn - 1. Any
// physics, cell type and order fits: a plane of the lattice can cut a part
// if no entry of K couples a node on one side of it to a node on the other
// side, and the pattern of K gives these planes.
//
// The dissection cuts each part of the lattice in two on the plane that
// is nearest to its middle, normal to the axis with the smallest
// cross-section, down to leaves of at most LEAF_POINTS points. Each part
// is a front. A front eliminates its nodes: the nodes of its cut plane, or
// all nodes of a leaf. Its border is the nodes on the sides of the part
// that are cut planes of the fronts above it. A front adds the rows of K
// of its nodes and the update matrices (the Schur complements on the
// border) of its two parts. Then it eliminates its nodes, keeps L and D,
// and gives its update matrix to the front above it. The fronts run
// children first (the reverse of the preorder of the tree), so the update
// matrices make a stack.
//
// A fixed dof keeps a 1 on the diagonal and no other entry. A load change
// with the same K needs only nested_solve: two triangular solves.
import std;
import math;
namespace solvers;

let LEAF_POINTS = 64;        // Reason: 32 and 128 give the same time on interposer_warpage (14.1 s)
let NESTED_WORK = 1.0e12;   // the multiplies and adds above which a solve takes CG with AMG
let ND_PANEL = 48;          // the columns of a panel of the dense elimination
let ND_MAX_DEPTH = 64;      // the depth of the tree: about log2 of the points

struct NdFront {
    lo: i64[3],         // the lattice points lo to hi (both in) along each axis
    hi: i64[3],
    nodes: i64,         // the nodes of the front at nd_nodes[nodes..]: ne, then nb
    ne: i64,            // the nodes that the front eliminates
    nb: i64,            // the nodes of its border
    factor: i64,        // the first value of the front in nd_factor
    leaf: bool,
}

var nd_lat: i64[3] = [0, 0, 0];
var nd_dpn = 0;
var nd_at: i64[>..] = [];       // the node at each lattice point, or -1
var nd_place: i64[>..] = [];    // the lattice point of each node
var nd_cuts: bool[>..] = [];    // whether plane p of axis a can cut: at p, lat[0] + p or lat[0] + lat[1] + p
var nd_fronts: NdFront[>..] = [];
var nd_nodes: i64[>..] = [];
var nd_factor: f64[>..] = [];   // for each front: the rows of L and D of its nodes, then L of its border
var nd_buf: f64[>..] = [];      // the dense front, n × n row-major: its lower triangle
var nd_panel: f64[>..] = [];    // L times D of the columns of a panel: ND_PANEL for each row of a front
var nd_stack: f64[>..] = [];    // the update matrices that wait for their front
var nd_map: i64[>..] = [];      // the place of each border dof of a child in its parent
var nd_pos: i64[>..] = [];      // the place of each node in the present front, or -1
var nd_diag: f64[>..] = [];     // the diagonal of K, for the singular check
var nd_values = 0;              // the values of nd_factor
var nd_entries = 0;             // the values of the arrays of the factorization, for the memory check
var nd_flops = 0.0;             // the multiplies and adds of the factorization
var nd_singular = -1;           // the dof of the last pivot with no stiffness, or -1

fn nested_entries() -> i64 { return nd_entries; }
fn nested_flops() -> f64 { return nd_flops; }
fn nested_singular() -> i64 { return nd_singular; }

// Whether a solve takes the direct path on the plan of nested_plan: its
// work is at most NESTED_WORK, and it fits in the memory with `vectors`
// vectors of the dofs more. Else the solve takes CG with AMG, which needs
// less memory and, on a large cube, less work.
fn nested_pays(vectors: i64) -> bool {
    guard nd_flops <= NESTED_WORK else { return false; }
    return memory_fits("lattice ldlt", nd_entries, nd_dpn * nd_place.len, vectors, 0);
}

// The lattice point of node n along each axis.
fn nd_point(n: i64) -> i64[3] {
    guard n >= 0 && n < nd_place.len else { abort("nd_point: a node of the plan"); }
    var p = nd_place[n];
    var q: i64[3] = [0, 0, 0];
    for a in 3 {
        q[a] = p % nd_lat[a];
        p /= nd_lat[a];
    }
    return q;
}

fn nd_cut_at(a: i64, p: i64) -> i64 { return if a == 0 { p } else if a == 1 { nd_lat[0] + p } else { nd_lat[0] + nd_lat[1] + p }; }

// The plan of the dissection of the pattern (row_start, col) with `dpn`
// dofs for each node on the lattice `lat`, node n at point place[n]: the
// tree, the nodes of each front and the sizes. It allocates no value of
// the factorization: nested_entries and nested_flops give its size and
// its work before nested_factor. False if a node is not on the lattice or
// two nodes are on one point.
fn nested_plan(lat: i64[3], place: i64[:], dpn: i64, row_start: i64[:], col: i32[:]) -> bool {
    guard dpn > 0 && lat[0] > 0 && lat[1] > 0 && lat[2] > 0 && row_start.len == dpn * place.len + 1 else { abort("nested_plan: dpn > 0, a lattice, and the rows of dpn dofs of each node"); }
    nd_lat = lat;
    nd_dpn = dpn;
    nd_at.clear();
    for p in lat[0] * lat[1] * lat[2] { nd_at.push(-1); }
    nd_place.clear();
    for p, n in place {
        guard p >= 0 && p < nd_at.len && nd_at[p] < 0 else { return false; }
        nd_at[p] = n;
        nd_place.push(p);
    }
    plan_cuts(row_start, col);
    plan_tree();
    nd_nodes.clear();
    for t in nd_fronts.len { front_list(t); }
    plan_sizes();
    return true;
}

// The planes that can cut: no entry of K couples nodes on both sides.
fn plan_cuts(row_start: i64[:], col: i32[:]) {
    var crossed: i64[>..] = [];   // the couplings across each plane, as differences
    for k in nd_lat[0] + nd_lat[1] + nd_lat[2] + 1 { crossed.push(0); }
    for n in nd_place.len {
        let p = nd_point(n);
        for e in row_start[nd_dpn * n]..row_start[nd_dpn * n + 1] {
            let q = nd_point(col[e] / nd_dpn);
            for a in 3 {
                let lo = min(p[a], q[a]);
                let hi = max(p[a], q[a]);
                if hi - lo < 2 { continue; }
                crossed[nd_cut_at(a, lo + 1)]++;
                crossed[nd_cut_at(a, hi)]--;
            }
        }
    }
    nd_cuts.clear();
    var now = 0;
    for k in crossed.len - 1 {
        if k == nd_lat[0] || k == nd_lat[0] + nd_lat[1] { now = 0; }
        now += crossed[k];
        nd_cuts.push(now == 0);
    }
}

// The cut of a part lo to hi: the axis and the plane, or -1, 0 if no plane
// can cut it.
fn part_cut(lo: i64[3], hi: i64[3]) -> i64, i64 {
    var axis = -1;
    var plane = 0;
    var cross = 0;
    for a in 3 {
        let mid = (lo[a] + hi[a]) / 2;
        var p = -1;
        for s in hi[a] - lo[a] {
            let near = if s % 2 == 0 { mid - s / 2 } else { mid + (s + 1) / 2 };
            if near > lo[a] && near < hi[a] && nd_cuts[nd_cut_at(a, near)] { p = near; break; }
        }
        if p < 0 { continue; }
        var c = 1;
        for b in 3 { if b != a { c *= hi[b] - lo[b] + 1; } }
        if axis < 0 || c < cross {
            axis = a;
            plane = p;
            cross = c;
        }
    }
    return axis, plane;
}

fn part_points(lo: i64[3], hi: i64[3]) -> i64 { return (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1); }

// The tree of fronts in preorder, from a stack of the parts.
fn plan_tree() {
    nd_fronts.clear();
    var parts: i64[>..] = [0, 0, 0, nd_lat[0] - 1, nd_lat[1] - 1, nd_lat[2] - 1];   // lo and hi of each part on the stack
    while parts.len > 0 {
        let top = parts.len - 6;
        var f = NdFront { lo: [parts[top], parts[top + 1], parts[top + 2]], hi: [parts[top + 3], parts[top + 4], parts[top + 5]],
                          nodes: 0, ne: 0, nb: 0, factor: 0, leaf: true };
        parts.resize(top);
        let axis, m = part_cut(f.lo, f.hi);
        f.leaf = axis < 0 || part_points(f.lo, f.hi) <= LEAF_POINTS;
        nd_fronts.push(f);
        if f.leaf { continue; }
        guard parts.len < 6 * ND_MAX_DEPTH else { abort("plan_tree: the tree is deeper than ND_MAX_DEPTH"); }
        var hi0 = f.hi;
        hi0[axis] = m;
        var lo1 = f.lo;
        lo1[axis] = m;
        parts.append(f.lo);
        parts.append(hi0);
        parts.append(lo1);
        parts.append(f.hi);
    }
}

// True if lattice point q of front f is on a cut side of the part: a side
// that is not on the edge of the lattice.
fn on_cut(f: NdFront, q: i64[3]) -> bool {
    for a in 3 {
        if (f.lo[a] > 0 && q[a] == f.lo[a]) || (f.hi[a] < nd_lat[a] - 1 && q[a] == f.hi[a]) { return true; }
    }
    return false;
}

// The nodes of front t: the nodes that it eliminates, then its border.
fn front_list(t: i64) {
    let f = nd_fronts[t];
    let start = nd_nodes.len;
    let axis, m = part_cut(f.lo, f.hi);
    for pass in 2 {
        for k in f.lo[2]..f.hi[2] + 1 {
            for j in f.lo[1]..f.hi[1] + 1 {
                for i in f.lo[0]..f.hi[0] + 1 {
                    let n = nd_at[i + nd_lat[0] * (j + nd_lat[1] * k)];
                    if n < 0 { continue; }
                    let q: i64[3] = [i, j, k];
                    let cut = on_cut(f, q);
                    let mine = !cut && (f.leaf || q[axis] == m);
                    if (pass == 0 && mine) || (pass == 1 && cut) { nd_nodes.push(n); }
                }
            }
        }
        if pass == 0 { nd_fronts[t].ne = nd_nodes.len - start; }
    }
    nd_fronts[t].nodes = start;
    nd_fronts[t].nb = nd_nodes.len - start - nd_fronts[t].ne;
}

// The place of each front in nd_factor, the sizes of the arrays and the
// work of the factorization.
fn plan_sizes() {
    var buf = 0;
    var rows = 0;
    nd_values = 0;
    nd_flops = 0.0;
    for t in nd_fronts.len {
        let e = nd_dpn * nd_fronts[t].ne;
        let b = nd_dpn * nd_fronts[t].nb;
        nd_fronts[t].factor = nd_values;
        nd_values += e * (e + 1) / 2 + b * e;
        buf = max(buf, (e + b) * (e + b));
        rows = max(rows, e + b);
        let ef = e * 1.0;
        nd_flops += ef * ef * ef / 3.0 + ef * ef * b + ef * b * b;
    }
    nd_entries = nd_values + buf + ND_PANEL * rows + stack_peak() + nd_place.len * nd_dpn;
}

// The largest size of the stack of update matrices, in the run order.
fn stack_peak() -> i64 {
    var sizes: i64[ND_MAX_DEPTH] = [0; ND_MAX_DEPTH];
    var depth = 0;
    var now = 0;
    var peak = 0;
    for tt in nd_fronts.len {
        let f = nd_fronts[nd_fronts.len - 1 - tt];
        let b = nd_dpn * f.nb;
        if !f.leaf {
            guard depth >= 2 else { abort("stack_peak: two children for each cut"); }
            depth -= 2;
            now -= sizes[depth] + sizes[depth + 1];
        }
        guard depth < ND_MAX_DEPTH else { abort("stack_peak: the depth < ND_MAX_DEPTH"); }
        sizes[depth] = b * b;
        depth++;
        now += b * b;
        peak = max(peak, now);
    }
    return peak;
}

// The global dof of place k of front f.
fn nd_dof(f: NdFront, k: i64) -> i64 { return nd_dpn * nd_nodes[f.nodes + k / nd_dpn] + k % nd_dpn; }

// The place of global dof d in the present front, or -1.
fn nd_local(d: i64) -> i64 {
    let p = nd_pos[d / nd_dpn];
    return if p < 0 { -1 } else { nd_dpn * p + d % nd_dpn };
}

// Sizes a list to n values of v.
fn nd_sized(v: f64[>..]&, n: i64) {
    v.clear();
    for k in n { v.push(0.0); }
}

// Frees the arrays of the factorization (nested_alloc), for a solve that
// takes another path after a nested_factor that returned false.
fn nested_free() {
    nd_factor = [];
    nd_buf = [];
    nd_panel = [];
    nd_stack = [];
    nd_diag = [];
    nd_map = [];
    nd_pos = [];
}

// The arrays of the factorization of the plan, once for each plan.
fn nested_alloc() {
    guard nd_fronts.len > 0 else { abort("nested_alloc: nested_plan"); }
    var buf = 0;
    var rows = 0;
    var border = 0;
    for f in nd_fronts {
        let n = nd_dpn * (f.ne + f.nb);
        buf = max(buf, n * n);
        rows = max(rows, n);
        border = max(border, nd_dpn * f.nb);
    }
    if nd_factor.len != nd_values { nd_sized(nd_factor, nd_values); }
    nd_sized(nd_buf, buf);
    nd_sized(nd_panel, ND_PANEL * rows);
    nd_sized(nd_stack, stack_peak());
    nd_sized(nd_diag, nd_dpn * nd_place.len);
    nd_map.clear();
    for k in border { nd_map.push(0); }
    nd_pos.clear();
    for n in nd_place.len { nd_pos.push(-1); }
}

// Adds the rows of K of the nodes that front f (size n) eliminates: the
// entries of the columns of its nodes, in the lower triangle. A column of
// a front below is in the update matrices already.
fn assemble_front(f: NdFront, n: i64, row_start: i64[:], col: i32[:], val: f64[:], fixed: u8[:]) {
    for k in nd_dpn * f.ne {
        let r = nd_dof(f, k);
        if fixed[r] != 0 { continue; }
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            let lc = nd_local(c);
            if lc < 0 || fixed[c] != 0 { continue; }
            if lc <= k { nd_buf[n * k + lc] += val[e]; }
            else if lc >= nd_dpn * f.ne { nd_buf[n * lc + k] += val[e]; }
        }
    }
}

// Adds the update matrix of child front c, at nd_stack[at..], to the front
// of size n.
fn extend_add(c: i64, at: i64, n: i64) {
    let f = nd_fronts[c];
    let e = nd_dpn * f.ne;
    let b = nd_dpn * f.nb;
    for r in b { nd_map[r] = nd_local(nd_dof(f, e + r)); }
    for r in b {
        let lr = nd_map[r];
        for s in r + 1 {
            let ls = nd_map[s];
            let v = nd_stack[at + b * r + s];
            if ls <= lr { nd_buf[n * lr + ls] += v; } else { nd_buf[n * ls + lr] += v; }
        }
    }
}

// Factors the columns k0..k1 of the front f of size n (a panel): the
// columns before k0 are done and their update is in the front. Sets
// nd_panel[ND_PANEL i + m - k0] = L[i][m] D[m]. False at a pivot with no
// stiffness: nd_singular is then its dof.
fn factor_panel(f: NdFront, n: i64, k0: i64, k1: i64, fixed: u8[:]) -> bool {
    guard k0 >= 0 && k1 > k0 && k1 - k0 <= ND_PANEL && k1 <= nd_dpn * f.ne && n * n <= nd_buf.len else { abort("factor_panel: 0 <= k0 < k1 <= the dofs of the front, k1 - k0 <= ND_PANEL, and the front fits"); }
    for k in k0..k1 {
        let rk = n * k;
        let pk = ND_PANEL * k - k0;
        var d = nd_buf[rk + k];
        for m in k0..k { d -= nd_buf[rk + m] * nd_panel[pk + m]; }
        let g = nd_dof(f, k);
        let scale = if fixed[g] != 0 { 1.0 } else { abs(nd_diag[g]) };
        if !(d > 1.0e-12 * scale) {
            nd_singular = g;
            return false;
        }
        nd_buf[rk + k] = d;
        for i in k + 1..n {
            let ri = n * i;
            var s = nd_buf[ri + k];
            for m in k0..k { s -= nd_buf[ri + m] * nd_panel[pk + m]; }
            let l = s / d;
            nd_buf[ri + k] = l;
            nd_panel[ND_PANEL * i + k - k0] = l * d;
        }
    }
    return true;
}

// Takes the panel k0..k1 from the rest of the front of size n:
// F[i][j] -= Σ L[i][m] D[m] L[j][m] for k1 <= j <= i, on tiles of 4 rows
// and 4 columns (update_tile), and the last rows one entry at a time.
// Reason: in the page (emcc, 1 thread) a tile does 4.0 Gflop/s and a row
// with 4 columns 2.1 Gflop/s; native, 15 and 12.
fn update_rest(n: i64, k0: i64, k1: i64) {
    guard k0 >= 0 && k1 > k0 && k1 - k0 <= ND_PANEL && k1 <= n && n * n <= nd_buf.len else { abort("update_rest: 0 <= k0 < k1 <= n, and k1 - k0 <= ND_PANEL"); }
    var c: f64[16] = [0.0; 16];
    var i = k1;
    while i + 4 <= n {
        var j = k1;
        while j <= i {
            update_tile(l_row(n, i, k0, k1), l_row(n, i + 1, k0, k1), l_row(n, i + 2, k0, k1), l_row(n, i + 3, k0, k1),
                        ld_row(j, k1 - k0), ld_row(j + 1, k1 - k0), ld_row(j + 2, k1 - k0), ld_row(j + 3, k1 - k0), c);
            for a in 4 {
                for q in 4 { if j + q <= i + a { nd_buf[n * (i + a) + j + q] -= c[4 * a + q]; } }
            }
            j += 4;
        }
        i += 4;
    }
    while i < n {
        let l = l_row(n, i, k0, k1);
        for j in k1..i + 1 {
            let ld = ld_row(j, k1 - k0);
            assert(ld.len == l.len);
            var s = 0.0;
            for m in l.len { s += l[m] * ld[m]; }
            nd_buf[n * i + j] -= s;
        }
        i++;
    }
}

// L of row i of the panel k0..k1 of the front of size n.
fn l_row(n: i64, i: i64, k0: i64, k1: i64) -> f64[:] { return nd_buf[n * i + k0..n * i + k1]; }

// L D of row j of the panel of w columns.
fn ld_row(j: i64, w: i64) -> f64[:] { return nd_panel[ND_PANEL * j..ND_PANEL * j + w]; }

// c[4 a + q] = Σ a_a[m] b_q[m]: 4 rows of L times 4 rows of L D. Each
// index is m, so no bounds check is left in the loop.
checkfree fn update_tile(a0: f64[:], a1: f64[:], a2: f64[:], a3: f64[:], b0: f64[:], b1: f64[:], b2: f64[:], b3: f64[:], c: f64[16]&) {
    assert(a1.len == a0.len && a2.len == a0.len && a3.len == a0.len && b0.len == a0.len && b1.len == a0.len && b2.len == a0.len && b3.len == a0.len);
    var s00 = 0.0;
    var s01 = 0.0;
    var s02 = 0.0;
    var s03 = 0.0;
    var s10 = 0.0;
    var s11 = 0.0;
    var s12 = 0.0;
    var s13 = 0.0;
    var s20 = 0.0;
    var s21 = 0.0;
    var s22 = 0.0;
    var s23 = 0.0;
    var s30 = 0.0;
    var s31 = 0.0;
    var s32 = 0.0;
    var s33 = 0.0;
    for m in a0.len {
        let x0 = a0[m];
        let x1 = a1[m];
        let x2 = a2[m];
        let x3 = a3[m];
        let y0 = b0[m];
        let y1 = b1[m];
        let y2 = b2[m];
        let y3 = b3[m];
        s00 += x0 * y0;
        s01 += x0 * y1;
        s02 += x0 * y2;
        s03 += x0 * y3;
        s10 += x1 * y0;
        s11 += x1 * y1;
        s12 += x1 * y2;
        s13 += x1 * y3;
        s20 += x2 * y0;
        s21 += x2 * y1;
        s22 += x2 * y2;
        s23 += x2 * y3;
        s30 += x3 * y0;
        s31 += x3 * y1;
        s32 += x3 * y2;
        s33 += x3 * y3;
    }
    c = [s00, s01, s02, s03, s10, s11, s12, s13, s20, s21, s22, s23, s30, s31, s32, s33];
}

// Eliminates the nodes of front f (size n), by panels. The rows of the
// border then hold L of the border, and the update matrix after them.
fn eliminate(f: NdFront, n: i64, fixed: u8[:]) -> bool {
    let e = nd_dpn * f.ne;
    for k in e { if fixed[nd_dof(f, k)] != 0 { nd_buf[n * k + k] = 1.0; } }
    var k0 = 0;
    while k0 < e {
        let k1 = min(k0 + ND_PANEL, e);
        if !factor_panel(f, n, k0, k1, fixed) { return false; }
        update_rest(n, k0, k1);
        k0 = k1;
    }
    return true;
}

// Keeps L and D of front f (size n) in nd_factor, and puts its update
// matrix at nd_stack[at..].
fn keep_front(f: NdFront, n: i64, at: i64) {
    let e = nd_dpn * f.ne;
    let b = nd_dpn * f.nb;
    var q = f.factor;
    for k in e {
        for m in k + 1 { nd_factor[q + m] = nd_buf[n * k + m]; }
        q += k + 1;
    }
    for r in b {
        let rr = n * (e + r);
        for m in e { nd_factor[q + m] = nd_buf[rr + m]; }
        q += e;
        for s in r + 1 { nd_stack[at + b * r + s] = nd_buf[rr + e + s]; }
    }
}

// The multifrontal factorization of K = (row_start, col, val) on the plan
// of nested_plan, with the fixed dofs as rows and columns of the identity.
// False if a pivot has no stiffness (nested_singular).
fn nested_factor(row_start: i64[:], col: i32[:], val: f64[:], fixed: u8[:]) -> bool {
    guard fixed.len == nd_dpn * nd_place.len && row_start.len == fixed.len + 1 else { abort("nested_factor: nested_plan of K, and a fixed value for each dof"); }
    nested_alloc();
    for d in fixed.len { nd_diag[d] = val[entry(row_start, col, d, d)]; }
    nd_singular = -1;
    var waiting: i64[ND_MAX_DEPTH] = [0; ND_MAX_DEPTH];   // the fronts of the update matrices on the stack
    var place: i64[ND_MAX_DEPTH] = [0; ND_MAX_DEPTH];     // and their place in nd_stack
    var depth = 0;
    var top = 0;
    for tt in nd_fronts.len {
        let t = nd_fronts.len - 1 - tt;
        let f = nd_fronts[t];
        let n = nd_dpn * (f.ne + f.nb);
        for k in n * n { nd_buf[k] = 0.0; }
        for k in f.ne + f.nb { nd_pos[nd_nodes[f.nodes + k]] = k; }
        if !f.leaf {
            for c in 2 {
                depth--;
                extend_add(waiting[depth], place[depth], n);
                top = place[depth];
            }
        }
        assemble_front(f, n, row_start, col, val, fixed);
        if !eliminate(f, n, fixed) { return false; }
        keep_front(f, n, top);
        waiting[depth] = t;
        place[depth] = top;
        depth++;
        top += nd_dpn * nd_dpn * f.nb * f.nb;
        for k in f.ne + f.nb { nd_pos[nd_nodes[f.nodes + k]] = -1; }
    }
    return true;
}

// Solves L D Lᵀ x = b with the factor of nested_factor: x is 0 on the fixed
// dofs.
fn nested_solve(b: f64[:], x: f64[:], fixed: u8[:]) { substitute(nd_factor, b, x, fixed); }

// Solves L D Lᵀ x = b with `factor`, the values of nd_factor of a factor of
// the plan. Each front copies its dofs of x to y, so that the loops run on
// contiguous rows (nd_dot, nd_take). Reason: on interposer_warpage the
// loops on the dofs of x took 2.0 s, and the loops on y take 0.2 s.
fn substitute(factor: f64[:], b: f64[:], x: f64[:], fixed: u8[:]) {
    guard b.len == x.len && fixed.len == x.len && x.len == nd_dpn * nd_place.len && factor.len == nd_values else { abort("substitute: a factor of the plan, and b and x of each dof"); }
    var y: f64[>..] = [];
    for f in nd_fronts { for k in nd_dpn * (f.ne + f.nb) - y.len { y.push(0.0); } }
    for j in x.len { x[j] = if fixed[j] != 0 { 0.0 } else { b[j] }; }
    for tt in nd_fronts.len { forward_front(factor, nd_fronts[nd_fronts.len - 1 - tt], x, y); }
    for f in nd_fronts { backward_front(factor, f, x, y); }
}

// The forward step and the diagonal step of front f: L y = b, then y / D.
fn forward_front(factor: f64[:], f: NdFront, x: f64[:], y: f64[:]) {
    let e = nd_dpn * f.ne;
    let n = e + nd_dpn * f.nb;
    for k in n { y[k] = x[nd_dof(f, k)]; }
    var q = f.factor;
    for k in e {
        y[k] -= nd_dot(factor[q..q + k], y[0..k]);
        q += k + 1;
    }
    for r in e..n {
        y[r] -= nd_dot(factor[q..q + e], y[0..e]);
        q += e;
    }
    q = f.factor;
    for k in e {
        y[k] /= factor[q + k];
        q += k + 1;
    }
    for k in n { x[nd_dof(f, k)] = y[k]; }
}

// The backward step of front f: Lᵀ x = y, with x of the border known.
fn backward_front(factor: f64[:], f: NdFront, x: f64[:], y: f64[:]) {
    let e = nd_dpn * f.ne;
    let n = e + nd_dpn * f.nb;
    for k in n { y[k] = x[nd_dof(f, k)]; }
    var q = f.factor + e * (e + 1) / 2;
    for r in e..n {
        nd_take(factor[q..q + e], y[r], y[0..e]);
        q += e;
    }
    for kk in e {
        let k = e - 1 - kk;
        let at = f.factor + k * (k + 1) / 2;
        nd_take(factor[at..at + k], y[k], y[0..k]);
    }
    for k in e { x[nd_dof(f, k)] = y[k]; }
}

// Σ a[m] b[m].
checkfree fn nd_dot(a: f64[:], b: f64[:]) -> f64 {
    assert(b.len == a.len);
    var s = 0.0;
    for m in a.len { s += a[m] * b[m]; }
    return s;
}

// y -= v a.
checkfree fn nd_take(a: f64[:], v: f64, y: f64[:]) {
    assert(y.len == a.len);
    for m in a.len { y[m] -= v * a[m]; }
}
