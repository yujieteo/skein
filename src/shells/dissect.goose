// A direct solve for a model of one grid (mesh.goose) with less memory
// than the skyline: a multifrontal LDLᵀ factorization on a nested
// dissection of the grid. References: A. George (1973), "Nested dissection
// of a regular finite element mesh", SIAM J. Numer. Anal. 10, 345-363, and
// I. S. Duff and J. K. Reid (1983), "The multifrontal solution of
// indefinite sparse symmetric linear equations", ACM TOMS 9, 302-325.
//
// The skyline of a grid of order 2 has about 6 (2 m + 3) entries in each
// column, with m the cells on the shorter side. At m = 130 and 720,000 dofs
// that is about 18 GB. The dissection needs about O(n log n) entries.
//
// The dissection cuts the cells of the grid in two along the longer side,
// on a line of the lattice between two cells, again and again down to
// leaves of at most LEAF_CELLS × LEAF_CELLS cells. A cut line separates the
// two halves: no cell has nodes on both sides of it. Each part of the tree
// is a front. A front eliminates its nodes: the nodes of its cut line, or
// all nodes of a leaf. Its border is the nodes on the sides of its cells
// that are cut lines of the fronts above it. The front of a leaf adds the
// stiffness of its cells and of its beams. The front of a cut adds the
// update matrices (the Schur complements on the border) of its two parts.
// Then it eliminates its nodes, keeps L and D, and gives its update matrix
// to the front above it. The fronts run children first (the reverse of
// the preorder of the tree), so the update matrices make a stack.
//
// A fixed dof keeps a 1 on the diagonal and no other entry, as in the
// skyline. The fronts are dense, with all 6 dofs of each node.
import std;
import math;
namespace shells;
import .laminate;
import .model;
import .mitc4;
import .mitc9;
import .beam;
import .mesh;
import solvers from "../solvers";

let LEAF_CELLS = 2;
let PANEL = 48;         // the columns of a panel of the dense elimination
let MAX_DEPTH = 64;     // the depth of the tree: about 2 log2 of the cells

struct Front {
    i0: i64,            // the cells [i0, i1) × [j0, j1) of the grid
    i1: i64,
    j0: i64,
    j1: i64,
    nodes: i64,         // the nodes of the front at front_nodes[nodes..]: ne, then nb
    ne: i64,            // the nodes that the front eliminates
    nb: i64,            // the nodes of its border
    factor: i64,        // the first value of the front in lfactor
    leaf: bool,
}

var fronts: Front[>..] = [];
var front_nodes: i64[>..] = [];
var lfactor: f64[>..] = [];        // for each front: the rows of L and D of its nodes, then L of its border
var front_buf: f64[>..] = [];      // the dense front, n × n row-major: its lower triangle
var front_panel: f64[>..] = [];    // L times D of the columns of a panel: PANEL for each row of a front
var front_stack: f64[>..] = [];    // the update matrices that wait for their front
var front_map: i64[>..] = [];      // the place of each border dof of a child in its parent
var front_pos: i64[>..] = [];      // the place of each node in the present front, or -1
var kdiag: f64[>..] = [];          // the diagonal of K, for the singular check
var leaf_beam_start: i64[>..] = [];    // the beams of front t are leaf_beams[leaf_beam_start[t]..leaf_beam_start[t + 1]]
var leaf_beams: i64[>..] = [];
var front_entries = 0;             // the values of all arrays above, for the memory check

// The lattice point (a, b) of node n of grid g.
fn node_lattice(g: i64, n: i64) -> i64, i64 {
    guard g >= 0 && g < grids.len && n >= grids[g].first_node else { abort("node_lattice: 0 <= g < grids.len, and n is a node of grid g"); }
    let gr = grids[g];
    let k = n - gr.first_node;
    let la = gr.order * gr.nx + 1;
    let lb = gr.order * gr.nt + 1;
    if gr.x_fast { return k % la, k / la; }
    return k / lb, k % lb;
}

// The cell (i, j) of grid g that holds both nodes of beam b, or -1, -1.
fn beam_cell(g: i64, b: i64) -> i64, i64 {
    guard g >= 0 && g < grids.len && b >= 0 && b < beams.len else { abort("beam_cell: 0 <= g < grids.len, and 0 <= b < beams.len"); }
    let gr = grids[g];
    let o = gr.order;
    let a0, b0 = node_lattice(g, beams[b].n[0]);
    let a1, b1 = node_lattice(g, beams[b].n[1]);
    let i = min(min(a0, a1) / o, gr.nx - 1);
    let j = min(min(b0, b1) / o, gr.nt - 1);
    if max(a0, a1) > o * (i + 1) || max(b0, b1) > o * (j + 1) { return -1, -1; }
    return i, j;
}

// The grid of the model if the dissection can solve it, or -1: one grid
// that has each node and each shell, and each beam in a cell.
fn dissect_grid() -> i64 {
    if grids.len != 1 || shells.len + beams.len == 0 { return -1; }
    let gr = grids[0];
    if gr.first_node != 0 || node_count() != (gr.order * gr.nx + 1) * (gr.order * gr.nt + 1) { return -1; }
    var cells = 0;
    for k in gr.nx * gr.nt { if grid_cells[gr.cells + k] >= 0 { cells++; } }
    if cells != shells.len { return -1; }
    for b in beams.len {
        let i, j = beam_cell(0, b);
        if i < 0 { return -1; }
    }
    return 0;
}

// The cut of front f: the axis (0 along x, 1 along θ) and the cell line.
fn front_cut(f: Front) -> i64, i64 {
    if f.i1 - f.i0 >= f.j1 - f.j0 { return 0, (f.i0 + f.i1) / 2; }
    return 1, (f.j0 + f.j1) / 2;
}

// The tree of fronts of grid g in preorder, from a stack of the parts.
fn build_fronts(g: i64) {
    let gr = grids[g];
    fronts.clear();
    var parts: i64[256] = [0; 256];     // 4 values for each part on the stack
    parts[1] = gr.nx;
    parts[3] = gr.nt;
    var top = 1;
    while top > 0 {
        top--;
        var f = Front { i0: parts[4 * top], i1: parts[4 * top + 1], j0: parts[4 * top + 2], j1: parts[4 * top + 3],
                        nodes: 0, ne: 0, nb: 0, factor: 0, leaf: false };
        f.leaf = f.i1 - f.i0 <= LEAF_CELLS && f.j1 - f.j0 <= LEAF_CELLS;
        fronts.push(f);
        if f.leaf { continue; }
        guard top + 2 <= MAX_DEPTH else { abort("build_fronts: the tree is deeper than MAX_DEPTH"); }
        let axis, m = front_cut(f);
        let parts2: i64[8] = if axis == 0 { [f.i0, m, f.j0, f.j1, m, f.i1, f.j0, f.j1] }
                             else { [f.i0, f.i1, f.j0, m, f.i0, f.i1, m, f.j1] };
        for k in 8 { parts[4 * top + k] = parts2[k]; }
        top += 2;
    }
}

// True if lattice point (a, b) of front f is on a cut side of its cells:
// a side that is not on the edge of grid g.
fn on_cut(g: i64, f: Front, a: i64, b: i64) -> bool {
    let gr = grids[g];
    let o = gr.order;
    return (f.i0 > 0 && a == o * f.i0) || (f.i1 < gr.nx && a == o * f.i1)
        || (f.j0 > 0 && b == o * f.j0) || (f.j1 < gr.nt && b == o * f.j1);
}

// The nodes of front t: the nodes that it eliminates, then its border.
fn front_lists(g: i64, t: i64) {
    let f = fronts[t];
    let o = grids[g].order;
    let start = front_nodes.len;
    let axis, m = front_cut(f);
    for pass in 2 {
        for b in o * f.j0..o * f.j1 + 1 {
            for a in o * f.i0..o * f.i1 + 1 {
                let cut = on_cut(g, f, a, b);
                let mine = if f.leaf { !cut } else if axis == 0 { a == o * m && !cut } else { b == o * m && !cut };
                if (pass == 0 && mine) || (pass == 1 && cut) { front_nodes.push(grid_node(g, a, b)); }
            }
        }
        if pass == 0 { fronts[t].ne = front_nodes.len - start; }
    }
    fronts[t].nodes = start;
    fronts[t].nb = front_nodes.len - start - fronts[t].ne;
}

// The beams of each leaf, from the leaf of each cell.
fn leaf_beam_lists(g: i64) {
    let gr = grids[g];
    var cell_leaf: i64[>..] = [];
    for k in gr.nx * gr.nt { cell_leaf.push(-1); }
    for t in fronts.len {
        let f = fronts[t];
        if !f.leaf { continue; }
        for j in f.j0..f.j1 { for i in f.i0..f.i1 { cell_leaf[i + gr.nx * j] = t; } }
    }
    var leaf: i64[>..] = [];
    leaf_beam_start.clear();
    for t in fronts.len + 1 { leaf_beam_start.push(0); }
    for b in beams.len {
        let i, j = beam_cell(g, b);
        leaf.push(cell_leaf[i + gr.nx * j]);
        leaf_beam_start[leaf[b] + 1]++;
    }
    for t in fronts.len { leaf_beam_start[t + 1] += leaf_beam_start[t]; }
    var next: i64[>..] = [];
    for t in fronts.len { next.push(leaf_beam_start[t]); }
    leaf_beams.clear();
    for b in beams.len { leaf_beams.push(0); }
    for b in beams.len {
        leaf_beams[next[leaf[b]]] = b;
        next[leaf[b]]++;
    }
}

// Sizes a list to n values of v.
fn sized_f64(v: f64[>..]&, n: i64) {
    v.clear();
    for k in n { v.push(0.0); }
}

// Makes the tree, the lists and the arrays of the dissection of grid g.
// The memory check runs after the tree and before the large arrays.
fn dissect_prepare(g: i64) {
    guard g >= 0 && g < grids.len && fixed.len == dof_count() else { abort("dissect_prepare: 0 <= g < grids.len, and a fixed value for each dof"); }
    build_fronts(g);
    front_nodes.clear();
    for t in fronts.len { front_lists(g, t); }
    leaf_beam_lists(g);
    var values = 0;
    var buf = 0;
    var rows = 0;
    var border = 0;
    for t in fronts.len {
        let ne6 = 6 * fronts[t].ne;
        let nb6 = 6 * fronts[t].nb;
        fronts[t].factor = values;
        values += ne6 * (ne6 + 1) / 2 + nb6 * ne6;
        buf = max(buf, (ne6 + nb6) * (ne6 + nb6));
        rows = max(rows, ne6 + nb6);
        border = max(border, nb6);
    }
    let stack = stack_peak();
    front_entries = values + buf + PANEL * rows + stack;
    solvers::memory_check("shell skyline ldlt", front_entries, dof_count(), 1, 0, "");
    sized_f64(lfactor, values);
    sized_f64(front_buf, buf);
    sized_f64(front_panel, PANEL * rows);
    sized_f64(front_stack, stack);
    sized_f64(kdiag, dof_count());
    front_map.clear();
    for k in border { front_map.push(0); }
    front_pos.clear();
    for n in node_count() { front_pos.push(-1); }
}

// The largest size of the stack of update matrices, in the run order.
fn stack_peak() -> i64 {
    var sizes: i64[64] = [0; 64];
    var depth = 0;
    var now = 0;
    var peak = 0;
    for tt in fronts.len {
        let f = fronts[fronts.len - 1 - tt];
        if !f.leaf {
            guard depth >= 2 else { abort("stack_peak: two children for each cut"); }
            depth -= 2;
            now -= sizes[depth] + sizes[depth + 1];
        }
        guard depth < MAX_DEPTH else { abort("stack_peak: the depth < MAX_DEPTH"); }
        sizes[depth] = 36 * f.nb * f.nb;
        depth++;
        now += 36 * f.nb * f.nb;
        peak = max(peak, now);
    }
    return peak;
}

// The global dof of place k of front f.
fn front_dof(f: Front, k: i64) -> i64 { return 6 * front_nodes[f.nodes + k / 6] + k % 6; }

// The place of global dof d in the present front.
fn local_dof(d: i64) -> i64 { return 6 * front_pos[d / 6] + d % 6; }

// Adds the m × m cell matrix k on the dofs d to the front of size n, and
// its diagonal to kdiag. A fixed dof gets no entry.
fn add_cell(n: i64, d: i64[54], m: i64, k: f64[:]) {
    guard (m == 12 || m == 24 || m == 54) && k.len >= m * m else { abort("add_cell: m is 12, 24 or 54, and k has m x m entries"); }
    for a in m {
        let ga = d[a];
        kdiag[ga] += k[m * a + a];
        if fixed[ga] != 0 { continue; }
        let la = local_dof(ga);
        for b in m {
            let gb = d[b];
            if fixed[gb] != 0 { continue; }
            let lb = local_dof(gb);
            if lb <= la { front_buf[n * la + lb] += k[m * a + b]; }
        }
    }
}

// Adds the cells and the beams of leaf front t of grid g to the front.
fn assemble_leaf(g: i64, t: i64, n: i64) {
    let f = fronts[t];
    var kbeam: f64[144] = [0.0; 144];
    for j in f.j0..f.j1 {
        for i in f.i0..f.i1 {
            let c = grid_cell(g, i, j);
            if c < 0 { continue; }
            cell_stiffness(c);
            let d, m = cell_dofs(c, false);
            add_cell(n, d, m, if m == 54 { ke9[..] } else { ke[..] });
        }
    }
    for k in leaf_beam_start[t]..leaf_beam_start[t + 1] {
        beam_stiffness(leaf_beams[k], kbeam);
        let d, m = cell_dofs(leaf_beams[k], true);
        add_cell(n, d, m, kbeam);
    }
}

// Adds the update matrix of child front c, at front_stack[at..], to the
// front of size n.
fn extend_add(c: i64, at: i64, n: i64) {
    let f = fronts[c];
    let nb6 = 6 * f.nb;
    for r in nb6 { front_map[r] = local_dof(front_dof(f, 6 * f.ne + r)); }
    for r in nb6 {
        let lr = front_map[r];
        for s in r + 1 {
            let ls = front_map[s];
            let v = front_stack[at + nb6 * r + s];
            if ls <= lr { front_buf[n * lr + ls] += v; } else { front_buf[n * ls + lr] += v; }
        }
    }
}

// Factors the columns k0..k1 of the front of size n (a panel): the
// columns before k0 are done and their update is in the front. Sets
// front_panel[PANEL i + m - k0] = L[i][m] D[m]. False at a pivot with no
// stiffness: singular_dof is then its dof.
fn factor_panel(f: Front, n: i64, k0: i64, k1: i64) -> bool {
    guard k0 >= 0 && k1 > k0 && k1 - k0 <= PANEL && k1 <= 6 * f.ne && n * n <= front_buf.len else { abort("factor_panel: 0 <= k0 < k1 <= 6 ne, k1 - k0 <= PANEL, and the front fits"); }
    for k in k0..k1 {
        let rk = n * k;
        let pk = PANEL * k - k0;
        var d = front_buf[rk + k];
        for m in k0..k { d -= front_buf[rk + m] * front_panel[pk + m]; }
        let g = front_dof(f, k);
        let scale = if fixed[g] != 0 { 1.0 } else { abs(kdiag[g]) };
        if !(d > 1.0e-12 * scale) {
            singular_dof = g;
            return false;
        }
        front_buf[rk + k] = d;
        for i in k + 1..n {
            let ri = n * i;
            var s = front_buf[ri + k];
            for m in k0..k { s -= front_buf[ri + m] * front_panel[pk + m]; }
            let l = s / d;
            front_buf[ri + k] = l;
            front_panel[PANEL * i + k - k0] = l * d;
        }
    }
    return true;
}

// Takes the panel k0..k1 from the rest of the front of size n:
// F[i][j] -= Σ L[i][m] D[m] L[j][m] for k1 <= j <= i, 4 columns j at once.
fn update_rest(n: i64, k0: i64, k1: i64) {
    guard k0 >= 0 && k1 > k0 && k1 - k0 <= PANEL && k1 <= n && n * n <= front_buf.len else { abort("update_rest: 0 <= k0 < k1 <= n, and k1 - k0 <= PANEL"); }
    for i in k1..n {
        let ri = n * i;
        var j = k1;
        while j + 4 <= i + 1 {
            let p0 = PANEL * j - k0;
            var s0 = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;
            var s3 = 0.0;
            for m in k0..k1 {
                let l = front_buf[ri + m];
                s0 += l * front_panel[p0 + m];
                s1 += l * front_panel[p0 + PANEL + m];
                s2 += l * front_panel[p0 + 2 * PANEL + m];
                s3 += l * front_panel[p0 + 3 * PANEL + m];
            }
            front_buf[ri + j] -= s0;
            front_buf[ri + j + 1] -= s1;
            front_buf[ri + j + 2] -= s2;
            front_buf[ri + j + 3] -= s3;
            j += 4;
        }
        while j <= i {
            var s = 0.0;
            for m in k0..k1 { s += front_buf[ri + m] * front_panel[PANEL * j - k0 + m]; }
            front_buf[ri + j] -= s;
            j++;
        }
    }
}

// Eliminates the nodes of front f (size n), by panels. The rows of the
// border then hold L of the border, and the update matrix after them.
fn eliminate(f: Front, n: i64) -> bool {
    let ne6 = 6 * f.ne;
    for k in ne6 { if fixed[front_dof(f, k)] != 0 { front_buf[n * k + k] = 1.0; } }
    var k0 = 0;
    while k0 < ne6 {
        let k1 = min(k0 + PANEL, ne6);
        if !factor_panel(f, n, k0, k1) { return false; }
        update_rest(n, k0, k1);
        k0 = k1;
    }
    return true;
}

// Keeps L and D of front f (size n) in lfactor, and puts its update matrix
// at front_stack[at..].
fn keep_front(f: Front, n: i64, at: i64) {
    let ne6 = 6 * f.ne;
    let nb6 = 6 * f.nb;
    var q = f.factor;
    for k in ne6 {
        for m in k + 1 { lfactor[q + m] = front_buf[n * k + m]; }
        q += k + 1;
    }
    for r in nb6 {
        let rr = n * (ne6 + r);
        for m in ne6 { lfactor[q + m] = front_buf[rr + m]; }
        q += ne6;
        for s in r + 1 { front_stack[at + nb6 * r + s] = front_buf[rr + ne6 + s]; }
    }
}

// The multifrontal factorization of grid g. False if a pivot has no
// stiffness (singular_dof).
fn dissect_factor(g: i64) -> bool {
    guard front_pos.len == node_count() && kdiag.len == dof_count() else { abort("dissect_factor: dissect_prepare made the arrays"); }
    singular_dof = -1;
    for &v in kdiag { v = 0.0; }
    var waiting: i64[64] = [0; 64];     // the fronts of the update matrices on the stack
    var place: i64[64] = [0; 64];       // and their place in front_stack
    var depth = 0;
    var top = 0;
    for tt in fronts.len {
        let t = fronts.len - 1 - tt;
        let f = fronts[t];
        let n = 6 * (f.ne + f.nb);
        for k in n * n { front_buf[k] = 0.0; }
        for k in f.ne + f.nb { front_pos[front_nodes[f.nodes + k]] = k; }
        if f.leaf { assemble_leaf(g, t, n); }
        else {
            for c in 2 {
                depth--;
                extend_add(waiting[depth], place[depth], n);
                top = place[depth];
            }
        }
        if !eliminate(f, n) { return false; }
        keep_front(f, n, top);
        waiting[depth] = t;
        place[depth] = top;
        depth++;
        top += 36 * f.nb * f.nb;
        for k in f.ne + f.nb { front_pos[front_nodes[f.nodes + k]] = -1; }
    }
    return true;
}

// Solves L D Lᵀ disp = load with the factor of dissect_factor.
fn dissect_substitute() {
    guard fixed.len == dof_count() && load.len == fixed.len && disp.len == fixed.len else { abort("dissect_substitute: the loads and the solution have each dof"); }
    for j in dof_count() { disp[j] = if fixed[j] != 0 { 0.0 } else { load[j] }; }
    for tt in fronts.len { forward_front(fronts[fronts.len - 1 - tt]); }
    for t in fronts.len { backward_front(fronts[t]); }
}

// The forward step and the diagonal step of front f: L y = b, then y / D.
fn forward_front(f: Front) {
    let ne6 = 6 * f.ne;
    var q = f.factor;
    for k in ne6 {
        var s = disp[front_dof(f, k)];
        for m in k { s -= lfactor[q + m] * disp[front_dof(f, m)]; }
        disp[front_dof(f, k)] = s;
        q += k + 1;
    }
    for r in 6 * f.nb {
        var s = 0.0;
        for m in ne6 { s += lfactor[q + m] * disp[front_dof(f, m)]; }
        disp[front_dof(f, ne6 + r)] -= s;
        q += ne6;
    }
    q = f.factor;
    for k in ne6 {
        disp[front_dof(f, k)] /= lfactor[q + k];
        q += k + 1;
    }
}

// The backward step of front f: Lᵀ x = y, with x of the border known.
fn backward_front(f: Front) {
    let ne6 = 6 * f.ne;
    let tri = ne6 * (ne6 + 1) / 2;
    for r in 6 * f.nb {
        let x = disp[front_dof(f, ne6 + r)];
        let q = f.factor + tri + ne6 * r;
        for m in ne6 { disp[front_dof(f, m)] -= lfactor[q + m] * x; }
    }
    for kk in ne6 {
        let k = ne6 - 1 - kk;
        let q = f.factor + k * (k + 1) / 2;
        let x = disp[front_dof(f, k)];
        for m in k { disp[front_dof(f, m)] -= lfactor[q + m] * x; }
    }
}
