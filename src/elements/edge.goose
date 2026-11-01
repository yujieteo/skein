// The edge elements of `family: nedelec` (spec §13): the Nédélec element
// of the first kind and lowest order on P1 triangles and tetrahedra, for
//
//     curl(ν curl(u)) + c u = s
//
// with u a `Complex` vector (or a real one, with no imaginary part), ν > 0
// and c complex, constant on the cells of a material. The weak form is
// ∫ ν curl u · curl v + c u · v dV + ∮ (n × h) · v dA = ∫ s · v dV, with
// h = ν curl u. A natural condition `cross(n, h) = z u` adds ∮ z u_t · v_t
// over its facets (u_t is the part of u along the facet).
//
// The dof of an edge is the line integral of u along it, from its lower
// vertex to the higher one. The basis function of edge (i, j) of a cell is
// w = λi ∇λj - λj ∇λi, with curl w = 2 ∇λi × ∇λj, times the sign of the
// local edge. In 2D, curl is the z part of that vector. A complex dof is two
// dofs, as in wave.goose: component 0 is the real part and 1 the imaginary
// part, and a complex matrix entry is a 2 × 2 block of the real matrix.
//
// The essential condition `cross(n, u) = cross(n, g)` holds each edge of the
// facets of the region at ∫ g · dl (spec §12.1, gap G4): on a curved wall
// only the part of u along the wall is held. The tree gauge (`gauge: tree`,
// gap G27) holds the edges of a spanning tree of the cells with c = 0 at 0,
// so that the gradients there do not make the matrix singular.
import std;
import math;
import dictionary;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .solve;
import .wave;
import .post;
import .eigen;
import solvers from "../solvers";

let EDGE_TOLERANCE = 1.0e-10;   // the GMRES residual relative to the right side
let EDGE_ITERATIONS = 20000;

var edge_v: i64[>..] = [];      // the two vertices of each edge, the lower index first
var cell_edge: i64[>..] = [];   // the edge of each local edge of each cell
var cell_sign: f64[>..] = [];   // +1 if the local edge runs from the lower vertex, else -1
var edge_nu: f64[>..] = [];     // ν of each cell
var edge_c: f64[>..] = [];      // c of each cell, as re and im
var gauge_count = 0;            // the edges of the tree of the last gauge
var facet_open: u8[>..] = [];   // 1 for each facet of each cell on the boundary of the mesh
var node_held: u8[>..] = [];    // 1 for each node of a held edge

fn edge_total() -> i64 { return edge_v.len / 2; }

fn local_edges() -> i64 { return edge_count(dim); }

// Starts the edge physics on the P1 mesh, with ν and c on each cell. In 2D,
// t is the thickness.
fn set_edge(nu: f64, c: f64[2], t: f64) {
    guard order == 1 && per_cell == dim + 1 && nu > 0.0 && t > 0.0 else { abort("set_edge: a P1 mesh, nu > 0 and t > 0"); }
    physics = EDGE;
    ncomp = 2;
    thickness = if dim == 3 { 1.0 } else { t };
    clear_stages();
    build_edges();
    other_dofs = 2 * edge_total();
    fixed.clear();
    prescribed.clear();
    load.clear();
    sol.clear();
    for d in other_dofs {
        fixed.push(0);
        prescribed.push(0.0);
        load.push(0.0);
        sol.push(0.0);
    }
    clear_matrix();
    edge_nu.clear();
    edge_c.clear();
    for e in cell_count() {
        edge_nu.push(nu);
        edge_c.append(c);
    }
    gauge_count = 0;
}

// The edges of the mesh, each once, and the edge and the sign of each
// local edge of each cell.
fn build_edges() {
    edge_v.clear();
    cell_edge.clear();
    cell_sign.clear();
    let nv = node_count();
    var at = dictionary<i64, i64> {};
    for c in cell_count() {
        for e in local_edges() {
            let a = cell_node(c, edge_node(dim, e, 0));
            let b = cell_node(c, edge_node(dim, e, 1));
            let k = at.get_or_insert(min(a, b) * nv + max(a, b), edge_total());
            if k == edge_total() {
                edge_v.push(min(a, b));
                edge_v.push(max(a, b));
            }
            cell_edge.push(k);
            cell_sign.push(if a < b { 1.0 } else { -1.0 });
        }
    }
    find_open_facets();
}

// Gives the cells in the shape g the coefficients ν and c. It gives the
// number of these cells.
fn edge_material(g: i64, nu: f64, c: f64[2]) -> i64 {
    guard physics == EDGE && g >= 0 && nu > 0.0 else { abort("edge_material: an edge physics, a shape and nu > 0"); }
    var count = 0;
    for e in cell_count() {
        if !cell_in(e, g) { continue; }
        edge_nu[e] = nu;
        edge_c[2 * e] = c[0];
        edge_c[2 * e + 1] = c[1];
        count++;
    }
    return count;
}

// The basis of cell c at the barycentric coordinates l: w[3 e + t] of each
// local edge e with its sign, and curl w[3 e + t] (in 2D, only t = 2). It
// gives the determinant of the map of the cell.
fn edge_basis(c: i64, l: f64[4], w: f64[18]&, cw: f64[18]&) -> f64 {
    let det = cell_geometry(c, [0.0; 3]);
    for e in local_edges() {
        let i = edge_node(dim, e, 0);
        let j = edge_node(dim, e, 1);
        let s = cell_sign[local_edges() * c + e];
        for t in 3 { w[3 * e + t] = s * (l[i] * grad[3 * j + t] - l[j] * grad[3 * i + t]); }
        let gi: f64[3] = [grad[3 * i], grad[3 * i + 1], grad[3 * i + 2]];
        let gj: f64[3] = [grad[3 * j], grad[3 * j + 1], grad[3 * j + 2]];
        cw[3 * e] = 2.0 * s * (gi[1] * gj[2] - gi[2] * gj[1]);
        cw[3 * e + 1] = 2.0 * s * (gi[2] * gj[0] - gi[0] * gj[2]);
        cw[3 * e + 2] = 2.0 * s * (gi[0] * gj[1] - gi[1] * gj[0]);
    }
    return det;
}

// The complex matrix me of cell c, edge by edge: the curl term and the
// term of c, by the rule of degree 2.
fn edge_cell(c: i64, me: f64[:]) {
    for &v in me { v = 0.0; }
    let ne = local_edges();
    var w: f64[18] = [0.0; 18];
    var cw: f64[18] = [0.0; 18];
    for q in quad_count(dim) {
        let xi, wt = quad_point(dim, q);
        var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
        for r in dim {
            l[r + 1] = xi[r];
            l[0] -= xi[r];
        }
        let dv = wt * edge_basis(c, l, w, cw) * thickness;
        for a in ne {
            for b in ne {
                var k = 0.0;
                var m = 0.0;
                for t in 3 {
                    k += cw[3 * a + t] * cw[3 * b + t];
                    m += w[3 * a + t] * w[3 * b + t];
                }
                me[2 * (ne * a + b)] += (edge_nu[c] * k + edge_c[2 * c] * m) * dv;
                me[2 * (ne * a + b) + 1] += edge_c[2 * c + 1] * m * dv;
            }
        }
    }
}

// The pattern of the edges: two edges are neighbours if they share a cell.
fn edge_pattern() {
    clear_matrix();
    let ne = edge_total();
    let le = local_edges();
    var keys: i64[>..] = [];
    for c in cell_count() {
        for a in le {
            for b in le { keys.push(cell_edge[le * c + a] * ne + cell_edge[le * c + b]); }
        }
    }
    sort(keys);
    row_start.push(0);
    var first = 0;
    for a in ne {
        var last = first;
        while last < keys.len && keys[last] / ne == a { last++; }
        for p in 2 {
            var prev = -1;
            for e in first..last {
                if keys[e] == prev { continue; }
                prev = keys[e];
                for q in 2 {
                    col.push((2 * (keys[e] % ne) + q) as i32);
                    val.push(0.0);
                }
            }
            row_start.push(col.len);
        }
        first = last;
    }
}

// The matrix of the cells in val, after the pattern. The source and the
// natural conditions are in `load` already.
fn edge_assemble() {
    edge_pattern();
    let m = wave_restart(dof_total());
    solvers::memory_check("simplex wave gmres ilu0", col.len, dof_total(), m + 9 + (m + 1) * m / max(dof_total(), 1), 0, "");
    assembled_source.clear();
    for d in dof_total() { assembled_source.push(0.0); }
    let le = local_edges();
    var me: f64[72] = [0.0; 72];
    for c in cell_count() {
        edge_cell(c, me[..2 * le * le]);
        for a in le {
            for b in le {
                let k = le * a + b;
                add_block(cell_edge[le * c + a], cell_edge[le * c + b], [me[2 * k], me[2 * k + 1]]);
            }
        }
    }
    add_impedance();
}

// --- the source and the conditions ---------------------------------------------

// Adds ∫ s · v dV to the load, with s = S(p) at each point p: the real
// parts of the components, then the imaginary parts.
fn edge_source<S>() {
    guard physics == EDGE && dofs_sized() else { abort("edge_source: an edge physics and the arrays of the dofs"); }
    let le = local_edges();
    var w: f64[18] = [0.0; 18];
    var cw: f64[18] = [0.0; 18];
    for c in cell_count() {
        for q in rich_count() {
            let xi, wt = rich_point(q);
            var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
            for r in dim {
                l[r + 1] = xi[r];
                l[0] -= xi[r];
            }
            let dv = wt * edge_basis(c, l, w, cw) * thickness;
            var p: f64[3] = [0.0; 3];
            for v in dim + 1 {
                for t in 3 { p[t] += l[v] * node_x(cell_node(c, v), t); }
            }
            let s: f64[6] = S(p);
            for a in le {
                let d = 2 * cell_edge[le * c + a];
                for t in 3 {
                    load[d] += s[t] * w[3 * a + t] * dv;
                    load[d + 1] += s[3 + t] * w[3 * a + t] * dv;
                }
            }
        }
    }
}

// Whether facet f of cell c has its vertices on one part of the region.
fn facet_on_region(c: i64, f: i64) -> bool {
    for k in region.len {
        var all = true;
        for l in dim { if !on_part(cell_node(c, facet_node(dim, f, l)), k) { all = false; } }
        if all { return true; }
    }
    return false;
}

// Whether facet f of cell c is on the boundary of the mesh and on the region:
// a facet inside the mesh can have its vertices on a curved wall, but its
// edges are chords through the cells.
fn edge_facet_on(c: i64, f: i64) -> bool {
    return facet_open[facet_count(dim) * c + f] != 0 && facet_on_region(c, f);
}

// The nodes of the held edges in node_held.
fn find_held_nodes() {
    node_held.clear();
    for v in node_count() { node_held.push(0); }
    for e in edge_total() {
        if fixed[2 * e] != 0 {
            node_held[edge_v[2 * e]] = 1;
            node_held[edge_v[2 * e + 1]] = 1;
        }
    }
}

// The key of facet f of cell c: its edge in 2D, else its two edges of the
// lowest indices, which give the third.
fn facet_key(c: i64, f: i64) -> i64 {
    let le = local_edges();
    var lo = edge_total();
    var mid = edge_total();
    for e in le {
        if !facet_has_edge(f, e) { continue; }
        let k = cell_edge[le * c + e];
        if k < lo {
            mid = lo;
            lo = k;
        } else if k < mid {
            mid = k;
        }
    }
    return if dim == 2 { lo } else { lo * edge_total() + mid };
}

// facet_open: 1 for each facet of each cell that no other cell has.
fn find_open_facets() {
    var keys: i64[>..] = [];
    for c in cell_count() {
        for f in facet_count(dim) { keys.push(facet_key(c, f)); }
    }
    var sorted: i64[>..] = [];
    sorted.append(keys);
    sort(sorted);
    facet_open.clear();
    for k in keys {
        var lo = 0;
        var hi = sorted.len - 1;
        while lo < hi {
            let m = (lo + hi) / 2;
            if sorted[m] < k { lo = m + 1; } else { hi = m; }
        }
        let shared = (lo > 0 && sorted[lo - 1] == k) || (lo + 1 < sorted.len && sorted[lo + 1] == k);
        facet_open.push(if shared { 0 } else { 1 });
    }
}

// The local edges of facet f of a cell: the edges between its vertices.
fn facet_has_edge(f: i64, e: i64) -> bool {
    var count = 0;
    for l in dim {
        let v = facet_node(dim, f, l);
        if v == edge_node(dim, e, 0) || v == edge_node(dim, e, 1) { count++; }
    }
    return count == 2;
}

// Holds each edge of the facets of the region at ∫ g · dl, with g = G(p)
// as in edge_source: the condition cross(n, u) = cross(n, g). It gives the
// number of edges that it holds.
fn fix_tangent<G>() -> i64 {
    guard physics == EDGE && dofs_sized() else { abort("fix_tangent: an edge physics and the arrays of the dofs"); }
    var count = 0;
    let le = local_edges();
    for c in cell_count() {
        for f in facet_count(dim) {
            if !edge_facet_on(c, f) { continue; }
            for e in le {
                let k = cell_edge[le * c + e];
                if !facet_has_edge(f, e) || fixed[2 * k] != 0 { continue; }
                let a = edge_v[2 * k];
                let b = edge_v[2 * k + 1];
                var v: f64[2] = [0.0, 0.0];
                for g in 2 {
                    let s = if g == 0 { LINE_G } else { 1.0 - LINE_G };
                    var p: f64[3] = [0.0; 3];
                    for t in 3 { p[t] = (1.0 - s) * node_x(a, t) + s * node_x(b, t); }
                    let gv: f64[6] = G(p);
                    for t in 3 {
                        v[0] += 0.5 * gv[t] * (node_x(b, t) - node_x(a, t));
                        v[1] += 0.5 * gv[3 + t] * (node_x(b, t) - node_x(a, t));
                    }
                }
                for p in 2 {
                    fixed[2 * k + p] = 1;
                    prescribed[2 * k + p] = v[p];
                }
                count++;
            }
        }
    }
    return count;
}

// The facets of the impedance conditions of the step: cell, facet and z.
var imp_facets: i64[>..] = [];
var imp_z: f64[>..] = [];

// Adds the facets of the region to the condition cross(n, h) = z u.
fn edge_impedance(z: f64[2]) -> i64 {
    guard physics == EDGE else { abort("edge_impedance: an edge physics"); }
    var count = 0;
    for c in cell_count() {
        for f in facet_count(dim) {
            if !edge_facet_on(c, f) { continue; }
            imp_facets.push(c);
            imp_facets.push(f);
            imp_z.append(z);
            count++;
        }
    }
    return count;
}

// The unit normal and the size of facet f of cell c (in 2D, the length
// times the thickness).
fn facet_normal(c: i64, f: i64) -> f64[3], f64 {
    let a = cell_node(c, facet_node(dim, f, 0));
    let b = cell_node(c, facet_node(dim, f, 1));
    var t1: f64[3] = [0.0; 3];
    for t in 3 { t1[t] = node_x(b, t) - node_x(a, t); }
    var n: f64[3] = [t1[1], -t1[0], 0.0];
    var size = sqrt(t1[0] * t1[0] + t1[1] * t1[1]) * thickness;
    if dim == 3 {
        let e = cell_node(c, facet_node(dim, f, 2));
        var t2: f64[3] = [0.0; 3];
        for t in 3 { t2[t] = node_x(e, t) - node_x(a, t); }
        n = [t1[1] * t2[2] - t1[2] * t2[1], t1[2] * t2[0] - t1[0] * t2[2], t1[0] * t2[1] - t1[1] * t2[0]];
        size = 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    }
    let len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for t in 3 { n[t] /= len; }
    return n, size;
}

// The barycentric coordinates in its cell of point q of the rule of degree
// 2 on facet f, and the weight of the point on the unit facet.
fn facet_point(f: i64, q: i64) -> f64[4], f64 {
    let xi, w = quad_point(dim - 1, q);
    var fl: f64[3] = [1.0 - xi[0] - xi[1], xi[0], xi[1]];
    if dim == 2 { fl = [1.0 - xi[0], xi[0], 0.0]; }
    var l: f64[4] = [0.0; 4];
    for v in dim { l[facet_node(dim, f, v)] = fl[v]; }
    return l, w * (if dim == 2 { 1.0 } else { 2.0 });
}

// Adds ∮ z u_t · v_t over the facets of the impedance conditions.
fn add_impedance() {
    let le = local_edges();
    var w: f64[18] = [0.0; 18];
    var cw: f64[18] = [0.0; 18];
    for i in imp_facets.len / 2 {
        let c = imp_facets[2 * i];
        let f = imp_facets[2 * i + 1];
        let n, size = facet_normal(c, f);
        for q in quad_count(dim - 1) {
            let l, wt = facet_point(f, q);
            edge_basis(c, l, w, cw);
            for a in le {
                for b in le {
                    var dn_a = 0.0;
                    var dn_b = 0.0;
                    var m = 0.0;
                    for t in 3 {
                        dn_a += w[3 * a + t] * n[t];
                        dn_b += w[3 * b + t] * n[t];
                        m += w[3 * a + t] * w[3 * b + t];
                    }
                    let tt = (m - dn_a * dn_b) * wt * size;
                    add_block(cell_edge[le * c + a], cell_edge[le * c + b], [imp_z[2 * i] * tt, imp_z[2 * i + 1] * tt]);
                }
            }
        }
    }
}

fn clear_impedance() {
    imp_facets.clear();
    imp_z.clear();
}

// --- the tree gauge --------------------------------------------------------------

// The root of vertex v in the forest `up`, with the path halved.
fn tree_root(up: i64[:], v: i64) -> i64 {
    var r = v;
    var steps = 0;
    while up[r] != r && steps < up.len {
        up[r] = up[up[r]];
        r = up[r];
        steps++;
    }
    return r;
}

// The tree gauge: holds at 0 the edges of a spanning tree of the vertices
// of the cells with c = 0. The vertices of a cell with c ≠ 0 and of a held
// edge are one root first, so a gradient that the tree removes is 0 there.
// Call it after the essential conditions. It gives the edges of the tree.
fn tree_gauge() -> i64 {
    guard physics == EDGE && dofs_sized() else { abort("tree_gauge: an edge physics and the arrays of the dofs"); }
    var up: i64[>..] = [];
    for v in node_count() { up.push(v); }
    var root = -1;
    var free_edge: u8[>..] = [];
    for k in edge_total() { free_edge.push(0); }
    let le = local_edges();
    for c in cell_count() {
        let conductor = edge_c[2 * c] != 0.0 || edge_c[2 * c + 1] != 0.0;
        for e in le {
            let k = cell_edge[le * c + e];
            if !conductor && fixed[2 * k] == 0 && free_edge[k] == 0 { free_edge[k] = 1; }
            if conductor { free_edge[k] = 2; }
            if !conductor && fixed[2 * k] == 0 { continue; }
            for s in 2 {
                let v = tree_root(up, edge_v[2 * k + s]);
                if root < 0 { root = v; }
                up[v] = tree_root(up, root);
            }
        }
    }
    var count = 0;
    for k in edge_total() {
        if free_edge[k] != 1 { continue; }
        let a = tree_root(up, edge_v[2 * k]);
        let b = tree_root(up, edge_v[2 * k + 1]);
        if a == b { continue; }
        up[a] = b;
        for p in 2 {
            fixed[2 * k + p] = 1;
            prescribed[2 * k + p] = 0.0;
        }
        count++;
    }
    gauge_count = count;
    return count;
}

// --- the solve -------------------------------------------------------------------

// The right side with the fixed values, and the fixed rows the identity.
fn edge_rhs() -> f64[>..] {
    let n = sol.len;
    var b: f64[>..] = [];
    for d in n {
        sol[d] = if fixed[d] != 0 { prescribed[d] } else { 0.0 };
        b.push(if fixed[d] != 0 { prescribed[d] } else { load[d] });
    }
    for r in n {
        for e in row_start[r]..row_start[r + 1] {
            let c = col[e] as i64;
            if fixed[r] != 0 {
                val[e] = if c == r { 1.0 } else { 0.0 };
            } else if fixed[c] != 0 {
                b[r] -= val[e] * sol[c];
                val[e] = 0.0;
            }
        }
    }
    return b;
}

// --- the gradient correction ----------------------------------------------------
//
// The gradient ∇φ of each P1 nodal function φ is in the space of the
// edges, and on the gradients the matrix has only the term of c. When
// |c| h² is small to ν (0.02 for the coax line of gap G4), these are many
// small eigenvalues that ILU(0) does not hold, and GMRES stalls. When c is
// not 0 on any cell, the preconditioner adds a correction in the space of
// the gradients (Hiptmair 1998): after ILU(0), z += G S⁻¹ Gᵀ (r - A z).
// G gives each node the gradient of its φ on the edges, and on them A is
// close to the phase of c times S = ∫ |c| ∇φi · ∇φj dV, the matrix of
// grad_rs, grad_col and grad_val. One V-cycle of AMG solves S, once for the
// real parts and once for the imaginary parts. A node on a held edge is out
// of the correction: its row of S is the identity.

var grad_rs: i64[>..] = [];
var grad_col: i32[>..] = [];
var grad_val: f64[>..] = [];
var grad_phase: f64[2] = [1.0, 0.0];   // the phase of the sum of c over the cells
var grad_t: f64[>..] = [];             // r - A z
var grad_y: f64[>..] = [];             // Gᵀ (r - A z): the real parts of the nodes, then the imaginary parts
var grad_w: f64[>..] = [];             // S⁻¹ of grad_y

// The pattern of S: each node and its neighbours along the edges.
fn gradient_pattern() {
    let nv = node_count();
    var keys: i64[>..] = [];
    for v in nv { keys.push(v * nv + v); }
    for e in edge_total() {
        keys.push(edge_v[2 * e] * nv + edge_v[2 * e + 1]);
        keys.push(edge_v[2 * e + 1] * nv + edge_v[2 * e]);
    }
    sort(keys);
    grad_rs.clear();
    grad_col.clear();
    grad_val.clear();
    grad_rs.push(0);
    var k = 0;
    for v in nv {
        while k < keys.len && keys[k] / nv == v {
            grad_col.push((keys[k] % nv) as i32);
            grad_val.push(0.0);
            k++;
        }
        grad_rs.push(grad_col.len);
    }
}

// The values of S, the identity in the rows of the held nodes, and the
// phase of c. False if c is 0 on a cell.
fn gradient_values() -> bool {
    var vol = 0.0;
    for q in quad_count(dim) {
        let xi, wt = quad_point(dim, q);
        vol += wt + 0.0 * xi[0];
    }
    var sum: f64[2] = [0.0, 0.0];
    for c in cell_count() {
        let a = sqrt(edge_c[2 * c] * edge_c[2 * c] + edge_c[2 * c + 1] * edge_c[2 * c + 1]);
        guard a > 0.0 else { return false; }
        sum[0] += edge_c[2 * c] / a;
        sum[1] += edge_c[2 * c + 1] / a;
        let s = a * vol * cell_geometry(c, [0.0; 3]) * thickness;
        for i in dim + 1 {
            for j in dim + 1 {
                var g = 0.0;
                for t in 3 { g += grad[3 * i + t] * grad[3 * j + t]; }
                grad_val[solvers::entry(grad_rs, grad_col, cell_node(c, i), cell_node(c, j))] += s * g;
            }
        }
    }
    let m = sqrt(sum[0] * sum[0] + sum[1] * sum[1]);
    guard m > 0.0 else { return false; }
    grad_phase = [sum[0] / m, sum[1] / m];
    for r in node_count() {
        for e in grad_rs[r]..grad_rs[r + 1] {
            let c = grad_col[e] as i64;
            if node_held[r] != 0 || node_held[c] != 0 { grad_val[e] = if r == c { 1.0 } else { 0.0 }; }
        }
    }
    return true;
}

// Builds S and its AMG levels. False if there is no correction: c is 0 on
// a cell, or S is not positive definite.
fn gradient_setup() -> bool {
    find_held_nodes();
    gradient_pattern();
    guard gradient_values() && solvers::amg_setup(grad_rs, grad_col, grad_val, 1) else { return false; }
    grad_t.clear();
    for d in sol.len { grad_t.push(0.0); }
    grad_y.clear();
    grad_w.clear();
    for v in 2 * node_count() {
        grad_y.push(0.0);
        grad_w.push(0.0);
    }
    return true;
}

// z += G S⁻¹ Gᵀ (r - A z), with S divided by the phase of c.
fn gradient_apply(r: f64[:], z: f64[:]) {
    let nv = node_count();
    solvers::rows_multiply(row_start, col, val, z, grad_t);
    for d in r.len { grad_t[d] = r[d] - grad_t[d]; }
    for &y in grad_y { y = 0.0; }
    for e in edge_total() {
        for q in 2 {
            grad_y[q * nv + edge_v[2 * e + 1]] += grad_t[2 * e + q];
            grad_y[q * nv + edge_v[2 * e]] -= grad_t[2 * e + q];
        }
    }
    for v in nv {
        if node_held[v] != 0 {
            grad_y[v] = 0.0;
            grad_y[nv + v] = 0.0;
        }
    }
    solvers::amg_apply(grad_y[..nv], grad_w[..nv]);
    solvers::amg_apply(grad_y[nv..], grad_w[nv..]);
    let p = grad_phase;
    for e in edge_total() {
        if fixed[2 * e] != 0 { continue; }
        let a = edge_v[2 * e];
        let b = edge_v[2 * e + 1];
        let wr = grad_w[b] - grad_w[a];
        let wi = grad_w[nv + b] - grad_w[nv + a];
        z[2 * e] += wr * p[0] + wi * p[1];
        z[2 * e + 1] += wi * p[0] - wr * p[1];
    }
}

// Assembles and solves by GMRES with ILU(0) and the gradient correction.
// False if GMRES does not converge.
fn edge_solve() -> bool {
    guard physics == EDGE && dofs_sized() && edge_nu.len == cell_count() else { abort("edge_solve: an edge physics and the arrays of the dofs"); }
    let start = clock();
    edge_assemble();
    let sizes = wave_sizes();
    let b = edge_rhs();
    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    var converged = solvers::ilu0(row_start, col, val, lu, diag);
    var its = 0;
    let corrected = converged && gradient_setup();
    if converged {
        solvers::rows_start(row_start, col, val);
        its, converged = solvers::gmres(row_start, col, val, b, sol, EDGE_TOLERANCE, wave_restart(b.len), EDGE_ITERATIONS) { r, z =>
            solvers::ilu0_apply(row_start, col, lu, diag, r, z);
            if corrected { gradient_apply(r, z); }
        };
        solvers::rows_stop();
    }
    guard wave_sizes() == sizes else { abort("edge_solve: a solve grows no array of the state (rule 3)"); }
    last_iterations = its;
    last_seconds = clock() - start;
    return converged;
}

// --- values ----------------------------------------------------------------------

// u (curl = false) or curl u in cell c at the barycentric coordinates l:
// the real parts of the 3 components, then the imaginary parts.
fn edge_cell_value(c: i64, l: f64[4], curl: bool) -> f64[6] {
    var w: f64[18] = [0.0; 18];
    var cw: f64[18] = [0.0; 18];
    edge_basis(c, l, w, cw);
    var v: f64[6] = [0.0; 6];
    let le = local_edges();
    for e in le {
        let k = cell_edge[le * c + e];
        for t in 3 {
            let b = if curl { cw[3 * e + t] } else { w[3 * e + t] };
            v[t] += sol[2 * k] * b;
            v[3 + t] += sol[2 * k + 1] * b;
        }
    }
    return v;
}

// u or curl u at point p: the mean of the cells that hold p, else the
// value of the nearest cell. The bool says whether a cell is near p.
fn edge_value_at(p: f64[3], curl: bool) -> f64[6], bool {
    guard physics == EDGE && sol.len == dof_total() else { abort("edge_value_at: an edge physics and a solution"); }
    var sum: f64[6] = [0.0; 6];
    var count = 0;
    var best = -1;
    var best_l = -1.0e300;
    for c in cell_count() {
        let l = barycentric(c, p);
        var lowest = l[0];
        for a in 1..dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
        }
        if lowest < -1.0e-9 { continue; }
        let v = edge_cell_value(c, l, curl);
        for t in 6 { sum[t] += v[t]; }
        count++;
    }
    if count > 0 {
        for t in 6 { sum[t] /= count; }
        return sum, true;
    }
    guard best >= 0 else { return sum, false; };
    return edge_cell_value(best, barycentric(best, p), curl), best_l >= -0.1;
}

// The line integral of component t (or of all, with t < 0, along the line)
// of u or curl u along the line from a to b, by the midpoint rule on n
// pieces. Re and im, and whether each point is near a cell.
fn edge_line(a: f64[3], b: f64[3], t: i64, curl: bool, n: i64) -> f64[2], bool {
    guard t >= -1 && t < 3 && n >= 1 else { abort("edge_line: -1 <= t < 3 and n >= 1"); }
    var s: f64[2] = [0.0, 0.0];
    var all = true;
    var len = 0.0;
    for j in 3 { len += (b[j] - a[j]) * (b[j] - a[j]); }
    len = sqrt(len);
    for i in n {
        let f = (i + 0.5) / n;
        var p: f64[3] = [0.0; 3];
        for j in 3 { p[j] = a[j] + f * (b[j] - a[j]); }
        let v, found = edge_value_at(p, curl);
        if !found { all = false; }
        for part in 2 {
            if t >= 0 {
                s[part] += v[3 * part + t] * len / n;
            } else {
                for j in 3 { s[part] += v[3 * part + j] * (b[j] - a[j]) / n; }
            }
        }
    }
    return s, all;
}

// The nodal mean of u (curl = false) or curl u of the cells at each node:
// 6 values for each node, as in edge_cell_value.
fn edge_nodal(curl: bool) -> f64[>..] {
    var v: f64[>..] = [];
    var count: f64[>..] = [];
    for n in node_count() {
        for t in 6 { v.push(0.0); }
        count.push(0.0);
    }
    for c in cell_count() {
        for a in dim + 1 {
            var l: f64[4] = [0.0; 4];
            l[a] = 1.0;
            let x = edge_cell_value(c, l, curl);
            let n = cell_node(c, a);
            for t in 6 { v[6 * n + t] += x[t]; }
            count[n] += 1.0;
        }
    }
    for n in node_count() {
        for t in 6 { v[6 * n + t] /= max(count[n], 1.0); }
    }
    return v;
}

// Writes the mesh with the nodal means of u and curl u: the real and the
// imaginary parts of each, as 4 vectors.
fn write_edge_vtu(path: const u8[:], field: const u8[:]) -> bool {
    guard physics == EDGE && sol.len == dof_total() else { abort("write_edge_vtu: an edge physics and a solution"); }
    var out: u8[>..] = [];
    vtu_mesh(out);
    for curl in [false, true] {
        let v = edge_nodal(curl);
        for part in 2 {
            format(out, "<DataArray type=\"Float64\" Name=\"", if curl { "curl_" } else { "" }, field, if part == 1 { "_imag" } else { "" },
                   "\" NumberOfComponents=\"3\" format=\"ascii\">\n");
            for n in node_count() { format(out, v[6 * n + 3 * part], " ", v[6 * n + 3 * part + 1], " ", v[6 * n + 3 * part + 2], "\n"); }
            out.append("</DataArray>\n");
        }
    }
    out.append(VTU_END);
    return write_file(path, out);
}

// --- the far field (gap G13) --------------------------------------------------------

// The samples of a closed surface for the far field: at each point, its
// position, the unit normal out of the surface, the weight (an area), and
// E and H as in edge_cell_value.
var far_p: f64[>..] = [];
var far_n: f64[>..] = [];
var far_w: f64[>..] = [];
var far_e: f64[>..] = [];
var far_h: f64[>..] = [];

fn clear_samples() {
    far_p.clear();
    far_n.clear();
    far_w.clear();
    far_e.clear();
    far_h.clear();
}

// The samples of the facets of the region, with E = u and H = hf curl u.
// Each facet is once, and E and H are the means of the cells on its two
// sides. The normal points away from the centroid of the nodes of the
// region, so the surface must be star-shaped about that point. 3D only.
fn surface_samples(hf: f64[2]) {
    guard physics == EDGE && dim == 3 && sol.len == dof_total() && node_count() < 2000000 else { abort("surface_samples: a 3d edge physics with a solution and fewer than 2e6 nodes"); }
    var mid: f64[3] = [0.0; 3];
    var on = 0.0;
    for v in node_count() {
        if !on_region(v) { continue; }
        for t in 3 { mid[t] += node_x(v, t); }
        on += 1.0;
    }
    for t in 3 { mid[t] /= max(on, 1.0); }
    let nv = node_count();
    var at = dictionary<i64, i64> {};
    var sides: i64[>..] = [];   // two (cell, facet) pairs for each facet, -1 for none
    for c in cell_count() {
        for f in facet_count(dim) {
            if !facet_on_region(c, f) { continue; }
            var v: i64[3] = [cell_node(c, facet_node(dim, f, 0)), cell_node(c, facet_node(dim, f, 1)), cell_node(c, facet_node(dim, f, 2))];
            sort(v);
            let k = at.get_or_insert((v[0] * nv + v[1]) * nv + v[2], sides.len / 4);
            if 4 * k == sides.len { sides.append([c, f, -1, -1]); } else { sides[4 * k + 2] = c; sides[4 * k + 3] = f; }
        }
    }
    clear_samples();
    for k in sides.len / 4 { facet_samples(sides[4 * k..4 * k + 4], mid, hf); }
}

// The samples of one facet, from the cells of `side` (cell, facet, cell,
// facet, with -1 for no second cell).
fn facet_samples(side: i64[:], mid: f64[3], hf: f64[2]) {
    let n0, size = facet_normal(side[0], side[1]);
    var n = n0;
    var centre: f64[3] = [0.0; 3];
    for l in 3 {
        for t in 3 { centre[t] += node_x(cell_node(side[0], facet_node(dim, side[1], l)), t) / 3.0; }
    }
    var away = 0.0;
    for t in 3 { away += n[t] * (centre[t] - mid[t]); }
    if away < 0.0 { for t in 3 { n[t] = -n[t]; } }
    let cells = if side[2] >= 0 { 2 } else { 1 };
    for q in quad_count(2) {
        let l, wt = facet_point(side[1], q);
        var p: f64[3] = [0.0; 3];
        for v in 4 {
            for t in 3 { p[t] += l[v] * node_x(cell_node(side[0], v), t); }
        }
        var e: f64[6] = [0.0; 6];
        var h: f64[6] = [0.0; 6];
        for i in cells {
            let c = side[2 * i];
            let b = barycentric(c, p);
            let u = edge_cell_value(c, b, false);
            let cu = edge_cell_value(c, b, true);
            for t in 3 {
                e[t] += u[t] / cells;
                e[3 + t] += u[3 + t] / cells;
                h[t] += (hf[0] * cu[t] - hf[1] * cu[3 + t]) / cells;
                h[3 + t] += (hf[0] * cu[3 + t] + hf[1] * cu[t]) / cells;
            }
        }
        far_p.append(p);
        far_n.append(n);
        far_w.push(wt * size);
        far_e.append(e);
        far_h.append(h);
    }
}

// The radiation intensity U (power per solid angle) in the direction d of
// the samples of a closed surface, by the transform of the surface currents
// J = n × H and M = -n × E (Balanis, Antenna Theory, §12.2-12.3) into a
// medium of wave number k and impedance eta, with the time factor
// exp(i omega t).
fn far_intensity(d: f64[3], k: f64, eta: f64) -> f64 {
    guard k > 0.0 && eta > 0.0 && far_w.len * 3 == far_p.len else { abort("far_intensity: k > 0, eta > 0 and samples"); }
    var r = d;
    let len = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    guard len > 0.0 else { abort("far_intensity: a direction"); }
    for t in 3 { r[t] /= len; }
    // phi = z × r, theta = phi × r; phi = y on the z axis.
    var ph: f64[3] = [-r[1], r[0], 0.0];
    let pl = sqrt(ph[0] * ph[0] + ph[1] * ph[1]);
    ph = if pl > 1.0e-12 { [ph[0] / pl, ph[1] / pl, 0.0] } else { [0.0, 1.0, 0.0] };
    let th: f64[3] = [ph[1] * r[2] - ph[2] * r[1], ph[2] * r[0] - ph[0] * r[2], ph[0] * r[1] - ph[1] * r[0]];
    // N and L projected on theta and phi: re and im of N_th, N_ph, L_th, L_ph.
    var a: f64[8] = [0.0; 8];
    for i in far_w.len {
        var phase = 0.0;
        for t in 3 { phase += k * r[t] * far_p[3 * i + t]; }
        let ex: f64[2] = [cos(phase) * far_w[i], sin(phase) * far_w[i]];
        for part in 2 {
            let off = 6 * i + 3 * part;
            let n = far_n[3 * i..3 * i + 3];
            let j: f64[3] = [n[1] * far_h[off + 2] - n[2] * far_h[off + 1], n[2] * far_h[off] - n[0] * far_h[off + 2], n[0] * far_h[off + 1] - n[1] * far_h[off]];
            let m: f64[3] = [n[2] * far_e[off + 1] - n[1] * far_e[off + 2], n[0] * far_e[off + 2] - n[2] * far_e[off], n[1] * far_e[off] - n[0] * far_e[off + 1]];
            let proj: f64[4] = [dot3(j, th), dot3(j, ph), dot3(m, th), dot3(m, ph)];
            for q in 4 {
                // (part of the current) × exp(i phase) × w
                let z = if part == 0 { ex } else { [-ex[1], ex[0]] };
                a[2 * q] += proj[q] * z[0];
                a[2 * q + 1] += proj[q] * z[1];
            }
        }
    }
    // |L_ph + eta N_th|^2 + |L_th - eta N_ph|^2
    let u1 = (a[6] + eta * a[0]) * (a[6] + eta * a[0]) + (a[7] + eta * a[1]) * (a[7] + eta * a[1]);
    let u2 = (a[4] - eta * a[2]) * (a[4] - eta * a[2]) + (a[5] - eta * a[3]) * (a[5] - eta * a[3]);
    return k * k / (32.0 * PI * PI * eta) * (u1 + u2);
}

fn dot3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// The power out of the closed surface of the samples: ∮ ½ Re(E × H*) · n dA.
fn surface_power() -> f64 {
    var p = 0.0;
    for i in far_w.len {
        let e = far_e[6 * i..6 * i + 6];
        let h = far_h[6 * i..6 * i + 6];
        let n = far_n[3 * i..3 * i + 3];
        // Re(E × conj(H)) = Er × Hr + Ei × Hi
        for part in 2 {
            let o = 3 * part;
            let x: f64[3] = [e[o + 1] * h[o + 2] - e[o + 2] * h[o + 1], e[o + 2] * h[o] - e[o] * h[o + 2], e[o] * h[o + 1] - e[o + 1] * h[o]];
            p += 0.5 * dot3(x, [n[0], n[1], n[2]]) * far_w[i];
        }
    }
    return p;
}
