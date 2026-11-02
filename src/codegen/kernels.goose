// The emitter of the cell kernels (spec §18.3, rule 4 of the Power of
// Ten). A cell kernel computes the matrix of one cell. The emitter writes
// the kernels of each lane into src/<lane>/kernels.goose: the hex27 and
// hex8i solids of the runtime, the simplex cells of src/elements/ and the
// MITC4 shell and the beam of src/shells/.
//
// The emitter writes the integrand of each block of the cell matrix, one
// block for each pair of node groups, with the loops over the components
// unrolled. If the kernel with all its blocks is longer than 60 lines, each
// block becomes a function, and the kernel calls them in order.
//
// The kernels are the same for each model, because they depend only on the
// cell type and the physics. So the emitter writes them into the packages
// once, and not into each generated program. The codegen test fails while
// a kernels.goose file differs from what the emitter writes, and
// `build/test_codegen write` writes the files again.
import std;
namespace codegen;

// The kernel that the emitter writes now: the blocks of its cell matrix.
var block_text: u8[>..] = [];   // the statements of each block, one after the other
var block_ends: i64[>..] = [];  // the end of each block in block_text
var block_rows: i64[>..] = [];  // the node group of the rows of each block
var block_cols: i64[>..] = [];  // the node group of the columns of each block

fn start_kernel() {
    block_text.clear();
    block_ends.clear();
    block_rows.clear();
    block_cols.clear();
}

// Ends block (a, b), whose statements are in block_text from the end of
// the block before it.
fn end_block(a: i64, b: i64) {
    guard a >= 0 && b >= 0 else { abort("end_block: node groups a and b >= 0"); }
    block_ends.push(block_text.len);
    block_rows.push(a);
    block_cols.push(b);
}

fn block_start(k: i64) -> i64 { return if k == 0 { 0 } else { block_ends[k - 1] }; }

fn line_count(s: u8[:]) -> i64 {
    var n = 0;
    for c in s { if c == '\n' { n++; } }
    return n;
}

// Appends the lines of `s` with `depth` more levels of indent.
fn indented(out: u8[>..]&, s: u8[:], depth: i64) {
    var at_start = true;
    for c in s {
        if at_start && c != '\n' { for i in 4 * depth { out.push(' '); } }
        out.push(c);
        at_start = c == '\n';
    }
}

// Writes the kernel `name`: `comment`, then `fn name(params) {`, the
// statements `head`, the blocks at `depth` more levels of indent than the
// body, the statements `tail` and `}`. If that is longer than MAX_FN_LINES,
// block (a, b) becomes the function name_a_b(block_params), with the guard
// `block_guard` (a condition), and the kernel calls it with `block_args`.
fn write_kernel(out: u8[>..]&, comment: u8[:], name: u8[:], params: u8[:], head: u8[:], depth: i64, tail: u8[:],
                block_params: u8[:], block_guard: u8[:], block_args: u8[:]) {
    guard block_ends.len > 0 && depth >= 0 else { abort("write_kernel: a block, and depth >= 0"); }
    let inline_lines = 2 + line_count(head) + line_count(block_text) + line_count(tail);
    if inline_lines > MAX_FN_LINES {
        for k in block_ends.len {
            let a = block_rows[k];
            let b = block_cols[k];
            let body = block_text[block_start(k)..block_ends[k]];
            guard 3 + line_count(body) <= MAX_FN_LINES else { abort("write_kernel: a block has more than 60 lines"); }
            format(out, "// Block (", a, ", ", b, ") of ", name, ": the rows of node group ", a, ", the columns of node group ", b, ".\n");
            format(out, "fn ", name, "_", a, "_", b, "(", block_params, ") {\n");
            format(out, "    guard ", block_guard, " else { abort(\"", name, "_", a, "_", b, ": ", block_guard, "\"); }\n");
            format(out, body, "}\n\n");
        }
    }
    format(out, comment, "fn ", name, "(", params, ") {\n", head);
    if inline_lines > MAX_FN_LINES {
        var calls: u8[>..] = [];
        for k in block_ends.len { format(calls, "    ", name, "_", block_rows[k], "_", block_cols[k], "(", block_args, ");\n"); }
        indented(out, calls, depth);
    } else {
        indented(out, block_text, depth);
    }
    format(out, tail, "}\n");
    guard line_count(head) + line_count(tail) + 2 + (if inline_lines > MAX_FN_LINES { block_ends.len } else { 0 }) <= MAX_FN_LINES else { abort("write_kernel: the kernel has more than 60 lines"); }
}

// --- the blocks of linear elasticity and of heat conduction ---------------------

// The block of K(i p, j q) = λ ∂p Ni ∂q Nj + μ (∂q Ni ∂p Nj + δpq ∇Ni · ∇Nj)
// times dv, for the nodes i in a_lo..a_hi and j in b_lo..b_hi, with `dim`
// components and `size` columns of the matrix `k`. `lambda` and `mu` name
// the Lamé constants. With `sum`, ∇Ni · ∇Nj is a sum of `dim` terms from
// 0.0; else it is one expression of 3 terms.
fn elastic_block(lo_a: i64, hi_a: i64, lo_b: i64, hi_b: i64, dim: i64, size: i64, k: u8[:], lambda: u8[:], mu: u8[:],
                 sum: bool) {
    guard lo_a < hi_a && lo_b < hi_b && (dim == 2 || dim == 3) && size >= dim * hi_a && size >= dim * hi_b else { abort("elastic_block: node ranges, dim 2 or 3, and size >= dim nodes"); }
    format(block_text, "    for i in ", lo_a, "..", hi_a, " {\n");
    format(block_text, "        for j in ", lo_b, "..", hi_b, " {\n");
    gradient_dot(dim, sum);
    format(block_text, "            let row = ", size * dim, " * i + ", dim, " * j;\n");
    for p in dim {
        for q in dim {
            let v = if p == 0 && q == 0 { "var v" } else { "v" };
            format(block_text, "            ", v, " = ", lambda, " * grad[", component(p, "i"), "] * grad[", component(q, "j"), "] + ", mu,
                   " * grad[", component(q, "i"), "] * grad[", component(p, "j"), "];\n");
            if p == q { format(block_text, "            v += ", mu, " * dot;\n"); }
            let at: u8[>..] = if p == 0 && q == 0 { "row" } else { str("row + ", size * p + q) };
            format(block_text, "            ", k, "[", at, "] += v * dv;\n");
        }
    }
    block_text.append("        }\n    }\n");
}

// The block of K(i, j) = k ∇Ni · ∇Nj times dv, for the nodes i in
// a_lo..a_hi and j in b_lo..b_hi of a matrix with `size` columns.
fn heat_block(lo_a: i64, hi_a: i64, lo_b: i64, hi_b: i64, dim: i64, size: i64) {
    guard lo_a < hi_a && lo_b < hi_b && (dim == 2 || dim == 3) && size >= hi_a && size >= hi_b else { abort("heat_block: node ranges, dim 2 or 3, and size >= the nodes"); }
    format(block_text, "    for i in ", lo_a, "..", hi_a, " {\n");
    format(block_text, "        for j in ", lo_b, "..", hi_b, " {\n");
    gradient_dot(dim, true);
    format(block_text, "            ke[", size, " * i + j] += kq * dot * dv;\n");
    block_text.append("        }\n    }\n");
}

// `dot` = ∇Ni · ∇Nj from grad[3 n + t] = ∂Nn/∂x_t.
fn gradient_dot(dim: i64, sum: bool) {
    guard dim == 2 || dim == 3 else { abort("gradient_dot: dim 2 or 3"); }
    if !sum {
        block_text.append("            let dot = grad[3 * i] * grad[3 * j] + grad[3 * i + 1] * grad[3 * j + 1]\n");
        block_text.append("                + grad[3 * i + 2] * grad[3 * j + 2];\n");
        return;
    }
    block_text.append("            var dot = 0.0;\n");
    for t in dim { format(block_text, "            dot += grad[", component(t, "i"), "] * grad[", component(t, "j"), "];\n"); }
}

// The index 3 n + t of grad, for the node `n` (a name) and the component t.
fn component(t: i64, n: u8[:]) -> u8[>..] {
    guard t >= 0 && t < 3 else { abort("component: 0 <= t < 3"); }
    return if t == 0 { str("3 * ", n) } else { str("3 * ", n, " + ", t) };
}

// The first node and the end of node group g, for `nodes` nodes in groups
// that start at `starts`.
fn group_bounds(starts: i64[:], nodes: i64, g: i64) -> i64, i64 {
    guard g >= 0 && g < starts.len else { abort("group_bounds: 0 <= g < starts.len"); }
    return starts[g], if g + 1 < starts.len { starts[g + 1] } else { nodes };
}

// --- the runtime ----------------------------------------------------------------

let RUNTIME_KERNELS = "src/runtime/kernels.goose";

fn file_header(out: u8[>..]&, what: u8[:], ns: u8[:], imports: u8[:]) {
    format(out, "// ", what, ", written by the emitter of the cell kernels\n");
    out.append("// (src/codegen/kernels.goose). Do not edit this file: edit the emitter,\n");
    out.append("// then run `build/test_codegen write` from the root of the repository.\n");
    format(out, "import std;\nnamespace ", ns, ";\n", imports, "\n");
}

fn runtime_kernels() -> u8[>..] {
    var out: u8[>..] = [];
    file_header(out, "The cell kernels of the runtime", "runtime", "import .fem;\nimport .solid;\nimport .hex8i;\n");
    hex27_kernel(out);
    out.push('\n');
    hex8i_kernel(out);
    return out;
}

// The gradients of the shape functions of `nodes` nodes at a quadrature
// point, after the shape function derivatives of the cell into shape_dn.
fn solid_gradients(out: u8[>..]&, nodes: i64) {
    out.append("        var jac: f64[9] = [0.0; 9];\n");
    out.append("        cell_jacobian(c, shape_dn, jac);\n");
    out.append("        var inv: f64[9] = [0.0; 9];\n");
    out.append("        let det = inverse3(jac, inv);\n");
    out.append("        assert(det > 0.0);\n");
    format(out, "        for l in ", nodes, " {\n");
    out.append("            for s in 3 {\n");
    out.append("                grad[3 * l + s] = inv[3 * s] * shape_dn[3 * l] + inv[3 * s + 1] * shape_dn[3 * l + 1]\n");
    out.append("                    + inv[3 * s + 2] * shape_dn[3 * l + 2];\n");
    out.append("            }\n        }\n");
}

// The Lamé constants of the material of cell c (solid.goose).
let LAME = "    let lambda, mu = cell_lame(c);\n";

let ELASTIC_PARAMS = "grad: f64[:], lambda: f64, mu: f64, dv: f64, ke: f64[>..]&";

// The groups of the 27 nodes of hex27: the layers ξ_2 = -1, 0 and 1.
fn hex27_kernel(out: u8[>..]&) {
    start_kernel();
    let starts = [0, 9, 18];
    for a in 3 {
        for b in 3 {
            let lo_a, hi_a = group_bounds(starts, 27, a);
            let lo_b, hi_b = group_bounds(starts, 27, b);
            elastic_block(lo_a, hi_a, lo_b, hi_b, 3, 81, "ke", "lambda", "mu", false);
            end_block(a, b);
        }
    }
    var head: u8[>..] = [];
    head.append("    guard c >= 0 && c < cell_count() && cell_type[c] as i64 == HEX27 && cell_start[c + 1] - cell_start[c] == 27 else { abort(\"hex27_stiffness: a hex27 cell with 27 nodes\"); }\n");
    head.append("    guard cell_quadrature(HEX27, quad_xi, quad_w) else { return false; }\n");
    head.append("    ke.clear();\n    for i in 81 * 81 { ke.push(0.0); }\n");
    head.append("    grad.clear();\n    for i in 81 { grad.push(0.0); }\n");
    head.append(LAME);
    head.append("    for w, g in quad_w {\n");
    head.append("        let xi: f64[3] = [quad_xi[3 * g], quad_xi[3 * g + 1], quad_xi[3 * g + 2]];\n");
    head.append("        guard cell_shape(HEX27, xi, shape_n, shape_dn) else { return false; }\n");
    solid_gradients(head, 27);
    head.append("        let dv = w * det;\n");
    write_kernel(out, "// The stiffness of the hex27 cell c, with 3 dofs for each node:\n// K(i p, j q) = ∫ λ ∂p Ni ∂q Nj + μ (∂q Ni ∂p Nj + δpq ∇Ni · ∇Nj) dV.\n",
                 "hex27_stiffness", "c: i64, ke: f64[>..]&", head, 1, "    }\n    return true;\n",
                 ELASTIC_PARAMS, "grad.len == 81 && ke.len == 81 * 81", "grad, lambda, mu, dv, ke");
}

// The groups of the 11 shape functions of hex8i: the 8 nodes and the 3
// modes.
fn hex8i_kernel(out: u8[>..]&) {
    start_kernel();
    let starts = [0, 8];
    for a in 2 {
        for b in 2 {
            let lo_a, hi_a = group_bounds(starts, 11, a);
            let lo_b, hi_b = group_bounds(starts, 11, b);
            elastic_block(lo_a, hi_a, lo_b, hi_b, 3, 33, "ke", "lambda", "mu", false);
            end_block(a, b);
        }
    }
    var head: u8[>..] = [];
    head.append("    guard c >= 0 && c < cell_count() && cell_type[c] as i64 == HEX8I && cell_start[c + 1] - cell_start[c] == 8 else { abort(\"hex8i_full: a hex8i cell with 8 nodes\"); }\n");
    head.append("    quad_xi.clear();\n    quad_w.clear();\n    gauss_hex2(quad_xi, quad_w);\n");
    head.append("    full.clear();\n    for i in ALL_DOFS * ALL_DOFS { full.push(0.0); }\n");
    head.append("    grad.clear();\n    for i in ALL_DOFS { grad.push(0.0); }\n");
    head.append(LAME);
    head.append("    for w, g in quad_w {\n");
    head.append("        let xi: f64[3] = [quad_xi[3 * g], quad_xi[3 * g + 1], quad_xi[3 * g + 2]];\n");
    head.append("        shape_n.clear();\n        shape_dn.clear();\n        hex8_shape(xi, shape_n, shape_dn);\n");
    solid_gradients(head, 8);
    head.append("        mode_gradient(c, xi, det, grad, NODE_DOFS);\n");
    head.append("        let dv = w * det;\n");
    write_kernel(out, "// The 33 × 33 matrix `full` of the 8 nodes and the 3 modes of the hex8i\n// cell c, with the integrand of hex27_stiffness.\n",
                 "hex8i_full", "c: i64", head, 1, "    }\n",
                 ELASTIC_PARAMS, "grad.len == 33 && ke.len == 33 * 33", "grad, lambda, mu, dv, full");
}

// --- the elements lane ------------------------------------------------------------

let ELEMENTS_KERNELS = "src/elements/kernels.goose";

let SIMPLEX_NAMES = ["tri3", "tri6", "tet4", "tet10"];

fn elements_kernels() -> u8[>..] {
    var out: u8[>..] = [];
    file_header(out, "The cell kernels of the simplex cells", "elements", "import .shape;\nimport .mesh;\nimport .assemble;\n");
    out.append("// The stiffness or the conductivity of cell c, and its source vector:\n");
    out.append("// elastic K(i p, j q) = ∫ λ ∂p Ni ∂q Nj + μ (∂q Ni ∂p Nj + δpq ∇Ni · ∇Nj) dV,\n");
    out.append("// f(i p) = ∫ b_p Ni dV, heat K(i, j) = ∫ k ∇Ni · ∇Nj dV and f(i) = ∫ q_v Ni dV.\n");
    out.append("fn cell_matrix(c: i64) {\n");
    out.append("    guard c >= 0 && c < cell_count() && (physics == ELASTIC || physics == HEAT) && (dim == 2 || dim == 3) && (order == 1 || order == 2) else { abort(\"cell_matrix: c < cell_count(), a physics, dim 2 or 3 and order 1 or 2\"); }\n");
    out.append("    if hex_cells {\n        if physics == ELASTIC { elastic_hex8(c); } else { heat_hex8(c); }\n        return;\n    }\n");
    for physics in 2 {
        let word = if physics == 0 { "ELASTIC" } else { "HEAT" };
        for dim in 2..4 {
            for order in 1..3 {
                format(out, "    if physics == ", word, " && dim == ", dim, " && order == ", order, " { ", simplex_kernel_name(physics, dim, order), "(c); }\n");
            }
        }
    }
    out.append("}\n");
    for physics in 2 {
        out.push('\n');
        cube_kernel(out, physics);
    }
    for physics in 2 {
        for dim in 2..4 {
            for order in 1..3 {
                out.push('\n');
                simplex_kernel(out, physics, dim, order);
            }
        }
    }
    return out;
}

// The 8-node hexahedron of order 1: one group of nodes, and the integrand of a simplex.
fn cube_kernel(out: u8[>..]&, physics: i64) {
    simplex_cell_kernel(out, physics, 3, 1, true);
}

// physics 0 is elastic, 1 heat.
fn simplex_kernel_name(physics: i64, dim: i64, order: i64) -> u8[>..] {
    guard (physics == 0 || physics == 1) && (dim == 2 || dim == 3) && (order == 1 || order == 2) else { abort("simplex_kernel_name: a physics, dim 2 or 3 and order 1 or 2"); }
    return str(if physics == 0 { "elastic_" } else { "heat_" }, SIMPLEX_NAMES[2 * (dim - 2) + order - 1]);
}

// The groups of the nodes of a simplex: its vertices, then the nodes at
// the midpoints of its edges (order 2).
fn simplex_kernel(out: u8[>..]&, physics: i64, dim: i64, order: i64) {
    simplex_cell_kernel(out, physics, dim, order, false);
}

fn simplex_cell_kernel(out: u8[>..]&, physics: i64, dim: i64, order: i64, cube: bool) {
    let name = if cube { str(if physics == 0 { "elastic_" } else { "heat_" }, "hex8") } else { simplex_kernel_name(physics, dim, order) };
    let nodes = if cube { 8 } else if order == 1 { dim + 1 } else { (dim + 1) * (dim + 2) / 2 };
    let comps = if physics == 0 { dim } else { 1 };
    let size = nodes * comps;
    let groups = if order == 1 { 1 } else { 2 };
    let starts = [0, dim + 1];
    start_kernel();
    for a in groups {
        for b in groups {
            let lo_a, hi_a = group_bounds(starts[..groups], nodes, a);
            let lo_b, hi_b = group_bounds(starts[..groups], nodes, b);
            if physics == 0 {
                elastic_block(lo_a, hi_a, lo_b, hi_b, dim, size, "ke", "lame", "shear", true);
            } else {
                heat_block(lo_a, hi_a, lo_b, hi_b, dim, size);
            }
            end_block(a, b);
        }
    }
    var head: u8[>..] = [];
    let word = if physics == 0 { "ELASTIC" } else { "HEAT" };
    format(head, "    guard c >= 0 && c < cell_count() && physics == ", word, " && dim == ", dim, " && per_cell == ", nodes, " && ncomp == ", comps,
           " else { abort(\"", name, ": a cell, the physics, dim ", dim, ", ", nodes, " nodes and ", comps, " components\"); }\n");
    head.append("    for &v in ke { v = 0.0; }\n    for &v in fe { v = 0.0; }\n");
    format(head, "    for qp in quad_count(", dim, ") {\n");
    format(head, "        let xi, w = quad_point(", dim, ", qp);\n");
    head.append("        let dv = w * cell_geometry(c, xi) * thickness;\n");
    if physics == 1 { head.append("        let kq = conductivity_at(c);\n"); }
    format(head, "        for i in ", nodes, " {\n");
    if physics == 1 {
        head.append("            fe[i] += source * shp[i] * dv;\n");
    } else {
        for p in dim { format(head, "            fe[", if p == 0 { str(dim, " * i") } else { str(dim, " * i + ", p) }, "] += volume_force[", p, "] * shp[i] * dv;\n"); }
    }
    head.append("        }\n");
    write_kernel(out, "", name, "c: i64", head, 1, "    }\n", "dv: f64",
                 str("per_cell == ", nodes, " && ncomp == ", comps), "dv");
}

// --- the shells lane --------------------------------------------------------------

let SHELLS_KERNELS = "src/shells/kernels.goose";

fn shells_kernels() -> u8[>..] {
    var out: u8[>..] = [];
    file_header(out, "The cell kernels of the shells and the beams", "shells", "import .model;\nimport .mitc4;\nimport .beam;\n");
    shell_kernel(out);
    out.push('\n');
    beam_kernel(out);
    return out;
}

// The groups of the 4 nodes of MITC4: the sides η = -1 (nodes 0 and 1) and
// η = 1 (nodes 2 and 3). The block of K(i, j) = (kd Di Dj + Σk rows_k,i
// gen_k,j) area, for the 6 dofs of each node of the groups.
fn shell_kernel(out: u8[>..]&) {
    start_kernel();
    for a in 2 {
        for b in 2 {
            format(block_text, "    for i in ", 12 * a, "..", 12 * a + 12, " {\n");
            format(block_text, "        for j in ", 12 * b, "..", 12 * b + 12, " {\n");
            block_text.append("            var v = kd * rows[192 + i] * rows[192 + j];\n");
            for k in 8 {
                let at: u8[>..] = if k == 0 { "" } else { str(24 * k, " + ") };
                format(block_text, "            v += rows[", at, "i] * gen[", at, "j];\n");
            }
            block_text.append("            ke[24 * i + j] += v * area;\n");
            block_text.append("        }\n    }\n");
            end_block(a, b);
        }
    }
    var head: u8[>..] = [];
    head.append("    guard c >= 0 && c < shells.len && shells[c].section >= 0 && shells[c].section < sections.len else { abort(\"shell_stiffness: 0 <= c < shells.len, and its section < sections.len\"); }\n");
    head.append("    load_cell(c);\n");
    head.append("    let s = sections[shells[c].section];\n");
    head.append("    let kd = DRILL * s.a[8];\n");
    head.append("    for &v in ke { v = 0.0; }\n");
    head.append("    var gen: f64[192] = [0.0; 192];      // N, M and Q of each dof column j: 8 rows of 24\n");
    head.append("    for gp in 4 {\n");
    head.append("        let xi = if gp % 2 == 0 { -GAUSS } else { GAUSS };\n");
    head.append("        let eta = if gp < 2 { -GAUSS } else { GAUSS };\n");
    head.append("        local_rows(c, xi, eta);\n");
    head.append("        for j in 24 {\n");
    head.append("            let am = mat3_row(s.a, 0, j);\n");
    head.append("            let bk = mat3_row(s.b, 72, j);\n");
    head.append("            let bm = mat3_row(s.b, 0, j);\n");
    head.append("            let dk = mat3_row(s.d, 72, j);\n");
    head.append("            for k in 3 {\n");
    head.append("                gen[24 * k + j] = am[k] + bk[k];\n");
    head.append("                gen[72 + 24 * k + j] = bm[k] + dk[k];\n");
    head.append("            }\n");
    head.append("            gen[144 + j] = s.s[0] * rows[144 + j] + s.s[1] * rows[168 + j];\n");
    head.append("            gen[168 + j] = s.s[2] * rows[144 + j] + s.s[3] * rows[168 + j];\n");
    head.append("        }\n");
    write_kernel(out, "// The stiffness ke of shell cell c, by 2 × 2 Gauss points.\n", "shell_stiffness", "c: i64", head, 1, "    }\n",
                 "kd: f64, gen: f64[:]", "gen.len == 192", "kd, gen");
}

// The beam has one group of nodes: its 2 nodes. The block is out = tbᵀ kt.
fn beam_kernel(out: u8[>..]&) {
    start_kernel();
    block_text.append("    for i in 12 {\n");
    block_text.append("        for j in 12 {\n");
    block_text.append("            var v = 0.0;\n");
    block_text.append("            for k in 12 { v += tb[12 * k + i] * kt[12 * k + j]; }\n");
    block_text.append("            out[12 * i + j] = v;\n");
    block_text.append("        }\n    }\n");
    end_block(0, 0);
    var head: u8[>..] = [];
    head.append("    guard b >= 0 && b < beams.len && out.len == 144 else { abort(\"beam_stiffness: 0 <= b < beams.len, and out has 144 entries\"); }\n");
    head.append("    beam_local(b);\n");
    head.append("    var kt: f64[144] = [0.0; 144];\n");
    head.append("    for i in 12 {\n");
    head.append("        for j in 12 {\n");
    head.append("            var v = 0.0;\n");
    head.append("            for k in 12 { v += kb[12 * i + k] * tb[12 * k + j]; }\n");
    head.append("            kt[12 * i + j] = v;\n");
    head.append("        }\n    }\n");
    write_kernel(out, "// The stiffness of beam b in global dofs: tbᵀ kb tb, into out[0..144].\n", "beam_stiffness", "b: i64, out: f64[:]", head, 0, "",
                 "kt: f64[:], out: f64[:]", "kt.len == 144 && out.len == 144", "kt, out");
}

// --- the files ----------------------------------------------------------------------

let KERNEL_FILES = [RUNTIME_KERNELS, ELEMENTS_KERNELS, SHELLS_KERNELS];

// The text that the emitter writes into KERNEL_FILES[f].
fn kernel_file(f: i64) -> u8[>..] {
    guard f >= 0 && f < KERNEL_FILES.len else { abort("kernel_file: 0 <= f < KERNEL_FILES.len"); }
    if f == 0 { return runtime_kernels(); }
    if f == 1 { return elements_kernels(); }
    return shells_kernels();
}
