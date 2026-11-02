// Field terms: the solution of one physics as a source, or a strain, of the
// next physics on the same mesh.
//
// After a solve, `keep(slot)` holds the nodal solution as the field of
// `slot`. A later physics adds a field term: an expression of the kept
// fields, written as the postfix code of a small stack machine that the
// assembly runs at each quadrature point:
//
//     OP_CONST a      push consts[a]
//     OP_VAL a b      push component b of the field of slot a
//     OP_GRAD a b     push d(component b / 4)/dx_(b % 4) of the field of slot a
//     OP_COORD a      push the coordinate a of the point
//     OP_ADD .. OP_SQRT   pop one or two values and push the result
//
// A term has a kind: FT_SOURCE adds ∫ s N_i to the heat source vector, and
// FT_STRAIN is an isotropic strain θ that adds ∫ (3 λ + 2 μ) θ ∂p N_i to the
// force of an elastic physics, and that the stress subtracts. The terms
// belong to the physics: the next `heat` or `elastic` call clears them.
import std;
import math;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;

let FT_SOURCE = 0;
let FT_STRAIN = 1;

let OP_CONST = 0;
let OP_VAL = 1;
let OP_GRAD = 2;
let OP_COORD = 3;
let OP_ADD = 10;
let OP_SUB = 11;
let OP_MUL = 12;
let OP_DIV = 13;
let OP_NEG = 14;
let OP_SQRT = 15;

let KEEP_SLOTS = 8;

var ft_kind: i64[>..] = [];
var ft_code_at: i64[>..] = [0];     // the code of term t is ft_code[ft_code_at[t]..ft_code_at[t + 1]], 3 values for each instruction
var ft_const_at: i64[>..] = [0];
var ft_code: i64[>..] = [];
var ft_consts: f64[>..] = [];

var kept_vals: f64[>..] = [];       // the field of slot s at node n, component p: kept_vals[s * 3 * nodes + 3 n + p]
var kept_comps: i64[>..] = [];

var ev_stack: f64[16];

fn clear_field_terms() {
    ft_kind.clear();
    ft_code_at.clear();
    ft_code_at.push(0);
    ft_const_at.clear();
    ft_const_at.push(0);
    ft_code.clear();
    ft_consts.clear();
}

// Keeps the solution of the physics now as the field of `slot`.
fn keep_field(slot: i64) {
    guard slot >= 0 && slot < KEEP_SLOTS && (physics == ELASTIC || physics == HEAT) && sol.len == dof_total() else { abort("keep_field: a slot, a heat or elastic physics and its solution"); }
    let nn = node_count();
    if kept_vals.len != KEEP_SLOTS * 3 * nn {
        kept_vals.clear();
        for i in KEEP_SLOTS * 3 * nn { kept_vals.push(0.0); }
        kept_comps.clear();
        for s in KEEP_SLOTS { kept_comps.push(0); }
    }
    kept_comps[slot] = ncomp;
    for n in nn {
        for p in ncomp { kept_vals[slot * 3 * nn + 3 * n + p] = sol[ncomp * n + p]; }
    }
}

// Adds a term of kind `kind` with the instructions `code` and the
// constants `consts`.
fn add_field_term(kind: i64, code: i64[:], consts: f64[:]) {
    guard (kind == FT_SOURCE || kind == FT_STRAIN) && code.len % 3 == 0 && code.len > 0 else { abort("add_field_term: a kind and a code of 3 values for each instruction"); }
    guard kind != FT_SOURCE || physics == HEAT else { abort("add_field_term: a source needs a heat physics"); }
    guard kind != FT_STRAIN || (physics == ELASTIC && (dim == 3 || plane == PLANE_STRESS)) else { abort("add_field_term: a strain needs a 3D or plane stress elastic physics"); }
    ft_kind.push(kind);
    ft_code.append(code);
    ft_consts.append(consts);
    ft_code_at.push(ft_code.len);
    ft_const_at.push(ft_consts.len);
}

fn kept_at(slot: i64, n: i64, p: i64) -> f64 {
    guard slot >= 0 && slot < KEEP_SLOTS && kept_vals.len == KEEP_SLOTS * 3 * node_count() && p >= 0 && p < kept_comps[slot] else { abort("kept_at: a kept slot and a component of it"); }
    return kept_vals[slot * 3 * node_count() + 3 * n + p];
}

// The value of term t at the point of the last cell_geometry of cell c.
fn eval_term(t: i64, c: i64) -> f64 {
    var sp = 0;
    var i = ft_code_at[t];
    let stop = ft_code_at[t + 1];
    let k0 = ft_const_at[t];
    while i < stop {
        let op = ft_code[i];
        let a = ft_code[i + 1];
        let b = ft_code[i + 2];
        i += 3;
        if op == OP_CONST {
            ev_stack[sp] = ft_consts[k0 + a];
            sp++;
        } else if op == OP_VAL {
            var v = 0.0;
            for l in per_cell { v += shp[l] * kept_at(a, cell_node(c, l), b); }
            ev_stack[sp] = v;
            sp++;
        } else if op == OP_GRAD {
            var g = 0.0;
            for l in per_cell { g += grad[3 * l + b % 4] * kept_at(a, cell_node(c, l), b / 4); }
            ev_stack[sp] = g;
            sp++;
        } else if op == OP_COORD {
            var x = 0.0;
            for l in per_cell { x += shp[l] * node_x(cell_node(c, l), a); }
            ev_stack[sp] = x;
            sp++;
        } else if op == OP_NEG {
            ev_stack[sp - 1] = -ev_stack[sp - 1];
        } else if op == OP_SQRT {
            ev_stack[sp - 1] = sqrt(max(ev_stack[sp - 1], 0.0));
        } else {
            guard sp >= 2 && sp < ev_stack.len else { abort("eval_term: two values for a binary instruction"); }
            let y = ev_stack[sp - 1];
            let x = ev_stack[sp - 2];
            sp--;
            ev_stack[sp - 1] = if op == OP_ADD { x + y } else if op == OP_SUB { x - y } else if op == OP_MUL { x * y } else { x / y };
        }
        guard sp >= 1 && sp <= ev_stack.len else { abort("eval_term: a stack of 1 to 16 values"); }
    }
    return ev_stack[0];
}

// The pressure-like factor of the isotropic strain θ: the stress of
// -θ I is -(3 λ + 2 μ) θ in 3D, and -2 (λ + μ) θ in plane stress.
fn strain_factor() -> f64 { return if dim == 3 { 3.0 * lame + 2.0 * shear } else { 2.0 * (lame + shear) }; }

// The stress of the strain terms of cell c at the point of the last cell_geometry.
fn strain_stress(c: i64) -> f64 {
    var s = 0.0;
    for t in ft_kind.len { if ft_kind[t] == FT_STRAIN { s += strain_factor() * eval_term(t, c); } }
    return s;
}

// Adds the field terms to the source fe of cell c, after cell_matrix.
fn add_field_terms(c: i64) {
    for t in ft_kind.len {
        for qp in quad_count(dim) {
            let xi, w = quad_point(dim, qp);
            let dv = w * cell_geometry(c, xi) * thickness;
            let s = eval_term(t, c);
            if ft_kind[t] == FT_SOURCE {
                for i in per_cell { fe[i] += s * shp[i] * dv; }
            } else {
                let k = strain_factor() * s;
                for i in per_cell {
                    for p in dim { fe[dim * i + p] += k * grad[3 * i + p] * dv; }
                }
            }
        }
    }
}
