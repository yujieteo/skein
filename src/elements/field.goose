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
let OP_NEG = 14;
let OP_SQRT = 15;

let KEEP_SLOTS = 8;

var ft_kind: i64[>..] = [];
var ft_law: i64[>..] = [];
var ft_code_at: i64[>..] = [0];
var ft_const_at: i64[>..] = [0];
var ft_code: i64[>..] = [];
var ft_consts: f64[>..] = [];

var kept_vals: f64[>..] = [];
var kept_comps: i64[>..] = [];

var ev_stack: f64[16];

fn clear_field_terms() {
    ft_kind.clear();
    ft_law.clear();
    ft_code_at.clear();
    ft_code_at.push(0);
    ft_const_at.clear();
    ft_const_at.push(0);
    ft_code.clear();
    ft_consts.clear();
}

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

fn add_field_term(kind: i64, code: i64[:], consts: f64[:]) { add_field_term_in(-1, kind, code, consts); }

fn add_field_term_in(law: i64, kind: i64, code: i64[:], consts: f64[:]) {
    guard (kind == FT_SOURCE || kind == FT_STRAIN) && code.len % 3 == 0 && code.len > 0 else { abort("add_field_term: a kind and a code of 3 values for each instruction"); }
    guard kind != FT_SOURCE || physics == HEAT else { abort("add_field_term: a source needs a heat physics"); }
    guard kind != FT_STRAIN || (physics == ELASTIC && (dim == 3 || plane == PLANE_STRESS)) else { abort("add_field_term: a strain needs a 3D or plane stress elastic physics"); }
    ft_kind.push(kind);
    ft_law.push(law);
    ft_code.append(code);
    ft_consts.append(consts);
    ft_code_at.push(ft_code.len);
    ft_const_at.push(ft_consts.len);
}

fn kept_at(slot: i64, n: i64, p: i64) -> f64 {
    guard slot >= 0 && slot < KEEP_SLOTS && kept_vals.len == KEEP_SLOTS * 3 * node_count() && p >= 0 && p < kept_comps[slot] else { abort("kept_at: a kept slot and a component of it"); }
    return kept_vals[slot * 3 * node_count() + 3 * n + p];
}

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

fn eval_facet_code(code: i64[:], consts: f64[:], x: f64[3], u: f64[3]) -> f64 {
    var sp = 0;
    var i = 0;
    while i < code.len {
        let op = code[i];
        let a = code[i + 1];
        let b = code[i + 2];
        i += 3;
        if op == OP_CONST {
            ev_stack[sp] = consts[a];
            sp++;
        } else if op == OP_VAL {
            guard a == KEEP_SLOTS && b >= 0 && b < 3 else { abort("eval_facet_code: the solution now"); }
            ev_stack[sp] = u[b];
            sp++;
        } else if op == OP_COORD {
            guard a >= 0 && a < 3 else { abort("eval_facet_code: a coordinate 0 to 2"); }
            ev_stack[sp] = x[a];
            sp++;
        } else if op == OP_NEG {
            ev_stack[sp - 1] = -ev_stack[sp - 1];
        } else if op == OP_SQRT {
            ev_stack[sp - 1] = sqrt(max(ev_stack[sp - 1], 0.0));
        } else {
            guard op != OP_GRAD && sp >= 2 && sp < ev_stack.len else { abort("eval_facet_code: two values for a binary instruction"); }
            let y = ev_stack[sp - 1];
            let z = ev_stack[sp - 2];
            sp--;
            ev_stack[sp - 1] = if op == OP_ADD { z + y } else if op == OP_SUB { z - y } else if op == OP_MUL { z * y } else { z / y };
        }
        guard sp >= 1 && sp <= ev_stack.len else { abort("eval_facet_code: a stack of 1 to 16 values"); }
    }
    return ev_stack[0];
}

fn strain_factor() -> f64 { return if dim == 3 { 3.0 * lame + 2.0 * shear } else { 2.0 * (lame + shear) }; }

fn strain_stress(c: i64) -> f64 {
    var s = 0.0;
    for t in ft_kind.len { if ft_kind[t] == FT_STRAIN && term_acts(t, c) { s += strain_factor() * eval_term(t, c); } }
    return s;
}

fn term_acts(t: i64, c: i64) -> bool {
    guard ft_law[t] >= 0 else { return true; };
    return ft_law[t] == (if cell_law.len == 0 { 0 } else { cell_law[c] });
}

fn add_field_terms(c: i64) {
    for t in ft_kind.len {
        guard term_acts(t, c) else { continue; };
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
