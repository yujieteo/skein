// The rules of the items (spec §3 to §22, in the order of §19): names,
// declarations, materials, physics, conditions, methods, steps, results,
// expectations, parameters, scales and studies.
import std;
import dictionary;
import math;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .expr;
import .builtins;
import .typeexpr;

// --- names (spec §4) ---------------------------------------------------------------------

// Adds a declaration to the namespace, or reports `duplicate_name`.
fn declare(kind: i64, tok: i64, node: i64, ty: Ty) -> i64 {
    let name = tok_text(tok);
    let old = visible.get(key(name));
    if old {
        let o = syms[old];
        if o.kind == S_BUILTIN {
            error(tok, "duplicate_name", str("'", name, "' is a built-in name (Appendix A)"),
                  str("rename '", name, "'"));
        } else {
            error(tok, "duplicate_name",
                  str("'", name, "' is already declared on line ", syntax::tokens[o.tok].line),
                  str("rename one of the two declarations of '", name, "'"));
        }
    }
    let s = add_sym(kind, tok, node, ty);
    visible.insert(key(name), s);
    return s;
}

fn add_builtins() {
    let consts = BUILTIN_CONSTS;
    for name in consts {
        var t = any_ty();
        if name == "I" { t = num_ty(NO_DIM, 2); }
        if name == "pi" { t = base_ty(B_REAL); }
        if name == "x" || name == "y" || name == "z" { t = num_ty(LENGTH, 0); t.varying = true; }
        if name == "t" { t = num_ty(TIME, 0); t.varying = true; }
        if name == "n" { t = num_ty(NO_DIM, 1); t.varying = true; }
        if name == "im" { t = base_ty(B_REAL); t.num = NS_COMPLEX; }
        if name == "dual_eps" { t = base_ty(B_REAL); t.num = NS_DUAL; }
        syms.push(Sym { kind: S_BUILTIN, tok: -1, node: -1, ty: t, has_val: false, val: 0.0, phys: -1 });
        visible.insert(key(name), syms.len - 1);
    }
    for name in BUILTIN_FNS {
        syms.push(Sym { kind: S_BUILTIN, tok: -1, node: -1, ty: any_ty(), has_val: false, val: 0.0, phys: -1 });
        visible.insert(key(name), syms.len - 1);
    }
}

// Each name that the model declares, for `use_before_decl` (spec §3.1).
fn note(tok: i64) {
    if !declared.contains(key(tok_text(tok))) { declared.insert(key(tok_text(tok)), tok); }
}

fn note_declarations(prog: i64) {
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        let k = syntax::nodes[item];
        match k.kind {
            syntax::N_PARAM, syntax::N_LET, syntax::N_FN, syntax::N_TABLE, syntax::N_MATERIAL,
            syntax::N_MESH, syntax::N_POINTS, syntax::N_SYMMETRY, syntax::N_STEP,
            syntax::N_EXPECT => { note(k.tok); },
            syntax::N_PHYSICS => {
                note(k.tok);
                var c = k.first;
                while c >= 0 {
                    let ck = syntax::nodes[c].kind;
                    if ck == syntax::N_UNKNOWN || ck == syntax::N_LET || ck == syntax::N_STATE {
                        note(syntax::nodes[c].tok);
                    }
                    c = syntax::nodes[c].next;
                }
            },
            _ => {},
        }
        item = k.next;
    }
}

// --- helpers -------------------------------------------------------------------------------

fn child(n: i64, i: i64) -> i64 { syntax::child(n, i) }

// Records the dimension of a declaration, for `scales_incomplete`.
fn use_dim(t: Ty) {
    if t.kind != K_NUM || is_dimless(t.dim) { return; }
    for d in used_dims {
        if d == t.dim { return; }
    }
    used_dims.push(t.dim);
}

// The type of a region expression: a geometry value.
fn region_ty(n: i64) -> Ty {
    let t = expr_ty(n);
    if !is_any(t) && t.kind != K_GEOM {
        error(first_tok(n), "type_mismatch", str("the region is ", ty_str(t)),
              "give a Solid, Surface, Curve or Point");
        return any_ty();
    }
    return t;
}

// The type of a setting value: a word, an expression, a list or a range.
fn setting_ty(n: i64) -> Ty {
    let saved = setting_mode;
    setting_mode = true;
    let t = expr_ty(n);
    setting_mode = saved;
    return t;
}

// The value of a setting with name `name` in block `b`, or -1.
fn setting(b: i64, name: u8[:]) -> i64 {
    if b < 0 { return -1; }
    var c = syntax::nodes[b].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_SETTING && tok_text(syntax::nodes[c].tok) == name {
            return syntax::nodes[c].first;
        }
        c = syntax::nodes[c].next;
    }
    return -1;
}

// Checks the setting names of block `b`: a name outside `valid` is the
// error `id`. Then checks the values.
fn check_block(b: i64, valid: const u8[:][:], id: u8[:], what: u8[:]) {
    if b < 0 { return; }
    var c = syntax::nodes[b].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        if k.kind == syntax::N_SETTING && !contains_name(valid, tok_text(k.tok)) {
            error(k.tok, id, str("'", tok_text(k.tok), "' is not a setting of ", what),
                  str("the settings of ", what, ": ", join(valid)));
        }
        c = k.next;
    }
    check_values(b);
}

// Checks the values of the settings of block `b`.
fn check_values(b: i64) {
    if b < 0 { return; }
    var c = syntax::nodes[b].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        if k.kind == syntax::N_SETTING {
            setting_ty(k.first);
        } else if k.kind == syntax::N_REFINE || k.kind == syntax::N_LAYERS {
            region_ty(k.first);
            setting_ty(syntax::nodes[k.first].next);
        }
        c = k.next;
    }
}

fn join(names: const u8[:][:]) -> u8[>..] {
    var out: u8[>..] = [];
    for n, i in names {
        if i > 0 { out.append(", "); }
        out.append(n);
    }
    return out;
}

// Reports `not_constant` for a varying value (spec §7.6).
fn want_constant(n: i64, t: Ty, what: u8[:]) {
    if t.varying {
        error(first_tok(n), "not_constant", str(what, " must be a constant, but it varies in space or time"),
              "use a value that does not depend on x, y, z, t or a field");
    }
}

// --- the header (spec §3, §6.5) ------------------------------------------------------------

fn check_header(h: i64) {
    header_tok = syntax::nodes[h].tok;
    let sp = syntax::child_of_kind(h, syntax::N_SPACE);
    if sp >= 0 {
        let s = tok_text(syntax::nodes[sp].tok);
        if s == "1d" { space = 1; }
        if s == "2d" { space = 2; }
        if s == "axi" { space = 2; axi = true; thickness = true; }
    }
    let us = syntax::child_of_kind(h, syntax::N_UNITSYS);
    if us >= 0 && tok_text(syntax::nodes[us].tok) == "nondim" {
        nondim = true;
        set_nondim(true);
    }
}

// --- param, let, fn, table, material (spec §5, §7.4, §10.2, §20.1) ---------------------------

fn check_value_decl(n: i64, kind: i64) -> i64 {
    let k = syntax::nodes[n];
    let decl = type_of(k.first, true);
    let value = syntax::nodes[k.first].next;
    let val = expr_ty(value);
    assign(first_tok(value), decl, val, str("'", tok_text(k.tok), "'"));
    var ty = decl;
    if is_any(ty) || (ty.kind == K_RESULTS && val.kind == K_RESULTS) { ty = val; }
    ty.varying = val.varying;
    ty.list = val.list;
    use_dim(decl);
    let s = declare(kind, k.tok, n, ty);
    let v, ok = const_value(value);
    if ok && !val.varying {
        syms[s].has_val = true;
        syms[s].val = v;
    }
    if kind == S_PARAM {
        want_constant(value, val, str("the param '", tok_text(k.tok), "'"));
        var c = syntax::nodes[value].next;
        while c >= 0 {
            if syntax::nodes[c].kind == syntax::N_TOL {
                assign(first_tok(syntax::nodes[c].first), decl, expr_ty(syntax::nodes[c].first), "the tolerance");
            }
            if syntax::nodes[c].kind == syntax::N_RANGE {
                let lo = syntax::nodes[c].first;
                let hi = syntax::nodes[lo].next;
                assign(first_tok(lo), decl, expr_ty(lo), "the low end of the range");
                assign(first_tok(hi), decl, expr_ty(hi), "the high end of the range");
                if ok { in_range(first_tok(value), v, c, tok_text(k.tok)); }
            }
            c = syntax::nodes[c].next;
        }
    }
    return s;
}

// `param_range`: the value `v` (SI) is outside the range node `r`.
fn in_range(tok: i64, v: f64, r: i64, name: u8[:]) {
    let lo, oklo = const_value(syntax::nodes[r].first);
    let hi, okhi = const_value(syntax::nodes[syntax::nodes[r].first].next);
    if (oklo && v < lo) || (okhi && v > hi) {
        error(tok, "param_range", str("the value of '", name, "' is outside its range"),
              str("give a value in ", node_text(r)));
    }
}

// The range node of param symbol `s`, or -1.
fn param_range_node(s: i64) -> i64 {
    var c = syntax::nodes[syntax::nodes[syms[s].node].first].next;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_RANGE { return c; }
        c = syntax::nodes[c].next;
    }
    return -1;
}

fn check_fn(n: i64) {
    let k = syntax::nodes[n];
    let mark = locals.len;
    var c = k.first;
    while c >= 0 && syntax::nodes[c].kind == syntax::N_FNPARAM {
        locals.push(Local { name: key(tok_text(syntax::nodes[c].tok)), ty: type_of(syntax::nodes[c].first, true) });
        c = syntax::nodes[c].next;
    }
    let result = if c >= 0 { type_of(c, true) } else { any_ty() };
    let body = if c >= 0 { syntax::nodes[c].next } else { -1 };
    if body >= 0 {
        let t = expr_ty(body);
        assign(first_tok(body), result, t, str("the result of '", tok_text(k.tok), "'"));
    }
    locals.resize(mark);
    declare(S_FN, k.tok, n, result);
}

fn check_table(n: i64) {
    let k = syntax::nodes[n];
    var cols = 0;
    var names: u8[>..] = [];
    var c = k.first;
    while c >= 0 && syntax::nodes[c].kind == syntax::N_COLUMN {
        let ct = type_of(syntax::nodes[c].first, true);
        use_dim(ct);
        let u = syntax::child_of_kind(c, syntax::N_UNIT);
        if u >= 0 {
            let ut, s = unit_ty(u);
            if is_num(ut) && is_num(ct) && ut.dim != ct.dim {
                error(syntax::nodes[u].tok, "unit_mismatch",
                      str("the column '", tok_text(syntax::nodes[c].tok), "' is ", ty_str(ct), ", but its unit is ", dim_str(ut.dim)),
                      str("give a unit of ", ty_str(ct)));
            }
        }
        if cols > 0 { names.append(", "); }
        names.append(tok_text(syntax::nodes[c].tok));
        cols++;
        c = syntax::nodes[c].next;
    }
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_ROW {
            let w = syntax::child_count(c);
            if w != cols {
                error(syntax::nodes[c].tok, "table_row_width", str("the row has ", w, " values, but the table has ", cols, " columns"),
                      str("give one value for each column: ", names));
            }
            var v = syntax::nodes[c].first;
            while v >= 0 {
                expr_ty(v);
                v = syntax::nodes[v].next;
            }
        }
        c = syntax::nodes[c].next;
    }
    declare(S_TABLE, k.tok, n, kind_ty(K_DECL));
}

fn check_material(n: i64) {
    let k = syntax::nodes[n];
    let m = declare(S_MATERIAL, k.tok, n, kind_ty(K_DECL));
    let mark = locals.len;
    var c = k.first;
    while c >= 0 {
        let p = syntax::nodes[c];
        if p.kind == syntax::N_PROP {
            let name = tok_text(p.tok);
            for l in locals[mark..] {
                if l.name == name {
                    error(p.tok, "duplicate_name", str("the material gives '", name, "' twice"), "remove one of the two");
                }
            }
            let decl = type_of(p.first, true);
            let value = syntax::nodes[p.first].next;
            assign(first_tok(value), decl, expr_ty(value), str("the property '", name, "'"));
            use_dim(decl);
            locals.push(Local { name: key(name), ty: decl });
            props.push(Prop { mat: m, tok: p.tok, ty: decl });
        }
        c = p.next;
    }
    locals.resize(mark);
}

// --- mesh and points (spec §9) -------------------------------------------------------------

fn check_discretization(n: i64, kind: i64) {
    let k = syntax::nodes[n];
    region_ty(k.first);
    let b = syntax::nodes[k.first].next;
    check_values(b);
    let size = setting(b, "size");
    if size >= 0 { want_constant(size, setting_ty(size), "the size"); }
    declare(kind, k.tok, n, kind_ty(if kind == S_MESH { K_MESH } else { K_POINTS }));
}
