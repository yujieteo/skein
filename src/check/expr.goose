// The types of expressions (spec §5 to §8): names, literals, operators,
// calls of built-in and user functions, members and `at`.
import std;
import dictionary;
import math;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .builtins;
import .typeexpr;

// The names of Appendix A, the geometry functions of spec §8.1 and the
// functions of the phase 1 gaps G6, G7 and G12. A declaration may not use
// them [`duplicate_name`].
let BUILTIN_FNS = [
    "sym", "skew", "tr", "det", "inv", "transpose", "dev", "dot", "cross", "outer", "norm",
    "ddot", "strain", "defgrad", "hooke", "elasticity", "interp", "sum", "von_mises",
    "principal", "fourier", "abs", "sqrt", "exp", "ln", "sin", "cos", "tan", "tanh", "atan",
    "atan2", "min", "max", "diff", "reaction", "integral", "mean", "results", "boundary",
    "grad", "div", "curl", "lap", "ddt", "ddt2", "re", "imag", "arg", "conj", "eigenvalue",
    "box", "cylinder", "sphere", "rect", "circle", "polygon", "fill", "extrude",
    "revolve", "sweep", "line", "arc", "spline", "point", "shell", "loft", "pattern",
    "pattern_polar", "fillet", "move", "rotate", "scale",
    // Shells and beams (spec §4.3), argmax and when (§14.5), rigid and moment (§12.1).
    "tangent", "ply_stiffness", "laminate", "along", "argmax", "when", "rigid", "moment",
    // rms of a harmonic step (§14.5, gap G2), radiation of a natural bc (§12.1),
    // the wall distance of a turbulence model (§4.2).
    "rms", "radiation", "wall_distance",
];
let BUILTIN_CONSTS = ["I", "pi", "x", "y", "z", "t", "n", "tau", "im", "omega"];

fn tok_text(t: i64) -> const u8[:] { syntax::text(t) }
fn node_tok(n: i64) -> i64 { syntax::nodes[n].tok }
fn node_kind(n: i64) -> i64 { syntax::nodes[n].kind }
fn first_child(n: i64) -> i64 { syntax::nodes[n].first }
fn next_sib(n: i64) -> i64 { syntax::nodes[n].next }

// The first token of the expression `n`, where an error about it goes.
recursive fn first_tok(n: i64) -> i64 {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_BINARY || k.kind == syntax::N_CALL || k.kind == syntax::N_MEMBER || k.kind == syntax::N_AT {
        if k.first >= 0 { return first_tok(k.first); }
    }
    return k.tok;
}

// --- names ---------------------------------------------------------------------------

// The edit distance of two names, for the suggested fix of a misspelt name.
fn edit_distance(a: u8[:], b: u8[:]) -> i64 {
    guard a.len < 64 && b.len < 64 else { return 99; }
    var prev: i64[65] = [0; 65];
    var cur: i64[65] = [0; 65];
    for j in b.len + 1 { prev[j] = j; }
    for i in 1..a.len + 1 {
        cur[0] = i;
        for j in 1..b.len + 1 {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            cur[j] = min(min(prev[j] + 1, cur[j - 1] + 1), prev[j - 1] + cost);
        }
        for j in b.len + 1 { prev[j] = cur[j]; }
    }
    return prev[b.len];
}

// The visible name nearest to `name`, or "".
fn nearest_name(name: u8[:]) -> u8[>..] {
    var best: u8[>..] = [];
    var best_d = max(2, name.len / 3) + 1;
    visible.each() { k, v =>
        let d = edit_distance(name, k);
        if d < best_d || (d == best_d && best.len > 0 && compare(k, best) < 0) {
            best_d = d;
            best.clear();
            best.append(k);
        }
    };
    for p in props {
        let k = tok_text(p.tok);
        let d = edit_distance(name, k);
        if d < best_d { best_d = d; best.clear(); best.append(k); }
    }
    return best;
}

fn unknown_name(tok: i64) {
    let name = tok_text(tok);
    if declared.contains(key(name)) {
        let at = declared.get(key(name));
        var line = 0;
        if at { line = syntax::tokens[at].line; }
        error(tok, "use_before_decl", str("'", name, "' is used above its declaration on line ", line),
              str("move the declaration of '", name, "' above this line"));
        return;
    }
    let near = nearest_name(name);
    if near.len > 0 {
        error(tok, "unknown_name", str("'", name, "' is not declared"), str("did you mean '", near, "'?"));
    } else {
        error(tok, "unknown_name", str("'", name, "' is not declared"), "declare it with let, param or fn above this line");
    }
}

// The materials of physics `p` that give the property `name`: the type of
// the first, and the number that give it and that lack it.
fn prop_in_physics(p: i64, name: u8[:]) -> Ty, i64, i64 {
    var ty = any_ty();
    var have = 0;
    var lack = 0;
    for pm in phys_mats {
        if pm.phys != p { continue; }
        var found = false;
        for pr in props {
            if pr.mat == pm.mat && tok_text(pr.tok) == name {
                if have == 0 { ty = pr.ty; }
                found = true;
                break;
            }
        }
        if found { have++; } else { lack++; }
    }
    return ty, have, lack;
}

fn prop_reported(p: i64, name: u8[:]) -> bool {
    for s in prop_seen {
        if s.phys == p && s.name == name { return true; }
    }
    prop_seen.push(PropSeen { phys: p, name: key(name) });
    return false;
}

// The type of the name at token `tok`.
fn name_ty(tok: i64) -> Ty {
    let name = tok_text(tok);
    for i in 0..locals.len {
        let l = locals[locals.len - 1 - i];
        if l.name == name { return l.ty; }
    }
    if cur_phys >= 0 {
        let ty, have, lack = prop_in_physics(cur_phys, name);
        if have > 0 {
            if lack > 0 && !prop_reported(cur_phys, name) {
                error(tok, "material_prop_missing",
                      str("the physics uses the property '", name, "', but a material of it does not give it"),
                      str("add '", name, "' to each material of the physics"));
            }
            check_prop_types(tok, cur_phys, name);
            return ty;
        }
    }
    let s = visible.get(key(name));
    if s { return syms[s].ty; }
    for p in props {
        if tok_text(p.tok) != name { continue; }
        if cur_phys >= 0 {
            if !prop_reported(cur_phys, name) {
                error(tok, "material_prop_missing",
                      str("the physics uses the property '", name, "', but no material of it gives it"),
                      str("add '", name, "' to each material of the physics"));
            }
            return p.ty;
        }
        // A property outside a physics: in a `bc`, `output`, `report` or `expect`.
        if cur_step >= 0 || in_expect || in_bc { return p.ty; }
        break;
    }
    if setting_mode { return kind_ty(K_WORD); }
    unknown_name(tok);
    return any_ty();
}

// `material_prop_type`: two materials of physics `p` give `name` with two types.
fn check_prop_types(tok: i64, p: i64, name: u8[:]) {
    var first = any_ty();
    var have = false;
    for pm in phys_mats {
        if pm.phys != p { continue; }
        for pr in props {
            if pr.mat != pm.mat || tok_text(pr.tok) != name { continue; }
            if !have {
                first = pr.ty;
                have = true;
            } else if pr.ty.dim != first.dim || pr.ty.rank != first.rank || pr.ty.num != first.num {
                if !prop_reported(-2 - p, name) {
                    error(tok, "material_prop_type",
                          str("two materials give '", name, "' as ", ty_str(first), " and ", ty_str(pr.ty)),
                          "give the property one type in each material");
                }
                return;
            }
        }
    }
}

// --- literals --------------------------------------------------------------------------

// The value of a number token.
fn number_value(tok: i64) -> f64 {
    let v, ok = parse_flt(tok_text(tok));
    return if ok { v } else { 0.0 };
}

fn is_int_text(s: u8[:]) -> bool {
    for c in s {
        if c == '.' || c == 'e' || c == 'E' { return false; }
    }
    return true;
}

// The dimension and the SI factor of the unit expression `n` (spec §6.3).
// `affine` counts the `degC` names in it.
recursive fn unit_value(n: i64, affine: i64&) -> i64[7], f64 {
    let k = syntax::nodes[n];
    match k.kind {
        syntax::N_UNIT => { return unit_value(k.first, affine); },
        syntax::N_UNIT_NAME => {
            let u, ok = find_unit(tok_text(k.tok));
            if !ok {
                error(k.tok, "unknown_unit", str("'", tok_text(k.tok), "' is not a unit"),
                      "use a unit of spec §6.3, for example mm, kN or MPa");
                return NO_DIM, 1.0;
            }
            if u.affine { affine++; }
            return u.dim, u.scale;
        },
        syntax::N_UNIT_MUL, syntax::N_UNIT_DIV => {
            let da, sa = unit_value(k.first, affine);
            let db, sb = unit_value(syntax::nodes[k.first].next, affine);
            if k.kind == syntax::N_UNIT_MUL { return da + db, sa * sb; }
            return da - db, sa / sb;
        },
        syntax::N_UNIT_POW => {
            let d, s = unit_value(k.first, affine);
            // The base is a unit name with whole exponents, so each half power is valid.
            let r, _ok = dim_pow_halves(d, k.aux);
            var f = 1.0;
            for i in abs(k.aux) / 2 { f *= s; }
            if k.aux % 2 != 0 { f *= sqrt(s); }
            if k.aux < 0 { f = 1.0 / f; }
            return r, f;
        },
        _ => { return NO_DIM, 1.0; },
    }
}

// The type of the unit `u` (an N_UNIT node) after a literal, and its SI factor.
fn unit_ty(u: i64) -> Ty, f64 {
    var affine = 0;
    let before = diags.len;
    let d, s = unit_value(u, affine);
    // After an unknown unit the type is unknown, so no error follows from it.
    if diags.len > before { return any_ty(), s; }
    let inner = syntax::nodes[u].first;
    if affine > 0 && syntax::nodes[inner].kind != syntax::N_UNIT_NAME {
        error(syntax::nodes[u].tok, "unit_affine", "degC is valid only alone after a literal (spec §6.3)",
              "write temperature differences and coefficients in K, for example W/(m*K)");
    }
    if nondim {
        error(syntax::nodes[u].tok, "nondim_unit", "a literal has a unit in a nondim model",
              "remove the unit");
    }
    var t = num_ty(d, 0);
    if syntax::nodes[inner].kind == syntax::N_UNIT_NAME {
        let name = tok_text(syntax::nodes[inner].tok);
        if name == "rad" || name == "deg" { t.base = B_ANGLE; }
    }
    return t, s;
}

// The type of a number node.
fn number_ty(n: i64) -> Ty {
    let k = syntax::nodes[n];
    let u = syntax::child_of_kind(n, syntax::N_UNIT);
    if u >= 0 {
        let t, s = unit_ty(u);
        return t;
    }
    let text = tok_text(k.tok);
    var t = base_ty(if is_int_text(text) { B_INT } else { B_REAL });
    t.zero = number_value(k.tok) == 0.0 && text == "0";
    return t;
}

// --- constant values ------------------------------------------------------------------

// The value of the constant expression `n` in SI, and true, or false when it
// is not a constant that the checker can compute.
recursive fn const_value(n: i64) -> f64, bool {
    let k = syntax::nodes[n];
    match k.kind {
        syntax::N_NUMBER => {
            var v = number_value(k.tok);
            let u = syntax::child_of_kind(n, syntax::N_UNIT);
            if u >= 0 {
                var affine = 0;
                let d, s = unit_value(u, affine);
                if affine > 0 {
                    if k.aux >= 0 { v = -v; }
                    return v + 273.15, true;
                }
                v *= s;
            }
            return v, true;
        },
        syntax::N_PERCENT => {
            let t = tok_text(k.tok);
            let v, ok = parse_flt(t[..t.len - 1]);
            return v / 100.0, ok;
        },
        syntax::N_NAME => {
            let name = tok_text(k.tok);
            if name == "pi" { return 3.141592653589793, true; }
            for i in 0..locals.len {
                if locals[locals.len - 1 - i].name == name { return 0.0, false; }
            }
            let s = visible.get(key(name));
            if s && syms[s].has_val { return syms[s].val, true; }
            return 0.0, false;
        },
        syntax::N_UNARY => {
            let v, ok = const_value(k.first);
            return -v, ok && syntax::tokens[k.tok].kind == syntax::TK_MINUS;
        },
        syntax::N_BINARY => {
            let a, oka = const_value(k.first);
            guard oka else { return 0.0, false; }
            let b, okb = const_value(syntax::nodes[k.first].next);
            guard okb else { return 0.0, false; }
            match syntax::tokens[k.tok].kind {
                syntax::TK_PLUS => { return a + b, true; },
                syntax::TK_MINUS => { return a - b, true; },
                syntax::TK_STAR => { return a * b, true; },
                syntax::TK_SLASH => { return a / b, b != 0.0; },
                syntax::TK_CARET => { return pow(a, b), true; },
                _ => { return 0.0, false; },
            }
        },
        _ => { return 0.0, false; },
    }
}

// --- the rules of the operators (spec §6.1, §7.2) ----------------------------------------

fn op_text(tok: i64) -> const u8[:] { tok_text(tok) }

// `a + b`, `a - b` and the two sides of an equation or a comparison: the
// same rank and the same dimension. `what` names the place in a message.
fn sum_ty(tok: i64, a: Ty, b: Ty, what: u8[:]) -> Ty {
    if is_any(a) || is_any(b) { return joined(any_ty(), a, b); }
    if !is_num(a) || !is_num(b) {
        error(tok, "type_mismatch", str(ty_str(a), " ", what, " ", ty_str(b)),
              "both sides must be quantities");
        return joined(any_ty(), a, b);
    }
    if a.zero { return joined(b, a, b); }
    if b.zero { return joined(a, a, b); }
    if a.rank != b.rank {
        error(tok, "rank_mismatch", str(ty_str(a), " ", what, " ", ty_str(b)),
              str("the left side is ", rank_name(a.rank), " and the right side is ", rank_name(b.rank)));
        return joined(any_ty(), a, b);
    }
    if a.dim != b.dim {
        let la = dim_str(a.dim);
        let lb = dim_str(b.dim);
        let q = dim_str(a.dim - b.dim);
        error(tok, "unit_mismatch", str(ty_str(a), " ", what, " ", ty_str(b)),
              str("the dimensions are ", la, " and ", lb, ". The right side needs a factor of ", q));
        return joined(any_ty(), a, b);
    }
    var r = a;
    r.num = max(a.num, b.num);
    r.base = if a.base == B_INT && b.base == B_INT { B_INT } else if a.base == B_ANGLE || b.base == B_ANGLE { B_ANGLE } else { B_REAL };
    r.zero = false;
    r.percent = false;
    return joined(r, a, b);
}

fn product_base(a: Ty, b: Ty) -> i64 {
    if a.base == B_INT && b.base == B_INT { return B_INT; }
    return B_REAL;
}

// `a*b` (spec §7.2).
fn mul_ty(tok: i64, a: Ty, b: Ty) -> Ty {
    if is_any(a) || is_any(b) { return joined(any_ty(), a, b); }
    if !is_num(a) || !is_num(b) {
        error(tok, "type_mismatch", str(ty_str(a), " * ", ty_str(b)), "both factors must be quantities");
        return joined(any_ty(), a, b);
    }
    var rank = 0;
    if a.rank == 0 {
        rank = b.rank;
    } else if b.rank == 0 {
        rank = a.rank;
    } else if a.rank == 1 && b.rank == 1 {
        error(tok, "ambiguous_product", "'*' between two vectors (spec §7.2)", "use dot, outer or cross");
        return joined(any_ty(), a, b);
    } else if a.rank == 2 && (b.rank == 1 || b.rank == 2) {
        rank = b.rank;
    } else if a.rank == 4 && b.rank == 2 {
        rank = 2;
    } else {
        error(tok, "rank_mismatch", str(rank_name(a.rank), " times ", rank_name(b.rank), " is not defined (spec §7.2)"),
              "use dot, ddot, outer or transpose");
        return joined(any_ty(), a, b);
    }
    var r = num_ty(a.dim + b.dim, rank);
    r.num = max(a.num, b.num);
    r.base = product_base(a, b);
    return joined(r, a, b);
}

// `a/b`.
fn div_ty(tok: i64, a: Ty, b: Ty) -> Ty {
    if is_any(a) || is_any(b) { return joined(any_ty(), a, b); }
    if !is_num(a) || !is_num(b) {
        error(tok, "type_mismatch", str(ty_str(a), " / ", ty_str(b)), "both sides must be quantities");
        return joined(any_ty(), a, b);
    }
    if b.rank != 0 {
        error(tok, "rank_mismatch", str("the divisor is ", rank_name(b.rank)), "divide by a scalar only");
        return joined(any_ty(), a, b);
    }
    var r = num_ty(a.dim - b.dim, a.rank);
    r.num = max(a.num, b.num);
    return joined(r, a, b);
}

// `a^e` with the exponent node `e` of type `b` (spec §6.1). A dimensionless
// base takes any constant exponent (gap G30).
fn pow_ty(tok: i64, a: Ty, b: Ty, e: i64) -> Ty {
    if is_any(a) || is_any(b) { return joined(any_ty(), a, b); }
    if !is_num(a) || !is_num(b) || a.rank != 0 || b.rank != 0 {
        error(tok, "type_mismatch", str(ty_str(a), " ^ ", ty_str(b)), "the base and the exponent must be scalars");
        return joined(any_ty(), a, b);
    }
    if !is_dimless(b.dim) {
        error(tok, "unit_mismatch", str("the exponent has the dimension ", dim_str(b.dim)), "the exponent must be a Real");
        return joined(any_ty(), a, b);
    }
    if is_dimless(a.dim) {
        var r = a;
        r.base = if a.base == B_INT && b.base == B_INT { B_INT } else { B_REAL };
        r.zero = false;
        return joined(r, a, b);
    }
    let v, ok = const_value(e);
    let whole = ok && 2.0 * v == floor(2.0 * v) && abs(v) < 64.0;
    let d, halves = dim_pow_halves(a.dim, if whole { (2.0 * v) as i64 } else { 0 });
    if !whole || !halves {
        error(tok, "unit_fractional",
              str("a quantity of dimension ", dim_str(a.dim), " needs a constant exponent that gives exponents in halves"),
              "use a constant exponent such as 2 or 0.5, or divide the base by a reference value");
        return joined(any_ty(), a, b);
    }
    var r = num_ty(d, 0);
    r.num = a.num;
    return joined(r, a, b);
}

// `a == b` and the other comparisons.
fn compare_ty(tok: i64, a: Ty, b: Ty) -> Ty {
    let s = sum_ty(tok, a, b, op_text(tok));
    var r = kind_ty(K_BOOL);
    return joined(r, a, s);
}

fn bool_ty(tok: i64, a: Ty, b: Ty) -> Ty {
    for t in [a, b] {
        if !is_any(t) && t.kind != K_BOOL {
            error(tok, "type_mismatch", str("'", op_text(tok), "' needs Bool values, not ", ty_str(t)),
                  "compare the values with == or < first");
        }
    }
    return joined(kind_ty(K_BOOL), a, b);
}

// `a | b`, `a & b` and `a - b` of two geometry values (spec §8.1).
fn geom_op_ty(tok: i64, a: Ty, b: Ty) -> Ty {
    if is_any(a) || is_any(b) { return joined(any_ty(), a, b); }
    if a.kind == K_BOOL && b.kind == K_BOOL {
        error(tok, "type_mismatch", str("'", op_text(tok), "' joins geometry values. Booleans take '",
                                        op_text(tok), op_text(tok), "'"),
              str("write '", op_text(tok), op_text(tok), "'"));
        return joined(kind_ty(K_BOOL), a, b);
    }
    if a.kind != K_GEOM || b.kind != K_GEOM {
        error(tok, "type_mismatch", str(ty_str(a), " ", op_text(tok), " ", ty_str(b)),
              "both sides must be geometry values");
        return joined(any_ty(), a, b);
    }
    var r = geom_ty(max(a.geo, b.geo));
    if syntax::tokens[tok].kind == syntax::TK_MINUS { r.geo = a.geo; }
    return joined(r, a, b);
}

// --- expressions ----------------------------------------------------------------------

fn is_builtin_fn(name: u8[:]) -> bool { contains_name(BUILTIN_FNS, name) }

// The type of a vector literal: all components of one dimension (spec §6.4).
fn vector_join(tok: i64, acc: Ty, c: Ty, i: i64) -> Ty {
    if is_any(acc) || is_any(c) { return joined(any_ty(), acc, c); }
    if !is_num(c) {
        error(tok, "type_mismatch", str("component ", i + 1, " of a vector is ", ty_str(c)),
              "each component must be a quantity");
        return any_ty();
    }
    if i == 0 { return c; }
    if c.rank != acc.rank {
        error(tok, "rank_mismatch", str("component ", i + 1, " of a vector is ", rank_name(c.rank)),
              "give each component the same rank");
        return any_ty();
    }
    if acc.zero { return joined(c, acc, c); }
    if c.zero { return joined(acc, acc, c); }
    if c.dim != acc.dim {
        error(tok, "unit_mismatch",
              str("component ", i + 1, " of a vector has the dimension ", dim_str(c.dim), ", not ", dim_str(acc.dim)),
              "give all components one dimension (spec §6.4)");
        return any_ty();
    }
    var r = joined(acc, acc, c);
    r.num = max(acc.num, c.num);
    if c.base != acc.base { r.base = B_REAL; }
    return r;
}

fn vector_size(tok: i64, count: i64, rank: i64) {
    if rank == 0 && count != space {
        error(tok, "type_mismatch", str("a vector of ", count, " components in a ", space, "D model"),
              str("give ", space, " components"));
    }
}

// The member `name` of a value of type `t` at token `tok` (spec §5.2, §5.6, §15.1).
fn member_ty(tok: i64, obj: i64, t: Ty) -> Ty {
    let name = tok_text(tok);
    if is_any(t) { return t; }
    if syntax::nodes[obj].kind == syntax::N_NAME {
        let s = visible.get(key(tok_text(syntax::nodes[obj].tok)));
        if s && syms[s].kind == S_TABLE {
            var c = syntax::nodes[syms[s].node].first;
            while c >= 0 {
                if syntax::nodes[c].kind == syntax::N_COLUMN && tok_text(syntax::nodes[c].tok) == name {
                    var ct = type_of(syntax::nodes[c].first, false);
                    ct.list = true;
                    ct.table = s;
                    return ct;
                }
                c = syntax::nodes[c].next;
            }
            error(tok, "unknown_name", str("the table has no column '", name, "'"), "use a column of the table");
            return any_ty();
        }
    }
    if t.kind == K_RESULTS {
        for f in result_fields {
            if f.res == t.aux && f.name == name { return f.ty; }
        }
        error(tok, "unknown_name", str("the results have no field '", name, "'"), "use an unknown of the other model");
        return any_ty();
    }
    if t.kind == K_MESH && (name == "size") { return num_ty(LENGTH, 0); }
    if t.kind == K_NUM && t.rank > 0 {
        let comps = if t.rank == 1 { name.len == 1 } else { name.len == 2 };
        var ok = comps && t.rank <= 2;
        for c in name {
            let limit = if axi && t.rank == 2 { "xyt" } else if space == 3 { "xyz" } else if space == 2 { "xy" } else { "x" };
            if !contains(limit, c) { ok = false; }
        }
        if ok {
            var r = t;
            r.rank = 0;
            return r;
        }
    }
    error(tok, "unknown_name", str(ty_str(t), " has no member '", name, "'"),
          "use .x, .y or .z of a vector, or .xx, .xy, … of a Tensor2");
    return any_ty();
}

// The type of a call of a user function (spec §7.4).
fn user_call_ty(tok: i64, s: i64, a: Args) -> Ty {
    let f = syms[s].node;
    var params = 0;
    var c = syntax::nodes[f].first;
    while c >= 0 && syntax::nodes[c].kind == syntax::N_FNPARAM {
        if params < a.count && params < MAX_ARGS {
            let pt = type_of(syntax::nodes[c].first, false);
            assign(first_tok(a.node[params]), pt, a.ty[params], str("argument ", params + 1, " of ", tok_text(syntax::nodes[f].tok)));
        }
        params++;
        c = syntax::nodes[c].next;
    }
    if params != a.count {
        error(tok, "type_mismatch", str(tok_text(syntax::nodes[f].tok), " takes ", params, " arguments, not ", a.count),
              "give one argument for each parameter");
    }
    var r = syms[s].ty;
    r.varying = args_varying(a);
    for i in min(a.count, MAX_ARGS) { r.list = r.list || a.ty[i].list; r.table = rows(r.table, a.ty[i].table); }
    return r;
}

// The type of expression `n`.
recursive fn expr_ty(n: i64) -> Ty {
    let k = syntax::nodes[n];
    match k.kind {
        syntax::N_NUMBER => { return number_ty(n); },
        syntax::N_PERCENT => {
            var t = base_ty(B_REAL);
            t.percent = true;
            return t;
        },
        syntax::N_STRING => { return kind_ty(K_STR); },
        syntax::N_BOOL => { return kind_ty(K_BOOL); },
        syntax::N_NAME => { return name_ty(k.tok); },
        syntax::N_VECTOR => {
            let unit = syntax::child_of_kind(n, syntax::N_UNIT);
            var ut = any_ty();
            if unit >= 0 {
                let t, f = unit_ty(unit);
                ut = t;
            }
            var acc = any_ty();
            var count = 0;
            var c = k.first;
            while c >= 0 {
                if syntax::nodes[c].kind != syntax::N_UNIT {
                    let ct = if unit >= 0 && is_bare(c) { ut } else { expr_ty(c) };
                    acc = if count == 0 { ct } else { vector_join(first_tok(c), acc, ct, count) };
                    count++;
                }
                c = syntax::nodes[c].next;
            }
            if is_any(acc) { return acc; }
            vector_size(k.tok, count, acc.rank);
            acc.rank = acc.rank + 1;
            acc.zero = false;
            if acc.base == B_INT { acc.base = B_REAL; }
            return acc;
        },
        syntax::N_RECORD => {
            var c = k.first;
            while c >= 0 {
                expr_ty(syntax::nodes[c].first);
                c = syntax::nodes[c].next;
            }
            return kind_ty(K_RECORD);
        },
        syntax::N_LIST => {
            var c = k.first;
            while c >= 0 {
                if syntax::nodes[c].kind != syntax::N_UNIT { expr_ty(c); }
                c = syntax::nodes[c].next;
            }
            return kind_ty(K_LIST);
        },
        syntax::N_RANGE => {
            let a = expr_ty(k.first);
            let b = expr_ty(syntax::nodes[k.first].next);
            return sum_ty(k.tok, a, b, "..");
        },
        syntax::N_UNARY => {
            let a = expr_ty(k.first);
            if syntax::tokens[k.tok].kind == syntax::TK_BANG { return bool_ty(k.tok, a, a); }
            if !is_any(a) && !is_num(a) { return not_number(k.tok, a); }
            var r = a;
            r.zero = false;
            return r;
        },
        syntax::N_BINARY => {
            let op = syntax::tokens[k.tok].kind;
            let a = expr_ty(k.first);
            let rhs = syntax::nodes[k.first].next;
            let b = expr_ty(rhs);
            match op {
                syntax::TK_PLUS => { return sum_ty(k.tok, a, b, "+"); },
                syntax::TK_MINUS => {
                    if a.kind == K_GEOM || b.kind == K_GEOM { return geom_op_ty(k.tok, a, b); }
                    return sum_ty(k.tok, a, b, "-");
                },
                syntax::TK_STAR => { return mul_ty(k.tok, a, b); },
                syntax::TK_SLASH => { return div_ty(k.tok, a, b); },
                syntax::TK_CARET => { return pow_ty(k.tok, a, b, rhs); },
                syntax::TK_BAR, syntax::TK_AMP => { return geom_op_ty(k.tok, a, b); },
                syntax::TK_OROR, syntax::TK_ANDAND => { return bool_ty(k.tok, a, b); },
                _ => { return compare_ty(k.tok, a, b); },
            }
        },
        syntax::N_CALL => {
            let callee = k.first;
            var a = Args { count: 0, ty: [any_ty(); 8], node: [-1; 8], name: [-1; 8] };
            let saved = setting_mode;
            // The profile of extrude, revolve, sweep and loft is a 2D Surface
            // (spec §8), and loft has 2.
            let cn = call_name(callee);
            let profile = space == 3 && (cn == "extrude" || cn == "revolve" || cn == "sweep" || cn == "loft");
            let profiles = if cn == "loft" { 2 } else { 1 };
            var c = syntax::nodes[callee].next;
            while c >= 0 {
                let v = syntax::nodes[c].first;
                if profile && a.count < profiles { space = 2; }
                let t = expr_ty(v);
                if profile { space = 3; }
                if a.count < MAX_ARGS {
                    a.ty[a.count] = t;
                    a.node[a.count] = v;
                    a.name[a.count] = syntax::nodes[c].aux;
                } else {
                    extra_point(call_name(callee), t, v, a.count);
                }
                a.count++;
                c = syntax::nodes[c].next;
            }
            setting_mode = saved;
            return call_ty(k.tok, callee, a);
        },
        syntax::N_MEMBER => {
            let t = expr_ty(k.first);
            return member_ty(k.tok, k.first, t);
        },
        syntax::N_AT => {
            let f = expr_ty(k.first);
            let p = expr_ty(syntax::nodes[k.first].next);
            return at_ty(k.tok, f, p);
        },
        _ => { return any_ty(); },
    }
}

fn not_number(tok: i64, a: Ty) -> Ty {
    error(tok, "type_mismatch", str("'-' of ", ty_str(a)), "negate a quantity only");
    return any_ty();
}

fn call_name(callee: i64) -> const u8[:] {
    if syntax::nodes[callee].kind == syntax::N_NAME { return tok_text(syntax::nodes[callee].tok); }
    return "";
}

// The type of a call: of a user function, a built-in function, or a word
// with arguments in a setting, for example `generalized_alpha(0.8)`.
fn call_ty(tok: i64, callee: i64, a: Args) -> Ty {
    let ct = syntax::nodes[callee].tok;
    guard syntax::nodes[callee].kind == syntax::N_NAME else {
        error(ct, "type_mismatch", "only a function can be called", "call a function by its name");
        return any_ty();
    }
    let name = tok_text(ct);
    let s = visible.get(key(name));
    if s && syms[s].kind == S_FN { return user_call_ty(ct, s, a); }
    if s && syms[s].kind != S_BUILTIN {
        error(ct, "type_mismatch", str("'", name, "' is not a function"), "call a function, or remove the brackets");
        return any_ty();
    }
    if is_builtin_fn(name) { return builtin_ty(ct, name, a); }
    if setting_mode { return kind_ty(K_WORD); }
    unknown_name(ct);
    return any_ty();
}

// `f at p` (spec §7.5).
fn at_ty(tok: i64, f: Ty, p: Ty) -> Ty {
    if !is_any(p) && !(p.kind == K_GEOM && p.geo == 0) && !(is_num(p) && p.rank == 1 && (p.dim == LENGTH || p.zero)) {
        error(tok, "type_mismatch", str("the point of 'at' is ", ty_str(p)), "give a Point or a vector of Length");
    }
    var r = f;
    r.varying = false;
    r.list = f.list || p.list;
    r.table = rows(f.table, p.table);
    return r;
}

// A bare number inside a literal with a unit takes that unit.
fn bare_in_unit(t: Ty, unit: i64) -> Ty {
    if unit < 0 { return t; }
    let ut, s = unit_ty(unit);
    var r = ut;
    r.list = t.list;
    r.table = t.table;
    return r;
}

// True for a number with no unit, or `-` and one.
recursive fn is_bare(n: i64) -> bool {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_NUMBER { return syntax::child_of_kind(n, syntax::N_UNIT) < 0; }
    if k.kind == syntax::N_UNARY { return is_bare(k.first); }
    return false;
}

// --- assignability (spec §5.7) ---------------------------------------------------------

// Checks that a value of type `val` may initialize a declaration of type
// `decl`. `what` names the declaration in a message.
fn assign(tok: i64, decl: Ty, val: Ty, what: u8[:]) -> bool {
    if is_any(decl) || is_any(val) { return true; }
    if decl.kind == K_GEOM && val.kind == K_GEOM {
        if decl.geo != val.geo {
            error(tok, "geometry_dim", str(what, " is a ", GEO_NAMES[decl.geo], ", but the value is a ", GEO_NAMES[val.geo]),
                  str("declare it as ", GEO_NAMES[val.geo]));
            return false;
        }
        return true;
    }
    if decl.kind != val.kind {
        error(tok, "type_mismatch", str(what, " is ", ty_str(decl), ", but the value is ", ty_str(val)),
              str("the value has the type ", ty_str(val)));
        return false;
    }
    if decl.kind != K_NUM { return true; }
    if val.zero { return true; }
    if decl.rank != val.rank {
        error(tok, "type_mismatch", str(what, " is ", ty_str(decl), ", but the value is ", ty_str(val)),
              str("the value has the type ", ty_str(val)));
        return false;
    }
    if val.num > decl.num || (val.num != NS_REAL && val.num != decl.num) {
        error(tok, "type_mismatch", str(what, " is ", ty_str(decl), ", but the value is ", ty_str(val)),
              str("the value has the type ", ty_str(val)));
        return false;
    }
    if decl.dim != val.dim {
        error(tok, "unit_mismatch", str(what, " is ", ty_str(decl), ", but the value is ", ty_str(val)),
              str("the value has the dimension ", dim_str(val.dim), ". It needs a factor of ", dim_str(decl.dim - val.dim)));
        return false;
    }
    if decl.base == B_INT && val.base != B_INT && is_dimless(decl.dim) && decl.rank == 0 {
        error(tok, "type_mismatch", str(what, " is Int, but the value is ", ty_str(val)),
              "write an integer with no '.' and no exponent");
        return false;
    }
    return true;
}
