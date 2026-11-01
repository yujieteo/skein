// The runtime expressions of the shell generator: a value of the model as
// Goose code of the generated program.
//
// A field is a value at a node of the mesh, for example a let of the
// physics. Its code uses the names of its field function in the program:
// p (the point), w and th (the unknowns), tn and te (N and eps_m), and mat
// (the material of the node). A global value is the value of a report or
// an expect: `f at p`, `max(f, r)`, `min(f, r)`, `reaction(r)`,
// `eigenvalue(k)`, and numbers, + - * / ^ and functions of these. A part
// of an expression that is the same at each node is a number of the code.
// `global.f`, with global the results value, is the field f on the model
// of the results value.
//
// Each expression also has a dimension (the halves of length, mass and time
// of codegen::Unit), for the unit of a report.
import std;
import math;
namespace shells_gen;
import syntax from "../syntax";
import codegen.codegen;

struct Dim {
    l: i64,
    m: i64,
    t: i64,
}

fn no_dim() -> Dim { return Dim { l: 0, m: 0, t: 0 }; }
fn dim_add(a: Dim, b: Dim, s: i64) -> Dim { return Dim { l: a.l + s * b.l, m: a.m + s * b.m, t: a.t + s * b.t }; }
fn dim_scale(a: Dim, k: f64) -> Dim {
    return Dim { l: round(a.l * k) as! i64, m: round(a.m * k) as! i64, t: round(a.t * k) as! i64 };
}
fn same_dim(a: Dim, b: Dim) -> bool { return a.l == b.l && a.m == b.m && a.t == b.t; }

let D_LENGTH = Dim { l: 2, m: 0, t: 0 };
let D_FORCE = Dim { l: 2, m: 2, t: -4 };
let D_PRESSURE = Dim { l: -2, m: 2, t: -4 };
let D_FREQUENCY = Dim { l: 0, m: 0, t: -2 };

// The unit of dimension d in the unit system of the header.
fn unit_text(d: Dim) -> u8[>..] {
    let sys = codegen::system;
    if same_dim(d, no_dim()) { return str(""); }
    if same_dim(d, D_LENGTH) { return str(sys.length); }
    if same_dim(d, D_FORCE) { return str(sys.force); }
    if same_dim(d, D_PRESSURE) { return str(sys.pressure); }
    if same_dim(d, dim_add(D_FORCE, D_LENGTH, -1)) { return str(sys.force, "/", sys.length); }
    if same_dim(d, dim_add(D_FORCE, D_LENGTH, 1)) { return str(sys.force, "*", sys.length); }
    if same_dim(d, D_FREQUENCY) { return str("Hz"); }
    return str(sys.length, "^", d.l / 2.0, "*kg^", d.m / 2.0, "*s^", d.t / 2.0);
}

// The dimension of a declared type.
recursive fn type_dim(t: i64) -> Dim {
    let k = kind(t);
    if k == syntax::N_DIM_MUL || k == syntax::N_DIM_DIV {
        return dim_add(type_dim(child(t, 0)), type_dim(child(t, 1)), if k == syntax::N_DIM_MUL { 1 } else { -1 });
    }
    if k != syntax::N_TYPE_NAME { return no_dim(); }
    let w = name(t);
    if w == "Length" { return D_LENGTH; }
    if w == "Area" { return Dim { l: 4, m: 0, t: 0 }; }
    if w == "Force" { return D_FORCE; }
    if w == "Pressure" { return D_PRESSURE; }
    if w == "Density" { return Dim { l: -6, m: 2, t: 0 }; }
    if w == "Mass" { return Dim { l: 0, m: 2, t: 0 }; }
    if w == "Time" { return Dim { l: 0, m: 0, t: 2 }; }
    if w == "Frequency" { return D_FREQUENCY; }
    if w == "Vector" || w == "Tensor2" {
        let a = first(t);
        return if a >= 0 { type_dim(a) } else { no_dim() };
    }
    return no_dim();
}

// The declared type of name `word`, or -1: a let of the physics, a
// property of a shell material, or a let or param item.
fn type_of_name(word: u8[:]) -> i64 {
    let i = find(let_names, word);
    if i >= 0 { return let_types[i]; }
    let p = prop_type(word);
    if p >= 0 { return p; }
    let v = find(value_names, word);
    if v >= 0 { return value_types[v]; }
    return -1;
}

var eigen_modal = false;    // eigenvalue(k) is a frequency

// --- what depends on the node ------------------------------------------------------

// Whether expression v changes from node to node, or from solve to solve.
recursive fn varies(v: i64) -> bool {
    let k = kind(v);
    if k == syntax::N_NAME {
        let w = name(v);
        if bound(w) >= 0 { return false; }
        if is_w(w) || is_theta(w) || w == "N" || w == "eps_m" || w == "Nx" || is_source(v) { return true; }
        if find(let_names, w) >= 0 || prop_type(w) >= 0 { return true; }
        let i = find(value_names, w);
        if i >= 0 {
            if depth > 100 { return false; }
            depth++;
            let r = varies(value_nodes[i]);
            depth--;
            return r;
        }
        return axis_of(w) >= 0 || w == "n";
    }
    if k == syntax::N_AT || is_call(v, "reaction") || is_call(v, "max") || is_call(v, "min") || is_call(v, "eigenvalue") {
        return true;
    }
    var c = first(v);
    while c >= 0 {
        if kind(c) != syntax::N_UNIT && varies(if kind(c) == syntax::N_ARG { first(c) } else { c }) { return true; }
        c = next(c);
    }
    return false;
}

// --- the code ----------------------------------------------------------------------

var node_mode = true;

// The fields of the program: field k is the expression fields[k] at a node,
// as an f64[9] (a scalar in [0], a vector in [0..3]).
var fields: i64[>..] = [];
var field_ranks: i64[>..] = [];

// The part of out from `start`, which it removes from out.
fn cut_tail(out: u8[>..]&, start: i64) -> u8[>..] {
    var t: u8[>..] = [];
    t.append(out[start..]);
    while out.len > start { out.pop(); }
    return t;
}

fn field_index(v: i64) -> i64 {
    for f, i in fields { if f == v { return i; } }
    fields.push(v);
    field_ranks.push(rank_of(v));
    return fields.len - 1;
}

// The rank of the code of emit_rt for v: 1, 3 or 9. emit_rt reports the
// faults.
recursive fn rank_of(v: i64) -> i64 {
    if !varies(v) { return value(v).n; }
    let k = kind(v);
    if k == syntax::N_NAME || source_member(v) {
        let w = name(v);
        if axis_of(w) >= 0 { return 1; }
        if w == "n" || is_w(w) || is_theta(w) { return 3; }
        if w == "N" || w == "eps_m" { return 9; }
        var node = -1;
        if find(let_names, w) >= 0 { node = let_nodes[find(let_names, w)]; }
        if node < 0 && find(value_names, w) >= 0 { node = value_nodes[find(value_names, w)]; }
        if node < 0 || depth > 100 { return 1; }
        depth++;
        let r = rank_of(node);
        depth--;
        return r;
    }
    if k == syntax::N_UNARY { return rank_of(first(v)); }
    if k == syntax::N_MEMBER { return 1; }
    if k == syntax::N_VECTOR { return 3; }
    if k == syntax::N_BINARY {
        let o = op(v);
        let ra = rank_of(child(v, 0));
        if o == "+" || o == "-" || o == "/" { return ra; }
        if o == "^" { return 1; }
        let rb = rank_of(child(v, 1));
        if ra == 9 && rb == 3 { return 3; }
        return if ra == 1 { rb } else { ra };
    }
    if k == syntax::N_AT {
        let f = child(v, 0);
        if kind(f) == syntax::N_NAME && name(f) == "Nx" { return 1; }
        return rank_of(f);
    }
    if is_call(v, "reaction") || is_call(v, "cross") { return 3; }
    return 1;
}

fn lit_val(out: u8[>..]&, x: Val) -> i64 {
    if x.n == 1 {
        lit(out, x.v[0]);
        return 1;
    }
    out.push('[');
    for k in x.n {
        if k > 0 { out.append(", "); }
        lit(out, x.v[k]);
    }
    out.push(']');
    return x.n;
}

let RT_FIX = "write numbers, names, x, y, z, n, the unknowns, N, eps_m, lets of the physics, dot, cross, norm, sqrt, abs, sin, cos and + - * / ^ of these";

// Appends the code of expression v. Gives its rank (1, 3 or 9) and its
// dimension.
recursive fn emit_rt(v: i64, out: u8[>..]&) -> i64, Dim {
    let d = dim_of(v);
    if !varies(v) {
        let x = value(v);
        return lit_val(out, x), d;
    }
    let k = kind(v);
    if k == syntax::N_NAME { return emit_word(name(v), tok(v), out), d; }
    if source_member(v) {
        if !node_mode { unsupported(tok(v), str("'", source_of(v), "' is a value at a node: write it with at"), RT_FIX); }
        return emit_word(name(v), tok(v), out), d;
    }
    if k == syntax::N_UNARY && op(v) == "-" {
        let start = out.len;
        let r, e = emit_rt(first(v), out);
        let a = cut_tail(out, start);
        format(out, "neg", r, "(", a, ")");
        return r, d;
    }
    if k == syntax::N_MEMBER {
        let a = axis_of(name(v));
        out.append("comp(");
        let r, e = emit_rt(first(v), out);
        format(out, ", ", a, ")");
        if r != 3 || a < 0 { unsupported(tok(v), "a member .x, .y or .z of a value that is not a vector", RT_FIX); }
        return 1, d;
    }
    if k == syntax::N_VECTOR {
        let u = syntax::child_of_kind(v, syntax::N_UNIT);
        let f = if u >= 0 { unit_factor(u) } else { 1.0 };
        if f != 1.0 {
            out.append("scale3(");
            lit(out, f);
            out.append(", ");
        }
        out.push('[');
        var count = 0;
        var c = first(v);
        while c >= 0 {
            if kind(c) != syntax::N_UNIT {
                if count > 0 { out.append(", "); }
                let r, e = emit_rt(c, out);
                if r != 1 { unsupported(tok(c), "a component of a vector is a vector", "write a scalar"); }
                count++;
            }
            c = next(c);
        }
        if count != 3 { unsupported(tok(v), "a vector in this slice has 3 components", "write (x, y, z)"); }
        out.append(if f != 1.0 { "])" } else { "]" });
        return 3, d;
    }
    if k == syntax::N_BINARY { return emit_binary(v, out), d; }
    if k == syntax::N_AT || k == syntax::N_CALL { return emit_call(v, out), d; }
    unsupported(tok(v), str("the expression '", source_of(v), "' is not in this slice"), RT_FIX);
    return 1, d;
}

// The code of name `w` at a node.
recursive fn emit_word(w: u8[:], at: i64, out: u8[>..]&) -> i64 {
    if node_mode && axis_of(w) >= 0 {
        format(out, "p[", axis_of(w), "]");
        return 1;
    }
    if node_mode && w == "n" {
        out.append("radial(p)");
        return 3;
    }
    if node_mode && is_w(w) {
        out.append("w");
        return 3;
    }
    if node_mode && is_theta(w) {
        out.append("th");
        return 3;
    }
    if node_mode && (w == "N" || w == "eps_m") {
        out.append(if w == "N" { "tn" } else { "te" });
        return 9;
    }
    if node_mode && prop_type(w) >= 0 && find(let_names, w) < 0 { return emit_prop(w, at, out); }
    var node = -1;
    if node_mode && find(let_names, w) >= 0 { node = let_nodes[find(let_names, w)]; }
    if node < 0 && find(value_names, w) >= 0 { node = value_nodes[find(value_names, w)]; }
    if node >= 0 && depth < 100 {
        depth++;
        out.push('(');
        let r, e = emit_rt(node, out);
        out.push(')');
        depth--;
        return r;
    }
    unsupported(at, str("'", w, "' has no value here"), RT_FIX);
    return 1;
}

// A property of the materials of the shell: the value of the material of
// the node.
fn emit_prop(w: u8[:], at: i64, out: u8[>..]&) -> i64 {
    out.append("pick(mat, [");
    for m in shell_mats.len {
        if m > 0 { out.append(", "); }
        let x = prop_value(m, w, at);
        if x.n != 1 { unsupported(at, "a property of the material at a node is a scalar", RT_FIX); }
        lit(out, x.v[0]);
    }
    out.append("])");
    return 1;
}

recursive fn emit_binary(v: i64, out: u8[>..]&) -> i64 {
    let o = op(v);
    let s0 = out.len;
    let ra, da = emit_rt(child(v, 0), out);
    let s1 = out.len;
    let rb, db = emit_rt(child(v, 1), out);
    let b = cut_tail(out, s1);
    let a = cut_tail(out, s0);
    if (o == "+" || o == "-") && ra == rb {
        if ra == 1 {
            format(out, "(", a, " ", o, " ", b, ")");
        } else {
            format(out, if o == "+" { "add" } else { "sub" }, ra, "(", a, ", ", b, ")");
        }
        return ra;
    }
    if o == "*" && ra == 1 && rb == 1 {
        format(out, "(", a, " * ", b, ")");
        return 1;
    }
    if o == "*" && (ra == 1 || rb == 1) {
        let r = if ra == 1 { rb } else { ra };
        if ra == 1 { format(out, "scale", r, "(", a, ", ", b, ")"); } else { format(out, "scale", r, "(", b, ", ", a, ")"); }
        return r;
    }
    if o == "*" && ra == 9 && rb == 3 {
        format(out, "tdotv(", a, ", ", b, ")");
        return 3;
    }
    if o == "/" && rb == 1 {
        if ra == 1 { format(out, "(", a, " / ", b, ")"); } else { format(out, "scale", ra, "(1.0 / ", b, ", ", a, ")"); }
        return ra;
    }
    if o == "^" && ra == 1 && rb == 1 {
        format(out, "pow(", a, ", ", b, ")");
        return 1;
    }
    unsupported(tok(v), str("the operator ", o, " of these values is not in this slice"), RT_FIX);
    return 1;
}

let SCALAR_FNS = ["sqrt", "abs", "sin", "cos", "tan"];

recursive fn emit_call(v: i64, out: u8[>..]&) -> i64 {
    if kind(v) == syntax::N_AT || is_call(v, "reaction") || is_call(v, "max") || is_call(v, "min") || is_call(v, "eigenvalue") {
        if node_mode {
            unsupported(tok(v), str("'", source_of(v), "' is a value of the whole model, not of a node"), RT_FIX);
            return 1;
        }
        return emit_global(v, out);
    }
    let f = if kind(child(v, 0)) == syntax::N_NAME { name(child(v, 0)) } else { "" };
    let count = syntax::child_count(v) - 1;
    if count == 1 && position(SCALAR_FNS, f) >= 0 {
        format(out, f, "(");
        let r, d = emit_rt(arg(v, 0), out);
        out.push(')');
        if r != 1 { unsupported(tok(v), str(f, " of a value that is not a scalar"), RT_FIX); }
        return 1;
    }
    if count == 2 && (f == "dot" || f == "cross") {
        out.append(if f == "dot" { "dot3(" } else { "cross3(" });
        let ra, da = emit_rt(arg(v, 0), out);
        out.append(", ");
        let rb, db = emit_rt(arg(v, 1), out);
        out.push(')');
        if ra != 3 || rb != 3 { unsupported(tok(v), str(f, " of values that are not vectors"), RT_FIX); }
        return if f == "dot" { 1 } else { 3 };
    }
    if count == 1 && f == "norm" {
        out.append("norm3(");
        let r, d = emit_rt(arg(v, 0), out);
        out.push(')');
        return 1;
    }
    unsupported(tok(v), str("the function '", f, "' is not in this slice"), RT_FIX);
    return 1;
}

// --- the global values ---------------------------------------------------------------

let GLOBAL_FIX = "write f at p, max(f, r), min(f, r), reaction(r) or eigenvalue(k), and numbers and + - * / of these";

// `f at p`, `max(f, r)`, `min(f, r)`, `reaction(r)` and `eigenvalue(k)`.
recursive fn emit_global(v: i64, out: u8[>..]&) -> i64 {
    if kind(v) == syntax::N_AT {
        let f = child(v, 0);
        let at = vec3(value(child(v, 1)));
        if kind(f) == syntax::N_NAME && name(f) == "Nx" && beam_tok >= 0 {
            out.append("axial_at(");
            lits(out, at);
            out.append(")");
            return 1;
        }
        let k = field_index(f);
        out.append(if field_ranks[k] == 1 { "at1(" } else if field_ranks[k] == 3 { "at3(" } else { "(" });
        format(out, "field_at(", k, ", [");
        lits(out, at);
        out.append("], ");
        quoted(out, source_of(v));
        out.append("))");
        return field_ranks[k];
    }
    if (is_call(v, "max") || is_call(v, "min")) && syntax::child_count(v) == 3 {
        let f = arg(v, 0);
        let rv = arg(v, 1);
        let r = find_region(rv);
        let whole = physics_tok >= 0 && surface_shape(rv) == mesh_g;
        let big = is_call(v, "max");
        if kind(f) == syntax::N_NAME && name(f) == "Nx" && beam_tok >= 0 {
            if r < 0 || regions[r].piece_count == 0 {
                unsupported(tok(rv), "max and min of Nx are on a beam Curve", GLOBAL_FIX);
                return 1;
            }
            let rc = region_runs_code(r);
            format(out, "beam_extreme(", rc, ", ", if big { "true" } else { "false" }, ")");
            return 1;
        }
        let k = field_index(f);
        if field_ranks[k] != 1 { unsupported(tok(v), "max and min of a value that is not a scalar", GLOBAL_FIX); }
        format(out, "node_extreme(", k, ", ");
        if whole {
            out.append("surface_nodes()");
        } else if r >= 0 && (regions[r].node >= 0 || regions[r].point) {
            let rc = region_code(r);
            out.append(rc);
        } else {
            unsupported(tok(rv), "the region of max or min is the meshed surface, or a Curve of the shell", GLOBAL_FIX);
        }
        format(out, ", ", if big { "true" } else { "false" }, ")");
        return 1;
    }
    if is_call(v, "reaction") && syntax::child_count(v) == 2 {
        let r = find_region(arg(v, 0));
        if r < 0 || !(regions[r].node >= 0 || regions[r].point) {
            unsupported(tok(v), "reaction needs a Curve of the shell or a Point that a let names", GLOBAL_FIX);
            return 3;
        }
        let rc = region_code(r);
        format(out, "reaction_of(", rc, ")");
        return 3;
    }
    if is_call(v, "eigenvalue") && syntax::child_count(v) == 2 {
        format(out, "shells::eigenvalue(", round(scalar_value(arg(v, 0))) as! i64, ")");
        return 1;
    }
    unsupported(tok(v), str("'", source_of(v), "' is not in this slice"), GLOBAL_FIX);
    return 1;
}

// Appends the code of global expression v. Gives its rank and dimension.
fn emit_global_expr(v: i64, out: u8[>..]&) -> i64, Dim {
    node_mode = false;
    let r, d = emit_rt(v, out);
    node_mode = true;
    return r, d;
}

// The function of field k in the program: the field at node k of view v.
fn emit_field_fn(out: u8[>..]&, k: i64) {
    var code: u8[>..] = [];
    let rank, d = emit_rt(fields[k], code);
    format(out, "// field ", k, ": ", source_of(fields[k]), "\nfn field_", k, "(v: i64, k: i64) -> f64[9] {\n");
    out.append("    let p = v_point(v, k);\n    let w = v_w(v, k);\n    let th = v_th(v, k);\n");
    out.append("    let tn = v_tn(v, k);\n    let te = v_te(v, k);\n    let mat = v_mat(v, k);\n");
    format(out, "    return pack", rank, "(", code, ");\n}\n\n");
}

// --- dimensions ----------------------------------------------------------------------

recursive fn dim_of(v: i64) -> Dim {
    let k = kind(v);
    if k == syntax::N_NUMBER {
        let u = first(v);
        if u < 0 || angle_unit(u) > 0.0 { return no_dim(); }
        let x = codegen::unit_value(u);
        return Dim { l: x.l, m: x.m, t: x.t };
    }
    if k == syntax::N_VECTOR {
        var d = no_dim();
        var c = first(v);
        while c >= 0 {
            if kind(c) == syntax::N_UNIT {
                if angle_unit(c) < 0.0 {
                    let x = codegen::unit_value(c);
                    d = Dim { l: x.l, m: x.m, t: x.t };
                }
            } else if same_dim(d, no_dim()) {
                d = dim_of(c);
            }
            c = next(c);
        }
        return d;
    }
    if k == syntax::N_NAME { return name_dim(name(v)); }
    if source_member(v) { return name_dim(name(v)); }
    if k == syntax::N_UNARY || k == syntax::N_MEMBER { return dim_of(first(v)); }
    if k == syntax::N_AT { return dim_of(child(v, 0)); }
    if k == syntax::N_BINARY {
        let o = op(v);
        let a = dim_of(child(v, 0));
        if o == "*" { return dim_add(a, dim_of(child(v, 1)), 1); }
        if o == "/" { return dim_add(a, dim_of(child(v, 1)), -1); }
        if o == "^" { return dim_scale(a, scalar_value(child(v, 1))); }
        return if same_dim(a, no_dim()) { dim_of(child(v, 1)) } else { a };
    }
    if k == syntax::N_CALL {
        let f = if kind(child(v, 0)) == syntax::N_NAME { name(child(v, 0)) } else { "" };
        if f == "reaction" { return D_FORCE; }
        if f == "eigenvalue" { return if eigen_modal { D_FREQUENCY } else { no_dim() }; }
        if f == "max" || f == "min" || f == "abs" || f == "norm" { return dim_of(arg(v, 0)); }
        if f == "sqrt" { return dim_scale(dim_of(arg(v, 0)), 0.5); }
        if f == "dot" || f == "cross" { return dim_add(dim_of(arg(v, 0)), dim_of(arg(v, 1)), 1); }
    }
    return no_dim();
}

fn name_dim(w: u8[:]) -> Dim {
    if axis_of(w) >= 0 || is_w(w) { return D_LENGTH; }
    if w == "N" { return dim_add(D_FORCE, D_LENGTH, -1); }
    if w == "Nx" { return D_FORCE; }
    if w == "eps_m" || w == "n" || is_theta(w) { return no_dim(); }
    let t = type_of_name(w);
    if t >= 0 { return type_dim(t); }
    return no_dim();
}
