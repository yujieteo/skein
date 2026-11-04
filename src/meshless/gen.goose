import std;
import math;
namespace meshless_gen;
import syntax from "../syntax";
import codegen.codegen;

fn unsupported(tok: i64, msg: u8[:], fix: u8[:]) { codegen::unsupported(tok, msg, fix); }
fn ok() -> bool { syntax::errors.len == 0 }
fn kind(n: i64) -> i64 { syntax::nodes[n].kind }
fn tok(n: i64) -> i64 { syntax::nodes[n].tok }
fn name(n: i64) -> const u8[:] { syntax::text(syntax::nodes[n].tok) }
fn first(n: i64) -> i64 { syntax::nodes[n].first }
fn next(n: i64) -> i64 { syntax::nodes[n].next }
fn child(n: i64, i: i64) -> i64 { syntax::child(n, i) }
fn source_of(n: i64) -> const u8[:] { codegen::source_of(n) }
fn set64(into: u8[..64]&, s: u8[:]) { codegen::set64(into, s); }
fn lit(out: u8[>..]&, x: f64) { codegen::lit(out, x); }
fn quoted(out: u8[>..]&, s: u8[:]) { codegen::quoted(out, s); }
fn is_call(n: i64, callee: u8[:]) -> bool {
    return n >= 0 && kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == callee;
}
fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }
fn arg_count(n: i64) -> i64 { syntax::child_count(n) - 1 }
fn binary(n: i64, op: u8[:]) -> bool { kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == op }

fn compact(n: i64) -> u8[>..] {
    var s: u8[>..] = [];
    for c in source_of(n) { if c != ' ' && c != '\n' { s.push(c); } }
    return s;
}

fn central_step(c: i64) -> bool {
    guard kind(c) == syntax::N_STEP else { return false; };
    let proc = child(c, 0);
    guard name(proc) == "transient" else { return false; };
    var a = first(proc);
    while a >= 0 {
        if syntax::text(syntax::nodes[a].aux) == "scheme" && source_of(first(a)) == "central_difference" { return true; }
        a = next(a);
    }
    return false;
}

fn claims(prog: i64) -> bool {
    var rkpm = false;
    var central = false;
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_METHOD && syntax::text(syntax::nodes[c].aux) == "rkpm" { rkpm = true; }
        if central_step(c) { central = true; }
        c = next(c);
    }
    return rkpm && central;
}

struct Cond {
    region: i64,
    comp: i64,
    value: i64,
}

struct Step {
    name: u8[..64],
    end: f64,
    dt: f64,
    first_bc: i64,
    bc_count: i64,
    every: i64,
    fields: i64,
}

struct Expect {
    node: i64,
    step: i64,
}

var shape_names: u8[..64][>..];
var shape_nodes: i64[>..];
var curve_names: u8[..64][>..];
var curve_nodes: i64[>..];
var local_names: u8[..64][>..];
var local_nodes: i64[>..];
var local_kinds: i64[>..];
var extra_names: u8[..64][>..];
var material_items: i64[>..];
var conds: Cond[>..];
var global_conds: i64[>..];
var active: i64[>..];
var reports: i64[>..];
var report_steps: i64[>..];
var steps: Step[>..];
var expects: Expect[>..];
var value_types: i64[>..];
var aggregates: i64[>..];
var physics_name: u8[..64];
var u_name: u8[..64];
var points_name: u8[..64];
var start_node = -1;
var stress_name: u8[..64];
var rho_node = -1;
var gravity: f64[2] = [0.0, 0.0];
var domain = -1;
var lo: f64[2] = [0.0, 0.0];
var hi: f64[2] = [0.0, 0.0];
var spacing = 0.0;
var cells: i64[2] = [0, 0];
var support = 0.0;
var row_table = -1;
var row_index = 0;
var in_chain = false;

fn reset() {
    syntax::errors.clear();
    codegen::reset();
    shape_names.clear();
    shape_nodes.clear();
    curve_names.clear();
    curve_nodes.clear();
    local_names.clear();
    local_nodes.clear();
    local_kinds.clear();
    extra_names.clear();
    material_items.clear();
    conds.clear();
    global_conds.clear();
    active.clear();
    reports.clear();
    report_steps.clear();
    steps.clear();
    expects.clear();
    aggregates.clear();
    value_types.clear();
    physics_name = [];
    u_name = [];
    points_name = [];
    stress_name = [];
    start_node = -1;
    rho_node = -1;
    gravity = [0.0, 0.0];
    domain = -1;
    spacing = 0.0;
    cells = [0, 0];
    support = 0.0;
    row_table = -1;
    in_chain = false;
}

struct Vec {
    n: i64,
    v: f64[2],
}

fn vector_value(n: i64) -> Vec {
    var r = Vec { n: 0, v: [0.0, 0.0] };
    var f = 1.0;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNIT {
            f = codegen::to_system(c);
        } else {
            if r.n < 2 { r.v[r.n] = codegen::eval(c).v[0]; }
            r.n++;
        }
        c = next(c);
    }
    if r.n != 2 { unsupported(tok(n), "a vector of a 2d model has 2 components", "write (x, y)"); }
    return Vec { n: 2, v: [r.v[0] * f, r.v[1] * f] };
}

recursive fn value(n: i64) -> Vec {
    if kind(n) == syntax::N_NAME {
        let v = codegen::find_value(name(n));
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let r = value(v);
            codegen::depth--;
            return r;
        }
    }
    if kind(n) == syntax::N_VECTOR { return vector_value(n); }
    if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        let a = value(first(n));
        return Vec { n: a.n, v: [-a.v[0], -a.v[1]] };
    }
    let r = codegen::eval(n);
    return Vec { n: 1, v: [r.v[0], 0.0] };
}

fn scalar_value(n: i64) -> f64 {
    let v = value(n);
    if v.n != 1 { unsupported(tok(n), "the value is a vector", "write a scalar"); }
    return v.v[0];
}

fn find_name(names: u8[..64][:], word: u8[:]) -> i64 {
    for s, i in names { if s == word { return i; } }
    return -1;
}

fn resolve(word: u8[:]) -> i64 {
    let v = codegen::find_value(word);
    if v >= 0 { return v; }
    return codegen::property(word);
}

fn is_field_member(n: i64) -> bool {
    return kind(n) == syntax::N_MEMBER && kind(first(n)) == syntax::N_NAME && name(first(n)) == u_name && (name(n) == "x" || name(n) == "y");
}

recursive fn varying(n: i64) -> bool {
    if kind(n) == syntax::N_NAME {
        let w = name(n);
        if w == "x" || w == "y" || w == "t" || w == u_name || w == "I" || find_name(local_names, w) >= 0 { return true; }
        let v = resolve(w);
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let r = varying(v);
            codegen::depth--;
            return r;
        }
        return false;
    }
    var c = first(n);
    while c >= 0 {
        if varying(c) { return true; }
        c = next(c);
    }
    return false;
}

let SCALAR = 0;
let TENSOR = 1;
let FUNCTIONS = ["sin", "cos", "tan", "exp", "sqrt", "abs", "tanh", "atan"];
let EXPR_FIX = "write + - * / ^ of numbers, lets, x, y, t, u.x, u.y, and sin, cos, tan, atan, exp, ln, sqrt, abs, tanh, min, max or interp";
let CHAIN_FIX = "write lets of defgrad(u), det, inv, transpose, tr, I, + - * / ^ and sin, cos, tan, atan, exp, ln, sqrt, abs, tanh, min, max, for example let F: Tensor2<Real> = defgrad(u)";

fn table_of(n: i64) -> i64 {
    let o = if kind(n) == syntax::N_MEMBER { first(n) } else { -1 };
    return if o >= 0 && kind(o) == syntax::N_NAME { find_name(codegen::table_names, name(o)) } else { -1 };
}

recursive fn has_table(n: i64) -> bool {
    if table_of(n) >= 0 { return true; }
    var c = first(n);
    while c >= 0 {
        if has_table(c) { return true; }
        c = next(c);
    }
    return false;
}

recursive fn kind_of(n: i64) -> i64 {
    let k = kind(n);
    if k == syntax::N_NAME {
        let w = name(n);
        if w == "I" && find_name(local_names, w) < 0 { return TENSOR; }
        let l = find_name(local_names, w);
        if l >= 0 { return local_kind(l); }
        return SCALAR;
    }
    if is_call(n, "defgrad") || is_call(n, "inv") || is_call(n, "transpose") { return TENSOR; }
    if k == syntax::N_UNARY { return kind_of(first(n)); }
    if k == syntax::N_BINARY {
        let a = kind_of(child(n, 0));
        let b = kind_of(child(n, 1));
        let op = syntax::text(tok(n));
        if op == "+" || op == "-" { return a; }
        if op == "*" { return if a == TENSOR || b == TENSOR { TENSOR } else { SCALAR }; }
        if op == "/" { return a; }
    }
    return SCALAR;
}

recursive fn local_kind(l: i64) -> i64 {
    if local_kinds[l] < 0 {
        local_kinds[l] = SCALAR;
        local_kinds[l] = kind_of(local_nodes[l]);
    }
    return local_kinds[l];
}

recursive fn emit(out: u8[>..]&, n: i64) {
    let k = kind(n);
    if k == syntax::N_MEMBER && table_of(n) >= 0 {
        emit_table_cell(out, n);
    } else if !varying(n) && !has_table(n) {
        out.push('(');
        lit(out, scalar_value(n));
        out.push(')');
    } else if k == syntax::N_NAME && (name(n) == "x" || name(n) == "y" || name(n) == "t") {
        out.append(name(n));
    } else if is_field_member(n) {
        out.append(if name(n) == "x" { "ux" } else { "uy" });
    } else if reduction_of(n) >= 0 {
        format(out, "aggregate_", aggregate_of(n), "()");
    } else if k == syntax::N_NAME && find_name(local_names, name(n)) >= 0 {
        if !in_chain {
            unsupported(tok(n), "a let of the physics is only in the let chain of the stress and in output", EXPR_FIX);
            return;
        }
        format(out, "l_", name(n));
    } else if k == syntax::N_NAME && name(n) == "I" {
        out.append("[1.0, 0.0, 0.0, 1.0]");
    } else if k == syntax::N_NAME && codegen::depth < 100 && resolve(name(n)) >= 0 {
        codegen::depth++;
        emit(out, resolve(name(n)));
        codegen::depth--;
    } else {
        emit_operation(out, n);
    }
}

recursive fn emit_operation(out: u8[>..]&, n: i64) {
    let k = kind(n);
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(n))) >= 0 {
        emit_arithmetic(out, n);
    } else if binary(n, "^") {
        if kind_of(child(n, 0)) == TENSOR { unsupported(tok(n), "a power of a tensor is not in the meshless generator", CHAIN_FIX); }
        out.append("pow(");
        emit(out, child(n, 0));
        out.append(", ");
        emit(out, child(n, 1));
        out.push(')');
    } else if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        if kind_of(first(n)) == TENSOR {
            out.append("t_scale(-1.0, ");
            emit(out, first(n));
            out.push(')');
        } else {
            out.append("(-");
            emit(out, first(n));
            out.push(')');
        }
    } else if is_call(n, "defgrad") && arg_count(n) == 1 && kind(arg(n, 0)) == syntax::N_NAME && name(arg(n, 0)) == u_name {
        out.append("f");
    } else if (is_call(n, "det") || is_call(n, "tr") || is_call(n, "inv") || is_call(n, "transpose")) && arg_count(n) == 1 && kind_of(arg(n, 0)) == TENSOR {
        format(out, "t_", name(child(n, 0)), "(");
        emit(out, arg(n, 0));
        out.push(')');
    } else if is_call(n, "interp") && arg_count(n) == 3 {
        emit_interp(out, n);
    } else if is_call(n, "when") && arg_count(n) == 2 {
        format(out, "when_time(series_", aggregate_of(arg(n, 0)), ", ");
        emit(out, arg(n, 1));
        out.push(')');
    } else if (is_call(n, "min") || is_call(n, "max")) && arg_count(n) == 2 {
        format(out, name(child(n, 0)), "(");
        emit(out, arg(n, 0));
        out.append(", ");
        emit(out, arg(n, 1));
        out.push(')');
    } else if is_call(n, "ln") && arg_count(n) == 1 {
        out.append("log(");
        emit(out, arg(n, 0));
        out.push(')');
    } else if k == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && position(FUNCTIONS, name(child(n, 0))) >= 0 && arg_count(n) == 1 {
        format(out, name(child(n, 0)), "(");
        emit(out, arg(n, 0));
        out.push(')');
    } else {
        unsupported(tok(n), str("the expression '", source_of(n), "' is not in the meshless generator"), CHAIN_FIX);
    }
}

fn emit_arithmetic(out: u8[>..]&, n: i64) {
    let op = syntax::text(tok(n));
    let a = child(n, 0);
    let b = child(n, 1);
    let ka = kind_of(a);
    let kb = kind_of(b);
    if ka == SCALAR && kb == SCALAR {
        out.push('(');
        emit(out, a);
        format(out, " ", op, " ");
        emit(out, b);
        out.push(')');
        return;
    }
    if op == "*" && ka == SCALAR {
        out.append("t_scale(");
        emit(out, a);
        out.append(", ");
        emit(out, b);
        out.push(')');
    } else if op == "*" && kb == SCALAR {
        out.append("t_scale(");
        emit(out, b);
        out.append(", ");
        emit(out, a);
        out.push(')');
    } else if op == "/" && kb == SCALAR {
        out.append("t_scale(1.0 / ");
        emit(out, b);
        out.append(", ");
        emit(out, a);
        out.push(')');
    } else if op == "*" && ka == TENSOR && kb == TENSOR {
        out.append("t_mul(");
        emit(out, a);
        out.append(", ");
        emit(out, b);
        out.push(')');
    } else if (op == "+" || op == "-") && ka == kb {
        format(out, "t_", if op == "+" { "add" } else { "sub" }, "(");
        emit(out, a);
        out.append(", ");
        emit(out, b);
        out.push(')');
    } else {
        unsupported(tok(n), str("the operation '", source_of(n), "' mixes a scalar and a tensor in a way that is not defined"), CHAIN_FIX);
    }
}

fn emit_interp(out: u8[>..]&, n: i64) {
    out.append("interp1(");
    for i in 2 {
        var col: f64[>..] = [];
        if !codegen::column_of(arg(n, i), col) { return; }
        out.push('[');
        for v, j in col {
            if j > 0 { out.append(", "); }
            lit(out, v);
        }
        out.append("], ");
    }
    emit(out, arg(n, 2));
    out.push(')');
}

fn emit_table_cell(out: u8[>..]&, n: i64) {
    let t = table_of(n);
    var col: f64[>..] = [];
    codegen::table_column(codegen::table_nodes[t], name(n), tok(n), col);
    if row_table != t || row_index >= col.len {
        unsupported(tok(n), "a column of a table is only in an expect that makes a test for each row", "write the column in the expect");
        return;
    }
    out.push('(');
    lit(out, col[row_index]);
    out.push(')');
}

let COND_FIX = "write a condition of x and y with == < <= > >= && || and !";

recursive fn emit_cond(out: u8[>..]&, n: i64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY { syntax::text(tok(n)) } else { "" };
    if op == "&&" || op == "||" {
        out.push('(');
        emit_cond(out, child(n, 0));
        format(out, ") ", op, " (");
        emit_cond(out, child(n, 1));
        out.push(')');
        return;
    }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "!" {
        out.append("!(");
        emit_cond(out, first(n));
        out.push(')');
        return;
    }
    let i = position(["==", "!=", "<", "<=", ">", ">="], op);
    if i < 0 {
        unsupported(tok(n), str("the condition '", source_of(n), "' is not in the meshless generator"), COND_FIX);
        return;
    }
    out.append(if i <= 1 { "abs(" } else { "(" });
    emit(out, child(n, 0));
    out.append(" - ");
    emit(out, child(n, 1));
    out.append([") <= tol", ") > tol", ") < -tol", ") <= tol", ") > tol", ") >= -tol"][i]);
}

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    if syntax::text(tok(child(n, 0))) != "2d" { unsupported(tok(child(n, 0)), "the meshless generator has 2d models", "write 2d in the header"); }
    let sys = child(n, 1);
    for s, i in codegen::SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            codegen::system = codegen::SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

fn value_item(n: i64) {
    let v = child(n, 1);
    var w: u8[..64] = [];
    set64(w, name(n));
    if is_call(v, "rect") {
        shape_names.push(w);
        shape_nodes.push(v);
    } else if is_call(v, "boundary") {
        curve_names.push(w);
        curve_nodes.push(v);
    } else {
        codegen::value_names.push(w);
        codegen::value_nodes.push(v);
        value_types.push(child(n, 0));
    }
}

let PHYSICS_FIX = "write physics water on column solves u: Vector<Length> = (0, 0) { material m, assume finite_strain, let F: Tensor2<Real> = defgrad(u), let P: Tensor2<Pressure> = ..., eq rho*ddt2(u) - div(P) = rho*(0, -g) }";

fn physics_item(n: i64) {
    set64(physics_name, name(n));
    var eqs: i64[>..] = [];
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            domain = find_name(shape_names, name(child(c, 0)));
            if domain < 0 || kind(child(c, 0)) != syntax::N_NAME { unsupported(tok(c), "the region of the physics is not a rect that a let names", "write physics water on column, with let column: Surface = rect(...)"); }
        } else if k == syntax::N_UNKNOWN {
            unknown_item(c);
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            var found = -1;
            for i in material_items { if name(i) == name(m) { found = i; } }
            if found < 0 { unsupported(tok(m), "the material is not a material item", "write material fluid { rho: Density = 1 kg/m^3 }"); }
            codegen::scope_material = found;
        } else if k == syntax::N_LET {
            local_names.push(default<u8[..64]>());
            set64(local_names[local_names.len - 1], name(c));
            local_nodes.push(child(c, 1));
            local_kinds.push(-1);
        } else if k == syntax::N_EQN {
            eqs.push(c);
        } else if k == syntax::N_ASSUME {
            assumptions(c);
        } else {
            unsupported(tok(c), "the physics statement is not in the meshless generator", PHYSICS_FIX);
        }
        c = next(c);
    }
    if domain < 0 && ok() { unsupported(tok(n), "the physics has no region", "write physics water on column"); }
    if !ok() { return; }
    equations(n, eqs);
}

fn assumptions(c: i64) {
    var a = first(c);
    while a >= 0 {
        if name(a) != "finite_strain" {
            unsupported(tok(a), str("the assumption '", name(a), "' is not in the meshless generator"), "write assume finite_strain: the stress is the first Piola-Kirchhoff stress of the start shape");
        }
        a = next(a);
    }
}

fn unknown_item(c: i64) {
    let t = child(c, 0);
    if kind(t) == syntax::N_TYPE_APPLY && name(t) == "Vector" && u_name.len == 0 {
        set64(u_name, name(c));
        start_node = child(c, 1);
    } else {
        unsupported(tok(c), "the unknown of the meshless generator is one Vector of Length", PHYSICS_FIX);
    }
}

fn is_zero(n: i64) -> bool { return kind(n) == syntax::N_NUMBER && codegen::eval(n).v[0] == 0.0; }

fn equations(n: i64, eqs: i64[:]) {
    if u_name.len == 0 { unsupported(tok(n), "the physics needs a Vector unknown", PHYSICS_FIX); return; }
    if eqs.len != 1 { unsupported(tok(n), "the physics has one eq", PHYSICS_FIX); return; }
    if !momentum_equation(eqs[0]) { unsupported(tok(eqs[0]), "the eq is not rho*ddt2(u) - div(P) = rho*(bx, by)", PHYSICS_FIX); }
}

fn momentum_equation(e: i64) -> bool {
    let lhs = child(e, 0);
    guard binary(lhs, "-") && binary(child(lhs, 0), "*") && kind(child(child(lhs, 0), 0)) == syntax::N_NAME && is_call(child(lhs, 1), "div") else { return false; };
    let rho = name(child(child(lhs, 0), 0));
    guard compact(child(child(lhs, 0), 1)) == str("ddt2(", u_name, ")") else { return false; };
    guard kind(arg(child(lhs, 1), 0)) == syntax::N_NAME else { return false; };
    let p = find_name(local_names, name(arg(child(lhs, 1), 0)));
    guard p >= 0 && resolve(rho) >= 0 else { return false; };
    set64(stress_name, name(arg(child(lhs, 1), 0)));
    if local_kind(p) != TENSOR { return false; }
    rho_node = resolve(rho);
    return is_zero(child(e, 1)) || body_force(child(e, 1), rho);
}

recursive fn body_force(n: i64, rho: u8[:]) -> bool {
    if binary(n, "+") { return body_force(child(n, 0), rho) && body_force(child(n, 1), rho); }
    guard binary(n, "*") && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == rho && kind(child(n, 1)) == syntax::N_VECTOR else { return false; };
    let g = vector_value(child(n, 1));
    gravity = [gravity[0] + g.v[0], gravity[1] + g.v[1]];
    return true;
}

fn material(n: i64) { material_items.push(n); }

let POINTS_FIX = "write points c on column { spacing: a/40 }";

fn points_item(n: i64) {
    if spacing > 0.0 {
        unsupported(tok(n), "a second points: the meshless generator has one", "remove the second points");
        return;
    }
    set64(points_name, name(n));
    let region = child(n, 0);
    if kind(region) != syntax::N_NAME || find_name(shape_names, name(region)) != domain {
        unsupported(tok(region), "the points are not on the region of the physics", POINTS_FIX);
    }
    var s = first(child(n, 1));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "spacing" {
            spacing = codegen::eval(first(s)).v[0];
        } else {
            unsupported(tok(s), "the points setting is not in the meshless generator", "use only spacing");
        }
        s = next(s);
    }
    if spacing <= 0.0 { unsupported(tok(n), "the points need a spacing > 0", POINTS_FIX); }
}

fn lattice_cells(n: i64) {
    let d = shape_nodes[domain];
    let a = vector_value(arg(d, 0));
    let b = vector_value(arg(d, 1));
    lo = a.v;
    hi = b.v;
    for k in 2 {
        let count = round((hi[k] - lo[k]) / spacing) as! i64;
        cells[k] = count;
        if count < 2 || abs(count * spacing - (hi[k] - lo[k])) > 1.0e-6 * spacing {
            unsupported(tok(n), "the spacing does not divide the sides of the rect", "write a spacing that divides both sides, for example a/40");
        }
    }
}

let METHOD_FIX = "write method water: rkpm { points: c, basis: linear, support: 3, integrate: scni, bc: nitsche }";

fn method_item(n: i64) {
    if name(n) != physics_name {
        unsupported(tok(n), "the method is not of the physics", METHOD_FIX);
        return;
    }
    if syntax::text(syntax::nodes[n].aux) != "rkpm" {
        unsupported(tok(n), "the meshless generator has the method rkpm", METHOD_FIX);
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        let v = first(s);
        let w = if kind(v) == syntax::N_NAME { name(v) } else { "" };
        if name(s) == "support" {
            support = codegen::eval(v).v[0];
        } else if name(s) == "points" {
            if w != points_name { unsupported(tok(v), "the points are not the points of the model", METHOD_FIX); }
        } else if name(s) == "basis" && w == "linear" {
        } else if name(s) == "integrate" && w == "scni" {
        } else if name(s) == "bc" && w == "nitsche" {
        } else {
            unsupported(tok(s), "the method setting is not in the meshless generator: it has points, basis: linear, support, integrate: scni and bc: nitsche", METHOD_FIX);
        }
        s = next(s);
    }
    if support <= 1.0 { unsupported(tok(n), "the method needs support > 1 (times the spacing)", METHOD_FIX); }
}

let STEP_FIX = "write step collapse: transient(end: 0.17 s, dt: 5.36 us, scheme: central_difference) { ... }";

fn procedure(proc: i64, s: Step&) {
    if name(proc) != "transient" {
        unsupported(tok(proc), "the meshless generator has transient steps", STEP_FIX);
        return;
    }
    var scheme = false;
    var a = first(proc);
    while a >= 0 {
        let w = syntax::text(syntax::nodes[a].aux);
        let v = first(a);
        if w == "end" {
            s.end = codegen::eval(v).v[0];
        } else if w == "dt" {
            s.dt = codegen::eval(v).v[0];
        } else if w == "scheme" && source_of(v) == "central_difference" {
            scheme = true;
        } else {
            unsupported(tok(a), "the meshless generator has transient(end:, dt:, scheme: central_difference)", STEP_FIX);
        }
        a = next(a);
    }
    if !scheme { unsupported(tok(proc), "the meshless generator has the scheme central_difference", "add scheme: central_difference"); }
    if s.end <= 0.0 || s.dt <= 0.0 { unsupported(tok(proc), "the step needs an end after the start and dt > 0 (the step is not adaptive)", STEP_FIX); }
}

let BC_FIX = "write bc left: u.x = 0, with left a let of boundary(column, condition)";

fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let r = if kind(region) == syntax::N_NAME { find_name(curve_names, name(region)) } else { -1 };
    if r < 0 {
        unsupported(tok(region), "the region of the bc is not a boundary that a let names", BC_FIX);
        return -1;
    }
    let target = child(n, 1);
    let value = child(n, 2);
    let first_cond = conds.len;
    if is_field_member(target) {
        conds.push(Cond { region: r, comp: if name(target) == "x" { 0 } else { 1 }, value: value });
    } else if kind(target) == syntax::N_NAME && name(target) == u_name && kind(value) == syntax::N_VECTOR && syntax::child_count(value) == 2 {
        conds.push(Cond { region: r, comp: 0, value: child(value, 0) });
        conds.push(Cond { region: r, comp: 1, value: child(value, 1) });
    } else {
        unsupported(tok(target), "the bc of the meshless generator is on u.x, u.y or u = (x, y)", BC_FIX);
        return -1;
    }
    if conds.len >= 62 {
        unsupported(tok(n), "the meshless generator has at most 62 bcs", "join bcs on the same region");
        return -1;
    }
    return first_cond;
}

fn output_item(n: i64, s: Step&) {
    var o = first(n);
    while o >= 0 {
        if kind(o) == syntax::N_EVERY {
            s.every = round(codegen::eval(first(o)).v[0]) as! i64;
        } else if kind(o) == syntax::N_NAME && compact(o) == u_name {
            s.fields |= 1;
        } else if kind(o) == syntax::N_NAME && find_name(local_names, name(o)) >= 0 && local_kind(find_name(local_names, name(o))) == SCALAR {
            var at = find_name(extra_names, name(o));
            if at < 0 {
                extra_names.push(default<u8[..64]>());
                set64(extra_names[extra_names.len - 1], name(o));
                at = extra_names.len - 1;
            }
            s.fields |= 2 << at;
        } else {
            unsupported(tok(o), "the output of the meshless generator is u and the scalar lets of the physics", "write output u, p every 4000");
        }
        o = next(o);
    }
    if s.every <= 0 { s.every = 1000000000; }
}

fn step(n: i64) {
    var s = Step { name: [], end: 0.0, dt: 0.0, first_bc: 0, bc_count: 0, every: 0, fields: 0 };
    set64(s.name, name(n));
    procedure(child(n, 0), s);
    var own: i64[>..] = [];
    var c = next(child(n, 0));
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            let b = conds.len;
            bc(c);
            for i in b..conds.len { own.push(i); }
        } else if k == syntax::N_OUTPUT {
            output_item(c, s);
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                reports.push(r);
                report_steps.push(steps.len);
                r = next(r);
            }
        } else if k == syntax::N_EXPECT {
            expects.push(Expect { node: c, step: steps.len });
        } else {
            unsupported(tok(c), "the step statement is not in the meshless generator", "use only bc, output, report and expect");
        }
        c = next(c);
    }
    s.first_bc = active.len;
    for g in global_conds {
        var replaced = false;
        for b in own { if conds[b].region == conds[g].region && conds[b].comp == conds[g].comp { replaced = true; } }
        if !replaced { active.push(g); }
    }
    active.append(own);
    s.bc_count = active.len - s.first_bc;
    steps.push(s);
}

fn reduction_of(n: i64) -> i64 {
    var r = n;
    var guard_depth = 0;
    while kind(r) == syntax::N_NAME && resolve(name(r)) >= 0 && guard_depth < 100 {
        r = resolve(name(r));
        guard_depth++;
    }
    let reduction = (is_call(r, "max") || is_call(r, "min")) && arg_count(r) == 2;
    if !reduction || kind(arg(r, 1)) != syntax::N_NAME || find_name(shape_names, name(arg(r, 1))) != domain { return -1; }
    return r;
}

fn aggregate_of(n: i64) -> i64 {
    let r = reduction_of(n);
    if r < 0 {
        unsupported(tok(n), str("'", source_of(n), "' is not a max or min of a field over the domain of the physics"),
                    "write let front: Length = max(x + u.x, column)");
        return -1;
    }
    for a, i in aggregates { if a == r { return i; } }
    aggregates.push(r);
    return aggregates.len - 1;
}

recursive fn register_reductions(n: i64) {
    if is_call(n, "when") && arg_count(n) == 2 {
        aggregate_of(arg(n, 0));
        register_reductions(arg(n, 1));
        return;
    }
    var c = first(n);
    while c >= 0 {
        register_reductions(c);
        c = next(c);
    }
}

fn reduction_unit(n: i64) -> u8[>..] {
    guard kind(n) == syntax::N_NAME else { return str(""); };
    let at = find_name(codegen::value_names, name(n));
    guard at >= 0 && at < value_types.len && kind(value_types[at]) == syntax::N_TYPE_NAME else { return str(""); };
    let t = name(value_types[at]);
    if t == "Length" { return str(codegen::system.length); }
    if t == "Time" { return str("s"); }
    if t == "Pressure" { return str(codegen::system.pressure); }
    return str("");
}

fn read(prog: i64) -> bool {
    reset();
    var physics_nodes: i64[>..] = [];
    var points_nodes: i64[>..] = [];
    var methods: i64[>..] = [];
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
        } else if k == syntax::N_TABLE {
            var w: u8[..64] = [];
            set64(w, name(c));
            codegen::table_names.push(w);
            codegen::table_nodes.push(c);
        } else if k == syntax::N_MATERIAL {
            material(c);
        } else if k == syntax::N_PHYSICS {
            physics_nodes.push(c);
        } else if k == syntax::N_POINTS {
            points_nodes.push(c);
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in the meshless generator"), "remove it");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_nodes.len != 1 { unsupported(at, "the meshless generator has one physics", PHYSICS_FIX); }
    if !ok() { return false; }
    physics_item(physics_nodes[0]);
    if !ok() { return false; }
    if points_nodes.len != 1 { unsupported(at, "the model has one points", POINTS_FIX); }
    if !ok() { return false; }
    points_item(points_nodes[0]);
    if ok() { lattice_cells(points_nodes[0]); }
    if methods.len != 1 { unsupported(at, "the model has one method", METHOD_FIX); }
    if !ok() { return false; }
    method_item(methods[0]);
    if !ok() { return false; }
    read_steps(prog, at);
    check_program();
    return ok();
}

fn read_steps(prog: i64, at: i64) {
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            let b = conds.len;
            bc(c);
            for i in b..conds.len { global_conds.push(i); }
        } else if k == syntax::N_STEP {
            step(c);
        } else if k == syntax::N_EXPECT {
            expects.push(Expect { node: c, step: -1 });
        }
        c = next(c);
    }
    if steps.len == 0 { unsupported(at, "the model has no step", STEP_FIX); }
    for x in expects { register_reductions(child(x.node, 0)); register_reductions(child(x.node, 1)); }
    for r in reports { aggregate_of(first(r)); }
}

fn check_program() {
    var scratch: u8[>..] = [];
    emit_functions(scratch);
    for x in expects { emit_expect(scratch, x); }
}

fn dof_estimate() -> i64 { return 2 * (cells[0] + 1) * (cells[1] + 1); }

fn emit_vector(out: u8[>..]&, n: i64) {
    if n < 0 || kind(n) != syntax::N_VECTOR {
        out.append("[0.0, 0.0]");
        return;
    }
    var f = 1.0;
    var parts: i64[>..] = [];
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNIT { f = codegen::to_system(c); } else { parts.push(c); }
        c = next(c);
    }
    if parts.len != 2 { unsupported(tok(n), "a vector of a 2d model has 2 components", "write (x, y)"); return; }
    out.push('[');
    for i in 2 {
        if i > 0 { out.append(", "); }
        lit(out, f);
        out.append(" * ");
        emit(out, parts[i]);
    }
    out.push(']');
}

fn emit_chain(out: u8[>..]&, fname: u8[:], ret: u8[:]) {
    let r = find_name(local_names, ret);
    let tensor = local_kind(r) == TENSOR;
    format(out, "fn ", fname, "(f: f64[4]) -> ", if tensor { "f64[4]" } else { "f64" }, " {\n");
    in_chain = true;
    for i in 0..r + 1 {
        format(out, "    let l_", local_names[i], if local_kind(i) == TENSOR { ": f64[4]" } else { "" }, " = ");
        emit(out, local_nodes[i]);
        out.append(";\n");
    }
    in_chain = false;
    format(out, "    return l_", ret, ";\n}\n\n");
}

fn emit_functions(out: u8[>..]&) {
    emit_chain(out, "stress", stress_name);
    for e in extra_names { emit_chain(out, str("chain_", e), e); }
    out.append("fn start_displacement(x: f64, y: f64) -> f64[2] {\n    let t = 0.0;\n    let ux = 0.0;\n    let uy = 0.0;\n    return ");
    emit_vector(out, start_node);
    out.append(";\n}\n\nfn boundary(k: i64, x: f64, y: f64, t: f64) -> f64 {\n    let ux = 0.0;\n    let uy = 0.0;\n");
    for b, k in conds {
        format(out, "    if k == ", k, " { return ");
        emit(out, b.value);
        out.append("; }\n");
    }
    out.append("    return 0.0;\n}\n\n");
    for c, i in curve_names {
        format(out, "fn select_", c, "() -> bool {\n    return meshless::select_faces() { x, y => ");
        if syntax::child_count(curve_nodes[i]) <= 2 { out.append("true"); } else { emit_cond(out, arg(curve_nodes[i], 1)); }
        out.append(" };\n}\n\n");
    }
    for r, i in aggregates { emit_reduction(out, r, i); }
    out.append("fn record(t: f64) {\n    times.push(t);\n");
    for r, i in aggregates { format(out, "    series_", i, ".push(aggregate_", i, "());\n"); }
    out.append("}\n\n");
}

fn emit_reduction(out: u8[>..]&, r: i64, i: i64) {
    let most = is_call(r, "max");
    format(out, "var series_", i, ": f64[>..] = [];\n\nfn field_", i, "(x: f64, y: f64, ux: f64, uy: f64) -> f64 {\n    let t = 0.0;\n    return ");
    emit(out, arg(r, 0));
    format(out, ";\n}\n\nfn aggregate_", i, "() -> f64 {\n    var m = ", if most { "-1.0e300" } else { "1.0e300" }, ";\n    for p in meshless::point_count() {\n");
    format(out, "        m = ", if most { "max" } else { "min" }, "(m, field_", i, "(meshless::px(p), meshless::py(p), meshless::node_u(0, p), meshless::node_u(1, p)));\n    }\n    return m;\n}\n\n");
}

fn emit_report(out: u8[>..]&, step: u8[:], r: i64) {
    let i = aggregate_of(first(r));
    format(out, "    // report ", source_of(r), "\n    {\n        let name = ");
    quoted(out, source_of(r));
    out.append(";\n        report(");
    quoted(out, step);
    format(out, ", name, quantity(aggregate_", i, "(), ");
    quoted(out, reduction_unit(first(r)));
    out.append("));\n    }\n");
}

recursive fn table_in(n: i64) -> i64 {
    if table_of(n) >= 0 { return table_of(n); }
    var c = first(n);
    while c >= 0 {
        let t = table_in(c);
        if t >= 0 { return t; }
        c = next(c);
    }
    return -1;
}

fn row_count(t: i64) -> i64 {
    var rows = 0;
    var c = first(codegen::table_nodes[t]);
    while c >= 0 {
        if kind(c) == syntax::N_ROW { rows++; }
        c = next(c);
    }
    return rows;
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    var t = table_in(child(n, 0));
    if t < 0 { t = table_in(child(n, 1)); }
    let rows = if t >= 0 { row_count(t) } else { 1 };
    let within = child(n, 2);
    let relative = kind(within) == syntax::N_PERCENT;
    var tolerance = 0.0;
    if relative {
        let text = syntax::text(tok(within));
        let v, good = parse_flt(text[..text.len - 1]);
        tolerance = v / 100.0;
    } else {
        tolerance = codegen::eval(within).v[0];
    }
    for r in rows {
        row_table = t;
        row_index = r;
        format(out, "    // expect ", name(n), "\n    {\n        let name = ");
        let label: u8[>..] = if rows > 1 { str(name(n), "[", r + 1, "]") } else { str(name(n)) };
        quoted(out, label);
        out.append(";\n        when_missing = false;\n        expect(name, ");
        emit(out, child(n, 0));
        out.append(", ");
        emit(out, child(n, 1));
        out.append(", ");
        lit(out, tolerance);
        format(out, ", ", if relative { "true" } else { "false" }, ", ");
        quoted(out, str(source_of(child(n, 1)), " within ", source_of(within)));
        out.append(");\n    }\n");
    }
    row_table = -1;
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "fn after_", s.name, "(k: i64, t: f64) {\n    record(t);\n");
    if s.fields != 0 {
        format(out, "    if k % ", s.every, " == 0 {\n        let path = str(");
        quoted(out, str(results_dir, "/", s.name, "_"));
        out.append(", k, \".vtu\");\n        var names: u8[>..] = [];\n        var values: f64[>..] = [];\n");
        for e, i in extra_names {
            if (s.fields & (2 << i)) != 0 {
                format(out, "        if names.len > 0 { names.push(','); }\n        names.append(\"", e, "\");\n");
                format(out, "        for p in meshless::point_count() { values.push(chain_", e, "(meshless::node_gradient(p))); }\n");
            }
        }
        format(out, "        output(path, meshless::write_dynamic_vtu(path, ", if (s.fields & 1) != 0 { "true" } else { "false" }, ", names, values));\n    }\n");
    }
    format(out, "}\n\n// step ", s.name, ": transient, central_difference\nfn step_", s.name, "() {\n    meshless::clear_bcs();\n");
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        let c = curve_names[conds[b].region];
        format(out, "    region(select_", c, "(), \"", c, "\");\n    meshless::fix_component(", conds[b].comp, ", ", b, ");\n");
    }
    out.append("    times.clear();\n");
    for r, i in aggregates { format(out, "    series_", i, ".clear();\n"); }
    out.append("    record(meshless::time);\n    solved(");
    quoted(out, s.name);
    out.append(", meshless::solve_dynamic(");
    lit(out, s.end);
    out.append(", ");
    lit(out, s.dt);
    out.append(", gravity_x, gravity_y, stress, boundary, ");
    format(out, "after_", s.name, "));\n");
    for i in 0..reports.len { if report_steps[i] == si { emit_report(out, s.name, reports[i]); } }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    out.append("}\n\n");
}

let PRELUDE = """
    var reports: u8[>..] = [];
    var report_count = 0;
    var checks: u8[>..] = [];
    var check_count = 0;
    var solves: u8[>..] = [];
    var solve_count = 0;
    var results: u8[>..] = [];
    var result_count = 0;
    var failed = false;
    var tol = 0.0;
    var times: f64[>..] = [];
    var when_missing = false;

    fn t_mul(a: f64[4], b: f64[4]) -> f64[4] {
        return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]];
    }

    fn t_add(a: f64[4], b: f64[4]) -> f64[4] { return [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]; }

    fn t_sub(a: f64[4], b: f64[4]) -> f64[4] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]; }

    fn t_scale(s: f64, a: f64[4]) -> f64[4] { return [s * a[0], s * a[1], s * a[2], s * a[3]]; }

    fn t_det(a: f64[4]) -> f64 { return a[0] * a[3] - a[1] * a[2]; }

    fn t_tr(a: f64[4]) -> f64 { return a[0] + a[3]; }

    fn t_transpose(a: f64[4]) -> f64[4] { return [a[0], a[2], a[1], a[3]]; }

    fn t_inv(a: f64[4]) -> f64[4] {
        let d = a[0] * a[3] - a[1] * a[2];
        return [a[3] / d, -a[1] / d, -a[2] / d, a[0] / d];
    }

    fn toon(into: u8[>..]&, s: u8[:]) {
        var needs = s.len == 0;
        for c in s {
            if c == ',' || c == '"' || c == ':' || c == '\\' || c == ' ' { needs = true; }
        }
        if !needs { into.append(s); return; }
        into.push('"');
        for c in s {
            if c == '"' || c == '\\' { into.push('\\'); }
            into.push(c);
        }
        into.push('"');
    }

    fn quantity(x: f64, unit: u8[:]) -> u8[>..] { if unit.len == 0 { str(x) } else { str(x, " ", unit) } }

    fn report(step: u8[:], expr: u8[:], value: u8[:]) {
        reports.append("  ");
        toon(reports, step);
        reports.push(',');
        toon(reports, expr);
        reports.push(',');
        toon(reports, value);
        reports.push('\n');
        report_count++;
    }

    fn check(stage: u8[:], id: u8[:], pass: bool, value: u8[:], limit: u8[:]) {
        checks.append("  ");
        checks.append(stage);
        checks.push(',');
        toon(checks, id);
        checks.append(if pass { ",pass," } else { ",fail," });
        toon(checks, value);
        checks.push(',');
        toon(checks, limit);
        checks.push('\n');
        check_count++;
        if !pass { failed = true; }
    }

    fn region(found: bool, name: u8[:]) {
        if !found {
            write_stderr(str("skein run: no boundary edge of the points is on the region '", name, "'\n"));
            exit(1);
        }
    }

    fn solved(step: u8[:], finished: bool) {
        format(solves, "  ", step, ",", meshless::point_count(), ",", meshless::dof_count(), ",", meshless::dyn_steps, ",0,", meshless::dyn_seconds, ",", if finished { "yes" } else { "no" }, "\n");
        solve_count++;
        check("solve", "finite", finished, str(meshless::dyn_steps, " steps"), "each displacement finite and below 1000 sizes of the body");
    }

    fn interp1(xs: f64[:], ys: f64[:], x: f64) -> f64 {
        for i in 1..xs.len {
            if x >= xs[i - 1] && x <= xs[i] {
                if xs[i] == xs[i - 1] { return ys[i]; }
                return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
            }
        }
        write_stderr("skein run: interp is outside the first column\n");
        failed = true;
        return 0.0;
    }

    fn when_time(series: f64[:], value: f64) -> f64 {
        if series.len > 0 && series[0] >= value { return times[0]; }
        for i in 1..series.len {
            if series[i] >= value {
                return times[i - 1] + (times[i] - times[i - 1]) * (value - series[i - 1]) / (series[i] - series[i - 1]);
            }
        }
        when_missing = true;
        return 0.0;
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", !when_missing && abs(got - want) <= t, str(name, ": ", if when_missing { str("the value is not reached") } else { quantity(got, "") }), limit);
    }

    fn output(path: u8[:], written: bool) {
        if !written {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    fn print_records() {
        var s: u8[>..] = [];
        format(s, "report[", report_count, "]{step,expr,value}:\n", reports);
        format(s, "verify[", check_count, "]{stage,id,result,value,limit}:\n", checks);
        format(s, "solve[", solve_count, "]{step,cells,dofs,iterations,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

fn generate(meshless_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    out.append("// The values are in the unit system of the model header.\nimport std;\nimport os;\nimport math;\nimport meshless from ");
    quoted(out, meshless_path);
    out.append(";\nimport solvers from ");
    quoted(out, str(meshless_path[..meshless_path.len - 8], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    out.append("var gravity_x = 0.0;\nvar gravity_y = 0.0;\n\n");
    emit_functions(out);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n    // points: a lattice of ");
    format(out, cells[0] + 1, " x ", cells[1] + 1, ", the spacing ", spacing, "\n    meshless::points_lattice(");
    lit(out, lo[0]);
    out.append(", ");
    lit(out, lo[1]);
    out.append(", ");
    lit(out, hi[0]);
    out.append(", ");
    lit(out, hi[1]);
    format(out, ", ", cells[0], ", ", cells[1], ");\n    tol = meshless::tolerance();\n");
    format(out, "    // physics ", physics_name, ": total Lagrangian solid, the displacement ", u_name, "\n    gravity_x = ");
    lit(out, gravity[0]);
    out.append(";\n    gravity_y = ");
    lit(out, gravity[1]);
    out.append(";\n    if !meshless::explicit_setup(");
    lit(out, support);
    out.append(", ");
    lit(out, scalar_value(rho_node));
    out.append(") {\n        write_stderr(\"skein run: the support is too small for the points: a shape function is singular\\n\");\n        exit(1);\n    }\n    meshless::start_displacement(start_displacement);\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
