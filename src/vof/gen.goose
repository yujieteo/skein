import std;
import math;
namespace vof_gen;
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

fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS {
            var velocity = false;
            var pressure = false;
            var real = false;
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_UNKNOWN {
                    let t = child(u, 0);
                    if kind(t) == syntax::N_TYPE_NAME && name(t) == "Real" { real = true; }
                    if kind(t) == syntax::N_TYPE_NAME && name(t) == "Pressure" { pressure = true; }
                    if kind(t) == syntax::N_TYPE_APPLY && name(t) == "Vector" && kind(first(t)) == syntax::N_TYPE_NAME && name(first(t)) == "Velocity" { velocity = true; }
                }
                u = next(u);
            }
            if velocity && pressure && real { return true; }
        }
        c = next(c);
    }
    return false;
}

struct Bc {
    region: i64,
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

let FLUX_NAMES = ["upwind", "linear_upwind", "central"];
let FLUX_CONSTANTS = ["vof::UPWIND", "vof::LINEAR_UPWIND", "vof::CENTRAL"];
let LIMITER_NAMES = ["none", "minmod", "van_leer"];
let LIMITER_CONSTANTS = ["vof::NONE", "vof::MINMOD", "vof::VAN_LEER"];

var shape_names: u8[..64][>..];
var shape_nodes: i64[>..];
var curve_names: u8[..64][>..];
var curve_nodes: i64[>..];
var local_names: u8[..64][>..];
var local_nodes: i64[>..];
var material_items: i64[>..];
var bcs: Bc[>..];
var global_bcs: i64[>..];
var active: i64[>..];
var reports: i64[>..];
var report_steps: i64[>..];
var steps: Step[>..];
var expects: Expect[>..];
var value_types: i64[>..];
var aggregates: i64[>..];
var physics_name: u8[..64];
var v_name: u8[..64];
var p_name: u8[..64];
var a_name: u8[..64];
var initial_node = -1;
var velocity_node = -1;
var rho_node = -1;
var mu_node = -1;
var gravity: f64[2] = [0.0, 0.0];
var sigma_node = -1;
var domain = -1;
var lo: f64[2] = [0.0, 0.0];
var hi: f64[2] = [0.0, 0.0];
var mesh_size = 0.0;
var cells: i64[2] = [0, 0];
var flux = 1;
var limiter = 2;
var row_table = -1;
var row_index = 0;

fn reset() {
    syntax::errors.clear();
    codegen::reset();
    shape_names.clear();
    shape_nodes.clear();
    curve_names.clear();
    curve_nodes.clear();
    local_names.clear();
    local_nodes.clear();
    material_items.clear();
    bcs.clear();
    global_bcs.clear();
    active.clear();
    reports.clear();
    report_steps.clear();
    steps.clear();
    expects.clear();
    aggregates.clear();
    value_types.clear();
    physics_name = [];
    v_name = [];
    p_name = [];
    a_name = [];
    initial_node = -1;
    velocity_node = -1;
    rho_node = -1;
    mu_node = -1;
    gravity = [0.0, 0.0];
    sigma_node = -1;
    domain = -1;
    mesh_size = 0.0;
    cells = [0, 0];
    flux = 1;
    limiter = 2;
    row_table = -1;
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
    let l = find_name(local_names, word);
    if l >= 0 { return local_nodes[l]; }
    let v = codegen::find_value(word);
    if v >= 0 { return v; }
    return codegen::property(word);
}

recursive fn varying(n: i64) -> bool {
    if kind(n) == syntax::N_NAME {
        let w = name(n);
        if w == "x" || w == "y" || w == "t" || w == a_name || w == p_name || w == v_name { return true; }
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

let FUNCTIONS = ["sin", "cos", "tan", "exp", "sqrt", "abs", "tanh", "atan"];
let EXPR_FIX = "write + - * / ^ of numbers, lets, x, y, t, alpha, p, v.x, v.y, and sin, cos, tan, atan, exp, ln, sqrt, abs, tanh, min, max or interp";

fn table_of(n: i64) -> i64 {
    let o = if kind(n) == syntax::N_MEMBER { first(n) } else { -1 };
    return if o >= 0 && kind(o) == syntax::N_NAME { find_name(codegen::table_names, name(o)) } else { -1 };
}

recursive fn emit(out: u8[>..]&, n: i64) {
    let k = kind(n);
    if k == syntax::N_MEMBER && table_of(n) >= 0 {
        emit_table_cell(out, n);
    } else if !varying(n) && !has_table(n) && !(k == syntax::N_NAME && find_name(local_names, name(n)) >= 0) {
        out.push('(');
        lit(out, scalar_value(n));
        out.push(')');
    } else if k == syntax::N_NAME && (name(n) == "x" || name(n) == "y" || name(n) == "t") {
        out.append(name(n));
    } else if k == syntax::N_NAME && name(n) == a_name {
        out.push('a');
    } else if k == syntax::N_NAME && name(n) == p_name {
        out.push('p');
    } else if k == syntax::N_MEMBER && kind(first(n)) == syntax::N_NAME && name(first(n)) == v_name && (name(n) == "x" || name(n) == "y") {
        out.append(if name(n) == "x" { "vx" } else { "vy" });
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
        out.push('(');
        emit(out, child(n, 0));
        format(out, " ", syntax::text(tok(n)), " ");
        emit(out, child(n, 1));
        out.push(')');
    } else if binary(n, "^") {
        out.append("pow(");
        emit(out, child(n, 0));
        out.append(", ");
        emit(out, child(n, 1));
        out.push(')');
    } else if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        out.append("(-");
        emit(out, first(n));
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
        unsupported(tok(n), str("the expression '", source_of(n), "' is not in the two-phase generator"), EXPR_FIX);
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
        unsupported(tok(n), str("the condition '", source_of(n), "' is not in the two-phase generator"), COND_FIX);
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
    if syntax::text(tok(child(n, 0))) != "2d" { unsupported(tok(child(n, 0)), "the two-phase generator has 2d models", "write 2d in the header"); }
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

let PHYSICS_FIX = "write physics flow on tank solves v: Vector<Velocity>, p: Pressure, alpha: Real = 0 { material fluid, let rho, mu, kappa and sigma, eq rho*(ddt(v) + grad(v)*v) - div(sigma) = rho*(0, -g) + s*kappa*grad(alpha), eq div(v) = 0, eq ddt(alpha) + div(alpha*v) = 0 }";

fn physics_item(n: i64) {
    set64(physics_name, name(n));
    var eqs: i64[>..] = [];
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            domain = find_name(shape_names, name(child(c, 0)));
            if domain < 0 || kind(child(c, 0)) != syntax::N_NAME { unsupported(tok(c), "the region of the physics is not a rect that a let names", "write physics flow on tank, with let tank: Surface = rect(...)"); }
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
        } else if k == syntax::N_EQN {
            eqs.push(c);
        } else if k == syntax::N_ASSUME {
            assumptions(c);
        } else {
            unsupported(tok(c), "the physics statement is not in the two-phase generator", PHYSICS_FIX);
        }
        c = next(c);
    }
    if domain < 0 && ok() { unsupported(tok(n), "the physics has no region", "write physics flow on tank"); }
    if !ok() { return; }
    equations(n, eqs);
}

fn assumptions(c: i64) {
    var a = first(c);
    while a >= 0 {
        if position(["incompressible", "newtonian", "laminar"], name(a)) < 0 {
            unsupported(tok(a), str("the assumption '", name(a), "' is not in the two-phase generator"), "write assume incompressible, newtonian, laminar");
        }
        a = next(a);
    }
}

fn unknown_item(c: i64) {
    let t = child(c, 0);
    if kind(t) == syntax::N_TYPE_APPLY && name(t) == "Vector" && v_name.len == 0 {
        set64(v_name, name(c));
        velocity_node = child(c, 1);
    } else if kind(t) == syntax::N_TYPE_NAME && name(t) == "Pressure" && p_name.len == 0 {
        set64(p_name, name(c));
        if child(c, 1) >= 0 { unsupported(tok(c), "the pressure starts at 0 in the two-phase generator", "remove the start value of the pressure"); }
    } else if kind(t) == syntax::N_TYPE_NAME && name(t) == "Real" && a_name.len == 0 {
        set64(a_name, name(c));
        initial_node = child(c, 1);
    } else {
        unsupported(tok(c), "the unknowns of the two-phase generator are one Vector, one Pressure and one Real", PHYSICS_FIX);
    }
}

fn is_zero(n: i64) -> bool { return kind(n) == syntax::N_NUMBER && codegen::eval(n).v[0] == 0.0; }

fn equations(n: i64, eqs: i64[:]) {
    if v_name.len == 0 || p_name.len == 0 || a_name.len == 0 { unsupported(tok(n), "the physics needs a Vector, a Pressure and a Real unknown", PHYSICS_FIX); return; }
    var momentum = false;
    var divergence = false;
    var transport = false;
    for e in eqs {
        let l = compact(child(e, 0));
        if l == str("div(", v_name, ")") && is_zero(child(e, 1)) && !divergence {
            divergence = true;
        } else if (l == str("ddt(", a_name, ")+div(", a_name, "*", v_name, ")") || l == str("ddt(", a_name, ")+div(", v_name, "*", a_name, ")")) && is_zero(child(e, 1)) && !transport {
            transport = true;
        } else if !momentum && momentum_equation(e) {
            momentum = true;
        } else {
            unsupported(tok(e), "the eq is not an equation of the two-phase generator", PHYSICS_FIX);
        }
    }
    if ok() && !(momentum && divergence && transport) { unsupported(tok(n), "the physics needs the momentum eq, eq div(v) = 0 and the transport eq of alpha", PHYSICS_FIX); }
}

fn momentum_equation(e: i64) -> bool {
    let lhs = child(e, 0);
    guard binary(lhs, "-") && binary(child(lhs, 0), "*") && kind(child(child(lhs, 0), 0)) == syntax::N_NAME && is_call(child(lhs, 1), "div") else { return false; };
    let rho = name(child(child(lhs, 0), 0));
    let inertia = compact(child(child(lhs, 0), 1));
    let want = str("ddt(", v_name, ")+grad(", v_name, ")*", v_name);
    guard inertia == want || inertia == str("(", want, ")") else { return false; };
    guard kind(arg(child(lhs, 1), 0)) == syntax::N_NAME else { return false; };
    let sigma = resolve(name(arg(child(lhs, 1), 0)));
    guard sigma >= 0 && find_name(local_names, rho) >= 0 else { return false; };
    let mu = viscosity_name(compact(sigma));
    guard mu.len > 0 && find_name(local_names, mu) >= 0 else { return false; };
    rho_node = local_nodes[find_name(local_names, rho)];
    mu_node = local_nodes[find_name(local_names, mu)];
    return is_zero(child(e, 1)) || forces(child(e, 1), rho);
}

fn viscosity_name(text: u8[:]) -> u8[>..] {
    let head = str("-", p_name, "*I+2*");
    let tail = str("*sym(grad(", v_name, "))");
    guard starts_with(text, head) && ends_with(text, tail) && text.len > head.len + tail.len else { return str(""); };
    return str(text[head.len..text.len - tail.len]);
}

recursive fn forces(n: i64, rho: u8[:]) -> bool {
    if binary(n, "+") { return forces(child(n, 0), rho) && forces(child(n, 1), rho); }
    if binary(n, "*") && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == rho && kind(child(n, 1)) == syntax::N_VECTOR {
        let g = vector_value(child(n, 1));
        gravity = g.v;
        return true;
    }
    guard binary(n, "*") && binary(child(n, 0), "*") && is_call(child(n, 1), "grad") && kind(arg(child(n, 1), 0)) == syntax::N_NAME && name(arg(child(n, 1), 0)) == a_name else { return false; };
    let k = child(child(n, 0), 1);
    guard kind(k) == syntax::N_NAME && resolve(name(k)) >= 0 && sigma_node < 0 else { return false; };
    guard compact(resolve(name(k))) == str("-div(grad(", a_name, ")/norm(grad(", a_name, ")))") else { return false; };
    sigma_node = child(child(n, 0), 0);
    return true;
}

fn material(n: i64) { material_items.push(n); }

let MESH_FIX = "write mesh m on tank { size: a/40, cell: quad }";

fn mesh_item(n: i64) {
    if mesh_size > 0.0 {
        unsupported(tok(n), "a second mesh: the two-phase generator has one", "remove the second mesh");
        return;
    }
    let region = child(n, 0);
    if kind(region) != syntax::N_NAME || find_name(shape_names, name(region)) != domain {
        unsupported(tok(region), "the mesh is not on the region of the physics", MESH_FIX);
    }
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            mesh_size = codegen::eval(v).v[0];
        } else if kind(s) == syntax::N_SETTING && name(s) == "cell" && source_of(v) == "quad" {
        } else {
            unsupported(tok(s), "the mesh setting is not in the two-phase generator", "use only size and cell: quad");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", MESH_FIX); }
}

fn mesh_cells(n: i64) {
    let d = shape_nodes[domain];
    let a = vector_value(arg(d, 0));
    let b = vector_value(arg(d, 1));
    lo = a.v;
    hi = b.v;
    for k in 2 {
        let count = round((hi[k] - lo[k]) / mesh_size) as! i64;
        cells[k] = count;
        if count < 2 || abs(count * mesh_size - (hi[k] - lo[k])) > 1.0e-6 * mesh_size {
            unsupported(tok(n), "the size of the mesh does not divide the sides of the rect", "write a size that divides both sides, for example a/40");
        }
    }
}

fn method_item(n: i64) {
    if name(n) != physics_name {
        unsupported(tok(n), "the method is not of the physics", "write method flow: fv { mesh: m }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        let v = first(s);
        let w = if kind(v) == syntax::N_NAME { name(v) } else { "" };
        if name(s) == "limiter" {
            limiter = position(LIMITER_NAMES, w);
            if limiter < 0 {
                limiter = 2;
                unsupported(tok(v), "the limiter of the two-phase generator is none, minmod or van_leer", "write limiter: van_leer");
            }
        } else if name(s) == "flux" {
            flux = position(FLUX_NAMES, w);
            if flux < 0 {
                flux = 1;
                unsupported(tok(v), "the flux of the two-phase generator is upwind, linear_upwind or central", "write flux: linear_upwind");
            }
        } else if name(s) == "gradient" {
            if w != "green_gauss" { unsupported(tok(v), str("the gradient '", w, "' is not in the two-phase generator: least_squares is unsupported, the gradient is green_gauss"), "write gradient: green_gauss"); }
        } else if name(s) == "pv_coupling" {
            if w != "rhie_chow" { unsupported(tok(v), "the pv_coupling of the two-phase generator is rhie_chow (the face flux of the fractional step)", "write pv_coupling: rhie_chow"); }
        } else if name(s) != "mesh" {
            unsupported(tok(s), "the method setting is not in the two-phase generator", "use mesh, gradient, flux, limiter and pv_coupling");
        }
        s = next(s);
    }
}

let STEP_FIX = "write step collapse: transient(end: 0.18 s, dt: 0.4 ms, scheme: backward_euler) { ... }";

fn procedure(proc: i64, s: Step&) {
    if name(proc) != "transient" {
        unsupported(tok(proc), "the two-phase generator has transient steps", STEP_FIX);
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
            s.dt = if kind(v) == syntax::N_NAME && name(v) == "stable" { 0.0 } else { codegen::eval(v).v[0] };
        } else if w == "scheme" && source_of(v) == "backward_euler" {
            scheme = true;
        } else {
            unsupported(tok(a), "the two-phase generator has transient(end:, dt:, scheme: backward_euler)", STEP_FIX);
        }
        a = next(a);
    }
    if !scheme { unsupported(tok(proc), "the two-phase generator has the scheme backward_euler: the pressure and the divergence of the step are implicit", "add scheme: backward_euler"); }
    if s.end <= 0.0 || s.dt <= 0.0 { unsupported(tok(proc), "the step needs an end after the start and dt > 0 (the step is not adaptive)", STEP_FIX); }
}

let BC_FIX = "write bc walls: v = 0, with walls a let of boundary(tank, condition)";

fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let r = if kind(region) == syntax::N_NAME { find_name(curve_names, name(region)) } else { -1 };
    if r < 0 {
        unsupported(tok(region), "the region of the bc is not a boundary that a let names", BC_FIX);
        return -1;
    }
    if compact(child(n, 1)) != v_name {
        unsupported(tok(child(n, 1)), "the bc of the two-phase generator is on the velocity v (a side with no bc is open: p = 0)", BC_FIX);
        return -1;
    }
    let value = child(n, 2);
    if !(kind(value) == syntax::N_VECTOR || is_zero(value)) {
        unsupported(tok(value), "the value of the velocity is a vector or 0", "write bc walls: v = 0, or v = (1 m/s, 0 m/s)");
        return -1;
    }
    if bcs.len >= 62 {
        unsupported(tok(n), "the two-phase generator has at most 62 bcs", "join bcs on the same region");
        return -1;
    }
    bcs.push(Bc { region: r, value: value });
    return bcs.len - 1;
}

fn output_item(n: i64, s: Step&) {
    var o = first(n);
    while o >= 0 {
        if kind(o) == syntax::N_EVERY {
            s.every = round(codegen::eval(first(o)).v[0]) as! i64;
        } else if kind(o) == syntax::N_NAME && compact(o) == v_name {
            s.fields |= 1;
        } else if kind(o) == syntax::N_NAME && compact(o) == p_name {
            s.fields |= 2;
        } else if kind(o) == syntax::N_NAME && compact(o) == a_name {
            s.fields |= 4;
        } else {
            unsupported(tok(o), "the output of the two-phase generator is v, p and alpha", "write output v, p, alpha every 10");
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
            let b = bc(c);
            if b >= 0 { own.push(b); }
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
            unsupported(tok(c), "the step statement is not in the two-phase generator", "use only bc, output, report and expect");
        }
        c = next(c);
    }
    s.first_bc = active.len;
    for g in global_bcs {
        var replaced = false;
        for b in own { if bcs[b].region == bcs[g].region { replaced = true; } }
        if !replaced { active.push(g); }
    }
    active.append(own);
    s.bc_count = active.len - s.first_bc;
    steps.push(s);
}

fn aggregate_of(n: i64) -> i64 {
    var r = n;
    var guard_depth = 0;
    while kind(r) == syntax::N_NAME && resolve(name(r)) >= 0 && guard_depth < 100 {
        r = resolve(name(r));
        guard_depth++;
    }
    let reduction = (is_call(r, "max") || is_call(r, "min")) && arg_count(r) == 2;
    if !reduction || kind(arg(r, 1)) != syntax::N_NAME || find_name(shape_names, name(arg(r, 1))) != domain {
        unsupported(tok(n), str("'", source_of(n), "' is not a max or min of a field over the domain of the physics"),
                    "write let front: Length = max(x*interp(wet.f, wet.s, alpha), tank)");
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
    var mesh_nodes: i64[>..] = [];
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
        } else if k == syntax::N_MESH {
            mesh_nodes.push(c);
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in the two-phase generator"), "remove it");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_nodes.len != 1 { unsupported(at, "the two-phase generator has one physics", PHYSICS_FIX); }
    if !ok() { return false; }
    physics_item(physics_nodes[0]);
    if !ok() { return false; }
    if mesh_nodes.len != 1 { unsupported(at, "the model has one mesh", MESH_FIX); }
    if !ok() { return false; }
    mesh_item(mesh_nodes[0]);
    if ok() { mesh_cells(mesh_nodes[0]); }
    for m in methods { method_item(m); }
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
            let b = bc(c);
            if b >= 0 { global_bcs.push(b); }
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

fn dof_estimate() -> i64 { return 4 * cells[0] * cells[1]; }

recursive fn has_table(n: i64) -> bool {
    if table_of(n) >= 0 { return true; }
    var c = first(n);
    while c >= 0 {
        if has_table(c) { return true; }
        c = next(c);
    }
    return false;
}

fn emit_velocity(out: u8[>..]&, n: i64) {
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

fn emit_functions(out: u8[>..]&) {
    out.append("fn properties(a: f64) -> f64[2] {\n    return [");
    emit(out, rho_node);
    out.append(", ");
    emit(out, mu_node);
    out.append("];\n}\n\nfn initial(x: f64, y: f64) -> f64 {\n    return ");
    if initial_node >= 0 { emit(out, initial_node); } else { out.append("0.0"); }
    out.append(";\n}\n\nfn initial_velocity(x: f64, y: f64) -> f64[2] {\n    let t = 0.0;\n    return ");
    emit_velocity(out, velocity_node);
    out.append(";\n}\n\nfn boundary(k: i64, x: f64, y: f64, t: f64) -> f64[2] {\n");
    for b, k in bcs {
        format(out, "    if k == ", k, " { return ");
        emit_velocity(out, b.value);
        out.append("; }\n");
    }
    out.append("    return [0.0, 0.0];\n}\n\n");
    for c, i in curve_names {
        format(out, "fn select_", c, "() -> bool {\n    return vof::select_faces() { x, y => ");
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
    format(out, "var series_", i, ": f64[>..] = [];\n\nfn field_", i, "(x: f64, y: f64, a: f64, p: f64, vx: f64, vy: f64) -> f64 {\n    let t = 0.0;\n    return ");
    emit(out, arg(r, 0));
    format(out, ";\n}\n\nfn aggregate_", i, "() -> f64 {\n    var m = ", if most { "-1.0e300" } else { "1.0e300" }, ";\n    for c in vof::cell_count() {\n        let q = vof::centre(c);\n");
    format(out, "        m = ", if most { "max" } else { "min" }, "(m, field_", i, "(q[0], q[1], vof::alpha[c], vof::pres[c], vof::vel[2 * c], vof::vel[2 * c + 1]));\n    }\n    return m;\n}\n\n");
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
        format(out, ", k, \".vtu\");\n        output(path, vof::write_vtu(path, ", s.fields, "));\n    }\n");
    }
    format(out, "}\n\n// step ", s.name, ": transient, backward_euler\nfn step_", s.name, "() {\n    vof::clear_bcs();\n");
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        let c = curve_names[bcs[b].region];
        format(out, "    region(select_", c, "(), \"", c, "\");\n    vof::fix_selected(", b, ");\n");
    }
    out.append("    times.clear();\n");
    for r, i in aggregates { format(out, "    series_", i, ".clear();\n"); }
    out.append("    record(vof::time);\n    solved(");
    quoted(out, s.name);
    out.append(", vof::solve_transient(");
    lit(out, s.end);
    out.append(", ");
    lit(out, s.dt);
    format(out, ", properties, boundary, after_", s.name, "));\n");
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
            write_stderr(str("skein run: no boundary face of the mesh is on the region '", name, "'\n"));
            exit(1);
        }
    }

    fn solved(step: u8[:], finished: bool) {
        format(solves, "  ", step, ",", vof::cell_count(), ",", 3 * vof::cell_count(), ",", vof::steps_done, ",", vof::linear_iterations, ",", vof::last_seconds, ",", if finished { "yes" } else { "no" }, "\n");
        solve_count++;
        check("solve", "finite", finished, str(vof::steps_done, " steps, Courant number at most ", vof::max_courant), "each value finite, Courant number at most 1");
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

fn generate(vof_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    out.append("// The values are in the unit system of the model header.\nimport std;\nimport os;\nimport math;\nimport vof from ");
    quoted(out, vof_path);
    out.append(";\nimport solvers from ");
    quoted(out, str(vof_path[..vof_path.len - 3], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    emit_functions(out);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n    // mesh: ");
    format(out, cells[0], " x ", cells[1], " rectangles\n    vof::mesh_rect(");
    lit(out, lo[0]);
    out.append(", ");
    lit(out, lo[1]);
    out.append(", ");
    lit(out, hi[0]);
    out.append(", ");
    lit(out, hi[1]);
    format(out, ", ", cells[0], ", ", cells[1], ");\n    tol = vof::tolerance();\n");
    format(out, "    // physics ", physics_name, ": two fluids, the volume fraction ", a_name, "\n    vof::set_gravity(");
    lit(out, gravity[0]);
    out.append(", ");
    lit(out, gravity[1]);
    out.append(");\n    vof::set_tension(");
    lit(out, if sigma_node >= 0 { scalar_value(sigma_node) } else { 0.0 });
    format(out, ");\n    vof::set_scheme(", FLUX_CONSTANTS[flux], ", ", LIMITER_CONSTANTS[limiter], ");\n    vof::start(initial);\n    vof::start_velocity(initial_velocity);\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
