// The code generator of the ODE models: a Skein syntax tree to a Goose
// program that calls src/ode/ode.goose.
//
// It takes a physics with only global unknowns (spec §10.5) and a method
// `ode`: lets of scalars and vectors, and an `eq` for each unknown, with
// `ddt`, `ddt2`, `dot`, `sin`, `cos`, `sqrt`, `exp`, `ln`, a whole power,
// the time `t` and the unknowns. The `ddt` of any expression is exact: the
// program evaluates the equations on series in time (src/ode/jet.goose).
// One `transient(end:, dt:, scheme:)` step integrates from the initial
// values, and its reports and expects read the unknowns at the end. Each
// other construct gives an error `unsupported` at its token, with a fix.
//
// The values are in the unit system of the header.
import std;
import math;
namespace ode_gen;
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
    kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == callee
}
fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }
fn binary(n: i64, op: u8[:]) -> bool { kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == op }

// Whether this generator takes the model: a physics with no region and no
// template, and an unknown.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS {
            var global = false;
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_ON || kind(u) == syntax::N_USES { global = false; break; }
                if kind(u) == syntax::N_UNKNOWN { global = true; }
                u = next(u);
            }
            if global { return true; }
        }
        c = next(c);
    }
    return false;
}

struct Item {
    name: u8[..64],
    node: i64,        // the value of a let, or the type of an unknown
    count: i64,       // the components of a let
    start: i64,       // the first local of its components
}

struct Report {
    item: i64,
    unknown: i64,
}

struct Expect {
    node: i64,
    unknown: i64,
    step: i64,
}

var lets: Item[>..];            // the lets of the physics, in order
var unknowns: Item[>..];
var initials: f64[>..];         // the initial value of each unknown
var equations: i64[>..];        // the N_EQN nodes
var reports: Report[>..];
var expects: Expect[>..];
var steps: u8[..64][>..];       // the name of the step: this generator has one
var step_end = 0.0;
var step_dt = 0.0;
var step_points = 6;
var step_scheme: u8[..64] = [];
var locals = 0;                 // the components of all the lets

fn reset() {
    syntax::errors.clear();
    codegen::reset();
    lets.clear();
    unknowns.clear();
    initials.clear();
    equations.clear();
    reports.clear();
    expects.clear();
    steps.clear();
    step_end = 0.0;
    step_dt = 0.0;
    step_points = 6;
    step_scheme = [];
    locals = 0;
}

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    let sys = child(n, 1);
    for s, i in codegen::SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            codegen::system = codegen::SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

fn scalar_value(n: i64) -> f64 {
    let v = codegen::eval(n);
    if v.n != 1 { unsupported(tok(n), "the value is a vector", "write a scalar"); }
    return v.v[0];
}

fn find_item(items: Item[:], word: u8[:]) -> i64 {
    for it, i in items { if it.name == word { return i; } }
    return -1;
}

// --- the expressions ---------------------------------------------------------------

let EXPR_FIX = "write + - * / of the unknowns, t, lets and numbers, a whole power, ddt, ddt2, dot, sin, cos, sqrt, exp, ln, (a, b) and v.x";

// Whether the value of n depends on an unknown, on t, on a let of the
// physics or on ddt: it is not a constant of the model.
recursive fn varies(n: i64) -> bool {
    if kind(n) == syntax::N_NAME {
        return name(n) == "t" || find_item(unknowns, name(n)) >= 0 || find_item(lets, name(n)) >= 0;
    }
    if is_call(n, "ddt") || is_call(n, "ddt2") { return true; }
    var c = first(n);
    while c >= 0 {
        if varies(c) { return true; }
        c = next(c);
    }
    return false;
}

// The components of the value of n: 1 for a scalar, the length of a vector.
recursive fn components(n: i64) -> i64 {
    let k = kind(n);
    if k == syntax::N_VECTOR {
        var count = 0;
        var c = first(n);
        while c >= 0 {
            if kind(c) != syntax::N_UNIT { count++; }
            c = next(c);
        }
        return count;
    }
    if k == syntax::N_NAME {
        let l = find_item(lets, name(n));
        return if l >= 0 { lets[l].count } else { 1 };
    }
    if k == syntax::N_UNARY { return components(first(n)); }
    if binary(n, "+") || binary(n, "-") { return max(components(child(n, 0)), components(child(n, 1))); }
    if binary(n, "*") { return max(components(child(n, 0)), components(child(n, 1))); }
    if binary(n, "/") { return components(child(n, 0)); }
    if is_call(n, "ddt") || is_call(n, "ddt2") { return components(arg(n, 0)); }
    return 1;
}

fn axis_of(word: u8[:]) -> i64 { return if word == "x" { 0 } else if word == "y" { 1 } else if word == "z" { 2 } else { -1 }; }

// The Goose expression of component c of n, a series, in `out`.
recursive fn jet(out: u8[>..]&, n: i64, c: i64) {
    let k = kind(n);
    if k == syntax::N_NAME {
        let u = find_item(unknowns, name(n));
        let l = find_item(lets, name(n));
        if u >= 0 {
            format(out, "ode::uj[", u, "]");
        } else if l >= 0 {
            format(out, "L", lets[l].start + min(c, lets[l].count - 1));
        } else if name(n) == "t" {
            out.append("ode::tj");
        } else {
            out.append("ode::jconst(");
            lit(out, scalar_value(n));
            out.append(")");
        }
        return;
    }
    if !varies(n) && k != syntax::N_VECTOR {
        out.append("ode::jconst(");
        lit(out, scalar_value(n));
        out.append(")");
        return;
    }
    if k == syntax::N_VECTOR {
        var i = 0;
        var v = first(n);
        while v >= 0 {
            if kind(v) == syntax::N_UNIT {
                unsupported(tok(v), "a vector with a unit in the physics", "write the unit on each component");
            } else {
                if i == c { jet(out, v, 0); return; }
                i++;
            }
            v = next(v);
        }
        unsupported(tok(n), "the vector has fewer components", "give the vectors one length");
        out.append("ode::jconst(0.0)");
        return;
    }
    jet_operation(out, n, c);
}

// A unary, binary or call node whose value varies.
recursive fn jet_operation(out: u8[>..]&, n: i64, c: i64) {
    let k = kind(n);
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        out.append("ode::jneg(");
        jet(out, first(n), c);
        out.append(")");
    } else if binary(n, "+") || binary(n, "-") {
        out.append(if binary(n, "+") { "ode::jadd(" } else { "ode::jsub(" });
        jet(out, child(n, 0), c);
        out.append(", ");
        jet(out, child(n, 1), c);
        out.append(")");
    } else if binary(n, "*") {
        let a = components(child(n, 0));
        out.append("ode::jmul(");
        jet(out, child(n, 0), if a > 1 { c } else { 0 });
        out.append(", ");
        jet(out, child(n, 1), if a > 1 { 0 } else { c });
        out.append(")");
    } else if binary(n, "/") {
        if components(child(n, 1)) != 1 { unsupported(tok(n), "a division by a vector", "divide by a scalar"); }
        out.append("ode::jdiv(");
        jet(out, child(n, 0), c);
        out.append(", ");
        jet(out, child(n, 1), 0);
        out.append(")");
    } else if binary(n, "^") {
        jet_power(out, n);
    } else if k == syntax::N_MEMBER && axis_of(name(n)) >= 0 {
        jet(out, first(n), axis_of(name(n)));
    } else if k == syntax::N_CALL {
        jet_call(out, n, c);
    } else {
        unsupported(tok(n), str("the expression '", source_of(n), "' is not in this generator"), EXPR_FIX);
        out.append("ode::jconst(0.0)");
    }
}

recursive fn jet_power(out: u8[>..]&, n: i64) {
    let e = child(n, 1);
    let p = if varies(e) { 0.5 } else { scalar_value(e) };
    if p != floor(p) || abs(p) > 16.0 {
        unsupported(tok(e), "the power is not a whole number from -16 to 16", "write a whole exponent, or sqrt(x)");
    }
    out.append("ode::jpow(");
    jet(out, child(n, 0), 0);
    format(out, ", ", p as! i64, ")");
}

recursive fn jet_call(out: u8[>..]&, n: i64, c: i64) {
    for f in ["sin", "cos", "sqrt", "exp", "ln"] {
        if !codegen::call_args(n, f, 1) { continue; }
        format(out, "ode::j", if f == "ln" { "log" } else { f }, "(");
        jet(out, arg(n, 0), 0);
        out.append(")");
        return;
    }
    if codegen::call_args(n, "ddt", 1) || codegen::call_args(n, "ddt2", 1) {
        out.append(if is_call(n, "ddt2") { "ode::jddt(ode::jddt(" } else { "ode::jddt(" });
        jet(out, arg(n, 0), c);
        out.append(if is_call(n, "ddt2") { "))" } else { ")" });
    } else if codegen::call_args(n, "dot", 2) {
        let count = components(arg(n, 0));
        if components(arg(n, 1)) != count { unsupported(tok(n), "dot of two vectors of different lengths", "give the vectors one length"); }
        for i in count - 1 { out.append("ode::jadd("); }
        for i in count {
            if i > 0 { out.append(", "); }
            out.append("ode::jmul(");
            jet(out, arg(n, 0), i);
            out.append(", ");
            jet(out, arg(n, 1), i);
            out.append(")");
            if i > 0 { out.append(")"); }
        }
    } else {
        unsupported(tok(n), str("the function '", source_of(child(n, 0)), "' is not in this generator"), EXPR_FIX);
        out.append("ode::jconst(0.0)");
    }
}

// --- the model ---------------------------------------------------------------------

let PHYSICS_FIX = "write physics linkage solves q: Angle = 0 rad { let ... eq ... } and method linkage: ode {}";
let STEP_FIX = "write step run: transient(end: 1 s, dt: 1 ms, scheme: bdf2) { report q }";

fn physics_item(n: i64) {
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_UNKNOWN {
            var w: Item = Item { name: [], node: first(c), count: 1, start: 0 };
            set64(w.name, name(c));
            unknowns.push(w);
            initials.push(if next(first(c)) >= 0 { scalar_value(next(first(c))) } else { 0.0 });
        } else if k == syntax::N_LET {
            var w: Item = Item { name: [], node: child(c, 1), count: components(child(c, 1)), start: locals };
            set64(w.name, name(c));
            locals += w.count;
            lets.push(w);
        } else if k == syntax::N_EQN {
            equations.push(c);
        } else {
            unsupported(tok(c), str("the statement '", syntax::kind_name(k), "' of a physics of global unknowns is not in this generator"),
                        PHYSICS_FIX);
        }
        c = next(c);
    }
    if equations.len != unknowns.len {
        unsupported(tok(n), str("the physics has ", unknowns.len, " unknowns and ", equations.len, " eqs"), "write one eq for each unknown");
    }
    if unknowns.len > 16 { unsupported(tok(n), "more than 16 unknowns", "use at most 16 unknowns"); }
}

fn method_item(n: i64) {
    let kind_tok = syntax::nodes[n].aux;
    if kind_tok < 0 || syntax::text(kind_tok) != "ode" {
        unsupported(tok(n), "the method of a physics of global unknowns is ode", PHYSICS_FIX);
    }
}

// The points of the polynomial of the scheme (src/ode/ode.goose).
fn scheme(v: i64) {
    let word = if kind(v) == syntax::N_CALL { name(child(v, 0)) } else { source_of(v) };
    set64(step_scheme, word);
    if word == "backward_euler" {
        step_points = 2;
    } else if word == "bdf2" {
        step_points = 3;
    } else if word == "generalized_alpha" {
        step_points = 6;
    } else {
        unsupported(tok(v), "this generator has the schemes backward_euler, bdf2 and generalized_alpha", STEP_FIX);
    }
}

fn step(n: i64) {
    if steps.len > 0 { unsupported(tok(n), "a second step: this generator has one", "remove the second step"); return; }
    let proc = child(n, 0);
    if name(proc) != "transient" { unsupported(tok(proc), "this generator has a transient step", STEP_FIX); return; }
    var w: u8[..64] = [];
    set64(w, name(n));
    steps.push(w);
    var a = first(proc);
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        let word = if named >= 0 { syntax::text(named) } else { "" };
        if word == "end" {
            step_end = scalar_value(first(a));
        } else if word == "dt" {
            step_dt = scalar_value(first(a));
        } else if word == "scheme" {
            scheme(first(a));
        } else {
            unsupported(tok(a), "the argument of transient is not in this generator", STEP_FIX);
        }
        a = next(a);
    }
    let count = if step_dt > 0.0 { step_end / step_dt } else { 0.0 };
    if !(count >= 1.0) || abs(count - round(count)) > 1.0e-9 * count || count > 1.0e7 {
        unsupported(tok(proc), "the step needs dt > 0 and end an whole number of dt, at most 10 million", STEP_FIX);
    }
    var c = next(proc);
    while c >= 0 {
        if kind(c) == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                reports.push(Report { item: r, unknown: unknown_of(first(r)) });
                r = next(r);
            }
        } else if kind(c) == syntax::N_EXPECT {
            expect(c, 0);
        } else {
            unsupported(tok(c), "the step statement is not in this generator", "use only report and expect");
        }
        c = next(c);
    }
}

// The unknown that n names, or -1 and an error.
fn unknown_of(n: i64) -> i64 {
    let u = if kind(n) == syntax::N_NAME { find_item(unknowns, name(n)) } else { -1 };
    if u < 0 { unsupported(tok(n), str("the value '", source_of(n), "' is not an unknown"), "write the name of an unknown, for example report xb"); }
    return u;
}

fn expect(n: i64, step: i64) {
    expects.push(Expect { node: n, unknown: unknown_of(child(n, 0)), step: step });
}

fn read(prog: i64) -> bool {
    reset();
    var physics = false;
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            codegen::value_names.push(name_of(c));
            codegen::value_nodes.push(child(c, 1));
        } else if k == syntax::N_PHYSICS && !physics {
            physics = true;
            physics_item(c);
        } else if k == syntax::N_METHOD {
            method_item(c);
        } else if k != syntax::N_SCALES && k != syntax::N_STEP && k != syntax::N_EXPECT && k != syntax::N_ALLOW {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this generator"),
                        "remove it: this generator runs one physics of global unknowns");
        }
        c = next(c);
    }
    if !ok() { return false; }
    c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_STEP { step(c); }
        c = next(c);
    }
    c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_EXPECT { expect(c, -1); }
        c = next(c);
    }
    if steps.len == 0 { unsupported(tok(first(prog)), "the model has no step", STEP_FIX); }
    return ok();
}

fn name_of(n: i64) -> u8[..64] {
    var w: u8[..64] = [];
    set64(w, name(n));
    return w;
}

fn dof_estimate() -> i64 { return unknowns.len; }

// --- the Goose program -------------------------------------------------------------

let PRELUDE = """
    // The reports and the verification records, in TOON (command_line.md §6).
    var reports: u8[>..] = [];
    var report_count = 0;
    var checks: u8[>..] = [];
    var check_count = 0;
    var solves: u8[>..] = [];
    var solve_count = 0;
    // A fault of the model or of a solve.
    var failed = false;
    // A failed check.
    var missed = false;

    // Appends `s` as a TOON value, in quotes when it needs them.
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
        if !pass { missed = true; }
    }

    fn solved(step: u8[:], regular: bool) {
        format(solves, "  ", step, ",", ode::count, ",", ode::solves, ",", ode::newton_steps, ",", ode::last_seconds, ",", if regular { "yes" } else { "no" }, "\n");
        solve_count++;
        if !regular {
            write_stderr(str("skein run: the Newton iteration of step '", step, "' does not converge at t = ", ode::time(), "\n"));
            failed = true;
        }
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", abs(got - want) <= t, str(name, ": ", quantity(got, unit)), limit);
    }

    fn print_records() {
        var s: u8[>..] = [];
        format(s, "report[", report_count, "]{step,expr,value}:\n", reports);
        format(s, "verify[", check_count, "]{stage,id,result,value,limit}:\n", checks);
        format(s, "solve[", solve_count, "]{step,unknowns,steps,newton_steps,seconds,converged}:\n", solves);
        s.append("results[0]:\n");
        s.append(solvers::memory_records());
        write_stdout(s);
    }
""";

// The unit of the type name t: a base unit of the header, or "".
fn base_unit(t: i64) -> const u8[:] {
    guard kind(t) == syntax::N_TYPE_NAME else { return ""; };
    let w = name(t);
    if w == "Length" { return codegen::system.length[..]; }
    if w == "Force" { return codegen::system.force[..]; }
    if w == "Angle" { return "rad"; }
    return "";
}

// The unit of a type: a base unit, or a product of two, or "".
fn type_unit(t: i64) -> u8[>..] {
    if kind(t) == syntax::N_DIM_MUL {
        let a = base_unit(child(t, 0));
        let b = base_unit(child(t, 1));
        return if a.len > 0 && b.len > 0 { str(a, "*", b) } else { [] };
    }
    return str(base_unit(t));
}

// The locals of the lets, as globals.
fn emit_lets(out: u8[>..]&) {
    for i in locals { format(out, "var L", i, ": f64[5] = [0.0, 0.0, 0.0, 0.0, 0.0];\n"); }
    out.push('\n');
    out.append("fn eval_lets() {\n");
    for l in lets {
        for c in l.count {
            format(out, "    L", l.start + c, " = ");
            jet(out, l.node, c);
            out.append(";\n");
        }
    }
    out.append("}\n\n");
}

// The residual of each eq at the series of the candidate.
fn emit_residual(out: u8[>..]&) {
    out.append("fn residual(t: f64) {\n    eval_lets();\n");
    for e, k in equations {
        format(out, "    ode::resid[", k, "] = ode::jsub(");
        jet(out, child(e, 0), 0);
        out.append(", ");
        jet(out, child(e, 1), 0);
        out.append(")[0];\n");
    }
    out.append("}\n\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    let within = child(n, 2);
    let right = child(n, 1);
    let unit = if kind(right) == syntax::N_NUMBER { syntax::child_of_kind(right, syntax::N_UNIT) } else { -1 };
    format(out, "    // expect ", name(n), "\n    {\n        let name = ");
    quoted(out, name(n));
    format(out, ";\n        expect(name, ode::value(", x.unknown, "), ");
    lit(out, scalar_value(right));
    out.append(", ");
    lit(out, scalar_value(within));
    format(out, ", ", if kind(within) == syntax::N_PERCENT { "true" } else { "false" }, ", ");
    quoted(out, if unit >= 0 { source_of(unit) } else { "" });
    out.append(", ");
    quoted(out, str(source_of(right), " within ", source_of(within)));
    out.append(");\n    }\n");
}

fn emit_step(out: u8[>..]&) {
    let s = steps[0];
    let count = round(step_end / step_dt) as! i64;
    format(out, "// step ", s, "\nfn step_", s, "() {\n    ode::begin(", unknowns.len, ", ");
    lit(out, step_dt);
    format(out, ", ", step_points, ", ", count, ");\n");
    for u, i in unknowns {
        format(out, "    ode::initial(", i, ", ");
        lit(out, initials[i]);
        out.append(");\n");
    }
    format(out, "    var regular = true;\n    for n in 1..", count + 1, " {\n");
    out.append("        let ok = ode::advance(n) { t => residual(t); };\n        if !ok { regular = false; break; }\n    }\n");
    out.append("    solved(");
    quoted(out, s);
    out.append(", regular);\n");
    for r in reports {
        if r.unknown < 0 { continue; }
        out.append("    {\n        let name = ");
        quoted(out, source_of(r.item));
        format(out, ";\n        report(");
        quoted(out, s);
        format(out, ", name, quantity(ode::value(", r.unknown, "), ");
        quoted(out, type_unit(unknowns[r.unknown].node));
        out.append("));\n    }\n");
    }
    for x in expects { if x.step == 0 && x.unknown >= 0 { emit_expect(out, x); } }
    out.append("}\n\n");
}

// Writes the Goose program of the model that `read` read. The program
// imports the ode package from `ode_path`.
fn generate(ode_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    let src = ode_path[..ode_path.len - 3];
    out.append("// The values are in the unit system of the model header.\nimport std;\nimport os;\nimport math;\nimport ode from ");
    quoted(out, str(src, "ode"));
    out.append(";\nimport solvers from ");
    quoted(out, str(src, "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    out.push('\n');
    emit_lets(out);
    emit_residual(out);
    emit_step(out);
    format(out, "fn main() {\n    step_", steps[0], "();\n");
    for x in expects { if x.step < 0 && x.unknown >= 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed || missed { exit(1); }\n}\n");
}
