// The part of the simplex generator for the edge physics of
// `family: nedelec` and the Darcy flow of `family: raviart_thomas` (spec
// §13), on P1 cells (src/elements/edge.goose and src/elements/mixed.goose).
//
// The edge physics: `let h = curl(u)/X` (or `X*curl(u)`, or `curl(u)`) and
// `eq curl(h) + c*u = 0` for a vector u, Complex in `harmonic` steps. Its
// bcs: `cross(n, u) = cross(n, g)` or `u = 0` (essential), and
// `cross(n, h) = z*u` (natural). `gauge: tree` holds a spanning tree of the
// cells with c = 0.
//
// The Darcy flow: `eq q + k*grad(p) = 0` and `eq div(q) = f`, in `static`
// steps. Its bcs: `p = v` (natural in the mixed form), and `dot(q, n) = g`
// (essential).
//
// A report or an expect is a Complex expression of constants, of points
// (`hc.x at p`) and of line integrals (`integral(ef.x, c)` on a `line`), and
// of `far_field(u, H, s, d)` on a boundary, with H a constant times curl(u),
// and for the Darcy flow of `integral(dot(q, n), s)` on a boundary.
import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import .gen;

let EDGE = 3;       // as in src/elements/assemble.goose
let DARCY = 4;      // MIXED there

let EDGE_FIX = "write let h: Vector<...> = curl(u)/mu and eq curl(h) + c*u = 0, with family: nedelec";
let DARCY_FIX = "write eq q + k*grad(p) = 0 and eq div(q) = f, with family: raviart_thomas";
let PROBE_FIX = "write re, imag, abs or arg of + - * / of constants, F.x at a point, integral(F.x, c) on a line c, integral(dot(q, n), s), or far_field(u, hf*curl(u), s, d)";

var family_word: u8[..64];          // "nedelec", "raviart_thomas" or ""
var curl_name: u8[..64];            // h of the edge physics
var curl_coef = -1;                 // X of h = curl(u)/X or X*curl(u), or -1
var curl_divides = false;
var mass_terms: i64[>..];           // the coefficients of the terms c*u of the eq
var mass_signs: f64[>..];
var gauge_tree = false;
var darcy_k = -1;                   // k of q + k*grad(p) = 0
var darcy_f = -1;                   // f of div(q) = f, or -1 for 0
var physics_lets: i64[>..];         // the let items of the physics
var probes: EdgeExpect[>..];        // the expects of the edge and Darcy physics
var quiet = false;                  // a probe that is not in this generator sets probe_bad, and gives no error
var probe_bad = false;

struct EdgeExpect {
    name: u8[..64],
    expr: i64,
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
    step: i64,
}

fn edge_reset() {
    family_word = [];
    curl_name = [];
    curl_coef = -1;
    curl_divides = false;
    mass_terms.clear();
    mass_signs.clear();
    gauge_tree = false;
    darcy_k = -1;
    darcy_f = -1;
    physics_lets.clear();
    probes.clear();
    quiet = false;
    probe_bad = false;
}

fn probe_fail(at: i64, msg: u8[:], fix: u8[:]) {
    if quiet { probe_bad = true; } else { unsupported(at, msg, fix); }
}

// The family of the methods: nedelec and raviart_thomas need this part of
// the generator. A model with one goes to this generator (src/codegen/run.goose).
fn family_of(methods: i64[:]) -> const u8[:] {
    for m in methods {
        var s = first(child(m, 0));
        while s >= 0 {
            if kind(s) == syntax::N_SETTING && name(s) == "family" && kind(first(s)) == syntax::N_NAME { return name(first(s)); }
            s = next(s);
        }
    }
    return "";
}

fn family_claims(prog: i64) -> bool {
    var methods: i64[>..] = [];
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_METHOD { methods.push(c); }
        c = next(c);
    }
    let f = family_of(methods);
    return f == "nedelec" || f == "raviart_thomas";
}

fn is_vector_type(n: i64) -> bool { return kind(n) == syntax::N_TYPE_APPLY && syntax::text(tok(n)) == "Vector"; }

// The physics block of the edge physics or of the Darcy flow.
fn edge_physics(n: i64) {
    physics = if family_word == "nedelec" { EDGE } else { DARCY };
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            body(first(c));
        } else if k == syntax::N_UNKNOWN {
            edge_unknown(c);
        } else if k == syntax::N_MATASSIGN {
            materials_of(c);
        } else if k == syntax::N_LET {
            physics_lets.push(c);
            if physics == EDGE { curl_let(c); }
        } else if k == syntax::N_EQN {
            if physics == EDGE { curl_eq(c); } else { darcy_eq(c); }
        } else {
            unsupported(tok(c), "the physics statement is not in this generator", if physics == EDGE { EDGE_FIX } else { DARCY_FIX });
        }
        c = next(c);
    }
    if physics == EDGE && (curl_name.len == 0 || unknown_name.len == 0) { unsupported(tok(n), "the edge physics needs a vector unknown, h and its eq", EDGE_FIX); }
    if physics == DARCY && (darcy_k < 0 || flux_name.len == 0 || unknown_name.len == 0) {
        unsupported(tok(n), "the Darcy flow needs a vector q, a scalar p and the two eqs", DARCY_FIX);
    }
}

// The unknown u of the edge physics, or q (flux_name) and p (unknown_name)
// of the Darcy flow.
fn edge_unknown(c: i64) {
    let vector = is_vector_type(child(c, 0)) ||
        kind(child(c, 0)) == syntax::N_TYPE_APPLY && syntax::text(tok(child(c, 0))) == "Complex" && is_vector_type(first(child(c, 0)));
    if physics == DARCY && vector && flux_name.len == 0 {
        set64(flux_name, name(c));
    } else if physics == DARCY && !vector && unknown_name.len == 0 {
        set64(unknown_name, name(c));
    } else if physics == EDGE && vector && unknown_name.len == 0 {
        set64(unknown_name, name(c));
    } else {
        unsupported(tok(c), "the unknowns are one vector u of nedelec, or a vector q and a scalar p of raviart_thomas",
                    if physics == EDGE { EDGE_FIX } else { DARCY_FIX });
    }
}

// `let h = curl(u)/X`, `X*curl(u)` or `curl(u)`. Another let is a value of
// the reports.
fn curl_let(c: i64) {
    let v = child(c, 1);
    let cu = str("curl(", unknown_name, ")");
    guard curl_name.len == 0 else { return; };
    if source_of(v) == cu {
        set64(curl_name, name(c));
    } else if kind(v) == syntax::N_BINARY && (syntax::text(tok(v)) == "/" || syntax::text(tok(v)) == "*") {
        let d = syntax::text(tok(v)) == "/";
        if source_of(child(v, 0)) == cu {
            curl_coef = child(v, 1);
            curl_divides = d;
            set64(curl_name, name(c));
        } else if !d && source_of(child(v, 1)) == cu {
            curl_coef = child(v, 0);
            set64(curl_name, name(c));
        }
    }
}

// eq curl(h) + c1*u - c2*u ... = 0.
fn curl_eq(c: i64) {
    var left = child(c, 0);
    let fix = EDGE_FIX;
    if source_of(child(c, 1)) != "0" { unsupported(tok(child(c, 1)), "the edge physics has no source in this generator", fix); }
    while kind(left) == syntax::N_BINARY && (syntax::text(tok(left)) == "+" || syntax::text(tok(left)) == "-") {
        let t = child(left, 1);
        guard kind(t) == syntax::N_BINARY && syntax::text(tok(t)) == "*" && source_of(child(t, 1)) == unknown_name else {
            unsupported(tok(t), str("the term '", source_of(t), "' is not c*", unknown_name), fix);
            return;
        };
        mass_terms.push(child(t, 0));
        mass_signs.push(if syntax::text(tok(left)) == "-" { -1.0 } else { 1.0 });
        left = child(left, 0);
    }
    if curl_name.len == 0 || source_of(left) != str("curl(", curl_name, ")") {
        unsupported(tok(left), "the eq of the edge physics starts with curl(h) of the let h = curl(u)/X", fix);
    }
}

// eq q + k*grad(p) = 0, or eq div(q) = f.
fn darcy_eq(c: i64) {
    let left = child(c, 0);
    let right = child(c, 1);
    if source_of(left) == str("div(", flux_name, ")") {
        if source_of(right) != "0" { darcy_f = right; }
        return;
    }
    let g = str("grad(", unknown_name, ")");
    let ok = kind(left) == syntax::N_BINARY && syntax::text(tok(left)) == "+" && source_of(child(left, 0)) == flux_name &&
        kind(child(left, 1)) == syntax::N_BINARY && syntax::text(tok(child(left, 1))) == "*" && source_of(child(child(left, 1), 1)) == g &&
        source_of(right) == "0";
    if !ok {
        unsupported(tok(c), "this generator solves eq q + k*grad(p) = 0 and eq div(q) = f", DARCY_FIX);
        return;
    }
    darcy_k = child(child(left, 1), 0);
}

// The settings family, order and gauge of the method.
fn edge_method(s: i64) -> bool {
    guard kind(s) == syntax::N_SETTING else { return false; };
    if name(s) == "family" { return true; }
    if name(s) == "gauge" && physics == EDGE && name(first(s)) == "tree" {
        gauge_tree = true;
        return true;
    }
    return false;
}

// `cross(n, a)`: the node of a, else -1.
fn crossed(n: i64) -> i64 {
    guard is_call(n, "cross") && syntax::child_count(n) == 3 && source_of(arg(n, 0)) == "n" else { return -1; };
    return arg(n, 1);
}

// A bc of the edge physics or of the Darcy flow on boundary s.
fn edge_bc(n: i64, s: i64) -> i64 {
    let target = child(n, 1);
    let v = child(n, 2);
    var b = Bc { surf: s, kind: BC_FIX, comp: -1, v: [0.0, 0.0, 0.0], tok: tok(n), node: -1 };
    let t = source_of(target);
    let fix = if physics == EDGE { str("write bc wall: cross(n, ", unknown_name, ") = 0, cross(n, ", unknown_name, ") = cross(n, g), or cross(n, ",
                                       curl_name, ") = z*", unknown_name) }
              else { str("write bc inlet: ", unknown_name, " = 1 Pa, or bc wall: dot(", flux_name, ", n) = 0 m/s") };
    let tangent = crossed(target) >= 0 && source_of(crossed(target)) == unknown_name || t == unknown_name;
    if physics == EDGE && tangent && source_of(v) == "0" {
        b.node = -1;
    } else if physics == EDGE && tangent && crossed(v) >= 0 {
        b.node = crossed(v);
    } else if physics == EDGE && crossed(target) >= 0 && source_of(crossed(target)) == curl_name && kind(v) == syntax::N_BINARY &&
              syntax::text(tok(v)) == "*" && source_of(child(v, 1)) == unknown_name {
        b.kind = BC_LOAD;
        b.node = child(v, 0);
    } else if physics == DARCY && t == unknown_name {
        b.kind = BC_LOAD;
        b.node = v;
    } else if physics == DARCY && (t == str("dot(", flux_name, ", n)") || t == str(flux_name, "*n")) {
        b.node = v;
    } else {
        unsupported(tok(target), str("the bc target '", t, "' is not in this generator"), fix);
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

// --- reports and expects ---------------------------------------------------------

// Whether n holds a field: an `at`, an integral or far_field.
recursive fn has_field(n: i64) -> bool {
    if kind(n) == syntax::N_AT || is_call(n, "integral") || is_call(n, "far_field") { return true; }
    var c = first(n);
    while c >= 0 {
        if has_field(c) { return true; }
        c = next(c);
    }
    return false;
}

let CPARTS = ["re", "imag", "abs", "arg", "conj"];

// The physics let named `word`, or -1.
fn physics_let(word: u8[:]) -> i64 {
    for l in physics_lets { if name(l) == word { return child(l, 1); } }
    return -1;
}

// Writes the Goose expression of the Complex value of n, as f64[2]. The
// statements that select a boundary go to `pre`.
recursive fn emit_cx(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if !has_field(n) {
        let z = cval(n);
        out.push('[');
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.push(']');
    } else if k == syntax::N_BINARY && (op == "+" || op == "-" || op == "*" || op == "/") {
        out.append(if op == "+" { "cadd(" } else if op == "-" { "csub(" } else if op == "*" { "cmul(" } else { "cdiv(" });
        emit_cx(out, pre, child(n, 0));
        out.append(", ");
        emit_cx(out, pre, child(n, 1));
        out.push(')');
    } else if k == syntax::N_UNARY && op == "-" {
        out.append("csub([0.0, 0.0], ");
        emit_cx(out, pre, first(n));
        out.push(')');
    } else if k == syntax::N_CALL && syntax::child_count(n) == 2 && kind(child(n, 0)) == syntax::N_NAME && position(CPARTS, name(child(n, 0))) >= 0 {
        format(out, "cpart(", position(CPARTS, name(child(n, 0))), ", ");
        emit_cx(out, pre, arg(n, 0));
        out.push(')');
    } else if k == syntax::N_AT {
        let p = value(child(n, 1));
        if p.n != dim { probe_fail(tok(n), "the point is not a vector", "write at (x, y, z) mm"); }
        emit_field(out, child(n, 0), p.v, p.v);
    } else if is_call(n, "integral") && syntax::child_count(n) == 3 {
        emit_integral(out, pre, n);
    } else if is_call(n, "far_field") && syntax::child_count(n) == 5 {
        emit_far_field(out, pre, n);
    } else {
        probe_fail(tok(n), str("the expression '", source_of(n), "' is not in this generator"), PROBE_FIX);
    }
}

// integral(F.x, c) on a line c, or integral(dot(q, n), s) on a boundary s.
fn emit_integral(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let f = arg(n, 0);
    let r = arg(n, 1);
    let s = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
    if physics == DARCY && s >= 0 && source_of(f) == str("dot(", flux_name, ", n)") {
        emit_region(pre, s);
        format(pre, "    let flux_", s, " = elements::boundary_flux();\n");
        format(out, "[flux_", s, ", 0.0]");
        return;
    }
    let c = if kind(r) == syntax::N_NAME { codegen::find_value(name(r)) } else { -1 };
    guard c >= 0 && is_call(c, "line") && syntax::child_count(c) == 3 else {
        probe_fail(tok(r), "the integral is on a line that a let names", PROBE_FIX);
        return;
    };
    emit_field(out, f, value(arg(c, 0)).v, value(arg(c, 1)).v);
}

// far_field(u, H, s, d) of the edge physics on a boundary s, with
// H = hf*curl(u), into the medium of the physics: curl(curl(u)/mu) -
// omega^2*eps*u = 0 has nu = 1/mu and c = -omega^2*eps, so k = sqrt(-c/nu)
// and eta = omega*mu/k. The generator checks the medium only at an omega
// of a step.
fn emit_far_field(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let r = arg(n, 2);
    let s = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
    let hf, ok = curl_factor(arg(n, 1));
    let nu, c = edge_terms();
    guard physics == EDGE && source_of(arg(n, 0)) == unknown_name && ok && s >= 0 && (omega_now == 0.0 || c[0] < 0.0 && c[1] == 0.0) else {
        probe_fail(tok(n), "far_field needs E = u, H = a constant times curl(u), a boundary of the body, and a medium with no loss", PROBE_FIX);
        return;
    };
    let k = sqrt(max(-c[0] / nu, 0.0));
    emit_region(pre, s);
    format(pre, "    let far_", n, " = elements::far_field(");
    lits(pre, value(arg(n, 3)).v);
    pre.append(", ");
    lit(pre, k);
    pre.append(", ");
    lit(pre, if k > 0.0 { omega_now / (nu * k) } else { 0.0 });
    pre.append(", ");
    lit(pre, hf[0]);
    pre.append(", ");
    lit(pre, hf[1]);
    pre.append(");\n");
    format(out, "[far_", n, ", 0.0]");
}

// The constant hf of a field H = hf*curl(u): curl(u), a let of the physics,
// and - * / of these with constants. False for another field.
recursive fn curl_factor(n: i64) -> f64[2], bool {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if source_of(n) == str("curl(", unknown_name, ")") { return [1.0, 0.0], true; }
    if k == syntax::N_NAME && physics_let(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let f, ok = curl_factor(physics_let(name(n)));
        codegen::depth--;
        return f, ok;
    }
    if k == syntax::N_UNARY && op == "-" {
        let f, ok = curl_factor(first(n));
        return [-f[0], -f[1]], ok;
    }
    if op == "*" && !has_unknown(child(n, 0)) {
        let f, ok = curl_factor(child(n, 1));
        return cmul(cval(child(n, 0)), f), ok;
    }
    if (op == "*" || op == "/") && !has_unknown(child(n, 1)) {
        let f, ok = curl_factor(child(n, 0));
        let z = cval(child(n, 1));
        return (if op == "*" { cmul(f, z) } else { cdiv(f, z) }), ok;
    }
    return [0.0, 0.0], false;
}

// Writes component F.x of a field at the point a (a == b), or its integral
// on the line from a to b. F is u, curl(u), a let of the physics, or + - *
// and / of these with constants.
recursive fn emit_field(out: u8[>..]&, f: i64, a: f64[3], b: f64[3]) {
    if kind(f) == syntax::N_NAME && name(f) == unknown_name && physics == DARCY {
        emit_probe_sample(out, "darcy_", "", 3, a, b);
        return;
    }
    guard kind(f) == syntax::N_MEMBER && axis_of(name(f)) >= 0 else {
        probe_fail(tok(f), str("'", source_of(f), "' is not a component of a field of the physics"), PROBE_FIX);
        return;
    };
    emit_part(out, first(f), axis_of(name(f)), a, b);
}

recursive fn emit_part(out: u8[>..]&, f: i64, comp: i64, a: f64[3], b: f64[3]) {
    let k = kind(f);
    let op = if k == syntax::N_BINARY { syntax::text(tok(f)) } else { "" };
    if k == syntax::N_NAME && physics == EDGE && name(f) == unknown_name {
        emit_probe_sample(out, "edge_", "false, ", comp, a, b);
    } else if k == syntax::N_NAME && physics == DARCY && name(f) == flux_name {
        emit_probe_sample(out, "darcy_", "", comp, a, b);
    } else if physics == EDGE && source_of(f) == str("curl(", unknown_name, ")") {
        emit_probe_sample(out, "edge_", "true, ", comp, a, b);
    } else if k == syntax::N_NAME && physics_let(name(f)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        emit_part(out, physics_let(name(f)), comp, a, b);
        codegen::depth--;
    } else if (op == "+" || op == "-") || (op == "*" || op == "/") && !has_unknown(child(f, 1)) || op == "*" && !has_unknown(child(f, 0)) {
        out.append(if op == "+" { "cadd(" } else if op == "-" { "csub(" } else if op == "*" { "cmul(" } else { "cdiv(" });
        for i in 2 {
            if i == 1 { out.append(", "); }
            let c = child(f, i);
            if has_unknown(c) {
                emit_part(out, c, comp, a, b);
            } else {
                let z = cval(c);
                out.push('[');
                lit(out, z[0]);
                out.append(", ");
                lit(out, z[1]);
                out.push(']');
            }
        }
        out.push(')');
    } else {
        probe_fail(tok(f), str("'", source_of(f), "' is not a field of the physics in this generator"), PROBE_FIX);
    }
}

// Whether n names the unknown, the flux or a let of the physics.
recursive fn has_unknown(n: i64) -> bool {
    if kind(n) == syntax::N_NAME {
        return name(n) == unknown_name || name(n) == flux_name || physics_let(name(n)) >= 0;
    }
    var c = first(n);
    while c >= 0 {
        if has_unknown(c) { return true; }
        c = next(c);
    }
    return false;
}

// `edge_point(curl, comp, x, y, z)` or `edge_segment(curl, comp, a, b)`.
fn emit_probe_sample(out: u8[>..]&, kind_name: u8[:], curl: u8[:], comp: i64, a: f64[3], b: f64[3]) {
    let line = a[0] != b[0] || a[1] != b[1] || a[2] != b[2];
    format(out, kind_name, if line { "segment(" } else { "point(" }, curl, comp, ", ");
    lits(out, a);
    if line {
        out.append(", ");
        lits(out, b);
    }
    out.push(')');
}

fn edge_report(item: i64) {
    var probe: u8[>..] = [];
    var pre: u8[>..] = [];
    emit_cx(probe, pre, first(item));
    var r = Report { expr: first(item), label: [] };
    codegen::set256(r.label, one_line(item));
    reports.push(r);
}

fn edge_expect(n: i64, step: i64) {
    let left = child(n, 0);
    let right = child(n, 1);
    let within = child(n, 2);
    var probe: u8[>..] = [];
    var pre: u8[>..] = [];
    quiet = true;
    probe_bad = false;
    emit_cx(probe, pre, left);
    quiet = false;
    // An expect that this generator cannot evaluate stays, and the report
    // says that it did not run.
    if probe_bad {
        var s: u8[..256] = [];
        codegen::set256(s, str(name(n), ": ", one_line(left), " is not in this generator"));
        skipped.push(s);
        return;
    }
    var x = EdgeExpect { name: [], expr: left, want: scalar_value(right), within: scalar_value(within),
                         relative: kind(within) == syntax::N_PERCENT, label: [], step: step };
    set64(x.name, name(n));
    codegen::set256(x.label, str(source_of(right), " within ", source_of(within)));
    probes.push(x);
}

// --- the program -------------------------------------------------------------------

let EDGE_PRELUDE = """
    // The Complex values of the reports and expects of the edge and Darcy
    // physics, as [re, im], and whether each point is in the mesh.
    var found_all = true;

    fn cadd(a: f64[2], b: f64[2]) -> f64[2] { [a[0] + b[0], a[1] + b[1]] }
    fn csub(a: f64[2], b: f64[2]) -> f64[2] { [a[0] - b[0], a[1] - b[1]] }
    fn cmul(a: f64[2], b: f64[2]) -> f64[2] { [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]] }

    fn cdiv(a: f64[2], b: f64[2]) -> f64[2] {
        let d = b[0] * b[0] + b[1] * b[1];
        return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
    }

    // re, imag, abs, arg or conj.
    fn cpart(part: i64, a: f64[2]) -> f64[2] {
        if part == 0 { return [a[0], 0.0]; }
        if part == 1 { return [a[1], 0.0]; }
        if part == 2 { return [sqrt(a[0] * a[0] + a[1] * a[1]), 0.0]; }
        if part == 3 { return [atan2(a[1], a[0]), 0.0]; }
        return [a[0], -a[1]];
    }

    fn cquantity(v: f64[2]) -> u8[>..] {
        if v[1] == 0.0 { return str(v[0]); }
        return str(v[0], if v[1] < 0.0 { " - " } else { " + " }, abs(v[1]), "i");
    }

    // Component comp of u, or of curl(u) with curl, at a point.
    fn edge_point(curl: bool, comp: i64, x: f64, y: f64, z: f64) -> f64[2] {
        let re, found_re = elements::edge_at(curl, 0, comp, x, y, z);
        let im, found_im = elements::edge_at(curl, 1, comp, x, y, z);
        if !found_re || !found_im { found_all = false; }
        return [re, im];
    }

    // Its integral on the line from (x0, y0, z0) to (x1, y1, z1).
    fn edge_segment(curl: bool, comp: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> f64[2] {
        let re, im, found = elements::edge_line_integral(curl, comp, x0, y0, z0, x1, y1, z1, 64);
        if !found { found_all = false; }
        return [re, im];
    }

    // Component comp of q (0 to 2), or p (3), of the Darcy flow at a point.
    fn darcy_point(comp: i64, x: f64, y: f64, z: f64) -> f64[2] {
        let v, found = elements::darcy_at(comp, x, y, z);
        if !found { found_all = false; }
        return [v, 0.0];
    }

    // Its integral on a line, by the midpoint rule on 64 pieces.
    fn darcy_segment(comp: i64, x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> f64[2] {
        let length = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
        var s = 0.0;
        for i in 64 {
            let t = (i as f64 + 0.5) / 64.0;
            s += darcy_point(comp, x0 + t * (x1 - x0), y0 + t * (y1 - y0), z0 + t * (z1 - z0))[0];
        }
        return [s * length / 64.0, 0.0];
    }

    fn output_darcy(path: u8[:], p: u8[:], q: u8[:]) {
        if !elements::write_darcy(path, p, q) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    """;

// The Goose expression of a real function of the point p (f64[3]):
// component j of a vector n, or the scalar n with j = -1.
recursive fn emit_space(out: u8[>..]&, n: i64, j: i64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if !mentions_xyz(n) {
        let v = value(n);
        out.push('(');
        lit(out, if j >= 0 && v.n > 1 { v.v[j] } else { v.v[0] });
        out.push(')');
    } else if k == syntax::N_NAME && axis_of(name(n)) >= 0 && j < 0 {
        format(out, "p[", axis_of(name(n)), "]");
    } else if k == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        emit_space(out, codegen::find_value(name(n)), j);
        codegen::depth--;
    } else if k == syntax::N_VECTOR && j >= 0 {
        var i = 0;
        var c = first(n);
        while c >= 0 && (kind(c) == syntax::N_UNIT || i < j) {
            if kind(c) != syntax::N_UNIT { i++; }
            c = next(c);
        }
        let unit = syntax::child_of_kind(n, syntax::N_UNIT);
        out.push('(');
        lit(out, if unit >= 0 { codegen::to_system(unit) } else { 1.0 });
        out.append(" * ");
        emit_space(out, c, -1);
        out.push(')');
    } else if k == syntax::N_BINARY && (op == "+" || op == "-" || op == "*" || op == "/") {
        out.push('(');
        emit_space(out, child(n, 0), if j >= 0 && is_vector(child(n, 0)) { j } else { -1 });
        format(out, " ", op, " ");
        emit_space(out, child(n, 1), if j >= 0 && is_vector(child(n, 1)) { j } else { -1 });
        out.push(')');
    } else if op == "^" {
        out.append("pow(");
        emit_space(out, child(n, 0), -1);
        out.append(", ");
        emit_space(out, child(n, 1), -1);
        out.push(')');
    } else if k == syntax::N_UNARY && op == "-" {
        out.append("(-");
        emit_space(out, first(n), j);
        out.push(')');
    } else if k == syntax::N_CALL && syntax::child_count(n) == 2 && kind(child(n, 0)) == syntax::N_NAME &&
              position(["sin", "cos", "tan", "exp", "ln", "sqrt", "abs"], name(child(n, 0))) >= 0 {
        format(out, if name(child(n, 0)) == "ln" { "log" } else { name(child(n, 0)) }, "(");
        emit_space(out, arg(n, 0), -1);
        out.push(')');
    } else {
        unsupported(tok(n), str("the expression '", source_of(n), "' is not in this generator"),
                    "write + - * / ^ of numbers, lets, x, y, z, and sin, cos, tan, exp, ln, sqrt or abs");
    }
}

recursive fn mentions_xyz(n: i64) -> bool {
    if kind(n) == syntax::N_NAME && axis_of(name(n)) >= 0 { return true; }
    if kind(n) == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = mentions_xyz(codegen::find_value(name(n)));
        codegen::depth--;
        return r;
    }
    var c = first(n);
    while c >= 0 {
        if mentions_xyz(c) { return true; }
        c = next(c);
    }
    return false;
}

recursive fn is_vector(n: i64) -> bool {
    let k = kind(n);
    if k == syntax::N_VECTOR { return true; }
    if k == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = is_vector(codegen::find_value(name(n)));
        codegen::depth--;
        return r;
    }
    if k == syntax::N_BINARY && (syntax::text(tok(n)) == "*" || syntax::text(tok(n)) == "/") { return is_vector(child(n, 0)) || is_vector(child(n, 1)); }
    if k == syntax::N_BINARY || k == syntax::N_UNARY { return is_vector(first(n)); }
    return false;
}

// The coefficients of the physics at the omega of the frequency that the
// generator writes: ν, and c as re and im.
fn edge_terms() -> f64, f64[2] {
    var nu = 1.0;
    if curl_coef >= 0 { nu = if curl_divides { 1.0 / scalar_value(curl_coef) } else { scalar_value(curl_coef) }; }
    var c: f64[2] = [0.0, 0.0];
    for t, i in mass_terms {
        let z = cval(t);
        c[0] += mass_signs[i] * z[0];
        c[1] += mass_signs[i] * z[1];
    }
    return nu, c;
}

fn emit_edge_bc(out: u8[>..]&, b: Bc) {
    format(out, "    // bc ", surfs[b.surf].name, "\n    {\n    ");
    emit_region(out, b.surf);
    out.append("        let held = elements::");
    if physics == DARCY {
        out.append(if b.kind == BC_LOAD { "pressure(" } else { "normal_flux(" });
        lit(out, scalar_value(b.node));
        out.append(");\n");
    } else if b.kind == BC_LOAD {
        let z = cval(b.node);
        out.append("impedance(");
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.append(");\n");
    } else {
        out.append("tangent() { p => [");
        for j in 3 {
            if b.node >= 0 && j < dim { emit_space(out, b.node, j); } else { out.append("0.0"); }
            out.append(", ");
        }
        out.append("0.0, 0.0, 0.0] };\n");
    }
    out.append("        region(held, ");
    quoted(out, surfs[b.surf].name);
    out.append(");\n    }\n");
}

// The source text of n on one line.
fn one_line(n: i64) -> u8[>..] {
    var s: u8[>..] = [];
    var space = false;
    for c in source_of(n) {
        let blank = c == ' ' || c == '\n' || c == '\t';
        if blank && !space { s.push(' '); }
        if !blank { s.push(c); }
        space = blank;
    }
    return s;
}

fn emit_probe(out: u8[>..]&, e: i64, label: u8[:]) {
    var v: u8[>..] = [];
    var pre: u8[>..] = [];
    emit_cx(v, pre, e);
    format(out, "    found_all = true;\n", pre, "    let v = ", v, ";\n    point(found_all, ");
    quoted(out, label);
    out.append(");\n");
}

fn emit_edge_report(out: u8[>..]&, step_label: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n    ");
    emit_probe(out, r.expr, r.label);
    out.append("        report(");
    quoted(out, step_label);
    out.append(", ");
    quoted(out, r.label);
    out.append(", cquantity(v));\n    }\n");
}

fn emit_edge_expect(out: u8[>..]&, x: EdgeExpect) {
    format(out, "    // expect ", x.name, "\n    {\n    ");
    emit_probe(out, x.expr, x.name);
    out.append("        expect(");
    quoted(out, x.name);
    out.append(", v[0], ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"\", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}

// A step of the edge physics or of the Darcy flow: a function for each
// frequency of a harmonic step, or one for a static step (omega 0).
fn emit_edge_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    let count = if s.procedure == HARMONIC { s.freq_count } else { 1 };
    format(out, "// step ", s.name, ": ", if s.procedure == HARMONIC { "harmonic" } else { "static" }, "\nfn step_", s.name, "() { step_", s.name, "_1(); }\n\n");
    for i in count { emit_edge_case(out, si, i, results_dir); }
}

fn emit_edge_case(out: u8[>..]&, si: i64, i: i64, results_dir: u8[:]) {
    let s = steps[si];
    let harmonic = s.procedure == HARMONIC;
    let f = if harmonic { freqs[s.first_freq + i] } else { 0.0 };
    omega_now = 2.0 * PI * f;
    var label: u8[>..] = [];
    if harmonic { format(label, s.name, " at ", f, " ", freq_unit()); } else { label.append(s.name); }
    format(out, "fn step_", s.name, "_", i + 1, "() {\n");
    if physics == EDGE {
        let nu, c = edge_terms();
        out.append("    elements::edge(");
        lit(out, nu);
        out.append(", ");
        lit(out, c[0]);
        out.append(", ");
        lit(out, c[1]);
    } else {
        out.append("    elements::darcy(");
        lit(out, scalar_value(darcy_k));
        out.append(", ");
        lit(out, if darcy_f >= 0 { scalar_value(darcy_f) } else { 0.0 });
    }
    out.append(", ");
    lit(out, thickness);
    out.append(");\n");
    // The fixed values first, then the gauge, which keeps them.
    for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_FIX { emit_edge_bc(out, b); } }
    for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_LOAD { emit_edge_bc(out, b); } }
    if gauge_tree { out.append("    elements::gauge();\n"); }
    out.append("    solve(");
    quoted(out, label);
    out.append(");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_edge_report(out, label, r); }
    let last = i == (if harmonic { s.freq_count } else { 1 }) - 1;
    if last { for x in probes { if x.step == si { emit_edge_expect(out, x); } } }
    if last && s.output {
        if physics == DARCY { out.append("    output_darcy("); } else { out.append("    output("); }
        quoted(out, str(results_dir, "/", s.name, ".vtu"));
        out.append(", ");
        quoted(out, unknown_name);
        if physics == DARCY {
            out.append(", ");
            quoted(out, flux_name);
        }
        out.append(");\n");
    }
    if !last { format(out, "    step_", s.name, "_", i + 2, "();\n"); }
    out.append("}\n\n");
}

// The dofs of a mesh with `points` points: 2 for each edge of the edge
// physics, and a facet and a cell of the Darcy flow. A mesh has about 7
// edges, 12 facets and 6 cells for each point in 3d, and 3 edges and 2
// cells in 2d.
fn edge_dofs(points: i64) -> i64 {
    if physics == EDGE { return points * (if dim == 3 { 14 } else { 6 }); }
    return points * (if dim == 3 { 18 } else { 5 });
}
