// The part of the simplex generator for the edge physics of
// `family: nedelec` and the Darcy flow of `family: raviart_thomas` (spec
// §13), on P1 cells (src/elements/edge.goose and src/elements/mixed.goose).
//
// The edge physics: `let h = curl(u)/X` (or `X*curl(u)`, or `curl(u)`) and
// `eq curl(h) + c*u = 0` for a vector u, Complex in `harmonic` steps. Its
// bcs: `cross(n, u) = cross(n, g)` or `u = 0` (essential), and
// `cross(n, h) = z*u` (natural). `gauge: tree` holds a spanning tree of the
// cells with c = 0. `order: 2` gives the elements of degree 2 on the P1 cells
// (src/elements/edge.goose). A `modal(count:, shift:)` step gives the
// resonances of the physics (src/elements/edge_modal.goose), with a z of the
// bc that depends on omega.
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
var mass_ddt2: i64[>..];            // 1 for a term c*ddt2(u): -omega^2*c in a harmonic step
var skew_form = false;              // let h = -2*skew(grad(u))/X, eq div(h) + ...: h is curl(u)/X, and h*n is cross(n, h)
var port_bc: i64[>..];              // the token of the bc of each source term of a natural condition
var port_node: i64[>..];            // its vector field
var port_const: i64[>..];           // its constant factor, or -1
var port_sign: f64[>..];
var port_zbc: i64[>..];             // the token of a bc with a term z*(u - dot(u, n)*n), and the node of its z

var gauge_tree = false;
var edge_rhs = -1;                  // the source s of eq curl(h) + c*u = s, or -1 for 0
var edge_order = 1;                 // the order of the elements: 1, or 2 with `order: 2` (the mesh stays P1)
var port_zsign: f64[>..];           // the sign of the term z*u of each bc of port_zbc
var port_zf_bc: i64[>..];           // the factors of z: the token of the bc, the node and 1 for a divisor
var port_zf_node: i64[>..];
var port_zf_div: i64[>..];
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
    mass_ddt2.clear();
    skew_form = false;
    port_bc.clear();
    port_node.clear();
    port_const.clear();
    port_sign.clear();
    port_zbc.clear();
    gauge_tree = false;
    edge_rhs = -1;
    edge_order = 1;
    port_zsign.clear();
    port_zf_bc.clear();
    port_zf_node.clear();
    port_zf_div.clear();
    darcy_k = -1;
    darcy_f = -1;
    physics_lets.clear();
    probes.clear();
    region_names.clear();
    region_shapes.clear();
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

// The text of s with no space, newline or parenthesis.
fn compact(s: u8[:]) -> u8[>..] {
    var r: u8[>..] = [];
    for ch in s { if ch != ' ' && ch != '\n' && ch != '\t' && ch != '(' && ch != ')' { r.push(ch); } }
    return r;
}

// `let h = curl(u)/X`, `X*curl(u)` or `curl(u)`, or the form of a stress
// tensor, `-2*skew(grad(u))/X`, whose h*n is cross(n, curl(u)/X). Another
// let is a value of the reports.
fn curl_let(c: i64) {
    let v = child(c, 1);
    let cu = str("curl(", unknown_name, ")");
    guard curl_name.len == 0 else { return; };
    let skew = str("-2*skew(grad(", unknown_name, "))");
    if compact(source_of(v)) == compact(skew) {
        skew_form = true;
        set64(curl_name, name(c));
    } else if kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "/" && compact(source_of(child(v, 0))) == compact(skew) {
        skew_form = true;
        curl_coef = child(v, 1);
        curl_divides = true;
        set64(curl_name, name(c));
    } else if source_of(v) == cu {
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
    if source_of(child(c, 1)) != "0" { edge_rhs = child(c, 1); }
    while kind(left) == syntax::N_BINARY && (syntax::text(tok(left)) == "+" || syntax::text(tok(left)) == "-") {
        let t = child(left, 1);
        let ddt2 = kind(t) == syntax::N_BINARY && syntax::text(tok(t)) == "*" && source_of(child(t, 1)) == str("ddt2(", unknown_name, ")");
        guard kind(t) == syntax::N_BINARY && syntax::text(tok(t)) == "*" && (source_of(child(t, 1)) == unknown_name || ddt2) else {
            unsupported(tok(t), str("the term '", source_of(t), "' is not c*", unknown_name, " or c*ddt2(", unknown_name, ")"), fix);
            return;
        };
        mass_terms.push(child(t, 0));
        mass_ddt2.push(if ddt2 { 1 } else { 0 });
        mass_signs.push(if syntax::text(tok(left)) == "-" { -1.0 } else { 1.0 });
        left = child(left, 0);
    }
    if curl_name.len == 0 || source_of(left) != str(if skew_form { "div(" } else { "curl(" }, curl_name, ")") {
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
    if physics == EDGE && skew_form && t == str(curl_name, "*n") {
        // Hs*n is cross(n, h): z*(u - dot(u, n)*n) is the impedance, the other terms a field.
        b.kind = BC_LOAD;
        port_terms(b, v, 1.0);
        bcs.push(b);
        return bcs.len - 1;
    } else if physics == EDGE && kind(target) == syntax::N_MEMBER && kind(first(target)) == syntax::N_NAME &&
              name(first(target)) == unknown_name && axis_of(name(target)) >= 0 && source_of(v) == "0" {
        // u.y = 0 with u.z = 0 on a plane of x holds the tangent of u. emit_edge_bc checks the pair.
        b.comp = axis_of(name(target));
        b.node = -1;
    } else if physics == EDGE && tangent && source_of(v) == "0" {
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

// The terms of the value n of Hs*n = n, with sign: z*(u - dot(u, n)*n) gives
// the impedance z of the bc, and a vector field times a constant gives a
// source.
recursive fn port_terms(b: Bc, n: i64, sign: f64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if k == syntax::N_BINARY && (op == "+" || op == "-") {
        port_terms(b, child(n, 0), sign);
        port_terms(b, child(n, 1), if op == "-" { -sign } else { sign });
        return;
    }
    if k == syntax::N_UNARY && op == "-" {
        port_terms(b, first(n), -sign);
        return;
    }
    let tangential = compact(str(unknown_name, "-dot(", unknown_name, ", n)*n"));
    // z*(u - dot(u, n)*n), z*u, or a product or a quotient with cross(n, cross(n, u)) = -(u - dot(u, n)*n):
    // a product of factors that name no unknown and no x, y, z, and the tangent of u.
    if k == syntax::N_BINARY && (op == "*" || op == "/") && position(port_zbc, b.tok) < 0 {
        var nodes: i64[>..] = [];
        var divs: i64[>..] = [];
        flatten_product(n, 0, nodes, divs);
        var at = -1;
        var turn = 1.0;
        var clean = true;
        for m, i in nodes {
            let c = compact(source_of(m));
            let t = if divs[i] != 0 { 0.0 } else if c == tangential || c == compact(unknown_name) { 1.0 } else if c == compact(str("cross(n, cross(n, ", unknown_name, "))")) { -1.0 } else { 0.0 };
            if t != 0.0 && at < 0 {
                at = i;
                turn = t;
            } else if has_unknown(m) || mentions_xyz(m) {
                clean = false;
            }
        }
        if at >= 0 && clean && sign > 0.0 {
            port_zbc.push(b.tok);
            port_zsign.push(turn);
            for m, i in nodes {
                if i == at { continue; }
                port_zf_bc.push(b.tok);
                port_zf_node.push(m);
                port_zf_div.push(divs[i]);
            }
            return;
        }
    }
    var constant = -1;
    var field = -1;
    if k == syntax::N_BINARY && op == "*" {
        if has_unknown(child(n, 1)) || mentions_xyz(child(n, 1)) { field = child(n, 1); constant = child(n, 0); }
        else if has_unknown(child(n, 0)) || mentions_xyz(child(n, 0)) { field = child(n, 0); constant = child(n, 1); }
    }
    if field >= 0 && constant >= 0 && !has_unknown(field) && !has_unknown(constant) && !mentions_xyz(constant) {
        port_bc.push(b.tok);
        port_node.push(field);
        port_const.push(constant);
        port_sign.push(sign);
    } else {
        unsupported(tok(n), str("the term '", source_of(n), "' of Hs*n is not in this generator"),
                    str("write c*(", unknown_name, " - dot(", unknown_name, ", n)*n) or c*cross(n, cross(n, ", unknown_name, "))/d for the impedance and c*g for a source, with g a vector of x, y and z"));
    }
}

// The factors of a product and a quotient of nodes, with 1 in divs for a divisor.
recursive fn flatten_product(n: i64, div: i64, nodes: i64[>..]&, divs: i64[>..]&) {
    if kind(n) == syntax::N_BINARY && (syntax::text(tok(n)) == "*" || syntax::text(tok(n)) == "/") {
        flatten_product(child(n, 0), div, nodes, divs);
        flatten_product(child(n, 1), if syntax::text(tok(n)) == "/" { 1 - div } else { div }, nodes, divs);
        return;
    }
    nodes.push(n);
    divs.push(div);
}

// The z of the bc with token t at omega_now, for the harmonic step and the
// start value of the modal step: the product of its factors, with its sign.
fn zcoef(t: i64) -> f64[2] {
    var z: f64[2] = [1.0, 0.0];
    let i = position(port_zbc, t);
    if i < 0 { return [0.0, 0.0]; }
    for b, j in port_zf_bc {
        if b != t { continue; }
        z = if port_zf_div[j] != 0 { cdiv(z, cval(port_zf_node[j])) } else { cmul(z, cval(port_zf_node[j])) };
    }
    return [port_zsign[i] * z[0], port_zsign[i] * z[1]];
}

// Whether n has omega, in itself or in a let.
recursive fn mentions_omega(n: i64) -> bool {
    if kind(n) == syntax::N_NAME && name(n) == "omega" { return true; }
    if kind(n) == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = mentions_omega(codegen::find_value(name(n)));
        codegen::depth--;
        return r;
    }
    var c = first(n);
    while c >= 0 {
        if mentions_omega(c) { return true; }
        c = next(c);
    }
    return false;
}

// The Goose expression of the Complex value of n as f64[2], with omega the
// variable w.
recursive fn emit_cz(out: u8[>..]&, n: i64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if !mentions_omega(n) {
        let z = cval(n);
        out.push('[');
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.push(']');
    } else if k == syntax::N_NAME && name(n) == "omega" {
        out.append("[w, 0.0]");
    } else if k == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        emit_cz(out, codegen::find_value(name(n)));
        codegen::depth--;
    } else if k == syntax::N_BINARY && (op == "+" || op == "-" || op == "*" || op == "/") {
        out.append(if op == "+" { "cadd(" } else if op == "-" { "csub(" } else if op == "*" { "cmul(" } else { "cdiv(" });
        emit_cz(out, child(n, 0));
        out.append(", ");
        emit_cz(out, child(n, 1));
        out.push(')');
    } else if k == syntax::N_UNARY && op == "-" {
        out.append("csub([0.0, 0.0], ");
        emit_cz(out, first(n));
        out.push(')');
    } else if is_call(n, "sqrt") && syntax::child_count(n) == 2 {
        out.append("csqrt(");
        emit_cz(out, arg(n, 0));
        out.push(')');
    } else if op == "^" && (scalar_value(child(n, 1)) == 0.5 || scalar_value(child(n, 1)) == floor(scalar_value(child(n, 1))) && scalar_value(child(n, 1)) >= 1.0 && scalar_value(child(n, 1)) <= 8.0) {
        let e = scalar_value(child(n, 1));
        if e == 0.5 {
            out.append("csqrt(");
            emit_cz(out, child(n, 0));
            out.push(')');
        } else {
            out.append("cpow(");
            emit_cz(out, child(n, 0));
            format(out, ", ", e as! i64, ")");
        }
    } else {
        unsupported(tok(n), str("the expression '", source_of(n), "' of omega is not in this generator"),
                    "write + - * / of numbers, lets, im and omega, sqrt, ^ 2 or ^ 0.5");
    }
}

// The z of the bc b at the variable w (omega) of a modal step.
fn emit_zcoef(out: u8[>..]&, b: Bc) {
    let i = position(port_zbc, b.tok);
    if i < 0 && b.node >= 0 {
        // cross(n, h) = z*u
        emit_cz(out, b.node);
        return;
    }
    var acc: u8[>..] = [];
    acc.append("[1.0, 0.0]");
    for t, j in port_zf_bc {
        if t != b.tok { continue; }
        var next_acc: u8[>..] = [];
        next_acc.append(if port_zf_div[j] != 0 { "cdiv(" } else { "cmul(" });
        next_acc.append(acc);
        next_acc.append(", ");
        emit_cz(next_acc, port_zf_node[j]);
        next_acc.push(')');
        acc.clear();
        acc.append(next_acc);
    }
    out.append("cmul([");
    lit(out, if i >= 0 { port_zsign[i] } else { 0.0 });
    format(out, ", 0.0], ", acc, ")");
}

// --- reports and expects ---------------------------------------------------------

// Whether n holds a field: an `at`, an integral or far_field.
recursive fn has_field(n: i64) -> bool {
    if kind(n) == syntax::N_AT || is_call(n, "integral") || is_call(n, "far_field") || is_call(n, "rms") || is_call(n, "eigenvalue") { return true; }
    if complex_solid && is_call(n, "mean") { return true; }
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
    if is_call(n, "rms") && syntax::child_count(n) == 3 && complex_solid {
        format(out, "[rms_value(", rms_index(n), "), 0.0]");
    } else if !has_field(n) {
        let z = cval(n);
        out.push('[');
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.push(']');
    } else if is_call(n, "eigenvalue") && syntax::child_count(n) == 2 && physics == EDGE {
        format(out, "modal_value(", scalar_value(arg(n, 0)) as! i64, ")");
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
    } else if k == syntax::N_AT && wave_eddy {
        let p = value(child(n, 1));
        if p.n != dim || source_of(child(n, 0)) != unknown_name {
            probe_fail(tok(n), "the point is not a vector, or the field is not the unknown", "write u at (x, y) mm");
            return;
        }
        out.append("wave_point(");
        lit(out, p.v[0]);
        out.append(", ");
        lit(out, p.v[1]);
        out.append(", ");
        lit(out, p.v[2]);
        out.push(')');
    } else if k == syntax::N_AT {
        let p = value(child(n, 1));
        if p.n != dim { probe_fail(tok(n), "the point is not a vector", "write at (x, y, z) mm"); }
        emit_field(out, child(n, 0), p.v, p.v);
    } else if is_call(n, "mean") && syntax::child_count(n) == 3 && complex_solid {
        emit_mean(out, pre, n);
    } else if is_call(n, "integral") && syntax::child_count(n) == 3 {
        emit_integral(out, pre, n);
    } else if is_call(n, "far_field") && syntax::child_count(n) == 5 {
        emit_far_field(out, pre, n);
    } else {
        probe_fail(tok(n), str("the expression '", source_of(n), "' is not in this generator"), PROBE_FIX);
    }
}

// mean(u.x, s) of the Complex unknown on a boundary s that a let names.
fn emit_mean(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let f = arg(n, 0);
    let r = arg(n, 1);
    let s = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
    guard kind(f) == syntax::N_MEMBER && source_of(first(f)) == unknown_name && axis_of(name(f)) >= 0 && s >= 0 else {
        probe_fail(tok(n), "mean takes a component of the unknown and a boundary that a let names", "write mean(u.z, tip)");
        return;
    };
    emit_region(pre, s);
    format(pre, "    let mean_", n, " = solid_mean(", axis_of(name(f)), ");\n");
    format(out, "mean_", n);
}

// The rms calls of the expression n, in the order of the program.
var rms_nodes: i64[>..];

// The number of the rms call n.
fn rms_index(n: i64) -> i64 {
    for m, i in rms_nodes { if m == n { return i; } }
    if rms_nodes.len >= 16 { probe_fail(tok(n), "more than 16 rms calls", "use at most 16 rms calls in a model"); return 0; }
    rms_nodes.push(n);
    return rms_nodes.len - 1;
}

// The first rms call in the expression n, or -1.
recursive fn find_rms(n: i64) -> i64 {
    if is_call(n, "rms") && syntax::child_count(n) == 3 { return n; }
    var c = first(n);
    while c >= 0 {
        let m = find_rms(c);
        if m >= 0 { return m; }
        c = next(c);
    }
    return -1;
}

// integral(F.x, c) on a line c, or integral(dot(q, n), s) on a boundary s.
fn emit_integral(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let f = arg(n, 0);
    let r = arg(n, 1);
    if wave_eddy {
        emit_cell_integral(out, pre, n);
        return;
    }
    let s = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
    if physics == DARCY && s >= 0 && source_of(f) == str("dot(", flux_name, ", n)") {
        emit_region(pre, s);
        format(pre, "    let flux_", s, " = elements::boundary_flux();\n");
        format(out, "[flux_", s, ", 0.0]");
        return;
    }
    if physics == EDGE && s >= 0 && is_call(f, "dot") && syntax::child_count(f) == 3 {
        // The integral over a boundary of dot(u, g) or dot(g, g), with g a vector of x, y and z.
        let a = arg(f, 0);
        let b = arg(f, 1);
        let field = source_of(a) == unknown_name || source_of(b) == unknown_name;
        let g = if source_of(a) == unknown_name { b } else { a };
        if has_unknown(g) || (!field && source_of(a) != source_of(b)) {
            probe_fail(tok(f), "the integral on a boundary is dot(u, g) or dot(g, g) with g a vector of x, y and z", PROBE_FIX);
            return;
        }
        emit_region(pre, s);
        format(pre, "    let facet_", n, " = elements::facet_integral(", if field { "true" } else { "false" }, ") { p => [");
        if field {
            for j in 3 {
                if j < dim { emit_space(pre, g, j); } else { pre.append("0.0"); }
                pre.append(", ");
            }
            pre.append("0.0, 0.0, 0.0");
        } else {
            for j in dim {
                if j > 0 { pre.append(" + "); }
                emit_space(pre, g, j);
                pre.append(" * ");
                emit_space(pre, g, j);
            }
            pre.append(", 0.0, 0.0, 0.0, 0.0, 0.0");
        }
        pre.append("] };\n");
        format(out, "facet_", n);
        return;
    }
    let c = if kind(r) == syntax::N_NAME { codegen::find_value(name(r)) } else { -1 };
    guard c >= 0 && is_call(c, "line") && syntax::child_count(c) == 3 else {
        probe_fail(tok(r), "the integral is on a line that a let names", PROBE_FIX);
        return;
    };
    emit_field(out, f, value(arg(c, 0)).v, value(arg(c, 1)).v);
}

var region_names: u8[..64][>..];   // the regions of the integrals, and their shapes: a shape is made once
var region_shapes: i64[>..];

fn cached_region_shape(r: i64) -> i64 {
    for w, i in region_names { if w == name(r) { return region_shapes[i]; } }
    let g = region_shape(r);
    var w: u8[..64] = [];
    set64(w, name(r));
    region_names.push(w);
    region_shapes.push(g);
    return g;
}

// The material item of the region r: the one of the law `on r`, or the
// base material.
fn region_item(r: i64) -> i64 {
    for w in region_laws { if ends_with(w.label, str(" on ", source_of(r))) { return w.item; } }
    return material_node;
}

// integral(f, r) over the cells of a region r, with f linear in u and ddt(u) and
// in the properties of the material of r: f = a u + b gives a I + b A, with I
// the integral of u and A the area of r.
fn emit_cell_integral(out: u8[>..]&, pre: u8[>..]&, n: i64) {
    let f = arg(n, 0);
    let r = arg(n, 1);
    guard kind(r) == syntax::N_NAME && find_geometry(name(r)) >= 0 else {
        probe_fail(tok(r), "the integral is on a region that a let names", PROBE_FIX);
        return;
    };
    let base = material_node;
    material_node = region_item(r);
    eddy_on = true;
    eddy_u = [0.0, 0.0];
    let b = cval(f);
    eddy_u = [1.0, 0.0];
    let a1 = cval(f);
    eddy_u = [2.0, 0.0];
    let a2 = cval(f);
    eddy_on = false;
    material_node = base;
    if abs(a2[0] - 2.0 * a1[0] + b[0]) > 1.0e-9 * (abs(a2[0]) + abs(b[0]) + 1.0e-300) || abs(a2[1] - 2.0 * a1[1] + b[1]) > 1.0e-9 * (abs(a2[1]) + abs(b[1]) + 1.0e-300) {
        probe_fail(tok(f), "the integrand is not affine in the unknown", PROBE_FIX);
    }
    format(pre, "    let cells_", n, " = wave_cells(", cached_region_shape(r), ");\n");
    out.append("cadd(cmul(");
    format(out, "[", a1[0] - b[0], ", ", a1[1] - b[1], "]");
    format(out, ", [cells_", n, "[0], cells_", n, "[1]]), [");
    lit(out, b[0]);
    format(out, " * cells_", n, "[2], ");
    lit(out, b[1]);
    format(out, " * cells_", n, "[2]])");
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
    if k == syntax::N_NAME && complex_solid && name(f) == unknown_name {
        if a[0] != b[0] || a[1] != b[1] || a[2] != b[2] { probe_fail(tok(f), "the integral of u on a line is not in this generator", PROBE_FIX); }
        emit_probe_sample(out, "solid_", "", comp, a, b);
    } else if k == syntax::N_NAME && physics == EDGE && name(f) == unknown_name {
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
    var r = Report { expr: first(item), label: [], row: -1 };
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

// The part of the program of a Complex vector unknown of the elastic physics.
let HARMONIC_PRELUDE = """
    // Component comp of u at a point, as [re, im].
    fn solid_point(comp: i64, x: f64, y: f64, z: f64) -> f64[2] {
        let re, found_re = elements::harmonic_at(0, comp, x, y, z);
        let im, found_im = elements::harmonic_at(1, comp, x, y, z);
        if !found_re || !found_im { found_all = false; }
        return [re, im];
    }

    // Its mean on the region of the last region call.
    fn solid_mean(comp: i64) -> f64[2] {
        return [elements::harmonic_mean(0, comp), elements::harmonic_mean(1, comp)];
    }

    // The root mean square of a random response (spec §14.5): the integral of
    // |v|^2 psd on the frequencies, by the trapezoidal rule, one slot for
    // each rms call.
    var rms_sum: f64[16] = [0.0; 16];
    var rms_y: f64[16] = [0.0; 16];
    var rms_f: f64[16] = [0.0; 16];
    var rms_n: i64[16] = [0; 16];

    fn rms_add(k: i64, v: f64[2], psd: f64, f: f64) {
        let y = (v[0] * v[0] + v[1] * v[1]) * psd;
        if rms_n[k] > 0 { rms_sum[k] += 0.5 * (rms_y[k] + y) * (f - rms_f[k]); }
        rms_y[k] = y;
        rms_f[k] = f;
        rms_n[k]++;
    }

    fn rms_value(k: i64) -> f64 { return sqrt(rms_sum[k]); }

    """;

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

    fn csqrt(a: f64[2]) -> f64[2] {
        let r = sqrt(a[0] * a[0] + a[1] * a[1]);
        let im = sqrt(max((r - a[0]) / 2.0, 0.0));
        return [sqrt(max((r + a[0]) / 2.0, 0.0)), if a[1] < 0.0 { -im } else { im }];
    }

    fn cpow(a: f64[2], n: i64) -> f64[2] {
        var r: f64[2] = [1.0, 0.0];
        for i in n { r = cmul(r, a); }
        return r;
    }

    // The record of a modal step, which has no iterations of GMRES.
    fn resonance_done(step: u8[:]) {
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",0,", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // eigenvalue(k) of the last modal step: the frequency omega / (2 pi), as [re, im]. NaN if the step did not converge.
    fn modal_value(k: i64) -> f64[2] {
        if !converged { return [0.0 / 0.0, 0.0 / 0.0]; }
        let re, im = elements::modal_frequency(k);
        return [re, im];
    }

    fn output_resonance(path: u8[:], field: u8[:], k: i64) {
        if !converged { return; }
        if !elements::write_modal_mode(path, field, k) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
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
    } else if k == syntax::N_NAME && property(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        emit_space(out, property(name(n)), j);
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
    } else if (is_call(n, "min") || is_call(n, "max") || is_call(n, "atan2")) && syntax::child_count(n) == 3 {
        format(out, name(child(n, 0)), "(");
        emit_space(out, arg(n, 0), -1);
        out.append(", ");
        emit_space(out, arg(n, 1), -1);
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
                    "write + - * / ^ of numbers, lets, x, y, z, and sin, cos, tan, exp, ln, sqrt, abs, min, max or atan2");
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
    if kind(n) == syntax::N_NAME && property(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = mentions_xyz(property(name(n)));
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
        // ddt2 of a harmonic step is (im*omega)^2.
        let k = if mass_ddt2[i] != 0 { -omega_now * omega_now } else { 1.0 };
        c[0] += mass_signs[i] * k * z[0];
        c[1] += mass_signs[i] * k * z[1];
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
    } else if b.kind == BC_LOAD && physics == EDGE && skew_form {
        // Hs*n = z*(u - dot(u, n)*n) + g: the impedance, then the source.
        let zi = position(port_zbc, b.tok);
        let z = if zi >= 0 { zcoef(b.tok) } else { [0.0, 0.0] };
        out.append("impedance(");
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.append(");\n");
        for t, i in port_bc {
            if t != b.tok { continue; }
            let c = cval(port_const[i]);
            out.append("        let loaded = elements::traction() { p => [");
            for part in 2 {
                for j in 3 {
                    if j < dim {
                        lit(out, port_sign[i] * c[part]);
                        out.append(" * ");
                        emit_space(out, port_node[i], j);
                    } else {
                        out.append("0.0");
                    }
                    out.append(", ");
                }
            }
            out.append("] };\n        region(loaded, ");
            quoted(out, surfs[b.surf].name);
            out.append(");\n");
        }
    } else if b.kind == BC_LOAD {
        let z = cval(b.node);
        out.append("impedance(");
        lit(out, z[0]);
        out.append(", ");
        lit(out, z[1]);
        out.append(");\n");
    } else if b.comp >= 0 {
        out.append("tangent() { p => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] };\n");
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

// Whether bc k is the first of the bcs u.a = 0 of its boundary. The boundary
// must be planes of one axis, and the bcs must hold the two other components,
// so they hold the tangent of u.
fn tangent_pair(list: Bc[:], k: i64) -> bool {
    var mask = 0;
    for b, i in list {
        if b.comp < 0 || b.surf != list[k].surf { continue; }
        if i < k { return false; }
        mask |= 1 << b.comp;
    }
    let sf = surfs[list[k].surf];
    var plane = sf.axes.len > 0 && sf.shapes.len == 0 && !sf.hole;
    for a in sf.axes { if a != sf.axes[0] { plane = false; } }
    let want = if plane { 7 - (1 << sf.axes[0]) } else { -1 };
    if dim != 3 || !plane || mask != want {
        unsupported(list[k].tok, "the components of u that a bc sets to 0 are not the tangent of u on planes of one axis",
                    str("write bc ", sf.name, ": u.y = 0 and bc ", sf.name, ": u.z = 0 on a plane of x, or cross(n, u) = 0"));
        return false;
    }
    return true;
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

// The source s of eq curl(h) + c*u = s: for the base material and for each
// material on a region, the cells of that material get the integral of s
// with the properties of that material (the vector js of a coil, 0 in the
// air). s is a real vector function of x, y and z.
fn emit_edge_sources(out: u8[>..]&) {
    guard edge_rhs >= 0 else { return; }
    let base = material_node;
    emit_edge_source_in(out, -1);
    for w in region_laws {
        material_node = w.item;
        emit_edge_source_in(out, w.g);
    }
    material_node = base;
}

fn emit_edge_source_in(out: u8[>..]&, g: i64) {
    if !mentions_xyz(edge_rhs) {
        let v = value(edge_rhs);
        if v.v[0] == 0.0 && v.v[1] == 0.0 && v.v[2] == 0.0 { return; }
    }
    format(out, "    elements::current_in(", g, ") { p => [");
    for j in 3 {
        if j < dim { emit_space(out, edge_rhs, j); } else { out.append("0.0"); }
        out.append(", ");
    }
    out.append("0.0, 0.0, 0.0] };\n");
}

// A step of the edge physics or of the Darcy flow: a function for each
// frequency of a harmonic step, or one for a static step (omega 0).
fn emit_edge_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    if s.procedure == MODAL {
        emit_edge_modal(out, si, results_dir);
        return;
    }
    let count = if s.procedure == HARMONIC { s.freq_count } else { 1 };
    format(out, "// step ", s.name, ": ", if s.procedure == HARMONIC { "harmonic" } else { "static" }, "\nfn step_", s.name, "() { step_", s.name, "_1(); }\n\n");
    for i in count { emit_edge_case(out, si, i, results_dir); }
}

// The terms of the rms calls of step si at the frequency f: the value of
// their function and the density of the PSD, for the trapezoidal rule of
// rms_add.
fn emit_rms_adds(out: u8[>..]&, si: i64, f: f64) {
    var calls: i64[>..] = [];
    for r in reports[steps[si].first_report..steps[si].first_report + steps[si].report_count] {
        let m = find_rms(r.expr);
        if m >= 0 { calls.push(m); }
    }
    // An expect of no step reads the last harmonic step.
    var last = -1;
    for t, k in steps { if t.procedure == HARMONIC { last = k; } }
    for x in probes { if x.step == si || x.step < 0 && si == last { let m = find_rms(x.expr); if m >= 0 { calls.push(m); } } }
    for m in calls {
        var v: u8[>..] = [];
        var pre: u8[>..] = [];
        emit_cx(v, pre, arg(m, 0));
        let psd = cval(arg(m, 1));
        format(out, "    // rms: the term of this frequency\n    {\n        found_all = true;\n    ", pre, "        let v = ", v, ";\n        rms_add(", rms_index(m), ", v, ");
        lit(out, psd[0]);
        out.append(", ");
        lit(out, f);
        out.append(");\n        point(found_all, \"rms\");\n    }\n");
    }
}

fn emit_edge_case(out: u8[>..]&, si: i64, i: i64, results_dir: u8[:]) {
    let s = steps[si];
    let harmonic = s.procedure == HARMONIC;
    let f = if harmonic { freqs[s.first_freq + i] } else { 0.0 };
    omega_now = 2.0 * PI * f;
    var label: u8[>..] = [];
    if harmonic { format(label, s.name, " at ", f, " ", freq_unit()); } else { label.append(s.name); }
    format(out, "fn step_", s.name, "_", i + 1, "() {\n");
    if complex_solid {
        out.append("    elements::clear_bcs();\n    elements::set_harmonic(");
        lit(out, omega_now);
        out.append(", ");
        lit(out, s.damping);
        out.append(");\n");
    } else if physics == EDGE {
        let nu, c = edge_terms();
        if edge_order == 2 { out.append("    elements::edge_order(2);\n"); }
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
    if !complex_solid {
        out.append(", ");
        lit(out, thickness);
        out.append(");\n");
        if physics == EDGE {
            let base = material_node;
            for w in region_laws {
                material_node = w.item;
                let nu, c = edge_terms();
                format(out, "    cells(elements::edge_in(", w.g, ", ");
                lit(out, nu);
                out.append(", ");
                lit(out, c[0]);
                out.append(", ");
                lit(out, c[1]);
                out.append("), ");
                quoted(out, w.label);
                out.append(");\n");
            }
            material_node = base;
            emit_edge_sources(out);
        }
    }
    // The fixed values first, then the gauge, which keeps them.
    if complex_solid {
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if !is_load(b) { emit_bc(out, b); } }
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if is_load(b) { emit_bc(out, b); } }
    } else {
        let step_bcs = active[s.first_bc..s.first_bc + s.bc_count];
        for b, k in step_bcs {
            if b.kind != BC_FIX { continue; }
            if b.comp >= 0 && !tangent_pair(step_bcs, k) { continue; }
            emit_edge_bc(out, b);
        }
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_LOAD { emit_edge_bc(out, b); } }
    }
    if gauge_tree { out.append("    elements::gauge();\n"); }
    out.append("    solve(");
    quoted(out, label);
    out.append(");\n");
    if complex_solid { emit_rms_adds(out, si, f); }
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

// modal(count:, shift:) of the edge physics: the `count` resonances nearest
// the shift (src/elements/edge_modal.goose). The shift is a frequency, kept
// in freqs.
fn edge_modal_proc(proc: i64, s: Step&) {
    let fix = "write step resonance: modal(count: 1, shift: 8 GHz) { ... }, with a term c*ddt2(u) in the eq";
    s.procedure = MODAL;
    s.count = 0;
    s.first_freq = freqs.len;
    var a = first(proc);
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        let w = if named >= 0 { syntax::text(named) } else { "" };
        if w == "count" {
            s.count = scalar_value(first(a)) as! i64;
        } else if w == "shift" {
            freqs.push(scalar_value(first(a)));
        } else {
            unsupported(tok(a), "the argument of modal is not in this generator", fix);
        }
        a = next(a);
    }
    s.freq_count = freqs.len - s.first_freq;
    var has_mass = false;
    for d in mass_ddt2 { if d != 0 { has_mass = true; } }
    if s.count < 1 { unsupported(tok(proc), "a modal step needs count >= 1", fix); }
    if s.freq_count != 1 || !(freqs[s.first_freq] > 0.0) {
        unsupported(tok(proc), "a modal step of an edge physics needs a shift > 0, above the null space of the curl", fix);
    }
    if !has_mass { unsupported(tok(proc), "a modal step needs a term c*ddt2(u) in the eq", fix); }
}

// The coefficient m (re and im) of the terms c*ddt2(u) at the material of
// material_node: the matrix M of the eigenvalue problem.
fn edge_mass_terms() -> f64[2] {
    var m: f64[2] = [0.0, 0.0];
    for t, i in mass_terms {
        if mass_ddt2[i] == 0 { continue; }
        let z = cval(t);
        m[0] += mass_signs[i] * z[0];
        m[1] += mass_signs[i] * z[1];
    }
    return m;
}

// A modal step: the matrices at omega 0 (c without the term of ddt2), the
// mass m of ddt2 on the cells, the conditions, then the resonances.
fn emit_edge_modal(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    let step_bcs = active[s.first_bc..s.first_bc + s.bc_count];
    // The z of each condition cross(n, h) = z*u as a function of omega.
    format(out, "// step ", s.name, ": modal\nfn wall_z_", s.name, "(g: i64, w: f64) -> f64[2] {\n");
    var group = 0;
    for b in step_bcs {
        if b.kind != BC_LOAD { continue; }
        format(out, "    if g == ", group, " { return ");
        emit_zcoef(out, b);
        out.append("; }\n");
        group++;
    }
    out.append("    return [0.0, 0.0];\n}\n\n");
    omega_now = 0.0;
    format(out, "fn step_", s.name, "() {\n");
    if edge_order == 2 { out.append("    elements::edge_order(2);\n"); }
    if edge_rhs >= 0 { unsupported(tok(edge_rhs), "a modal step has no source", "write eq curl(h) + c*ddt2(u) = 0"); }
    let nu, c = edge_terms();
    let m = edge_mass_terms();
    out.append("    elements::edge(");
    lit(out, nu);
    out.append(", ");
    lit(out, c[0]);
    out.append(", ");
    lit(out, c[1]);
    out.append(", ");
    lit(out, thickness);
    out.append(");\n    elements::edge_mass_in(-1, ");
    lit(out, m[0]);
    out.append(", ");
    lit(out, m[1]);
    out.append(");\n");
    let base = material_node;
    for w in region_laws {
        material_node = w.item;
        let wnu, wc = edge_terms();
        let wm = edge_mass_terms();
        format(out, "    cells(elements::edge_in(", w.g, ", ");
        lit(out, wnu);
        out.append(", ");
        lit(out, wc[0]);
        out.append(", ");
        lit(out, wc[1]);
        out.append("), ");
        quoted(out, w.label);
        format(out, ");\n    cells(elements::edge_mass_in(", w.g, ", ");
        lit(out, wm[0]);
        out.append(", ");
        lit(out, wm[1]);
        out.append("), ");
        quoted(out, w.label);
        out.append(");\n");
    }
    material_node = base;
    for b, k in step_bcs {
        if b.kind != BC_FIX { continue; }
        if b.comp >= 0 && !tangent_pair(step_bcs, k) { continue; }
        emit_edge_bc(out, b);
    }
    group = 0;
    for b in step_bcs {
        if b.kind != BC_LOAD { continue; }
        for t in port_bc {
            if t == b.tok { unsupported(b.tok, "a modal step has no source in a condition", "write the condition with the terms of z*(u - dot(u, n)*n) only"); }
        }
        format(out, "    // bc ", surfs[b.surf].name, "\n    {\n    ");
        emit_region(out, b.surf);
        format(out, "        let held = elements::impedance_group(", group, ");\n        region(held, ");
        quoted(out, surfs[b.surf].name);
        out.append(");\n    }\n");
        group++;
    }
    var label: u8[>..] = [];
    label.append(s.name);
    format(out, "    converged = elements::resonance(", s.count, ", ");
    lit(out, 2.0 * PI * freqs[s.first_freq]);
    format(out, ") { g, w => wall_z_", s.name, "(g, w) };\n    resonance_done(");
    quoted(out, label);
    out.append(");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_edge_report(out, label, r); }
    for x in probes { if x.step == si { emit_edge_expect(out, x); } }
    if s.output {
        for k in 1..s.count + 1 {
            out.append("    output_resonance(");
            quoted(out, str(results_dir, "/", s.name, "_mode", k, ".vtu"));
            out.append(", ");
            quoted(out, unknown_name);
            format(out, ", ", k, ");\n");
        }
    }
    out.append("}\n\n");
}

// The dofs of a mesh with `points` points: 2 for each edge of the edge
// physics, and a facet and a cell of the Darcy flow. A mesh has about 7
// edges, 12 facets and 6 cells for each point in 3d, and 3 edges and 2
// cells in 2d.
fn edge_dofs(points: i64) -> i64 {
    if physics == EDGE && edge_order == 2 { return points * (if dim == 3 { 78 } else { 20 }); }
    if physics == EDGE { return points * (if dim == 3 { 14 } else { 6 }); }
    return points * (if dim == 3 { 18 } else { 5 });
}
