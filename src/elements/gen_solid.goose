import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import .gen;
import .gen_contact;
import .gen_aniso;

let SOLID_J2 = 1;
let SOLID_J2_FINITE = 2;
let SOLID_ENERGY = 3;

let PLASTIC_FIX = "write plastic { yield: von_mises(sigma) - (sy + Hk*alpha), strain: eps_p, hardening: alpha }, with state eps_p: Tensor2<Real> = 0 and state alpha: Real = 0";
let ENERGY_FIX = "write let S: Tensor2<Pressure> = 2*diff(W, C), with C = transpose(defgrad(u))*defgrad(u) and W of principal(C), det(C) and tr(C)";

var solid_law = 0;
var finite_strain = false;
var finite_eq = false;
var pk_name: u8[..64];
var plastic_tok = -1;
var state_names: u8[..64][>..];
var strain_state: u8[..64];
var yield_node = -1;
var hard_node = -1;
var solid_lets: i64[>..];
var c_name: u8[..64];
var principal_name: u8[..64];
var energy_node = -1;
var energy_fns: i64[>..];
var velocity_name: u8[..64];
var velocity_start = -1;
var dynamic_eq = false;

fn solid_reset() {
    solid_law = 0;
    finite_strain = false;
    finite_eq = false;
    pk_name = [];
    plastic_tok = -1;
    state_names.clear();
    strain_state = [];
    yield_node = -1;
    hard_node = -1;
    solid_lets.clear();
    c_name = [];
    principal_name = [];
    energy_node = -1;
    energy_fns.clear();
    velocity_name = [];
    velocity_start = -1;
    dynamic_eq = false;
    aniso_reset();
}

let DYNAMIC_FIX = "write solves u: Vector<Length>, v: Vector<Velocity> = (0, 0, -V0), with eq ddt(u) - v = 0 and eq rho*ddt(v) - div(defgrad(u)*S) = 0";

fn solid_unknown(c: i64) {
    if velocity_name.len > 0 { unsupported(tok(c), "a third unknown", DYNAMIC_FIX); }
    set64(velocity_name, name(c));
    velocity_start = child(c, 1);
}

fn solid_dynamic_eq(c: i64) -> bool {
    guard source_of(child(c, 1)) == "0" else { return false; }
    let left = child(c, 0);
    if kind(left) == syntax::N_BINARY && syntax::text(tok(left)) == "-" && source_of(child(left, 0)) == str("rho*ddt2(", unknown_name, ")") &&
       is_call(child(left, 1), "div") && source_of(arg(child(left, 1), 0)) == flux_name {
        dynamic_eq = true;
        return true;
    }
    guard velocity_name.len > 0 else { return false; }
    guard kind(left) == syntax::N_BINARY && syntax::text(tok(left)) == "-" else { return false; }
    let a = child(left, 0);
    let b = child(left, 1);
    if source_of(a) == str("ddt(", unknown_name, ")") && source_of(b) == velocity_name { return true; }
    let inertia = source_of(a) == str("rho*ddt(", velocity_name, ")");
    guard inertia && is_call(b, "div") && solid_eq(arg(b, 0)) else { return false; }
    dynamic_eq = true;
    return true;
}

fn traction_target() -> u8[>..] {
    if pk_name.len > 0 { return str(pk_name, "*n"); }
    if finite_eq { return str("defgrad(", unknown_name, ")*", flux_name, "*n"); }
    return str(flux_name, "*n");
}

fn solid_velocity_bc(target: i64, v: i64, b: Bc&) -> bool {
    guard velocity_name.len > 0 && kind(target) == syntax::N_MEMBER && source_of(first(target)) == velocity_name else { return false; }
    b.comp = axis_of(name(target));
    if b.comp < 0 || source_of(v) != "0" { unsupported(tok(v), "this generator holds a component of v at 0", "write bc wall: v.z = 0"); }
    return true;
}

fn solid_state(n: i64) {
    if source_of(child(n, 1)) != "0" { unsupported(tok(n), "a state starts at 0 in this generator", PLASTIC_FIX); }
    var w: u8[..64] = [];
    set64(w, name(n));
    state_names.push(w);
}

fn is_state(word: u8[:]) -> bool {
    for s in state_names { if s == word { return true; } }
    return false;
}

fn solid_let(n: i64) -> bool {
    let v = child(n, 1);
    let u = unknown_name;
    let s = if is_call(v, "hooke") && syntax::child_count(v) == 4 { arg(v, 2) } else { -1 };
    if s >= 0 && source_of(arg(v, 0)) == "E" && source_of(arg(v, 1)) == "nu" && kind(s) == syntax::N_BINARY &&
       syntax::text(tok(s)) == "-" && source_of(child(s, 0)) == str("strain(", u, ")") && is_state(source_of(child(s, 1))) {
        set64(strain_state, source_of(child(s, 1)));
        set_flux(n);
        return true;
    }
    if flux_name.len > 0 && source_of(v) == str("defgrad(", u, ")*", flux_name) && physics == ELASTIC {
        set64(pk_name, name(n));
        return true;
    }
    if source_of(v) == str("transpose(defgrad(", u, "))*defgrad(", u, ")") {
        set64(c_name, name(n));
        return true;
    }
    if c_name.len > 0 && source_of(v) == str("principal(", c_name, ")") {
        set64(principal_name, name(n));
        return true;
    }
    let d = if kind(v) == syntax::N_BINARY && source_of(child(v, 0)) == "2" { child(v, 1) } else { -1 };
    if d >= 0 && syntax::text(tok(v)) == "*" && is_call(d, "diff") && syntax::child_count(d) == 3 && c_name.len > 0 &&
       source_of(arg(d, 1)) == c_name {
        for l in solid_lets { if name(l) == source_of(arg(d, 0)) { energy_node = child(l, 1); } }
        if energy_node < 0 { unsupported(tok(d), "W of diff(W, C) is not a let of the physics before it", ENERGY_FIX); }
        solid_law = SOLID_ENERGY;
        set_flux(n);
        return true;
    }
    solid_lets.push(n);
    return true;
}

fn set_flux(n: i64) {
    if physics >= 0 && physics != ELASTIC { unsupported(tok(n), "a second flux", PHYSICS_FIX); }
    physics = ELASTIC;
    set64(flux_name, name(n));
}

fn solid_eq(x: i64) -> bool {
    if pk_name.len > 0 && source_of(x) == pk_name {
        finite_eq = true;
        return true;
    }
    guard source_of(x) == str("defgrad(", unknown_name, ")*", flux_name) else { return false; }
    finite_eq = true;
    return true;
}

fn solid_plastic(n: i64) {
    plastic_tok = tok(n);
    var s = first(first(n));
    while s >= 0 {
        let v = first(s);
        if name(s) == "yield" {
            plastic_yield(v);
        } else if name(s) == "strain" && source_of(v) != strain_state {
            unsupported(tok(v), "the strain of plastic is not the state in hooke", PLASTIC_FIX);
        } else if name(s) == "hardening" && !is_state(source_of(v)) {
            unsupported(tok(v), "the hardening of plastic is not a state", PLASTIC_FIX);
        } else if name(s) != "strain" && name(s) != "hardening" {
            unsupported(tok(s), "the setting of plastic is not in this generator", PLASTIC_FIX);
        }
        s = next(s);
    }
    solid_law = SOLID_J2;
}

fn plastic_yield(v: i64) {
    let ok = kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "-" && source_of(child(v, 0)) == str("von_mises(", flux_name, ")");
    let r = if ok { child(v, 1) } else { -1 };
    let sum = r >= 0 && kind(r) == syntax::N_BINARY && syntax::text(tok(r)) == "+";
    let h = if sum { child(r, 1) } else { -1 };
    yield_node = if sum { child(r, 0) } else { r };
    if h >= 0 && kind(h) == syntax::N_BINARY && syntax::text(tok(h)) == "*" && is_state(source_of(child(h, 1))) {
        hard_node = child(h, 0);
    } else if sum || !ok {
        unsupported(tok(v), "this generator has the yield von_mises(sigma) - (sy + Hk*alpha)", PLASTIC_FIX);
    }
}

fn solid_finish(n: i64) {
    if pk_name.len > 0 && solid_law == 0 && strain_state.len == 0 {
        if !finite_strain || !finite_eq {
            unsupported(tok(n), "P = defgrad(u)*S needs assume finite_strain and eq div(P) = 0", "write both");
        }
        finite_strain = false;
        finite_eq = false;
    }
    if solid_law == SOLID_ENERGY || strain_state.len > 0 {
        if finite_strain != finite_eq {
            unsupported(tok(n), "assume finite_strain goes with eq div(defgrad(u)*S) = 0", "write both, or neither");
        }
        if dim == 2 && plane != PLANE_STRAIN { unsupported(tok(n), "a 2d nonlinear solid needs plane_strain", "add assume plane_strain(1 mm)"); }
    }
    if strain_state.len > 0 && plastic_tok < 0 { unsupported(tok(n), "hooke(E, nu, strain(u) - eps_p) needs a plastic block", PLASTIC_FIX); }
    if solid_law == SOLID_J2 && finite_strain { solid_law = SOLID_J2_FINITE; }
    if velocity_name.len > 0 && (!dynamic_eq || solid_law != SOLID_J2_FINITE) {
        unsupported(tok(n), "the first-order form is in this generator for J2 at finite strain only", DYNAMIC_FIX);
    }
    if solid_law == SOLID_ENERGY && !finite_strain { unsupported(tok(n), "an energy W needs assume finite_strain", ENERGY_FIX); }
    if finite_strain && solid_law == 0 { unsupported(tok(n), "finite_strain needs a plastic block or an energy", ENERGY_FIX); }
    if solid_law == SOLID_ENERGY { energy_calls(energy_node); }
    for l in solid_lets {
        if solid_law != SOLID_ENERGY && aniso_node < 0 {
            unsupported(tok(l), str("the let '", name(l), "' is not a flux of this generator"),
                        "write let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u)), or let q: Vector<HeatFlux> = -k*grad(T)");
        }
    }
}

recursive fn energy_calls(n: i64) {
    guard n >= 0 else { return; }
    if kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME {
        let f = find_fn(name(child(n, 0)));
        var known = false;
        for g in energy_fns { if g == f { known = true; } }
        if f >= 0 && !known {
            energy_fns.push(f);
            energy_calls(syntax::child(f, syntax::child_count(f) - 1));
        }
    }
    if kind(n) == syntax::N_NAME {
        for l in solid_lets { if name(l) == name(n) { energy_calls(child(l, 1)); } }
    }
    var c = first(n);
    while c >= 0 {
        energy_calls(c);
        c = next(c);
    }
}

fn find_fn(word: u8[:]) -> i64 {
    var c = first(program_node);
    while c >= 0 {
        if kind(c) == syntax::N_FN && name(c) == word { return c; }
        c = next(c);
    }
    return -1;
}

fn solid_property(word: u8[:]) -> bool {
    if solid_law == SOLID_ENERGY { return true; }
    return solid_law > 0 && (word == source_of(yield_node) || hard_node >= 0 && word == source_of(hard_node));
}

fn solid_ramp(proc: i64, s: Step&) -> bool {
    let a = first(proc);
    guard (solid_law > 0 || has_ties()) && a >= 0 && next(a) < 0 && syntax::nodes[a].aux >= 0 && syntax::text(syntax::nodes[a].aux) == "ramp" else { return false; }
    s.count = scalar_value(first(a)) as! i64;
    if s.count < 1 { unsupported(tok(a), "the ramp needs at least 1 increment", "write static(ramp: 10)"); }
    return true;
}

let SOLID_PRELUDE = """
    // A static step of the nonlinear solid. Its cg_iterations are those of
    // all its Newton iterations.
    fn nonlinear(step: u8[:], done: bool) {
        converged = done;
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // An explicit step of the nonlinear solid from the start velocity
    // (vx, vy, vz). Its cg_iterations are its time steps.
    fn explicit_solid(step: u8[:], end: f64, dt: f64, vx: f64, vy: f64, vz: f64) {
        converged = elements::explicit_solid(end, dt, vx, vy, vz);
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    """;

fn emit_solid_law(out: u8[>..]&) {
    if solid_law == SOLID_ENERGY {
        out.append("    elements::hyperelastic(");
        lit(out, thickness);
        out.append(");\n");
        return;
    }
    guard solid_law > 0 else { return; }
    out.append("    elements::plastic(");
    lit(out, scalar_value(yield_node));
    out.append(", ");
    lit(out, if hard_node >= 0 { scalar_value(hard_node) } else { 0.0 });
    format(out, ", ", if solid_law == SOLID_J2_FINITE { "true" } else { "false" }, ");\n");
}

fn emit_solid_solve(out: u8[>..]&, s: Step) {
    if s.procedure == TRANSIENT {
        let v0 = if velocity_start >= 0 { value(velocity_start).v } else { [0.0; 3] };
        out.append("    explicit_solid(");
        quoted(out, s.name);
        for x in [s.end, s.dt, v0[0], v0[1], v0[2]] {
            out.append(", ");
            lit(out, x);
        }
        out.append(");\n");
        return;
    }
    out.append("    nonlinear(");
    quoted(out, s.name);
    if solid_law == SOLID_ENERGY {
        format(out, ", elements::solve_energy(", s.count, ", energy));\n");
    } else {
        format(out, ", elements::solve_plastic(", s.count, "));\n");
    }
}

fn emit_energy(out: u8[>..]&) {
    out.append(SOLID_PRELUDE);
    guard solid_law == SOLID_ENERGY else { return; }
    for f in energy_fns {
        format(out, "// fn ", name(f), " of the model, in hyper-dual numbers\nfn energy_", name(f), "(");
        let count = syntax::child_count(f);
        for i in count - 2 {
            if i > 0 { out.append(", "); }
            format(out, name(child(f, i)), ": elements::Hd");
        }
        out.append(") -> elements::Hd {\n    return ");
        hd_expr(out, child(f, count - 1), f);
        out.append(";\n}\n\n");
    }
    let c = if principal_name.len > 0 { principal_name[..] } else { "c" };
    format(out, "// W of physics ", physics_name, ": c holds the principal values of C\nfn energy(", c, ": elements::Hd[3]) -> elements::Hd {\n");
    for l in solid_lets {
        if child(l, 1) == energy_node { continue; }
        format(out, "    let ", name(l), " = ");
        hd_expr(out, child(l, 1), -1);
        out.append(";\n");
    }
    out.append("    return ");
    hd_expr(out, energy_node, -1);
    out.append(";\n}\n\n");
}

fn is_param(f: i64, word: u8[:]) -> bool {
    guard f >= 0 else { return false; }
    for i in syntax::child_count(f) - 2 { if name(child(f, i)) == word { return true; } }
    return false;
}

recursive fn hd_expr(out: u8[>..]&, n: i64, f: i64) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY { syntax::text(tok(n)) } else { "" };
    if op == "+" || op == "-" {
        out.append("(");
        hd_expr(out, child(n, 0), f);
        format(out, " ", op, " ");
        hd_expr(out, child(n, 1), f);
        out.append(")");
    } else if op == "*" || op == "/" || op == "^" {
        out.append(if op == "*" { "elements::hd_mul(" } else if op == "/" { "elements::hd_div(" } else { "elements::hd_pow_hd(" });
        hd_expr(out, child(n, 0), f);
        out.append(", ");
        hd_expr(out, child(n, 1), f);
        out.append(")");
    } else if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        out.append("(");
        hd_expr(out, first(n), f);
        out.append(" * -1.0)");
    } else if k == syntax::N_CALL {
        hd_call(out, n, f);
    } else if k == syntax::N_MEMBER && f < 0 && principal_name.len > 0 && source_of(first(n)) == principal_name && axis_of(name(n)) >= 0 {
        format(out, principal_name, "[", axis_of(name(n)), "]");
    } else if k == syntax::N_NAME && (is_param(f, name(n)) || f < 0 && energy_let(name(n))) {
        out.append(name(n));
    } else if k == syntax::N_NUMBER || k == syntax::N_NAME {
        out.append("elements::hd(");
        lit(out, scalar_value(n));
        out.append(")");
    } else {
        unsupported(tok(n), str("'", source_of(n), "' is not in the energy of this generator"), ENERGY_FIX);
    }
}

fn energy_let(word: u8[:]) -> bool {
    for l in solid_lets { if name(l) == word && child(l, 1) != energy_node { return true; } }
    return false;
}

let HD_CALLS = ["sqrt", "ln", "log", "exp"];
let HD_FNS = ["elements::hd_sqrt(", "elements::hd_ln(", "elements::hd_ln(", "elements::hd_exp("];

recursive fn hd_call(out: u8[>..]&, n: i64, f: i64) {
    let callee = source_of(child(n, 0));
    let count = syntax::child_count(n) - 1;
    let c = if principal_name.len > 0 { principal_name[..] } else { "c" };
    if count == 1 && c_name.len > 0 && source_of(arg(n, 0)) == c_name && (callee == "det" || callee == "tr") {
        if callee == "det" {
            format(out, "elements::hd_mul(elements::hd_mul(", c, "[0], ", c, "[1]), ", c, "[2])");
        } else {
            format(out, "(", c, "[0] + ", c, "[1] + ", c, "[2])");
        }
        return;
    }
    for h, i in HD_CALLS {
        if callee == h && count == 1 {
            out.append(HD_FNS[i]);
            hd_expr(out, arg(n, 0), f);
            out.append(")");
            return;
        }
    }
    let g = find_fn(callee);
    guard g >= 0 && syntax::child_count(g) - 2 == count else {
        unsupported(tok(n), str("the call '", callee, "' is not in the energy of this generator"),
                    "use det(C), tr(C), sqrt, ln, exp, and fn items of the model");
        return;
    };
    format(out, "energy_", callee, "(");
    for i in count {
        if i > 0 { out.append(", "); }
        hd_expr(out, arg(n, i), f);
    }
    out.append(")");
}
