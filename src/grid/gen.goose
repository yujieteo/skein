// The code generator of the grid models: a Skein syntax tree to a Goose
// program that calls src/grid/grid.goose.
//
// It takes the models with a `method` of the kind `fv` or `dg`. The model
// is 2d. It has one physics with one scalar unknown u, and the physics
// block has one equation of the form of src/grid/grid.goose:
//
//     ct ddt(u) + div(a u) - d div(grad(u)) = f(x, y, t)
//
// with the terms in any sequence and on any side, `lap(u)` for
// div(grad(u)), and ct, a and d constant. Each term with no u is a part of
// f, and can use x, y, t and the varying lets. The mesh is a `rect` with
// `cell: quad`, and the conditions are `bc r: u = g(x, y, t)` on a
// `boundary(rect)` or a `boundary(rect, condition)`. The steps are `static`
// and `transient(end:, dt:, scheme: ssprk3)`. An `expect u == g within e`
// compares the L2 norm of u - g with e (spec §15.3). Each other construct
// gives an error `unsupported` at its token, with a fix.
//
// The values are in the unit system of the header, and the reports give
// u in the unit of its type.
import std;
import math;
namespace grid_gen;
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

// Whether this generator takes the model: a method of the kind fv or dg.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_METHOD && (syntax::text(syntax::nodes[c].aux) == "fv" || syntax::text(syntax::nodes[c].aux) == "dg") { return true; }
        c = next(c);
    }
    return false;
}

// --- the model ----------------------------------------------------------------

let FV = 0;
let DG = 1;
let FLUX_NAMES = ["upwind", "linear_upwind", "central", "lax_friedrichs"];
let FLUX_CONSTANTS = ["grid::UPWIND", "grid::LINEAR_UPWIND", "grid::CENTRAL", "grid::LAX_FRIEDRICHS"];
let GRADIENT_NAMES = ["green_gauss", "least_squares"];
let GRADIENT_CONSTANTS = ["grid::GREEN_GAUSS", "grid::LEAST_SQUARES"];
let LIMITER_NAMES = ["none", "minmod", "van_leer", "venkatakrishnan"];
let LIMITER_CONSTANTS = ["grid::NONE", "grid::MINMOD", "grid::VAN_LEER", "grid::VENKATAKRISHNAN"];

// The essential condition `bc region: u = value`.
struct Bc {
    region: i64,     // an index in curve_names
    value: i64,      // the expression
}

struct Step {
    name: u8[..64],
    transient: bool,
    end: f64,
    dt: f64,           // 0 for `dt: stable`
    first_bc: i64,     // the conditions of the step in `active`
    bc_count: i64,
    first_report: i64,
    report_count: i64,
    output: bool,
}

struct Expect {
    node: i64,
    step: i64,         // the step, or -1 for an item
}

var level = 0;   // skein verify convergence: the mesh has 2^level times the cells on each axis

var shape_names: u8[..64][>..];   // the lets of a rect
var shape_nodes: i64[>..];
var curve_names: u8[..64][>..];   // the lets of a boundary
var curve_nodes: i64[>..];
var bcs: Bc[>..];
var global_bcs: i64[>..];
var active: i64[>..];             // the conditions of each step, in `bcs`
var reports: i64[>..];            // the N_REPORTITEM of each report
var steps: Step[>..];
var expects: Expect[>..];
var material_node = -1;
var physics_name: u8[..64];
var unknown: u8[..64];            // the name of u
var unknown_type: u8[..64];
var initial_node = -1;
var source_terms: i64[>..];       // the terms of f
var source_signs: f64[>..];
var ct = 0.0;
var adv: f64[2] = [0.0, 0.0];
var lap = 0.0;                    // the factor of div(grad(u)) in lhs - rhs
var flip = 1.0;                   // -1 if the equation is the negative of the form
var domain = -1;                  // the rect of the mesh, in shape_nodes
var lo: f64[2] = [0.0, 0.0];
var hi: f64[2] = [0.0, 0.0];
var cells: i64[2] = [0, 0];
var method = -1;
var flux = 1;                     // linear_upwind
var gradient = 0;                 // green_gauss
var limiter = 0;                  // none
var order = 2;
var penalty = 2.0;

fn reset() {
    syntax::errors.clear();
    codegen::reset();
    shape_names.clear();
    shape_nodes.clear();
    curve_names.clear();
    curve_nodes.clear();
    bcs.clear();
    global_bcs.clear();
    active.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    material_node = -1;
    physics_name = [];
    unknown = [];
    unknown_type = [];
    initial_node = -1;
    source_terms.clear();
    source_signs.clear();
    ct = 0.0;
    adv = [0.0, 0.0];
    lap = 0.0;
    flip = 1.0;
    domain = -1;
    method = -1;
    flux = 1;
    gradient = 0;
    limiter = 0;
    order = 2;
    penalty = 2.0;
}

// --- values ---------------------------------------------------------------------

// A scalar or a 2D vector in the unit system of the header.
struct Vec {
    n: i64,
    v: f64[2],
}

recursive fn value(n: i64) -> Vec {
    let k = kind(n);
    if k == syntax::N_NAME {
        var v = codegen::find_value(name(n));
        if v < 0 { v = codegen::property(name(n)); }
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let r = value(v);
            codegen::depth--;
            return r;
        }
    }
    if k == syntax::N_VECTOR {
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
    if (binary(n, "*") || binary(n, "/")) && kind(child(n, 0)) != syntax::N_NUMBER {
        let a = value(child(n, 0));
        let b = value(child(n, 1));
        if a.n == 2 && b.n == 1 {
            let f = if binary(n, "*") { b.v[0] } else { 1.0 / b.v[0] };
            return Vec { n: 2, v: [a.v[0] * f, a.v[1] * f] };
        }
        if a.n == 1 && b.n == 2 && binary(n, "*") { return Vec { n: 2, v: [a.v[0] * b.v[0], a.v[0] * b.v[1]] }; }
    }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        let a = value(first(n));
        if a.n == 2 { return Vec { n: 2, v: [-a.v[0], -a.v[1]] }; }
    }
    let r = codegen::eval(n);
    return Vec { n: 1, v: [r.v[0], 0.0] };
}

fn scalar_value(n: i64) -> f64 {
    let v = value(n);
    if v.n != 1 { unsupported(tok(n), "the value is a vector", "write a scalar"); }
    return v.v[0];
}

// Whether expression n depends on x, y or t, through the lets too.
recursive fn varying(n: i64) -> bool {
    if kind(n) == syntax::N_NAME {
        let w = name(n);
        if w == "x" || w == "y" || w == "t" { return true; }
        let v = codegen::find_value(w);
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

// Whether expression n holds the unknown.
recursive fn has_unknown(n: i64) -> bool {
    if kind(n) == syntax::N_NAME { return name(n) == unknown; }
    var c = first(n);
    while c >= 0 {
        if has_unknown(c) { return true; }
        c = next(c);
    }
    return false;
}

// --- the items ------------------------------------------------------------------

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    let space = syntax::text(tok(child(n, 0)));
    if space != "2d" { unsupported(tok(child(n, 0)), "the fv and dg generator has 2d models", "write 2d in the header"); }
    let sys = child(n, 1);
    for s, i in codegen::SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            codegen::system = codegen::SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

// `let` and `param`: a rect, a boundary, or a value.
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
    }
}

fn find_name(names: u8[..64][:], word: u8[:]) -> i64 {
    for s, i in names { if s == word { return i; } }
    return -1;
}

let PHYSICS_FIX = "write physics heat on plate solves u: Temperature { eq ddt(u) + div(a*u) - d*div(grad(u)) = f }, with a, d constant";

fn physics_item(n: i64) {
    set64(physics_name, name(n));
    var eq_count = 0;
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_UNKNOWN {
            if unknown.len > 0 { unsupported(tok(c), "a second unknown: the fv and dg generator has one", PHYSICS_FIX); }
            set64(unknown, name(c));
            if kind(child(c, 0)) != syntax::N_TYPE_NAME { unsupported(tok(c), "the unknown of the fv and dg generator is a scalar", PHYSICS_FIX); }
            set64(unknown_type, name(child(c, 0)));
            initial_node = child(c, 1);
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            material_node = find_material(name(m));
            if next(m) >= 0 { unsupported(tok(next(m)), "a second material: the fv and dg generator has one", "remove it"); }
        } else if k == syntax::N_EQN {
            eq_count++;
            equation(c);
        } else if k != syntax::N_ON {
            unsupported(tok(c), "the physics statement is not in the fv and dg generator", PHYSICS_FIX);
        }
        c = next(c);
    }
    if eq_count != 1 { unsupported(tok(n), "the fv and dg generator needs one eq", PHYSICS_FIX); }
}

var material_items: i64[>..];

fn find_material(word: u8[:]) -> i64 {
    for m in material_items { if name(m) == word { return m; } }
    return -1;
}

// The terms of `eq left = right`, as the terms of left - right.
fn equation(n: i64) {
    terms(child(n, 0), 1.0);
    terms(child(n, 1), -1.0);
    // The form has -d div(grad(u)) with d >= 0, and ct >= 0.
    if lap > 0.0 || lap == 0.0 && ct < 0.0 { flip = -1.0; }
    ct *= flip;
    lap *= flip;
    adv = [adv[0] * flip, adv[1] * flip];
    if ct < 0.0 { unsupported(tok(n), "ddt(u) and div(grad(u)) have the same sign: the equation is not stable forward in time", PHYSICS_FIX); }
}

// The value of a constant factor, or an error.
fn factor(n: i64) -> f64 {
    if varying(n) { unsupported(tok(n), "a factor of the unknown is not constant in the fv and dg generator", PHYSICS_FIX); }
    return scalar_value(n);
}

// Adds the terms of n times s.
recursive fn terms(n: i64, s: f64) {
    if !has_unknown(n) {
        source_terms.push(n);
        source_signs.push(s);
        return;
    }
    if binary(n, "+") || binary(n, "-") {
        terms(child(n, 0), s);
        terms(child(n, 1), if binary(n, "-") { -s } else { s });
    } else if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        terms(first(n), -s);
    } else if binary(n, "*") && !has_unknown(child(n, 0)) {
        terms(child(n, 1), s * factor(child(n, 0)));
    } else if binary(n, "*") && !has_unknown(child(n, 1)) {
        terms(child(n, 0), s * factor(child(n, 1)));
    } else if binary(n, "/") && !has_unknown(child(n, 1)) {
        terms(child(n, 0), s / factor(child(n, 1)));
    } else if is_call(n, "ddt") && source_of(arg(n, 0)) == unknown {
        ct += s;
    } else if is_call(n, "lap") && source_of(arg(n, 0)) == unknown {
        lap += s;
    } else if is_call(n, "div") {
        flux_terms(arg(n, 0), s);
    } else {
        unsupported(tok(n), str("the term '", source_of(n), "' is not in the fv and dg generator"), PHYSICS_FIX);
    }
}

// Adds the terms of div(n) times s: grad(u) and a*u.
recursive fn flux_terms(n: i64, s: f64) {
    if binary(n, "+") || binary(n, "-") {
        flux_terms(child(n, 0), s);
        flux_terms(child(n, 1), if binary(n, "-") { -s } else { s });
    } else if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        flux_terms(first(n), -s);
    } else if is_call(n, "grad") && source_of(arg(n, 0)) == unknown {
        lap += s;
    } else if (binary(n, "*") || binary(n, "/")) && has_unknown(child(n, 0)) != has_unknown(child(n, 1)) {
        let ui = if has_unknown(child(n, 0)) { 0 } else { 1 };
        let other = child(n, 1 - ui);
        if varying(other) { unsupported(tok(other), "a factor of the unknown is not constant in the fv and dg generator", PHYSICS_FIX); }
        let v = value(other);
        if v.n == 2 && source_of(child(n, ui)) == unknown && binary(n, "*") {
            adv = [adv[0] + s * v.v[0], adv[1] + s * v.v[1]];
        } else if v.n == 1 {
            flux_terms(child(n, ui), if binary(n, "*") { s * v.v[0] } else { s / v.v[0] });
        } else {
            unsupported(tok(n), str("the flux '", source_of(n), "' is not in the fv and dg generator"), PHYSICS_FIX);
        }
    } else {
        unsupported(tok(n), str("the flux '", source_of(n), "' is not in the fv and dg generator"), PHYSICS_FIX);
    }
}

// mesh m on r { size: h, cell: quad }, with r a rect.
fn mesh(n: i64) {
    let r = child(n, 0);
    domain = if kind(r) == syntax::N_NAME { find_name(shape_names, name(r)) } else { -1 };
    if domain < 0 {
        unsupported(tok(r), "the fv and dg generator meshes a rect", "write let plate: Surface = rect((0, 0) m, (1, 1) m), and mesh m on plate");
        return;
    }
    let shape = shape_nodes[domain];
    let a = value(arg(shape, 0));
    let b = value(arg(shape, 1));
    lo = [min(a.v[0], b.v[0]), min(a.v[1], b.v[1])];
    hi = [max(a.v[0], b.v[0]), max(a.v[1], b.v[1])];
    var size = 0.0;
    var quad = false;
    var s = first(child(n, 1));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            size = scalar_value(first(s));
        } else if kind(s) == syntax::N_SETTING && name(s) == "cell" && source_of(first(s)) == "quad" {
            quad = true;
        } else {
            unsupported(tok(s), "the mesh setting is not in the fv and dg generator", "write mesh m on plate { size: 0.05 m, cell: quad }");
        }
        s = next(s);
    }
    if !quad { unsupported(tok(n), "the fv and dg generator meshes with cell: quad", "add cell: quad"); }
    if size <= 0.0 { unsupported(tok(n), "the mesh needs a size", "add size: 0.05 m"); }
    for j in 2 { cells[j] = max(1, round((hi[j] - lo[j]) / max(size, 1.0e-300)) as! i64) << level; }
}

fn method_item(n: i64) {
    let word = syntax::text(syntax::nodes[n].aux);
    if name(n) != physics_name || word != "fv" && word != "dg" {
        unsupported(tok(n), "the fv and dg generator has one method, fv or dg, of the physics", "write method heat: fv { mesh: m }");
        return;
    }
    method = if word == "fv" { FV } else { DG };
    if method == DG { flux = 0; }
    var s = first(child(n, 0));
    while s >= 0 {
        let v = first(s);
        let w = if kind(v) == syntax::N_NAME { name(v) } else { "" };
        if name(s) == "flux" {
            flux = position(FLUX_NAMES, w);
            let good = if method == FV { flux >= 0 && flux <= 2 } else { flux == 0 || flux == 3 || w == "roe" };
            if w == "roe" { flux = 0; }   // Roe is upwind for one scalar
            if !good { unsupported(tok(v), "the flux is not one of the method", if method == FV { "write flux: upwind, linear_upwind or central" } else { "write flux: upwind, lax_friedrichs or roe" }); }
            flux = max(flux, 0);
        } else if name(s) == "gradient" && method == FV {
            gradient = max(0, position(GRADIENT_NAMES, w));
            if position(GRADIENT_NAMES, w) < 0 { unsupported(tok(v), "the gradient is not one of fv", "write gradient: green_gauss or least_squares"); }
        } else if name(s) == "limiter" && method == FV {
            limiter = max(0, position(LIMITER_NAMES, w));
            if position(LIMITER_NAMES, w) < 0 { unsupported(tok(v), "the limiter is not one of fv", "write limiter: none, minmod, van_leer or venkatakrishnan"); }
        } else if name(s) == "order" && method == DG {
            order = scalar_value(v) as! i64;
            if order < 1 || order > 3 { unsupported(tok(v), "the dg generator has order 1, 2 and 3", "write order: 2"); }
        } else if name(s) == "penalty" && method == DG {
            penalty = scalar_value(v);
            if penalty <= 0.0 { unsupported(tok(v), "the penalty is more than 0", "write penalty: 2"); }
        } else if name(s) != "mesh" {
            unsupported(tok(s), "the method setting is not in the fv and dg generator",
                        "use mesh, flux, gradient and limiter of fv, or mesh, order, flux and penalty of dg");
        }
        s = next(s);
    }
}

let BC_FIX = "write bc left: u = 0 K, with let left: Curve = boundary(plate, x == 0 m)";

// `bc r: u = g`. It gives the index in bcs, or -1.
fn bc(n: i64) -> i64 {
    let r = child(n, 0);
    let target = child(n, 1);
    let ci = if kind(r) == syntax::N_NAME { find_name(curve_names, name(r)) } else { -1 };
    if ci < 0 {
        unsupported(tok(r), "the region of a bc is a boundary that a let names", BC_FIX);
        return -1;
    }
    if kind(target) != syntax::N_NAME || name(target) != unknown {
        unsupported(tok(target), "the fv and dg generator has the essential condition u = g", BC_FIX);
        return -1;
    }
    let shape = arg(curve_nodes[ci], 0);
    if kind(shape) != syntax::N_NAME || find_name(shape_names, name(shape)) != domain {
        unsupported(tok(shape), "the boundary is not on the rect of the mesh", BC_FIX);
        return -1;
    }
    bcs.push(Bc { region: ci, value: child(n, 2) });
    return bcs.len - 1;
}

let STEP_FIX = "write step run: static { ... } or step run: transient(end: 1 s, dt: stable, scheme: ssprk3) { ... }";

fn procedure(proc: i64, s: Step&) {
    if name(proc) == "static" && first(proc) < 0 { return; }
    if name(proc) != "transient" {
        unsupported(tok(proc), "the fv and dg generator has static and transient steps", STEP_FIX);
        return;
    }
    s.transient = true;
    var scheme = false;
    var a = first(proc);
    while a >= 0 {
        let w = syntax::text(syntax::nodes[a].aux);
        let v = first(a);
        if w == "end" {
            s.end = scalar_value(v);
        } else if w == "dt" {
            s.dt = if kind(v) == syntax::N_NAME && name(v) == "stable" { 0.0 } else { scalar_value(v) };
        } else if w == "scheme" && source_of(v) == "ssprk3" {
            scheme = true;
        } else {
            unsupported(tok(a), "the fv and dg generator has transient(end:, dt:, scheme: ssprk3)", STEP_FIX);
        }
        a = next(a);
    }
    if !scheme { unsupported(tok(proc), "the fv and dg generator has the explicit scheme ssprk3", "add scheme: ssprk3"); }
    if s.end <= 0.0 || s.dt < 0.0 { unsupported(tok(proc), "the step needs an end after the start, and dt >= 0", STEP_FIX); }
}

fn step(n: i64) {
    var s = Step { name: [], transient: false, end: 0.0, dt: 0.0, first_bc: 0, bc_count: 0, first_report: reports.len,
                   report_count: 0, output: false };
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
            s.output = true;
            if kind(first(c)) != syntax::N_NAME || name(first(c)) != unknown || next(first(c)) >= 0 {
                unsupported(tok(c), "the fv and dg generator writes the unknown", str("write output ", unknown));
            }
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                report(r);
                r = next(r);
            }
        } else if k == syntax::N_EXPECT {
            expect(c, steps.len);
        } else {
            unsupported(tok(c), "the step statement is not in the fv and dg generator", "use only bc, output, report and expect");
        }
        c = next(c);
    }
    // The conditions of the step replace the global ones on the same region.
    s.first_bc = active.len;
    for g in global_bcs {
        var replaced = false;
        for b in own { if bcs[b].region == bcs[g].region { replaced = true; } }
        if !replaced { active.push(g); }
    }
    active.append(own);
    s.bc_count = active.len - s.first_bc;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

let REPORT_FIX = "write report u at (0.5, 0.5) m";

// `u at p`.
fn point_ok(e: i64) -> bool {
    guard kind(e) == syntax::N_AT && kind(child(e, 0)) == syntax::N_NAME && name(child(e, 0)) == unknown else { return false; }
    return value(child(e, 1)).n == 2;
}

fn report(item: i64) {
    if !point_ok(first(item)) {
        unsupported(tok(first(item)), str("the report '", source_of(first(item)), "' is not in the fv and dg generator"), REPORT_FIX);
        return;
    }
    reports.push(item);
}

// `expect e: u at p == v within w`, or `expect e: u == g within w`.
fn expect(n: i64, step: i64) {
    let left = child(n, 0);
    let field = kind(left) == syntax::N_NAME && name(left) == unknown;
    if !field && !point_ok(left) {
        unsupported(tok(left), "the fv and dg generator compares u at a point, or u with a field (L2)", "write expect e: u == sin(pi*x/L)*1 K within 1%");
        return;
    }
    expects.push(Expect { node: n, step: step });
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// it.
fn read(prog: i64) -> bool {
    reset();
    material_items.clear();
    var physics_node = -1;
    var mesh_node = -1;
    var method_node = -1;
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
        } else if k == syntax::N_MATERIAL {
            material_items.push(c);
        } else if k == syntax::N_PHYSICS && physics_node < 0 {
            physics_node = c;
        } else if k == syntax::N_MESH && mesh_node < 0 {
            mesh_node = c;
        } else if k == syntax::N_METHOD && method_node < 0 {
            method_node = c;
        } else if k != syntax::N_SCALES && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in the fv and dg generator"),
                        "remove it: this generator runs one physics on one mesh with one method");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_node < 0 || mesh_node < 0 || method_node < 0 {
        unsupported(at, "the model needs a physics, a mesh and a method", PHYSICS_FIX);
        return false;
    }
    physics_item(physics_node);
    codegen::scope_material = material_node;
    if ok() { method_item(method_node); }
    if ok() { mesh(mesh_node); }
    if !ok() { return false; }
    c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_BC {
            let b = bc(c);
            if b >= 0 { global_bcs.push(b); }
        } else if kind(c) == syntax::N_STEP {
            step(c);
        } else if kind(c) == syntax::N_EXPECT {
            expect(c, -1);
        }
        c = next(c);
    }
    if steps.len == 0 { unsupported(at, "the model has no step", STEP_FIX); }
    for s in steps { if s.transient && ct <= 0.0 { unsupported(at, "a transient step needs ddt(u) in the equation", PHYSICS_FIX); } }
    return ok();
}

// The dofs of the mesh.
fn dof_estimate() -> i64 { return cells[0] * cells[1] * (if method == DG { (order + 1) * (order + 1) } else { 1 }); }

// The order of the L2 error of the method for a smooth solution
// (docs/checks.md §4.3), and its reason.
fn expected_order() -> f64, u8[>..] {
    if method == DG { return order + 1.0, str("dg, order ", order); }
    if (adv[0] != 0.0 || adv[1] != 0.0) && flux == 0 { return 1.0, str("fv, flux upwind"); }
    return 2.0, str("fv, flux ", FLUX_NAMES[flux]);
}

// --- the Goose program ----------------------------------------------------------

let FUNCTIONS = ["sin", "cos", "tan", "exp", "log", "sqrt", "abs", "tanh"];
let COND_FIX = "write a condition of x and y with == < <= > >= && || and !";

// The Goose expression of the field n of x, y and t.
recursive fn emit(out: u8[>..]&, n: i64) {
    let k = kind(n);
    if !varying(n) {
        out.push('(');
        lit(out, scalar_value(n));
        out.push(')');
    } else if k == syntax::N_NAME && (name(n) == "x" || name(n) == "y" || name(n) == "t") {
        out.append(name(n));
    } else if k == syntax::N_NAME && codegen::depth < 100 {
        codegen::depth++;
        emit(out, codegen::find_value(name(n)));
        codegen::depth--;
    } else if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(n))) >= 0 {
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
    } else if k == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && position(FUNCTIONS, name(child(n, 0))) >= 0 && syntax::child_count(n) == 2 {
        format(out, name(child(n, 0)), "(");
        emit(out, arg(n, 0));
        out.push(')');
    } else {
        unsupported(tok(n), str("the expression '", source_of(n), "' is not in the fv and dg generator"),
                    "write + - * / ^ of numbers, lets, x, y, t, and sin, cos, tan, exp, log, sqrt, abs or tanh");
    }
}

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
        unsupported(tok(n), str("the condition '", source_of(n), "' is not in the fv and dg generator"), COND_FIX);
        return;
    }
    // Each comparison has the tolerance of spec §8.2.
    out.append(if i <= 1 { "abs(" } else { "(" });
    emit(out, child(n, 0));
    out.append(" - ");
    emit(out, child(n, 1));
    out.append([") <= tol", ") > tol", ") < -tol", ") <= tol", ") > tol", ") >= -tol"][i]);
}

let PRELUDE = """
    // The reports and the verification records, in TOON (command_line.md §6).
    var reports: u8[>..] = [];
    var report_count = 0;
    var checks: u8[>..] = [];
    var check_count = 0;
    var solves: u8[>..] = [];
    var solve_count = 0;
    var results: u8[>..] = [];
    var result_count = 0;
    var errors: u8[>..] = [];
    var error_count = 0;
    var failed = false;
    var tol = 0.0;

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
        if !pass { failed = true; }
    }

    // A bc on a region that holds no face is a fault of the model.
    fn region(found: bool, name: u8[:]) {
        if !found {
            write_stderr(str("skein run: no boundary face of the mesh is on the region '", name, "'\n"));
            exit(1);
        }
    }

    // u at (x, y). A point outside the mesh is a fault of the model.
    fn at(x: f64, y: f64, name: u8[:]) -> f64 {
        let v, found = grid::value_at(x, y);
        if !found {
            write_stderr(str("skein run: the point of '", name, "' is not in the mesh\n"));
            failed = true;
        }
        return v;
    }

    fn solved(step: u8[:], converged: bool, rtol: f64) {
        format(solves, "  ", step, ",", grid::cell_count(), ",", grid::dof_count(), ",", grid::iterations, ",",
               grid::linear_iterations, ",", grid::last_seconds, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if rtol > 0.0 {
            let r = grid::last_residual / max(grid::first_residual, 1.0e-300);
            check("solve", "residual", converged, str(r), str("<= ", rtol));
        } else {
            check("solve", "finite", converged, str(grid::iterations, " steps"), "each value finite");
        }
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", abs(got - want) <= t, str(name, ": ", quantity(got, unit)), limit);
    }

    // The L2 error of u against a field, its record for skein verify
    // convergence, and its check.
    fn l2(name: u8[:], error: f64, norm: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        format(errors, "  ", name, ",", grid::cell_count(), ",", error, ",", norm, "\n");
        error_count++;
        let t = if relative { within * norm } else { within };
        check("post", "expect", error <= t, str(name, ": L2 error ", quantity(error, unit)), limit);
    }

    fn output(path: u8[:], name: u8[:]) {
        if !grid::write_vtu(path, name) {
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
        format(s, "solve[", solve_count, "]{step,cells,dofs,iterations,gmres_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        format(s, "l2[", error_count, "]{expect,cells,error,norm}:\n", errors);
        s.append(solvers::memory_records());
        write_stdout(s);
    }
""";

// The unit of u in the reports.
fn unit() -> const u8[:] {
    if unknown_type == "Temperature" { return "K"; }
    if unknown_type == "Length" { return codegen::system.length; }
    if unknown_type == "Pressure" { return codegen::system.pressure; }
    return "";
}

// The functions of the fields: f, the values of the conditions, the
// initial value, and the field of each L2 expect.
fn emit_fields(out: u8[>..]&) {
    out.append("fn source(x: f64, y: f64, t: f64) -> f64 {\n    return 0.0");
    for n, i in source_terms {
        out.append(" + ");
        lit(out, -flip * source_signs[i]);
        out.append(" * ");
        emit(out, n);
    }
    out.append(";\n}\n\nfn boundary(k: i64, x: f64, y: f64, t: f64) -> f64 {\n");
    for b, k in bcs {
        format(out, "    if k == ", k, " { return ");
        emit(out, b.value);
        out.append("; }\n");
    }
    out.append("    return 0.0;\n}\n\nfn initial(x: f64, y: f64) -> f64 {\n    let t = 0.0;\n    return ");
    if initial_node >= 0 { emit(out, initial_node); } else { out.append("0.0"); }
    out.append(";\n}\n\n");
    for x in expects {
        if kind(child(x.node, 0)) != syntax::N_NAME { continue; }
        format(out, "fn exact_", name(x.node), "(x: f64, y: f64, t: f64) -> f64 {\n    return ");
        emit(out, child(x.node, 1));
        out.append(";\n}\n\n");
    }
    for c, i in curve_names {
        format(out, "fn select_", c, "() -> bool {\n    return grid::select_faces() { x, y => ");
        if syntax::child_count(curve_nodes[i]) <= 2 { out.append("true"); } else { emit_cond(out, arg(curve_nodes[i], 1)); }
        out.append(" };\n}\n\n");
    }
}

fn emit_point(out: u8[>..]&, e: i64) {
    let p = value(child(e, 1));
    out.append("at(");
    lit(out, p.v[0]);
    out.append(", ");
    lit(out, p.v[1]);
    out.append(", name)");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    let within = child(n, 2);
    let relative = kind(within) == syntax::N_PERCENT;
    let tolerance = scalar_value(within);
    format(out, "    // expect ", name(n), "\n    {\n        let name = ");
    quoted(out, name(n));
    out.append(";\n");
    if kind(child(n, 0)) == syntax::N_NAME {
        format(out, "        let error, norm = grid::l2_error(exact_", name(n), ");\n        l2(name, error, norm, ");
    } else {
        out.append("        expect(name, ");
        emit_point(out, child(n, 0));
        out.append(", ");
        lit(out, scalar_value(child(n, 1)));
        out.append(", ");
    }
    lit(out, tolerance);
    format(out, ", ", if relative { "true" } else { "false" }, ", \"", unit(), "\", ");
    quoted(out, str(source_of(child(n, 1)), " within ", source_of(within)));
    out.append(");\n    }\n");
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, "\nfn step_", s.name, "() {\n    grid::clear_bcs();\n");
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        let c = curve_names[bcs[b].region];
        format(out, "    region(select_", c, "(), \"", c, "\");\n    grid::fix_selected(", b, ");\n");
    }
    if s.transient {
        out.append("    solved(");
        quoted(out, s.name);
        out.append(", grid::solve_transient(");
        lit(out, s.end);
        out.append(", ");
        lit(out, s.dt);
        out.append(", source, boundary), 0.0);\n");
    } else {
        out.append("    solved(");
        quoted(out, s.name);
        out.append(", grid::solve_static(1.0e-10, source, boundary), 1.0e-10);\n");
    }
    for r in reports[s.first_report..s.first_report + s.report_count] {
        format(out, "    {\n        let name = ");
        quoted(out, source_of(r));
        out.append(";\n        report(");
        quoted(out, s.name);
        out.append(", name, quantity(");
        emit_point(out, first(r));
        format(out, ", \"", unit(), "\"));\n    }\n");
    }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        out.append("    output(");
        quoted(out, str(results_dir, "/", s.name, ".vtu"));
        out.append(", ");
        quoted(out, unknown);
        out.append(");\n");
    }
    out.append("}\n\n");
}

// Writes the Goose program of the model that `read` read. The program
// imports the grid package from `grid_path`, and writes its VTU files in
// `results_dir`.
fn generate(grid_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    out.append("// The values are in the unit system of the model header.\nimport std;\nimport os;\nimport math;\nimport grid from ");
    quoted(out, grid_path);
    // The grid package imports the package solvers, beside it.
    out.append(";\nimport solvers from ");
    quoted(out, str(grid_path[..grid_path.len - 4], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    out.push('\n');
    emit_fields(out);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    format(out, "fn main() {\n    // mesh: ", cells[0], " x ", cells[1], " rectangles\n    grid::mesh_rect(");
    lit(out, lo[0]);
    out.append(", ");
    lit(out, lo[1]);
    out.append(", ");
    lit(out, hi[0]);
    out.append(", ");
    lit(out, hi[1]);
    format(out, ", ", cells[0], ", ", cells[1], ");\n    tol = grid::tolerance();\n");
    format(out, "    // physics ", physics_name, ": ct ddt(u) + div(a u) - d div(grad(u)) = f\n    grid::equation(");
    lit(out, ct);
    out.append(", ");
    lit(out, adv[0]);
    out.append(", ");
    lit(out, adv[1]);
    out.append(", ");
    lit(out, -lap);
    out.append(");\n");
    if method == FV {
        format(out, "    grid::fv(", FLUX_CONSTANTS[flux], ", ", GRADIENT_CONSTANTS[gradient], ", ", LIMITER_CONSTANTS[limiter], ");\n");
    } else {
        format(out, "    grid::dg(", order, ", ", FLUX_CONSTANTS[flux], ", ");
        lit(out, penalty);
        out.append(");\n");
    }
    out.append("    grid::start(initial);\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
