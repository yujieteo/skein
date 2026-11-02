// The code generator of the models of conservation laws: a Skein syntax
// tree to a Goose program that calls src/cons/cons.goose.
//
// It takes the models with a method fv on a physics of more than one
// unknown, or of a vector unknown, whose equations are the conservation
// form `eq ddt(U) + div(F) = S`. U is a scalar, a vector or a tensor of 2
// components on each axis. F, S and the conditions are written with the
// lets of the physics: + - * / ^ of scalars, vectors (x, y) and tensors,
// dot, outer, norm, transpose, tr, sym, ddot, I, n, x, y, t, and the
// functions sqrt, exp, ln, sin, cos, tan, atan, abs, min and max. The
// generator knows no physics: it writes the flux of the model as Goose.
// The domain is a CSG shape of rect, circle and polygon, meshed with
// triangles by src/mesher/. The steps are `transient(end:, dt:, scheme:
// ssprk3)`. Each other construct gives an error `unsupported` at its token,
// with a fix.
//
// A `bc region: U = value` is a condition on an unknown (the state outside
// the side), and a `bc region: dot(F, n) = value` or `bc region: F*n =
// value` is a condition on the flux of an equation, the flux through the
// side. Its value is read at the state in the cell. A side with no
// condition has the flux of the state in the cell. A condition that names
// an unknown on a side of the flux gives the Rusanov flux of the state in
// the cell and the state outside (src/cons/solve.goose).
import std;
import math;
namespace cons_gen;
import syntax from "../syntax";
import codegen.codegen;
import mesher from "../mesher";

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

// The value of argument i of call n.
fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }
fn arg_count(n: i64) -> i64 { syntax::child_count(n) - 1 }

fn binary(n: i64, op: u8[:]) -> bool { kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == op }

// The source of n with no spaces.
fn compact(n: i64) -> u8[>..] {
    var s: u8[>..] = [];
    for c in source_of(n) { if c != ' ' { s.push(c); } }
    return s;
}

// --- the claim ---------------------------------------------------------------

// Whether this generator takes the model: a method fv on a physics that
// solves more than one unknown or a vector unknown, with no template.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_METHOD && syntax::text(syntax::nodes[c].aux) == "fv" && system_physics(prog, name(c)) { return true; }
        c = next(c);
    }
    return false;
}

fn system_physics(prog: i64, word: u8[:]) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS && name(c) == word {
            var unknowns = 0;
            var vector = false;
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_USES { return false; }
                if kind(u) == syntax::N_UNKNOWN {
                    unknowns++;
                    if kind(child(u, 0)) != syntax::N_TYPE_NAME { vector = true; }
                }
                u = next(u);
            }
            return unknowns > 1 || vector;
        }
        c = next(c);
    }
    return false;
}

// --- the model ----------------------------------------------------------------

let SHAPE = 0;
let BOUNDARY = 1;

struct Named {
    name: u8[..64],
    kind: i64,
    value: i64,     // the expression node
}

struct Unknown {
    name: u8[..64],
    rank: i64,
    offset: i64,
    initial: i64,   // the node of the start value, or -1
}

// ddt(U) + div(F) = S: the unknown, the flux node (or -1) and the node of S (or -1).
struct Eq {
    unknown: i64,
    flux: i64,
    rhs: i64,
}

let BC_STATE = 0;   // a value of an unknown, or of a component of it
let BC_FLUX = 1;    // the flux of the equation through the side

struct Bc {
    region: i64,    // the node of the region
    kind: i64,
    unknown: i64,
    comp: i64,      // the component, or -1 for all
    value: i64,
}

struct Step {
    name: u8[..64],
    end: f64,
    dt: f64,
    first_bc: i64,      // the bcs that hold in the step: a range of `active`
    bc_count: i64,
    output: bool,
}

struct Expect {
    node: i64,
    step: i64,
}

var named: Named[>..];
var unknowns: Unknown[>..];
var eqs: Eq[>..];
var local_names: u8[..64][>..];
var local_nodes: i64[>..];
var material_items: i64[>..];
var bcs: Bc[>..];
var active: i64[>..];
var global_bcs: i64[>..];
var steps: Step[>..];
var output_nodes: i64[>..];       // the fields of the outputs
var output_steps: i64[>..];       // the step of each
var expects: Expect[>..];
var reports: i64[>..];
var report_steps: i64[>..];
var aggregates: i64[>..];        // the integral and mean calls of the reports and expects
var cell_regions: i64[>..];      // the indices in `named` of the cells that an aggregate names
var physics_name: u8[..64];
var domain = -1;
var mesh_size = 0.0;
var comps = 0;
var limiter = 2;                 // cons::VAN_LEER
var material_node = -1;

fn reset() {
    codegen::reset();
    named.clear();
    unknowns.clear();
    eqs.clear();
    local_names.clear();
    local_nodes.clear();
    material_items.clear();
    bcs.clear();
    active.clear();
    global_bcs.clear();
    steps.clear();
    output_nodes.clear();
    output_steps.clear();
    expects.clear();
    reports.clear();
    report_steps.clear();
    aggregates.clear();
    cell_regions.clear();
    physics_name = [];
    domain = -1;
    mesh_size = 0.0;
    comps = 0;
    limiter = 2;
    material_node = -1;
    body.clear();
    memo_names.clear();
    memo_values.clear();
}

fn find_named(word: u8[:]) -> i64 {
    for g, i in named { if g.name == word { return i; } }
    return -1;
}

fn named_of(n: i64, k: i64) -> i64 {
    guard kind(n) == syntax::N_NAME else { return -1; }
    let i = find_named(name(n));
    return if i >= 0 && named[i].kind == k { i } else { -1 };
}

fn find_unknown(word: u8[:]) -> i64 {
    for u, i in unknowns { if u.name == word { return i; } }
    return -1;
}

recursive fn is_shape(n: i64) -> bool {
    if is_call(n, "rect") || is_call(n, "circle") || is_call(n, "polygon") { return true; }
    if kind(n) == syntax::N_NAME { return named_of(n, SHAPE) >= 0; }
    let op = if kind(n) == syntax::N_BINARY { syntax::text(tok(n)) } else { "" };
    return (op == "-" || op == "|" || op == "&") && is_shape(child(n, 0)) && is_shape(child(n, 1));
}

// `a | b` of boundaries.
recursive fn is_union(n: i64) -> bool {
    guard kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == "|" else { return false; }
    return (named_of(child(n, 0), BOUNDARY) >= 0 || is_union(child(n, 0))) && (named_of(child(n, 1), BOUNDARY) >= 0 || is_union(child(n, 1)));
}

fn value_item(n: i64) {
    let v = child(n, 1);
    var g = Named { name: [], kind: SHAPE, value: v };
    set64(g.name, name(n));
    if is_shape(v) {
        named.push(g);
    } else if is_call(v, "boundary") || is_union(v) {
        g.kind = BOUNDARY;
        named.push(g);
    } else {
        codegen::value_names.push(g.name);
        codegen::value_nodes.push(v);
    }
}

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    if syntax::text(tok(child(n, 0))) != "2d" { unsupported(tok(child(n, 0)), "the generator of conservation laws has 2d models", "write 2d in the header"); }
    let sys = child(n, 1);
    for s, i in codegen::SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            codegen::system = codegen::SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

let PHYSICS_FIX = "write physics flow on dom solves rho: Density, mom: Vector<Density*Velocity> { eq ddt(rho) + div(mom) = 0 ... }";

fn physics_item(n: i64) {
    if physics_name.len > 0 {
        unsupported(tok(n), "a second physics: the generator of conservation laws has one", "remove it");
        return;
    }
    set64(physics_name, name(n));
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            let r = child(c, 0);
            domain = named_of(r, SHAPE);
            if domain < 0 { unsupported(tok(r), "the region of the physics is not a shape that a let names", "write physics flow on dom, with let dom: Surface = rect(...)"); }
        } else if k == syntax::N_UNKNOWN {
            unknown_item(c);
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            material_node = -1;
            for i in material_items { if name(i) == name(m) { material_node = i; } }
            if material_node < 0 { unsupported(tok(m), "the material is not a material item", "write material gas { gamma: Real = 1.4 }"); }
            codegen::scope_material = material_node;
        } else if k == syntax::N_LET {
            local_names.push(default<u8[..64]>());
            set64(local_names[local_names.len - 1], name(c));
            local_nodes.push(child(c, 1));
        } else if k == syntax::N_EQN {
            equation(c);
        } else if k != syntax::N_ASSUME {
            unsupported(tok(c), "the physics statement is not in the generator of conservation laws", PHYSICS_FIX);
        }
        c = next(c);
    }
    if eqs.len != unknowns.len { unsupported(tok(n), "each unknown needs one eq with ddt of it", PHYSICS_FIX); }
}

// The rank of the unknown type: a scalar is 0, and Vector<..> is 1.
fn unknown_item(c: i64) {
    let t = child(c, 0);
    var u = Unknown { name: [], rank: 0, offset: comps, initial: child(c, 1) };
    set64(u.name, name(c));
    if kind(t) == syntax::N_TYPE_APPLY && name(t) == "Vector" {
        u.rank = 1;
    } else if kind(t) == syntax::N_TYPE_APPLY {
        unsupported(tok(t), "the unknown is a scalar or a Vector (a tensor has a flux of rank 3)", PHYSICS_FIX);
    }
    comps += if u.rank == 0 { 1 } else { 2 };
    unknowns.push(u);
}

// `eq ddt(U) + div(F) = S`.
fn equation(c: i64) {
    var e = Eq { unknown: -1, flux: -1, rhs: child(c, 1) };
    terms(child(c, 0), c, e);
    if e.unknown < 0 {
        unsupported(tok(c), "the eq has no ddt of an unknown", "write eq ddt(rho) + div(mom) = 0");
        return;
    }
    for q in eqs { if q.unknown == e.unknown { unsupported(tok(c), "a second eq for the same unknown", "write one eq with ddt of it"); } }
    eqs.push(e);
}

recursive fn terms(n: i64, eq_node: i64, e: Eq&) {
    if binary(n, "+") {
        terms(child(n, 0), eq_node, e);
        terms(child(n, 1), eq_node, e);
    } else if is_call(n, "ddt") && kind(arg(n, 0)) == syntax::N_NAME && find_unknown(name(arg(n, 0))) >= 0 && e.unknown < 0 {
        e.unknown = find_unknown(name(arg(n, 0)));
    } else if is_call(n, "div") && e.flux < 0 {
        e.flux = arg(n, 0);
    } else {
        unsupported(tok(n), str("the term '", source_of(n), "' is not in the generator of conservation laws"),
                    "write eq ddt(U) + div(F) = S, with the flux F and the source S on the right side");
    }
}

fn material(n: i64) { material_items.push(n); }

let MESH_FIX = "write mesh m on dom { size: 0.01 m, cell: tri }";

fn mesh_item(n: i64) {
    if mesh_size > 0.0 {
        unsupported(tok(n), "a second mesh: the generator of conservation laws has one", "remove the second mesh");
        return;
    }
    let region = child(n, 0);
    if domain < 0 || kind(region) != syntax::N_NAME || name(region) != named[domain].name {
        unsupported(tok(region), "the mesh is not on the region of the physics", MESH_FIX);
    }
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            mesh_size = codegen::eval(v).v[0];
        } else if kind(s) == syntax::N_SETTING && name(s) == "cell" && source_of(v) == "tri" {
            // the cell
        } else {
            unsupported(tok(s), "the mesh setting is not in the generator of conservation laws", "use only size and cell: tri");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", MESH_FIX); }
}

let LIMITER_NAMES = ["none", "minmod", "van_leer"];

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
                unsupported(tok(v), "the limiter of a system is none, minmod or van_leer", "write limiter: van_leer");
            }
        } else if name(s) == "flux" {
            if w != "upwind" { unsupported(tok(v), "the flux of a system is upwind (the Rusanov flux with the spectral radius of the Jacobian)", "write flux: upwind"); }
        } else if name(s) == "gradient" {
            if w != "green_gauss" { unsupported(tok(v), "the gradient of a system is green_gauss", "write gradient: green_gauss"); }
        } else if name(s) != "mesh" {
            unsupported(tok(s), "the method setting is not in the generator of conservation laws", "use mesh, flux: upwind, gradient: green_gauss and limiter");
        }
        s = next(s);
    }
}

// --- the steps ----------------------------------------------------------------------

let STEP_FIX = "write step march: transient(end: 2 s, dt: 0.6 ms, scheme: ssprk3) { ... }";

fn procedure(proc: i64, s: Step&) {
    if name(proc) != "transient" {
        unsupported(tok(proc), "the generator of conservation laws has transient steps", STEP_FIX);
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
        } else if w == "scheme" && source_of(v) == "ssprk3" {
            scheme = true;
        } else {
            unsupported(tok(a), "the generator of conservation laws has transient(end:, dt:, scheme: ssprk3)", STEP_FIX);
        }
        a = next(a);
    }
    if !scheme { unsupported(tok(proc), "the generator of conservation laws has the explicit scheme ssprk3", "add scheme: ssprk3"); }
    if s.end <= 0.0 || s.dt < 0.0 { unsupported(tok(proc), "the step needs an end after the start, and dt >= 0 or stable (the Rusanov step is not stable above a limit that skein run tells)", STEP_FIX); }
}

// The region of a bc: a boundary that a let names, or `a | b` of them.
recursive fn region_ok(r: i64) -> bool {
    if named_of(r, BOUNDARY) >= 0 { return true; }
    return binary(r, "|") && region_ok(child(r, 0)) && region_ok(child(r, 1));
}

let BC_FIX = "write bc inlet: rho = 1.4 kg/m^3, or bc wall: dot(mom, n) = 0, with the flux of the eq on the left";

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let target = child(n, 1);
    if !region_ok(region) {
        unsupported(tok(region), "the region of the bc is not a boundary that a let names", BC_FIX);
        return -1;
    }
    var b = Bc { region: region, kind: BC_STATE, unknown: -1, comp: -1, value: child(n, 2) };
    let t = compact(target);
    for u, i in unknowns {
        if t == u.name { b.unknown = i; }
        if u.rank == 1 && t == str(u.name, ".x") { b.unknown = i; b.comp = 0; }
        if u.rank == 1 && t == str(u.name, ".y") { b.unknown = i; b.comp = 1; }
    }
    if b.unknown < 0 {
        b.kind = BC_FLUX;
        for e in eqs {
            if e.flux < 0 { continue; }
            let f = compact(e.flux);
            if t == str("dot(", f, ",n)") || t == str(f, "*n") || (unknowns[e.unknown].rank == 0 && t == str("dot(n,", f, ")")) { b.unknown = e.unknown; }
        }
    }
    if b.unknown < 0 {
        unsupported(tok(target), str("the bc target '", t, "' is not an unknown or the flux of an eq"), BC_FIX);
        return -1;
    }
    if bcs.len >= 62 {
        unsupported(tok(n), "the generator of conservation laws has at most 62 bcs", "join bcs on the same region");
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

fn output_item(n: i64, s: Step&, si: i64) {
    s.output = true;
    var o = first(n);
    while o >= 0 {
        if kind(o) == syntax::N_EVERY {
            unsupported(tok(o), "a transient step writes its output at the end", "remove every");
        } else {
            output_nodes.push(o);
            output_steps.push(si);
        }
        o = next(o);
    }
}

fn step(n: i64) {
    var s = Step { name: [], end: 0.0, dt: 0.0, first_bc: 0, bc_count: 0, output: false };
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
            output_item(c, s, steps.len);
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
            unsupported(tok(c), "the step statement is not in the generator of conservation laws", "use only bc, output, report and expect");
        }
        c = next(c);
    }
    // The conditions of the step replace the global ones with the same region, target and kind.
    s.first_bc = active.len;
    for g in global_bcs {
        var replaced = false;
        for b in own { if same_target(bcs[b], bcs[g]) { replaced = true; } }
        if !replaced { active.push(g); }
    }
    active.append(own);
    s.bc_count = active.len - s.first_bc;
    steps.push(s);
}

fn same_target(a: Bc, b: Bc) -> bool {
    return compact(a.region) == compact(b.region) && a.kind == b.kind && a.unknown == b.unknown && a.comp == b.comp;
}

// --- the mini tensor algebra ---------------------------------------------------------

// A value of the model written as Goose: its rank (0 scalar, 1 vector, 2
// tensor) and the temporaries that hold its components.
struct T {
    rank: i64,
    v: i64[4],
}

var body: u8[>..];                    // the lets of the function that the emitter writes
var temps = 0;
var memo_names: u8[..64][>..];
var memo_values: T[>..];
var in_flux = false;                  // the function has no x, y and t
var depth = 0;

fn count_of(rank: i64) -> i64 { return if rank == 0 { 1 } else if rank == 1 { 2 } else { 4 }; }

// A new temporary with the Goose expression text.
fn temp(text: u8[:]) -> i64 {
    temps++;
    format(body, "    let t", temps, " = ", text, ";\n");
    return temps;
}

fn tn(id: i64) -> u8[>..] { return str("t", id); }

fn make(rank: i64) -> T { return T { rank: rank, v: [0, 0, 0, 0] }; }

fn scalar_temp(text: u8[:]) -> T {
    var r = make(0);
    r.v[0] = temp(text);
    return r;
}

fn constant(x: f64) -> T {
    var s: u8[>..] = [];
    lit(s, x);
    return scalar_temp(s);
}

fn rank_error(n: i64, what: u8[:]) -> T {
    unsupported(tok(n), str("'", source_of(n), "': ", what), "check the ranks: scalar, vector (x, y) and tensor");
    return constant(0.0);
}

fn combine(a: T, b: T, op: u8[:], n: i64) -> T {
    if a.rank != b.rank { return rank_error(n, str("the sides of '", op, "' have different ranks")); }
    var r = make(a.rank);
    for i in count_of(a.rank) { r.v[i] = temp(str(tn(a.v[i]), " ", op, " ", tn(b.v[i]))); }
    return r;
}

// a times b: a scalar by anything, a tensor by a vector or a tensor.
fn times(a: T, b: T, n: i64) -> T {
    if a.rank == 0 || b.rank == 0 {
        let s = if a.rank == 0 { a } else { b };
        let o = if a.rank == 0 { b } else { a };
        var r = make(o.rank);
        for i in count_of(o.rank) { r.v[i] = temp(str(tn(s.v[0]), " * ", tn(o.v[i]))); }
        return r;
    }
    return contract(a, b, n);
}

// The matrix product, the matrix by a vector, and the vector by a matrix. The dot of 2 vectors is a scalar.
fn contract(a: T, b: T, n: i64) -> T {
    if a.rank == 1 && b.rank == 1 {
        var r = make(0);
        r.v[0] = temp(str(tn(a.v[0]), " * ", tn(b.v[0]), " + ", tn(a.v[1]), " * ", tn(b.v[1])));
        return r;
    }
    if a.rank == 2 && b.rank == 1 {
        var r = make(1);
        for i in 2 { r.v[i] = temp(str(tn(a.v[2 * i]), " * ", tn(b.v[0]), " + ", tn(a.v[2 * i + 1]), " * ", tn(b.v[1]))); }
        return r;
    }
    if a.rank == 1 && b.rank == 2 {
        var r = make(1);
        for j in 2 { r.v[j] = temp(str(tn(a.v[0]), " * ", tn(b.v[j]), " + ", tn(a.v[1]), " * ", tn(b.v[2 + j]))); }
        return r;
    }
    if a.rank == 2 && b.rank == 2 {
        var r = make(2);
        for i in 2 {
            for j in 2 { r.v[2 * i + j] = temp(str(tn(a.v[2 * i]), " * ", tn(b.v[j]), " + ", tn(a.v[2 * i + 1]), " * ", tn(b.v[2 + j]))); }
        }
        return r;
    }
    return rank_error(n, "the product has a scalar or a rank that is not 1 or 2");
}

fn negate(a: T) -> T {
    var r = make(a.rank);
    for i in count_of(a.rank) { r.v[i] = temp(str("-", tn(a.v[i]))); }
    return r;
}

fn memo_find(word: u8[:]) -> i64 {
    for m, i in memo_names { if m == word { return i; } }
    return -1;
}

// The state u[k] of an unknown.
fn state_of(i: i64) -> T {
    let u = unknowns[i];
    var r = make(u.rank);
    for k in count_of(u.rank) { r.v[k] = temp(str("u[", u.offset + k, "]")); }
    return r;
}

fn name_value(n: i64) -> T {
    let w = name(n);
    let u = find_unknown(w);
    if u >= 0 { return state_of(u); }
    if w == "I" {
        var r = make(2);
        r.v = [temp("1.0"), temp("0.0"), temp("0.0"), temp("1.0")];
        return r;
    }
    if w == "n" {
        var r = make(1);
        r.v[0] = temp("nx");
        r.v[1] = temp("ny");
        return r;
    }
    if w == "x" || w == "y" || w == "t" {
        if in_flux { return rank_error(n, "the flux and the conditions on a flux do not read x, y and t"); }
        return scalar_temp(w);
    }
    if w == "pi" { return constant(PI); }
    let m = memo_find(w);
    if m >= 0 { return memo_values[m]; }
    var node = -1;
    for l, i in local_names { if l == w { node = local_nodes[i]; } }
    if node < 0 { node = codegen::find_value(w); }
    if node < 0 { node = codegen::property(w); }
    if node < 0 || depth > 100 {
        unsupported(tok(n), str("'", w, "' is not a value that the generator of conservation laws knows"), "declare it before, with let or param, or in the physics");
        return constant(0.0);
    }
    depth++;
    let r = emit(node);
    depth--;
    memo_names.push(default<u8[..64]>());
    set64(memo_names[memo_names.len - 1], w);
    memo_values.push(r);
    return r;
}

recursive fn emit(n: i64) -> T {
    let k = kind(n);
    if k == syntax::N_NUMBER || k == syntax::N_PERCENT { return constant(codegen::eval(n).v[0]); }
    if k == syntax::N_NAME { return name_value(n); }
    if k == syntax::N_VECTOR { return emit_vector(n); }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" { return negate(emit(first(n))); }
    if k == syntax::N_BINARY { return emit_binary(n); }
    if k == syntax::N_CALL { return emit_call(n); }
    if k == syntax::N_MEMBER { return emit_member(n); }
    return rank_error(n, "the expression is not in the generator of conservation laws");
}

recursive fn emit_vector(n: i64) -> T {
    var r = make(1);
    var f = 1.0;
    var i = 0;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNIT {
            f = codegen::to_system(c);
        } else if i < 2 {
            let x = emit(c);
            if x.rank != 0 { return rank_error(n, "a component of a vector is a scalar"); }
            r.v[i] = x.v[0];
            i++;
        } else {
            i++;
        }
        c = next(c);
    }
    if i != 2 { return rank_error(n, "a vector of this 2d model has 2 components"); }
    if f != 1.0 {
        var s: u8[>..] = [];
        lit(s, f);
        for j in 2 { r.v[j] = temp(str(tn(r.v[j]), " * ", s)); }
    }
    return r;
}

// `v.x` and `v.y` of a vector, and `T.xx` to `T.yy` of a tensor.
recursive fn emit_member(n: i64) -> T {
    let a = emit(first(n));
    let w = name(n);
    var comp = -1;
    if a.rank == 1 && w == "x" { comp = 0; }
    if a.rank == 1 && w == "y" { comp = 1; }
    if a.rank == 2 && w.len == 2 && (w[0] == 'x' || w[0] == 'y') && (w[1] == 'x' || w[1] == 'y') { comp = (if w[0] == 'x' { 0 } else { 2 }) + (if w[1] == 'x' { 0 } else { 1 }); }
    if comp < 0 { return rank_error(n, "the member is x or y of a vector, or xx to yy of a tensor"); }
    var r = make(0);
    r.v[0] = a.v[comp];
    return r;
}

recursive fn emit_binary(n: i64) -> T {
    let op = syntax::text(tok(n));
    let a = emit(child(n, 0));
    if op == "^" {
        let b = emit(child(n, 1));
        if a.rank != 0 || b.rank != 0 { return rank_error(n, "the power has scalars"); }
        let two = kind(child(n, 1)) == syntax::N_NUMBER && source_of(child(n, 1)) == "2";
        let text: u8[>..] = if two { str(tn(a.v[0]), " * ", tn(a.v[0])) } else { str("pow(", tn(a.v[0]), ", ", tn(b.v[0]), ")") };
        return scalar_temp(text);
    }
    let b = emit(child(n, 1));
    if op == "+" || op == "-" { return combine(a, b, op, n); }
    if op == "*" { return times(a, b, n); }
    if op == "/" {
        if b.rank != 0 { return rank_error(n, "the divisor is a scalar"); }
        var r = make(a.rank);
        for i in count_of(a.rank) { r.v[i] = temp(str(tn(a.v[i]), " / ", tn(b.v[0]))); }
        return r;
    }
    return rank_error(n, str("the operator '", op, "' is not in the generator of conservation laws"));
}

let SCALAR_FUNCTIONS = ["sqrt", "exp", "ln", "sin", "cos", "tan", "atan", "abs"];
let SCALAR_GOOSE = ["sqrt", "exp", "log", "sin", "cos", "tan", "atan", "abs"];

recursive fn emit_call(n: i64) -> T {
    let callee = child(n, 0);
    guard kind(callee) == syntax::N_NAME else { return rank_error(n, "the call is a function name"); }
    let w = name(callee);
    let count = arg_count(n);
    guard count >= 1 && count <= 2 else { return rank_error(n, "the function has 1 or 2 arguments"); }
    var args: T[2] = [make(0), make(0)];
    for i in count { args[i] = emit(arg(n, i)); }
    let f = position(SCALAR_FUNCTIONS, w);
    if f >= 0 && count == 1 {
        if args[0].rank != 0 { return rank_error(n, str(w, " has a scalar")); }
        return scalar_temp(str(SCALAR_GOOSE[f], "(", tn(args[0].v[0]), ")"));
    }
    if (w == "min" || w == "max") && count == 2 && args[0].rank == 0 && args[1].rank == 0 {
        return scalar_temp(str(w, "(", tn(args[0].v[0]), ", ", tn(args[1].v[0]), ")"));
    }
    if w == "dot" && count == 2 { return contract(args[0], args[1], n); }
    if w == "outer" && count == 2 && args[0].rank == 1 && args[1].rank == 1 {
        var r = make(2);
        for i in 2 { for j in 2 { r.v[2 * i + j] = temp(str(tn(args[0].v[i]), " * ", tn(args[1].v[j]))); } }
        return r;
    }
    if w == "norm" && count == 1 {
        var s: u8[>..] = [];
        for i in count_of(args[0].rank) { format(s, if i > 0 { " + " } else { "" }, tn(args[0].v[i]), " * ", tn(args[0].v[i])); }
        return scalar_temp(str("sqrt(", s, ")"));
    }
    if w == "transpose" && count == 1 && args[0].rank == 2 {
        var r = make(2);
        r.v = [args[0].v[0], args[0].v[2], args[0].v[1], args[0].v[3]];
        return r;
    }
    if w == "tr" && count == 1 && args[0].rank == 2 { return scalar_temp(str(tn(args[0].v[0]), " + ", tn(args[0].v[3]))); }
    if w == "sym" && count == 1 && args[0].rank == 2 {
        var r = make(2);
        let a = args[0];
        r.v = [a.v[0], temp(str("0.5 * (", tn(a.v[1]), " + ", tn(a.v[2]), ")")), 0, a.v[3]];
        r.v[2] = r.v[1];
        return r;
    }
    if w == "ddot" && count == 2 && args[0].rank == 2 && args[1].rank == 2 {
        var s: u8[>..] = [];
        for i in 4 { format(s, if i > 0 { " + " } else { "" }, tn(args[0].v[i]), " * ", tn(args[1].v[i])); }
        return scalar_temp(s);
    }
    return rank_error(n, "the function is not in the generator of conservation laws");
}

// --- the mesh file ------------------------------------------------------------------

// The point (x, y) of a vector node with its unit, in SI.
fn point2(n: i64) -> f64[2] {
    var p: f64[2] = [0.0, 0.0];
    var f = 1.0;
    var i = 0;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNIT { f = codegen::to_system(c); } else if i < 2 { p[i] = codegen::eval(c).v[0]; i++; }
        c = next(c);
    }
    return [p[0] * f, p[1] * f];
}

recursive fn mesher_shape(n: i64) -> i64 {
    if kind(n) == syntax::N_NAME { return mesher_shape(named[named_of(n, SHAPE)].value); }
    if is_call(n, "rect") {
        let a = point2(arg(n, 0));
        let b = point2(arg(n, 1));
        return mesher::rect(a[0], a[1], b[0], b[1]);
    }
    if is_call(n, "circle") {
        let c = point2(arg(n, 0));
        return mesher::circle(c[0], c[1], codegen::eval(arg(n, 1)).v[0]);
    }
    if is_call(n, "polygon") {
        var xy: f64[>..] = [];
        for i in arg_count(n) {
            let p = point2(arg(n, i));
            xy.push(p[0]);
            xy.push(p[1]);
        }
        return mesher::polygon(xy[..]);
    }
    let a = mesher_shape(child(n, 0));
    let b = mesher_shape(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { mesher::subtract(a, b) } else if op == "|" { mesher::union(a, b) } else { mesher::intersect(a, b) };
}

recursive fn mentions_xy(n: i64) -> bool {
    if kind(n) == syntax::N_NAME { return name(n) == "x" || name(n) == "y"; }
    var c = first(n);
    while c >= 0 {
        if mentions_xy(c) { return true; }
        c = next(c);
    }
    return false;
}

// The value of a scalar expression of x and y at a point.
recursive fn number_at(n: i64, x: f64, y: f64) -> f64 {
    if !mentions_xy(n) { return codegen::eval(n).v[0]; }
    let k = kind(n);
    if k == syntax::N_NAME { return if name(n) == "x" { x } else { y }; }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" { return -number_at(first(n), x, y); }
    if k == syntax::N_BINARY {
        let op = syntax::text(tok(n));
        let a = number_at(child(n, 0), x, y);
        let b = number_at(child(n, 1), x, y);
        if op == "+" { return a + b; }
        if op == "-" { return a - b; }
        if op == "*" { return a * b; }
        if op == "/" { return a / b; }
        if op == "^" { return pow(a, b); }
    }
    if is_call(n, "abs") { return abs(number_at(arg(n, 0), x, y)); }
    if is_call(n, "sqrt") { return sqrt(number_at(arg(n, 0), x, y)); }
    unsupported(tok(n), str("the expression '", source_of(n), "' of x and y is not in a region"), "write + - * / ^, abs and sqrt of x, y and values");
    return 0.0;
}

// Whether the condition n holds at (x, y), with the tolerance tol for ==, <= and >=.
recursive fn holds(n: i64, x: f64, y: f64, tol: f64) -> bool {
    if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "!" { return !holds(first(n), x, y, tol); }
    if kind(n) == syntax::N_BINARY {
        let op = syntax::text(tok(n));
        if op == "&&" { return holds(child(n, 0), x, y, tol) && holds(child(n, 1), x, y, tol); }
        if op == "||" { return holds(child(n, 0), x, y, tol) || holds(child(n, 1), x, y, tol); }
        let a = number_at(child(n, 0), x, y);
        let b = number_at(child(n, 1), x, y);
        if op == "==" { return abs(a - b) <= tol; }
        if op == "<" { return a < b - tol; }
        if op == "<=" { return a <= b + tol; }
        if op == ">" { return a > b + tol; }
        if op == ">=" { return a >= b - tol; }
    }
    unsupported(tok(n), str("the condition '", source_of(n), "' is not in a region"), "write a condition of x and y with == < <= > >= && || !");
    return false;
}

// Whether the point (x, y) is on the region r (a boundary, or `a | b`).
recursive fn on_region(r: i64, x: f64, y: f64, tol: f64) -> bool {
    if binary(r, "|") { return on_region(child(r, 0), x, y, tol) || on_region(child(r, 1), x, y, tol); }
    var v = named[named_of(r, BOUNDARY)].value;
    if binary(v, "|") { return on_region(v, x, y, tol); }
    if kind(v) == syntax::N_NAME { return on_region(v, x, y, tol); }
    let shape = arg(v, 0);
    let s = named_of(shape, SHAPE);
    if s >= 0 && s != domain {
        if abs(mesher::phi(mesher_shape(shape), x, y, 0.0)) > tol { return false; }
    }
    guard arg_count(v) >= 2 else { return true; }
    let cond = arg(v, 1);
    // A shape: the side is on its edge. Else a condition of x and y.
    let edge = named_of(cond, SHAPE);
    if edge >= 0 { return abs(mesher::phi(mesher_shape(cond), x, y, 0.0)) <= tol; }
    return holds(cond, x, y, tol);
}

var cells = 0;
var vertices = 0;

// Meshes the domain of the model that `read` read, and writes the mesh in
// the directory dir: the vertices, the triangles, the bcs that hold on each
// side, and the cells of each region. It gives an error, and false, if the
// mesher does not finish or a file cannot be written.
fn write_mesh(dir: u8[:]) -> bool {
    mesher::clear_shapes();
    mesher::clear_refines();
    mesher::clear_interfaces();
    let shape = mesher_shape(named[domain].value);
    guard mesher::mesh_2d(shape, mesh_size) else {
        write_stderr(str("skein run: the mesher cannot mesh '", named[domain].name, "' with the size of the mesh\n"));
        return false;
    };
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    for n in mesher::point_count() {
        let p = mesher::vertex(n);
        xy.push(p[0]);
        xy.push(p[1]);
    }
    for c in mesher::cell_count() {
        let a = mesher::cell_node(c, 0);
        let b = mesher::cell_node(c, 1);
        let d = mesher::cell_node(c, 2);
        let signed = (xy[2 * b] - xy[2 * a]) * (xy[2 * d + 1] - xy[2 * a + 1]) - (xy[2 * d] - xy[2 * a]) * (xy[2 * b + 1] - xy[2 * a + 1]);
        tri.push(a);
        tri.push(if signed >= 0.0 { b } else { d });
        tri.push(if signed >= 0.0 { d } else { b });
    }
    cells = tri.len / 3;
    vertices = xy.len / 2;
    let tol = 0.02 * mesh_size;
    var count = 3 * cells;
    var mask: i64[>..] = [];
    for t in count { mask.push(0); }
    var inner: u8[>..] = [];
    for t in count { inner.push(0); }
    mark_inner(tri, inner);
    for c in cells {
        for k in 3 {
            if inner[3 * c + k] != 0 { continue; }
            let a = tri[3 * c + k];
            let b = tri[3 * c + (k + 1) % 3];
            let mx = 0.5 * (xy[2 * a] + xy[2 * b]);
            let my = 0.5 * (xy[2 * a + 1] + xy[2 * b + 1]);
            for bi, i in bcs { if on_region(bi.region, mx, my, tol) { mask[3 * c + k] |= 1 << i; } }
        }
    }
    var flags: u8[>..] = [];
    for r in cell_regions {
        let region = mesher_shape(named[r].value);
        for c in cells {
            let x = (xy[2 * tri[3 * c]] + xy[2 * tri[3 * c + 1]] + xy[2 * tri[3 * c + 2]]) / 3.0;
            let y = (xy[2 * tri[3 * c] + 1] + xy[2 * tri[3 * c + 1] + 1] + xy[2 * tri[3 * c + 2] + 1]) / 3.0;
            flags.push(if mesher::phi(region, x, y, 0.0) < 0.0 { 1 } else { 0 });
        }
    }
    mesher::clear_shapes();
    if !ok() { return false; }
    if !write_file(str(dir, "/cons.xy"), to_bytes(xy[..])) || !write_file(str(dir, "/cons.tri"), to_bytes(tri[..])) ||
       !write_file(str(dir, "/cons.sides"), to_bytes(mask[..])) || !write_file(str(dir, "/cons.cells"), to_bytes(flags[..])) {
        write_stderr(str("skein run: cannot write the mesh in '", dir, "'\n"));
        return false;
    }
    return true;
}

// Marks the sides that two cells share.
fn mark_inner(tri: i64[:], inner: u8[>..]&) {
    var keys: i64[>..] = [];
    let nv = vertices;
    for c in tri.len / 3 {
        for k in 3 {
            let a = tri[3 * c + k];
            let b = tri[3 * c + (k + 1) % 3];
            keys.push((min(a, b) * nv + max(a, b)) * tri.len + 3 * c + k);
        }
    }
    sort(keys);
    for i in 1..keys.len {
        if keys[i] / tri.len == keys[i - 1] / tri.len {
            inner[keys[i] % tri.len] = 1;
            inner[keys[i - 1] % tri.len] = 1;
        }
    }
}

// The cells of the mesh, from the area of the domain and the size of the
// mesh: the mesher makes about 2 / (√3 h²/2) cells for each unit of area. The
// area is that of the bounding box of the shape, which is more.
fn dof_estimate() -> i64 {
    mesher::clear_shapes();
    let area = box_area(named[domain].value);
    mesher::clear_shapes();
    return (comps * area / (0.4330127 * mesh_size * mesh_size)) as! i64;
}

recursive fn box_area(n: i64) -> f64 {
    if kind(n) == syntax::N_NAME { return box_area(named[named_of(n, SHAPE)].value); }
    if is_call(n, "rect") {
        let a = point2(arg(n, 0));
        let b = point2(arg(n, 1));
        return abs((b[0] - a[0]) * (b[1] - a[1]));
    }
    if is_call(n, "circle") { return PI * codegen::eval(arg(n, 1)).v[0] * codegen::eval(arg(n, 1)).v[0]; }
    if is_call(n, "polygon") {
        var lo: f64[2] = [1.0e300, 1.0e300];
        var hi: f64[2] = [-1.0e300, -1.0e300];
        for i in arg_count(n) {
            let p = point2(arg(n, i));
            for j in 2 { lo[j] = min(lo[j], p[j]); hi[j] = max(hi[j], p[j]); }
        }
        return (hi[0] - lo[0]) * (hi[1] - lo[1]);
    }
    let a = box_area(child(n, 0));
    let b = box_area(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { a } else if op == "|" { a + b } else { min(a, b) };
}

// --- reading the model ----------------------------------------------------------------

// Whether a result expression holds integral or mean, and registers each.
recursive fn register_results(e: i64) {
    if is_call(e, "integral") || is_call(e, "mean") {
        let r = named_of(arg(e, 1), SHAPE);
        if arg_count(e) != 2 || r < 0 {
            unsupported(tok(e), "integral and mean have a field and a cell region (a Surface that a let names)", "write mean(p, behind), with let behind: Surface = polygon(...)");
            return;
        }
        if position(cell_regions, r) < 0 { cell_regions.push(r); }
        aggregates.push(e);
        return;
    }
    var c = first(e);
    while c >= 0 {
        register_results(c);
        c = next(c);
    }
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
        } else if k == syntax::N_MATERIAL {
            material(c);
        } else if k == syntax::N_PHYSICS {
            physics_nodes.push(c);
        } else if k == syntax::N_MESH {
            mesh_nodes.push(c);
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in the generator of conservation laws"),
                        "remove it: the generator runs one transient model of conservation laws");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    for n in physics_nodes { physics_item(n); }
    if physics_nodes.len == 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if !ok() { return false; }
    for m in mesh_nodes { mesh_item(m); }
    if mesh_size <= 0.0 && ok() { unsupported(at, "the model has no mesh", MESH_FIX); }
    for m in methods { method_item(m); }
    if !ok() { return false; }
    read_steps(prog, at);
    check_functions();
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
    for x in expects { register_results(child(x.node, 0)); }
    for r in reports { register_results(first(r)); }
}

// Writes the functions into a scratch, and evaluates the regions at a point,
// so that the errors of the model come from `read`.
fn check_functions() {
    var scratch: u8[>..] = [];
    emit_functions(scratch, "");
    guard ok() else { return; }
    mesher::clear_shapes();
    for b in bcs { on_region(b.region, 0.0, 0.0, 1.0e-9); }
    mesher::clear_shapes();
}

// --- the emitters --------------------------------------------------------------------

var normal_ok = true;     // n is the unit normal of a side in the function

fn begin(flux_only: bool, normal: bool) {
    body.clear();
    temps = 0;
    memo_names.clear();
    memo_values.clear();
    depth = 0;
    in_flux = flux_only;
    normal_ok = normal;
}

// Starts a block of a function with its own lets.
fn begin_block() {
    body.clear();
    memo_names.clear();
    memo_values.clear();
    depth = 0;
}

fn is_zero(n: i64) -> bool { return kind(n) == syntax::N_NUMBER && codegen::eval(n).v[0] == 0.0; }

// The value of node n with the rank wanted. A literal 0 is a zero of any rank.
fn emit_ranked(n: i64, rank: i64) -> T {
    if n < 0 || is_zero(n) {
        var z = make(rank);
        for i in count_of(rank) { z.v[i] = temp("0.0"); }
        return z;
    }
    let r = emit(n);
    if r.rank != rank {
        unsupported(tok(n), str("'", source_of(n), "' has the rank ", r.rank, ", and the unknown has ", rank), "give a scalar, a vector (x, y), with the rank of the unknown");
        return make(rank);
    }
    return r;
}

fn assign(out: u8[>..]&, target: u8[:], offset: i64, r: T) {
    for i in count_of(r.rank) { format(out, "    ", target, "[", offset + i, "] = ", tn(r.v[i]), ";\n"); }
}

fn emit_flux(out: u8[>..]&) {
    begin(true, true);
    var tail: u8[>..] = [];
    for e in eqs {
        let u = unknowns[e.unknown];
        if e.flux < 0 {
            for i in count_of(u.rank) { format(tail, "    f[", u.offset + i, "] = 0.0;\n"); }
            continue;
        }
        let t = emit(e.flux);
        if t.rank != u.rank + 1 {
            unsupported(tok(e.flux), str("the flux '", source_of(e.flux), "' of ", u.name, " has the rank ", t.rank), "give a vector flux to a scalar unknown, and a tensor flux to a vector unknown");
            continue;
        }
        for i in count_of(u.rank) { format(tail, "    f[", u.offset + i, "] = ", tn(t.v[2 * i]), " * nx + ", tn(t.v[2 * i + 1]), " * ny;\n"); }
    }
    format(out, "fn flux(u: f64[:], nx: f64, ny: f64, f: f64[:]) {\n", body, tail, "}\n\n");
}

fn emit_source(out: u8[>..]&) {
    begin(false, false);
    var tail: u8[>..] = [];
    for e in eqs {
        let u = unknowns[e.unknown];
        assign(tail, "s", u.offset, emit_ranked(e.rhs, u.rank));
    }
    format(out, "fn source(u: f64[:], x: f64, y: f64, t: f64, s: f64[:]) {\n", body, tail, "}\n\n");
}

fn emit_initial(out: u8[>..]&) {
    begin(false, false);
    var tail: u8[>..] = [];
    for u in unknowns {
        let t = emit_ranked(u.initial, u.rank);
        for i in count_of(u.rank) { format(tail, "    if k == ", u.offset + i, " { return ", tn(t.v[i]), "; }\n"); }
    }
    format(out, "fn initial(x: f64, y: f64, k: i64) -> f64 {\n    let t = 0.0;\n", body, tail, "    return 0.0;\n}\n\n");
}

// The state outside the sides and the flux of a step. Each condition is a
// function of its own, so that no generated function passes 60 lines (rule 4
// of the Power of Ten), and ghost_<step> and wall_<step> call those of the
// bits that the mask of the side sets.
fn emit_conditions(out: u8[>..]&, si: i64) {
    let s = steps[si];
    var ghost: u8[>..] = [];
    var wall: u8[>..] = [];
    begin(false, true);
    for ai in s.first_bc..s.first_bc + s.bc_count {
        let b = bcs[active[ai]];
        let u = unknowns[b.unknown];
        begin_block();
        let state = b.kind == BC_STATE;
        let rank = if !state { u.rank } else if b.comp >= 0 { 0 } else { u.rank };
        let t = emit_ranked(b.value, rank);
        let target = if state { "g" } else { "f" };
        format(out, "fn bc_", s.name, "_", active[ai], "(u: f64[:], x: f64, y: f64, t: f64, nx: f64, ny: f64, ", target, ": f64[:]) {\n", body);
        assign(out, target, u.offset + (if state { max(b.comp, 0) } else { 0 }), t);
        out.append("}\n\n");
        let call = str("    if ((mask >> ", active[ai], ") & 1) != 0 { bc_", s.name, "_", active[ai], "(u, x, y, t, nx, ny, ", target, "); }\n");
        if state { ghost.append(call); } else { wall.append(call); }
    }
    format(out, "fn ghost_", s.name, "(mask: i64, u: f64[:], x: f64, y: f64, t: f64, nx: f64, ny: f64, g: f64[:]) {\n    for k in ", comps, " { g[k] = u[k]; }\n", ghost, "}\n\n");
    format(out, "fn wall_", s.name, "(mask: i64, u: f64[:], x: f64, y: f64, t: f64, nx: f64, ny: f64, f: f64[:]) {\n", wall, "}\n\n");
}

fn emit_aggregates(out: u8[>..]&) {
    for e, i in aggregates {
        begin(false, false);
        let t = emit(arg(e, 0));
        if t.rank != 0 { unsupported(tok(e), "the field of integral and mean is a scalar", "write mean(norm(mom)/rho, behind)"); }
        format(out, "fn agg", i, "(u: f64[:], x: f64, y: f64) -> f64 {\n", body, "    return ", tn(t.v[0]), ";\n}\n\n");
    }
}

// The fields of the output of step si, as one function of the field and the component.
fn emit_output(out: u8[>..]&, si: i64) {
    let s = steps[si];
    begin(false, false);
    var tail: u8[>..] = [];
    var i = 0;
    for q, o in output_nodes {
        if output_steps[o] != si { continue; }
        if kind(q) != syntax::N_NAME {
            unsupported(tok(q), "the output is a name of an unknown or of a let of the physics", "write output rho, mom, p");
            continue;
        }
        let t = emit(q);
        for j in count_of(t.rank) { format(tail, "    if i == ", i, " && k == ", j, " { return ", tn(t.v[j]), "; }\n"); }
        i++;
    }
    format(out, "fn out_", s.name, "(u: f64[:], x: f64, y: f64, i: i64, k: i64) -> f64 {\n", body, tail, "    return 0.0;\n}\n\n");
}

// The names of the fields of the output, and their component counts.
fn output_names(si: i64, names: u8[>..]&, counts: u8[>..]&) {
    begin(false, false);
    var i = 0;
    for q, o in output_nodes {
        if output_steps[o] != si { continue; }
        if i > 0 {
            names.push(';');
            counts.append(", ");
        }
        names.append(name(q));
        let t = emit(q);
        format(counts, count_of(t.rank));
        i++;
    }
}

fn emit_functions(out: u8[>..]&, results_dir: u8[:]) {
    emit_flux(out);
    emit_source(out);
    emit_initial(out);
    emit_aggregates(out);
    for si in steps.len {
        emit_conditions(out, si);
        if steps[si].output { emit_output(out, si); }
    }
}

// --- results ------------------------------------------------------------------------------

// The Goose expression of the SI value of result e: a number, integral and
// mean of a field on a region, + - * / ^ and atan, tan, sqrt and abs of them.
recursive fn emit_result(out: u8[>..]&, e: i64) {
    if is_call(e, "integral") || is_call(e, "mean") {
        for a, i in aggregates {
            if a == e {
                format(out, "cons::", name(child(e, 0)), "(", position(cell_regions, named_of(arg(e, 1), SHAPE)), ", agg", i, ")");
                return;
            }
        }
    }
    if kind(e) == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(e))) >= 0 {
        out.push('(');
        emit_result(out, child(e, 0));
        format(out, " ", syntax::text(tok(e)), " ");
        emit_result(out, child(e, 1));
        out.push(')');
        return;
    }
    if kind(e) == syntax::N_UNARY && syntax::text(tok(e)) == "-" {
        out.append("(-");
        emit_result(out, first(e));
        out.push(')');
        return;
    }
    for f, i in SCALAR_FUNCTIONS {
        if is_call(e, f) && arg_count(e) == 1 {
            format(out, SCALAR_GOOSE[i], "(");
            emit_result(out, arg(e, 0));
            out.push(')');
            return;
        }
    }
    if !has_aggregate(e) {
        out.push('(');
        lit(out, codegen::eval(e).v[0]);
        out.push(')');
        return;
    }
    unsupported(tok(e), str("the result '", source_of(e), "' is not in the generator of conservation laws"), "write + - * / of integral(f, region), mean(f, region), values, and atan, tan, sqrt");
    out.append("0.0");
}

recursive fn has_aggregate(e: i64) -> bool {
    if is_call(e, "integral") || is_call(e, "mean") { return true; }
    var c = first(e);
    while c >= 0 {
        if has_aggregate(c) { return true; }
        c = next(c);
    }
    return false;
}

fn emit_report(out: u8[>..]&, step: u8[:], r: i64) {
    format(out, "    // report ", source_of(r), "\n    {\n        let name = ");
    quoted(out, source_of(r));
    out.append(";\n        report(");
    quoted(out, step);
    out.append(", name, quantity(");
    emit_result(out, first(r));
    out.append(", \"\"));\n    }\n");
}

// One check for an expect. The value is in SI, and the unit of the right
// side is for the text only.
fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    let want = codegen::eval(child(n, 1)).v[0];
    let within = child(n, 2);
    let relative = kind(within) == syntax::N_PERCENT;
    var tolerance = 0.0;
    if relative {
        let t = syntax::text(tok(within));
        let v, good = parse_flt(t[..t.len - 1]);
        tolerance = v / 100.0;
    } else {
        tolerance = codegen::eval(within).v[0];
    }
    let rhs = child(n, 1);
    let unit = if kind(rhs) == syntax::N_NUMBER && first(rhs) >= 0 { source_of(first(rhs)) } else { "" };
    let factor = if unit.len > 0 { codegen::to_system(first(rhs)) } else { 1.0 };
    var got: u8[>..] = [];
    emit_result(got, child(n, 0));
    format(out, "    // expect ", name(n), "\n    {\n        let name = ");
    quoted(out, name(n));
    format(out, ";\n        expect(name, ", got, ", ");
    lit(out, want);
    out.append(", ");
    lit(out, tolerance);
    format(out, ", ", if relative { "true" } else { "false" }, ", ");
    lit(out, factor);
    out.append(", ");
    quoted(out, unit);
    out.append(", ");
    quoted(out, str(source_of(rhs), " within ", source_of(within)));
    out.append(");\n    }\n");
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": transient, ssprk3\nfn step_", s.name, "() {\n    solved(");
    quoted(out, s.name);
    out.append(", cons::solve_transient(");
    lit(out, s.end);
    out.append(", ");
    lit(out, s.dt);
    format(out, ", flux, source, ghost_", s.name, ", wall_", s.name, "));\n");
    for i in 0..reports.len { if report_steps[i] == si { emit_report(out, s.name, reports[i]); } }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        var names: u8[>..] = [];
        var counts: u8[>..] = [];
        output_names(si, names, counts);
        var path: u8[>..] = [];
        quoted(path, str(results_dir, "/", s.name, ".vtu"));
        format(out, "    output(", path, ", cons::write_vtu(", path, ", ");
        quoted(out, names);
        format(out, ", [", counts, "], out_", s.name, "));\n");
    }
    out.append("}\n\n");
}

// The part of the generated program that is the same for each model.
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
    var failed = false;

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

    fn solved(step: u8[:], finished: bool) {
        format(solves, "  ", step, ",", cons::cell_count(), ",", cons::dof_count(), ",", cons::iterations, ",0,", cons::last_seconds, ",", if finished { "yes" } else { "no" }, "\n");
        solve_count++;
        check("solve", "finite", finished, str(cons::iterations, " steps"), "each value finite, dt below the stable step");
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, factor: f64, unit: u8[:], limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", abs(got - want) <= t, str(name, ": ", quantity(got / factor, unit)), limit);
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
        format(s, "solve[", solve_count, "]{step,cells,dofs,iterations,gmres_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

let LIMITER_CONSTANTS = ["cons::NONE", "cons::MINMOD", "cons::VAN_LEER"];

// Writes the Goose program of the model that `read` read. The program
// imports the cons package from `cons_path`, reads the mesh that
// `write_mesh` wrote in `results_dir`, and writes its VTU files there.
fn generate(cons_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    out.append("// The values are in SI.\nimport std;\nimport os;\nimport math;\nimport cons from ");
    quoted(out, cons_path);
    // The cons package imports the package solvers, beside it.
    out.append(";\nimport solvers from ");
    quoted(out, str(cons_path[..cons_path.len - 4], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    emit_functions(out, results_dir);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    format(out, "fn main() {\n    // mesh on ", named[domain].name, ": triangles, written by skein run\n    if !cons::read_mesh(");
    quoted(out, results_dir);
    format(out, ", ", cell_regions.len, ") {\n        write_stderr(\"skein run: cannot read the mesh in ");
    out.append(results_dir);
    format(out, "\\n\");\n        exit(3);\n    }\n    // physics ", physics_name, ": ", unknowns.len, " unknowns, ", comps, " components\n    cons::setup(", comps, ", ", LIMITER_CONSTANTS[limiter], ");\n    cons::start(initial);\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
