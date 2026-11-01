// The code generator of the optim models: a Skein syntax tree to a Goose
// program that calls src/optim/optim.goose.
//
// It takes the 2d models of the template linear_elastic with a mesh of
// `cell: quad` on a `rect`: bilinear quadrilaterals on a grid of equal
// rectangles, in plane stress or plane strain, with one material and one
// static step. The conditions are `u = 0`, `u = (a, b)`, `u.x = a` and
// `u.y = b` on a grid line `boundary(rect, x == a)`, on `boundary(rect)`
// or on a `point(p)` that is a node, and `sigma*n = (tx, ty)` on a side
// or at a point (a force). A report or an expect is a linear form of the
// values of u at points: `k*(u.y at p) + c`, with k and c constant. A
// `topology` study runs the loop of femahoy gap G1 (src/optim/topology.goose)
// on the objective, a report of this form. Each other construct gives an
// error `unsupported` at its token, with a fix.
//
// The values are in the unit system of the header.
import std;
import math;
namespace optim_gen;
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

// Whether this generator takes the model: a 2d model of linear_elastic
// with a mesh of cell: quad.
fn claims(prog: i64) -> bool {
    var quad = false;
    var elastic = false;
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_HEADER && syntax::text(tok(child(c, 0))) != "2d" { return false; }
        if kind(c) == syntax::N_MESH {
            var s = first(child(c, 1));
            while s >= 0 {
                if kind(s) == syntax::N_SETTING && name(s) == "cell" && source_of(first(s)) == "quad" { quad = true; }
                s = next(s);
            }
        }
        if kind(c) == syntax::N_PHYSICS {
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_USES && name(u) == "linear_elastic" { elastic = true; }
                u = next(u);
            }
        }
        c = next(c);
    }
    return quad && elastic;
}

let CURVE = 0;
let POINT = 1;

// A region of a bc: a grid line axis == value, the whole boundary
// (axis -1), or a point.
struct Region {
    name: u8[..64],
    kind: i64,
    axis: i64,
    at: f64[2],
}

// A bc: a fixed value of component comp (-1 for both), or a load.
struct Bc {
    region: i64,
    load: bool,
    comp: i64,
    v: f64[2],
}

// A linear form a u_comp(p) + b, and the unit of its value.
struct Form {
    a: f64,
    comp: i64,
    p: f64[2],
    b: f64,
    unit: u8[..64],
}

struct Report {
    item: i64,       // the N_REPORTITEM
    form: Form,
}

struct Expect {
    node: i64,
    form: Form,
    step: i64,       // -1 for an expect after the steps
}

var level = 0;   // skein verify convergence: the mesh has 2^level times the cells on each axis
var shape_names: u8[..64][>..];
var shape_nodes: i64[>..];
var regions: Region[>..];
var bcs: Bc[>..];
var global_bcs: i64[>..];
var step_bcs: i64[>..];
var reports: Report[>..];
var expects: Expect[>..];
var steps: u8[..64][>..];      // the name of the step: this generator has one
var output = false;
var material_items: i64[>..];
var material_node = -1;
var unknown: u8[..64];
var young = 0.0;
var poisson = 0.0;
var plane = 0;                 // optim::PLANE_STRESS
var thickness = 1.0;
var lo: f64[2] = [0.0, 0.0];
var hi: f64[2] = [0.0, 0.0];
var cells: i64[2] = [0, 0];
var domain = -1;

// The topology study of the run (spec §20.2), or -1: run.goose sets it.
var topology_study = -1;
var objective = -1;            // its report, in `reports`
var volume = 0.0;
var penalty = 3.0;
var filter = 0.0;
var most = 200;
var threshold = 0.0;           // 0: no solid

fn reset() {
    syntax::errors.clear();
    codegen::reset();
    shape_names.clear();
    shape_nodes.clear();
    regions.clear();
    bcs.clear();
    global_bcs.clear();
    step_bcs.clear();
    reports.clear();
    expects.clear();
    steps.clear();
    output = false;
    material_items.clear();
    material_node = -1;
    unknown = [];
    young = 0.0;
    poisson = 0.0;
    plane = 0;
    thickness = 1.0;
    domain = -1;
    objective = -1;
    volume = 0.0;
    penalty = 3.0;
    filter = 0.0;
    most = 200;
    threshold = 0.0;
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

fn point_value(n: i64) -> f64[2] {
    let v = value(n);
    if v.n != 2 { unsupported(tok(n), "a point of a 2d model has 2 components", "write (x, y) mm"); }
    return v.v;
}

fn find_name(names: u8[..64][:], word: u8[:]) -> i64 {
    for s, i in names { if s == word { return i; } }
    return -1;
}

// The exponents of length, mass and time of expression n, in halves as
// codegen::Unit, with u a Length. known is false for a construct that
// this generator does not follow.
recursive fn dims(n: i64) -> i64[3], bool {
    let k = kind(n);
    if k == syntax::N_NUMBER {
        guard first(n) >= 0 else { return [0, 0, 0], true; };
        let u = codegen::unit_value(first(n));
        return [u.l, u.m, u.t], true;
    }
    if k == syntax::N_AT { return [2, 0, 0], true; }
    if k == syntax::N_NAME {
        var v = codegen::find_value(name(n));
        if v < 0 { v = codegen::property(name(n)); }
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let d, known = dims(v);
            codegen::depth--;
            return d, known;
        }
        return [0, 0, 0], false;
    }
    if k == syntax::N_UNARY { return dims(first(n)); }
    if binary(n, "+") || binary(n, "-") { return dims(child(n, 0)); }
    if binary(n, "*") || binary(n, "/") {
        let a, ka = dims(child(n, 0));
        let b, kb = dims(child(n, 1));
        let s = if binary(n, "*") { 1 } else { -1 };
        return [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]], ka && kb;
    }
    return [0, 0, 0], false;
}

// The unit of the dimension d in the system of the header: the force to
// the power of the mass, times the length to the rest. "" if it is not so.
fn unit_text(d: i64[3], known: bool) -> u8[>..] {
    guard known && d[0] % 2 == 0 && d[1] % 2 == 0 && d[2] == -2 * d[1] else { return []; };
    let force = d[1] / 2;
    let length = (d[0] - d[1]) / 2;
    if force == 1 && length == -2 { return str(codegen::system.pressure); }
    var s: u8[>..] = [];
    for part in 2 {
        let e = if part == 0 { force } else { length };
        guard e != 0 else { continue; }
        if s.len > 0 { s.push('*'); }
        s.append(if part == 0 { codegen::system.force[..] } else { codegen::system.length[..] });
        if e != 1 { format(s, "^", e); }
    }
    return s;
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

// The linear form of expression n: a u_c(p) + b. False if n is not one.
recursive fn linear(n: i64, f: Form&) -> bool {
    if !has_unknown(n) {
        f.a = 0.0;
        f.comp = -1;
        f.b = scalar_value(n);
        return ok();
    }
    if kind(n) == syntax::N_AT {
        let m = child(n, 0);
        guard kind(m) == syntax::N_MEMBER && kind(first(m)) == syntax::N_NAME && name(first(m)) == unknown &&
              (name(m) == "x" || name(m) == "y") else { return false; };
        f.a = 1.0;
        f.comp = if name(m) == "x" { 0 } else { 1 };
        f.p = point_value(child(n, 1));
        f.b = 0.0;
        return true;
    }
    if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        guard linear(first(n), f) else { return false; };
        f.a = -f.a;
        f.b = -f.b;
        return true;
    }
    guard binary(n, "*") || binary(n, "/") || binary(n, "+") || binary(n, "-") else { return false; };
    var g = Form { a: 0.0, comp: -1, p: [0.0, 0.0], b: 0.0, unit: [] };
    guard linear(child(n, 0), f) && linear(child(n, 1), g) else { return false; };
    if binary(n, "+") || binary(n, "-") {
        let s = if binary(n, "-") { -1.0 } else { 1.0 };
        guard f.comp < 0 || g.comp < 0 || (f.comp == g.comp && f.p == g.p) else { return false; };
        if f.comp < 0 { f.comp = g.comp; f.p = g.p; }
        f.a += s * g.a;
        f.b += s * g.b;
        return true;
    }
    // One side is a constant: g, or f of a product.
    if f.comp < 0 {
        guard binary(n, "*") else { return false; };
        let k = f.b;
        f = g;
        f.a *= k;
        f.b *= k;
        return true;
    }
    guard g.comp < 0 else { return false; };
    let k = if binary(n, "/") { 1.0 / g.b } else { g.b };
    f.a *= k;
    f.b *= k;
    return true;
}

let FORM_FIX = "write a value of u at a point times a constant, for example -(u.y at (0 mm, 20 mm))*1 N";

// The form of expression n, with its unit, or an error at n.
fn form_of(n: i64) -> Form {
    var f = Form { a: 0.0, comp: -1, p: [0.0, 0.0], b: 0.0, unit: [] };
    if !linear(n, f) || f.comp < 0 {
        unsupported(tok(n), str("'", source_of(n), "' is not a linear form of u at a point"), FORM_FIX);
    }
    let d, known = dims(n);
    f.unit.append(unit_text(d, known));
    return f;
}

// --- the items ------------------------------------------------------------------

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

let REGION_FIX = "write let left: Curve = boundary(beam, x == 0 mm), or let corner: Point = point((60, 0) mm)";

// `let` and `param`: a rect, a boundary, a point, or a value.
fn value_item(n: i64) {
    let v = child(n, 1);
    var w: u8[..64] = [];
    set64(w, name(n));
    if is_call(v, "rect") {
        shape_names.push(w);
        shape_nodes.push(v);
    } else if is_call(v, "boundary") {
        boundary_region(w, v);
    } else if is_call(v, "point") {
        regions.push(Region { name: w, kind: POINT, axis: -1, at: point_value(arg(v, 0)) });
    } else {
        codegen::value_names.push(w);
        codegen::value_nodes.push(v);
    }
}

// `boundary(r)`, or `boundary(r, x == a)` or `y == a`: a grid line.
fn boundary_region(w: u8[..64], v: i64) {
    var r = Region { name: w, kind: CURVE, axis: -1, at: [0.0, 0.0] };
    if syntax::child_count(v) > 2 {
        let c = arg(v, 1);
        let side = if binary(c, "==") && kind(child(c, 0)) == syntax::N_NAME { name(child(c, 0)) } else { "" };
        if side != "x" && side != "y" {
            unsupported(tok(c), "a boundary of this generator is a grid line x == a or y == a", REGION_FIX);
            return;
        }
        r.axis = if side == "x" { 0 } else { 1 };
        r.at[0] = scalar_value(child(c, 1));
    }
    regions.push(r);
}

fn find_region(word: u8[:]) -> i64 {
    for r, i in regions { if r.name == word { return i; } }
    return -1;
}

let MESH_FIX = "write mesh m on beam { size: 1 mm, cell: quad }, with let beam: Surface = rect((0, 0) mm, (60, 20) mm)";

fn mesh(n: i64) {
    let r = child(n, 0);
    domain = if kind(r) == syntax::N_NAME { find_name(shape_names, name(r)) } else { -1 };
    if domain < 0 || kind(n) != syntax::N_MESH {
        unsupported(tok(r), "this generator meshes a rect", MESH_FIX);
        return;
    }
    let shape = shape_nodes[domain];
    let a = point_value(arg(shape, 0));
    let b = point_value(arg(shape, 1));
    lo = [min(a[0], b[0]), min(a[1], b[1])];
    hi = [max(a[0], b[0]), max(a[1], b[1])];
    var size = 0.0;
    var s = first(child(n, 1));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            size = scalar_value(first(s));
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "cell") {
            unsupported(tok(s), "the mesh setting is not in this generator", MESH_FIX);
        }
        s = next(s);
    }
    if size <= 0.0 { unsupported(tok(n), "the mesh needs a size", "add size: 1 mm"); }
    for j in 2 { cells[j] = max(1, round((hi[j] - lo[j]) / max(size, 1.0e-300)) as! i64) << level; }
}

let PHYSICS_FIX = "write physics plate on beam: linear_elastic solves u: Vector<Length> { material solid, assume plane_stress(1 mm) }";

fn physics_item(n: i64) {
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_UNKNOWN {
            set64(unknown, name(c));
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            for item in material_items { if name(item) == name(m) { material_node = item; } }
            if next(m) >= 0 || first(m) >= 0 { unsupported(tok(m), "this generator has one material on the whole body", "remove the second material"); }
        } else if k == syntax::N_ASSUME {
            var a = first(c);
            while a >= 0 {
                let w = name(a);
                if w == "plane_stress" || w == "plane_strain" {
                    plane = if w == "plane_stress" { 0 } else { 1 };
                    if first(a) >= 0 { thickness = scalar_value(first(a)); }
                } else if w != "small_strain" && w != "isotropic" {
                    unsupported(tok(a), str("the assumption '", w, "' is not in this generator"), PHYSICS_FIX);
                }
                a = next(a);
            }
        } else if k == syntax::N_USES && (name(c) != "linear_elastic" || first(c) >= 0) {
            unsupported(tok(c), "this generator has the template linear_elastic with no body force", PHYSICS_FIX);
        } else if k != syntax::N_ON && k != syntax::N_USES {
            unsupported(tok(c), "the physics statement is not in this generator", PHYSICS_FIX);
        }
        c = next(c);
    }
    if material_node < 0 { unsupported(tok(n), "the physics needs a material", PHYSICS_FIX); return; }
    codegen::scope_material = material_node;
    let e = codegen::property("E");
    let nu = codegen::property("nu");
    if e < 0 || nu < 0 { unsupported(tok(material_node), "the material needs E and nu", "give E: Pressure = 1 MPa and nu: Real = 0.3"); return; }
    young = scalar_value(e);
    poisson = scalar_value(nu);
    if young <= 0.0 || poisson <= -1.0 || poisson >= 0.5 { unsupported(tok(material_node), "the material needs E > 0 and -1 < nu < 0.5", "give E: Pressure = 1 MPa and nu: Real = 0.3"); }
}

let BC_FIX = "write bc left: u.x = 0, bc corner: u = 0, or bc tip: sigma*n = (0, -1) N";

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let r = child(n, 0);
    let target = child(n, 1);
    let v = child(n, 2);
    let ri = if kind(r) == syntax::N_NAME { find_region(name(r)) } else { -1 };
    if ri < 0 {
        unsupported(tok(r), "the region of a bc is a boundary or a point that a let names", REGION_FIX);
        return -1;
    }
    var b = Bc { region: ri, load: false, comp: -1, v: [0.0, 0.0] };
    let t = source_of(target);
    let member = kind(target) == syntax::N_MEMBER && kind(first(target)) == syntax::N_NAME && name(first(target)) == unknown;
    if t == unknown {
        let x = value(v);
        if x.n == 2 { b.v = x.v; } else if x.v[0] != 0.0 { unsupported(tok(v), "u is a vector", "write u = 0, or a vector"); }
    } else if member && (name(target) == "x" || name(target) == "y") {
        b.comp = if name(target) == "x" { 0 } else { 1 };
        b.v[0] = scalar_value(v);
    } else if t == "sigma*n" {
        b.load = true;
        b.v = point_value(v);
        let reg = regions[ri];
        if reg.kind == CURVE && reg.axis < 0 { unsupported(tok(r), "a traction of this generator is on one side of the rect", BC_FIX); }
    } else {
        unsupported(tok(target), str("the bc target '", t, "' is not in this generator"), BC_FIX);
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

let STEP_FIX = "write step load: static { bc ..., output u, report c = -(u.y at (0 mm, 20 mm))*1 N }";

fn step(n: i64) {
    if steps.len > 0 { unsupported(tok(n), "a second step: this generator has one static step", "remove the second step"); return; }
    let proc = child(n, 0);
    if name(proc) != "static" || first(proc) >= 0 { unsupported(tok(proc), "this generator has a static step", STEP_FIX); }
    var w: u8[..64] = [];
    set64(w, name(n));
    steps.push(w);
    var c = next(proc);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            let b = bc(c);
            if b >= 0 { step_bcs.push(b); }
        } else if k == syntax::N_OUTPUT {
            output = true;
            if kind(first(c)) != syntax::N_NAME || name(first(c)) != unknown || next(first(c)) >= 0 {
                unsupported(tok(c), "this generator writes the unknown, with the density of each cell", str("write output ", unknown));
            }
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                reports.push(Report { item: r, form: form_of(first(r)) });
                r = next(r);
            }
        } else if k == syntax::N_EXPECT {
            expect(c, 0);
        } else {
            unsupported(tok(c), "the step statement is not in this generator", "use only bc, output, report and expect");
        }
        c = next(c);
    }
}

// `expect e: form == v within w`.
fn expect(n: i64, step: i64) {
    expects.push(Expect { node: n, form: form_of(child(n, 0)), step: step });
}

// The arguments of the topology study (spec §20.2).
fn topology_args() {
    let pr = first(topology_study);
    var a = first(pr);
    while a >= 0 {
        let word = if syntax::nodes[a].aux >= 0 { syntax::text(syntax::nodes[a].aux) } else { "objective" };
        let v = first(a);
        if word == "objective" {
            for r, i in reports { if syntax::nodes[r.item].aux >= 0 && syntax::text(syntax::nodes[r.item].aux) == name(v) { objective = i; } }
            if objective < 0 { unsupported(tok(v), "the objective is not a named report of the step", "write report compliance = ... in the step"); }
        } else if word == "region" {
            if kind(v) != syntax::N_NAME || find_name(shape_names, name(v)) != domain {
                unsupported(tok(v), "the design region of this generator is the rect of the mesh", "write region: beam, the region of the mesh");
            }
        } else if word == "volume" {
            volume = scalar_value(v);
        } else if word == "penalty" {
            penalty = scalar_value(v);
        } else if word == "filter" {
            filter = scalar_value(v);
        } else if word == "max_iterations" {
            most = scalar_value(v) as! i64;
        } else if word == "threshold" {
            threshold = scalar_value(v);
        }
        a = next(a);
    }
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// it.
fn read(prog: i64) -> bool {
    reset();
    var physics_node = -1;
    var mesh_node = -1;
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
        } else if (k == syntax::N_MESH || k == syntax::N_POINTS) && mesh_node < 0 {
            mesh_node = c;
        } else if k != syntax::N_SCALES && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT && k != syntax::N_STUDY {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this generator"),
                        "remove it: this generator runs one physics on one mesh");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_node < 0 || mesh_node < 0 {
        unsupported(at, "the model needs a physics and a mesh", PHYSICS_FIX);
        return false;
    }
    physics_item(physics_node);
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
    if topology_study >= 0 && ok() { topology_args(); }
    return ok();
}

// The dofs of the mesh.
fn dof_estimate() -> i64 { return 2 * (cells[0] + 1) * (cells[1] + 1); }

// --- the Goose program ----------------------------------------------------------

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

    // A bc on a region that holds no node is a fault of the model.
    fn region(found: bool, name: u8[:]) {
        if !found {
            write_stderr(str("skein run: no node of the grid is on the region '", name, "'\n"));
            exit(1);
        }
    }

    // Component comp of u at (x, y). A point outside the grid is a fault of the model.
    fn at(comp: i64, x: f64, y: f64, name: u8[:]) -> f64 {
        let v, found = optim::value_at(comp, x, y);
        if !found {
            write_stderr(str("skein run: the point of '", name, "' is not in the grid\n"));
            failed = true;
        }
        return v;
    }

    fn solved(step: u8[:], regular: bool) {
        format(solves, "  ", step, ",", optim::cell_count(), ",", optim::dof_count(), ",", optim::solves, ",",
               optim::last_seconds, ",", if regular { "yes" } else { "no" }, "\n");
        solve_count++;
        if !regular {
            write_stderr(str("skein run: the conditions of step '", step, "' do not hold the body: the matrix is singular\n"));
            failed = true;
        }
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", abs(got - want) <= t, str(name, ": ", quantity(got, unit)), limit);
    }

    fn output(path: u8[:], name: u8[:]) {
        if !optim::write_vtu(path, name) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    // The history of a topology loop, one row for each iteration.
    fn history() -> u8[>..] {
        var s = str("topology[", optim::history_objective.len, "]{iteration,objective,volume,change}:\n");
        for v, i in optim::history_objective {
            format(s, "  ", i + 1, ",", v, ",", optim::history_volume[i], ",", optim::history_change[i], "\n");
        }
        return s;
    }

    fn print_records() {
        var s: u8[>..] = [];
        format(s, "report[", report_count, "]{step,expr,value}:\n", reports);
        format(s, "verify[", check_count, "]{stage,id,result,value,limit}:\n", checks);
        format(s, "solve[", solve_count, "]{step,cells,dofs,solves,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        if optim::history_objective.len > 0 { s.append(history()); }
        s.append(solvers::memory_records());
        write_stdout(s);
    }
""";

// The value of form f, as Goose: a at(comp, x, y, name) + b.
fn emit_form(out: u8[>..]&, f: Form) {
    lit(out, f.a);
    format(out, " * at(", f.comp, ", ");
    lit(out, f.p[0]);
    out.append(", ");
    lit(out, f.p[1]);
    out.append(", name) + ");
    lit(out, f.b);
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    let within = child(n, 2);
    format(out, "    // expect ", name(n), "\n    {\n        let name = ");
    quoted(out, name(n));
    out.append(";\n        expect(name, ");
    emit_form(out, x.form);
    out.append(", ");
    lit(out, scalar_value(child(n, 1)));
    out.append(", ");
    lit(out, scalar_value(within));
    format(out, ", ", if kind(within) == syntax::N_PERCENT { "true" } else { "false" }, ", ");
    quoted(out, x.form.unit);
    out.append(", ");
    quoted(out, str(source_of(child(n, 1)), " within ", source_of(within)));
    out.append(");\n    }\n");
}

// The calls of bc b.
fn emit_bc(out: u8[>..]&, b: Bc) {
    let r = regions[b.region];
    for comp in 2 {
        guard !b.load && (b.comp < 0 || b.comp == comp) else { continue; }
        let v = if b.comp < 0 { b.v[comp] } else { b.v[0] };
        if r.kind == POINT {
            format(out, "    region(optim::fix_point(");
            lit(out, r.at[0]);
            out.append(", ");
            lit(out, r.at[1]);
            format(out, ", ", comp, ", ");
        } else {
            for axis in 2 {
                guard r.axis < 0 || r.axis == axis else { continue; }
                for side in 2 {
                    guard r.axis < 0 || side == 0 else { continue; }
                    let at = if r.axis >= 0 { r.at[0] } else if side == 0 { lo[axis] } else { hi[axis] };
                    format(out, "    region(optim::fix_line(", axis, ", ");
                    lit(out, at);
                    format(out, ", ", comp, ", ");
                    lit(out, v);
                    format(out, "), \"", r.name, "\");\n");
                }
            }
            continue;
        }
        lit(out, v);
        format(out, "), \"", r.name, "\");\n");
    }
    guard b.load else { return; }
    format(out, "    region(optim::", if r.kind == POINT { str("point_load(") } else { str("edge_load(", r.axis, ", ") });
    lit(out, r.at[0]);
    out.append(", ");
    lit(out, if r.kind == POINT { r.at[1] } else { b.v[0] });
    out.append(", ");
    lit(out, if r.kind == POINT { b.v[0] } else { b.v[1] });
    if r.kind == POINT {
        out.append(", ");
        lit(out, b.v[1]);
    }
    format(out, "), \"", r.name, "\");\n");
}

fn emit_step(out: u8[>..]&, results_dir: u8[:]) {
    let s = steps[0];
    format(out, "// step ", s, "\nfn step_", s, "() {\n    optim::clear_bcs();\n");
    for b in global_bcs { emit_bc(out, bcs[b]); }
    for b in step_bcs { emit_bc(out, bcs[b]); }
    if topology_study >= 0 {
        let f = reports[objective].form;
        format(out, "    // the topology study ", name(topology_study), ": the objective ", source_of(reports[objective].item), "\n");
        format(out, "    if !optim::objective(", f.comp, ", ");
        lit(out, f.p[0]);
        out.append(", ");
        lit(out, f.p[1]);
        out.append(", ");
        lit(out, f.a);
        out.append(") { region(false, \"objective\"); }\n    solved(");
        quoted(out, s);
        out.append(", optim::topology(");
        lit(out, volume);
        out.append(", ");
        lit(out, penalty);
        out.append(", ");
        lit(out, filter);
        format(out, ", ", most, "));\n");
        emit_volume_check(out);
        if threshold > 0.0 {
            out.append("    // the solid: the cells at the threshold or more\n    solved(");
            quoted(out, s);
            out.append(", optim::solid(");
            lit(out, threshold);
            out.append("));\n");
        }
    } else {
        out.append("    solved(");
        quoted(out, s);
        out.append(", optim::solve());\n");
    }
    for r in reports {
        out.append("    {\n        let name = ");
        quoted(out, source_of(r.item));
        out.append(";\n        report(");
        quoted(out, s);
        out.append(", name, quantity(");
        emit_form(out, r.form);
        out.append(", ");
        quoted(out, r.form.unit);
        out.append("));\n    }\n");
    }
    for x in expects { if x.step == 0 { emit_expect(out, x); } }
    if output {
        out.append("    output(");
        quoted(out, str(results_dir, "/", s, ".vtu"));
        out.append(", ");
        quoted(out, unknown);
        out.append(");\n");
    }
    out.append("}\n\n");
}

// The check of the volume of the last design: within 1 % of the limit.
// Reason: the bisection of the update holds the volume to 1e-3 of the
// multiplier, and the last design is the one before the last update.
fn emit_volume_check(out: u8[>..]&) {
    out.append("    {\n        let k = optim::history_volume.len;\n        let v = if k > 0 { optim::history_volume[k - 1] } else { 0.0 };\n");
    out.append("        check(\"post\", \"volume\", abs(v - ");
    lit(out, volume);
    out.append(") <= 0.01 * ");
    lit(out, volume);
    out.append(", str(v), ");
    quoted(out, str(volume, " within 1%"));
    out.append(");\n    }\n");
}

// Writes the Goose program of the model that `read` read. The program
// imports the optim package from `optim_path`, and writes its VTU file in
// `results_dir`.
fn generate(optim_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    let src = optim_path[..optim_path.len - 5];
    out.append("// The values are in the unit system of the model header.\nimport std;\nimport os;\nimport math;\nimport optim from ");
    quoted(out, str(src, "optim"));
    // The package imports the package solvers, beside it.
    out.append(";\nimport solvers from ");
    quoted(out, str(src, "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    out.push('\n');
    emit_step(out, results_dir);
    format(out, "fn main() {\n    // mesh: ", cells[0], " x ", cells[1], " rectangles\n    optim::mesh_rect(");
    lit(out, lo[0]);
    out.append(", ");
    lit(out, lo[1]);
    out.append(", ");
    lit(out, hi[0]);
    out.append(", ");
    lit(out, hi[1]);
    format(out, ", ", cells[0], ", ", cells[1], ");\n    optim::elastic(");
    lit(out, young);
    out.append(", ");
    lit(out, poisson);
    out.append(", ");
    lit(out, thickness);
    format(out, ", ", if plane == 0 { "optim::PLANE_STRESS" } else { "optim::PLANE_STRAIN" }, ");\n");
    format(out, "    step_", steps[0], "();\n");
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
