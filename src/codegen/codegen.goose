// The code generator of the first slice (docs/implementation.md §6): a
// Skein syntax tree to a Goose program that calls src/runtime/api.goose.
//
// The slice is one linear static solid: a CSG box, a structured mesh of
// 27-node hexahedra, a linear elastic material, `u = 0` and `sigma*n`
// boundary conditions on the faces of the box, and `static` steps. Each
// other construct gives an error `unsupported` at its token, with a fix.
//
// The generator evaluates each value at compile time, in the unit system
// of the header, because the runtime has no units. The generated program
// prints the reports and the verification records in TOON
// (command_line.md §6.3, §6.4).
import std;
import math;
namespace codegen;
import syntax from "../syntax";

// --- errors ------------------------------------------------------------------

fn unsupported(tok: i64, msg: u8[:], fix: u8[:]) {
    syntax::error_at(tok, "unsupported", msg, fix);
}

fn ok() -> bool { syntax::errors.len == 0 }

fn kind(n: i64) -> i64 { syntax::nodes[n].kind }
fn tok(n: i64) -> i64 { syntax::nodes[n].tok }
fn name(n: i64) -> const u8[:] { syntax::text(syntax::nodes[n].tok) }

// The source text of node `n` and its children.
fn source_of(n: i64) -> const u8[:] {
    var lo = syntax::tokens[tok(n)].start;
    var hi = lo;
    span(n, lo, hi);
    return syntax::source[lo..hi];
}

recursive fn span(n: i64, lo: i64&, hi: i64&) {
    let t = syntax::tokens[tok(n)];
    if t.start < lo { lo = t.start; }
    if t.start + t.len > hi { hi = t.start + t.len; }
    if syntax::nodes[n].aux >= 0 && kind(n) != syntax::N_DIM_POW && kind(n) != syntax::N_UNIT_POW {
        let a = syntax::tokens[syntax::nodes[n].aux];
        if a.start + a.len > hi { hi = a.start + a.len; }
    }
    var c = syntax::nodes[n].first;
    while c >= 0 {
        span(c, lo, hi);
        c = syntax::nodes[c].next;
    }
    // A closing bracket ends a vector or a call.
    if kind(n) == syntax::N_VECTOR || kind(n) == syntax::N_CALL {
        var e = hi;
        while e < syntax::source.len && syntax::source[e] == ' ' { e++; }
        if e < syntax::source.len && syntax::source[e] == ')' { hi = e + 1; }
    }
}

// --- units (spec §6.3, §6.5) ---------------------------------------------------

// A factor to SI and the exponents of length, mass and time.
struct Unit {
    f: f64,
    l: i64,
    m: i64,
    t: i64,
}

struct BaseUnit {
    name: u8[..4],
    u: Unit,
}

let BASE_UNITS = [
    BaseUnit { name: "m", u: Unit { f: 1.0, l: 1, m: 0, t: 0 } },
    BaseUnit { name: "g", u: Unit { f: 1.0e-3, l: 0, m: 1, t: 0 } },
    BaseUnit { name: "s", u: Unit { f: 1.0, l: 0, m: 0, t: 1 } },
    BaseUnit { name: "N", u: Unit { f: 1.0, l: 1, m: 1, t: -2 } },
    BaseUnit { name: "Pa", u: Unit { f: 1.0, l: -1, m: 1, t: -2 } },
    BaseUnit { name: "t", u: Unit { f: 1000.0, l: 0, m: 1, t: 0 } },
    BaseUnit { name: "bar", u: Unit { f: 1.0e5, l: -1, m: 1, t: -2 } },
];

let PREFIX_CHARS = "numckMG";
let PREFIX_FACTORS = [1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2, 1.0e3, 1.0e6, 1.0e9];

// The unit system of the header: the size of its base units in SI, and the
// names that the reports print.
struct System {
    l: f64,
    m: f64,
    t: f64,
    length: u8[..4],
    force: u8[..4],
    pressure: u8[..4],
}

let SYSTEM_NAMES = ["SI", "mm_N_MPa_s", "mm_kg_ms"];
let SYSTEMS = [
    System { l: 1.0, m: 1.0, t: 1.0, length: "m", force: "N", pressure: "Pa" },
    System { l: 1.0e-3, m: 1000.0, t: 1.0, length: "mm", force: "N", pressure: "MPa" },
    System { l: 1.0e-3, m: 1.0, t: 1.0e-3, length: "mm", force: "kN", pressure: "GPa" },
];

var system = SYSTEMS[0];

fn find_base(word: u8[:]) -> i64 {
    for b, i in BASE_UNITS {
        if b.name == word { return i; }
    }
    return -1;
}

fn unit_name(t: i64) -> Unit {
    let word = syntax::text(t);
    let b = find_base(word);
    if b >= 0 { return BASE_UNITS[b].u; }
    if word.len > 1 {
        let p = position(PREFIX_CHARS, word[0]);
        let rest = find_base(word[1..]);
        if p >= 0 && rest >= 0 {
            var u = BASE_UNITS[rest].u;
            u.f *= PREFIX_FACTORS[p];
            return u;
        }
    }
    unsupported(t, str("the unit '", word, "' is not in this slice"),
                "use m, g, t, s, N, Pa or bar, with a prefix n, u, m, c, k, M or G, for example mm, kN or GPa");
    return Unit { f: 1.0, l: 0, m: 0, t: 0 };
}

recursive fn unit_value(n: i64) -> Unit {
    let k = kind(n);
    if k == syntax::N_UNIT { return unit_value(syntax::nodes[n].first); }
    if k == syntax::N_UNIT_NAME { return unit_name(tok(n)); }
    if k == syntax::N_UNIT_POW {
        let b = unit_value(syntax::nodes[n].first);
        let e = syntax::nodes[n].aux;
        return Unit { f: pow(b.f, e), l: b.l * e, m: b.m * e, t: b.t * e };
    }
    let a = unit_value(syntax::child(n, 0));
    let b = unit_value(syntax::child(n, 1));
    if k == syntax::N_UNIT_MUL { return Unit { f: a.f * b.f, l: a.l + b.l, m: a.m + b.m, t: a.t + b.t }; }
    return Unit { f: a.f / b.f, l: a.l - b.l, m: a.m - b.m, t: a.t - b.t };
}

// The factor from a unit to the unit system of the header.
fn to_system(n: i64) -> f64 {
    let u = unit_value(n);
    return u.f / (pow(system.l, u.l) * pow(system.m, u.m) * pow(system.t, u.t));
}

// --- values --------------------------------------------------------------------

// A scalar (`n` = 1) or a vector of 3 (`n` = 3), in the unit system of the
// header.
struct Val {
    n: i64,
    v: f64[3],
}

fn scalar(x: f64) -> Val { Val { n: 1, v: [x, 0.0, 0.0] } }

// The `let` and `param` items with a value, by name.
var value_names: u8[..64][>..];
var value_nodes: i64[>..];
var depth = 0;

fn find_value(word: u8[:]) -> i64 {
    for v, i in value_names {
        if v == word { return value_nodes[i]; }
    }
    return -1;
}

fn number(t: i64) -> f64 {
    var digits: u8[>..] = [];
    for c in syntax::text(t) {
        if c != '_' { digits.push(c); }
    }
    let x, good = parse_flt(digits);
    if !good { unsupported(t, "this number cannot be read", "write the number as digits, for example 2.5e3"); }
    return x;
}

// The value of expression `n`.
recursive fn eval(n: i64) -> Val {
    let k = kind(n);
    if k == syntax::N_NUMBER {
        if syntax::nodes[n].aux >= 0 {
            unsupported(tok(n), "temperatures are not in this slice", "remove the degC value");
        }
        var x = number(tok(n));
        if syntax::nodes[n].first >= 0 { x *= to_system(syntax::nodes[n].first); }
        return scalar(x);
    }
    if k == syntax::N_PERCENT {
        let t = syntax::text(tok(n));
        let x, good = parse_flt(t[..t.len - 1]);
        return scalar(x / 100.0);
    }
    if k == syntax::N_NAME {
        let v = find_value(name(n));
        if v < 0 {
            unsupported(tok(n), str("'", name(n), "' is not a value that this slice knows"),
                        "declare it before, with let or param and a number");
            return scalar(0.0);
        }
        if depth > 100 {
            unsupported(tok(n), str("the value of '", name(n), "' depends on itself"), "give it a number");
            return scalar(0.0);
        }
        depth++;
        let r = eval(v);
        depth--;
        return r;
    }
    if k == syntax::N_VECTOR {
        var r = Val { n: 3, v: [0.0, 0.0, 0.0] };
        var f = 1.0;
        var i = 0;
        var c = syntax::nodes[n].first;
        while c >= 0 {
            if kind(c) == syntax::N_UNIT {
                f = to_system(c);
            } else if i < 3 {
                let x = eval(c);
                if x.n != 1 { unsupported(tok(c), "a component of a vector is a vector", "write a scalar"); }
                r.v[i] = x.v[0];
                i++;
            } else {
                i++;
            }
            c = syntax::nodes[c].next;
        }
        if i != 3 { unsupported(tok(n), "a vector in this slice has 3 components", "write (x, y, z)"); }
        for j in 0..3 { r.v[j] *= f; }
        return r;
    }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        var r = eval(syntax::nodes[n].first);
        for j in 0..3 { r.v[j] = -r.v[j]; }
        return r;
    }
    if k == syntax::N_BINARY {
        let op = syntax::text(tok(n));
        let a = eval(syntax::child(n, 0));
        let b = eval(syntax::child(n, 1));
        if op == "+" || op == "-" {
            if a.n != b.n { unsupported(tok(n), "a scalar and a vector are added", "add two scalars or two vectors"); }
            var r = a;
            for j in 0..3 { r.v[j] = if op == "+" { a.v[j] + b.v[j] } else { a.v[j] - b.v[j] }; }
            return r;
        }
        if (op == "*" || op == "/") && b.n == 3 && (op == "/" || a.n == 3) {
            unsupported(tok(n), str("'", op, "' of two vectors, or of a scalar by a vector"),
                        "multiply or divide a vector by a scalar");
            return a;
        }
        if op == "*" {
            var r = if a.n == 3 { a } else { b };
            let s = if a.n == 3 { b.v[0] } else { a.v[0] };
            for j in 0..3 { r.v[j] *= s; }
            return r;
        }
        if op == "/" {
            var r = a;
            for j in 0..3 { r.v[j] /= b.v[0]; }
            return r;
        }
        if op == "^" && a.n == 1 && b.n == 1 { return scalar(pow(a.v[0], b.v[0])); }
    }
    unsupported(tok(n), str("the expression '", source_of(n), "' is not in this slice"),
                "write a number with a unit, a name, a vector (x, y, z), or + - * / ^ of these");
    return scalar(0.0);
}

// --- the model -------------------------------------------------------------------

// A face of the box: the part of its boundary where coordinate `axis` is
// `value`.
struct Face {
    name: u8[..64],
    axis: i64,
    value: f64,
}

// A boundary condition: `u = 0` (`fix`), or `sigma*n = t` on a face.
struct Bc {
    face: i64,      // index in `faces`
    fix: bool,
    t: f64[3],
    tok: i64,
}

struct Report {
    expr: i64,      // the expression node
    label: u8[..256],
}

struct Step {
    name: u8[..64],
    first_own: i64, // the bcs of the step itself: a range of `bcs`
    own_count: i64,
    first_bc: i64,  // all bcs that hold in the step: a range of `active`
    bc_count: i64,
    first_report: i64,
    report_count: i64,
    output: bool,
}

// An `expect` on a component of the displacement at a point.
struct Expect {
    name: u8[..64],
    comp: i64,      // 0, 1 or 2
    at: f64[3],
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
}

var model_name: u8[..64];
var solid_name: u8[..64];
var solid_tok = -1;
var lo: f64[3];
var hi: f64[3];
var faces: Face[>..];
var bcs: Bc[>..];
var active: Bc[>..];
var reports: Report[>..];
var steps: Step[>..];
var expects: Expect[>..];
// The skipped `expect` items: a name and why.
var skipped: u8[..256][>..];
var mesh_size = 0.0;
var mesh_tok = -1;
var order = 2;
var young = 0.0;
var poisson = -1.0;
var material_name: u8[..64];
var physics_name: u8[..64];
var unknown_name: u8[..64];
var physics_tok = -1;
// The `bc` items outside the steps, which hold in each step.
var global_bcs: i64[>..];

fn reset() {
    value_names.clear();
    value_nodes.clear();
    faces.clear();
    bcs.clear();
    active.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    skipped.clear();
    global_bcs.clear();
    model_name = [];
    solid_name = [];
    solid_tok = -1;
    material_name = [];
    physics_name = [];
    unknown_name = [];
    physics_tok = -1;
    mesh_size = 0.0;
    mesh_tok = -1;
    order = 2;
    young = 0.0;
    poisson = -1.0;
    depth = 0;
    system = SYSTEMS[0];
}

fn set64(into: u8[..64]&, s: u8[:]) {
    into = [];
    into.append(s[..min(s.len, 64)]);
}

fn set256(into: u8[..256]&, s: u8[:]) {
    into = [];
    into.append(s[..min(s.len, 256)]);
}

fn find_face(word: u8[:]) -> i64 {
    for f, i in faces {
        if f.name == word { return i; }
    }
    return -1;
}

fn axis_of(word: u8[:]) -> i64 {
    if word == "x" { return 0; }
    if word == "y" { return 1; }
    if word == "z" { return 2; }
    return -1;
}

let BOX_FIX = "write the solid as box((x0, y0, z0) mm, (x1, y1, z1) mm)";

// `let beam: Solid = box(a, b)`.
fn solid(n: i64, value: i64) {
    if kind(value) != syntax::N_CALL || name(syntax::child(value, 0)) != "box" || syntax::child_count(value) != 3 {
        unsupported(tok(value), str("the solid '", name(n), "' is not one box: this slice meshes only a box"), BOX_FIX);
        return;
    }
    if solid_tok >= 0 {
        unsupported(tok(n), "a second solid: this slice has one box", "remove the second solid");
        return;
    }
    let a = eval(syntax::nodes[syntax::child(value, 1)].first);
    let b = eval(syntax::nodes[syntax::child(value, 2)].first);
    if a.n != 3 || b.n != 3 { unsupported(tok(value), "box needs two corner points", BOX_FIX); return; }
    for j in 0..3 {
        lo[j] = min(a.v[j], b.v[j]);
        hi[j] = max(a.v[j], b.v[j]);
        if lo[j] == hi[j] { unsupported(tok(value), "the box has no volume", BOX_FIX); }
    }
    set64(solid_name, name(n));
    solid_tok = tok(n);
}

let FACE_FIX = "write the surface as boundary(beam, x == 0 mm), with x, y or z equal to a face of the box";

// `let root: Surface = boundary(beam, x == value)`.
fn surface(n: i64, value: i64) {
    let good = kind(value) == syntax::N_CALL && name(syntax::child(value, 0)) == "boundary" &&
               syntax::child_count(value) == 3;
    if !good {
        unsupported(tok(value), str("the surface '", name(n), "' is not one face of the box"), FACE_FIX);
        return;
    }
    let region = syntax::nodes[syntax::child(value, 1)].first;
    let cond = syntax::nodes[syntax::child(value, 2)].first;
    if kind(region) != syntax::N_NAME || name(region) != solid_name {
        unsupported(tok(region), "the surface is not on the box", FACE_FIX);
        return;
    }
    let is_face = kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" &&
                  kind(syntax::child(cond, 0)) == syntax::N_NAME;
    let axis = if is_face { axis_of(name(syntax::child(cond, 0))) } else { -1 };
    if axis < 0 {
        unsupported(tok(cond), "the surface is not one face of the box", FACE_FIX);
        return;
    }
    let v = eval(syntax::child(cond, 1));
    let tol = 1.0e-9 * (hi[axis] - lo[axis]);
    if abs(v.v[0] - lo[axis]) > tol && abs(v.v[0] - hi[axis]) > tol {
        unsupported(tok(cond), "the surface is not on a face of the box", FACE_FIX);
        return;
    }
    var f = Face { name: [], axis: axis, value: v.v[0] };
    set64(f.name, name(n));
    faces.push(f);
}

// `let` and `param`.
fn value_item(n: i64) {
    let ty = syntax::child(n, 0);
    let value = syntax::child(n, 1);
    let t = if kind(ty) == syntax::N_TYPE_NAME { name(ty) } else { "" };
    if t == "Solid" { solid(n, value); return; }
    if t == "Surface" { surface(n, value); return; }
    if t == "Volume" || t == "Curve" || t == "Region" {
        unsupported(tok(n), str("a ", t, " is not in this slice"), "use a Solid box and Surface faces of it");
        return;
    }
    var w: u8[..64] = [];
    set64(w, name(n));
    value_names.push(w);
    value_nodes.push(value);
}

fn material(n: i64) {
    if material_name.len > 0 {
        unsupported(tok(n), "a second material: this slice has one", "remove the second material");
        return;
    }
    set64(material_name, name(n));
    var p = syntax::nodes[n].first;
    while p >= 0 {
        let v = eval(syntax::child(p, 1));
        if name(p) == "E" {
            young = v.v[0];
        } else if name(p) == "nu" {
            poisson = v.v[0];
        } else if name(p) != "rho" {
            unsupported(tok(p), str("the property '", name(p), "' is not in this slice"),
                        "give the material only E and nu (and rho, which a static step does not use)");
        }
        p = syntax::nodes[p].next;
    }
    if young <= 0.0 || poisson <= -1.0 || poisson >= 0.5 {
        unsupported(tok(n), "the material needs E > 0 and -1 < nu < 0.5", "give E: Pressure = 210 GPa and nu: Real = 0.3");
    }
}

fn mesh(n: i64) {
    if mesh_tok >= 0 {
        unsupported(tok(n), "a second mesh: this slice has one", "remove the second mesh");
        return;
    }
    mesh_tok = tok(n);
    let region = syntax::child(n, 0);
    if kind(region) != syntax::N_NAME || name(region) != solid_name {
        unsupported(tok(region), "the mesh is not on the box", str("write mesh ", name(n), " on ", solid_name));
    }
    var s = syntax::nodes[syntax::child(n, 1)].first;
    while s >= 0 {
        let value = syntax::nodes[s].first;
        if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "refine and layers are not in this slice", "remove the line");
        } else if name(s) == "size" {
            mesh_size = eval(value).v[0];
        } else if name(s) == "cell" {
            if kind(value) != syntax::N_NAME || name(value) != "hex" {
                unsupported(tok(value), "this slice meshes a box with hex cells only", "write cell: hex");
            }
        } else {
            unsupported(tok(s), str("the mesh setting '", name(s), "' is not in this slice"), "use only size and cell");
        }
        s = syntax::nodes[s].next;
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 5 mm"); }
}

let PHYSICS_FIX = "write physics solid on beam: linear_elastic solves u: Vector<Length> { material steel }";

fn physics(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: this slice has one", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    var uses = false;
    var c = syntax::nodes[n].first;
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            let r = syntax::nodes[c].first;
            if kind(r) != syntax::N_NAME || name(r) != solid_name {
                unsupported(tok(r), "the physics is not on the box", PHYSICS_FIX);
            }
        } else if k == syntax::N_USES {
            uses = name(c) == "linear_elastic" && syntax::nodes[c].first < 0;
        } else if k == syntax::N_UNKNOWN {
            if unknown_name.len > 0 { unsupported(tok(c), "a second unknown", PHYSICS_FIX); }
            set64(unknown_name, name(c));
        } else if k == syntax::N_MATASSIGN {
            let m = syntax::nodes[c].first;
            if syntax::child_count(c) != 1 || syntax::nodes[m].first >= 0 || name(m) != material_name {
                unsupported(tok(c), "this slice gives one material to the whole physics", PHYSICS_FIX);
            }
        } else {
            unsupported(tok(c), "this slice has only the linear_elastic template, with no statements but material",
                        PHYSICS_FIX);
        }
        c = syntax::nodes[c].next;
    }
    if !uses { unsupported(tok(n), "this slice solves only the template linear_elastic", PHYSICS_FIX); }
}

fn method(n: i64) {
    if name(n) != physics_name || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this slice has only the cg method of the physics", "write method solid: cg { order: 2 }");
        return;
    }
    var s = syntax::nodes[syntax::child(n, 0)].first;
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            order = eval(syntax::nodes[s].first).v[0] as! i64;
            if order != 2 {
                unsupported(tok(s), "this slice has only order 2, because linear cells lock in bending",
                            "write order: 2");
            }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this slice", "use only mesh and order");
        }
        s = syntax::nodes[s].next;
    }
}

let BC_FIX = "write bc root: u = 0, or bc tip: sigma*n = (0, 0, -1) MPa";

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let region = syntax::child(n, 0);
    let target = syntax::child(n, 1);
    let value = syntax::child(n, 2);
    let face = if kind(region) == syntax::N_NAME { find_face(name(region)) } else { -1 };
    if face < 0 {
        unsupported(tok(region), "the region of the bc is not a face of the box that a let names", FACE_FIX);
        return -1;
    }
    var b = Bc { face: face, fix: false, t: [0.0, 0.0, 0.0], tok: tok(n) };
    if kind(target) == syntax::N_NAME && name(target) == unknown_name {
        let v = eval(value);
        b.fix = true;
        for j in 0..3 {
            if v.v[j] != 0.0 {
                unsupported(tok(value), "this slice fixes u = 0 only", BC_FIX);
            }
        }
    } else if source_of(target) == "sigma*n" {
        let v = eval(value);
        if v.n != 3 { unsupported(tok(value), "the traction is a vector", BC_FIX); }
        b.t = v.v;
    } else {
        unsupported(tok(target), str("the bc target '", source_of(target), "' is not in this slice"), BC_FIX);
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

fn step(n: i64) {
    let proc = syntax::child(n, 0);
    if name(proc) != "static" || syntax::nodes[proc].first >= 0 {
        unsupported(tok(proc), "this slice has only static steps with no arguments", "write step load: static { ... }");
    }
    var s = Step { name: [], first_own: bcs.len, own_count: 0, first_bc: 0, bc_count: 0,
                   first_report: reports.len, report_count: 0, output: false };
    set64(s.name, name(n));
    var c = syntax::nodes[proc].next;
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            bc(c);
        } else if k == syntax::N_OUTPUT {
            s.output = true;
            var o = syntax::nodes[c].first;
            while o >= 0 {
                if kind(o) == syntax::N_EVERY { unsupported(tok(o), "a static step has no every", "remove it"); }
                o = syntax::nodes[o].next;
            }
        } else if k == syntax::N_REPORT {
            var r = syntax::nodes[c].first;
            while r >= 0 {
                report(r);
                r = syntax::nodes[r].next;
            }
        } else {
            unsupported(tok(c), "the step statement is not in this slice", "use only bc, output and report");
        }
        c = syntax::nodes[c].next;
    }
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

let REPORT_FIX = "write report reaction(root), or report u.z at (100, 5, 5) mm";

// The kinds of a report: the reaction on a face, or u or one component of u
// at a point.
fn report(item: i64) {
    let e = syntax::nodes[item].first;
    if kind(e) == syntax::N_CALL && name(syntax::child(e, 0)) == "reaction" && syntax::child_count(e) == 2 {
        let r = syntax::nodes[syntax::child(e, 1)].first;
        let face = if kind(r) == syntax::N_NAME { find_face(name(r)) } else { -1 };
        if face < 0 { unsupported(tok(r), "reaction needs a face of the box that a let names", REPORT_FIX); return; }
    } else if kind(e) == syntax::N_AT {
        point_field(e);
    } else {
        unsupported(tok(e), str("the report '", source_of(e), "' is not in this slice"), REPORT_FIX);
        return;
    }
    var r = Report { expr: e, label: [] };
    set256(r.label, source_of(item));
    reports.push(r);
}

// The component of `u at p`: 0, 1 or 2, or 3 for the whole vector, or -1.
fn point_field(e: i64) -> i64 {
    let f = syntax::child(e, 0);
    var comp = -1;
    if kind(f) == syntax::N_NAME && name(f) == unknown_name { comp = 3; }
    if kind(f) == syntax::N_MEMBER && kind(syntax::nodes[f].first) == syntax::N_NAME &&
       name(syntax::nodes[f].first) == unknown_name {
        comp = axis_of(name(f));
    }
    if comp < 0 {
        unsupported(tok(f), str("'", source_of(f), "' at a point is not in this slice"),
                    str("use ", unknown_name, " or ", unknown_name, ".x, .y or .z at a point"));
        return -1;
    }
    let p = eval(syntax::child(e, 1));
    if p.n != 3 { unsupported(tok(e), "the point is not a vector", "write at (x, y, z) mm"); }
    return comp;
}

fn expect(n: i64) {
    let left = syntax::child(n, 0);
    let right = syntax::child(n, 1);
    let within = syntax::child(n, 2);
    var label: u8[>..] = [];
    format(label, source_of(right), " within ", source_of(within));
    // An expect on a field this slice cannot sample stays, and the report
    // says that it did not run.
    if kind(left) == syntax::N_AT && kind(syntax::child(left, 0)) == syntax::N_MEMBER &&
       name(syntax::nodes[syntax::child(left, 0)].first) != unknown_name {
        var s: u8[..256] = [];
        set256(s, str(name(n), ": ", source_of(left), " is not in this slice"));
        skipped.push(s);
        return;
    }
    if kind(left) != syntax::N_AT {
        unsupported(tok(left), "this slice checks only u.x, u.y or u.z at a point",
                    "write expect tip: u.z at (100, 5, 5) mm == -0.19 mm within 2%");
        return;
    }
    let comp = point_field(left);
    if comp == 3 {
        unsupported(tok(left), "an expect compares one component", str("write ", unknown_name, ".z"));
        return;
    }
    if comp < 0 { return; }
    var x = Expect { name: [], comp: comp, at: eval(syntax::child(left, 1)).v, want: eval(right).v[0],
                     within: eval(within).v[0], relative: kind(within) == syntax::N_PERCENT, label: [] };
    set64(x.name, name(n));
    set256(x.label, label);
    expects.push(x);
}

fn header(n: i64) {
    set64(model_name, name(n));
    let space = syntax::child(n, 0);
    if syntax::text(tok(space)) != "3d" {
        unsupported(tok(space), "this slice is 3d only", "write 3d in the header");
    }
    let sys = syntax::child(n, 1);
    for s, i in SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            system = SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

// Reads the program `prog` into the model above. False on an error, and
// `syntax::errors` holds the errors.
fn read(prog: i64) -> bool {
    reset();
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
            // The items after a wrong solid, surface or value give only
            // follow-on errors.
            if !ok() { return false; }
        } else if k == syntax::N_MATERIAL {
            material(c);
        } else if k == syntax::N_MESH {
            mesh(c);
        } else if k == syntax::N_PHYSICS {
            physics(c);
        } else if k == syntax::N_METHOD {
            method(c);
        } else if k == syntax::N_BC {
            let b = bc(c);
            if b >= 0 { global_bcs.push(b); }
        } else if k == syntax::N_STEP {
            step(c);
        } else if k == syntax::N_EXPECT {
            expect(c);
        } else if k != syntax::N_SCALES {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this slice"),
                        "remove it: this slice runs one linear static box");
        }
        c = syntax::nodes[c].next;
    }
    let first = syntax::nodes[prog].first;
    let at = if first >= 0 { tok(first) } else { tok(prog) };
    if solid_tok < 0 { unsupported(at, "the model has no box solid", BOX_FIX); }
    if mesh_tok < 0 { unsupported(at, "the model has no mesh", "add mesh m on beam { size: 5 mm, cell: hex }"); }
    if physics_tok < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if material_name.len == 0 { unsupported(at, "the model has no material", "add material steel { E: Pressure = 210 GPa, nu: Real = 0.3 }"); }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { bc root: u = 0 ... }"); }
    if ok() {
        step_bcs();
        supports();
    }
    return ok();
}

// Puts `b` in `set`, in place of a bc on the same face and target.
fn merge(set: Bc[>..]&, b: Bc) {
    for i in 0..set.len {
        if set[i].face == b.face && set[i].fix == b.fix {
            set[i] = b;
            return;
        }
    }
    set.push(b);
}

// The bcs that hold in each step (spec §12.2): the top-level bcs, wherever
// they are in the file, then the bcs of each step from that step on. A
// later bc replaces an earlier one on the same face and target.
fn step_bcs() {
    var set: Bc[>..] = [];
    for g in global_bcs { merge(set, bcs[g]); }
    for i in 0..steps.len {
        let s = steps[i];
        for b in bcs[s.first_own..s.first_own + s.own_count] { merge(set, b); }
        steps[i].first_bc = active.len;
        steps[i].bc_count = set.len;
        for b in set { active.push(b); }
    }
}

fn same_plane(a: i64, b: i64) -> bool {
    faces[a].axis == faces[b].axis && faces[a].value == faces[b].value
}

// The fixed faces of step `si`, one for each plane.
fn fixed_faces(si: i64) -> i64[>..] {
    let s = steps[si];
    var out: i64[>..] = [];
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if !b.fix { continue; }
        var seen = false;
        for f in out { if same_plane(f, b.face) { seen = true; } }
        if !seen { out.push(b.face); }
    }
    return out;
}

// The runtime keeps the supports from step to step, so each step must fix
// the same faces. The reaction sums each fixed face once, so two fixed
// faces must not share an edge.
fn supports() {
    let a = fixed_faces(0);
    if a.len == 0 {
        unsupported(physics_tok, "no face is fixed, so the box can move as a rigid body", "add bc root: u = 0");
        return;
    }
    for f in a[1..] {
        if faces[f].axis != faces[a[0]].axis {
            unsupported(bc_tok(0, f), str("the fixed faces ", faces[a[0]].name, " and ", faces[f].name,
                                          " share an edge, and this slice cannot sum their reaction"),
                        "fix one face, or two opposite faces");
        }
    }
    for si in 1..steps.len {
        let b = fixed_faces(si);
        if a != b {
            unsupported(bc_tok(si, b[b.len - 1]), str("step ", steps[si].name, " fixes a face that step ",
                                                      steps[0].name, " does not fix"),
                        "fix each face in the first step, or put the bc u = 0 outside the steps");
        }
    }
    for si in 0..steps.len {
        let s = steps[si];
        let fixed = fixed_faces(si);
        for r in reports[s.first_report..s.first_report + s.report_count] {
            if kind(r.expr) != syntax::N_CALL { continue; }
            let arg = syntax::nodes[syntax::child(r.expr, 1)].first;
            var on = false;
            for f in fixed { if same_plane(f, find_face(name(arg))) { on = true; } }
            if !on {
                unsupported(tok(arg), str("the face ", name(arg), " is not fixed in step ", s.name,
                                          ", so it has no reaction"),
                            "report the reaction of a face with bc u = 0");
            }
        }
    }
}

// The token of the fixed bc on face `f` in step `si`.
fn bc_tok(si: i64, f: i64) -> i64 {
    let s = steps[si];
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix && b.face == f { return b.tok; }
    }
    return physics_tok;
}

// --- the Goose program -------------------------------------------------------------

let AXES = ["runtime::X", "runtime::Y", "runtime::Z"];
let AXIS_NAMES = ["x", "y", "z"];

// The text of a Goose float literal.
fn lit(out: u8[>..]&, x: f64) {
    if x == 0.0 { out.append("0.0"); return; }
    format(out, x);
}

fn lits(out: u8[>..]&, v: f64[3]) {
    lit(out, v[0]);
    out.append(", ");
    lit(out, v[1]);
    out.append(", ");
    lit(out, v[2]);
}

// A Goose string literal.
fn quoted(out: u8[>..]&, s: u8[:]) {
    out.push('"');
    for c in s {
        if c == '"' || c == '\\' { out.push('\\'); }
        out.push(c);
    }
    out.push('"');
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

    fn quantity(x: f64, unit: u8[:]) -> u8[>..] { str(x, " ", unit) }

    fn vector(x: f64, y: f64, z: f64, unit: u8[:]) -> u8[>..] { str("(", x, ", ", y, ", ", z, ") ", unit) }

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

    fn check(id: u8[:], pass: bool, value: u8[:], limit: u8[:]) {
        checks.append("  post,");
        toon(checks, id);
        checks.append(if pass { ",pass," } else { ",fail," });
        toon(checks, value);
        checks.push(',');
        toon(checks, limit);
        checks.push('\n');
        check_count++;
        if !pass { failed = true; }
    }

    // A bc on a face that holds no node is a fault of the generator.
    fn face(found: bool, region: u8[:]) {
        if !found {
            write_stderr(str("skein run: internal fault: no node on the face '", region, "'. Please report it.\n"));
            exit(3);
        }
    }

    fn solve(step: u8[:]) {
        let converged = runtime::solve_static();
        format(solves, "  ", step, ",", runtime::dof_count(), ",", runtime::cg_iterations(), ",",
               runtime::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of the largest component of the applied load.
    fn equilibrium(rx: f64, ry: f64, rz: f64, unit: u8[:]) {
        let fx, fy, fz = runtime::applied_load();
        let r = sqrt((rx + fx) * (rx + fx) + (ry + fy) * (ry + fy) + (rz + fz) * (rz + fz));
        let limit = 1.0e-6 * max(abs(fx), max(abs(fy), abs(fz)));
        check("equilibrium", r <= limit, quantity(r, unit), quantity(limit, unit));
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let tol = if relative { within * abs(want) } else { within };
        check("expect", abs(got - want) <= tol, str(name, ": ", got, " ", unit), limit);
    }

    // An `expect` that this slice cannot evaluate.
    fn not_run(why: u8[:]) {
        checks.append("  post,expect,not_run,");
        toon(checks, why);
        checks.append(",\"\"\n");
        check_count++;
    }

    fn output(path: u8[:]) {
        if !runtime::write_vtu(path) {
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
        format(s, "solve[", solve_count, "]{step,dofs,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        write_stdout(s);
    }

    """;

fn emit_face_call(out: u8[>..]&, b: Bc) {
    let f = faces[b.face];
    out.append("    face(runtime::");
    out.append(if b.fix { "fix_face(" } else { "traction_face(" });
    out.append(AXES[f.axis]);
    out.append(", ");
    lit(out, f.value);
    if !b.fix {
        out.append(", ");
        lits(out, b.t);
    }
    out.append("), ");
    quoted(out, f.name);
    out.append(");\n");
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": static\n");
    format(out, "fn step_", s.name, "() {\n");
    if si > 0 {
        out.append("    // The supports of the first step stay.\n");
        out.append("    runtime::clear_loads();\n");
    }
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix && si > 0 { continue; }
        let f = faces[b.face];
        if b.fix {
            format(out, "    // bc ", f.name, ": ", unknown_name, " = 0\n");
        } else {
            format(out, "    // bc ", f.name, ": sigma*n = (", b.t[0], ", ", b.t[1], ", ", b.t[2], ") ",
                   system.pressure, "\n");
        }
        emit_face_call(out, b);
    }
    out.append("    solve(");
    quoted(out, s.name);
    out.append(");\n");
    // The reaction on all fixed faces, for the check `equilibrium`.
    out.append("    var rx, ry, rz = 0.0, 0.0, 0.0;\n");
    for fi in fixed_faces(si) {
        let f = faces[fi];
        format(out, "    {\n        let x, y, z = runtime::reaction_face(", AXES[f.axis], ", ");
        lit(out, f.value);
        out.append(");\n        rx += x;\n        ry += y;\n        rz += z;\n    }\n");
    }
    format(out, "    equilibrium(rx, ry, rz, \"", system.force, "\");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] {
        emit_report(out, s.name, r);
    }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(");\n");
    }
    out.append("}\n\n");
}

fn emit_point(out: u8[>..]&, at: f64[3]) {
    out.append("    let ux, uy, uz = runtime::displacement_at(");
    lits(out, at);
    out.append(");\n");
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n");
    out.append("    {\n    ");
    let e = r.expr;
    var value: u8[>..] = [];
    if kind(e) == syntax::N_CALL {
        let f = faces[find_face(name(syntax::nodes[syntax::child(e, 1)].first))];
        format(out, "    let x, y, z = runtime::reaction_face(", AXES[f.axis], ", ");
        lit(out, f.value);
        out.append(");\n");
        format(value, "vector(x, y, z, \"", system.force, "\")");
    } else {
        let comp = point_field(e);
        out.append("    ");
        emit_point(out, eval(syntax::child(e, 1)).v);
        if comp == 3 {
            format(value, "vector(ux, uy, uz, \"", system.length, "\")");
        } else {
            format(value, "quantity(u", AXIS_NAMES[comp], ", \"", system.length, "\")");
        }
    }
    out.append("        report(");
    quoted(out, step);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", value, ");\n    }\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n    ");
    emit_point(out, x.at);
    format(out, "        expect(");
    quoted(out, x.name);
    format(out, ", u", AXIS_NAMES[x.comp], ", ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"", system.length, "\", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}

// Writes the Goose program of the model that `read` read. The program
// imports the runtime from `runtime_path` and writes its VTU files under
// `results_dir`.
fn generate(runtime_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", system.length, ", force in ",
           system.force, ".\n");
    out.append("import std;\nimport os;\nimport runtime from ");
    quoted(out, runtime_path);
    out.append(";\n");
    out.append(PRELUDE);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n");
    format(out, "    // mesh on ", solid_name, ": 27-node hexahedra of at most ", mesh_size, " ",
           system.length, "\n");
    out.append("    runtime::mesh_box(");
    lits(out, lo);
    out.append(", ");
    lits(out, hi);
    for j in 0..3 {
        out.append(",\n                     runtime::cells_for(");
        lit(out, hi[j] - lo[j]);
        out.append(", ");
        lit(out, mesh_size);
        out.append(")");
    }
    out.append(");\n");
    format(out, "    // material ", material_name, "\n    runtime::material(");
    lit(out, young);
    out.append(", ");
    lit(out, poisson);
    out.append(");\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { emit_expect(out, x); }
    for s in skipped {
        out.append("    not_run(");
        quoted(out, s);
        out.append(");\n");
    }
    out.append("    print_records();\n");
    out.append("    if failed { exit(1); }\n");
    out.append("}\n");
}
