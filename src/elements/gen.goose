// The code generator of the simplex models: a Skein syntax tree to a Goose
// program that calls src/elements/api.goose.
//
// It takes the models that the box slice of src/codegen/ does not: a 2d
// model, a mesh of `tri` or `tet` cells, or a physics with a scalar
// unknown. The physics block decides the equations: the template
// `linear_elastic` or `heat(q_v)`, or the same statements written out
// (`hooke(E, nu, strain(u))` with `div(sigma) = 0`, `-k*grad(T)` or
// `fourier(k, T)` with `div(q) = q_v`). The geometry is a box or a rect,
// less a cylinder or a circle on its corner (lo.x, lo.y): the quarter
// model of a plate with a hole. Each other construct gives an error
// `unsupported` at its token, with a fix.
//
// This generator is a stand-in until the codegen dispatch tables of task
// skein-p2-interface land. It uses the units and the values of
// src/codegen/codegen.goose.
import std;
import math;
namespace elements_gen;
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
fn is_call(n: i64, callee: u8[:]) -> bool {
    kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == callee
}

// The value of argument i of call n.
fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }

// --- the claim ---------------------------------------------------------------

// Whether this generator takes the model: a 2d header, a mesh of tri or
// tet cells, or a physics with a scalar unknown.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER && syntax::text(tok(child(c, 0))) == "2d" { return true; }
        if k == syntax::N_MESH {
            var s = first(child(c, 1));
            while s >= 0 {
                if kind(s) == syntax::N_SETTING && name(s) == "cell" && kind(first(s)) == syntax::N_NAME &&
                   (name(first(s)) == "tri" || name(first(s)) == "tet") {
                    return true;
                }
                s = next(s);
            }
        }
        if k == syntax::N_PHYSICS {
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_UNKNOWN && kind(child(u, 0)) == syntax::N_TYPE_NAME { return true; }
                if kind(u) == syntax::N_USES && name(u) == "heat" { return true; }
                u = next(u);
            }
        }
        c = next(c);
    }
    return false;
}

// --- values --------------------------------------------------------------------

// A scalar or a vector of `dim` components, in the unit system of the
// header. A temperature in degC is an absolute temperature in K.
var dim = 3;

struct Vec {
    n: i64,
    v: f64[3],
}

fn scalar(x: f64) -> Vec { Vec { n: 1, v: [x, 0.0, 0.0] } }

recursive fn value(n: i64) -> Vec {
    let k = kind(n);
    if k == syntax::N_NUMBER && first(n) >= 0 {
        let u = first(first(n));
        if kind(u) == syntax::N_UNIT_NAME && syntax::text(tok(u)) == "degC" {
            let x = codegen::number(tok(n));
            let sign = if syntax::nodes[n].aux >= 0 { -1.0 } else { 1.0 };
            return scalar(sign * x + 273.15);
        }
    }
    if k == syntax::N_NAME {
        let v = codegen::find_value(name(n));
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let r = value(v);
            codegen::depth--;
            return r;
        }
    }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        var r = value(first(n));
        for j in 3 { r.v[j] = -r.v[j]; }
        return r;
    }
    if k == syntax::N_VECTOR {
        var r = Vec { n: 0, v: [0.0, 0.0, 0.0] };
        var f = 1.0;
        var c = first(n);
        while c >= 0 {
            if kind(c) == syntax::N_UNIT {
                f = codegen::to_system(c);
            } else {
                let x = value(c);
                if x.n != 1 { unsupported(tok(c), "a component of a vector is a vector", "write a scalar"); }
                if r.n < 3 { r.v[r.n] = x.v[0]; }
                r.n++;
            }
            c = next(c);
        }
        if r.n != dim {
            unsupported(tok(n), str("a vector of this ", dim, "d model has ", dim, " components"),
                        if dim == 2 { "write (x, y)" } else { "write (x, y, z)" });
        }
        for j in 3 { r.v[j] *= f; }
        return r;
    }
    let r = codegen::eval(n);
    return Vec { n: r.n, v: r.v };
}

fn scalar_value(n: i64) -> f64 {
    let v = value(n);
    if v.n != 1 { unsupported(tok(n), "the value is a vector", "write a scalar"); }
    return v.v[0];
}

// --- the model -------------------------------------------------------------------

let GEO_BOX = 0;
let GEO_HOLE = 1;

let ELASTIC = 0;
let HEAT = 1;

let SOLID_3D = 0;
let PLANE_STRESS = 1;
let PLANE_STRAIN = 2;

// The bc kinds: a fixed value of the unknown or a component of it, a
// traction `sigma*n` or a heat flux `dot(q, n)`.
let BC_FIX = 0;
let BC_LOAD = 1;

// A region of the boundary: the union of planes `axis == value`, or the
// surface of the hole.
struct Surf {
    name: u8[..64],
    hole: bool,
    axes: i64[..6],
    values: f64[..6],
}

struct Bc {
    surf: i64,      // index in `surfs`
    kind: i64,
    comp: i64,      // -1 for all components
    v: f64[3],
    tok: i64,
}

struct Report {
    expr: i64,
    label: u8[..256],
}

struct Step {
    name: u8[..64],
    first_own: i64,
    own_count: i64,
    first_bc: i64,
    bc_count: i64,
    first_report: i64,
    report_count: i64,
    output: bool,
}

// What a point field samples (api.goose: value_at, stress_at, flux_at).
let FIELD_VALUE = 0;
let FIELD_STRESS = 1;
let FIELD_FLUX = 2;

struct Expect {
    name: u8[..64],
    what: i64,
    comp: i64,
    at: f64[3],
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
}

var geometry_names: u8[..64][>..];
var geometry_nodes: i64[>..];
var surface_lets: i64[>..];
var geo = GEO_BOX;
var lo: f64[3];
var hi: f64[3];
var hole_r = 0.0;
var hole_name: u8[..64];
var body_name: u8[..64];
var surfs: Surf[>..];
var bcs: Bc[>..];
var active: Bc[>..];
var global_bcs: i64[>..];
var reports: Report[>..];
var steps: Step[>..];
var expects: Expect[>..];
var skipped: u8[..256][>..];
var mesh_size = 0.0;
var mesh_fine = 0.0;
var mesh_tok = -1;
var order = 2;
var physics = -1;
var plane = SOLID_3D;
var thickness = 1.0;
var physics_tok = -1;
var physics_name: u8[..64];
var unknown_name: u8[..64];
var flux_name: u8[..64];        // sigma or q
var material_name: u8[..64];
var material_node = -1;
var young = 0.0;
var poisson = -1.0;
var conductivity = 0.0;
var source = 0.0;
var source_node = -1;           // the right side of the eq, or -1 for 0

fn reset() {
    codegen::reset();
    dim = 3;
    geometry_names.clear();
    geometry_nodes.clear();
    surface_lets.clear();
    geo = GEO_BOX;
    hole_r = 0.0;
    hole_name = [];
    body_name = [];
    surfs.clear();
    bcs.clear();
    active.clear();
    global_bcs.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    skipped.clear();
    mesh_size = 0.0;
    mesh_fine = 0.0;
    mesh_tok = -1;
    order = 2;
    physics = -1;
    plane = SOLID_3D;
    thickness = 1.0;
    physics_tok = -1;
    physics_name = [];
    unknown_name = [];
    flux_name = [];
    material_name = [];
    material_node = -1;
    young = 0.0;
    poisson = -1.0;
    conductivity = 0.0;
    source = 0.0;
    source_node = -1;
}

fn set64(into: u8[..64]&, s: u8[:]) { codegen::set64(into, s); }

fn find_geometry(word: u8[:]) -> i64 {
    for g, i in geometry_names {
        if g == word { return geometry_nodes[i]; }
    }
    return -1;
}

fn find_surf(word: u8[:]) -> i64 {
    for s, i in surfs {
        if s.name == word { return i; }
    }
    return -1;
}

fn axis_of(word: u8[:]) -> i64 {
    if word == "x" { return 0; }
    if word == "y" { return 1; }
    if word == "z" && dim == 3 { return 2; }
    return -1;
}

// The geometry call that `n` is or names.
fn geometry_call(n: i64) -> i64 {
    if kind(n) == syntax::N_NAME { return find_geometry(name(n)); }
    return n;
}

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    let space = syntax::text(tok(child(n, 0)));
    if space != "2d" && space != "3d" {
        unsupported(tok(child(n, 0)), "this generator has 2d and 3d models", "write 2d or 3d in the header");
    }
    dim = if space == "2d" { 2 } else { 3 };
    let sys = child(n, 1);
    for s, i in codegen::SYSTEM_NAMES {
        if s == syntax::text(tok(sys)) {
            codegen::system = codegen::SYSTEMS[i];
            return;
        }
    }
    unsupported(tok(sys), "the unit system is not in this slice", "write SI, mm_N_MPa_s or mm_kg_ms");
}

// `let` and `param`: a geometry value, a boundary, or a value.
fn value_item(n: i64) {
    let v = child(n, 1);
    let geometric = kind(v) == syntax::N_CALL &&
        (is_call(v, "box") || is_call(v, "rect") || is_call(v, "cylinder") || is_call(v, "circle")) ||
        kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "-";
    var w: u8[..64] = [];
    set64(w, name(n));
    if geometric {
        geometry_names.push(w);
        geometry_nodes.push(v);
    } else if is_call(v, "boundary") {
        surface_lets.push(n);
    } else {
        codegen::value_names.push(w);
        codegen::value_nodes.push(v);
    }
}

let BODY_FIX = "write the body as box(p0, p1), or box(p0, p1) - cylinder(p0, (0, 0, 1), r, h) with the axis of the hole on the corner p0 (in 2d: rect and circle)";

fn corner(n: i64, i: i64) -> f64[3] {
    let v = value(arg(n, i));
    return v.v;
}

// The body of the physics: a box, or a box less a hole on its corner.
fn body(region: i64) {
    let n = if kind(region) == syntax::N_NAME { find_geometry(name(region)) } else { -1 };
    if n < 0 {
        unsupported(tok(region), "the region of the physics is not a solid that a let names", BODY_FIX);
        return;
    }
    set64(body_name, name(region));
    let shape = if dim == 3 { "box" } else { "rect" };
    var outer = n;
    var cut = -1;
    if kind(n) == syntax::N_BINARY {
        outer = geometry_call(child(n, 0));
        cut = child(n, 1);
    }
    if outer < 0 || !is_call(outer, shape) || syntax::child_count(outer) != 3 {
        unsupported(tok(n), str("the body is not one ", shape, ", or one ", shape, " less a hole"), BODY_FIX);
        return;
    }
    let a = corner(outer, 0);
    let b = corner(outer, 1);
    for j in 3 {
        lo[j] = min(a[j], b[j]);
        hi[j] = max(a[j], b[j]);
        if j < dim && lo[j] == hi[j] { unsupported(tok(outer), str("the ", shape, " has no size"), BODY_FIX); }
    }
    guard cut >= 0 else { return; }
    let c = geometry_call(cut);
    let round = if dim == 3 { "cylinder" } else { "circle" };
    if c < 0 || !is_call(c, round) {
        unsupported(tok(cut), str("the hole is not a ", round), BODY_FIX);
        return;
    }
    if kind(cut) == syntax::N_NAME { set64(hole_name, name(cut)); }
    let centre = corner(c, 0);
    let tol = 1.0e-6 * max(hi[0] - lo[0], hi[1] - lo[1]);
    if abs(centre[0] - lo[0]) > tol || abs(centre[1] - lo[1]) > tol {
        unsupported(tok(c), "the axis of the hole is not on the corner (lo.x, lo.y) of the body", BODY_FIX);
        return;
    }
    if dim == 3 {
        let axis = corner(c, 1);
        let r = scalar_value(arg(c, 2));
        let h = scalar_value(arg(c, 3));
        if abs(axis[0]) > 0.0 || abs(axis[1]) > 0.0 || axis[2] <= 0.0 || centre[2] > lo[2] + tol ||
           centre[2] + h < hi[2] - tol {
            unsupported(tok(c), "the hole does not go along z through the body", BODY_FIX);
            return;
        }
        hole_r = r;
    } else {
        hole_r = scalar_value(arg(c, 1));
    }
    if hole_r <= 0.0 || hole_r >= min(hi[0] - lo[0], hi[1] - lo[1]) {
        unsupported(tok(c), "the radius of the hole is not between 0 and the size of the body", BODY_FIX);
        return;
    }
    geo = GEO_HOLE;
}

let SURF_FIX = "write boundary(plate, x == 0 mm), with x, y or z on a face of the body, an || of these, or boundary(plate, bore) for the hole";

// Adds the planes of condition `cond` to `s`.
recursive fn surf_parts(cond: i64, s: Surf&) -> bool {
    if kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "||" {
        return surf_parts(child(cond, 0), s) && surf_parts(child(cond, 1), s);
    }
    if kind(cond) == syntax::N_NAME && hole_name.len > 0 && name(cond) == hole_name && geo == GEO_HOLE {
        s.hole = true;
        return true;
    }
    let is_plane = kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" &&
        kind(child(cond, 0)) == syntax::N_NAME;
    let axis = if is_plane { axis_of(name(child(cond, 0))) } else { -1 };
    if axis < 0 || s.axes.len == s.axes.cap {
        unsupported(tok(cond), "the boundary is not a face of the body", SURF_FIX);
        return false;
    }
    let v = scalar_value(child(cond, 1));
    let tol = 1.0e-6 * (hi[axis] - lo[axis]);
    if abs(v - lo[axis]) > tol && abs(v - hi[axis]) > tol {
        unsupported(tok(cond), "the plane is not a face of the body", SURF_FIX);
        return false;
    }
    s.axes.push(axis);
    s.values.push(v);
    return true;
}

// `let root: Surface = boundary(plate, x == 0 mm)`.
fn surface(n: i64) {
    let v = child(n, 1);
    if syntax::child_count(v) != 3 || kind(arg(v, 0)) != syntax::N_NAME || name(arg(v, 0)) != body_name {
        unsupported(tok(v), str("the boundary '", name(n), "' is not on the body of the physics"), SURF_FIX);
        return;
    }
    var s = Surf { name: [], hole: false, axes: [], values: [] };
    set64(s.name, name(n));
    if !surf_parts(arg(v, 1), s) { return; }
    if s.hole && s.axes.len > 0 {
        unsupported(tok(v), "a boundary is the hole or planes, not both", "make two boundaries");
        return;
    }
    surfs.push(s);
}

fn material(n: i64) {
    if material_node >= 0 {
        unsupported(tok(n), "a second material: this generator has one", "remove the second material");
        return;
    }
    material_node = n;
    set64(material_name, name(n));
}

// The material property `word`, or -1.
fn property(word: u8[:]) -> i64 {
    guard material_node >= 0 else { return -1; }
    var p = first(material_node);
    while p >= 0 {
        if name(p) == word { return child(p, 1); }
        p = next(p);
    }
    return -1;
}

fn read_material() {
    guard material_node >= 0 else { return; }
    var p = first(material_node);
    while p >= 0 {
        let w = name(p);
        let known = if physics == ELASTIC { w == "E" || w == "nu" } else { w == "k" || w == "q_v" };
        if !known && w != "rho" {
            unsupported(tok(p), str("the property '", w, "' is not in this generator"),
                        if physics == ELASTIC { "give the material E and nu" } else { "give the material k, and q_v for a source" });
        }
        p = next(p);
    }
    if physics == ELASTIC {
        let e = property("E");
        let nu = property("nu");
        if e >= 0 { young = scalar_value(e); }
        if nu >= 0 { poisson = scalar_value(nu); }
        if young <= 0.0 || poisson <= -1.0 || poisson >= 0.5 {
            unsupported(tok(material_node), "the material needs E > 0 and -1 < nu < 0.5",
                        "give E: Pressure = 210 GPa and nu: Real = 0.3");
        }
    } else {
        let k = property("k");
        if k >= 0 { conductivity = scalar_value(k); }
        if conductivity <= 0.0 {
            unsupported(tok(material_node), "the material needs k > 0", "give k: Conductivity = 50 W/(m*K)");
        }
    }
}

fn mesh(n: i64) {
    if mesh_tok >= 0 {
        unsupported(tok(n), "a second mesh: this generator has one", "remove the second mesh");
        return;
    }
    mesh_tok = tok(n);
    let region = child(n, 0);
    if kind(region) != syntax::N_NAME || name(region) != body_name {
        unsupported(tok(region), "the mesh is not on the body of the physics", str("write mesh ", name(n), " on ", body_name));
    }
    let cell = if dim == 3 { "tet" } else { "tri" };
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_REFINE {
            let r = first(s);
            let si = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
            if si < 0 || !surfs[si].hole {
                unsupported(tok(r), "this generator refines only the hole", "write refine hole: 0.1 mm");
            } else {
                mesh_fine = scalar_value(child(s, 1));
            }
        } else if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "layers are not in this generator", "remove the line");
        } else if name(s) == "size" {
            mesh_size = scalar_value(v);
        } else if name(s) == "cell" {
            if kind(v) != syntax::N_NAME || name(v) != cell {
                unsupported(tok(v), str("a ", dim, "d model of this generator has ", cell, " cells"),
                            str("write cell: ", cell));
            }
        } else {
            unsupported(tok(s), str("the mesh setting '", name(s), "' is not in this generator"),
                        "use only size, cell and refine");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 5 mm"); }
    if mesh_fine <= 0.0 { mesh_fine = mesh_size; }
}

let PHYSICS_FIX = "write physics solid on plate: linear_elastic solves u: Vector<Length> { material steel }, or physics thermal on bar: heat(q_v) solves T: Temperature { material steel }";

// A `let` of the flux: sigma = hooke(E, nu, strain(u)), or q = -k*grad(T)
// or fourier(k, T). It gives the physics, or -1.
fn flux_let(n: i64) -> i64 {
    let v = child(n, 1);
    let u = unknown_name;
    let hooke = is_call(v, "hooke") && syntax::child_count(v) == 4 && source_of(arg(v, 0)) == "E" &&
        source_of(arg(v, 1)) == "nu" && is_call(arg(v, 2), "strain") && source_of(arg(arg(v, 2), 0)) == u;
    if hooke { return ELASTIC; }
    let fourier = is_call(v, "fourier") && syntax::child_count(v) == 3 && source_of(arg(v, 0)) == "k" &&
        source_of(arg(v, 1)) == u;
    if fourier || source_of(v) == str("-k*grad(", u, ")") { return HEAT; }
    return -1;
}

fn physics_item(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: this generator has one", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    var eq_count = 0;
    var uses_template = false;
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            body(first(c));
        } else if k == syntax::N_USES {
            if name(c) == "linear_elastic" && first(c) < 0 {
                physics = ELASTIC;
                uses_template = true;
                set64(flux_name, "sigma");
            } else if name(c) == "heat" && syntax::child_count(c) <= 1 {
                physics = HEAT;
                uses_template = true;
                set64(flux_name, "q");
                if first(c) >= 0 { source_node = first(first(c)); }
            } else {
                unsupported(tok(c), "this generator has the templates linear_elastic and heat(q_v)", PHYSICS_FIX);
            }
        } else if k == syntax::N_UNKNOWN {
            if unknown_name.len > 0 { unsupported(tok(c), "a second unknown", PHYSICS_FIX); }
            set64(unknown_name, name(c));
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            if syntax::child_count(c) != 1 || first(m) >= 0 {
                unsupported(tok(c), "this generator gives one material to the whole physics", PHYSICS_FIX);
            }
        } else if k == syntax::N_ASSUME {
            assumptions(c);
        } else if k == syntax::N_LET {
            let p = flux_let(c);
            if p < 0 || physics >= 0 && p != physics {
                unsupported(tok(c), str("the let '", name(c), "' is not a flux of this generator"),
                            "write let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u)), or let q: Vector<HeatFlux> = -k*grad(T)");
            } else {
                physics = p;
                set64(flux_name, name(c));
            }
        } else if k == syntax::N_EQN {
            eq_count++;
            let left = child(c, 0);
            let right = child(c, 1);
            let div = is_call(left, "div") && flux_name.len > 0 && source_of(arg(left, 0)) == flux_name;
            if !div || physics == ELASTIC && source_of(right) != "0" {
                unsupported(tok(c), "this generator solves eq div(sigma) = 0 or eq div(q) = q_v", PHYSICS_FIX);
            } else if physics == HEAT {
                source_node = right;
            }
        } else {
            unsupported(tok(c), "the physics statement is not in this generator", PHYSICS_FIX);
        }
        c = next(c);
    }
    if physics < 0 { unsupported(tok(n), "the physics is neither linear elastic nor heat", PHYSICS_FIX); }
    if physics >= 0 && !uses_template && eq_count != 1 {
        unsupported(tok(n), "a written-out physics needs one eq", "add eq div(sigma) = 0 or eq div(q) = q_v");
    }
    if dim == 2 && physics == ELASTIC && plane == SOLID_3D {
        unsupported(tok(n), "a 2d elastic physics needs plane_stress(t) or plane_strain", "add assume plane_stress(1 mm)");
    }
}

fn assumptions(n: i64) {
    var a = first(n);
    while a >= 0 {
        let w = name(a);
        if w == "plane_stress" || w == "plane_strain" {
            if dim != 2 { unsupported(tok(a), str(w, " is for a 2d model"), "remove it"); }
            plane = if w == "plane_stress" { PLANE_STRESS } else { PLANE_STRAIN };
            if first(a) >= 0 { thickness = scalar_value(first(a)); }
            if first(a) < 0 && w == "plane_stress" {
                unsupported(tok(a), "plane_stress needs the thickness", "write plane_stress(1 mm)");
            }
        } else if w != "small_strain" && w != "isotropic" {
            unsupported(tok(a), str("the assumption '", w, "' is not in this generator"),
                        "use small_strain, isotropic, plane_stress(t) or plane_strain");
        }
        a = next(a);
    }
}

// The source of the heat equation: a value, or the material property q_v.
fn heat_source() {
    guard physics == HEAT else { return; }
    if source_node < 0 {
        let p = property("q_v");
        if p >= 0 { source = scalar_value(p); }
        return;
    }
    if kind(source_node) == syntax::N_NAME && codegen::find_value(name(source_node)) < 0 {
        let p = property(name(source_node));
        if p < 0 {
            unsupported(tok(source_node), str("'", name(source_node), "' is not a value or a property of the material"),
                        "give the source as a let or as q_v in the material");
            return;
        }
        source = scalar_value(p);
        return;
    }
    source = scalar_value(source_node);
}

fn method(n: i64) {
    if name(n) != physics_name || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this generator has only the cg method of the physics", "write method solid: cg { order: 2 }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            order = scalar_value(first(s)) as! i64;
            if order != 1 && order != 2 { unsupported(tok(s), "this generator has order 1 and order 2", "write order: 2"); }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this generator", "use only mesh and order");
        }
        s = next(s);
    }
}

let BC_FIX_TEXT = "write bc root: u = 0, bc sym: u.x = 0, bc tip: sigma*n = (0, 0, -1) MPa, bc wall: T = 20 degC or bc side: dot(q, n) = 1000 W/m^2";

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let target = child(n, 1);
    let v = child(n, 2);
    let s = if kind(region) == syntax::N_NAME { find_surf(name(region)) } else { -1 };
    if s < 0 {
        unsupported(tok(region), "the region of the bc is not a boundary that a let names", SURF_FIX);
        return -1;
    }
    var b = Bc { surf: s, kind: BC_FIX, comp: -1, v: [0.0, 0.0, 0.0], tok: tok(n) };
    let t = source_of(target);
    let member = kind(target) == syntax::N_MEMBER && kind(first(target)) == syntax::N_NAME &&
        name(first(target)) == unknown_name;
    if t == unknown_name {
        let x = value(v);
        if physics == ELASTIC && x.n == 1 && x.v[0] != 0.0 {
            unsupported(tok(v), "u is a vector", "write u = 0, or a vector");
        }
        if physics == HEAT && x.n != 1 { unsupported(tok(v), "T is a scalar", "write T = 20 degC"); }
        b.v = x.v;
    } else if member && physics == ELASTIC && axis_of(name(target)) >= 0 {
        b.comp = axis_of(name(target));
        b.v[0] = scalar_value(v);
    } else if physics == ELASTIC && t == str(flux_name, "*n") {
        let x = value(v);
        if x.n != dim { unsupported(tok(v), "the traction is a vector", BC_FIX_TEXT); }
        b.kind = BC_LOAD;
        b.v = x.v;
    } else if physics == HEAT && (t == str("dot(", flux_name, ", n)") || t == str(flux_name, "*n")) {
        b.kind = BC_LOAD;
        b.v[0] = scalar_value(v);
    } else {
        unsupported(tok(target), str("the bc target '", t, "' is not in this generator"), BC_FIX_TEXT);
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

fn step(n: i64) {
    let proc = child(n, 0);
    if name(proc) != "static" || first(proc) >= 0 {
        unsupported(tok(proc), "this generator has only static steps with no arguments", "write step load: static { ... }");
    }
    var s = Step { name: [], first_own: bcs.len, own_count: 0, first_bc: 0, bc_count: 0,
                   first_report: reports.len, report_count: 0, output: false };
    set64(s.name, name(n));
    var c = next(proc);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            bc(c);
        } else if k == syntax::N_OUTPUT {
            s.output = true;
            var o = first(c);
            while o >= 0 {
                if kind(o) == syntax::N_EVERY { unsupported(tok(o), "a static step has no every", "remove it"); }
                o = next(o);
            }
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                report(r);
                r = next(r);
            }
        } else {
            unsupported(tok(c), "the step statement is not in this generator", "use only bc, output and report");
        }
        c = next(c);
    }
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

let REPORT_FIX = "write report reaction(root), or report u.x at (1, 0, 0) mm, sigma.xx at (0, 1, 0) mm or T at (0, 0, 0) m";

fn report(item: i64) {
    let e = first(item);
    if is_call(e, "reaction") && syntax::child_count(e) == 2 {
        let r = arg(e, 0);
        if kind(r) != syntax::N_NAME || find_surf(name(r)) < 0 {
            unsupported(tok(r), "reaction needs a boundary that a let names", REPORT_FIX);
            return;
        }
    } else if kind(e) == syntax::N_AT {
        let what, comp = point_field(e, true);
        if what < 0 { return; }
    } else {
        unsupported(tok(e), str("the report '", source_of(e), "' is not in this generator"), REPORT_FIX);
        return;
    }
    var r = Report { expr: e, label: [] };
    codegen::set256(r.label, source_of(item));
    reports.push(r);
}

// The field of `f at p`: what it samples and the component, with -1 for
// the whole vector u. It gives -1, -1 for a field it cannot sample, and an
// error if `errors`.
fn point_field(e: i64, errors: bool) -> i64, i64 {
    let f = child(e, 0);
    var what = -1;
    var comp = -1;
    if kind(f) == syntax::N_NAME && name(f) == unknown_name {
        what = FIELD_VALUE;
        comp = if physics == HEAT { 0 } else { -1 };
    }
    if kind(f) == syntax::N_MEMBER && kind(first(f)) == syntax::N_NAME {
        let object = name(first(f));
        let m = name(f);
        if object == unknown_name && physics == ELASTIC && axis_of(m) >= 0 {
            what = FIELD_VALUE;
            comp = axis_of(m);
        } else if object == flux_name && physics == HEAT && axis_of(m) >= 0 {
            what = FIELD_FLUX;
            comp = axis_of(m);
        } else if object == flux_name && physics == ELASTIC && m.len == 2 {
            let i = if m[0] == 'z' { 2 } else { axis_of(m[..1]) };
            let j = if m[1] == 'z' { 2 } else { axis_of(m[1..]) };
            if i >= 0 && j >= 0 {
                what = FIELD_STRESS;
                comp = 3 * i + j;
            }
        }
    }
    if what < 0 {
        if errors {
            unsupported(tok(f), str("'", source_of(f), "' at a point is not in this generator"), REPORT_FIX);
        }
        return -1, -1;
    }
    let p = value(child(e, 1));
    if p.n != dim && errors { unsupported(tok(e), "the point is not a vector", "write at (x, y, z) mm"); }
    return what, comp;
}

fn expect(n: i64) {
    let left = child(n, 0);
    let right = child(n, 1);
    let within = child(n, 2);
    var label: u8[>..] = [];
    format(label, source_of(right), " within ", source_of(within));
    var what = -1;
    var comp = -1;
    if kind(left) == syntax::N_AT {
        what, comp = point_field(left, false);
    }
    // An expect on a field this generator cannot sample stays, and the
    // report says that it did not run.
    if what < 0 || comp < 0 {
        var s: u8[..256] = [];
        codegen::set256(s, str(name(n), ": ", source_of(left), " is not in this generator"));
        skipped.push(s);
        return;
    }
    var x = Expect { name: [], what: what, comp: comp, at: value(child(left, 1)).v, want: scalar_value(right),
                     within: scalar_value(within), relative: kind(within) == syntax::N_PERCENT, label: [] };
    // A tolerance of a temperature is a difference, not a degC value.
    if kind(within) == syntax::N_NUMBER && first(within) >= 0 && source_of(first(within)) == "degC" {
        x.within -= 273.15;
    }
    set64(x.name, name(n));
    codegen::set256(x.label, label);
    expects.push(x);
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// the errors.
fn read(prog: i64) -> bool {
    reset();
    var physics_node = -1;
    var mesh_node = -1;
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
            if physics_node >= 0 { unsupported(tok(c), "a second physics: this generator has one", "remove the second physics"); }
            physics_node = c;
        } else if k == syntax::N_MESH {
            if mesh_node >= 0 { unsupported(tok(c), "a second mesh: this generator has one", "remove the second mesh"); }
            mesh_node = c;
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k != syntax::N_SCALES && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this generator"),
                        "remove it: this generator runs one linear static or steady heat model");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_node < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if !ok() { return false; }
    physics_item(physics_node);
    if !ok() { return false; }
    for s in surface_lets { surface(s); }
    if material_node < 0 { unsupported(at, "the model has no material", "add material steel { E: Pressure = 210 GPa, nu: Real = 0.3 }"); }
    read_material();
    heat_source();
    if mesh_node < 0 { unsupported(at, "the model has no mesh", "add mesh m on plate { size: 5 mm }"); } else { mesh(mesh_node); }
    for m in methods { method(m); }
    if !ok() { return false; }
    c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            let b = bc(c);
            if b >= 0 { global_bcs.push(b); }
        } else if k == syntax::N_STEP {
            step(c);
        } else if k == syntax::N_EXPECT {
            expect(c);
        }
        c = next(c);
    }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { bc root: u = 0 ... }"); }
    if ok() { step_bcs(); }
    if ok() { supports(); }
    return ok();
}

// Puts `b` in `set`, in place of a bc on the same region and target.
fn merge(set: Bc[>..]&, b: Bc) {
    for i in 0..set.len {
        if set[i].surf == b.surf && set[i].kind == b.kind && set[i].comp == b.comp {
            set[i] = b;
            return;
        }
    }
    set.push(b);
}

// The bcs that hold in each step (spec §12.2), as in src/codegen.
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

// Each step needs a fixed value, and each reaction report a fixed region.
fn supports() {
    for s in steps {
        var fixes = 0;
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_FIX { fixes++; } }
        if fixes == 0 {
            unsupported(physics_tok, str("step ", s.name, " fixes no value, so the solution is not unique"),
                        if physics == ELASTIC { "add bc root: u = 0" } else { "add bc wall: T = 20 degC" });
        }
        for r in reports[s.first_report..s.first_report + s.report_count] {
            if kind(r.expr) != syntax::N_CALL { continue; }
            let arg0 = arg(r.expr, 0);
            let si = find_surf(name(arg0));
            var on = false;
            for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_FIX && b.surf == si { on = true; } }
            if !on {
                unsupported(tok(arg0), str("the boundary ", name(arg0), " is not fixed in step ", s.name, ", so it has no reaction"),
                            "report the reaction of a boundary with a fixed value");
            }
        }
    }
}

// --- the Goose program -------------------------------------------------------------

let AXES = ["elements::X", "elements::Y", "elements::Z"];
let AXIS_NAMES = ["x", "y", "z"];

fn lit(out: u8[>..]&, x: f64) { codegen::lit(out, x); }

fn lits(out: u8[>..]&, v: f64[3]) {
    lit(out, v[0]);
    out.append(", ");
    lit(out, v[1]);
    out.append(", ");
    lit(out, v[2]);
}

fn quoted(out: u8[>..]&, s: u8[:]) { codegen::quoted(out, s); }

// The unit of the unknown, of its reaction, and of the flux.
fn value_unit() -> const u8[:] { if physics == HEAT { "K" } else { codegen::system.length } }

fn reaction_unit() -> const u8[:] {
    guard physics == HEAT else { return codegen::system.force; }
    if codegen::system.pressure == "Pa" { return "W"; }
    if codegen::system.pressure == "MPa" { return "mW"; }
    return "kW";
}

fn field_unit(what: i64) -> const u8[:] {
    if what == FIELD_VALUE { return value_unit(); }
    if what == FIELD_STRESS { return codegen::system.pressure; }
    if codegen::system.pressure == "Pa" { return "W/m^2"; }
    if codegen::system.pressure == "MPa" { return "mW/mm^2"; }
    return "kW/mm^2";
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

    // A bc on a region that holds no node is a fault of the generator.
    fn region(found: bool, name: u8[:]) {
        if !found {
            write_stderr(str("skein run: internal fault: no node on the boundary '", name, "'. Please report it.\n"));
            exit(3);
        }
    }

    // A point outside the mesh is a fault of the model.
    fn point(found: bool, expr: u8[:]) {
        if !found {
            write_stderr(str("skein run: the point of '", expr, "' is not in the mesh\n"));
            failed = true;
        }
    }

    fn solve(step: u8[:]) {
        let converged = elements::solve();
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of the sum of the sizes of the nodal loads
    // and reactions.
    fn equilibrium(unit: u8[:]) {
        let rx, ry, rz = elements::total_reaction();
        let fx, fy, fz = elements::applied_load();
        let r = sqrt((rx + fx) * (rx + fx) + (ry + fy) * (ry + fy) + (rz + fz) * (rz + fz));
        let limit = 1.0e-6 * elements::load_scale();
        check("equilibrium", r <= limit, quantity(r, unit), quantity(limit, unit));
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let tol = if relative { within * abs(want) } else { within };
        check("expect", abs(got - want) <= tol, str(name, ": ", got, " ", unit), limit);
    }

    // An `expect` that this generator cannot evaluate.
    fn not_run(why: u8[:]) {
        checks.append("  post,expect,not_run,");
        toon(checks, why);
        checks.append(",\"\"\n");
        check_count++;
    }

    fn output(path: u8[:], field: u8[:]) {
        if !elements::write_vtu(path, field) {
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
        format(s, "solve[", solve_count, "]{step,cells,dofs,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        write_stdout(s);
    }

    """;

// Selects the region of boundary `si`.
fn emit_region(out: u8[>..]&, si: i64) {
    let s = surfs[si];
    out.append("    region(");
    if s.hole {
        out.append("elements::hole(");
        lit(out, lo[0]);
        out.append(", ");
        lit(out, lo[1]);
        out.append(", ");
        lit(out, hole_r);
        out.append(")");
    }
    for a, i in s.axes {
        format(out, if i == 0 { "elements::plane(" } else { " && elements::also_plane(" }, AXES[a], ", ");
        lit(out, s.values[i]);
        out.append(")");
    }
    out.append(", ");
    quoted(out, s.name);
    out.append(");\n");
}

fn emit_bc(out: u8[>..]&, b: Bc) {
    let s = surfs[b.surf];
    format(out, "    // bc ", s.name, "\n");
    emit_region(out, b.surf);
    if b.kind == BC_LOAD && physics == HEAT {
        out.append("    elements::flux(");
        lit(out, b.v[0]);
        out.append(");\n");
    } else if b.kind == BC_LOAD {
        out.append("    elements::traction(");
        lits(out, b.v);
        out.append(");\n");
    } else if physics == HEAT {
        out.append("    elements::fix(0, ");
        lit(out, b.v[0]);
        out.append(");\n");
    } else if b.comp >= 0 {
        format(out, "    elements::fix(", AXES[b.comp], ", ");
        lit(out, b.v[0]);
        out.append(");\n");
    } else {
        for j in dim {
            format(out, "    elements::fix(", AXES[j], ", ");
            lit(out, b.v[j]);
            out.append(");\n");
        }
    }
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": static\n");
    format(out, "fn step_", s.name, "() {\n");
    out.append("    elements::clear_bcs();\n");
    // The fixed values first: a load on a node that a later fix holds still
    // counts in the reaction.
    for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_FIX { emit_bc(out, b); } }
    for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_LOAD { emit_bc(out, b); } }
    out.append("    solve(");
    quoted(out, s.name);
    out.append(");\n");
    format(out, "    equilibrium(\"", reaction_unit(), "\");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s.name, r); }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(", ");
        quoted(out, unknown_name);
        out.append(");\n");
    }
    out.append("}\n\n");
}

// `let v, found = <sample at p>;`
fn emit_sample(out: u8[>..]&, var_name: u8[:], what: i64, comp: i64, at: f64[3]) {
    format(out, "        let ", var_name, ", found_", var_name, " = elements::");
    if what == FIELD_VALUE { format(out, "value_at(", comp, ", "); }
    if what == FIELD_FLUX { format(out, "flux_at(", AXES[comp], ", "); }
    if what == FIELD_STRESS { format(out, "stress_at(", AXES[comp / 3], ", ", AXES[comp % 3], ", "); }
    lits(out, at);
    out.append(");\n");
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n");
    let e = r.expr;
    var v: u8[>..] = [];
    if kind(e) == syntax::N_CALL {
        out.append("    ");
        emit_region(out, find_surf(name(arg(e, 0))));
        out.append("        let x, y, z = elements::reaction();\n");
        if physics == HEAT {
            format(v, "quantity(x, \"", reaction_unit(), "\")");
        } else {
            format(v, "vector(x, y, z, \"", reaction_unit(), "\")");
        }
    } else {
        let what, comp = point_field(e, true);
        let at = value(child(e, 1)).v;
        if comp >= 0 {
            emit_sample(out, "v", what, comp, at);
            format(out, "        point(found_v, ");
            quoted(out, r.label);
            out.append(");\n");
            format(v, "quantity(v, \"", field_unit(what), "\")");
        } else {
            for j in 3 { emit_sample(out, str("u", AXIS_NAMES[j]), what, j, at); }
            format(out, "        point(found_ux, ");
            quoted(out, r.label);
            out.append(");\n");
            format(v, "vector(ux, uy, uz, \"", value_unit(), "\")");
        }
    }
    out.append("        report(");
    quoted(out, step);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", v, ");\n    }\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n");
    emit_sample(out, "v", x.what, x.comp, x.at);
    out.append("        point(found_v, ");
    quoted(out, x.name);
    out.append(");\n        expect(");
    quoted(out, x.name);
    out.append(", v, ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"", field_unit(x.what), "\", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}

// Writes the Goose program of the model that `read` read. The program
// imports the elements from `elements_path` and writes its VTU files under
// `results_dir`.
fn generate(elements_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    let sys = codegen::system;
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", sys.length, ", force in ", sys.force,
           ", temperature in K.\n");
    out.append("import std;\nimport os;\nimport elements from ");
    quoted(out, elements_path);
    out.append(";\n");
    out.append(PRELUDE);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n");
    let cell = if dim == 3 { "tetrahedra" } else { "triangles" };
    format(out, "    // mesh on ", body_name, ": P", order, " ", cell, " of at most ", mesh_size, " ", sys.length);
    if geo == GEO_HOLE {
        format(out, ", ", mesh_fine, " ", sys.length, " at the hole\n    elements::mesh_hole(", dim, ", ", order, ", ");
        lits(out, lo);
        out.append(", ");
        lits(out, hi);
        out.append(", ");
        lit(out, hole_r);
        out.append(", ");
        lit(out, mesh_size);
        out.append(", ");
        lit(out, mesh_fine);
        out.append(");\n");
    } else {
        format(out, "\n    elements::mesh_box(", dim, ", ", order, ", ");
        lits(out, lo);
        out.append(", ");
        lits(out, hi);
        for j in 3 {
            out.append(",\n                       elements::cells_for(");
            lit(out, if j < dim { hi[j] - lo[j] } else { 1.0 });
            out.append(", ");
            lit(out, mesh_size);
            out.append(")");
        }
        out.append(");\n");
    }
    format(out, "    // physics ", physics_name, ", material ", material_name, "\n");
    if physics == HEAT {
        out.append("    elements::heat(");
        lit(out, conductivity);
        out.append(", ");
        lit(out, source);
        out.append(", ");
        lit(out, thickness);
        out.append(");\n");
    } else {
        out.append(if plane == PLANE_STRESS { "    elements::elastic_plane_stress(" }
                   else if plane == PLANE_STRAIN { "    elements::elastic_plane_strain(" } else { "    elements::elastic(" });
        lit(out, young);
        out.append(", ");
        lit(out, poisson);
        if plane != SOLID_3D {
            out.append(", ");
            lit(out, thickness);
        }
        out.append(");\n");
    }
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
