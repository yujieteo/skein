// The code generator of the simplex models: a Skein syntax tree to a Goose
// program that calls src/elements/api.goose.
//
// It takes the models that the box slice of src/codegen/ does not: a 2d
// model, a mesh of `tri` or `tet` cells, or a physics with a scalar
// unknown. The physics block decides the equations: the template
// `linear_elastic(b)` or `heat(q_v)`, or the same statements written out
// (`hooke(E, nu, strain(u))` with `div(sigma) = 0`, `-k*grad(T)` or
// `fourier(k, T)` with `div(q) = q_v`), or `div(grad(u)) + c*u = s` for a
// `Complex` scalar unknown u in `harmonic` steps, with `pml` and `bloch`
// (src/elements/wave.goose). A box or a rect, and one less a
// cylinder or a circle on its corner (lo.x, lo.y) (the quarter model of a
// plate with a hole), get the structured meshes of src/elements/mesh.goose.
// Each other shape of the mesher reader (src/mesher/verify.goose) gets the
// mesh of src/mesher/, with its refine lines: the generated program builds
// the shapes of the reader again and meshes them. A 2d elastic model with
// `points` and `method rkpm` solves with src/meshless/ on the nodes of P1
// triangles of the same mesh (RKPM_PRELUDE). Each other construct gives an
// error `unsupported` at its token, with a fix.
//
// The dispatch tables of src/codegen/dispatch.goose add a template or a
// procedure to the box model of src/codegen/codegen.goose. A model of this
// generator differs in the dimension, the geometry and the cells, so
// src/codegen/run.goose sends the whole model here. It uses the units and
// the values of src/codegen/codegen.goose.
import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import mesher from "../mesher";   // the reader of the geometry, and its shapes
import .gen_edge;                  // the edge physics and the Darcy flow
import .gen_solid;                 // plasticity, finite strain and hyperelastic energies
import .gen_contact;               // contact and the rigid face

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
                if kind(u) == syntax::N_UNKNOWN && complex_scalar(child(u, 0)) { return true; }
                if kind(u) == syntax::N_USES && name(u) == "heat" { return true; }
                u = next(u);
            }
        }
        c = next(c);
    }
    return false;
}

// Whether type n is Complex<D> of a scalar D.
fn complex_scalar(n: i64) -> bool {
    return kind(n) == syntax::N_TYPE_APPLY && syntax::text(tok(n)) == "Complex" && kind(first(n)) == syntax::N_TYPE_NAME;
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
        var v = codegen::find_value(name(n));
        if v < 0 { v = property(name(n)); }
        if v >= 0 && codegen::depth < 100 {
            codegen::depth++;
            let r = value(v);
            codegen::depth--;
            return r;
        }
    }
    if k == syntax::N_BINARY && (syntax::text(tok(n)) == "*" || syntax::text(tok(n)) == "/") {
        let a = value(child(n, 0));
        let b = value(child(n, 1));
        if b.n != 1 && (syntax::text(tok(n)) == "/" || a.n != 1) {
            unsupported(tok(n), "a product of two vectors, or a scalar divided by a vector", "multiply or divide a vector by a scalar");
            return a;
        }
        var r = if b.n == 1 { a } else { b };
        let f = if b.n == 1 { b.v[0] } else { a.v[0] };
        for j in 3 { r.v[j] = if syntax::text(tok(n)) == "/" { r.v[j] / f } else { r.v[j] * f }; }
        return r;
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

fn cmul(a: f64[2], b: f64[2]) -> f64[2] { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }

fn cdiv(a: f64[2], b: f64[2]) -> f64[2] {
    let d = b[0] * b[0] + b[1] * b[1];
    if d == 0.0 { return [0.0, 0.0]; }
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
}

// A Complex constant: `im`, and `omega` of the frequency that the
// generator writes, with +, -, *, / and integer powers, over the real
// values of value().
recursive fn cval(n: i64) -> f64[2] {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if k == syntax::N_NAME && name(n) == "im" { return [0.0, 1.0]; }
    if k == syntax::N_NAME && name(n) == "omega" { return [omega_now, 0.0]; }
    if k == syntax::N_NAME && codegen::find_value(name(n)) < 0 && property(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = cval(property(name(n)));
        codegen::depth--;
        return r;
    }
    if k == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = cval(codegen::find_value(name(n)));
        codegen::depth--;
        return r;
    }
    if op == "-" && k == syntax::N_UNARY {
        let a = cval(first(n));
        return [-a[0], -a[1]];
    }
    if op == "+" || op == "-" || op == "*" || op == "/" {
        let a = cval(child(n, 0));
        let b = cval(child(n, 1));
        if op == "+" { return [a[0] + b[0], a[1] + b[1]]; }
        if op == "-" { return [a[0] - b[0], a[1] - b[1]]; }
        return if op == "*" { cmul(a, b) } else { cdiv(a, b) };
    }
    if op == "^" {
        let a = cval(child(n, 0));
        let e = scalar_value(child(n, 1));
        if a[1] == 0.0 && a[0] > 0.0 { return [exp(e * log(a[0])), 0.0]; }
        guard e == floor(e) && abs(e) <= 16.0 else {
            unsupported(tok(n), "a power of a Complex value needs an integer exponent", "write z^2 or z*z");
            return a;
        };
        var r: f64[2] = [1.0, 0.0];
        for i in abs(e) as! i64 { r = cmul(r, a); }
        return if e < 0.0 { cdiv([1.0, 0.0], r) } else { r };
    }
    return [scalar_value(n), 0.0];
}

// --- the model -------------------------------------------------------------------

let GEO_BOX = 0;
let GEO_HOLE = 1;
let GEO_SHAPE = 2;

let ELASTIC = 0;
let HEAT = 1;
let WAVE = 2;       // as in src/elements/assemble.goose

let SOLID_3D = 0;
let PLANE_STRESS = 1;
let PLANE_STRAIN = 2;

// The bc kinds: a fixed value of the unknown or a component of it, a
// traction `sigma*n` or a heat flux `dot(q, n)`, or the unknown of another
// model on the cut of a submodel.
let BC_FIX = 0;
let BC_LOAD = 1;
let BC_SOURCE = 2;  // `bc cut: u = global.u`: the unknown of a `results`
let BC_RADIATE = 3; // `bc r: dot(q, n) = radiation(eps)`: a face of the enclosure of the step
let BC_PRESSURE = 4; // `bc inner: sigma*n = -p*n`: v[0] holds p
let BC_FILM = 5;    // `bc skin: dot(q, n) = h*(T - T_air)`: v holds h and T_air

fn is_load(b: Bc) -> bool { return b.kind == BC_LOAD || b.kind == BC_PRESSURE || b.kind == BC_MOMENT; }

// The pressure p of the traction value -p*n, or NaN.
fn pressure_of(v: i64) -> f64 {
    let neg = kind(v) == syntax::N_UNARY && syntax::text(tok(v)) == "-";
    let m = if neg { first(v) } else { v };
    guard kind(m) == syntax::N_BINARY && syntax::text(tok(m)) == "*" && source_of(child(m, 1)) == "n" else { return 0.0 / 0.0; }
    let p = scalar_value(child(m, 0));
    return if neg { p } else { -p };
}

// A region of the boundary: the union of planes `axis == value` and of
// the surfaces of the shapes `shapes` of the mesher, or the surface of the
// hole.
struct Surf {
    name: u8[..64],
    hole: bool,
    axes: i64[..6],
    values: f64[..6],
    shapes: i64[..6],
    points: f64[..48],  // x, y and z of each point of a Point let: the node nearest to each
    minus: i64,     // the boundary whose nodes a - b removes, or -1
    body: i64,      // the body of a contact that holds the boundary, or -1 (gen_contact.goose)
}

// A refine line of the mesh of a shape: the shape, the size, and whether
// the target is the surface of the shape (src/mesher/refine.goose).
struct Refine {
    g: i64,
    size: f64,
    surface: bool,
}

struct Bc {
    surf: i64,      // index in `surfs`
    kind: i64,
    comp: i64,      // -1 for all components, or the index in `sources` of BC_SOURCE
    v: f64[3],
    tok: i64,
    node: i64,      // the value of a bc of the wave physics: a Complex constant, which can hold omega
}

// A `results` of another model (spec §15.1): the let, the file, its model
// and unit system, the step, and the --set word of `case:`, or "".
struct Source {
    name: u8[..64],
    path: u8[..256],
    model: u8[..64],
    system: i64,
    step: u8[..64],
    set: u8[..256],
}

// A material on a region (`material backfill on lift1`): the shape of the
// region, and the law of the material in the units of the model.
struct Law {
    g: i64,
    item: i64,      // the material item
    young: f64,
    poisson: f64,
    force: f64[3],
    label: u8[..256],
}

// `birth r` or `death r` of a step: the shape of r.
struct Change {
    g: i64,
    birth: bool,
    label: u8[..256],
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
    first_change: i64,
    change_count: i64,
    output: bool,
    procedure: i64,     // STATIC, MODAL, BUCKLING or HARMONIC
    count: i64,         // the modes of a modal or buckling step
    first_freq: i64,    // the frequencies of a harmonic step in `freqs`
    freq_count: i64,
    damping: f64,       // the damping ratio of a harmonic step
    end: f64,           // the end time of a transient step
    dt: f64,            // its time step, or 0 for dt: stable
}

let STATIC = 0;
let MODAL = 1;
let BUCKLING = 2;
let HARMONIC = 3;
let TRANSIENT = 4;  // central_difference (src/elements/explicit.goose)

// What a point field samples (api.goose: value_at, stress_at, flux_at).
let FIELD_VALUE = 0;
let FIELD_STRESS = 1;
let FIELD_FLUX = 2;
let FIELD_EIGEN = 3;    // eigenvalue(comp + 1)
let FIELD_WAVE = 4;     // a part of the wave unknown: comp is the index in PARTS
let FIELD_REACTION = 5; // reaction(r).x, or the heat into the body by integral(dot(q, n), r) on a fixed boundary r
let FIELD_MEAN = 6;     // mean(T, r) or mean(u.x, r) on the boundary r
let FIELD_BODY_MEAN = 7; // mean(sigma.xx, body): comp is 3 i + j
let FIELD_MAX = 8;      // max(sigma.yy, r) or max(-sigma.yy, r) on the boundary r: comp is 3 i + j, + 9 for the minus

struct Expect {
    name: u8[..64],
    what: i64,
    comp: i64,
    at: f64[3],
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
    step: i64,      // the step that checks it at its end, or -1 for an item
    modal: bool,    // an eigenvalue of a modal step: a frequency
    surf: i64,      // the boundary of FIELD_MEAN and FIELD_REACTION, or -1
    scale: f64,     // the factor of the value: -1 for a minus sign or for the heat out of the body
    inverse: bool,  // the value is scale / the field, as -V0/integral(dot(J, n), spot)
    unit: u8[..64], // the unit of the value, for a reaction of a heat physics with no k; else ""
}

var geometry_names: u8[..64][>..];
var geometry_nodes: i64[>..];
var surface_lets: i64[>..];
var point_lets: i64[>..];
var conductivity_name: u8[..64] = [];  // the property of the flux of the heat physics: k, or eps of -eps*grad(phi)
var geo = GEO_BOX;
var lo: f64[3];
var hi: f64[3];
var hole_r = 0.0;
var hole_name: u8[..64];
var body_shape = -1;            // the shape of the mesher of GEO_SHAPE
var refines: Refine[>..];
var body_name: u8[..64];
var surfs: Surf[>..];
var bcs: Bc[>..];
var active: Bc[>..];
var global_bcs: i64[>..];
var reports: Report[>..];
var steps: Step[>..];
var expects: Expect[>..];
var skipped: u8[..256][>..];
var sources: Source[>..];
var cut_source = -1;            // the source of the BC_SOURCE bcs, or -1
var mesh_size = 0.0;
var mesh_fine = 0.0;
var mesh_tok = -1;
var order = 2;
var rkpm_support = 0.0;         // the support of the rkpm method, or 0 for cg
var physics = -1;
var program_node = -1;          // the program that read() reads
var plane = SOLID_3D;
var thickness = 1.0;
var physics_tok = -1;
var physics_name: u8[..64];
var unknown_name: u8[..64];
var flux_name: u8[..64];        // sigma or q
var material_name: u8[..64];
var material_node = -1;         // the material of law 0: of the cells that no other material holds
var material_items: i64[>..];   // each material item
var region_laws: Law[>..];      // the other materials, on regions
var changes: Change[>..];       // the birth and death statements of the steps
var interfaces: i64[>..];       // the shapes whose boundary the mesh keeps
var young = 0.0;
var poisson = -1.0;
var conductivity = 0.0;
var source = 0.0;
var source_node = -1;           // the right side of the eq, or -1 for 0
var body_node = -1;             // the argument b of linear_elastic, or -1
var body_force: f64[3] = [0.0; 3];
var density = 0.0;
var wave_c_node = -1;           // c of eq div(grad(u)) + c*u = s, or -1 for 0
var wave_s_node = -1;           // s, or -1 for 0
var omega_now = 0.0;            // omega of the frequency that the generator writes
var freqs: f64[>..];            // the frequencies of the harmonic steps
var pml_g = -1;                 // the shape of the pml layer, or -1
var pml_tok = -1;
var pml_lo: f64[3] = [0.0; 3];  // the box inside the layer
var pml_hi: f64[3] = [0.0; 3];
var pml_h = 0.0;
var bloch_low = -1;             // the low and high boundaries of a bloch group, or -1
var bloch_high = -1;
var bloch_period: f64[3] = [0.0; 3];
var bloch_theta = 0.0;          // dot(k, period)
var step_now = Step { name: [], first_own: 0, own_count: 0, first_bc: 0, bc_count: 0, first_report: 0, report_count: 0,
                      first_change: 0, change_count: 0, output: false, procedure: STATIC, count: 0, first_freq: 0, freq_count: 0,
                      damping: 0.0, end: 0.0, dt: 0.0 };   // the step that step() reads

fn reset() {
    density = 0.0;
    codegen::reset();
    dim = 3;
    geometry_names.clear();
    geometry_nodes.clear();
    surface_lets.clear();
    point_lets.clear();
    set64(conductivity_name, "k");
    geo = GEO_BOX;
    hole_r = 0.0;
    hole_name = [];
    body_shape = -1;
    refines.clear();
    body_name = [];
    surfs.clear();
    bcs.clear();
    active.clear();
    global_bcs.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    skipped.clear();
    sources.clear();
    cut_source = -1;
    mesh_size = 0.0;
    mesh_fine = 0.0;
    mesh_tok = -1;
    order = 2;
    rkpm_support = 0.0;
    physics = -1;
    plane = SOLID_3D;
    thickness = 1.0;
    physics_tok = -1;
    physics_name = [];
    unknown_name = [];
    flux_name = [];
    material_name = [];
    material_node = -1;
    material_items.clear();
    region_laws.clear();
    changes.clear();
    interfaces.clear();
    young = 0.0;
    poisson = -1.0;
    conductivity = 0.0;
    source = 0.0;
    source_node = -1;
    body_node = -1;
    body_force = [0.0; 3];
    wave_c_node = -1;
    wave_s_node = -1;
    omega_now = 0.0;
    freqs.clear();
    pml_g = -1;
    pml_tok = -1;
    bloch_low = -1;
    bloch_high = -1;
    bloch_theta = 0.0;
    edge_reset();
    solid_reset();
    contact_reset();
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

let SHAPE_NAMES = ["box", "rect", "cylinder", "sphere", "circle", "polygon", "fill", "extrude", "revolve", "loft", "pattern", "pattern_polar"];

fn is_shape(n: i64) -> bool {
    guard n >= 0 && kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME else { return false; }
    for s in SHAPE_NAMES { if name(child(n, 0)) == s { return true; } }
    return false;
}

// Whether `n` is a shape: a shape call, the name of one, or |, & and - of
// shapes.
recursive fn geometric(n: i64) -> bool {
    if kind(n) == syntax::N_NAME { return find_geometry(name(n)) >= 0; }
    if kind(n) == syntax::N_BINARY {
        let op = syntax::text(tok(n));
        return (op == "|" || op == "&" || op == "-") && geometric(child(n, 0)) && geometric(child(n, 1));
    }
    return is_shape(n);
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
// `point(p)`, or a `pattern` of one: a Point let.
recursive fn point_set(v: i64) -> bool {
    return is_call(v, "point") || is_call(v, "pattern") && syntax::child_count(v) == 4 && point_set(arg(v, 0));
}

// Adds the points of the point set v to s.
recursive fn add_points(v: i64, s: Surf&) {
    if is_call(v, "point") {
        let p = value(arg(v, 0));
        if p.n != dim || s.points.len + 3 > s.points.cap {
            unsupported(tok(v), "a Point let has at most 16 points of the dimension of the model", "write point((x, y, z) mm)");
            return;
        }
        for j in 3 { s.points.push(p.v[j]); }
        return;
    }
    let first_point = s.points.len;
    add_points(arg(v, 0), s);
    let step = value(arg(v, 1)).v;
    let count = scalar_value(arg(v, 2)) as! i64;
    let last = s.points.len;
    for k in 1..count {
        for i in first_point / 3..last / 3 {
            if s.points.len + 3 > s.points.cap {
                unsupported(tok(v), "a Point let has at most 16 points", "use fewer points");
                return;
            }
            for j in 3 { s.points.push(s.points[3 * i + j] + k * step[j]); }
        }
    }
}

// `let tip: Point = pattern(point((L, 0 mm, 0 mm)), (0 mm, W/2, 0 mm), 3)`.
fn point_surf(n: i64) {
    var s = Surf { name: [], hole: false, axes: [], values: [], shapes: [], points: [], minus: -1, body: -1 };
    set64(s.name, name(n));
    add_points(child(n, 1), s);
    surfs.push(s);
}

fn value_item(n: i64) {
    let v = child(n, 1);
    var w: u8[..64] = [];
    set64(w, name(n));
    if point_set(v) {
        point_lets.push(n);
    } else if geometric(v) {
        geometry_names.push(w);
        geometry_nodes.push(v);
    } else if is_call(v, "boundary") {
        surface_lets.push(n);
    } else if is_call(v, "results") {
        results_item(n, v);
    } else {
        codegen::value_names.push(w);
        codegen::value_nodes.push(v);
    }
}

// `let global: Results = results(file, step: s, case: p)`. The file is
// relative to the model file.
fn results_item(n: i64, v: i64) {
    var src = Source { name: [], path: [], model: [], system: -1, step: [], set: [] };
    set64(src.name, name(n));
    let file = arg(v, 0);
    guard kind(file) == syntax::N_STRING else { return; };
    let model_path = syntax::files[syntax::tokens[tok(file)].file].path;
    var cut = 0;
    for c, i in model_path { if c == '/' { cut = i + 1; } }
    let quoted_file = name(file);
    codegen::set256(src.path, str(model_path[..cut], quoted_file[1..quoted_file.len - 1]));
    for i in 1..syntax::child_count(v) - 1 {
        let a = child(v, i + 1);
        let named = syntax::nodes[a].aux;
        let x = first(a);
        if named >= 0 && syntax::text(named) == "case" {
            let p = if kind(x) == syntax::N_NAME { codegen::find_value(name(x)) } else { -1 };
            guard p >= 0 else {
                unsupported(tok(x), "the case of results is not a param", "write case: p, with p a param of this model");
                return;
            };
            var word: u8[>..] = [];
            format(word, name(x), "=", source_of(p));
            codegen::set256(src.set, word);
        } else if kind(x) == syntax::N_STRING {
            let t = name(x);
            set64(src.step, t[1..t.len - 1]);
        }
    }
    var text: u8[>..] = [];
    guard read_file(src.path, text) else {
        unsupported(tok(file), str("cannot read '", src.path, "'"), "give the path of a model, relative to this file");
        return;
    };
    let f = syntax::add_file(src.path, text);
    syntax::lex(f);
    let prog = syntax::parse(f);
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_HEADER {
            set64(src.model, name(c));
            for sys, i in codegen::SYSTEM_NAMES { if sys == syntax::text(tok(child(c, 1))) { src.system = i; } }
        }
        c = next(c);
    }
    if src.model.len == 0 || src.system < 0 || src.step.len == 0 {
        unsupported(tok(v), str("the model, the unit system or the step of '", src.path, "' is not in this generator"),
                    "give a model file with a header such as model beam: 3d, SI, and a step");
        return;
    }
    sources.push(src);
}

fn find_source(word: u8[:]) -> i64 {
    for s, i in sources { if s.name == word { return i; } }
    return -1;
}

let BODY_FIX = "write the body as a shape of box, cylinder, sphere, rect, circle, polygon, fill, extrude, revolve, loft, pattern or pattern_polar, with |, & and -";

fn corner(n: i64, i: i64) -> f64[3] {
    let v = value(arg(n, i));
    return v.v;
}

// Whether shape `n` is a box (a rect in 2d), or one less a hole along z
// with its axis on the corner (lo.x, lo.y): the structured meshes of
// mesh_box and mesh_hole.
fn structured(n: i64) -> bool {
    let shape = if dim == 3 { "box" } else { "rect" };
    let round = if dim == 3 { "cylinder" } else { "circle" };
    let cut = kind(n) == syntax::N_BINARY;
    let outer = if cut { geometry_call(child(n, 0)) } else { n };
    guard outer >= 0 && is_call(outer, shape) && syntax::child_count(outer) == 3 else { return false; }
    guard cut else { return true; }
    let c = geometry_call(child(n, 1));
    guard syntax::text(tok(n)) == "-" && c >= 0 && is_call(c, round) else { return false; }
    let a = corner(outer, 0);
    let b = corner(outer, 1);
    let centre = corner(c, 0);
    let tol = 1.0e-6 * max(abs(b[0] - a[0]), abs(b[1] - a[1]));
    guard abs(centre[0] - min(a[0], b[0])) <= tol && abs(centre[1] - min(a[1], b[1])) <= tol else { return false; }
    let r = scalar_value(arg(c, if dim == 3 { 2 } else { 1 }));
    guard r > 0.0 && r < min(abs(b[0] - a[0]), abs(b[1] - a[1])) else { return false; }
    guard dim == 3 else { return true; }
    let axis = corner(c, 1);
    let h = scalar_value(arg(c, 3));
    return axis[0] == 0.0 && axis[1] == 0.0 && axis[2] > 0.0 && centre[2] <= min(a[2], b[2]) + tol &&
        centre[2] + h >= max(a[2], b[2]) - tol;
}

// The body of the physics: a box, or a box less a hole on its corner, for
// the structured meshes, else a shape of the mesher.
fn body(region: i64) {
    let n = if kind(region) == syntax::N_NAME { find_geometry(name(region)) } else { -1 };
    if n < 0 {
        unsupported(tok(region), "the region of the physics is not a solid that a let names", BODY_FIX);
        return;
    }
    set64(body_name, name(region));
    if !structured(n) {
        shape_body(region);
        return;
    }
    let outer = if kind(n) == syntax::N_BINARY { geometry_call(child(n, 0)) } else { n };
    let a = corner(outer, 0);
    let b = corner(outer, 1);
    for j in 3 {
        lo[j] = min(a[j], b[j]);
        hi[j] = max(a[j], b[j]);
        if j < dim && lo[j] == hi[j] { unsupported(tok(outer), "the body has no size", BODY_FIX); }
    }
    guard kind(n) == syntax::N_BINARY else { return; }
    let cut = child(n, 1);
    let c = geometry_call(cut);
    if kind(cut) == syntax::N_NAME { set64(hole_name, name(cut)); }
    hole_r = scalar_value(arg(c, if dim == 3 { 2 } else { 1 }));
    geo = GEO_HOLE;
}

// The body as a shape of the mesher reader: its bounds give lo and hi.
fn shape_body(region: i64) {
    let v = mesher::eval_value(region);
    guard ok() else { return; }
    guard v.n == 0 else {
        unsupported(tok(region), "the region of the physics is not a shape", BODY_FIX);
        return;
    }
    let b = mesher::bounds(v.g);
    let three = b[5] > b[2];
    guard three == (dim == 3) else {
        unsupported(tok(region), str("the body is not a ", dim, "d shape"), BODY_FIX);
        return;
    }
    for j in 3 {
        lo[j] = b[j];
        hi[j] = b[3 + j];
    }
    body_shape = v.g;
    geo = GEO_SHAPE;
}

let SURF_FIX = "write boundary(plate, x == 0 mm), with x, y or z on a face of the body, abs(x - 40 mm) == 15 mm or abs(x) == 15 mm for two faces, an || of these, or boundary(plate, bore) for the hole";

// Adds the planes of condition `cond` to `s`.
recursive fn surf_parts(cond: i64, s: Surf&) -> bool {
    if kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "||" {
        return surf_parts(child(cond, 0), s) && surf_parts(child(cond, 1), s);
    }
    if kind(cond) == syntax::N_NAME && hole_name.len > 0 && name(cond) == hole_name && geo == GEO_HOLE {
        s.hole = true;
        return true;
    }
    if geo == GEO_SHAPE && geometric(cond) && s.shapes.len < s.shapes.cap {
        let v = mesher::eval_value(cond);
        guard v.n == 0 else { return false; }
        s.shapes.push(v.g);
        return true;
    }
    // abs(x - c) == d: the planes x == c - d and x == c + d.
    let is_abs = kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" && is_call(child(cond, 0), "abs") &&
        syntax::child_count(child(cond, 0)) == 2;
    let inner = if is_abs { arg(child(cond, 0), 0) } else { -1 };
    if inner >= 0 && kind(inner) == syntax::N_BINARY && syntax::text(tok(inner)) == "-" && kind(child(inner, 0)) == syntax::N_NAME &&
       axis_of(name(child(inner, 0))) >= 0 && s.axes.len + 2 <= s.axes.cap {
        let axis = axis_of(name(child(inner, 0)));
        let c = scalar_value(child(inner, 1));
        let d = scalar_value(child(cond, 1));
        return face_plane(cond, axis, c - d, s) && face_plane(cond, axis, c + d, s);
    }
    // abs(x) == d: the planes x == -d and x == d.
    if inner >= 0 && kind(inner) == syntax::N_NAME && axis_of(name(inner)) >= 0 && s.axes.len + 2 <= s.axes.cap {
        let d = scalar_value(child(cond, 1));
        return face_plane(cond, axis_of(name(inner)), -d, s) && face_plane(cond, axis_of(name(inner)), d, s);
    }
    if half_space(cond, s) { return ok(); }
    if mantle(cond, s) { return ok(); }
    let is_plane = kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" &&
        kind(child(cond, 0)) == syntax::N_NAME;
    let axis = if is_plane { axis_of(name(child(cond, 0))) } else { -1 };
    if axis < 0 || s.axes.len == s.axes.cap {
        unsupported(tok(cond), "the boundary is not a face of the body", SURF_FIX);
        return false;
    }
    return face_plane(cond, axis, scalar_value(child(cond, 1)), s);
}

// `x > 0 mm` on a box, with 0 mm on its lo face (or `x < a` with a on its
// hi face): each face of the box but that one. It gives false when `cond`
// is no such condition.
fn half_space(cond: i64, s: Surf&) -> bool {
    guard kind(cond) == syntax::N_BINARY && kind(child(cond, 0)) == syntax::N_NAME else { return false; }
    let op = syntax::text(tok(cond));
    let above = op == ">" || op == ">=";
    guard above || op == "<" || op == "<=" else { return false; }
    let axis = axis_of(name(child(cond, 0)));
    guard axis >= 0 else { return false; }
    let v = scalar_value(child(cond, 1));
    let tol = 1.0e-6 * (hi[axis] - lo[axis]);
    let open = if above { lo[axis] } else { hi[axis] };
    if geo != GEO_BOX || abs(v - open) > tol || s.axes.len + 2 * dim - 1 > s.axes.cap {
        unsupported(tok(cond), "a half space is a boundary only on a box, from the plane of a face", SURF_FIX);
        return true;
    }
    for j in dim {
        for end in [lo[j], hi[j]] {
            if j == axis && end == open { continue; }
            s.axes.push(j);
            s.values.push(end);
        }
    }
    return true;
}

// The axis and the centre of `x^2` or `(x - c)^2`, or -1.
fn square(n: i64) -> i64, f64 {
    guard kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == "^" && source_of(child(n, 1)) == "2" else { return -1, 0.0; }
    let b = child(n, 0);
    if kind(b) == syntax::N_NAME { return axis_of(name(b)), 0.0; }
    guard kind(b) == syntax::N_BINARY && syntax::text(tok(b)) == "-" && kind(child(b, 0)) == syntax::N_NAME else { return -1, 0.0; }
    return axis_of(name(child(b, 0))), scalar_value(child(b, 1));
}

// `y^2 + z^2 == R^2` or `sqrt(y^2 + z^2) == R` on a 3d shape: the surface at distance R from the
// line along the third axis, as a cylinder of the mesher past the ends of
// the body. It gives false when `cond` is no such condition.
fn mantle(cond: i64, s: Surf&) -> bool {
    guard kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" else { return false; }
    // sqrt(y^2 + z^2) == R is the same surface.
    let root = is_call(child(cond, 0), "sqrt") && syntax::child_count(child(cond, 0)) == 2;
    let sum = if root { arg(child(cond, 0), 0) } else { child(cond, 0) };
    guard kind(sum) == syntax::N_BINARY && syntax::text(tok(sum)) == "+" else { return false; }
    let a, ca = square(child(sum, 0));
    let b, cb = square(child(sum, 1));
    guard a >= 0 && b >= 0 && a != b else { return false; }
    let r = scalar_value(child(cond, 1));
    let r2 = if root { r * r } else { r };
    if geo != GEO_SHAPE || dim != 3 || r2 <= 0.0 || s.shapes.len == s.shapes.cap {
        unsupported(tok(cond), "a round surface is a boundary only on a 3d shape of the mesher", SURF_FIX);
        return true;
    }
    let axis = 3 - a - b;
    let span = hi[axis] - lo[axis];
    var c: f64[3] = [0.0; 3];
    var d: f64[3] = [0.0; 3];
    c[a] = ca;
    c[b] = cb;
    c[axis] = lo[axis] - span;
    d[axis] = 1.0;
    s.shapes.push(mesher::cylinder(c[0], c[1], c[2], d[0], d[1], d[2], sqrt(r2), 3.0 * span));
    return true;
}

// Adds the plane `axis == v` of condition `cond` to `s`. On a box, it is a
// face. On a shape, it holds a part of the surface, so it is in the bounds.
fn face_plane(cond: i64, axis: i64, v: f64, s: Surf&) -> bool {
    let tol = 1.0e-6 * (hi[axis] - lo[axis]);
    let face = abs(v - lo[axis]) <= tol || abs(v - hi[axis]) <= tol;
    let inside = v >= lo[axis] - tol && v <= hi[axis] + tol;
    if !face && (geo != GEO_SHAPE || !inside) {
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
    let of = if syntax::child_count(v) == 3 && kind(arg(v, 0)) == syntax::N_NAME { name(arg(v, 0)) } else { "" };
    let k = find_body(of);
    if of != body_name && k < 0 {
        unsupported(tok(v), str("the boundary '", name(n), "' is not on the body of the physics"), SURF_FIX);
        return;
    }
    var s = Surf { name: [], hole: false, axes: [], values: [], shapes: [], points: [], minus: -1, body: k };
    set64(s.name, name(n));
    if !surf_parts(arg(v, 1), s) { return; }
    if s.hole && s.axes.len > 0 || s.shapes.len > 0 && s.axes.len > 0 {
        unsupported(tok(v), "a boundary is a surface or planes, not both", "make two boundaries");
        return;
    }
    surfs.push(s);
}

fn material(n: i64) {
    material_items.push(n);
    if material_node < 0 {
        material_node = n;
        set64(material_name, name(n));
    }
}

fn find_material(word: u8[:]) -> i64 {
    for m in material_items { if name(m) == word { return m; } }
    return -1;
}

// The shape of a region of cells: of a material, a birth or a death. The
// mesh keeps its boundary (mesher::add_interface), so a cell is in it or
// out of it as a whole.
fn region_shape(r: i64) -> i64 {
    let v = mesher::eval_value(r);
    guard ok() else { return -1; }
    guard v.n == 0 && geo == GEO_SHAPE && dim == 2 else {
        unsupported(tok(r), "a region of cells needs a 2d body of the mesher: the 3d mesher keeps no interfaces",
                    "use one material and no birth or death in a 3d model");
        return -1;
    };
    interfaces.push(v.g);
    return v.g;
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

// Reads the material of law 0, then the one of each region law.
fn read_material() {
    guard material_node >= 0 else { return; }
    let base = material_node;
    if physics != ELASTIC && region_laws.len > 0 {
        unsupported(tok(base), "a material on a region of a heat or wave physics is not in this generator", "use one material");
    }
    for &w in region_laws {
        material_node = w.item;
        read_law();
        w.young = young;
        w.poisson = poisson;
    }
    material_node = base;
    read_law();
}

// Reads the material `material_node` into young and poisson, or
// conductivity.
fn read_law() {
    var p = first(material_node);
    while p >= 0 {
        let w = name(p);
        // The eq of the wave physics names its properties.
        let known = if physics == ELASTIC { w == "E" || w == "nu" } else { physics == WAVE || physics == EDGE || physics == DARCY || w == conductivity_name || w == "q_v" };
        if !known && w != "rho" && !solid_property(w) {
            unsupported(tok(p), str("the property '", w, "' is not in this generator"),
                        if physics == ELASTIC { "give the material E and nu" } else { "give the material k, and q_v for a source" });
        }
        p = next(p);
    }
    if physics == ELASTIC {
        let e = property("E");
        let nu = property("nu");
        let rho = property("rho");
        if e >= 0 { young = scalar_value(e); }
        if nu >= 0 { poisson = scalar_value(nu); }
        if rho >= 0 { density = scalar_value(rho); }
        if solid_law != SOLID_ENERGY && (young <= 0.0 || poisson <= -1.0 || poisson >= 0.5) {
            unsupported(tok(material_node), "the material needs E > 0 and -1 < nu < 0.5",
                        "give E: Pressure = 210 GPa and nu: Real = 0.3");
        }
    } else if physics == HEAT {
        let k = property(conductivity_name);
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
    let size_name = if kind(n) == syntax::N_POINTS { "spacing" } else { "size" };
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_REFINE && geo == GEO_SHAPE {
            // A boundary by a condition of one surface (a mantle) refines that surface.
            let si = if kind(first(s)) == syntax::N_NAME { find_surf(name(first(s))) } else { -1 };
            let one = si >= 0 && surfs[si].shapes.len == 1 && surfs[si].axes.len == 0;
            var g = if one { surfs[si].shapes[0] } else { -1 };
            var surface = true;
            if !one { g, surface = mesher::refine_target(first(s)); }
            let size = scalar_value(child(s, 1));
            if g >= 0 && size > 0.0 { refines.push(Refine { g: g, size: size, surface: surface }); }
        } else if kind(s) == syntax::N_REFINE {
            let r = first(s);
            let si = if kind(r) == syntax::N_NAME { find_surf(name(r)) } else { -1 };
            if si < 0 || !surfs[si].hole {
                unsupported(tok(r), "this generator refines only the hole", "write refine hole: 0.1 mm");
            } else {
                mesh_fine = scalar_value(child(s, 1));
            }
        } else if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "layers and crack_tip are not in this generator", "remove the line");
        } else if name(s) == size_name {
            mesh_size = scalar_value(v);
        } else if name(s) == "cell" && kind(n) == syntax::N_MESH {
            if kind(v) != syntax::N_NAME || name(v) != cell {
                unsupported(tok(v), str("a ", dim, "d model of this generator has ", cell, " cells"),
                            str("write cell: ", cell));
            }
        } else {
            unsupported(tok(s), str("the setting '", name(s), "' is not in this generator"),
                        if kind(n) == syntax::N_MESH { "use only size, cell and refine" } else { "use only spacing and refine" });
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), str("the ", size_name, " must be > 0"), str("add ", size_name, ": 5 mm")); }
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
    // -eps*grad(phi) or -sigma_e*grad(phi): the heat physics with the
    // property eps or sigma_e in place of k.
    // The parser reads it as (-eps)*grad(phi).
    let product = kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "*" && source_of(child(v, 1)) == str("grad(", u, ")");
    let c = if product { child(v, 0) } else { -1 };
    let scaled = c >= 0 && kind(c) == syntax::N_UNARY && syntax::text(tok(c)) == "-" && kind(first(c)) == syntax::N_NAME;
    guard scaled else { return -1; }
    set64(conductivity_name, name(first(c)));
    return HEAT;
}

fn physics_item(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: this generator has one", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    if family_word.len > 0 {
        edge_physics(n);
        return;
    }
    var eq_count = 0;
    var uses_template = false;
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            body(first(c));
        } else if k == syntax::N_USES {
            if name(c) == "linear_elastic" && syntax::child_count(c) <= 1 {
                physics = ELASTIC;
                uses_template = true;
                set64(flux_name, "sigma");
                if first(c) >= 0 { body_node = first(first(c)); }
            } else if name(c) == "heat" && syntax::child_count(c) <= 1 {
                physics = HEAT;
                uses_template = true;
                set64(flux_name, "q");
                if first(c) >= 0 { source_node = first(first(c)); }
            } else {
                unsupported(tok(c), "this generator has the templates linear_elastic(b) and heat(q_v)", PHYSICS_FIX);
            }
        } else if k == syntax::N_UNKNOWN && unknown_name.len > 0 {
            solid_unknown(c);
        } else if k == syntax::N_UNKNOWN {
            set64(unknown_name, name(c));
            if complex_scalar(child(c, 0)) { physics = WAVE; }
        } else if k == syntax::N_MATASSIGN {
            materials_of(c);
        } else if k == syntax::N_ASSUME {
            assumptions(c);
        } else if k == syntax::N_STATE {
            solid_state(c);
        } else if k == syntax::N_PLASTIC {
            solid_plastic(c);
        } else if k == syntax::N_LET && flux_let(c) < 0 {
            solid_let(c);
        } else if k == syntax::N_LET {
            let p = flux_let(c);
            if physics >= 0 && p != physics {
                unsupported(tok(c), str("the let '", name(c), "' is not a flux of this generator"),
                            "write let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u)), or let q: Vector<HeatFlux> = -k*grad(T)");
            } else {
                physics = p;
                set64(flux_name, name(c));
            }
        } else if k == syntax::N_EQN && solid_dynamic_eq(c) {
            if dynamic_eq && eq_count == 0 { eq_count++; }
        } else if k == syntax::N_EQN && physics == WAVE {
            eq_count++;
            wave_eq(c);
        } else if k == syntax::N_EQN {
            eq_count++;
            let left = child(c, 0);
            let right = child(c, 1);
            let div = is_call(left, "div") && flux_name.len > 0 && (source_of(arg(left, 0)) == flux_name || solid_eq(arg(left, 0)));
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
    solid_finish(n);
    if physics < 0 { unsupported(tok(n), "the physics is neither linear elastic nor heat", PHYSICS_FIX); }
    if physics >= 0 && !uses_template && eq_count != 1 {
        unsupported(tok(n), "a written-out physics needs one eq", "add eq div(sigma) = 0 or eq div(q) = q_v");
    }
    if dim == 2 && physics == ELASTIC && plane == SOLID_3D {
        unsupported(tok(n), "a 2d elastic physics needs plane_stress(t) or plane_strain", "add assume plane_stress(1 mm)");
    }
}

// `material steel, copper on pin`: the first material is law 0. Each
// other one holds its region.
fn materials_of(c: i64) {
    var m = first(c);
    while m >= 0 {
        let item = find_material(name(m));
        if item < 0 {
            unsupported(tok(m), str("'", name(m), "' is not a material of the model"), PHYSICS_FIX);
        } else if m == first(c) {
            material_node = item;
            set64(material_name, name(m));
        } else if first(m) < 0 {
            unsupported(tok(m), "a second material needs a region", str("write ", name(m), " on <region>"));
        } else {
            var w = Law { g: region_shape(first(m)), item: item, young: 0.0, poisson: 0.0, force: [0.0; 3], label: [] };
            codegen::set256(w.label, str(name(m), " on ", source_of(first(m))));
            region_laws.push(w);
        }
        m = next(m);
    }
}

let WAVE_FIX = "write physics sound on duct solves p: Complex<Pressure>, with eq div(grad(p)) + (omega/c)^2*p = 0 in its block";

// eq div(grad(u)) + c*u = s, or div(grad(u)) = s, of a Complex unknown u.
fn wave_eq(n: i64) {
    let left = child(n, 0);
    let lap = str("div(grad(", unknown_name, "))");
    var form = source_of(left) == lap;
    if !form && kind(left) == syntax::N_BINARY && syntax::text(tok(left)) == "+" && source_of(child(left, 0)) == lap {
        let t = child(left, 1);
        if kind(t) == syntax::N_BINARY && syntax::text(tok(t)) == "*" && source_of(child(t, 1)) == unknown_name {
            wave_c_node = child(t, 0);
            form = true;
        }
    }
    if !form { unsupported(tok(n), "this generator solves eq div(grad(u)) + c*u = s for a Complex scalar u", WAVE_FIX); }
    if source_of(child(n, 1)) != "0" { wave_s_node = child(n, 1); }
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
        } else if w == "finite_strain" {
            finite_strain = true;
        } else if w != "small_strain" && w != "isotropic" {
            unsupported(tok(a), str("the assumption '", w, "' is not in this generator"),
                        "use small_strain, finite_strain, isotropic, plane_stress(t) or plane_strain");
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

// The body force b of linear_elastic: a vector, for example rho*(0, 0, -g)
// with rho of the material.
// The body force b of law 0 and of each region law, with the rho of each.
fn body_force_arg() {
    guard body_node >= 0 else { return; }
    let base = material_node;
    for &w in region_laws {
        material_node = w.item;
        w.force = value(body_node).v;
    }
    material_node = base;
    let v = value(body_node);
    if v.n != dim {
        unsupported(tok(body_node), str("the body force b of this ", dim, "d model is a vector of ", dim, " components"),
                    if dim == 2 { "write b: rho*(0, -g)" } else { "write b: rho*(0, 0, -g)" });
        return;
    }
    body_force = v.v;
}

fn method(n: i64) {
    if name(n) == physics_name && syntax::text(syntax::nodes[n].aux) == "rkpm" {
        rkpm_method(n);
        return;
    }
    if name(n) != physics_name || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this generator has only the cg method of the physics", "write method solid: cg { order: 2 }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            order = scalar_value(first(s)) as! i64;
            if order != 1 && order != 2 { unsupported(tok(s), "this generator has order 1 and order 2", "write order: 2"); }
        } else if kind(s) == syntax::N_SETTING && name(s) == "pml" && physics == WAVE {
            pml_setting(first(s));
        } else if (physics == EDGE || physics == DARCY) && edge_method(s) {
            // family and gauge: gen_edge.goose
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this generator",
                        "use only mesh, order, the pml of a wave physics, and the family and gauge: tree of edge elements");
        }
        s = next(s);
    }
    if (physics == EDGE || physics == DARCY) && order != 1 {
        unsupported(tok(n), str("the family ", family_word, " has order 1 in this generator"), "write order: 1");
    }
}

// method solid: rkpm { points: c, basis: linear, support: 2.5,
// integrate: scni, bc: nitsche } of a 2D elastic physics: the meshless
// package on the nodes of P1 triangles of the points item (spec §17.5).
fn rkpm_method(n: i64) {
    order = 1;
    rkpm_support = 2.5;
    var s = first(child(n, 0));
    while s >= 0 {
        let v = first(s);
        let word = if kind(s) == syntax::N_SETTING && kind(v) == syntax::N_NAME { name(v) } else { "" };
        if kind(s) == syntax::N_SETTING && name(s) == "support" {
            rkpm_support = scalar_value(v);
            if rkpm_support <= 1.0 { unsupported(tok(v), "the support must be > 1", "write support: 2.5"); }
        } else if !(kind(s) == syntax::N_SETTING && (name(s) == "points" || (name(s) == "basis" && word == "linear") ||
                    (name(s) == "integrate" && word == "scni") || (name(s) == "bc" && word == "nitsche"))) {
            unsupported(tok(s), "this generator has rkpm with the linear basis, scni and nitsche", "write method solid: rkpm { points: c, basis: linear, support: 2.5, integrate: scni, bc: nitsche }");
        }
        s = next(s);
    }
    if dim != 2 || physics != ELASTIC {
        unsupported(tok(n), "this generator has rkpm on a 2d linear_elastic physics", "use method solid: cg");
    }
}

// The parts of a model that the rkpm lane does not have: it solves static
// steps of one material with loads on the boundary.
fn rkpm_limits(prog: i64) {
    let fix = "use method solid: cg for this model";
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_STUDY { unsupported(tok(c), "a study of an rkpm model is not in this generator", fix); }
        c = next(c);
    }
    for st in steps {
        if st.procedure != STATIC || st.change_count > 0 { unsupported(mesh_tok, "rkpm has only static steps with no birth or death", fix); }
        for r in reports[st.first_report..st.first_report + st.report_count] {
            if max_von_mises(r.expr) { unsupported(tok(r.expr), "rkpm has no max(von_mises(sigma), region)", "report sigma at a point"); }
        }
    }
    if region_laws.len > 0 || body_node >= 0 || cut_source >= 0 {
        unsupported(mesh_tok, "rkpm has one material, no body force and no results of another model", fix);
    }
}

let PML_FIX = "write pml: { region: layer, thickness: 1 m }, with let layer: Surface = world - rect(a, b)";

// pml: { region: layer, thickness: h }. The layer is the body less a box
// (a rect in 2d): the stretch grows with the distance past that box.
fn pml_setting(v: i64) {
    pml_tok = tok(v);
    var region = -1;
    var c = first(v);
    while c >= 0 {
        let field = if syntax::nodes[c].aux >= 0 { syntax::text(syntax::nodes[c].aux) } else { "" };
        if field == "region" { region = first(c); }
        if field == "thickness" { pml_h = scalar_value(first(c)); }
        c = next(c);
    }
    let n = if region >= 0 { geometry_call(region) } else { -1 };
    let inner = if n >= 0 && kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == "-" { geometry_call(child(n, 1)) } else { -1 };
    guard inner >= 0 && is_call(inner, if dim == 3 { "box" } else { "rect" }) && syntax::child_count(inner) == 3 && pml_h > 0.0 else {
        unsupported(pml_tok, "the pml region is not the body less a box, or its thickness is not > 0", PML_FIX);
        return;
    };
    let a = corner(inner, 0);
    let b = corner(inner, 1);
    for j in 3 {
        pml_lo[j] = min(a[j], b[j]);
        pml_hi[j] = max(a[j], b[j]);
    }
    let g = mesher::eval_value(region);
    guard ok() && g.n == 0 else {
        unsupported(tok(region), "the pml region is not a shape", PML_FIX);
        return;
    };
    pml_g = g.g;
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
    if physics == EDGE || physics == DARCY { return edge_bc(n, s); }
    var b = Bc { surf: s, kind: BC_FIX, comp: -1, v: [0.0, 0.0, 0.0], tok: tok(n), node: -1 };
    if rigid_bc(target, v, b) {
        bcs.push(b);
        return bcs.len - 1;
    }
    let t = source_of(target);
    let member = kind(target) == syntax::N_MEMBER && kind(first(target)) == syntax::N_NAME &&
        name(first(target)) == unknown_name;
    let from = if kind(v) == syntax::N_MEMBER && kind(first(v)) == syntax::N_NAME { find_source(name(first(v))) } else { -1 };
    if t == unknown_name && from >= 0 {
        if name(v) != unknown_name {
            unsupported(tok(v), str("the cut takes '", name(v), "' of the other model"), str("write ", unknown_name, " = ", name(first(v)), ".", unknown_name));
        }
        if cut_source >= 0 && cut_source != from {
            unsupported(tok(v), "a second results on a cut: this generator has one", "take the cut values from one results");
        }
        cut_source = from;
        b.kind = BC_SOURCE;
        b.comp = from;
    } else if physics == WAVE && (t == unknown_name || t == str("dot(grad(", unknown_name, "), n)")) {
        b.kind = if t == unknown_name { BC_FIX } else { BC_LOAD };
        b.node = v;
    } else if physics == WAVE {
        unsupported(tok(target), str("the bc target '", t, "' is not in this generator"),
                    str("write bc inlet: ", unknown_name, " = 1 Pa, or bc wall: dot(grad(", unknown_name, "), n) = 0 Pa/m"));
        return -1;
    } else if t == unknown_name && physics == ELASTIC && mentions_xyz(v) {
        b.node = v;
    } else if t == unknown_name {
        let x = value(v);
        if physics == ELASTIC && x.n == 1 && x.v[0] != 0.0 {
            unsupported(tok(v), "u is a vector", "write u = 0, or a vector");
        }
        if physics == HEAT && x.n != 1 { unsupported(tok(v), "T is a scalar", "write T = 20 degC"); }
        b.v = x.v;
    } else if solid_velocity_bc(target, v, b) {
        b.v[0] = 0.0;
    } else if member && physics == ELASTIC && axis_of(name(target)) >= 0 {
        b.comp = axis_of(name(target));
        b.v[0] = scalar_value(v);
    } else if physics == ELASTIC && t == traction_target() && pressure_of(v) == pressure_of(v) {
        b.kind = BC_PRESSURE;
        b.v[0] = pressure_of(v);
    } else if physics == ELASTIC && t == traction_target() {
        let x = value(v);
        if x.n != dim { unsupported(tok(v), "the traction is a vector", BC_FIX_TEXT); }
        b.kind = BC_LOAD;
        b.v = x.v;
    } else if physics == HEAT && (t == str("dot(", flux_name, ", n)") || t == str(flux_name, "*n")) && is_call(v, "radiation") {
        if dim != 2 { unsupported(tok(v), "radiation(eps) runs in a 2d model only", "make the model 2d"); }
        b.kind = BC_RADIATE;
        b.v[0] = scalar_value(arg(v, 0));
    } else if physics == HEAT && (t == str("dot(", flux_name, ", n)") || t == str(flux_name, "*n")) && film(v, b) {
        b.kind = BC_FILM;
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

// `h*(T - T_air)`: h and T_air in b.v. False if v is no such value.
fn film(v: i64, b: Bc&) -> bool {
    guard kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "*" else { return false; }
    let d = child(v, 1);
    guard kind(d) == syntax::N_BINARY && syntax::text(tok(d)) == "-" && source_of(child(d, 0)) == unknown_name else { return false; }
    b.v[0] = scalar_value(child(v, 0));
    b.v[1] = scalar_value(child(d, 1));
    if b.v[0] <= 0.0 { unsupported(tok(v), "the film coefficient h must be > 0", "write dot(q, n) = h*(T - T_air) with h > 0"); }
    return true;
}

fn step(n: i64) {
    let proc = child(n, 0);
    var s = Step { name: [], first_own: bcs.len, own_count: 0, first_bc: 0, bc_count: 0,
                   first_report: reports.len, report_count: 0, first_change: changes.len, change_count: 0, output: false,
                   procedure: STATIC, count: 0, first_freq: 0, freq_count: 0, damping: 0.0, end: 0.0, dt: 0.0 };
    set64(s.name, name(n));
    procedure(proc, s);
    step_now = s;
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
        } else if k == syntax::N_EXPECT {
            expect(c, steps.len);
        } else if k == syntax::N_RELEASE {
            release(c);
        } else if k == syntax::N_BIRTH || k == syntax::N_DEATH {
            var r = first(c);
            while r >= 0 {
                if kind(r) == syntax::N_WHERE {
                    unsupported(tok(r), "birth and death with where are not in this generator", "select the cells with a region");
                } else {
                    var ch = Change { g: region_shape(r), birth: k == syntax::N_BIRTH, label: [] };
                    codegen::set256(ch.label, str(if ch.birth { "birth " } else { "death " }, source_of(r)));
                    changes.push(ch);
                }
                r = next(r);
            }
        } else {
            unsupported(tok(c), "the step statement is not in this generator", "use only bc, output, report, expect, release, birth and death");
        }
        c = next(c);
    }
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    s.change_count = changes.len - s.first_change;
    steps.push(s);
}

let REPORT_FIX = "write report reaction(root), reaction(root).x, mean(u.x, top), max(von_mises(sigma), hole), or report u.x at (1, 0, 0) mm, sigma.xx at (0, 1, 0) mm or T at (0, 0, 0) m";

// The procedure of a step: static with no arguments, or modal or buckling
// with count, on an elastic model, or harmonic on a wave physics.
fn procedure(proc: i64, s: Step&) {
    let p = name(proc);
    if p == "harmonic" || physics == WAVE {
        harmonic(proc, s);
        return;
    }
    if p == "transient" {
        transient(proc, s);
        return;
    }
    guard p == "modal" || p == "buckling" else {
        if p != "static" || first(proc) >= 0 && !solid_ramp(proc, s) {
            unsupported(tok(proc), "this generator has static steps with no arguments, modal and buckling steps with count, and transient steps with central_difference",
                        "write step load: static { ... } or step modes: modal(count: 6) { ... }");
        }
        return;
    };
    s.procedure = if p == "modal" { MODAL } else { BUCKLING };
    var a = first(proc);
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        if named >= 0 && syntax::text(named) == "count" {
            s.count = scalar_value(first(a)) as! i64;
        } else {
            unsupported(tok(a), str("the argument of ", p, " is not in this generator"), str("write ", p, "(count: 6)"));
        }
        a = next(a);
    }
    if s.count < 1 { unsupported(tok(proc), str("a ", p, " step needs count >= 1"), str("write ", p, "(count: 6)")); }
    if physics != ELASTIC { unsupported(tok(proc), str("a ", p, " step needs the linear_elastic template"), PHYSICS_FIX); }
    if s.procedure == MODAL && (density <= 0.0 || region_laws.len > 0) {
        unsupported(tok(proc), "a modal step needs the density rho of one material", "give rho: Density = 7850 kg/m^3 in the material");
    }
    if s.procedure == BUCKLING && steps.len == 0 {
        unsupported(tok(proc), "a buckling step needs a static step before it: its stress", "add step preload: static { ... } before it");
    }
}

let TRANSIENT_FIX = "write step pulse: transient(end: 1 ms, dt: stable, scheme: central_difference) { ... }";

// A transient step of linear_elastic: transient(end:, dt:, scheme:
// central_difference), from rest, with the density of one material.
fn transient(proc: i64, s: Step&) {
    s.procedure = TRANSIENT;
    var scheme = false;
    var a = first(proc);
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        let w = if named >= 0 { syntax::text(named) } else { "" };
        let v = first(a);
        if w == "end" {
            s.end = scalar_value(v);
        } else if w == "dt" {
            s.dt = if kind(v) == syntax::N_NAME && name(v) == "stable" { 0.0 } else { scalar_value(v) };
        } else if w == "scheme" && source_of(v) == "central_difference" {
            scheme = true;
        } else {
            unsupported(tok(a), "this generator has transient(end:, dt:, scheme: central_difference)", TRANSIENT_FIX);
        }
        a = next(a);
    }
    if !scheme { unsupported(tok(proc), "this generator has the explicit scheme central_difference", "add scheme: central_difference"); }
    if s.end <= 0.0 || s.dt < 0.0 { unsupported(tok(proc), "the step needs an end after the start, and dt >= 0", TRANSIENT_FIX); }
    if physics != ELASTIC { unsupported(tok(proc), "a transient step needs the linear_elastic template", PHYSICS_FIX); }
    if solid_law > 0 && !dynamic_eq { unsupported(tok(proc), "a transient step of a nonlinear solid needs the first-order form", DYNAMIC_FIX); }
    if density <= 0.0 || region_laws.len > 0 {
        unsupported(tok(proc), "a transient step needs the density rho of one material", "give rho: Density = 7850 kg/m^3 in the material");
    }
}

// A number with no unit, which takes the unit of its list.
fn bare(n: i64) -> bool {
    let m = if kind(n) == syntax::N_UNARY { first(n) } else { n };
    return kind(m) == syntax::N_NUMBER && syntax::child_of_kind(m, syntax::N_UNIT) < 0;
}

// harmonic(freq: [f1, f2, ...], damping: z) of the wave physics.
fn harmonic(proc: i64, s: Step&) {
    let fix = "write step sweep: harmonic(freq: [100, 200] Hz) { ... }, with a Complex unknown";
    if physics != WAVE && physics != EDGE || name(proc) != "harmonic" {
        unsupported(tok(proc), "this generator has harmonic steps of a Complex unknown, and no other step of it", fix);
        return;
    }
    s.procedure = HARMONIC;
    s.first_freq = freqs.len;
    var a = first(proc);
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        let w = if named >= 0 { syntax::text(named) } else { "" };
        if w == "freq" && kind(first(a)) == syntax::N_LIST {
            let unit = syntax::child_of_kind(first(a), syntax::N_UNIT);
            let f = if unit >= 0 { codegen::to_system(unit) } else { 1.0 };
            var c = first(first(a));
            while c >= 0 {
                if kind(c) != syntax::N_UNIT { freqs.push(scalar_value(c) * (if unit >= 0 && bare(c) { f } else { 1.0 })); }
                c = next(c);
            }
        } else if w == "damping" {
            s.damping = scalar_value(first(a));
        } else {
            unsupported(tok(a), "the argument of harmonic is not in this generator", fix);
        }
        a = next(a);
    }
    s.freq_count = freqs.len - s.first_freq;
    for f in freqs[s.first_freq..] { if !(f > 0.0) { unsupported(tok(proc), "a harmonic step needs frequencies > 0", fix); } }
    if s.freq_count == 0 || s.damping < 0.0 { unsupported(tok(proc), "a harmonic step needs a frequency and damping >= 0", fix); }
}

// The step of an eigenvalue(k): with `own`, the step `step_now` that step()
// reads, else the last modal or buckling step that read() has read.
fn eigen_step(own: bool) -> Step {
    var at = if own { step_now } else { Step { name: [], first_own: 0, own_count: 0, first_bc: 0, bc_count: 0, first_report: 0,
                                               report_count: 0, first_change: 0, change_count: 0, output: false, procedure: STATIC, count: 0,
                                               first_freq: 0, freq_count: 0, damping: 0.0, end: 0.0, dt: 0.0 } };
    if !own { for st in steps { if st.procedure == MODAL || st.procedure == BUCKLING { at = st; } } }
    return at;
}

// eigenvalue(k) of eigen_step(own): k - 1, or -1 if it is not a mode of it.
fn eigen_index(e: i64, own: bool) -> i64 {
    guard is_call(e, "eigenvalue") && syntax::child_count(e) == 2 else { return -1; }
    let count = eigen_step(own).count;
    let k = scalar_value(arg(e, 0)) as! i64;
    if k < 1 || k > count {
        unsupported(tok(e), str("eigenvalue(", k, ") is not a mode of the modal or buckling step"), str("give 1 <= k <= ", count));
    }
    return if k >= 1 && k <= count { k - 1 } else { -1 };
}

fn freq_unit() -> const u8[:] { return if codegen::system.t == 1.0 { "Hz" } else { "kHz" }; }

// `max(von_mises(sigma), r)` of an elastic model, with r a boundary that a
// let names.
fn max_von_mises(e: i64) -> bool {
    guard is_call(e, "max") && syntax::child_count(e) == 3 && physics == ELASTIC else { return false; }
    let f = arg(e, 0);
    let r = arg(e, 1);
    return is_call(f, "von_mises") && syntax::child_count(f) == 2 && source_of(arg(f, 0)) == flux_name &&
           kind(r) == syntax::N_NAME && find_surf(name(r)) >= 0;
}

// The parts of a Complex value that a report or an expect takes at a point.
let PARTS = ["re", "imag", "abs", "arg"];

// The index in PARTS of `re(u at p)` and the others, or -1.
fn wave_part(e: i64) -> i64 {
    guard physics == WAVE && kind(e) == syntax::N_CALL && syntax::child_count(e) == 2 && kind(child(e, 0)) == syntax::N_NAME else { return -1; }
    let f = arg(e, 0);
    guard kind(f) == syntax::N_AT && kind(child(f, 0)) == syntax::N_NAME && name(child(f, 0)) == unknown_name else { return -1; }
    for w, i in PARTS { if name(child(e, 0)) == w { return i; } }
    return -1;
}

fn report(item: i64) {
    let e = first(item);
    if physics == EDGE || physics == DARCY {
        edge_report(item);
        return;
    }
    if physics == WAVE {
        if wave_part(e) < 0 {
            unsupported(tok(e), str("the report '", source_of(e), "' is not in this generator"),
                        str("write report re(", unknown_name, " at (1, 0) m), or imag, abs or arg of it"));
            return;
        }
    } else if is_call(e, "reaction") && syntax::child_count(e) == 2 {
        let r = arg(e, 0);
        if kind(r) != syntax::N_NAME || find_surf(name(r)) < 0 {
            unsupported(tok(r), "reaction needs a boundary that a let names", REPORT_FIX);
            return;
        }
    } else if kind(e) == syntax::N_AT {
        let what, comp = point_field(e, true);
        if what < 0 { return; }
    } else if is_call(e, "eigenvalue") {
        if eigen_index(e, true) < 0 { return; }
    } else if physics == ELASTIC && (kind(e) == syntax::N_MEMBER || is_call(e, "mean")) {
        let what, comp, surf, scale = region_field(e);
        if what < 0 || comp < 0 || surf < 0 || scale != 1.0 {
            unsupported(tok(e), str("the report '", source_of(e), "' is not in this generator"), REPORT_FIX);
            return;
        }
    } else if !max_von_mises(e) {
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

// The boundary of region r: a boundary that a let names, or a - b of two
// of them, with the nodes of b removed (elements::except_plane). Or -1.
fn region_surf(r: i64) -> i64 {
    if kind(r) == syntax::N_NAME { return find_surf(name(r)); }
    guard kind(r) == syntax::N_BINARY && syntax::text(tok(r)) == "-" && kind(child(r, 0)) == syntax::N_NAME && kind(child(r, 1)) == syntax::N_NAME else { return -1; };
    let a = find_surf(name(child(r, 0)));
    let b = find_surf(name(child(r, 1)));
    let word = source_of(r);
    guard a >= 0 && b >= 0 && !surfs[b].hole && word.len <= 64 else { return -1; };
    for s, i in surfs { if s.name == word { return i; } }
    var s = surfs[a];
    set64(s.name, word);
    s.minus = b;
    surfs.push(s);
    return surfs.len - 1;
}

// An `expect` item (step -1), or an `expect` of step `step`.
fn expect(n: i64, step: i64) {
    if physics == EDGE || physics == DARCY {
        edge_expect(n, step);
        return;
    }
    let left = child(n, 0);
    let right = child(n, 1);
    let within = child(n, 2);
    var label: u8[>..] = [];
    format(label, source_of(right), " within ", source_of(within));
    var what = -1;
    var comp = -1;
    if kind(left) == syntax::N_AT {
        what, comp = point_field(left, false);
    } else if is_call(left, "eigenvalue") {
        comp = eigen_index(left, step >= 0);
        what = if comp >= 0 { FIELD_EIGEN } else { -1 };
    } else if wave_part(left) >= 0 {
        comp = wave_part(left);
        what = FIELD_WAVE;
    }
    var surf = -1;
    var scale = 1.0;
    var inverse = false;
    if what < 0 { what, comp, surf, scale, inverse = scaled_field(left); }
    // An expect on a field this generator cannot sample stays, and the
    // report says that it did not run.
    if what < 0 || comp < 0 {
        var s: u8[..256] = [];
        codegen::set256(s, str(name(n), ": ", source_of(left), " is not in this generator"));
        skipped.push(s);
        return;
    }
    let point = if what == FIELD_WAVE { arg(left, 0) } else { left };
    let no_point = what == FIELD_EIGEN || surf >= 0 || what == FIELD_BODY_MEAN;
    var x = Expect { name: [], what: what, comp: comp, at: if no_point { [0.0; 3] } else { value(child(point, 1)).v }, want: scalar_value(right),
                     within: scalar_value(within), relative: kind(within) == syntax::N_PERCENT, label: [], step: step,
                     modal: what == FIELD_EIGEN && eigen_step(step >= 0).procedure == MODAL, surf: surf, scale: scale,
                     inverse: inverse, unit: [] };
    // A reaction of -eps*grad(phi) is no heat: the value takes the unit of
    // the expected value.
    let unit = if kind(right) == syntax::N_NUMBER { syntax::child_of_kind(right, syntax::N_UNIT) } else { -1 };
    if what == FIELD_REACTION && conductivity_name != "k" && unit >= 0 {
        let f = codegen::to_system(unit);
        x.want /= f;
        if !x.relative { x.within /= f; }
        x.scale /= f;
        set64(x.unit, source_of(unit));
    }
    // A tolerance of a temperature is a difference, not a degC value.
    if kind(within) == syntax::N_NUMBER && first(within) >= 0 && source_of(first(within)) == "degC" {
        x.within -= 273.15;
    }
    set64(x.name, name(n));
    codegen::set256(x.label, label);
    expects.push(x);
}

// Whether a bc fixes the unknown on boundary s.
fn fixed_surf(s: i64) -> bool {
    for b in bcs { if b.surf == s && b.kind == BC_FIX { return true; } }
    return false;
}

// region_field(e), or e / c, c / e or c * e with a constant c: the
// field, the component, the boundary, the factor and whether the value is
// the factor / the field.
fn scaled_field(e: i64) -> i64, i64, i64, f64, bool {
    var f = e;
    var sign = 1.0;
    if kind(f) == syntax::N_UNARY && syntax::text(tok(f)) == "-" {
        sign = -1.0;
        f = first(f);
    }
    let op = if kind(f) == syntax::N_BINARY { syntax::text(tok(f)) } else { "" };
    if op == "/" || op == "*" {
        let what, comp, surf, scale = region_field(child(f, 0));
        if what >= 0 && op == "/" { return what, comp, surf, sign * scale / scalar_value(child(f, 1)), false; }
        if what >= 0 { return what, comp, surf, sign * scale * scalar_value(child(f, 1)), false; }
        let w2, c2, s2, k2 = region_field(child(f, 1));
        if w2 >= 0 && op == "/" { return w2, c2, s2, sign * scalar_value(child(f, 0)) / k2, true; }
        if w2 >= 0 { return w2, c2, s2, sign * k2 * scalar_value(child(f, 0)), false; }
        return -1, -1, -1, 1.0, false;
    }
    let what, comp, surf, scale = region_field(e);
    return what, comp, surf, scale, false;
}

// An expect of a value on a boundary: mean(T, r), mean(u.x, r),
// reaction(r).x, reaction(r) of a heat physics, or integral(dot(q, n), r)
// on a boundary r that a bc fixes, with a minus sign or not. It gives the
// field, the component, the boundary and the factor, or -1 for the field.
fn region_field(e: i64) -> i64, i64, i64, f64 {
    var scale = 1.0;
    var f = e;
    if kind(f) == syntax::N_UNARY && syntax::text(tok(f)) == "-" {
        scale = -1.0;
        f = first(f);
    }
    guard physics == ELASTIC || physics == HEAT else { return -1, -1, -1, 1.0; }
    if is_call(f, "mean") && syntax::child_count(f) == 3 && kind(arg(f, 1)) == syntax::N_NAME && name(arg(f, 1)) == body_name && physics == ELASTIC {
        let v = arg(f, 0);
        let m = if kind(v) == syntax::N_MEMBER && kind(first(v)) == syntax::N_NAME && name(first(v)) == flux_name { name(v) } else { "" };
        let i = if m.len == 2 { axis_of(m[..1]) } else { -1 };
        let j = if m.len == 2 { axis_of(m[1..]) } else { -1 };
        guard i >= 0 && j >= 0 else { return -1, -1, -1, 1.0; }
        return FIELD_BODY_MEAN, 3 * i + j, -1, scale;
    }
    if is_call(f, "max") && syntax::child_count(f) == 3 && region_surf(arg(f, 1)) >= 0 && physics == ELASTIC {
        var v = arg(f, 0);
        let minus = kind(v) == syntax::N_UNARY && syntax::text(tok(v)) == "-";
        if minus { v = first(v); }
        let m = if kind(v) == syntax::N_MEMBER && kind(first(v)) == syntax::N_NAME && name(first(v)) == flux_name { name(v) } else { "" };
        let i = if m.len == 2 { axis_of(m[..1]) } else { -1 };
        let j = if m.len == 2 { axis_of(m[1..]) } else { -1 };
        guard i >= 0 && j >= 0 else { return -1, -1, -1, 1.0; }
        return FIELD_MAX, 3 * i + j + (if minus { 9 } else { 0 }), region_surf(arg(f, 1)), scale;
    }
    if is_call(f, "mean") && syntax::child_count(f) == 3 && region_surf(arg(f, 1)) >= 0 {
        let v = arg(f, 0);
        if physics == HEAT && source_of(v) == unknown_name { return FIELD_MEAN, 0, region_surf(arg(f, 1)), scale; }
        let member = kind(v) == syntax::N_MEMBER && kind(first(v)) == syntax::N_NAME && name(first(v)) == unknown_name;
        if physics == ELASTIC && member && axis_of(name(v)) >= 0 { return FIELD_MEAN, axis_of(name(v)), region_surf(arg(f, 1)), scale; }
        return -1, -1, -1, 1.0;
    }
    // The heat out of the body through a fixed boundary is -reaction.
    if physics == HEAT && is_call(f, "integral") && syntax::child_count(f) == 3 && source_of(arg(f, 0)) == str("dot(", flux_name, ", n)") {
        let s = region_surf(arg(f, 1));
        guard s >= 0 && fixed_surf(s) else { return -1, -1, -1, 1.0; }
        return FIELD_REACTION, 0, s, -scale;
    }
    if physics == HEAT && is_call(f, "reaction") && syntax::child_count(f) == 2 && region_surf(arg(f, 0)) >= 0 {
        return FIELD_REACTION, 0, region_surf(arg(f, 0)), scale;
    }
    let call = if kind(f) == syntax::N_MEMBER { first(f) } else { -1 };
    if physics == ELASTIC && call >= 0 && is_call(call, "reaction") && syntax::child_count(call) == 2 && region_surf(arg(call, 0)) >= 0 &&
       axis_of(name(f)) >= 0 {
        return FIELD_REACTION, axis_of(name(f)), region_surf(arg(call, 0)), scale;
    }
    return -1, -1, -1, 1.0;
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// the errors.
fn read(prog: i64) -> bool {
    reset();
    program_node = prog;
    var physics_node = -1;
    var mesh_node = -1;
    var methods: i64[>..] = [];
    var symmetries: i64[>..] = [];
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
        } else if k == syntax::N_MESH || k == syntax::N_POINTS {
            if mesh_node >= 0 { unsupported(tok(c), "a second mesh: this generator has one", "remove the second mesh"); }
            mesh_node = c;
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k == syntax::N_SYMMETRY {
            symmetries.push(c);
        } else if k == syntax::N_CONTACT {
            contact_items.push(c);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT && k != syntax::N_STUDY && k != syntax::N_FN {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this generator"),
                        "remove it: this generator runs one linear static or steady heat model");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if physics_node < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if !ok() { return false; }
    let family = family_of(methods);
    if family == "nedelec" || family == "raviart_thomas" { set64(family_word, family); }
    mesher::read_geometry(prog);
    physics_item(physics_node);
    if !ok() { return false; }
    contact_bodies(prog);
    if !ok() { return false; }
    for s in surface_lets { surface(s); }
    for s in point_lets { point_surf(s); }
    for n in contact_items { contact_item(n); }
    if material_node < 0 { unsupported(at, "the model has no material", "add material steel { E: Pressure = 210 GPa, nu: Real = 0.3 }"); }
    read_material();
    heat_source();
    body_force_arg();
    if mesh_node < 0 { unsupported(at, "the model has no mesh", "add mesh m on plate { size: 5 mm }"); } else { mesh(mesh_node); }
    for m in methods { method(m); }
    for y in symmetries { symmetry(y); }
    if mesh_node >= 0 && (kind(mesh_node) == syntax::N_POINTS) != (rkpm_support > 0.0) {
        unsupported(tok(mesh_node), "points need the rkpm method, and the rkpm method needs points",
                    "write points c on plate { spacing: 2 mm } with method solid: rkpm { points: c }");
    }
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
            expect(c, -1);
        }
        c = next(c);
    }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { bc root: u = 0 ... }"); }
    if rkpm_support > 0.0 { rkpm_limits(prog); }
    if ok() { step_bcs(); }
    if ok() { harmonic_values(); }
    if ok() { supports(); }
    return ok();
}

// symmetry cell: bloch { period: p, low: a, high: b, k: kv } of the wave
// physics: the value on b is the value on a times exp(im*dot(kv, p)).
fn symmetry(n: i64) {
    let fix = "write symmetry cell: bloch { period: (1, 0) m, low: x0, high: x1, k: (3 m^-1, 0 m^-1) }";
    if physics != WAVE || bloch_low >= 0 {
        unsupported(tok(n), "this generator has one bloch group, of a Complex unknown", fix);
        return;
    }
    var kv: f64[3] = [0.0; 3];
    var s = first(syntax::nodes[child(n, 0)].next);
    while s >= 0 {
        let w = name(s);
        let v = first(s);
        if w == "period" { bloch_period = value(v).v; }
        if w == "k" { kv = value(v).v; }
        if w == "low" || w == "high" {
            let at = if kind(v) == syntax::N_NAME { find_surf(name(v)) } else { -1 };
            if at < 0 { unsupported(tok(v), "low and high are boundaries that a let names", fix); }
            if w == "low" { bloch_low = at; } else { bloch_high = at; }
        }
        s = next(s);
    }
    if bloch_low < 0 || bloch_high < 0 { unsupported(tok(n), "bloch needs low and high", fix); }
    bloch_theta = 0.0;
    for j in 3 { bloch_theta += kv[j] * bloch_period[j]; }
}

// Puts `b` in `set`, in place of the bcs that it replaces.
fn merge(set: Bc[>..]&, b: Bc) {
    var kept: Bc[>..] = [];
    var placed = false;
    for c in set {
        if !replaces(b, c) {
            kept.push(c);
        } else if !placed {
            kept.push(b);
            placed = true;
        }
    }
    if !placed { kept.push(b); }
    set.clear();
    set.append(kept);
}

// Whether b replaces c: the same region and target. A rigid face replaces
// each fixed component of u on its region, and a fixed u replaces a rigid
// face.
fn replaces(b: Bc, c: Bc) -> bool {
    guard b.surf == c.surf else { return false; }
    let held = (b.kind == BC_FIX || b.kind == BC_RIGID) && (c.kind == BC_FIX || c.kind == BC_RIGID);
    if held && (b.kind == BC_RIGID || c.kind == BC_RIGID) { return true; }
    return b.kind == c.kind && b.comp == c.comp;
}

// The bcs that hold in each step (spec §12.2), as in src/codegen.
fn step_bcs() {
    var set: Bc[>..] = [];
    for g in global_bcs { merge(set, bcs[g]); }
    for i in 0..steps.len {
        let s = steps[i];
        release_bcs(i, set);
        for b in bcs[s.first_own..s.first_own + s.own_count] { merge(set, b); }
        steps[i].first_bc = active.len;
        steps[i].bc_count = set.len;
        for b in set { active.push(b); }
    }
}

// The Complex values of each frequency of each harmonic step: c, s and the
// bcs at its omega, so that a fault in them stops the model here.
fn harmonic_values() {
    guard physics == WAVE else { return; }
    for s in steps {
        if s.procedure != HARMONIC { continue; }
        for f in freqs[s.first_freq..s.first_freq + s.freq_count] {
            omega_now = 2.0 * PI * f;
            let c = if wave_c_node >= 0 { cval(wave_c_node) } else { [0.0, 0.0] };
            if wave_s_node >= 0 { cval(wave_s_node); }
            for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind != BC_SOURCE { cval(b.node); } }
            if pml_g >= 0 && !(c[0] > 0.0) { unsupported(pml_tok, "the pml needs c > 0 in div(grad(u)) + c*u = s", PML_FIX); }
            if !ok() { return; }
        }
    }
}

// Each step needs a fixed value, and each reaction report a fixed region.
fn supports() {
    // A natural condition alone makes the edge physics and the Darcy flow unique.
    guard physics != EDGE && physics != DARCY else { return; }
    for s in steps {
        var fixes = 0;
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if b.kind == BC_FIX || b.kind == BC_SOURCE { fixes++; } }
        if fixes == 0 && physics != WAVE {
            unsupported(physics_tok, str("step ", s.name, " fixes no value, so the solution is not unique"),
                        if physics == ELASTIC { "add bc root: u = 0" } else { "add bc wall: T = 20 degC" });
        }
        for r in reports[s.first_report..s.first_report + s.report_count] {
            if !is_call(r.expr, "reaction") { continue; }
            let arg0 = arg(r.expr, 0);
            let si = find_surf(name(arg0));
            var on = false;
            for b in active[s.first_bc..s.first_bc + s.bc_count] { if (b.kind == BC_FIX || b.kind == BC_SOURCE) && b.surf == si { on = true; } }
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

// The prefix of the calls of a condition, a solve or a sample: the
// elements, or the rk_ functions of RKPM_PRELUDE.
fn api() -> const u8[:] { if rkpm_support > 0.0 { "rk_" } else { "elements::" } }

fn lits(out: u8[>..]&, v: f64[3]) {
    lit(out, v[0]);
    out.append(", ");
    lit(out, v[1]);
    out.append(", ");
    lit(out, v[2]);
}

fn quoted(out: u8[>..]&, s: u8[:]) { codegen::quoted(out, s); }

// The unit of the unknown, of its reaction, and of the flux.
// The wave physics prints its values with no unit: in the units of the
// model.
fn value_unit() -> const u8[:] { if physics == HEAT && conductivity_name == "k" { "K" } else if physics == HEAT || physics == WAVE { "" } else { codegen::system.length } }

fn reaction_unit() -> const u8[:] {
    guard physics == HEAT else { return codegen::system.force; }
    // The flux of -eps*grad(phi) is no heat: its unit is not in the system.
    if conductivity_name != "k" { return ""; }
    if codegen::system.pressure == "Pa" { return "W"; }
    if codegen::system.pressure == "MPa" { return "mW"; }
    return "kW";
}

fn field_unit(what: i64) -> const u8[:] {
    if what == FIELD_VALUE || what == FIELD_WAVE || what == FIELD_MEAN { return value_unit(); }
    if what == FIELD_STRESS || what == FIELD_BODY_MEAN || what == FIELD_MAX { return codegen::system.pressure; }
    if what == FIELD_REACTION { return reaction_unit(); }
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
    var converged = false;

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

    // A region of cells with no cell is a fault of the model.
    fn cells(found: bool, what: u8[:]) {
        if !found {
            write_stderr(str("skein run: no cell is in the region of '", what, "': give a finer mesh\n"));
            failed = true;
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
        converged = elements::solve();
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // A transient step: central_difference from rest to `end`. The
    // cg_iterations of its solve record are its time steps.
    fn explicit(step: u8[:], end: f64, dt: f64) {
        converged = elements::central_difference(end, dt);
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // A modal (or buckling) step: the count lowest modes, by LOBPCG.
    fn eigen(step: u8[:], count: i64, modal: bool) {
        converged = if modal { elements::modal(count) } else { elements::buckling(count) };
        format(solves, "  ", step, ",", elements::element_count(), ",", elements::dof_count(), ",",
               elements::cg_iterations(), ",", elements::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    fn mode(path: u8[:], field: u8[:], k: i64) {
        if !elements::write_mode(path, field, k) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of elements::load_scale(): the larger of the
    // sum of the sizes of the nodal loads and reactions, and of the terms
    // of K u on the fixed rows.
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

    // The state of a converged step, for a `results` of another model.
    fn state(dir: u8[:], step: u8[:]) {
        if converged && !elements::write_state(dir, step) {
            write_stderr(str("skein run: cannot write the state of step '", step, "' in '", dir, "'\n"));
            failed = true;
        }
    }

    fn print_records() {
        var s: u8[>..] = [];
        format(s, "report[", report_count, "]{step,expr,value}:\n", reports);
        format(s, "verify[", check_count, "]{stage,id,result,value,limit}:\n", checks);
        format(s, "solve[", solve_count, "]{step,cells,dofs,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

// The part of the program of the rkpm method: the meshless package on the
// nodes of the P1 mesh, with the regions of the elements. A condition on a
// region holds the boundary edges of the cloud with both points in it.
let RKPM_PRELUDE = """
    var fixed_values: f64[>..] = [];   // the value of each condition k of the step

    fn rk_value(k: i64, x: f64, y: f64, t: f64) -> f64 { return fixed_values[k]; }
    fn rk_no_source(x: f64, y: f64, t: f64) -> f64 { return 0.0; }

    fn rk_start(support: f64, E: f64, nu: f64, stress: bool, t: f64) {
        var points: f64[>..] = [];
        for n in elements::node_count() { points.append([elements::coords[3 * n], elements::coords[3 * n + 1]]); }
        meshless::set_cloud(points, elements::conn);
        meshless::elastic(E, nu, stress, t);
        if !meshless::rkpm(support) {
            write_stderr("skein run: an rkpm shape function is singular: give a larger support\n");
            exit(1);
        }
    }

    // The boundary edges of the region of the last region call.
    fn rk_select() {
        var mask: u8[>..] = [];
        for n in elements::node_count() { mask.push(if elements::on_region(n) { 1 } else { 0 }); }
        region(meshless::select_points(mask), "the region");
    }

    fn rk_clear_bcs() {
        meshless::clear_bcs();
        fixed_values.clear();
    }

    fn rk_fix(comp: i64, v: f64) {
        rk_select();
        fixed_values.push(v);
        meshless::fix_component(comp, fixed_values.len - 1);
    }

    fn rk_traction(tx: f64, ty: f64, tz: f64) {
        rk_select();
        meshless::load_selected(0, tx);
        meshless::load_selected(1, ty);
    }

    fn rk_reaction() -> f64, f64, f64 {
        rk_select();
        let r = meshless::reaction();
        return r[0], r[1], 0.0;
    }

    fn rk_solve() -> bool { return meshless::solve_static(rk_no_source, rk_value); }
    fn rk_element_count() -> i64 { return meshless::cell_count(); }
    fn rk_dof_count() -> i64 { return meshless::dof_count(); }
    fn rk_cg_iterations() -> i64 { return 0; }
    fn rk_solve_seconds() -> f64 { return meshless::last_seconds; }
    fn rk_value_at(c: i64, x: f64, y: f64, z: f64) -> f64, bool { return meshless::value_at(c, x, y); }
    fn rk_stress_at(i: i64, j: i64, x: f64, y: f64, z: f64) -> f64, bool { return meshless::stress_at(i, j, x, y); }
    fn rk_write_vtu(path: const u8[:], field: const u8[:]) -> bool { return meshless::write_vtu(path, field); }

    fn rk_total_reaction() -> f64, f64, f64 {
        let r, size, f = meshless::balance();
        return r[0], r[1], 0.0;
    }

    fn rk_applied_load() -> f64, f64, f64 {
        let r, size, f = meshless::balance();
        return f[0], f[1], 0.0;
    }

    fn rk_load_scale() -> f64 {
        let r, size, f = meshless::balance();
        return size;
    }

    """;

// The calls of PRELUDE that RKPM_PRELUDE takes.
let RKPM_CALLS = ["solve()", "element_count()", "dof_count()", "cg_iterations()", "solve_seconds()", "total_reaction()",
                  "applied_load()", "load_scale()", "write_vtu("];

// The part of the program of a wave physics.
let WAVE_PRELUDE = """
    // Part `part` of the Complex unknown at a point: re, imag, abs or arg.
    fn wave_at(part: i64, x: f64, y: f64, z: f64) -> f64, bool {
        let re, found_re = elements::value_at(0, x, y, z);
        let im, found_im = elements::value_at(1, x, y, z);
        let found = found_re && found_im;
        if part == 0 { return re, found; }
        if part == 1 { return im, found; }
        return if part == 2 { sqrt(re * re + im * im) } else { atan2(im, re) }, found;
    }

    """;

// The part of the program of a submodel (spec §15.2).
let SUBMODEL_PRELUDE = """
    extern "libc.so.6" {
        fn system(command: const u8&) -> i32;
    }

    // Reads the state of step `step` of the other model in `dir`. If it is
    // not there, `command` runs the other model first.
    fn results(dir: u8[:], step: u8[:], file: u8[:], command: u8[:], scales: f64[3]) {
        if elements::read_state(dir, step, scales[0], scales[1], scales[2]) { return; }
        var c: u8[>..] = [];
        c.append(command);
        c.push(0);
        if system(c[0]) != 0 { write_stderr(str("skein run: the run of '", file, "' for its results failed\n")); }
        if !elements::read_state(dir, step, scales[0], scales[1], scales[2]) {
            write_stderr(str("skein run: no state of step '", step, "' of '", file, "' in '", dir, "'\n"));
            exit(1);
        }
    }

    // Each node of the cut is in the mesh of the other model (checks.md,
    // `submodel_outside`).
    fn cut_outside(outside: i64, name: u8[:]) {
        check("submodel_outside", outside == 0, str(outside, " nodes of ", name), "0 nodes");
    }

    // On each part of the cut, the resultant force of the submodel is the
    // one of the other model, to 2 % of the largest (checks.md,
    // `submodel_cut`).
    fn cut_check(unit: u8[:]) {
        let d, scale = elements::cut_difference();
        check("submodel_cut", d <= 0.02 * scale, quantity(d, unit), quantity(0.02 * scale, unit));
    }

    """;

// The size of the unit of pressure of unit system `sys` in SI.
fn pressure_si(sys: codegen::System) -> f64 { return sys.m / (sys.l * sys.t * sys.t); }

// Reads the state of the cut source, and runs its model first if the state
// is not there. The run of the other model is locked, so the cases of a
// study that need the same state run it once.
fn emit_source(out: u8[>..]&) {
    let src = sources[cut_source];
    var sets: const u8[:][>..] = [];
    if src.set.len > 0 { sets.push(src.set); }
    let dir = codegen::results_path(src.model, sets);
    var inner: u8[>..] = [];
    format(inner, "test -e ", codegen::shell_quoted(str(dir, "/", src.step, ".state_sol")), " || exec ",
           codegen::shell_quoted(str(resource_dir(), "skein")), " run ", codegen::shell_quoted(src.path));
    if src.set.len > 0 { format(inner, " --set ", codegen::shell_quoted(src.set)); }
    inner.append(" 1>&2");
    var command: u8[>..] = [];
    format(command, "mkdir -p ", codegen::shell_quoted(dir), " && flock ", codegen::shell_quoted(str(dir, "/.lock")), " sh -c ",
           codegen::shell_quoted(inner));
    let other = codegen::SYSTEMS[src.system];
    let length = other.l / codegen::system.l;
    format(out, "    // ", src.name, ": the state of step ", src.step, " of ", src.path, "\n    results(");
    quoted(out, dir);
    out.append(", ");
    quoted(out, src.step);
    out.append(", ");
    quoted(out, src.path);
    out.append(", ");
    quoted(out, command);
    out.append(", [");
    lit(out, length);
    out.append(", ");
    lit(out, if physics == HEAT { 1.0 } else { length });
    out.append(", ");
    lit(out, pressure_si(other) / pressure_si(codegen::system));
    out.append("]);\n");
}

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
    for g, i in s.shapes { format(out, if i == 0 { "elements::surface(" } else { " && elements::also_surface(" }, g, ")"); }
    for i in s.points.len / 3 {
        out.append(if i == 0 { "elements::point(" } else { " && elements::also_point(" });
        lits(out, [s.points[3 * i], s.points[3 * i + 1], s.points[3 * i + 2]]);
        out.append(")");
    }
    if s.minus >= 0 {
        let m = surfs[s.minus];
        for a, i in m.axes {
            format(out, " && elements::except_plane(", AXES[a], ", ");
            lit(out, m.values[i]);
            out.append(")");
        }
        for g in m.shapes { format(out, " && elements::except_surface(", g, ")"); }
    }
    if s.body >= 0 { format(out, " && elements::in_body(", s.body, ")"); }
    out.append(", ");
    quoted(out, s.name);
    out.append(");\n");
}

fn emit_bc(out: u8[>..]&, b: Bc) {
    let s = surfs[b.surf];
    format(out, "    // bc ", s.name, "\n");
    emit_region(out, b.surf);
    if b.kind == BC_RIGID || b.kind == BC_MOMENT {
        emit_rigid(out, b);
    } else if b.kind == BC_SOURCE {
        format(out, "    cut_outside(elements::fix_state(), \"", s.name, "\");\n");
    } else if physics == WAVE {
        // The value at the omega of the frequency that the generator writes.
        let z = cval(b.node);
        out.append(if b.kind == BC_LOAD { "    elements::traction(" } else { "    elements::fix(0, " });
        lit(out, z[0]);
        out.append(if b.kind == BC_LOAD { ", " } else { ");\n    elements::fix(1, " });
        lit(out, z[1]);
        out.append(if b.kind == BC_LOAD { ", 0.0);\n" } else { ");\n" });
    } else if b.kind == BC_RADIATE {
        // σ = 5.670374419e-8 kg/(s^3 K^4) in the units of the model.
        out.append("    elements::radiate(");
        lit(out, b.v[0]);
        out.append(", ");
        lit(out, 5.670374419e-8 * codegen::system.t * codegen::system.t * codegen::system.t / codegen::system.m);
        out.append(");\n");
    } else if b.kind == BC_LOAD && s.points.len > 0 {
        out.append("    elements::node_force(");
        lits(out, b.v);
        out.append(");\n");
    } else if b.kind == BC_FILM {
        out.append("    elements::film(");
        lit(out, b.v[0]);
        out.append(", ");
        lit(out, b.v[1]);
        out.append(");\n");
    } else if b.kind == BC_LOAD && physics == HEAT {
        out.append("    elements::flux(");
        lit(out, b.v[0]);
        out.append(");\n");
    } else if b.kind == BC_PRESSURE {
        out.append("    elements::pressure_load(");
        lit(out, b.v[0]);
        out.append(");\n");
    } else if b.kind == BC_LOAD {
        format(out, "    ", api(), "traction(");
        lits(out, b.v);
        out.append(");\n");
    } else if physics == HEAT {
        out.append("    elements::fix(0, ");
        lit(out, b.v[0]);
        out.append(");\n");
    } else if b.node >= 0 {
        // A value of the point p of the node.
        for j in dim {
            format(out, "    elements::fix_with(", AXES[j], ") { p => ");
            emit_space(out, b.node, j);
            out.append(" };\n");
        }
    } else if b.comp >= 0 {
        format(out, "    ", api(), "fix(", AXES[b.comp], ", ");
        lit(out, b.v[0]);
        out.append(");\n");
    } else {
        for j in dim {
            format(out, "    ", api(), "fix(", AXES[j], ", ");
            lit(out, b.v[j]);
            out.append(");\n");
        }
    }
}

// A modal or buckling step: the bcs that fix, then LOBPCG. It writes a
// file for each mode. The stages of a study share no mode.
fn emit_eigen_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    let modal = s.procedure == MODAL;
    format(out, "// step ", s.name, ": ", if modal { "modal" } else { "buckling" }, "\nfn step_", s.name, "() {\n");
    out.append("    elements::clear_bcs();\n");
    for ch in changes[s.first_change..s.first_change + s.change_count] { emit_change(out, ch); }
    for b in active[s.first_bc..s.first_bc + s.bc_count] { if !is_load(b) { emit_bc(out, b); } }
    out.append("    eigen(");
    quoted(out, s.name);
    format(out, ", ", s.count, ", ", if modal { "true" } else { "false" }, ");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s, r); }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        for k in 1..s.count + 1 {
            var path: u8[>..] = [];
            format(path, results_dir, "/", s.name, "_mode", k, ".vtu");
            out.append("    mode(");
            quoted(out, path);
            out.append(", ");
            quoted(out, unknown_name);
            format(out, ", ", k, ");\n");
        }
    }
    out.append("}\n\n");
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:], share: codegen::Share) {
    let s = steps[si];
    if physics == EDGE || physics == DARCY {
        emit_edge_step(out, si, results_dir);
        return;
    }
    if s.procedure == HARMONIC {
        emit_harmonic_step(out, si, results_dir);
        return;
    }
    if s.procedure == MODAL || s.procedure == BUCKLING {
        emit_eigen_step(out, si, results_dir);
        return;
    }
    format(out, "// step ", s.name, ": ", if s.procedure == TRANSIENT { "transient, central_difference" } else { "static" }, "\n");
    format(out, "fn step_", s.name, "() {\n");
    // After the clear: clear_bcs keeps the solution of the last step,
    // which a birth reads.
    format(out, "    ", api(), "clear_bcs();\n");
    for ch in changes[s.first_change..s.first_change + s.change_count] { emit_change(out, ch); }
    if share.stages >= 2 && share.reads {
        codegen::emit_share(out, "", "elements::read_solution", share, s.name);
    } else {
        // The fixed values first: a load on a node that a later fix holds
        // still counts in the reaction.
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if !is_load(b) { emit_bc(out, b); } }
        for b in active[s.first_bc..s.first_bc + s.bc_count] { if is_load(b) { emit_bc(out, b); } }
        if s.procedure == TRANSIENT && solid_law == 0 {
            out.append("    explicit(");
            quoted(out, s.name);
            out.append(", ");
            lit(out, s.end);
            out.append(", ");
            lit(out, s.dt);
            out.append(");\n");
        } else if solid_law > 0 {
            emit_solid_solve(out, s);
        } else if has_ties() {
            out.append("    nonlinear(");
            quoted(out, s.name);
            format(out, ", elements::solve_contact(", max(1, s.count), "));\n");
        } else {
            out.append("    solve(");
            quoted(out, s.name);
            out.append(");\n");
        }
        if share.stages >= 2 { codegen::emit_share(out, "converged && ", "elements::write_solution", share, s.name); }
        // The state of the step, for a `results` of another model.
        if rkpm_support <= 0.0 {
            out.append("    state(");
            quoted(out, results_dir);
            out.append(", ");
            quoted(out, s.name);
            out.append(");\n");
        }
    }
    // A transient step has the inertia forces too: no equilibrium check.
    if s.procedure == STATIC { format(out, "    equilibrium(\"", reaction_unit(), "\");\n"); }
    if physics == ELASTIC {
        for b in active[s.first_bc..s.first_bc + s.bc_count] {
            if b.kind != BC_SOURCE { continue; }
            emit_region(out, b.surf);
            format(out, "    cut_check(\"", reaction_unit(), "\");\n");
        }
    }
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s, r); }
    for x in expects { if x.step == si { emit_expect(out, x); } }
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

// A harmonic step: a function for each frequency, which calls the next.
fn emit_harmonic_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": harmonic, ", s.freq_count, " frequencies\nfn step_", s.name, "() { step_", s.name, "_1(); }\n\n");
    for i in s.freq_count { emit_frequency(out, si, i, results_dir); }
}

// Frequency i of harmonic step si: the terms, the pml and the bcs at its
// omega, the solve and the reports. The last frequency checks the expects
// of the step and writes the output.
fn emit_frequency(out: u8[>..]&, si: i64, i: i64, results_dir: u8[:]) {
    let s = steps[si];
    let f = freqs[s.first_freq + i];
    omega_now = 2.0 * PI * f;
    let c = if wave_c_node >= 0 { cval(wave_c_node) } else { [0.0, 0.0] };
    let src = if wave_s_node >= 0 { cval(wave_s_node) } else { [0.0, 0.0] };
    var label: u8[>..] = [];
    format(label, s.name, " at ", f, " ", freq_unit());
    format(out, "// step ", s.name, ", frequency ", i + 1, " of ", s.freq_count, "\nfn step_", s.name, "_", i + 1, "() {\n");
    out.append("    elements::clear_bcs();\n    elements::wave_terms(");
    for v in [c[0], c[1], src[0], src[1]] {
        lit(out, v);
        out.append(", ");
    }
    lit(out, s.damping);
    out.append(");\n");
    if pml_g >= 0 {
        format(out, "    elements::pml(", pml_g, ", ");
        lits(out, pml_lo);
        out.append(", ");
        lits(out, pml_hi);
        out.append(", ");
        lit(out, pml_h);
        out.append(");\n");
    }
    for b in active[s.first_bc..s.first_bc + s.bc_count] { emit_bc(out, b); }
    out.append("    solve(");
    quoted(out, label);
    out.append(");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_wave_report(out, label, r); }
    let last = i == s.freq_count - 1;
    if last { for x in expects { if x.step == si { emit_expect(out, x); } } }
    if last && s.output {
        out.append("    output(");
        quoted(out, str(results_dir, "/", s.name, ".vtu"));
        out.append(", ");
        quoted(out, unknown_name);
        out.append(");\n");
    }
    if !last { format(out, "    step_", s.name, "_", i + 2, "();\n"); }
    out.append("}\n\n");
}

fn emit_wave_report(out: u8[>..]&, step_label: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n");
    emit_sample(out, "v", FIELD_WAVE, wave_part(r.expr), value(child(arg(r.expr, 0), 1)).v);
    out.append("        point(found_v, ");
    quoted(out, r.label);
    out.append(");\n        report(");
    quoted(out, step_label);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", quantity(v, \"", if wave_part(r.expr) == 3 { "rad" } else { value_unit() }, "\"));\n    }\n");
}

// The ties of the bloch group: the high boundary to the low one.
fn emit_bloch(out: u8[>..]&) {
    out.append("    // symmetry: bloch\n");
    emit_region(out, bloch_low);
    out.append("    elements::bloch_low();\n");
    emit_region(out, bloch_high);
    out.append("    if !elements::bloch(");
    lits(out, bloch_period);
    out.append(", ");
    lit(out, bloch_theta);
    out.append(") {\n        write_stderr(\"skein run: the period of the bloch group does not map low onto high\\n\");\n");
    out.append("        failed = true;\n    }\n");
}

fn emit_change(out: u8[>..]&, ch: Change) {
    format(out, "    cells(elements::", if ch.birth { "birth(" } else { "kill(" }, ch.g, "), ");
    quoted(out, ch.label);
    out.append(");\n");
}

// `let v, found = <sample at p>;`
fn emit_sample(out: u8[>..]&, var_name: u8[:], what: i64, comp: i64, at: f64[3]) {
    if what == FIELD_WAVE {
        format(out, "        let ", var_name, ", found_", var_name, " = wave_at(", comp, ", ");
        lits(out, at);
        out.append(");\n");
        return;
    }
    if what == FIELD_REACTION || what == FIELD_MEAN || what == FIELD_MAX {
        out.append("    ");
        emit_region(out, at[0] as! i64);
        if what == FIELD_MAX {
            let c = comp % 9;
            format(out, "        let ", var_name, " = elements::max_stress(", AXES[c / 3], ", ", AXES[c % 3], ", ", if comp >= 9 { "-1.0" } else { "1.0" }, ");\n");
        } else if what == FIELD_MEAN {
            format(out, "        let ", var_name, " = elements::mean_value(", AXES[comp], ");\n");
        } else {
            format(out, "        let r_x, r_y, r_z = ", api(), "reaction();\n        let ", var_name, " = r_", AXIS_NAMES[comp], ";\n");
        }
        format(out, "        let found_", var_name, " = ", var_name, " == ", var_name, ";\n");
        return;
    }
    format(out, "        let ", var_name, ", found_", var_name, " = ", api());
    if what == FIELD_VALUE { format(out, "value_at(", comp, ", "); }
    if what == FIELD_FLUX { format(out, "flux_at(", AXES[comp], ", "); }
    if what == FIELD_STRESS { format(out, "stress_at(", AXES[comp / 3], ", ", AXES[comp % 3], ", "); }
    if what == FIELD_EIGEN {
        format(out, "eigenvalue(", comp + 1, "), true;\n");
        return;
    }
    lits(out, at);
    out.append(");\n");
}

fn emit_report(out: u8[>..]&, step: Step, r: Report) {
    format(out, "    // report ", r.label, "\n    {\n");
    let e = r.expr;
    var v: u8[>..] = [];
    if is_call(e, "eigenvalue") {
        format(out, "        let x = elements::eigenvalue(", scalar_value(arg(e, 0)) as! i64, ");\n");
        format(v, "quantity(x, \"", if step.procedure == MODAL { freq_unit() } else { "" }, "\")");
    } else if max_von_mises(e) {
        out.append("    ");
        emit_region(out, find_surf(name(arg(e, 1))));
        out.append("        let x = elements::max_von_mises();\n");
        format(v, "quantity(x, \"", codegen::system.pressure, "\")");
    } else if physics == ELASTIC && (kind(e) == syntax::N_MEMBER || is_call(e, "mean")) {
        let what, comp, surf, scale = region_field(e);
        emit_sample(out, "x", what, comp, [surf * 1.0, 0.0, 0.0]);
        format(v, "quantity(x, \"", field_unit(what), "\")");
    } else if kind(e) == syntax::N_CALL {
        out.append("    ");
        emit_region(out, find_surf(name(arg(e, 0))));
        format(out, "        let x, y, z = ", api(), "reaction();\n");
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
    quoted(out, step.name);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", v, ");\n    }\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n");
    if x.what == FIELD_BODY_MEAN {
        out.append("        let v = ");
        lit(out, x.scale);
        format(out, " * elements::mean_stress(", x.comp / 3, ", ", x.comp % 3, ");\n        let found_v = true;\n");
    } else if x.what == FIELD_MAX {
        emit_sample(out, "m", x.what, x.comp, [x.surf * 1.0, 0.0, 0.0]);
        out.append("        let v = ");
        lit(out, x.scale);
        out.append(" * m;\n        let found_v = found_m;\n");
    } else if x.surf >= 0 {
        out.append("    ");
        emit_region(out, x.surf);
        if x.what == FIELD_MEAN {
            format(out, "        let v = ");
            lit(out, x.scale);
            format(out, " * elements::mean_value(", x.comp, ");\n");
        } else {
            format(out, "        let rx, ry, rz = ", api(), "reaction();\n        let r: f64[3] = [rx, ry, rz];\n        let v = ");
            lit(out, x.scale);
            format(out, if x.inverse { " / r[" } else { " * r[" }, x.comp, "];\n");
        }
        out.append("        let found_v = true;\n");
    } else {
        emit_sample(out, "v", x.what, x.comp, x.at);
    }
    out.append("        point(found_v, ");
    quoted(out, x.name);
    out.append(");\n        expect(");
    quoted(out, x.name);
    out.append(", v, ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    let unit = if x.unit.len > 0 { x.unit[..] } else if x.what == FIELD_WAVE && x.comp == 3 { "rad" } else if x.what == FIELD_REACTION { reaction_unit() } else if x.what != FIELD_EIGEN { field_unit(x.what) } else if x.modal { freq_unit() } else { "" };
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"", unit, "\", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}

// The dofs of the mesh of the model that `read` read: for the hole, the
// rows of cells that grow from `mesh_fine` by 1.25 (elements::grade_rows).
// A shape is meshed: a P2 node is at each point and at each edge, and a
// mesh has about 7 edges for each point in 3d, 3 in 2d.
fn dof_estimate() -> i64 {
    if geo == GEO_SHAPE {
        mesher::clear_refine();
        for r in refines { mesher::add_refine(r.g, r.size, r.surface); }
        mesher::clear_interfaces();
        for g in interfaces { mesher::add_interface(g); }
        let done = if dim == 3 { mesher::mesh_3d(body_shape, mesh_size) } else { mesher::mesh_2d(body_shape, mesh_size) };
        let points = if done { mesher::point_count() } else { 0 };
        if physics == EDGE || physics == DARCY { return edge_dofs(points); }
        let nodes = if order == 2 { points * (if dim == 3 { 8 } else { 4 }) } else { points };
        return nodes * per_node();
    }
    var cells: i64[3] = [1, 1, 1];
    if geo == GEO_HOLE {
        let w = hi[0] - lo[0];
        let h = hi[1] - lo[1];
        cells[0] = max(ceil(0.5 * PI * hole_r / mesh_fine) as! i64, ceil((w + h) / mesh_size) as! i64);
        cells[1] = ceil((sqrt(w * w + h * h) - hole_r) / mesh_size + log(mesh_size / mesh_fine) / log(1.25)) as! i64 + 1;
    } else {
        for j in 2 { cells[j] = max(ceil((hi[j] - lo[j]) / mesh_size - 1.0e-9) as! i64, 1); }
    }
    if dim == 3 { cells[2] = max(ceil((hi[2] - lo[2]) / mesh_size - 1.0e-9) as! i64, 1); }
    var nodes = 1;
    for j in dim { nodes *= order * cells[j] + 1; }
    if physics == EDGE || physics == DARCY { return edge_dofs(nodes); }
    return nodes * per_node();
}

// The dofs of a node.
fn per_node() -> i64 { return if physics == HEAT { 1 } else if physics == WAVE { 2 } else { dim }; }

// The shapes of the mesher reader and the refine lines, again in the
// program: mesh_shape meshes them, and a boundary can be the surface of one.
fn emit_shapes(out: u8[>..]&) {
    out.append("    // the shapes of the geometry, in the sequence of the mesher reader\n");
    out.append("    mesher::clear_shapes();\n    mesher::clear_refine();\n");
    if mesher::polygon_points.len > 0 {
        out.append("    mesher::polygon_points.append([");
        for v, i in mesher::polygon_points {
            if i > 0 { out.append(", "); }
            lit(out, v);
        }
        out.append("]);\n");
    }
    for g in mesher::shapes {
        format(out, "    mesher::add_shape(", g.kind, ", ", g.a, ", ", g.b, ", [");
        for v, i in g.p {
            if i > 0 { out.append(", "); }
            lit(out, v);
        }
        out.append("]);\n");
    }
    for r in refines {
        format(out, "    mesher::add_refine(", r.g, ", ");
        lit(out, r.size);
        format(out, ", ", if r.surface { "true" } else { "false" }, ");\n");
    }
    for g in interfaces { format(out, "    mesher::add_interface(", g, ");\n"); }
}

// Writes the Goose program of the model that `read` read. The program
// imports the elements from `elements_path` and writes its VTU files under
// `results_dir`. In a study, it writes or reads the stages of `share`.
fn generate(elements_path: u8[:], results_dir: u8[:], share: codegen::Share, out: u8[>..]&) {
    let sys = codegen::system;
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", sys.length, ", force in ", sys.force,
           ", temperature in K.\n");
    let src = elements_path[..elements_path.len - "elements".len];
    out.append("import std;\nimport os;\nimport elements from ");
    quoted(out, elements_path);
    out.append(";\nimport mesher from ");
    quoted(out, str(src, "mesher"));
    out.append(";   // the elements import it\nimport runtime from ");
    quoted(out, str(src, "runtime"));
    out.append(";   // the mesher imports it\nimport solvers from ");
    quoted(out, str(src, "solvers"));
    out.append(";   // the runtime imports it\n");
    if rkpm_support > 0.0 {
        out.append("import meshless from ");
        quoted(out, str(src, "meshless"));
        out.append(";\n");
        var prelude: u8[>..] = [];
        prelude.append(PRELUDE);
        for c in RKPM_CALLS {
            var next: u8[>..] = [];
            format_replaced(next, prelude, str("elements::", c), str("rk_", c));
            prelude.clear();
            prelude.append(next);
        }
        out.append(prelude);
        out.append(RKPM_PRELUDE);
    } else {
        out.append(PRELUDE);
    }
    if cut_source >= 0 { out.append(SUBMODEL_PRELUDE); }
    if physics == WAVE { out.append(WAVE_PRELUDE); }
    if physics == EDGE || physics == DARCY { out.append(EDGE_PRELUDE); }
    if solid_law > 0 || has_ties() { emit_energy(out); }
    for si in 0..steps.len { emit_step(out, si, results_dir, share); }
    out.append("fn main() {\n");
    let cell = if dim == 3 { "tetrahedra" } else { "triangles" };
    if geo == GEO_SHAPE || pml_g >= 0 { emit_shapes(out); }
    if share.stages >= 1 && share.reads {
        codegen::emit_share(out, "", "elements::read_mesh", share, "");
    } else {
        format(out, "    // mesh on ", body_name, ": P", order, " ", cell, " of at most ", mesh_size, " ", sys.length);
        if body_shapes.len > 0 {
            emit_bodies(out);
        } else if geo == GEO_SHAPE {
            format(out, ", by the mesher\n    if !elements::mesh_shape(", dim, ", ", order, ", ", body_shape, ", ");
            lit(out, mesh_size);
            out.append(") {\n        write_stderr(\"skein run: the mesher did not finish in its budget of points: give a larger size\\n\");\n");
            out.append("        exit(1);\n    }\n");
        } else if geo == GEO_HOLE {
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
        if share.stages >= 1 { codegen::emit_share(out, "", "elements::write_mesh", share, ""); }
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
    } else if physics == EDGE || physics == DARCY {
        out.append("    // each step sets the coefficients of the physics\n");
    } else if physics == WAVE {
        out.append("    elements::wave(");
        lit(out, thickness);
        out.append(");\n");
        if bloch_low >= 0 { emit_bloch(out); }
    } else if solid_law == SOLID_ENERGY {
        emit_solid_law(out);
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
        if rkpm_support > 0.0 {
            format(out, "    // method rkpm: support ", rkpm_support, " times the spacing, on the nodes of the mesh\n    rk_start(");
            lit(out, rkpm_support);
            out.append(", ");
            lit(out, young);
            out.append(", ");
            lit(out, poisson);
            format(out, ", ", if plane == PLANE_STRESS { "true" } else { "false" }, ", ");
            lit(out, thickness);
            out.append(");\n");
        }
        emit_solid_law(out);
        if body_node >= 0 {
            out.append("    elements::body_force(");
            lits(out, body_force);
            out.append(");\n");
        }
        if density > 0.0 {
            out.append("    elements::set_density(");
            lit(out, density);
            out.append(");\n");
        }
        for w in region_laws {
            format(out, "    cells(elements::material_in(", w.g, ", ");
            lit(out, w.young);
            out.append(", ");
            lit(out, w.poisson);
            out.append(", ");
            lits(out, w.force);
            out.append("), ");
            quoted(out, w.label);
            out.append(");\n");
        }
    }
    if cut_source >= 0 { emit_source(out); }
    emit_contacts(out);
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    for x in probes { if x.step < 0 { emit_edge_expect(out, x); } }
    for s in skipped {
        out.append("    not_run(");
        quoted(out, s);
        out.append(");\n");
    }
    out.append("    print_records();\n");
    out.append("    if failed { exit(1); }\n");
    out.append("}\n");
}
