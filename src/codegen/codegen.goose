// The code generator of the first slice (docs/implementation.md §6): a
// Skein syntax tree to a Goose program that calls src/runtime/api.goose.
//
// The slice is one linear static solid: a CSG box, a quarter of a tube or
// a union of boxes, a structured mesh of 27-node hexahedra, a linear
// elastic material on each region, with a thermal strain of constant
// temperatures (spec §10.4), `u = 0`, `u.x = 0` and `sigma*n` boundary
// conditions on the faces of the solid, `u.x = 0` at a Point, and `static`
// steps. Each other construct gives an error `unsupported` at its token,
// with a fix.
//
// The generator evaluates each value at compile time, in the unit system
// of the header, because the runtime has no units. The generated program
// prints the reports and the verification records in TOON
// (command_line.md §6.3, §6.4).
import std;
import math;
namespace codegen;
import syntax from "../syntax";
import .dispatch;

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
    let lo, hi = source_span(n);
    return syntax::source[lo..hi];
}

// The offsets in `syntax::source` of the text of `source_of(n)`.
fn source_span(n: i64) -> i64, i64 {
    var lo = syntax::tokens[tok(n)].start;
    var hi = lo;
    span(n, lo, hi);
    // Parentheses make no node, so add the ones that the span cuts.
    let src = syntax::source;
    var depth = 0;
    var least = 0;
    for c in src[lo..hi] {
        if c == '(' { depth++; }
        if c == ')' { depth--; }
        if depth < least { least = depth; }
    }
    var closes = -least;
    var opens = depth - least;
    while closes > 0 {
        var e = lo - 1;
        while e >= 0 && src[e] == ' ' { e--; }
        if e < 0 || src[e] != '(' { break; }
        lo = e;
        closes--;
    }
    while opens > 0 {
        var e = hi;
        while e < src.len && src[e] == ' ' { e++; }
        if e >= src.len || src[e] != ')' { break; }
        hi = e + 1;
        opens--;
    }
    return lo, hi;
}

recursive fn span(n: i64, lo: i64&, hi: i64&) {
    let t = syntax::tokens[tok(n)];
    if t.start < lo { lo = t.start; }
    if t.start + t.len > hi { hi = t.start + t.len; }
    if syntax::nodes[n].aux >= 0 && kind(n) != syntax::N_DIM_POW && kind(n) != syntax::N_UNIT_POW {
        let a = syntax::tokens[syntax::nodes[n].aux];
        if a.start + a.len > hi { hi = a.start + a.len; }
    }
    // The exponent of a power makes no node: a number, or `-` and one, after the `^`.
    if kind(n) == syntax::N_DIM_POW || kind(n) == syntax::N_UNIT_POW {
        let e = if syntax::tokens[tok(n) + 1].kind == syntax::TK_MINUS { tok(n) + 2 } else { tok(n) + 1 };
        let x = syntax::tokens[e];
        if x.start + x.len > hi { hi = x.start + x.len; }
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

// A factor to SI and the exponents of length, mass and time, in halves
// (spec §6.1).
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
    BaseUnit { name: "m", u: Unit { f: 1.0, l: 2, m: 0, t: 0 } },
    BaseUnit { name: "g", u: Unit { f: 1.0e-3, l: 0, m: 2, t: 0 } },
    BaseUnit { name: "s", u: Unit { f: 1.0, l: 0, m: 0, t: 2 } },
    BaseUnit { name: "N", u: Unit { f: 1.0, l: 2, m: 2, t: -4 } },
    BaseUnit { name: "Pa", u: Unit { f: 1.0, l: -2, m: 2, t: -4 } },
    BaseUnit { name: "t", u: Unit { f: 1000.0, l: 0, m: 2, t: 0 } },
    BaseUnit { name: "bar", u: Unit { f: 1.0e5, l: -2, m: 2, t: -4 } },
    // skein-elements: the heat units of src/elements/gen.goose. Each system uses the kelvin.
    BaseUnit { name: "W", u: Unit { f: 1.0, l: 4, m: 2, t: -6 } },
    BaseUnit { name: "K", u: Unit { f: 1.0, l: 0, m: 0, t: 0 } },
    // The frequency of eigenvalue(k) of a modal step (src/elements/eigen.goose).
    BaseUnit { name: "Hz", u: Unit { f: 1.0, l: 0, m: 0, t: -2 } },
    // The units of the edge physics (src/elements/gen_edge.goose). Each system uses the ampere,
    // so A has no exponent here, and V is W/A. deg is an angle in rad.
    BaseUnit { name: "A", u: Unit { f: 1.0, l: 0, m: 0, t: 0 } },
    BaseUnit { name: "V", u: Unit { f: 1.0, l: 4, m: 2, t: -6 } },
    BaseUnit { name: "deg", u: Unit { f: 0.017453292519943295, l: 0, m: 0, t: 0 } },
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
                "use m, g, t, s, N, Pa, bar, W, K, Hz, A, V or deg, with a prefix n, u, m, c, k, M or G, for example mm, kN or GPa");
    return Unit { f: 1.0, l: 0, m: 0, t: 0 };
}

recursive fn unit_value(n: i64) -> Unit {
    let k = kind(n);
    if k == syntax::N_UNIT { return unit_value(syntax::nodes[n].first); }
    if k == syntax::N_UNIT_NAME { return unit_name(tok(n)); }
    if k == syntax::N_UNIT_POW {
        let b = unit_value(syntax::nodes[n].first);
        let e = syntax::nodes[n].aux;
        return Unit { f: pow(b.f, e / 2.0), l: b.l * e / 2, m: b.m * e / 2, t: b.t * e / 2 };
    }
    let a = unit_value(syntax::child(n, 0));
    let b = unit_value(syntax::child(n, 1));
    if k == syntax::N_UNIT_MUL { return Unit { f: a.f * b.f, l: a.l + b.l, m: a.m + b.m, t: a.t + b.t }; }
    return Unit { f: a.f / b.f, l: a.l - b.l, m: a.m - b.m, t: a.t - b.t };
}

// The factor from a unit to the unit system of the header.
fn to_system(n: i64) -> f64 {
    let u = unit_value(n);
    return u.f / (pow(system.l, u.l / 2.0) * pow(system.m, u.m / 2.0) * pow(system.t, u.t / 2.0));
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
        // A temperature in degC is an absolute temperature in K. The aux of
        // a degC value is its minus sign.
        let u = syntax::nodes[n].first;
        if u >= 0 && kind(syntax::nodes[u].first) == syntax::N_UNIT_NAME && syntax::text(tok(syntax::nodes[u].first)) == "degC" {
            let sign = if syntax::nodes[n].aux >= 0 { -1.0 } else { 1.0 };
            return scalar(sign * number(tok(n)) + 273.15);
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
        var v = find_value(name(n));
        if v < 0 && name(n) == "pi" { return scalar(PI); }
        if v < 0 { v = property(name(n)); }
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
    // The functions of one scalar.
    for f in ["sqrt", "exp", "ln", "sin", "cos", "tan", "abs"] {
        if !call_args(n, f, 1) { continue; }
        let a = eval(syntax::nodes[syntax::child(n, 1)].first);
        if a.n != 1 { unsupported(tok(n), str(f, " of a vector"), str("write ", f, " of a scalar")); }
        let x = a.v[0];
        return scalar(if f == "sqrt" { sqrt(x) } else if f == "exp" { exp(x) } else if f == "ln" { log(x) } else if f == "sin" { sin(x) }
                      else if f == "cos" { cos(x) } else if f == "tan" { tan(x) } else { abs(x) });
    }
    // A Point is its position.
    if call_args(n, "point", 1) { return eval(syntax::nodes[syntax::child(n, 1)].first); }
    if call_args(n, "interp", 3) { return scalar(interp(n)); }
    unsupported(tok(n), str("the expression '", source_of(n), "' is not in this slice"),
                "write a number with a unit, a name, a vector (x, y, z), or + - * / ^ of these");
    return scalar(0.0);
}

// The tables by name.
var table_names: u8[..64][>..];
var table_nodes: i64[>..];

// The material item whose properties a name can be: the material of the
// law that the generator evaluates now, or -1.
var scope_material = -1;

// The value of property `word` of scope_material, or -1.
fn property(word: u8[:]) -> i64 {
    guard scope_material >= 0 else { return -1; };
    var p = syntax::nodes[scope_material].first;
    while p >= 0 {
        if name(p) == word { return syntax::child(p, 1); }
        p = syntax::nodes[p].next;
    }
    return -1;
}

// The value of a table cell: a number, or - and a number.
recursive fn row_number(n: i64) -> f64 {
    if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" { return -row_number(syntax::nodes[n].first); }
    if kind(n) != syntax::N_NUMBER || syntax::nodes[n].first >= 0 {
        unsupported(tok(n), "a value of a table in this slice is a number", "write the number, and the unit in the column");
        return 0.0;
    }
    return number(tok(n));
}

// The values of column `word` of table t, in the unit system of the header.
fn table_column(t: i64, word: u8[:], at: i64, out: f64[>..]&) {
    var col = -1;
    var factor = 1.0;
    var i = 0;
    var c = syntax::nodes[t].first;
    while c >= 0 {
        if kind(c) == syntax::N_COLUMN {
            if name(c) == word {
                col = i;
                let u = syntax::child_of_kind(c, syntax::N_UNIT);
                if u >= 0 { factor = to_system(u); }
            }
            i++;
        } else if kind(c) == syntax::N_ROW && col >= 0 {
            out.push(row_number(syntax::child(c, col)) * factor);
        }
        c = syntax::nodes[c].next;
    }
    if col < 0 { unsupported(at, str("the table has no column ", word), "name a column of the table"); }
}

// The values of the column `m`, t.col of a table t, in `out`. False after
// an error.
fn column_of(m: i64, out: f64[>..]&) -> bool {
    let o = if kind(m) == syntax::N_MEMBER { syntax::nodes[m].first } else { -1 };
    let t = if o >= 0 && kind(o) == syntax::N_NAME { position(table_names, name(o)) } else { -1 };
    if t < 0 {
        unsupported(tok(m), "interp takes a column of a table", INTERP_FIX);
        return false;
    }
    table_column(table_nodes[t], name(m), tok(m), out);
    return ok();
}

let INTERP_FIX = "write interp(props.T, props.k, x), with two columns of a table and x in the range of the first";

// `interp(t.xs, t.ys, x)`: the piecewise linear function of the column ys
// of the column xs, in increasing order, at x.
recursive fn interp(n: i64) -> f64 {
    let x = eval(syntax::nodes[syntax::child(n, 3)].first).v[0];
    var xs: f64[>..] = [];
    var ys: f64[>..] = [];
    guard column_of(syntax::nodes[syntax::child(n, 1)].first, xs) && column_of(syntax::nodes[syntax::child(n, 2)].first, ys) else { return 0.0; };
    for i in 1..xs.len {
        if xs[i] <= xs[i - 1] { unsupported(tok(n), "the first column of interp is not in increasing order", INTERP_FIX); return 0.0; }
        if x >= xs[i - 1] && x <= xs[i] { return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]); }
    }
    unsupported(tok(n), str("interp at ", x, " is outside the first column"), INTERP_FIX);
    return 0.0;
}

// --- the model -------------------------------------------------------------------

let BOX = 0;
let CYLINDER = 1;
let TUBE = 2;       // a cylinder minus a cylinder on the same axis
let SECTOR = 3;     // a tube & a box that cuts the quarter x >= cx, y >= cy
let UNION = 4;      // a union of boxes, a | b, on the grid of runtime::mesh_union

// A solid. A box, and a union of boxes, span lo to hi in x, y and z. The
// others span lo to hi in r, θ and z about the axis x = cx, y = cy, which
// is parallel to z. The boxes of a union are union_boxes[6 first..6 (first
// + count)]: x0, y0, z0, x1, y1, z1 for each.
struct Shape {
    kind: i64,
    lo: f64[3],
    hi: f64[3],
    cx: f64,
    cy: f64,
    first: i64,
    count: i64,
}

var union_boxes: f64[>..];

// A face of the meshed solid: the part of its boundary where lattice
// coordinate `axis` (x, y, z for a box, r, θ, z for a sector) is `value`.
// It is the node set and the face set `set` of the mesh
// (src/runtime/api.goose): x0, x1, y0, y1, z0 or z1 on a box, r0, r1, t0,
// t1, z0 or z1 on a sector. A surface of a lane's mesh has axis -1.
struct Face {
    name: u8[..64],
    axis: i64,
    value: f64,
    set: u8[..64],
}

// A boundary condition on a face: `u = 0` or `u.x = 0` (`fix` holds a bit
// for each fixed component), or `sigma*n = t + s n` (`fix` is 0).
struct Bc {
    face: i64,      // index in `faces`
    fix: i64,
    t: f64[3],
    s: f64,
    tok: i64,
}

struct Report {
    expr: i64,      // the expression node
    label: u8[..256],
}

struct Step {
    name: u8[..64],
    procedure: i64, // the lane of the procedure (dispatch.goose)
    solver: u8[..64],         // the linear solver and the preconditioner of `solve`
    preconditioner: u8[..64],
    first_own: i64, // the bcs of the step itself: a range of `bcs`
    own_count: i64,
    first_bc: i64,  // all bcs that hold in the step: a range of `active`
    bc_count: i64,
    first_report: i64,
    report_count: i64,
    output: bool,
}

let EXPECT_U = 0;
let EXPECT_SIGMA = 1;
let EXPECT_REACTION = 2;
let EXPECT_SUM = 3;      // a sum of components of u at points (point_sum)

// An `expect` on a component of u or sigma at a point, or of the reaction
// on a face.
struct Expect {
    name: u8[..64],
    kind: i64,
    comp: i64,      // 0 to 2 for u and the reaction, 0 to 5 for sigma
    at: f64[3],
    face: i64,
    tok: i64,
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
    step: i64,      // the step that checks it at its end, or -1 for an item
    expr: i64,      // the expression of EXPECT_SUM
}

var model_name: u8[..64];
// The solids by name, and the solid that the model meshes.
var shape_names: u8[..64][>..];
var shapes: Shape[>..];
var solid_name: u8[..64];
var solid_tok = -1;
var domain = Shape { kind: BOX, lo: [0.0, 0.0, 0.0], hi: [0.0, 0.0, 0.0], cx: 0.0, cy: 0.0, first: 0, count: 0 };
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
var template = CORE;            // the lane of the template (dispatch.goose)
var mesh_lane = CORE;           // the lane that meshes the geometry: CORE (a box or a sector), MESHER or SHELLS
var cell_name: u8[..64];        // the mesh setting `cell`
var cell_tok = -1;
var order_tok = -1;
var incompatible = false;       // the method setting `family: incompatible`
var family_tok = -1;
// The `bc` items outside the steps, which hold in each step.
var global_bcs: i64[>..];

// A `let p: Point`: the face of its node set in `faces`, its position, and
// whether a bc holds it. The grid of a union of boxes has a grid plane
// through each Point that a bc holds, so a node is there.
struct PointSet {
    face: i64,
    at: f64[3],
    held: bool,
}

var point_sets: PointSet[>..];

// The materials by name, as the items.
var material_names: u8[..64][>..];
var material_nodes: i64[>..];

// A material of the physics and its region: `material m on r`. The first
// law is the material of each cell, and the region of each law, in order,
// takes its material. The strain is the thermal strain alpha*(T - T0),
// with the properties of the material.
struct Law {
    material: i64,  // index in material_names
    region: i64,    // index in shapes, or -1 for the whole solid
    young: f64,
    poisson: f64,
    strain: f64,
    tok: i64,
}

var laws: Law[>..];
var thermal_node = -1;  // the expression X of hooke(E, nu, strain(u) - X*I), or -1

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
    point_sets.clear();
    material_names.clear();
    material_nodes.clear();
    laws.clear();
    thermal_node = -1;
    table_names.clear();
    table_nodes.clear();
    scope_material = -1;
    model_name = [];
    shape_names.clear();
    shapes.clear();
    union_boxes.clear();
    solid_name = [];
    solid_tok = -1;
    material_name = [];
    physics_name = [];
    unknown_name = [];
    physics_tok = -1;
    template = CORE;
    mesh_lane = CORE;
    set64(cell_name, "hex");
    cell_tok = -1;
    order_tok = -1;
    incompatible = false;
    family_tok = -1;
    mesh_size = 0.0;
    mesh_tok = -1;
    order = 2;
    young = 0.0;
    poisson = -1.0;
    depth = 0;
    system = SYSTEMS[0];
}

// `s` in single quotes for the shell.
fn shell_quoted(s: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    out.push('\'');
    for c in s {
        if c == '\'' { out.append("'\\''"); } else { out.push(c); }
    }
    out.push('\'');
    return out;
}

// The directory of the results of model `model`, run with the --set words
// `sets`: .skein/results/<model>, and with --set, set/<words> in it. The
// words are joined by ',', with no spaces and with '_' for '/'. So a run
// with other values does not write over the results of the model, and a
// submodel finds the state of the case that it needs (spec §15.1).
fn results_path(model: u8[:], sets: const u8[:][:]) -> u8[>..] { return set_path(".skein/results/", model, sets); }

// The directory `root`<model>, and with --set, set/<words> in it.
fn set_path(root: u8[:], model: u8[:], sets: const u8[:][:]) -> u8[>..] {
    var p = str(root, model);
    for w, i in sets {
        p.append(if i == 0 { "/set/" } else { "," });
        for c in w {
            if c == '/' { p.push('_'); } else if c != ' ' { p.push(c); }
        }
    }
    return p;
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

let AXIS_NAMES = ["x", "y", "z"];

fn axis_of(word: u8[:]) -> i64 {
    if word == "x" { return 0; }
    if word == "y" { return 1; }
    if word == "z" { return 2; }
    return -1;
}

let BOX_FIX = "write the solid as box((x0, y0, z0) mm, (x1, y1, z1) mm), or as a quarter of a tube: (outer - bore) & box((cx, cy, z0) mm, (x1, y1, z1))";
let TUBE_FIX = "write cylinder(c, (0, 0, 1), r, h) - cylinder(c, (0, 0, 1), r_bore, h)";
let SECTOR_FIX = "cut the tube with a box from its centre: (outer - bore) & box((cx, cy, z0) mm, (cx + r, cy + r, z1))";

fn find_shape(word: u8[:]) -> i64 {
    for v, i in shape_names {
        if v == word { return i; }
    }
    return -1;
}

fn call_args(value: i64, fn_name: u8[:], count: i64) -> bool {
    return kind(value) == syntax::N_CALL && name(syntax::child(value, 0)) == fn_name &&
           syntax::child_count(value) == count + 1;
}

fn arg(value: i64, i: i64) -> Val { eval(syntax::nodes[syntax::child(value, i + 1)].first) }

// The shape of a solid expression: a box, a cylinder along z, a cylinder
// minus a cylinder, or a tube & a box, with names of solids in place of
// each.
recursive fn shape_of(value: i64) -> Shape {
    var out = Shape { kind: -1, lo: [0.0, 0.0, 0.0], hi: [0.0, 0.0, 0.0], cx: 0.0, cy: 0.0, first: 0, count: 0 };
    let k = kind(value);
    if k == syntax::N_NAME {
        let i = find_shape(name(value));
        if i < 0 { unsupported(tok(value), str("'", name(value), "' is not a solid"), BOX_FIX); return out; }
        return shapes[i];
    }
    if call_args(value, "box", 2) {
        let a = arg(value, 0);
        let b = arg(value, 1);
        if a.n != 3 || b.n != 3 { unsupported(tok(value), "box needs two corner points", BOX_FIX); return out; }
        out.kind = BOX;
        for j in 0..3 {
            out.lo[j] = min(a.v[j], b.v[j]);
            out.hi[j] = max(a.v[j], b.v[j]);
            if out.lo[j] == out.hi[j] { unsupported(tok(value), "the box has no volume", BOX_FIX); }
        }
        return out;
    }
    if call_args(value, "cylinder", 4) {
        let c = arg(value, 0);
        let axis = arg(value, 1);
        let r = arg(value, 2).v[0];
        let h = arg(value, 3).v[0];
        if c.n != 3 || axis.n != 3 || axis.v[0] != 0.0 || axis.v[1] != 0.0 || axis.v[2] <= 0.0 {
            unsupported(tok(value), "this slice has cylinders on the axis (0, 0, 1) only", TUBE_FIX);
            return out;
        }
        if r <= 0.0 || h <= 0.0 { unsupported(tok(value), "the cylinder has no volume", TUBE_FIX); return out; }
        return Shape { kind: CYLINDER, lo: [0.0, 0.0, c.v[2]], hi: [r, 2.0 * PI, c.v[2] + h], cx: c.v[0], cy: c.v[1], first: 0, count: 0 };
    }
    if k == syntax::N_BINARY && syntax::text(tok(value)) == "-" {
        let a = shape_of(syntax::child(value, 0));
        let b = shape_of(syntax::child(value, 1));
        let good = a.kind == CYLINDER && b.kind == CYLINDER && a.cx == b.cx && a.cy == b.cy &&
                   a.lo[2] == b.lo[2] && a.hi[2] == b.hi[2] && b.hi[0] < a.hi[0];
        if !good {
            if ok() {
                unsupported(tok(value), "this slice subtracts only a smaller cylinder on the same axis and with the same ends",
                            TUBE_FIX);
            }
            return out;
        }
        out = a;
        out.kind = TUBE;
        out.lo[0] = b.hi[0];
        return out;
    }
    if k == syntax::N_BINARY && syntax::text(tok(value)) == "&" {
        let a = shape_of(syntax::child(value, 0));
        let b = shape_of(syntax::child(value, 1));
        let tube = if a.kind == TUBE { a } else { b };
        let cut = if a.kind == TUBE { b } else { a };
        let good = tube.kind == TUBE && cut.kind == BOX && cut.lo[0] == tube.cx && cut.lo[1] == tube.cy &&
                   cut.hi[0] >= tube.cx + tube.hi[0] && cut.hi[1] >= tube.cy + tube.hi[0] &&
                   max(cut.lo[2], tube.lo[2]) < min(cut.hi[2], tube.hi[2]);
        if !good {
            if ok() { unsupported(tok(value), "this slice cuts only the quarter x >= cx, y >= cy of a tube", SECTOR_FIX); }
            return out;
        }
        out = tube;
        out.kind = SECTOR;
        out.hi[1] = 0.5 * PI;
        out.lo[2] = max(cut.lo[2], tube.lo[2]);
        out.hi[2] = min(cut.hi[2], tube.hi[2]);
        return out;
    }
    if k == syntax::N_BINARY && syntax::text(tok(value)) == "|" { return box_union(value); }
    unsupported(tok(value), str("the solid '", source_of(value), "' is not in this slice"), BOX_FIX);
    return out;
}

let UNION_FIX = "write a union of boxes, for example box((0, 0, 0) mm, (25, 5, 0.73) mm) | box((0, 0, 0.73) mm, (25, 5, 0.75) mm)";

// The boxes of shape `sh`, a box or a union, at the end of union_boxes. It
// gives the first box and the count.
fn copy_boxes(sh: Shape) -> i64, i64 {
    let first = union_boxes.len / 6;
    if sh.kind == BOX {
        for j in 0..3 { union_boxes.push(sh.lo[j]); }
        for j in 0..3 { union_boxes.push(sh.hi[j]); }
        return first, 1;
    }
    for x in 6 * sh.first..6 * (sh.first + sh.count) { union_boxes.push(union_boxes[x]); }
    return first, sh.count;
}

// `a | b` of two boxes or unions of boxes.
recursive fn box_union(value: i64) -> Shape {
    let a = shape_of(syntax::child(value, 0));
    let b = shape_of(syntax::child(value, 1));
    var out = a;
    if (a.kind != BOX && a.kind != UNION) || (b.kind != BOX && b.kind != UNION) {
        if ok() { unsupported(tok(value), "this slice has a union of boxes only", UNION_FIX); }
        return out;
    }
    out.kind = UNION;
    let fa, ca = copy_boxes(a);
    let fb, cb = copy_boxes(b);
    out.first = fa;
    out.count = ca + cb;
    for j in 0..3 {
        out.lo[j] = min(a.lo[j], b.lo[j]);
        out.hi[j] = max(a.hi[j], b.hi[j]);
    }
    return out;
}

// Whether `value` has the form of a solid that this package meshes: a
// name, a box, a cylinder, or a -, & or | of them (shape_of).
fn core_solid(value: i64) -> bool {
    let k = kind(value);
    if k == syntax::N_NAME || call_args(value, "box", 2) || call_args(value, "cylinder", 4) { return true; }
    return k == syntax::N_BINARY && (syntax::text(tok(value)) == "-" || syntax::text(tok(value)) == "&" ||
                                     syntax::text(tok(value)) == "|");
}

// `let wall: Solid = ...`. A solid of another form goes to the mesher lane
// (mesher_read in lane_mesher.goose), which then meshes the model.
fn solid(n: i64, value: i64) {
    if !core_solid(value) && solid_tok < 0 && mesher_read(n, value) {
        mesh_lane = MESHER;
        set64(solid_name, name(n));
        solid_tok = tok(n);
        return;
    }
    let sh = shape_of(value);
    if !ok() { return; }
    var w: u8[..64] = [];
    set64(w, name(n));
    shape_names.push(w);
    shapes.push(sh);
}

// The solid that the mesh, the physics and the surfaces are on. The first
// of these to name a solid chooses it. It gives false, after an error, if
// `region` is not that solid or cannot be meshed.
fn on_domain(region: i64, what: u8[:]) -> bool {
    if mesh_lane != CORE {
        if kind(region) == syntax::N_NAME && name(region) == solid_name { return true; }
        unsupported(tok(region), str(what, " is not on ", solid_name, ": this slice meshes one solid"),
                    str("put it on ", solid_name));
        return false;
    }
    let i = if kind(region) == syntax::N_NAME { find_shape(name(region)) } else { -1 };
    if i < 0 {
        unsupported(tok(region), str(what, " is not on a solid that a let names"), BOX_FIX);
        return false;
    }
    if solid_tok >= 0 {
        if name(region) == solid_name { return true; }
        unsupported(tok(region), str(what, " is not on ", solid_name, ": this slice meshes one solid"),
                    str("put it on ", solid_name));
        return false;
    }
    let sh = shapes[i];
    if sh.kind != BOX && sh.kind != SECTOR && sh.kind != UNION {
        unsupported(tok(region), str("the solid '", name(region), "' is not a box, a union of boxes or a quarter of a tube"),
                    if sh.kind == TUBE { SECTOR_FIX } else { BOX_FIX });
        return false;
    }
    domain = sh;
    set64(solid_name, name(region));
    solid_tok = tok(region);
    return true;
}

let FACE_FIX = "write the surface as boundary(beam, x == 0 mm), with x, y or z equal to a face of the solid, or as boundary(wall, bore) for a cylinder face of a tube";

fn near(a: f64, b: f64, size: f64) -> bool { abs(a - b) <= 1.0e-9 * size }

// The face of the domain where `axis` == v, as a face of the lattice, or
// -1 for the axis if it is not a face.
fn plane_face(axis: i64, v: f64) -> i64, f64 {
    let d = domain;
    if d.kind == BOX || d.kind == UNION {
        let size = d.hi[axis] - d.lo[axis];
        if near(v, d.lo[axis], size) || near(v, d.hi[axis], size) { return axis, v; }
        return -1, 0.0;
    }
    let r = d.hi[0];
    if axis == 0 && near(v, d.cx, r) { return 1, d.hi[1]; }
    if axis == 1 && near(v, d.cy, r) { return 1, d.lo[1]; }
    if axis == 2 && (near(v, d.lo[2], r) || near(v, d.hi[2], r)) { return 2, v; }
    return -1, 0.0;
}

// The names of the lattice axes in the set names of the mesh.
let BOX_SETS = ["x", "y", "z"];
let SECTOR_SETS = ["r", "t", "z"];

// `let root: Surface = boundary(beam, x == value)`, or
// `let inner: Surface = boundary(wall, bore)` with bore a cylinder. A
// surface of a geometry of the mesher lane goes to mesher_surface
// (lane_mesher.goose), and of the shells lane to shells_surface
// (lane_shells.goose).
fn surface(n: i64, value: i64) {
    if mesh_lane != CORE {
        var g = Face { name: [], axis: -1, value: 0.0, set: [] };
        set64(g.name, name(n));
        let found = if mesh_lane == MESHER { mesher_surface(n, value, g.set) } else { shells_surface(n, value, g.set) };
        if !found {
            unsupported(tok(value), str("the surface '", name(n), "' is not a surface of the mesh"), FACE_FIX);
            return;
        }
        faces.push(g);
        return;
    }
    if !call_args(value, "boundary", 2) {
        unsupported(tok(value), str("the surface '", name(n), "' is not one face of the solid"), FACE_FIX);
        return;
    }
    let region = syntax::nodes[syntax::child(value, 1)].first;
    let cond = syntax::nodes[syntax::child(value, 2)].first;
    if !on_domain(region, "the surface") { return; }
    var axis = -1;
    var at = 0.0;
    let ci = if kind(cond) == syntax::N_NAME { find_shape(name(cond)) } else { -1 };
    if ci >= 0 && domain.kind == SECTOR {
        let c = shapes[ci];
        let r = c.hi[0];
        let good = c.kind == CYLINDER && c.cx == domain.cx && c.cy == domain.cy;
        if good && (near(r, domain.lo[0], r) || near(r, domain.hi[0], r)) { axis = 0; at = r; }
    } else if kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" &&
              kind(syntax::child(cond, 0)) == syntax::N_NAME {
        let a = axis_of(name(syntax::child(cond, 0)));
        if a >= 0 {
            let v = eval(syntax::child(cond, 1));
            let pa, pv = plane_face(a, v.v[0]);
            axis = pa;
            at = pv;
        }
    }
    if axis < 0 {
        unsupported(tok(cond), str("the surface is not on a face of ", solid_name), FACE_FIX);
        return;
    }
    var f = Face { name: [], axis: axis, value: at, set: [] };
    set64(f.name, name(n));
    let side = if near(at, domain.lo[axis], domain.hi[axis] - domain.lo[axis]) { 0 } else { 1 };
    set64(f.set, str(if domain.kind == SECTOR { SECTOR_SETS[axis] } else { BOX_SETS[axis] }, side));
    faces.push(f);
}

// `let p: Point = point((x, y, z) mm)`: a value, and the node set of a face
// that a bc can hold.
fn point(n: i64, value: i64) {
    if !call_args(value, "point", 1) {
        unsupported(tok(value), str("the point '", name(n), "' is not point((x, y, z) mm)"), "write let p: Point = point((0, 0, 0) mm)");
        return;
    }
    let at = eval(value);
    if at.n != 3 { unsupported(tok(value), "the point is not a vector", "write point((x, y, z) mm)"); }
    var f = Face { name: [], axis: -1, value: 0.0, set: [] };
    set64(f.name, name(n));
    set64(f.set, str("point_", name(n)));
    faces.push(f);
    point_sets.push(PointSet { face: faces.len - 1, at: at.v, held: false });
}

// `let` and `param`.
fn value_item(n: i64) {
    let ty = syntax::child(n, 0);
    let value = syntax::child(n, 1);
    let t = if kind(ty) == syntax::N_TYPE_NAME { name(ty) } else { "" };
    let shell_type = t == "Surface" || t == "Curve";
    // A model with no solid can have a shell or beam geometry of the shells lane.
    if shell_type && solid_tok < 0 && shells_read_geometry(n, value) {
        mesh_lane = SHELLS;
        set64(solid_name, name(n));
        solid_tok = tok(n);
        return;
    }
    if t == "Solid" { solid(n, value); return; }
    if t == "Surface" || (shell_type && mesh_lane == SHELLS) { surface(n, value); return; }
    if t == "Point" { point(n, value); }
    if t == "Volume" || t == "Curve" || t == "Region" {
        unsupported(tok(n), str("a ", t, " is not in this slice"), "use a Solid box and Surface faces of it");
        return;
    }
    var w: u8[..64] = [];
    set64(w, name(n));
    value_names.push(w);
    value_nodes.push(value);
}

// A material item. The physics evaluates its properties (laws).
fn material(n: i64) {
    var w: u8[..64] = [];
    set64(w, name(n));
    material_names.push(w);
    material_nodes.push(n);
    if material_name.len == 0 { set64(material_name, name(n)); }
    var p = syntax::nodes[n].first;
    while p >= 0 {
        if name(p) != "E" && name(p) != "nu" && name(p) != "alpha" && name(p) != "rho" {
            unsupported(tok(p), str("the property '", name(p), "' is not in this slice"),
                        "give the material only E, nu and alpha (and rho, which a static step does not use)");
        }
        p = syntax::nodes[p].next;
    }
}

// The E, nu and the thermal strain of each law, with the properties of its
// material. A model with one material and no material statement in the
// physics gives that material to the whole solid.
fn evaluate_laws() {
    if laws.len == 0 && material_nodes.len == 1 {
        laws.push(Law { material: 0, region: -1, young: 0.0, poisson: -1.0, strain: 0.0, tok: tok(material_nodes[0]) });
    }
    for &w in laws {
        scope_material = material_nodes[w.material];
        let e = property("E");
        let nu = property("nu");
        w.young = if e >= 0 { eval(e).v[0] } else { 0.0 };
        w.poisson = if nu >= 0 { eval(nu).v[0] } else { -1.0 };
        if thermal_node >= 0 { w.strain = eval(thermal_node).v[0]; }
        scope_material = -1;
        if w.young <= 0.0 || w.poisson <= -1.0 || w.poisson >= 0.5 {
            unsupported(w.tok, "the material needs E > 0 and -1 < nu < 0.5", "give E: Pressure = 210 GPa and nu: Real = 0.3");
        }
    }
    if laws.len > 0 {
        young = laws[0].young;
        poisson = laws[0].poisson;
        set64(material_name, material_names[laws[0].material]);
    }
}

// `material m, n on r, ...` of the physics.
fn material_laws(c: i64) {
    var m = syntax::nodes[c].first;
    while m >= 0 {
        let i = position(material_names, name(m));
        let r = syntax::nodes[m].first;
        let region = if r >= 0 && kind(r) == syntax::N_NAME { find_shape(name(r)) } else { -1 };
        if i < 0 {
            unsupported(tok(m), str("'", name(m), "' is not a material of the model"), PHYSICS_FIX);
        } else if r >= 0 && (region < 0 || (shapes[region].kind != BOX && shapes[region].kind != UNION)) {
            unsupported(tok(r), "the region of a material is not a box or a union of boxes that a let names", PHYSICS_FIX);
        } else {
            laws.push(Law { material: i, region: region, young: 0.0, poisson: -1.0, strain: 0.0, tok: tok(m) });
        }
        m = syntax::nodes[m].next;
    }
}

fn mesh(n: i64) {
    if mesh_tok >= 0 {
        unsupported(tok(n), "a second mesh: this slice has one", "remove the second mesh");
        return;
    }
    mesh_tok = tok(n);
    let region = syntax::child(n, 0);
    on_domain(region, "the mesh");
    var s = syntax::nodes[syntax::child(n, 1)].first;
    while s >= 0 {
        let value = syntax::nodes[s].first;
        if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "refine, layers and crack_tip are not in this slice", "remove the line");
        } else if name(s) == "size" {
            mesh_size = eval(value).v[0];
        } else if name(s) == "cell" {
            if kind(value) != syntax::N_NAME {
                unsupported(tok(value), "this slice meshes with hex cells only", "write cell: hex");
            }
            set64(cell_name, name(value));
            cell_tok = tok(value);
        } else {
            unsupported(tok(s), str("the mesh setting '", name(s), "' is not in this slice"), "use only size and cell");
        }
        s = syntax::nodes[s].next;
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 5 mm"); }
}

let PHYSICS_FIX = "write physics solid on beam: linear_elastic solves u: Vector<Length> { material steel }, or the statements of it: let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) - alpha*dT*I) and eq div(sigma) = 0";

fn physics(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: this slice has one", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    // A template of a lane goes to its generator (dispatch.goose).
    var c0 = syntax::nodes[n].first;
    while c0 >= 0 {
        if kind(c0) == syntax::N_USES && template_lane(name(c0)) > CORE {
            template = template_lane(name(c0));
            if !read_template(template, name(c0), n) {
                unsupported(tok(n), str("the template ", name(c0), " has no generator in this build"),
                            PHYSICS_FIX);
            }
            return;
        }
        c0 = syntax::nodes[c0].next;
    }
    var uses = false;
    var hooke = false;
    var eqs = 0;
    var c = syntax::nodes[n].first;
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            let r = syntax::nodes[c].first;
            on_domain(r, "the physics");
        } else if k == syntax::N_USES {
            uses = name(c) == "linear_elastic" && syntax::nodes[c].first < 0;
        } else if k == syntax::N_UNKNOWN {
            if unknown_name.len > 0 { unsupported(tok(c), "a second unknown", PHYSICS_FIX); }
            set64(unknown_name, name(c));
        } else if k == syntax::N_MATASSIGN {
            material_laws(c);
        } else if k == syntax::N_ASSUME {
            var a = syntax::nodes[c].first;
            while a >= 0 {
                if name(a) != "small_strain" && name(a) != "isotropic" {
                    unsupported(tok(a), str("the assumption '", name(a), "' is not in this slice"), "use small_strain and isotropic");
                }
                a = syntax::nodes[a].next;
            }
        } else if k == syntax::N_LET && name(c) == "sigma" {
            hooke = hooke_let(syntax::child(c, 1));
        } else if k == syntax::N_EQN && source_of(syntax::child(c, 0)) == "div(sigma)" && source_of(syntax::child(c, 1)) == "0" {
            eqs++;
        } else {
            unsupported(tok(c), "this slice has the template linear_elastic, or the statements of it: material, assume, let sigma and eq",
                        PHYSICS_FIX);
        }
        c = syntax::nodes[c].next;
    }
    if !uses && !(hooke && eqs == 1) {
        unsupported(tok(n), "this slice solves only the template linear_elastic, or let sigma = hooke(...) with eq div(sigma) = 0", PHYSICS_FIX);
    }
}

// `hooke(E, nu, strain(u))` or `hooke(E, nu, strain(u) - X*I)` with X the
// thermal strain, for example alpha*(T - T0) with constant temperatures
// (spec §10.4). It sets thermal_node to X. False after an error.
fn hooke_let(v: i64) -> bool {
    let strain = str("strain(", unknown_name, ")");
    let good = call_args(v, "hooke", 3) && source_of(syntax::nodes[syntax::child(v, 1)].first) == "E" &&
               source_of(syntax::nodes[syntax::child(v, 2)].first) == "nu";
    let e = if good { syntax::nodes[syntax::child(v, 3)].first } else { -1 };
    if e >= 0 && source_of(e) == strain { return true; }
    if e >= 0 && kind(e) == syntax::N_BINARY && syntax::text(tok(e)) == "-" && source_of(syntax::child(e, 0)) == strain {
        let t = syntax::child(e, 1);
        if kind(t) == syntax::N_BINARY && syntax::text(tok(t)) == "*" && kind(syntax::child(t, 1)) == syntax::N_NAME &&
           name(syntax::child(t, 1)) == "I" {
            thermal_node = syntax::child(t, 0);
            return true;
        }
    }
    unsupported(tok(v), "this slice has sigma = hooke(E, nu, strain(u)), with - alpha*(T - T0)*I for a thermal strain", PHYSICS_FIX);
    return false;
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
            order_tok = tok(s);
        } else if kind(s) == syntax::N_SETTING && name(s) == "family" {
            let f = syntax::nodes[s].first;
            family_tok = tok(f);
            incompatible = name(f) == "incompatible";
            if !incompatible && name(f) != "lagrange" {
                unsupported(family_tok, "this slice has only the families lagrange and incompatible",
                            "write family: incompatible, or remove the family");
            }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this slice", "use only mesh, order and family");
        }
        s = syntax::nodes[s].next;
    }
}

let BC_FIX = "write bc root: u = 0, bc sym: u.x = 0, bc tip: sigma*n = (0, 0, -1) MPa, or bc inner: sigma*n = -p*n";

// The scale s of `value` == s*n, with n the unit normal. False if `value`
// is not a scalar times n.
recursive fn normal_scale(value: i64) -> bool, f64 {
    let k = kind(value);
    if k == syntax::N_NAME && name(value) == "n" && find_value("n") < 0 { return true, 1.0; }
    if k == syntax::N_UNARY && syntax::text(tok(value)) == "-" {
        let good, s = normal_scale(syntax::nodes[value].first);
        return good, -s;
    }
    if k == syntax::N_BINARY && syntax::text(tok(value)) == "*" {
        for i in 0..2 {
            let good, s = normal_scale(syntax::child(value, i));
            if good {
                let other = eval(syntax::child(value, 1 - i));
                if other.n == 1 { return true, s * other.v[0]; }
            }
        }
    }
    return false, 0.0;
}

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let region = syntax::child(n, 0);
    let target = syntax::child(n, 1);
    let value = syntax::child(n, 2);
    let face = if kind(region) == syntax::N_NAME { find_face(name(region)) } else { -1 };
    if face < 0 {
        unsupported(tok(region), str("the region of the bc is not a face of ", solid_name, " that a let names"), FACE_FIX);
        return -1;
    }
    var b = Bc { face: face, fix: 0, t: [0.0, 0.0, 0.0], s: 0.0, tok: tok(n) };
    for &q in point_sets {
        if q.face == face { q.held = true; }
    }
    let comp = if kind(target) == syntax::N_MEMBER && kind(syntax::nodes[target].first) == syntax::N_NAME &&
                  name(syntax::nodes[target].first) == unknown_name { axis_of(name(target)) } else { -1 };
    if (kind(target) == syntax::N_NAME && name(target) == unknown_name) || comp >= 0 {
        let v = eval(value);
        b.fix = if comp >= 0 { 1 << comp } else { 7 };
        if comp >= 0 && v.n != 1 { unsupported(tok(value), "a component of u is a scalar", BC_FIX); }
        for j in 0..3 {
            if v.v[j] != 0.0 {
                unsupported(tok(value), "this slice fixes u = 0 only", BC_FIX);
            }
        }
    } else if source_of(target) == "sigma*n" && faces[face].axis < 0 && mesh_lane == CORE {
        unsupported(tok(region), "this slice has no load on a Point", "put the load on a Surface");
        return -1;
    } else if source_of(target) == "sigma*n" {
        let normal, scale = normal_scale(value);
        if normal {
            b.s = scale;
        } else {
            let v = eval(value);
            if v.n != 3 { unsupported(tok(value), "the traction is a vector", BC_FIX); }
            b.t = v.v;
        }
    } else {
        unsupported(tok(target), str("the bc target '", source_of(target), "' is not in this slice"), BC_FIX);
        return -1;
    }
    bcs.push(b);
    return bcs.len - 1;
}

fn step(n: i64) {
    let proc = syntax::child(n, 0);
    var s = Step { name: [], procedure: CORE, solver: [], preconditioner: [], first_own: bcs.len, own_count: 0,
                   first_bc: 0, bc_count: 0, first_report: reports.len, report_count: 0, output: false };
    set64(s.name, name(n));
    set64(s.solver, "cg");
    set64(s.preconditioner, "amg");   // spec §14.4
    // A procedure of a lane goes to its generator (dispatch.goose).
    let p = procedure_lane(name(proc));
    if p > CORE {
        s.procedure = p;
        steps.push(s);
        if !read_procedure(p, name(proc), n) {
            unsupported(tok(proc), str("the procedure ", name(proc), " has no generator in this build"),
                        "write step load: static { ... }");
        }
        return;
    }
    if name(proc) != "static" || syntax::nodes[proc].first >= 0 {
        unsupported(tok(proc), "this slice has only static steps with no arguments", "write step load: static { ... }");
    }
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
        } else if k == syntax::N_SOLVE {
            solve(c, s);
        } else if k == syntax::N_EXPECT {
            expect(c, steps.len);
        } else {
            unsupported(tok(c), "the step statement is not in this slice", "use only bc, output, report, expect and solve");
        }
        c = syntax::nodes[c].next;
    }
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

let SOLVE_FIX = "write solve { linear: cg, precond: amg }, or precond: jacobi";

// `solve { linear: <solver>, precond: <preconditioner> }`. A pair other
// than cg with jacobi goes to the solvers lane (lane_solvers.goose).
fn solve(n: i64, s: Step&) {
    var c = syntax::nodes[syntax::nodes[n].first].first;
    while c >= 0 {
        let value = syntax::nodes[c].first;
        let known = kind(c) == syntax::N_SETTING && (name(c) == "linear" || name(c) == "precond");
        if !known || kind(value) != syntax::N_NAME {
            unsupported(tok(c), "the solve setting is not in this slice", SOLVE_FIX);
        } else if name(c) == "linear" {
            set64(s.solver, name(value));
        } else {
            set64(s.preconditioner, name(value));
        }
        c = syntax::nodes[c].next;
    }
    if !(s.solver == "cg" && s.preconditioner == "jacobi") && !solvers_known(s.solver, s.preconditioner) {
        unsupported(tok(n), str("the linear solver ", s.solver, " with the preconditioner ", s.preconditioner,
                                " is not in this build"), SOLVE_FIX);
    }
}

let REPORT_FIX = "write report reaction(root), or report u.z at (100, 5, 5) mm";

// The kinds of a report: the reaction on a face, or u or one component of u
// at a point.
fn report(item: i64) {
    let e = syntax::nodes[item].first;
    if kind(e) == syntax::N_CALL && name(syntax::child(e, 0)) == "reaction" && syntax::child_count(e) == 2 {
        let r = syntax::nodes[syntax::child(e, 1)].first;
        let face = if kind(r) == syntax::N_NAME { find_face(name(r)) } else { -1 };
        if face < 0 { unsupported(tok(r), "reaction needs a face of the solid that a let names", REPORT_FIX); return; }
    } else if kind(e) == syntax::N_AT {
        point_field(e);
    } else if !point_sum(e) {
        if ok() { unsupported(tok(e), str("the report '", source_of(e), "' is not in this slice"), REPORT_FIX); }
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

// Whether `e` is a sum of components of u at points: u.z at p, and + and -
// of these, for example u.z at corner - u.z at centre.
recursive fn point_sum(e: i64) -> bool {
    let k = kind(e);
    if k == syntax::N_AT {
        let c = point_field(e);
        if c == 3 { unsupported(tok(e), "a sum of values at points adds components", str("write ", unknown_name, ".z")); }
        return c >= 0 && c < 3;
    }
    let op = syntax::text(tok(e));
    if k == syntax::N_BINARY && (op == "+" || op == "-") { return point_sum(syntax::child(e, 0)) && point_sum(syntax::child(e, 1)); }
    if k == syntax::N_UNARY && op == "-" { return point_sum(syntax::nodes[e].first); }
    return false;
}

// Emits the samples of the point sum `e` (point_sum), from sample `at` on,
// and appends its value in Goose to `expr`.
recursive fn emit_point_sum(out: u8[>..]&, e: i64, at: i64&, expr: u8[>..]&) {
    let k = kind(e);
    if k == syntax::N_AT {
        at++;
        // The block of the report or the expect indents the first line.
        format(out, if at == 1 { "    let ux" } else { "        let ux" }, at, ", uy", at, ", uz", at, " = runtime::displacement_at(");
        lits(out, eval(syntax::child(e, 1)).v);
        out.append(");\n");
        format(expr, "u", AXIS_NAMES[point_field(e)], at);
        return;
    }
    if k == syntax::N_UNARY {
        expr.append("-(");
        emit_point_sum(out, syntax::nodes[e].first, at, expr);
        expr.push(')');
        return;
    }
    expr.push('(');
    emit_point_sum(out, syntax::child(e, 0), at, expr);
    format(expr, " ", syntax::text(tok(e)), " ");
    emit_point_sum(out, syntax::child(e, 1), at, expr);
    expr.push(')');
}

let SIGMA_COMPONENTS = ["xx", "yy", "zz", "yz", "xz", "xy"];
let EXPECT_FIX = "write expect tip: u.z at (100, 5, 5) mm == -0.19 mm within 2%, or sigma.xx at a point, or reaction(root).z";

// An `expect` item (step -1), or an `expect` of step `step`.
fn expect(n: i64, step: i64) {
    let left = syntax::child(n, 0);
    let right = syntax::child(n, 1);
    let within = syntax::child(n, 2);
    var label: u8[>..] = [];
    format(label, source_of(right), " within ", source_of(within));
    var x = Expect { name: [], kind: EXPECT_U, comp: -1, at: [0.0, 0.0, 0.0], face: -1, tok: tok(n), want: 0.0,
                     within: 0.0, relative: kind(within) == syntax::N_PERCENT, label: [], step: step, expr: left };
    let f = if kind(left) == syntax::N_AT { syntax::child(left, 0) } else { left };
    let object = if kind(f) == syntax::N_MEMBER { syntax::nodes[f].first } else { -1 };
    if kind(left) == syntax::N_AT && object >= 0 && kind(object) == syntax::N_NAME && name(object) == "sigma" {
        x.kind = EXPECT_SIGMA;
        let c = position(SIGMA_COMPONENTS, name(f));
        if c < 0 {
            unsupported(tok(f), str("sigma has no component '", name(f), "'"), "use xx, yy, zz, yz, xz or xy");
            return;
        }
        x.comp = c;
        x.at = eval(syntax::child(left, 1)).v;
    } else if kind(left) == syntax::N_AT && object >= 0 && name(object) != unknown_name {
        // An expect on a field this slice cannot sample stays, and the
        // report says that it did not run.
        var s: u8[..256] = [];
        set256(s, str(name(n), ": ", source_of(left), " is not in this slice"));
        skipped.push(s);
        return;
    } else if kind(left) == syntax::N_AT {
        x.comp = point_field(left);
        if x.comp == 3 {
            unsupported(tok(left), "an expect compares one component", str("write ", unknown_name, ".z"));
            return;
        }
        if x.comp < 0 { return; }
        x.at = eval(syntax::child(left, 1)).v;
    } else if object >= 0 && kind(object) == syntax::N_CALL && name(syntax::child(object, 0)) == "reaction" &&
              syntax::child_count(object) == 2 {
        x.kind = EXPECT_REACTION;
        let r = syntax::nodes[syntax::child(object, 1)].first;
        x.face = if kind(r) == syntax::N_NAME { find_face(name(r)) } else { -1 };
        x.comp = axis_of(name(f));
        if x.face < 0 || x.comp < 0 {
            unsupported(tok(left), "this slice checks reaction(face).x, .y or .z of a face that a let names", EXPECT_FIX);
            return;
        }
    } else if kind(left) == syntax::N_BINARY || kind(left) == syntax::N_UNARY {
        x.kind = EXPECT_SUM;
        if !point_sum(left) {
            if ok() { unsupported(tok(left), "this slice checks a sum of components of u at points, with + and -", EXPECT_FIX); }
            return;
        }
    } else {
        unsupported(tok(left), "this slice checks only u or sigma components at a point, or a reaction component",
                    EXPECT_FIX);
        return;
    }
    x.want = eval(right).v[0];
    x.within = eval(within).v[0];
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
        } else if k == syntax::N_TABLE {
            var w: u8[..64] = [];
            set64(w, name(c));
            table_names.push(w);
            table_nodes.push(c);
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
            expect(c, -1);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_STUDY {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this slice"),
                        "remove it: this slice runs one linear static box");
        }
        c = syntax::nodes[c].next;
    }
    if ok() && template == CORE { evaluate_laws(); }
    let first = syntax::nodes[prog].first;
    let at = if first >= 0 { tok(first) } else { tok(prog) };
    if solid_tok < 0 { unsupported(at, "the model has no solid to mesh", BOX_FIX); }
    if mesh_tok < 0 { unsupported(at, "the model has no mesh", "add mesh m on beam { size: 5 mm, cell: hex }"); }
    if physics_tok < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if material_name.len == 0 { unsupported(at, "the model has no material", "add material steel { E: Pressure = 210 GPa, nu: Real = 0.3 }"); }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { bc root: u = 0 ... }"); }
    // A cell other than hex of order 2 goes to the lane that meshes the
    // geometry: the mesher, the shells or the elements lane (lane_*.goose).
    let lane_cell = if mesh_lane == MESHER { mesher_cell(cell_name, order) }
                    else if mesh_lane == SHELLS { shells_cell(cell_name, order) }
                    else { elements_cell(cell_name, order) };
    for q in point_sets {
        if q.held && domain.kind != UNION {
            unsupported(bc_tok_of(q.face), "a bc on a Point needs a union of boxes, whose grid has a node at the point",
                        "write the solid as a union of boxes, for example a | b");
        }
    }
    if incompatible {
        // The runtime has the incompatible modes on the box mesh only (hex8i).
        if thermal_node >= 0 {
            unsupported(family_tok, "this slice has no thermal strain on the incompatible modes", "remove the family");
        }
        if template != CORE || domain.kind != BOX || cell_name != "hex" || order != 1 {
            unsupported(family_tok, "this slice has the incompatible modes only for linear_elastic on hex cells of order 1 on a box",
                        "write cell: hex and order: 1, and mesh a box");
        }
    } else if (cell_name != "hex" || order != 2) && !lane_cell {
        if cell_name != "hex" {
            unsupported(cell_tok, "this slice meshes with hex cells only", "write cell: hex");
        } else {
            unsupported(order_tok, "this slice has only order 2, because linear cells lock in bending",
                        "write order: 2");
        }
    }
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

fn same_plane(a: i64, b: i64) -> bool { faces[a].set == faces[b].set }

// A fixed face and its fixed components, one bit for each.
struct Fixed {
    face: i64,
    mask: i64,
}

// The fixed faces of step `si`, one for each plane.
fn fixed_faces(si: i64) -> Fixed[>..] {
    let s = steps[si];
    var out: Fixed[>..] = [];
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix == 0 { continue; }
        var seen = false;
        for &f in out {
            if same_plane(f.face, b.face) {
                f.mask |= b.fix;
                seen = true;
            }
        }
        if !seen { out.push(Fixed { face: b.face, mask: b.fix }); }
    }
    return out;
}

// The fixed components of face `face` in step `si`.
fn fixed_mask(si: i64, face: i64) -> i64 {
    for f in fixed_faces(si) {
        if same_plane(f.face, face) { return f.mask; }
    }
    return 0;
}

fn same_fixed(a: Fixed[:], b: Fixed[:]) -> bool {
    guard a.len == b.len else { return false; }
    for f, i in a {
        if f.face != b[i].face || f.mask != b[i].mask { return false; }
    }
    return true;
}

// The runtime keeps the supports from step to step, so each step must fix
// the same faces. The check equilibrium uses the total reaction, which
// counts each fixed node once, so two fixed faces can share nodes. A
// reaction of a fixed face that shares nodes with another fixed face that
// fixes the same component fails when the program runs (emit_alone).
fn supports() {
    let a = fixed_faces(0);
    if a.len == 0 {
        unsupported(physics_tok, str("no face is fixed, so ", solid_name, " can move as a rigid body"),
                    "add bc root: u = 0");
        return;
    }
    for si in 1..steps.len {
        let b = fixed_faces(si);
        if !same_fixed(a, b) {
            unsupported(bc_tok(si, b[b.len - 1].face), str("step ", steps[si].name, " fixes a face that step ",
                                                           steps[0].name, " does not fix"),
                        "fix each face in the first step, or put the bc u = 0 outside the steps");
        }
    }
    for si in 0..steps.len {
        let s = steps[si];
        for r in reports[s.first_report..s.first_report + s.report_count] {
            if kind(r.expr) != syntax::N_CALL { continue; }
            let arg = syntax::nodes[syntax::child(r.expr, 1)].first;
            if fixed_mask(si, find_face(name(arg))) == 0 {
                unsupported(tok(arg), str("the face ", name(arg), " is not fixed in step ", s.name,
                                          ", so it has no reaction"),
                            "report the reaction of a face with bc u = 0");
            }
        }
    }
    for x in expects {
        if x.kind != EXPECT_REACTION { continue; }
        if (fixed_mask(0, x.face) & (1 << x.comp)) == 0 {
            unsupported(x.tok, str("expect ", x.name, ": the face ", faces[x.face].name, " does not fix u.",
                                         AXIS_NAMES[x.comp], ", so it has no reaction in ", AXIS_NAMES[x.comp]),
                        "expect a component of the reaction that a bc u = 0 or u.x = 0 holds");
        }
    }
}

// The token of the first bc on face `f`.
fn bc_tok_of(f: i64) -> i64 {
    for b in bcs {
        if b.face == f { return b.tok; }
    }
    return physics_tok;
}

// The token of the fixed bc on face `f` in step `si`.
fn bc_tok(si: i64, f: i64) -> i64 {
    let s = steps[si];
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix != 0 && b.face == f { return b.tok; }
    }
    return physics_tok;
}

// --- the Goose program -------------------------------------------------------------

let AXES = ["runtime::X", "runtime::Y", "runtime::Z"];

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

// What a case of a study shares with case 1 (spec §20.4). Case 1 writes
// its first `stages` stages in `dir`, and each other case (`reads`) reads
// them in place of their work. Stage 1 is the mesh. Stage 2 is the solve of
// each step, so a case that reads it runs only the post-processing.
struct Share {
    dir: u8[..1024],
    stages: i64,
    reads: bool,
}

fn no_share() -> Share { return Share { dir: [], stages: 0, reads: false }; }

// The call `call(dir)`, or `call(dir, step)` when `step` is not empty, of a
// case that shares stages, made only when the condition `when` holds. The
// case stops when the call fails.
fn emit_share(out: u8[>..]&, when: u8[:], call: u8[:], share: Share, step: u8[:]) {
    format(out, "    if ", when, "!", call, "(");
    quoted(out, share.dir);
    if step.len > 0 {
        out.append(", ");
        quoted(out, step);
    }
    out.append(") {\n        write_stderr(");
    quoted(out, str("skein run: cannot share the stages of case 1 in '", share.dir, "'"));
    out.append(");\n        write_stderr(\"\\n\");\n        exit(1);\n    }\n");
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

    // A bc on a face that holds no node is a fault of the generator.
    fn face(found: bool, region: u8[:]) {
        if !found {
            write_stderr(str("skein run: internal fault: no node on the face '", region, "'. Please report it.\n"));
            exit(3);
        }
    }

    // The reaction of a fixed face that shares nodes with another fixed
    // face that fixes the same component is not the load on that face
    // alone, so the run reports an error.
    fn shared(face: u8[:], set: u8[:], other: u8[:], other_set: u8[:]) -> bool {
        guard runtime::sets_share_nodes(set, other_set) else { return false; }
        write_stderr(str("skein run: the fixed faces '", face, "' and '", other, "' share nodes and fix the ",
                         "same component of u, so reaction(", face, ") is not the load on '", face, "' alone. ",
                         "Report the reaction of a fixed face that shares no node with another fixed face. The ",
                         "check equilibrium uses the total reaction of all fixed faces.\n"));
        failed = true;
        return true;
    }

    fn solve(step: u8[:], solver: u8[:], preconditioner: u8[:]) {
        runtime::linear_solver(solver, preconditioner);
        converged = runtime::solve_static();
        format(solves, "  ", step, ",", runtime::cell_count(), ",", runtime::dof_count(), ",", runtime::cg_iterations(), ",",
               runtime::solve_seconds(), ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of the largest applied load.
    fn equilibrium(rx: f64, ry: f64, rz: f64, largest: f64, unit: u8[:]) {
        let fx, fy, fz = runtime::applied_load();
        let r = sqrt((rx + fx) * (rx + fx) + (ry + fy) * (ry + fy) + (rz + fz) * (rz + fz));
        let limit = 1.0e-6 * largest;
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
        format(s, "solve[", solve_count, "]{step,cells,dofs,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

fn emit_face(out: u8[>..]&, call: u8[:], f: Face) {
    format(out, "    face(runtime::", call, "(");
    quoted(out, f.set);
}

fn emit_face_end(out: u8[>..]&, f: Face) {
    out.append("), ");
    quoted(out, f.name);
    out.append(");\n");
}

fn emit_face_call(out: u8[>..]&, b: Bc) {
    let f = faces[b.face];
    if b.fix == 7 {
        emit_face(out, "fix", f);
        emit_face_end(out, f);
    } else if b.fix != 0 {
        for c in 0..3 {
            if (b.fix & (1 << c)) == 0 { continue; }
            emit_face(out, "fix_component", f);
            format(out, ", ", AXES[c]);
            emit_face_end(out, f);
        }
    } else if b.s != 0.0 {
        emit_face(out, "normal_traction", f);
        out.append(", ");
        lit(out, b.s);
        emit_face_end(out, f);
    } else {
        emit_face(out, "traction", f);
        out.append(", ");
        lits(out, b.t);
        emit_face_end(out, f);
    }
}

// The largest load of step `si`, in the variable `largest`: the size of a
// traction times the area of its face (runtime::face_area).
fn emit_largest_load(out: u8[>..]&, si: i64) {
    let s = steps[si];
    out.append("    var largest = 0.0;\n");
    // A thermal load sums to 0: its scale is the sum of its positive loads.
    if thermal_node >= 0 { out.append("    largest = max(largest, runtime::load_spread());\n"); }
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix != 0 { continue; }
        out.append("    largest = max(largest, ");
        lit(out, sqrt(b.t[0] * b.t[0] + b.t[1] * b.t[1] + b.t[2] * b.t[2]) + abs(b.s));
        out.append(" * runtime::face_area(");
        quoted(out, faces[b.face].set);
        out.append("));\n");
    }
}

// The components of a reaction that the fixed components `mask` hold, with
// 0.0 for each other component.
fn reaction_parts(mask: i64) -> u8[>..] {
    var out: u8[>..] = [];
    for c in 0..3 {
        if c > 0 { out.append(", "); }
        out.append(if (mask & (1 << c)) != 0 { AXIS_NAMES[c] } else { "0.0" });
    }
    return out;
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:], share: Share) {
    let s = steps[si];
    format(out, "// step ", s.name, ": static\n");
    format(out, "fn step_", s.name, "() {\n");
    if share.stages >= 2 && share.reads {
        emit_share(out, "", "runtime::read_solution", share, s.name);
    } else {
        emit_solve(out, si, share);
    }
    // The reaction on all fixed dofs, each once, for the check `equilibrium`.
    out.append("    let rx, ry, rz = runtime::total_reaction();\n");
    emit_largest_load(out, si);
    format(out, "    equilibrium(rx, ry, rz, largest, \"", system.force, "\");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] {
        emit_report(out, si, r);
    }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(");\n");
    }
    out.append("}\n\n");
}

// The supports and the loads of step `si`, and its solve.
fn emit_solve(out: u8[>..]&, si: i64, share: Share) {
    let s = steps[si];
    if si > 0 {
        out.append("    // The supports of the first step stay.\n");
        out.append("    runtime::clear_loads();\n");
    }
    for b in active[s.first_bc..s.first_bc + s.bc_count] {
        if b.fix != 0 && si > 0 { continue; }
        let f = faces[b.face];
        if b.fix == 7 {
            format(out, "    // bc ", f.name, ": ", unknown_name, " = 0\n");
        } else if b.fix != 0 {
            for c in 0..3 {
                if (b.fix & (1 << c)) != 0 { format(out, "    // bc ", f.name, ": ", unknown_name, ".", AXIS_NAMES[c], " = 0\n"); }
            }
        } else if b.s != 0.0 {
            format(out, "    // bc ", f.name, ": sigma*n = ", b.s, "*n ", system.pressure, "\n");
        } else {
            format(out, "    // bc ", f.name, ": sigma*n = (", b.t[0], ", ", b.t[1], ", ", b.t[2], ") ",
                   system.pressure, "\n");
        }
        emit_face_call(out, b);
    }
    if thermal_node >= 0 {
        out.append("    // the thermal strain of each material: alpha*(T - T0)\n    runtime::thermal_load();\n");
    }
    out.append("    solve(");
    quoted(out, s.name);
    out.append(", ");
    quoted(out, s.solver);
    out.append(", ");
    quoted(out, s.preconditioner);
    out.append(");\n");
    if share.stages >= 2 { emit_share(out, "converged && ", "runtime::write_solution", share, s.name); }
}

fn emit_point(out: u8[>..]&, at: f64[3]) {
    out.append("    let ux, uy, uz = runtime::displacement_at(");
    lits(out, at);
    out.append(");\n");
}

// Emits `alone`: false if face `fi` shares nodes with another fixed face of
// step `si` that fixes a component in `mask`.
fn emit_alone(out: u8[>..]&, si: i64, fi: i64, mask: i64) {
    out.append("        var alone = true;\n");
    for g in fixed_faces(si) {
        if same_plane(g.face, fi) || (g.mask & mask) == 0 { continue; }
        out.append("        if shared(");
        quoted(out, faces[fi].name);
        out.append(", ");
        quoted(out, faces[fi].set);
        out.append(", ");
        quoted(out, faces[g.face].name);
        out.append(", ");
        quoted(out, faces[g.face].set);
        out.append(") { alone = false; }\n");
    }
}

fn emit_report(out: u8[>..]&, si: i64, r: Report) {
    let step = steps[si].name;
    format(out, "    // report ", r.label, "\n");
    out.append("    {\n    ");
    let e = r.expr;
    let reaction = kind(e) == syntax::N_CALL;
    var value: u8[>..] = [];
    if reaction {
        let fi = find_face(name(syntax::nodes[syntax::child(e, 1)].first));
        out.append("    let x, y, z = runtime::reaction(");
        quoted(out, faces[fi].set);
        out.append(");\n");
        emit_alone(out, si, fi, fixed_mask(si, fi));
        format(value, "vector(", reaction_parts(fixed_mask(si, fi)), ", \"", system.force, "\")");
    } else if kind(e) != syntax::N_AT {
        var at = 0;
        var sum: u8[>..] = [];
        emit_point_sum(out, e, at, sum);
        format(value, "quantity(", sum, ", \"", system.length, "\")");
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
    out.append(if reaction { "        if alone { report(" } else { "        report(" });
    quoted(out, step);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", value, if reaction { "); }\n    }\n" } else { ");\n    }\n" });
}

let SIGMA_NAMES = ["sxx", "syy", "szz", "syz", "sxz", "sxy"];

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n    ");
    var got: u8[>..] = [];
    var unit = system.length;
    if x.kind == EXPECT_U {
        emit_point(out, x.at);
        format(got, "u", AXIS_NAMES[x.comp]);
    } else if x.kind == EXPECT_SUM {
        var at = 0;
        emit_point_sum(out, x.expr, at, got);
    } else if x.kind == EXPECT_SIGMA {
        out.append("    let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(");
        lits(out, x.at);
        out.append(");\n");
        got.append(SIGMA_NAMES[x.comp]);
        unit = system.pressure;
    } else {
        out.append("    let x, y, z = runtime::reaction(");
        quoted(out, faces[x.face].set);
        out.append(");\n");
        emit_alone(out, 0, x.face, 1 << x.comp);
        got.append(AXIS_NAMES[x.comp]);
        unit = system.force;
    }
    let reaction = x.kind == EXPECT_REACTION;
    out.append(if reaction { "        if alone { expect(" } else { "        expect(" });
    quoted(out, x.name);
    format(out, ", ", got, ", ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"", unit, "\", ");
    quoted(out, x.label);
    out.append(if reaction { "); }\n    }\n" } else { ");\n    }\n" });
}

// Writes the Goose program of the model that `read` read. The program
// imports the runtime from `runtime_path` and writes its VTU files under
// `results_dir`. In a study, it writes or reads the stages of `share`.
fn generate(runtime_path: u8[:], results_dir: u8[:], share: Share, out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", system.length, ", force in ",
           system.force, ".\n");
    out.append("import std;\nimport os;\nimport runtime from ");
    quoted(out, runtime_path);
    out.append(";\n");
    // solvers lane: the runtime imports the package solvers (lane_solvers.goose).
    out.append("import solvers from ");
    quoted(out, str(runtime_path, "/../solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    for si in 0..steps.len {
        if steps[si].procedure == CORE {
            emit_step(out, si, results_dir, share);
        } else {
            emit_procedure(out, steps[si].procedure, si, results_dir);
        }
    }
    out.append("fn main() {\n");
    if share.stages >= 1 && share.reads {
        emit_share(out, "", "runtime::read_mesh", share, "");
    } else {
        if mesh_lane == MESHER {
            mesher_emit(out);
        } else if mesh_lane == SHELLS {
            shells_emit_mesh(out);
        } else if !elements_emit_mesh(out) {
            emit_core_mesh(out);
        }
        if share.stages >= 1 { emit_share(out, "", "runtime::write_mesh", share, ""); }
    }
    if template == CORE {
        format(out, "    // material ", material_name, "\n    runtime::material(");
        lit(out, young);
        out.append(", ");
        lit(out, poisson);
        out.append(");\n");
        emit_laws(out);
    } else {
        emit_template(out, template);
    }
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    for s in skipped {
        out.append("    not_run(");
        quoted(out, s);
        out.append(");\n");
    }
    out.append("    print_records();\n");
    out.append("    if failed { exit(1); }\n");
    out.append("}\n");
}

// The material of each region of the laws, box by box (runtime::material_box).
// A thermal strain of the first law without a region holds in the bounds of
// the solid.
fn emit_laws(out: u8[>..]&) {
    for w, i in laws {
        if w.region < 0 && (i > 0 || w.strain == 0.0) { continue; }
        let sh = if w.region >= 0 { shapes[w.region] } else { domain };
        format(out, "    // material ", material_names[w.material], " on ", if w.region >= 0 { shape_names[w.region] } else { solid_name },
               ": E, nu and the thermal strain\n");
        let first = union_boxes.len / 6;
        let f, count = copy_boxes(if sh.kind == UNION { sh } else { Shape { kind: BOX, lo: sh.lo, hi: sh.hi, cx: 0.0, cy: 0.0, first: 0, count: 0 } });
        for b in f..f + count {
            out.append("    runtime::material_box(");
            lit(out, w.young);
            out.append(", ");
            lit(out, w.poisson);
            out.append(", ");
            lit(out, w.strain);
            for j in 0..6 {
                out.append(", ");
                lit(out, union_boxes[6 * b + j]);
            }
            out.append(");\n");
        }
        union_boxes.resize(6 * first);
    }
}

// The grid lines along axis `a` of the union of boxes of the domain: each
// face of each box, and each Point that a bc holds, is a grid plane.
// Between two planes, the grid has the fewest equal cells of at most the
// mesh size. The program meshes the same grid (runtime::mesh_union): this
// copy gives only the dof estimate.
fn grid_lines(a: i64) -> f64[>..] {
    let lo = domain.lo[a];
    let hi = domain.hi[a];
    var cuts: f64[>..] = [lo, hi];
    for b in domain.first..domain.first + domain.count {
        for c in [union_boxes[6 * b + a], union_boxes[6 * b + 3 + a]] {
            if c > lo && c < hi { cuts.push(c); }
        }
    }
    for q in point_sets {
        if q.held && q.at[a] > lo && q.at[a] < hi { cuts.push(q.at[a]); }
    }
    sort(cuts);
    let tol = 1.0e-9 * (hi - lo);
    var out: f64[>..] = [lo];
    var last = lo;
    for c in cuts {
        if c - last <= tol { continue; }
        let n = max(1, ceil((c - last) / mesh_size - 1.0e-9) as! i64);
        for k in 1..n + 1 { out.push(last + (c - last) * k / n); }
        last = c;
    }
    return out;
}

// The grid of the union of boxes, and the node set of each Point that a
// bc holds.
fn emit_grid_mesh(out: u8[>..]&) {
    format(out, "    // mesh on ", solid_name, ": 27-node hexahedra of at most ", mesh_size, " ", system.length,
           ", with a grid plane on each face of each box and through each held Point\n");
    out.append("    runtime::mesh_union([");
    for x in 6 * domain.first..6 * (domain.first + domain.count) {
        if x > 6 * domain.first { out.append(", "); }
        lit(out, union_boxes[x]);
    }
    out.append("],\n                       [");
    var points = 0;
    for q in point_sets {
        if !q.held { continue; }
        if points > 0 { out.append(", "); }
        lits(out, q.at);
        points++;
    }
    out.append("], ");
    lit(out, mesh_size);
    out.append(");\n");
    for q in point_sets {
        if !q.held { continue; }
        format(out, "    face(runtime::point_set(");
        quoted(out, faces[q.face].set);
        out.append(", ");
        lits(out, q.at);
        out.append("), ");
        quoted(out, faces[q.face].name);
        out.append(");\n");
    }
}

// The dofs of the mesh of the model that `read` read: 3 at each node of
// the lattice of 27-node hexahedra (emit_core_mesh), at most.
fn dof_estimate() -> i64 {
    var nodes = 1;
    if domain.kind == UNION {
        for j in 0..3 { nodes *= 2 * grid_lines(j).len - 1; }
        return 3 * nodes;
    }
    for j in 0..3 {
        var length = domain.hi[j] - domain.lo[j];
        if domain.kind == SECTOR && j == 1 { length *= domain.hi[0]; }
        nodes *= 2 * max(ceil(length / mesh_size - 1.0e-9) as! i64, 1) + 1;
    }
    return 3 * nodes;
}

// The box or the sector of the domain, by the lattice mesher, or the union
// of boxes by the grid mesher.
fn emit_core_mesh(out: u8[>..]&) {
    if domain.kind == UNION { emit_grid_mesh(out); return; }
    let d = domain;
    if d.kind == SECTOR {
        format(out, "    // mesh on ", solid_name, ": 27-node hexahedra of at most ", mesh_size, " ",
               system.length, ", along r, the outer arc and z\n");
        out.append("    runtime::mesh_annulus(");
        lit(out, d.cx);
        out.append(", ");
        lit(out, d.cy);
        for j in 0..3 {
            out.append(", ");
            lit(out, d.lo[j]);
            out.append(", ");
            lit(out, d.hi[j]);
        }
    } else if incompatible {
        format(out, "    // mesh on ", solid_name, ": 8-node hexahedra with incompatible modes of at most ", mesh_size,
               " ", system.length, "\n");
        out.append("    runtime::mesh_box_incompatible(");
        lits(out, d.lo);
        out.append(", ");
        lits(out, d.hi);
    } else {
        format(out, "    // mesh on ", solid_name, ": 27-node hexahedra of at most ", mesh_size, " ",
               system.length, "\n");
        out.append("    runtime::mesh_box(");
        lits(out, d.lo);
        out.append(", ");
        lits(out, d.hi);
    }
    for j in 0..3 {
        out.append(",\n                     runtime::cells_for(");
        lit(out, d.hi[j] - d.lo[j]);
        if d.kind == SECTOR && j == 1 {
            out.append(" * ");
            lit(out, d.hi[0]);
        }
        out.append(", ");
        lit(out, mesh_size);
        out.append(")");
    }
    out.append(");\n");
}
