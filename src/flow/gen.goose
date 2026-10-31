// The code generator of the flow models: a Skein syntax tree to a Goose
// program that calls src/flow/flow.goose.
//
// It takes the models with a physics `navier_stokes(f)` or `stokes(f)`, or
// the same equations written out: an unknown Vector<Velocity> and an
// unknown Pressure. The model is 2d and SI. The domain is a CSG shape of
// rect and circle, meshed with triangles by src/mesher/ with its `refine`
// on segments and circles. The steps are `static`, solved by Newton's
// method with GMRES and ILU(0). Each other construct gives an error
// `unsupported` at its token, with a fix.
//
// The scaled form of spec §21.2: the three values of the `scales` line give
// the reference length, mass and time, so each quantity q of the dimension
// L^a M^b T^c has the reference value q_ref = L^a M^b T^c, and the generated
// program gives the runtime q / q_ref. Each report, expect and VTU file is
// dimensional again: the program multiplies each result by its q_ref.
import std;
import math;
namespace flow_gen;
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

// The source of n with no spaces.
fn compact(n: i64) -> u8[>..] {
    var s: u8[>..] = [];
    for c in source_of(n) { if c != ' ' { s.push(c); } }
    return s;
}

// --- the claim ---------------------------------------------------------------

// Whether this generator takes the model: a physics that uses stokes or
// navier_stokes, or that solves a Vector<Velocity> and a Pressure.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS {
            var velocity = false;
            var pressure = false;
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_USES && (name(u) == "navier_stokes" || name(u) == "stokes") { return true; }
                if kind(u) == syntax::N_UNKNOWN {
                    let t = child(u, 0);
                    if kind(t) == syntax::N_TYPE_NAME && name(t) == "Pressure" { pressure = true; }
                    if kind(t) == syntax::N_TYPE_APPLY && name(t) == "Vector" && kind(first(t)) == syntax::N_TYPE_NAME &&
                       name(first(t)) == "Velocity" {
                        velocity = true;
                    }
                }
                u = next(u);
            }
            if velocity && pressure { return true; }
        }
        c = next(c);
    }
    return false;
}

// --- quantities and the scales (spec §21) ------------------------------------

// A scalar or a 2D vector in SI, and its dimension: the exponents of L, M
// and T in halves.
struct Q {
    n: i64,
    v: f64[2],
    d: i64[3],
}

fn real(x: f64) -> Q { Q { n: 1, v: [x, 0.0], d: [0, 0, 0] } }

// The reference length, mass and time of the scales line.
var base: f64[3] = [1.0, 1.0, 1.0];

// q_ref of the dimension d.
fn reference(d: i64[3]) -> f64 {
    return pow(base[0], d[0] / 2.0) * pow(base[1], d[1] / 2.0) * pow(base[2], d[2] / 2.0);
}

fn scaled(q: Q) -> Q {
    var s = q;
    let r = reference(q.d);
    for j in 2 { s.v[j] = q.v[j] / r; }
    return s;
}

fn same_dim(a: i64[3], b: i64[3]) -> bool { a[0] == b[0] && a[1] == b[1] && a[2] == b[2] }

// The table row that a table column stands for, in an expect of a table.
var table_node = -1;
var table_row = -1;

recursive fn quantity(n: i64) -> Q {
    let k = kind(n);
    if k == syntax::N_NUMBER {
        if syntax::nodes[n].aux >= 0 { unsupported(tok(n), "a flow model has no temperatures", "remove the degC value"); }
        var q = real(codegen::number(tok(n)));
        if first(n) >= 0 {
            let u = codegen::unit_value(first(n));
            q.v[0] *= u.f;
            q.d = [u.l, u.m, u.t];
        }
        return q;
    }
    if k == syntax::N_NAME { return name_quantity(n); }
    if k == syntax::N_MEMBER && kind(first(n)) == syntax::N_NAME && table_node >= 0 && name(first(n)) == name(table_node) {
        return table_quantity(n);
    }
    if k == syntax::N_VECTOR { return vector_quantity(n); }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        var q = quantity(first(n));
        for j in 2 { q.v[j] = -q.v[j]; }
        return q;
    }
    if k == syntax::N_BINARY { return binary_quantity(n); }
    unsupported(tok(n), str("the expression '", source_of(n), "' is not a value of this generator"),
                "write a number with a unit, a name, a vector (x, y), or + - * / ^ of these");
    return real(0.0);
}

recursive fn name_quantity(n: i64) -> Q {
    var v = codegen::find_value(name(n));
    if v < 0 { v = property(name(n)); }
    if v < 0 && name(n) == "pi" { return real(PI); }
    if v < 0 || codegen::depth > 100 {
        unsupported(tok(n), str("'", name(n), "' is not a value that this generator knows"),
                    "declare it before, with let or param and a number");
        return real(0.0);
    }
    codegen::depth++;
    let q = quantity(v);
    codegen::depth--;
    return q;
}

// A vector: its dimension is that of its components, with the unit after
// them. A 0 with no unit fits each dimension.
recursive fn vector_quantity(n: i64) -> Q {
    var q = Q { n: 0, v: [0.0, 0.0], d: [0, 0, 0] };
    var f = 1.0;
    var unit_d: i64[3] = [0, 0, 0];
    var have = false;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNIT {
            let u = codegen::unit_value(c);
            f = u.f;
            unit_d = [u.l, u.m, u.t];
        } else {
            let x = quantity(c);
            if x.n != 1 { unsupported(tok(c), "a component of a vector is a vector", "write a scalar"); }
            if q.n < 2 { q.v[q.n] = x.v[0]; }
            if x.v[0] != 0.0 || !same_dim(x.d, [0, 0, 0]) {
                if have && !same_dim(x.d, q.d) {
                    unsupported(tok(c), "the components of the vector have different dimensions", "give each component the same unit");
                }
                q.d = x.d;
                have = true;
            }
            q.n++;
        }
        c = next(c);
    }
    if q.n != 2 { unsupported(tok(n), "a vector of this 2d model has 2 components", "write (x, y)"); }
    for j in 2 { q.v[j] *= f; }
    for j in 3 { q.d[j] += unit_d[j]; }
    return q;
}

recursive fn binary_quantity(n: i64) -> Q {
    let op = syntax::text(tok(n));
    let a = quantity(child(n, 0));
    let b = quantity(child(n, 1));
    if op == "+" || op == "-" {
        var r = if a.v[0] == 0.0 && a.v[1] == 0.0 && a.n == 1 { b } else { a };
        for j in 2 { r.v[j] = if op == "+" { a.v[j] + b.v[j] } else { a.v[j] - b.v[j] }; }
        return r;
    }
    if (op == "*" || op == "/") && b.n == 2 && (op == "/" || a.n == 2) {
        unsupported(tok(n), "a product of two vectors, or a scalar divided by a vector", "multiply or divide a vector by a scalar");
        return a;
    }
    if op == "*" || op == "/" {
        var r = if a.n == 2 { a } else { b };
        let s = if a.n == 2 { b.v[0] } else { a.v[0] };
        if op == "*" {
            for j in 2 { r.v[j] = if a.n == 2 { a.v[j] * s } else { s * b.v[j] }; }
            r.d = [a.d[0] + b.d[0], a.d[1] + b.d[1], a.d[2] + b.d[2]];
        } else {
            r = a;
            for j in 2 { r.v[j] = a.v[j] / b.v[0]; }
            r.d = [a.d[0] - b.d[0], a.d[1] - b.d[1], a.d[2] - b.d[2]];
        }
        return r;
    }
    if op == "^" && a.n == 1 && b.n == 1 && b.d == [0, 0, 0] && round(2.0 * b.v[0]) == 2.0 * b.v[0] {
        let h = round(2.0 * b.v[0]) as! i64;
        return Q { n: 1, v: [pow(a.v[0], b.v[0]), 0.0], d: [a.d[0] * h / 2, a.d[1] * h / 2, a.d[2] * h / 2] };
    }
    unsupported(tok(n), str("the operator '", op, "' is not in this generator"), "use + - * / and ^ with a number");
    return a;
}

// The value of column n.tok of the row table_row of table table_node.
fn table_quantity(n: i64) -> Q {
    var col = 0;
    var c = first(table_node);
    while c >= 0 && kind(c) == syntax::N_COLUMN {
        if name(c) == name(n) {
            var r = c;
            for i in table_row + 1 {
                r = next(r);
                while r >= 0 && kind(r) != syntax::N_ROW { r = next(r); }
            }
            var q = quantity(child(r, col));
            if child(c, 1) >= 0 {
                let u = codegen::unit_value(child(c, 1));
                q.v[0] *= u.f;
                q.d = [u.l, u.m, u.t];
            }
            return q;
        }
        col++;
        c = next(c);
    }
    unsupported(tok(n), str("the table has no column '", name(n), "'"), "name a column of the table");
    return real(0.0);
}

fn scalar_value(n: i64) -> f64 {
    let q = quantity(n);
    if q.n != 1 { unsupported(tok(n), "the value is a vector", "write a scalar"); }
    return q.v[0];
}

// The scaled value of n.
fn scaled_value(n: i64) -> f64 { return scaled(quantity(n)).v[0]; }

// `scales a, U, rho`: the reference length, mass and time from the three
// scales, by the solution of the 3 × 3 system of the logarithms.
fn scales(n: i64) {
    var m: f64[9] = [0.0; 9];
    var b: f64[3] = [0.0; 3];
    var count = 0;
    var c = first(n);
    while c >= 0 {
        let q = quantity(c);
        if count < 3 && q.n == 1 && q.v[0] > 0.0 {
            for j in 3 { m[3 * count + j] = q.d[j] / 2.0; }
            b[count] = log(q.v[0]);
        }
        if q.n != 1 || q.v[0] <= 0.0 { unsupported(tok(c), "a scale is a positive scalar", "write a length, a velocity and a density"); }
        count++;
        c = next(c);
    }
    let x, solved = solve3(m, b);
    if count != 3 || !solved {
        unsupported(tok(n), "the scales of a flow model are three values that give a length, a mass and a time",
                    "write scales L, U, rho with a length, a velocity and a density");
        return;
    }
    for j in 3 { base[j] = exp(x[j]); }
}

// The solution of m x = b, by Gauss elimination with partial pivots.
fn solve3(m: f64[9], b: f64[3]) -> f64[3], bool {
    var a = m;
    var r = b;
    for k in 3 {
        var p = k;
        for i in k + 1..3 { if abs(a[3 * i + k]) > abs(a[3 * p + k]) { p = i; } }
        if abs(a[3 * p + k]) < 1.0e-12 { return [0.0; 3], false; }
        for j in 3 {
            let t = a[3 * k + j];
            a[3 * k + j] = a[3 * p + j];
            a[3 * p + j] = t;
        }
        let t = r[k];
        r[k] = r[p];
        r[p] = t;
        for i in k + 1..3 {
            let f = a[3 * i + k] / a[3 * k + k];
            for j in 3 { a[3 * i + j] -= f * a[3 * k + j]; }
            r[i] -= f * r[k];
        }
    }
    var x: f64[3] = [0.0; 3];
    for kk in 3 {
        let k = 2 - kk;
        var s = r[k];
        for j in k + 1..3 { s -= a[3 * k + j] * x[j]; }
        x[k] = s / a[3 * k + k];
    }
    return x, true;
}

// --- the model -------------------------------------------------------------------

// A `let` of a shape (rect, circle, or - | & of shapes), a boundary
// region `boundary(shape, condition)`, or a point region `point(p)`.
let SHAPE = 0;
let BOUNDARY = 1;
let POINT = 2;

struct Named {
    name: u8[..64],
    kind: i64,
    value: i64,     // the expression node
}

// A `refine`: a segment from (p0, p1) to (p2, p3), or a circle of centre
// (p0, p1) and radius p2, scaled, with the scaled size h.
struct Refine {
    circle: bool,
    p: f64[4],
    h: f64,
}

let BC_V = 0;       // component comp of the velocity, or both for comp -1
let BC_P = 1;

// A `bc`: the region (an index in `named`), the target and the scaled value.
struct Bc {
    region: i64,
    kind: i64,
    comp: i64,
    v: f64[2],
}

struct Report {
    expr: i64,
    label: u8[..256],
}

struct Step {
    name: u8[..64],
    first_own: i64,     // the bcs of the step itself: a range of `bcs`
    own_count: i64,
    first_bc: i64,      // all bcs that hold in the step: a range of `active`
    bc_count: i64,
    first_report: i64,
    report_count: i64,
    rtol: f64,
    output: bool,
}

// An `expect`: for each row of a table, when its sides name a table.
struct Expect {
    node: i64,
    table: i64,     // the N_TABLE node, or -1
    step: i64,      // the step that checks it at its end, or -1 for an item
}

var named: Named[>..];
var tables: i64[>..];
var refines: Refine[>..];
var bcs: Bc[>..];
var active: Bc[>..];
var global_bcs: i64[>..];
var reports: Report[>..];
var steps: Step[>..];
var expects: Expect[>..];
var material_node = -1;
var material_name: u8[..64];
var physics_name: u8[..64];
var physics_tok = -1;
var domain = -1;            // the shape of the physics: an index in `named`
var velocity: u8[..64];     // the names of the unknowns
var pressure: u8[..64];
var stress: u8[..64];       // the let of the stress: sigma
var convective = true;
var force: f64[2] = [0.0, 0.0];
var force_node = -1;
var mesh_size = 0.0;
var grade = 1.3;
var order = 2;
var stabilize = 1;          // flow::SUPG
var stabilize_set = false;
var has_scales = false;

let STABILIZE_NAMES = ["none", "supg", "supg_pspg", "supg_pspg_lsic"];

fn reset() {
    codegen::reset();
    named.clear();
    tables.clear();
    refines.clear();
    bcs.clear();
    active.clear();
    global_bcs.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    material_node = -1;
    material_name = [];
    physics_name = [];
    physics_tok = -1;
    domain = -1;
    velocity = [];
    pressure = [];
    set64(stress, "sigma");
    convective = true;
    force = [0.0, 0.0];
    force_node = -1;
    mesh_size = 0.0;
    grade = 1.3;
    order = 2;
    stabilize = 1;
    stabilize_set = false;
    has_scales = false;
    base = [1.0, 1.0, 1.0];
    table_node = -1;
    table_row = -1;
    used_regions.clear();
    vector_comp = -1;
}

fn find_named(word: u8[:]) -> i64 {
    for g, i in named {
        if g.name == word { return i; }
    }
    return -1;
}

// The named shape, boundary or point of node n, or -1.
fn named_of(n: i64, k: i64) -> i64 {
    guard kind(n) == syntax::N_NAME else { return -1; }
    let i = find_named(name(n));
    return if i >= 0 && named[i].kind == k { i } else { -1 };
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

fn header(n: i64) {
    set64(codegen::model_name, name(n));
    if syntax::text(tok(child(n, 0))) != "2d" {
        unsupported(tok(child(n, 0)), "the flow generator has 2d models", "write 2d in the header");
    }
    if syntax::text(tok(child(n, 1))) != "SI" {
        unsupported(tok(child(n, 1)), "the flow generator has SI models", "write SI in the header");
    }
    codegen::system = codegen::SYSTEMS[0];
}

recursive fn is_shape(n: i64) -> bool {
    if is_call(n, "rect") || is_call(n, "circle") { return true; }
    if kind(n) == syntax::N_NAME { return named_of(n, SHAPE) >= 0; }
    let op = if kind(n) == syntax::N_BINARY { syntax::text(tok(n)) } else { "" };
    return (op == "-" || op == "|" || op == "&") && is_shape(child(n, 0)) && is_shape(child(n, 1));
}

// `let` and `param`: a shape, a region, or a value.
fn value_item(n: i64) {
    let v = child(n, 1);
    var g = Named { name: [], kind: SHAPE, value: v };
    set64(g.name, name(n));
    if is_shape(v) {
        named.push(g);
    } else if is_call(v, "boundary") || is_call(v, "point") {
        g.kind = if is_call(v, "point") { POINT } else { BOUNDARY };
        named.push(g);
    } else {
        codegen::value_names.push(g.name);
        codegen::value_nodes.push(v);
    }
}

fn material(n: i64) {
    if material_node >= 0 {
        unsupported(tok(n), "a second material: the flow generator has one", "remove the second material");
        return;
    }
    material_node = n;
    set64(material_name, name(n));
    var p = first(n);
    while p >= 0 {
        if name(p) != "rho" && name(p) != "mu" {
            unsupported(tok(p), str("the property '", name(p), "' is not in the flow generator"), "give the fluid rho and mu");
        }
        p = next(p);
    }
}

let PHYSICS_FIX = "write physics flow on dom: navier_stokes solves v: Vector<Velocity>, p: Pressure { material fluid }";

// The statements of a physics written out: `let sigma = -p*I + 2*mu*sym(grad(v))`,
// `eq rho*grad(v)*v - div(sigma) = 0` and `eq div(v) = 0`.
fn written_out(c: i64, eqs: i64&) {
    let k = kind(c);
    let v = velocity;
    let p = pressure;
    if k == syntax::N_LET {
        if compact(child(c, 1)) != str("-", p, "*I+2*mu*sym(grad(", v, "))") {
            unsupported(tok(c), str("the let '", name(c), "' is not the Newtonian stress"),
                        str("write let sigma: Tensor2<Pressure> = -", p, "*I + 2*mu*sym(grad(", v, "))"));
            return;
        }
        set64(stress, name(c));
        return;
    }
    if k != syntax::N_EQN {
        unsupported(tok(c), "the physics statement is not in the flow generator", PHYSICS_FIX);
        return;
    }
    let left = compact(child(c, 0));
    let right = compact(child(c, 1));
    let s = stress;
    eqs++;
    if left == str("div(", v, ")") && right == "0" { return; }
    if left == str("rho*grad(", v, ")*", v, "-div(", s, ")") && right == "0" { return; }
    unsupported(tok(c), "the eq is not an equation of the flow generator",
                str("write eq rho*grad(", v, ")*", v, " - div(", s, ") = 0 and eq div(", v, ") = 0"));
}

fn physics_item(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: the flow generator has one", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    var template = false;
    var eqs = 0;
    var statements: i64[>..] = [];
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            domain = named_of(first(c), SHAPE);
            if domain < 0 { unsupported(tok(first(c)), "the region of the physics is not a shape that a let names", "write physics flow on dom, with let dom: Surface = rect(...)"); }
        } else if k == syntax::N_USES {
            template = true;
            convective = name(c) == "navier_stokes";
            if name(c) != "navier_stokes" && name(c) != "stokes" {
                unsupported(tok(c), "the flow generator has the templates navier_stokes(f) and stokes(f)", PHYSICS_FIX);
            }
            if first(c) >= 0 { force_node = first(first(c)); }
        } else if k == syntax::N_UNKNOWN {
            let t = child(c, 0);
            if kind(t) == syntax::N_TYPE_NAME && name(t) == "Pressure" && pressure.len == 0 {
                set64(pressure, name(c));
            } else if kind(t) == syntax::N_TYPE_APPLY && velocity.len == 0 {
                set64(velocity, name(c));
            } else {
                unsupported(tok(c), "the flow generator solves one velocity and one pressure", PHYSICS_FIX);
            }
        } else if k == syntax::N_MATASSIGN {
            if syntax::child_count(c) != 1 || first(first(c)) >= 0 {
                unsupported(tok(c), "the flow generator gives one material to the whole physics", PHYSICS_FIX);
            }
        } else if k == syntax::N_ASSUME {
            assumptions(c);
        } else {
            statements.push(c);
        }
        c = next(c);
    }
    for s in statements { written_out(s, eqs); }
    if !template && eqs != 2 { unsupported(tok(n), "a written-out flow physics needs two eqs", PHYSICS_FIX); }
    if template && statements.len > 0 { unsupported(tok(statements[0]), "a physics with a template has no let or eq", "remove it"); }
    if domain < 0 && ok() { unsupported(tok(n), "the physics has no region", "write physics flow on dom"); }
}

fn assumptions(n: i64) {
    var a = first(n);
    while a >= 0 {
        let w = name(a);
        if w != "incompressible" && w != "newtonian" && w != "laminar" {
            unsupported(tok(a), str("the assumption '", w, "' is not in the flow generator"), "use incompressible, newtonian and laminar");
        }
        a = next(a);
    }
}

fn method(n: i64) {
    if name(n) != physics_name || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "the flow generator has only the cg method of the physics", "write method flow: cg { order: { v: 2, p: 1 } }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            order_setting(v);
        } else if kind(s) == syntax::N_SETTING && name(s) == "stabilize" {
            let i = if kind(v) == syntax::N_NAME { position(STABILIZE_NAMES, name(v)) } else { -1 };
            if i < 0 { unsupported(tok(v), "the flow generator has stabilize: none, supg, supg_pspg or supg_pspg_lsic", "write stabilize: supg_pspg"); }
            stabilize = max(i, 0);
            stabilize_set = true;
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in the flow generator", "use only mesh, order and stabilize");
        }
        s = next(s);
    }
}

// `order: 2` or `order: { v: 2, p: 1 }`: the velocity has order 1 or 2,
// and the pressure has order 1.
fn order_setting(v: i64) {
    if kind(v) != syntax::N_RECORD {
        order = scalar_value(v) as! i64;
        if order != 1 && order != 2 { unsupported(tok(v), "the flow generator has order 1 and 2", "write order: { v: 2, p: 1 }"); }
        return;
    }
    var a = first(v);
    while a >= 0 {
        let w = syntax::text(syntax::nodes[a].aux);
        let k = scalar_value(first(a)) as! i64;
        if w == velocity && (k == 1 || k == 2) {
            order = k;
        } else if w != pressure || k != 1 {
            unsupported(tok(a), "the flow generator has a velocity of order 1 or 2 and a pressure of order 1", "write order: { v: 2, p: 1 }");
        }
        a = next(a);
    }
}

// --- the mesh ----------------------------------------------------------------------

let REFINE_FIX = "refine a boundary on a side of a rect, boundary(box, y == 1 m), or on a circle, boundary(dom, hole)";

// The refines of boundary g with the scaled size h: a side x == c or
// y == c of the rect of the boundary, an || of these, or a circle.
recursive fn refine_parts(shape: i64, cond: i64, h: f64) {
    if kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "||" {
        refine_parts(shape, child(cond, 0), h);
        refine_parts(shape, child(cond, 1), h);
        return;
    }
    let circle = named_of(cond, SHAPE);
    if circle >= 0 && is_call(named[circle].value, "circle") {
        let c = named[circle].value;
        let centre = scaled(quantity(arg(c, 0)));
        refines.push(Refine { circle: true, p: [centre.v[0], centre.v[1], scaled_value(arg(c, 1)), 0.0], h: h });
        return;
    }
    let side = kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "==" && kind(child(cond, 0)) == syntax::N_NAME &&
        (name(child(cond, 0)) == "x" || name(child(cond, 0)) == "y");
    var rect = if kind(shape) == syntax::N_NAME && named_of(shape, SHAPE) >= 0 { named[named_of(shape, SHAPE)].value } else { shape };
    if !side || !is_call(rect, "rect") {
        unsupported(tok(cond), "the refine is not on a side of a rect or on a circle", REFINE_FIX);
        return;
    }
    let a = scaled(quantity(arg(rect, 0)));
    let b = scaled(quantity(arg(rect, 1)));
    let c = scaled_value(child(cond, 1));
    if name(child(cond, 0)) == "x" {
        refines.push(Refine { circle: false, p: [c, min(a.v[1], b.v[1]), c, max(a.v[1], b.v[1])], h: h });
    } else {
        refines.push(Refine { circle: false, p: [min(a.v[0], b.v[0]), c, max(a.v[0], b.v[0]), c], h: h });
    }
}

fn mesh(n: i64) {
    if mesh_size > 0.0 {
        unsupported(tok(n), "a second mesh: the flow generator has one", "remove the second mesh");
        return;
    }
    let region = child(n, 0);
    if domain < 0 || kind(region) != syntax::N_NAME || name(region) != named[domain].name {
        unsupported(tok(region), "the mesh is not on the region of the physics", "write mesh m on dom, with the region of the physics");
    }
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_REFINE {
            let r = named_of(v, BOUNDARY);
            if r < 0 {
                unsupported(tok(v), "the refine is not on a boundary that a let names", REFINE_FIX);
            } else {
                let b = named[r].value;
                refine_parts(arg(b, 0), arg(b, 1), scaled_value(child(s, 1)));
            }
        } else if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "layers and crack_tip are not in the flow generator", "remove the line");
        } else if name(s) == "size" {
            mesh_size = scaled_value(v);
        } else if name(s) == "grade" {
            grade = scalar_value(v);
            if grade <= 1.0 { unsupported(tok(v), "the grade is more than 1", "write grade: 1.3"); }
        } else if name(s) == "cell" {
            if kind(v) != syntax::N_NAME || name(v) != "tri" { unsupported(tok(v), "the flow generator has tri cells", "write cell: tri"); }
        } else {
            unsupported(tok(s), str("the mesh setting '", name(s), "' is not in the flow generator"), "use only size, cell, grade and refine");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 1 cm"); }
}

// --- the steps ---------------------------------------------------------------------

let BC_FIX_TEXT = "write bc wall: v = 0, bc lid: v = (1, 0) m/s, bc side: v.y = 0, bc corner: p = 0, or bc outlet: sigma*n = 0";

// A `bc` item or statement. It gives its index in `bcs`, or -1.
fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let target = child(n, 1);
    let v = child(n, 2);
    let r = if kind(region) == syntax::N_NAME { find_named(name(region)) } else { -1 };
    if r < 0 || named[r].kind == SHAPE {
        unsupported(tok(region), "the region of the bc is not a boundary or a point that a let names", BC_FIX_TEXT);
        return -1;
    }
    var b = Bc { region: r, kind: BC_V, comp: -1, v: [0.0, 0.0] };
    let t = compact(target);
    let q = scaled(quantity(v));
    if t == velocity {
        if q.n == 1 && q.v[0] != 0.0 { unsupported(tok(v), str(velocity, " is a vector"), "write v = 0, or a vector"); }
        b.v = q.v;
    } else if t == str(velocity, ".x") || t == str(velocity, ".y") {
        b.comp = if t[t.len - 1] == 'x' { 0 } else { 1 };
        b.v[0] = q.v[0];
    } else if t == pressure {
        b.kind = BC_P;
        b.v[0] = q.v[0];
    } else if t == str(stress, "*n") {
        // The natural condition of the weak form.
        if q.v[0] != 0.0 || q.v[1] != 0.0 { unsupported(tok(v), "the flow generator has the traction 0 only", "write sigma*n = 0"); }
        return -1;
    } else {
        unsupported(tok(target), str("the bc target '", t, "' is not in the flow generator"), BC_FIX_TEXT);
        return -1;
    }
    use_region(r);
    bcs.push(b);
    return bcs.len - 1;
}

fn solve_settings(n: i64, s: Step&) {
    var c = first(first(n));
    while c >= 0 {
        let v = first(c);
        let w = if kind(v) == syntax::N_NAME { name(v) } else { "" };
        if name(c) == "rtol" {
            s.rtol = scalar_value(v);
            if s.rtol <= 0.0 { unsupported(tok(v), "rtol is more than 0", "write rtol: 1e-8"); }
        } else if !(name(c) == "nonlinear" && w == "newton" || name(c) == "linear" && w == "gmres" || name(c) == "precond" && w == "ilu0") {
            unsupported(tok(c), "the flow generator solves with nonlinear: newton, linear: gmres and precond: ilu0",
                        "write solve { nonlinear: newton, linear: gmres, precond: ilu0, rtol: 1e-8 }");
        }
        c = next(c);
    }
}

fn step(n: i64) {
    let proc = child(n, 0);
    if name(proc) != "static" || first(proc) >= 0 {
        unsupported(tok(proc), "the flow generator has only static steps with no arguments", "write step steady: static { ... }");
    }
    var s = Step { name: [], first_own: bcs.len, own_count: 0, first_bc: 0, bc_count: 0,
                   first_report: reports.len, report_count: 0, rtol: 1.0e-8, output: false };
    set64(s.name, name(n));
    var c = next(proc);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            bc(c);
        } else if k == syntax::N_SOLVE {
            solve_settings(c, s);
        } else if k == syntax::N_OUTPUT {
            s.output = true;
            output_fields(c);
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                report(r);
                r = next(r);
            }
        } else if k == syntax::N_EXPECT {
            expect(c, steps.len);
        } else {
            unsupported(tok(c), "the step statement is not in the flow generator", "use only bc, solve, output, report and expect");
        }
        c = next(c);
    }
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

fn output_fields(n: i64) {
    var o = first(n);
    while o >= 0 {
        if kind(o) == syntax::N_EVERY {
            unsupported(tok(o), "a static step has no every", "remove it");
        } else if kind(o) != syntax::N_NAME || (name(o) != velocity && name(o) != pressure) {
            unsupported(tok(o), "the flow generator writes the velocity and the pressure", str("write output ", velocity, ", ", pressure));
        }
        o = next(o);
    }
}

let REPORT_FIX = "write v.x at (x, y), p at (x, y), integral((sigma*n).x, wall), or + - * / of these and values";

fn report(item: i64) {
    let e = first(item);
    let good = if whole_vector(e) { point_ok(e) } else { result_ok(e) };
    if !good { return; }
    var r = Report { expr: e, label: [] };
    codegen::set256(r.label, source_of(item));
    reports.push(r);
}

// `f at p` of the velocity, a component of it, or the pressure: the
// flow::VALUE_* of each component, with -1 for the second of a scalar.
fn point_field(e: i64) -> i64, i64 {
    let f = compact(child(e, 0));
    if f == velocity { return 0, 1; }
    if f == str(velocity, ".x") { return 0, -1; }
    if f == str(velocity, ".y") { return 1, -1; }
    if f == pressure { return 2, -1; }
    return -1, -1;
}

// Whether `e` is `v at p`, the velocity vector at a point. Only a report
// shows it: an expect or a + - * / uses one component.
fn whole_vector(e: i64) -> bool {
    guard kind(e) == syntax::N_AT else { return false; }
    let what, other = point_field(e);
    return other >= 0;
}

// The point of `f at p`: a vector, or a Point that a let names.
fn at_point(p: i64) -> Q {
    let i = named_of(p, POINT);
    return quantity(if i >= 0 { arg(named[i].value, 0) } else { p });
}

// Whether the point of `f at p` is a vector (x, y), with an error if not.
fn point_ok(e: i64) -> bool {
    let before = syntax::errors.len;
    let p = at_point(child(e, 1));
    if p.n != 2 && syntax::errors.len == before {
        unsupported(tok(child(e, 1)), "the point of 'at' is not a vector (x, y)", "write at (x, y) m, or the name of a point((x, y) m)");
    }
    return syntax::errors.len == before;
}

// `integral((sigma*n).x, wall)`: the component, or -1.
fn traction_integral(e: i64) -> i64 {
    guard is_call(e, "integral") && syntax::child_count(e) == 3 && named_of(arg(e, 1), BOUNDARY) >= 0 else { return -1; }
    let f = compact(arg(e, 0));
    if f == str("(", stress, "*n).x") { return 0; }
    if f == str("(", stress, "*n).y") { return 1; }
    return -1;
}

// Whether `e` is a result of this generator, with an error if not.
recursive fn result_ok(e: i64) -> bool {
    let k = kind(e);
    if k == syntax::N_AT {
        let what, other = point_field(e);
        if what < 0 { unsupported(tok(e), str("'", source_of(child(e, 0)), "' at a point is not in the flow generator"), REPORT_FIX); }
        if other >= 0 { unsupported(tok(e), str("the vector '", source_of(child(e, 0)), "' at a point is only in a report, not in an expect or a + - * /"), REPORT_FIX); }
        return what >= 0 && other < 0 && point_ok(e);
    }
    if is_call(e, "integral") {
        if traction_integral(e) >= 0 { use_region(named_of(arg(e, 1), BOUNDARY)); }
        if traction_integral(e) < 0 { unsupported(tok(e), "the flow generator integrates (sigma*n).x or (sigma*n).y on a boundary", REPORT_FIX); }
        return traction_integral(e) >= 0;
    }
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(e))) >= 0 {
        return result_ok(child(e, 0)) && result_ok(child(e, 1));
    }
    if k == syntax::N_UNARY && syntax::text(tok(e)) == "-" { return result_ok(first(e)); }
    let q = quantity(e);
    return ok() && q.n >= 1;
}

// The table that the sides of expect n name, or -1.
recursive fn table_of(n: i64) -> i64 {
    if kind(n) == syntax::N_MEMBER && kind(first(n)) == syntax::N_NAME {
        for t in tables { if name(t) == name(first(n)) { return t; } }
    }
    var c = first(n);
    while c >= 0 {
        let t = table_of(c);
        if t >= 0 { return t; }
        c = next(c);
    }
    return -1;
}

// An `expect` item (step -1), or an `expect` of step `step`.
fn expect(n: i64, step: i64) {
    var x = Expect { node: n, table: table_of(n), step: step };
    table_node = x.table;
    table_row = 0;
    let good = result_ok(child(n, 0));
    if good { quantity(child(n, 1)); }
    if good && kind(child(n, 2)) != syntax::N_PERCENT { quantity(child(n, 2)); }
    table_node = -1;
    if good && ok() { expects.push(x); }
}

// The rows of table t.
fn row_count(t: i64) -> i64 {
    var k = 0;
    var c = first(t);
    while c >= 0 {
        if kind(c) == syntax::N_ROW { k++; }
        c = next(c);
    }
    return k;
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
            physics_node = c;
        } else if k == syntax::N_MESH {
            mesh_node = c;
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k == syntax::N_TABLE {
            tables.push(c);
        } else if k == syntax::N_SCALES {
            has_scales = true;
        } else if k != syntax::N_USE && k != syntax::N_BC && k != syntax::N_STEP && k != syntax::N_EXPECT && k != syntax::N_STUDY {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in the flow generator"),
                        "remove it: the flow generator runs one steady flow model");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if !has_scales { unsupported(at, "the model has no scales line", "add scales L, U, rho"); }
    c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_SCALES { scales(c); }
        c = next(c);
    }
    if physics_node >= 0 { physics_item(physics_node); } else { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if !ok() { return false; }
    read_fluid(at);
    if mesh_node < 0 { unsupported(at, "the model has no mesh", "add mesh m on dom { size: 1 cm }"); } else { mesh(mesh_node); }
    for m in methods { method(m); }
    if !stabilize_set { stabilize = if convective { 1 } else { 0 }; }
    if order == 1 && stabilize < 2 { unsupported(at, "order 1 for the velocity needs pressure stabilization", "write stabilize: supg_pspg"); }
    if !ok() { return false; }
    read_steps(prog, at);
    return ok();
}

// rho and mu of the material, and the body force f of the template.
fn read_fluid(at: i64) {
    let r = property("rho");
    let m = property("mu");
    if material_node < 0 || r < 0 || m < 0 {
        unsupported(at, "the flow needs a material with rho and mu", "add material fluid { rho: Density = 1000 kg/m^3, mu: DynViscosity = 1e-3 Pa*s }");
        return;
    }
    if scalar_value(r) <= 0.0 || scalar_value(m) <= 0.0 { unsupported(tok(material_node), "rho and mu are more than 0", "give rho > 0 and mu > 0"); }
    if force_node >= 0 {
        let f = scaled(quantity(force_node));
        if f.n != 2 && !(f.n == 1 && f.v[0] == 0.0) { unsupported(tok(force_node), "the body force f is a vector", "write navier_stokes((0, -9.81) m/s^2*rho)"); }
        if f.n == 2 { force = f.v; }
    }
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
            expect(c, -1);
        }
        c = next(c);
    }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step steady: static { bc wall: v = 0 ... }"); }
    if ok() { step_bcs(); }
    // The conditions of the regions give their errors here, not in generate.
    var scratch: u8[>..] = [];
    emit_selects(scratch);
}

// Puts `b` in `set`, in place of a bc on the same region and target.
fn merge(set: Bc[>..]&, b: Bc) {
    for i in 0..set.len {
        if set[i].region == b.region && set[i].kind == b.kind && set[i].comp == b.comp {
            set[i] = b;
            return;
        }
    }
    set.push(b);
}

// The bcs that hold in each step (spec §12.2).
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

// --- the mesh file -----------------------------------------------------------------

// The mesher shape of the shape expression n, scaled.
recursive fn mesher_shape(n: i64) -> i64 {
    if kind(n) == syntax::N_NAME { return mesher_shape(named[named_of(n, SHAPE)].value); }
    if is_call(n, "rect") {
        let a = scaled(quantity(arg(n, 0)));
        let b = scaled(quantity(arg(n, 1)));
        return mesher::rect(a.v[0], a.v[1], b.v[0], b.v[1]);
    }
    if is_call(n, "circle") {
        let c = scaled(quantity(arg(n, 0)));
        return mesher::circle(c.v[0], c.v[1], scaled_value(arg(n, 1)));
    }
    let a = mesher_shape(child(n, 0));
    let b = mesher_shape(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { mesher::subtract(a, b) } else if op == "|" { mesher::union(a, b) } else { mesher::intersect(a, b) };
}

// Meshes the domain of the model that `read` read, and writes its
// vertices and triangles in the directory `dir` for flow::read_mesh. It
// gives an error, and false, if the mesher does not finish or a file
// cannot be written.
fn write_mesh(dir: u8[:]) -> bool {
    mesher::clear_shapes();
    mesher::clear_refines();
    mesher::set_grade(grade);
    for r in refines {
        if r.circle { mesher::refine_circle(r.p[0], r.p[1], r.p[2], r.h); } else { mesher::refine_segment(r.p[0], r.p[1], r.p[2], r.p[3], r.h); }
    }
    if !mesher::mesh_2d(mesher_shape(named[domain].value), mesh_size) {
        write_stderr(str("skein run: the mesher cannot mesh '", named[domain].name, "' with the size and the refines of the mesh\n"));
        return false;
    }
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    for n in mesher::point_count() {
        let p = mesher::vertex(n);
        xy.push(p[0]);
        xy.push(p[1]);
    }
    for c in mesher::cell_count() {
        for k in 3 { tri.push(mesher::cell_node(c, k)); }
    }
    cells = mesher::cell_count();
    vertices = mesher::point_count();
    if !write_file(str(dir, "/flow.xy"), to_bytes(xy[..])) || !write_file(str(dir, "/flow.tri"), to_bytes(tri[..])) {
        write_stderr(str("skein run: cannot write the mesh in '", dir, "'\n"));
        return false;
    }
    return true;
}

var cells = 0;
var vertices = 0;

// The dofs of the mesh, from the area of the domain and the sizes of the
// mesh: the mesher makes about 2 / (√3 h²/2) cells for each unit of area.
fn dof_estimate() -> i64 {
    let area = shape_area(named[domain].value);
    var cell_count = area / (0.4330127 * mesh_size * mesh_size);
    for r in refines {
        let length = if r.circle { TAU * r.p[2] } else { sqrt((r.p[2] - r.p[0]) * (r.p[2] - r.p[0]) + (r.p[3] - r.p[1]) * (r.p[3] - r.p[1])) };
        // The band where the size grows from h to the mesh size.
        let band = (mesh_size - r.h) / (grade - 1.0);
        cell_count += length * band / (0.4330127 * r.h * mesh_size);
    }
    let vertex_count = 0.5 * cell_count;
    return (if order == 1 { 3.0 * vertex_count } else { 3.0 * vertex_count + 2.0 * 1.5 * cell_count }) as! i64;
}

recursive fn shape_area(n: i64) -> f64 {
    if kind(n) == syntax::N_NAME { return shape_area(named[named_of(n, SHAPE)].value); }
    if is_call(n, "rect") {
        let a = scaled(quantity(arg(n, 0)));
        let b = scaled(quantity(arg(n, 1)));
        return abs((b.v[0] - a.v[0]) * (b.v[1] - a.v[1]));
    }
    if is_call(n, "circle") { return PI * scaled_value(arg(n, 1)) * scaled_value(arg(n, 1)); }
    let a = shape_area(child(n, 0));
    let b = shape_area(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { max(a - b, 0.0) } else if op == "|" { a + b } else { min(a, b) };
}

// --- the emitters ------------------------------------------------------------------

// The SI unit of the dimension d (in halves of L, M and T).
fn unit_of(d: i64[3]) -> u8[>..] {
    if same_dim(d, [0, 0, 0]) { return ""; }
    if same_dim(d, [2, 0, -2]) { return "m/s"; }
    if same_dim(d, [-2, 2, -4]) { return "Pa"; }
    if same_dim(d, [0, 2, -4]) { return "N/m"; }
    if same_dim(d, [2, 2, -4]) { return "N"; }
    var s: u8[>..] = [];
    let names = ["m", "kg", "s"];
    for j in 3 {
        if d[j] == 0 { continue; }
        if s.len > 0 { s.push('*'); }
        s.append(names[j]);
        if d[j] != 2 {
            s.push('^');
            if d[j] % 2 == 0 { format(s, d[j] / 2); } else { format(s, d[j] / 2.0); }
        }
    }
    return s;
}

// The component of a velocity vector in a report, or -1.
var vector_comp = -1;

let VALUE_NAMES = ["flow::VALUE_VX", "flow::VALUE_VY", "flow::VALUE_P"];
let VELOCITY_D: i64[3] = [2, 0, -2];
let PRESSURE_D: i64[3] = [-2, 2, -4];

// The Goose expression of the SI value of result e, and its dimension.
// The program evaluates it in a function with the variable `name`, the
// label of the report or expect.
recursive fn emit_result(out: u8[>..]&, e: i64) -> i64[3] {
    let k = kind(e);
    if k == syntax::N_AT {
        var what, other = point_field(e);
        if other >= 0 && vector_comp >= 0 { what = vector_comp; }
        let at = scaled(at_point(child(e, 1)));
        let d = if what == 2 { PRESSURE_D } else { VELOCITY_D };
        format(out, "at(", VALUE_NAMES[what], ", ");
        lit(out, at.v[0]);
        out.append(", ");
        lit(out, at.v[1]);
        out.append(", name) * ");
        lit(out, reference(d));
        return d;
    }
    if is_call(e, "integral") {
        // The force of the wall on the fluid: p_ref L for each unit of depth.
        let r = named_of(arg(e, 1), BOUNDARY);
        format(out, "wall_force(select_", named[r].name, "(), ", traction_integral(e), ", ");
        quoted(out, named[r].name);
        out.append(") * ");
        lit(out, reference([0, 2, -4]));
        return [0, 2, -4];
    }
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(e))) >= 0 {
        let op = syntax::text(tok(e));
        out.push('(');
        let a = emit_result(out, child(e, 0));
        format(out, " ", op, " ");
        let b = emit_result(out, child(e, 1));
        out.push(')');
        if op == "*" { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
        if op == "/" { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
        return a;
    }
    if k == syntax::N_UNARY {
        out.append("(-");
        let a = emit_result(out, first(e));
        out.push(')');
        return a;
    }
    let q = quantity(e);
    out.push('(');
    lit(out, q.v[0]);
    out.push(')');
    return q.d;
}

// The Goose condition of the region of a `boundary(shape, cond)`, in the
// coordinates x and y and the tolerance tol of the program.
fn emit_boundary(out: u8[>..]&, b: i64) {
    let shape = arg(b, 0);
    let s = named_of(shape, SHAPE);
    if s != domain {
        if s < 0 || !emit_shape_edge(out, named[s].value) {
            unsupported(tok(shape), "the shape of the boundary is not the domain, a rect or a circle", "write boundary(dom, ...)");
        }
        out.append(" && ");
    }
    out.push('(');
    emit_cond(out, arg(b, 1));
    out.push(')');
}

// The condition of the edge of a rect or a circle. False for another shape.
fn emit_shape_edge(out: u8[>..]&, n: i64) -> bool {
    if is_call(n, "circle") {
        let c = scaled(quantity(arg(n, 0)));
        out.append("abs(sqrt((x - ");
        lit(out, c.v[0]);
        out.append(") * (x - ");
        lit(out, c.v[0]);
        out.append(") + (y - ");
        lit(out, c.v[1]);
        out.append(") * (y - ");
        lit(out, c.v[1]);
        out.append(")) - ");
        lit(out, scaled_value(arg(n, 1)));
        out.append(") <= tol");
        return true;
    }
    guard is_call(n, "rect") else { return false; }
    let a = scaled(quantity(arg(n, 0)));
    let b = scaled(quantity(arg(n, 1)));
    out.append("on_rect(x, y, ");
    for j in 2 {
        lit(out, min(a.v[j], b.v[j]));
        out.append(", ");
    }
    lit(out, max(a.v[0], b.v[0]));
    out.append(", ");
    lit(out, max(a.v[1], b.v[1]));
    out.append(", tol)");
    return true;
}

let COND_FIX = "write a condition of x and y with == < <= > >= && || and abs, or the name of a rect or a circle";

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
    if k == syntax::N_NAME && named_of(n, SHAPE) >= 0 {
        if !emit_shape_edge(out, named[named_of(n, SHAPE)].value) { unsupported(tok(n), "the shape is not a rect or a circle", COND_FIX); }
        return;
    }
    let i = position(["==", "!=", "<", "<=", ">", ">="], op);
    if i < 0 {
        unsupported(tok(n), str("the condition '", source_of(n), "' is not in the flow generator"), COND_FIX);
        return;
    }
    // Each comparison has the tolerance of spec §8.2.
    out.append(if i <= 1 { "abs(" } else { "(" });
    emit_coordinate(out, child(n, 0));
    out.append(" - ");
    emit_coordinate(out, child(n, 1));
    out.append([") <= tol", ") > tol", ") < -tol", ") <= tol", ") > tol", ") >= -tol"][i]);
}

recursive fn emit_coordinate(out: u8[>..]&, n: i64) {
    let k = kind(n);
    if k == syntax::N_NAME && (name(n) == "x" || name(n) == "y") {
        out.append(name(n));
        return;
    }
    if is_call(n, "abs") {
        out.append("abs(");
        emit_coordinate(out, arg(n, 0));
        out.push(')');
        return;
    }
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(n))) >= 0 && mentions_xy(n) {
        out.push('(');
        emit_coordinate(out, child(n, 0));
        format(out, " ", syntax::text(tok(n)), " ");
        emit_coordinate(out, child(n, 1));
        out.push(')');
        return;
    }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" && mentions_xy(n) {
        out.append("(-");
        emit_coordinate(out, first(n));
        out.push(')');
        return;
    }
    out.push('(');
    lit(out, scaled_value(n));
    out.push(')');
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

    // A bc on a region that holds no node is a fault of the model.
    fn region(found: bool, name: u8[:]) {
        if !found {
            write_stderr(str("skein run: no node of the mesh is on the region '", name, "'\n"));
            exit(1);
        }
    }

    // On the edge of the rect from (x0, y0) to (x1, y1).
    fn on_rect(x: f64, y: f64, x0: f64, y0: f64, x1: f64, y1: f64, t: f64) -> bool {
        let inside = x >= x0 - t && x <= x1 + t && y >= y0 - t && y <= y1 + t;
        return inside && (abs(x - x0) <= t || abs(x - x1) <= t || abs(y - y0) <= t || abs(y - y1) <= t);
    }

    // The scaled value `what` at (x, y). A point outside the mesh is a
    // fault of the model.
    fn at(what: i64, x: f64, y: f64, name: u8[:]) -> f64 {
        let v, found = flow::value_at(what, x, y);
        if !found {
            write_stderr(str("skein run: the point of '", name, "' is not in the mesh\n"));
            failed = true;
        }
        return v;
    }

    // The scaled force of the wall on the fluid, component `comp`.
    fn wall_force(found: bool, comp: i64, name: u8[:]) -> f64 {
        region(found, name);
        return flow::region_force(comp);
    }

    fn solve(step: u8[:], rtol: f64) {
        let converged = flow::solve_steady(rtol);
        format(solves, "  ", step, ",", flow::cell_count(), ",", flow::dof_count(), ",", flow::newton_steps, ",",
               flow::linear_iterations, ",", flow::last_seconds, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        let r = flow::last_residual / max(flow::first_residual, 1.0e-300);
        check("solve", "residual", converged, str(r), str("<= ", rtol));
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let t = if relative { within * abs(want) } else { within };
        check("post", "expect", abs(got - want) <= t, str(name, ": ", quantity(got, unit)), limit);
    }

    fn output(path: u8[:], x_ref: f64, v_ref: f64, p_ref: f64) {
        if !flow::write_vtu(path, x_ref, v_ref, p_ref) {
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
        format(s, "solve[", solve_count, "]{step,cells,dofs,newton_steps,gmres_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        write_stdout(s);
    }

    """;

// The regions that the program selects: indices in `named`.
var used_regions: i64[>..];

fn use_region(r: i64) { if position(used_regions, r) < 0 { used_regions.push(r); } }

// `fn select_<name>() -> bool` for each region that the program selects.
fn emit_selects(out: u8[>..]&) {
    for r in used_regions {
        let g = named[r];
        format(out, "fn select_", g.name, "() -> bool {\n    return ");
        if g.kind == POINT {
            let p = scaled(quantity(arg(g.value, 0)));
            out.append("flow::select_point(");
            lit(out, p.v[0]);
            out.append(", ");
            lit(out, p.v[1]);
            out.append(");\n}\n\n");
        } else {
            out.append("flow::select_boundary() { x, y => ");
            emit_boundary(out, g.value);
            out.append(" };\n}\n\n");
        }
    }
}

fn emit_bc(out: u8[>..]&, b: Bc) {
    let g = named[b.region];
    format(out, "    region(select_", g.name, "(), ");
    quoted(out, g.name);
    out.append(");\n");
    if b.kind == BC_P {
        out.append("    flow::fix_pressure(");
        lit(out, b.v[0]);
        out.append(");\n");
        return;
    }
    for j in 2 {
        if b.comp >= 0 && b.comp != j { continue; }
        format(out, "    flow::fix_velocity(", j, ", ");
        lit(out, if b.comp >= 0 { b.v[0] } else { b.v[j] });
        out.append(");\n");
    }
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": static\nfn step_", s.name, "() {\n    flow::clear_bcs();\n");
    for b in active[s.first_bc..s.first_bc + s.bc_count] { emit_bc(out, b); }
    out.append("    solve(");
    quoted(out, s.name);
    out.append(", ");
    lit(out, s.rtol);
    out.append(");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s.name, r); }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(", ");
        lit(out, reference([2, 0, 0]));
        out.append(", ");
        lit(out, reference(VELOCITY_D));
        out.append(", ");
        lit(out, reference(PRESSURE_D));
        out.append(");\n");
    }
    out.append("}\n\n");
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n        let name = ");
    quoted(out, r.label);
    out.append(";\n        report(");
    quoted(out, step);
    out.append(", name, ");
    let e = r.expr;
    if whole_vector(e) {
        // The velocity vector at a point.
        var x: u8[>..] = [];
        vector_comp = 0;
        let d = emit_result(x, e);
        var y: u8[>..] = [];
        vector_comp = 1;
        emit_result(y, e);
        vector_comp = -1;
        format(out, "str(\"(\", ", x, ", \", \", ", y, ", \") ", unit_of(d), "\")");
    } else {
        var x: u8[>..] = [];
        let d = emit_result(x, e);
        format(out, "quantity(", x, ", \"", unit_of(d), "\")");
    }
    out.append(");\n    }\n");
}

// One check for an expect, or one for each row of its table.
fn emit_expect(out: u8[>..]&, x: Expect) {
    let n = x.node;
    let rows = if x.table >= 0 { row_count(x.table) } else { 1 };
    let within = child(n, 2);
    let relative = kind(within) == syntax::N_PERCENT;
    var tolerance = 0.0;
    if relative {
        let t = syntax::text(tok(within));
        let v, good = parse_flt(t[..t.len - 1]);
        tolerance = v / 100.0;
    } else {
        tolerance = quantity(within).v[0];
    }
    table_node = x.table;
    for row in rows {
        table_row = row;
        var id: u8[>..] = [];
        id.append(name(n));
        if x.table >= 0 { format(id, "[", row + 1, "]"); }
        var got: u8[>..] = [];
        let d = emit_result(got, child(n, 0));
        let want = quantity(child(n, 1)).v[0];
        format(out, "    // expect ", id, "\n    {\n        let name = ");
        quoted(out, id);
        format(out, ";\n        expect(name, ", got, ", ");
        lit(out, want);
        out.append(", ");
        lit(out, tolerance);
        format(out, ", ", if relative { "true" } else { "false" }, ", \"", unit_of(d), "\", ");
        quoted(out, str(want, " ", unit_of(d), " within ", source_of(within)));
        out.append(");\n    }\n");
    }
    table_node = -1;
}

// Writes the Goose program of the model that `read` read. The program
// imports the flow package from `flow_path`, reads the mesh that
// `write_mesh` wrote in `results_dir`, and writes its VTU files there.
fn generate(flow_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    out.append("// The flow package works in the scaled units of the scales line (spec §21.2): length ");
    lit(out, base[0]);
    out.append(" m, mass ");
    lit(out, base[1]);
    out.append(" kg, time ");
    lit(out, base[2]);
    out.append(" s.\nimport std;\nimport os;\nimport flow from ");
    quoted(out, flow_path);
    // The flow package imports the package solvers, beside it.
    out.append(";\nimport solvers from ");
    quoted(out, str(flow_path[..flow_path.len - 4], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    emit_selects(out);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n");
    format(out, "    // mesh on ", named[domain].name, ": P", order, "/P1 triangles, written by skein run\n    if !flow::read_mesh(");
    quoted(out, results_dir);
    format(out, ", ", order, ") {\n        write_stderr(\"skein run: cannot read the mesh in ");
    out.append(results_dir);
    out.append("\\n\");\n        exit(3);\n    }\n    tol = flow::tolerance();\n");
    format(out, "    // physics ", physics_name, ", material ", material_name, "\n    flow::fluid(");
    lit(out, scaled_value(property("rho")));
    out.append(", ");
    lit(out, scaled_value(property("mu")));
    format(out, ", ", if convective { "true" } else { "false" }, ", flow::", ["NONE", "SUPG", "SUPG_PSPG", "SUPG_PSPG_LSIC"][stabilize], ");\n");
    if force[0] != 0.0 || force[1] != 0.0 {
        out.append("    flow::body_force(");
        lit(out, force[0]);
        out.append(", ");
        lit(out, force[1]);
        out.append(");\n");
    }
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
