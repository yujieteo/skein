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

fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }

fn compact(n: i64) -> u8[>..] {
    var s: u8[>..] = [];
    for c in source_of(n) { if c != ' ' { s.push(c); } }
    return s;
}

fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS {
            var velocity = false;
            var pressure = false;
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_USES && (name(u) == "navier_stokes" || name(u) == "stokes" || name(u) == "k_omega_sst") { return true; }
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

struct Q {
    n: i64,
    v: f64[2],
    d: i64[3],
}

fn real(x: f64) -> Q { Q { n: 1, v: [x, 0.0], d: [0, 0, 0] } }

var base: f64[3] = [1.0, 1.0, 1.0];

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
    if k == syntax::N_CALL { return call_quantity(n); }
    unsupported(tok(n), str("the expression '", source_of(n), "' is not a value of this generator"),
                "write a number with a unit, a name, a vector (x, y), or + - * / ^ of these");
    return real(0.0);
}

var fn_names: u8[..64][>..];
var fn_nodes: i64[>..];

fn find_fn(word: u8[:]) -> i64 {
    for w, i in fn_names { if w == word { return fn_nodes[i]; } }
    return -1;
}

var bind_names: u8[..64][>..];
var bind_vals: Q[>..];

recursive fn call_quantity(n: i64) -> Q {
    let callee = child(n, 0);
    if kind(callee) != syntax::N_NAME {
        unsupported(tok(n), "the call is not of a name", "call a function by its name");
        return real(0.0);
    }
    let fi = position(["sin", "cos", "tan", "sqrt", "exp", "ln", "abs"], name(callee));
    if fi >= 0 && syntax::child_count(n) == 2 {
        let a = quantity(arg(n, 0));
        if a.n != 1 { unsupported(tok(n), "a function of a vector", "write a scalar"); return real(0.0); }
        let x = a.v[0];
        if fi == 6 { return Q { n: 1, v: [abs(x), 0.0], d: [a.d[0], a.d[1], a.d[2]] }; }
        if fi == 3 {
            if a.d[0] % 2 != 0 || a.d[1] % 2 != 0 || a.d[2] % 2 != 0 {
                unsupported(tok(n), "sqrt of a quantity with an odd dimension exponent", "take the square root of a quantity with even exponents, or divide by a unit first");
                return real(0.0);
            }
            return Q { n: 1, v: [sqrt(x), 0.0], d: [a.d[0] / 2, a.d[1] / 2, a.d[2] / 2] };
        }
        if a.d[0] != 0 || a.d[1] != 0 || a.d[2] != 0 {
            unsupported(tok(n), str(name(callee), " of a quantity with a dimension"), "divide the argument by a unit to make it dimensionless");
            return real(0.0);
        }
        return real(if fi == 0 { sin(x) } else if fi == 1 { cos(x) } else if fi == 2 { tan(x) } else if fi == 4 { exp(x) } else if fi == 5 { log(x) } else { abs(x) });
    }
    var f = find_fn(name(callee));
    if f < 0 || codegen::depth > 100 {
        unsupported(tok(n), str("the function '", name(callee), "' is not a fn of the model"), "call sin, cos, sqrt, or a fn that the model declares");
        return real(0.0);
    }
    var count = 0;
    var body = -1;
    var c = first(f);
    while c >= 0 {
        if kind(c) == syntax::N_FNPARAM { count++; }
        body = c;
        c = next(c);
    }
    if count != syntax::child_count(n) - 1 || count > 8 {
        unsupported(tok(n), "the call does not give one value for each parameter", "give each parameter of the function");
        return real(0.0);
    }
    var vals = default<Q[8]>();
    for k in count { vals[k] = quantity(arg(n, k)); }
    let base_len = bind_names.len;
    var k = 0;
    c = first(f);
    while c >= 0 {
        if kind(c) == syntax::N_FNPARAM {
            var w: u8[..64] = [];
            set64(w, name(c));
            bind_names.push(w);
            bind_vals.push(vals[k]);
            k++;
        }
        c = next(c);
    }
    codegen::depth++;
    let q = quantity(body);
    codegen::depth--;
    while bind_names.len > base_len {
        bind_names.pop();
        bind_vals.pop();
    }
    return q;
}

recursive fn name_quantity(n: i64) -> Q {
    var b = bind_names.len;
    while b > 0 {
        b--;
        if bind_names[b] == name(n) { return bind_vals[b]; }
    }
    if probe_on && thermal && name(n) == temp { return real(probe_value); }
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

fn scaled_value(n: i64) -> f64 { return scaled(quantity(n)).v[0]; }

fn scales(n: i64) {
    var m: f64[9] = [0.0; 9];
    var b: f64[3] = [0.0; 3];
    var count = 0;
    var c = first(n);
    while c >= 0 {
        let q = quantity(c);
        if count >= 3 && thermal_scale_ok(c) { count++; c = next(c); continue; }
        if count < 3 && q.n == 1 && q.v[0] > 0.0 {
            for j in 3 { m[3 * count + j] = q.d[j] / 2.0; }
            b[count] = log(q.v[0]);
        }
        if q.n != 1 || q.v[0] <= 0.0 { unsupported(tok(c), "a scale is a positive scalar", "write a length, a velocity and a density"); }
        count++;
        c = next(c);
    }
    let x, solved = solve3(m, b);
    if count < 3 || count > 4 || !solved {
        unsupported(tok(n), "the scales of a flow model are three values that give a length, a mass and a time",
                    "write scales L, U, rho with a length, a velocity and a density");
        return;
    }
    for j in 3 { base[j] = exp(x[j]); }
}

fn thermal_scale_ok(c: i64) -> bool {
    let q = quantity(c);
    return q.n == 1 && q.v[0] > 0.0 && same_dim(q.d, [0, 0, 0]);
}

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

let SHAPE = 0;
let BOUNDARY = 1;
let POINT = 2;
let LINE = 3;

struct Named {
    name: u8[..64],
    kind: i64,
    value: i64,
}

struct Refine {
    circle: bool,
    p: f64[4],
    h: f64,
}

let BC_V = 0;
let BC_P = 1;
let BC_U = 2;
let BC_T = 3;
let BC_FLUX = 4;
let BC_NT = 5;
let BC_D = 6;

struct Bc {
    region: i64,
    kind: i64,
    comp: i64,
    v: f64[2],
    expr: i64,
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
    rtol: f64,
    output: bool,
    solid_output: bool,
}

struct Expect {
    node: i64,
    table: i64,
    step: i64,
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
var domain = -1;
var velocity: u8[..64];
var pressure: u8[..64];
var stress: u8[..64];
var convective = true;
var turbulent = false;
var spalart = false;
var sa_nt: u8[..64];
var sa_d: u8[..64];
var sa_d_name: u8[..64];
var sa_d_tok = -1;
var sa_lets = 0;
var sa_eqs = 0;
var v_start: f64[2] = [0.0, 0.0];
var nt_start = 0.0;
var layer_region = -1;
var layer_first = 0.0;
var layer_ratio = 0.0;
var layer_count = 0;
var turb_k: u8[..64];
var turb_w: u8[..64];
var period: f64[2] = [0.0, 0.0];
var period_high = -1;
var force: f64[2] = [0.0, 0.0];
var force_node = -1;
var mesh_size = 0.0;
var grade = 1.3;
var order = 2;
var stabilize = 1;
var stabilize_set = false;
var has_scales = false;
var solid_tok = -1;
var solid_name: u8[..64];
var solid_domain = -1;
var displacement: u8[..64];
var solid_material = -1;
var solid_size = 0.0;
var interface_node = -1;

var thermal = false;
var temp: u8[..64];
var temp_start = 0.0;
var flux: u8[..64];
var energy_seen = false;
var heat_source = 0.0;
var buoyancy: f64[2] = [0.0, 0.0];
var probe_on = false;
var probe_value = 0.0;
var th_tok = -1;
var th_name: u8[..64];
var th_domain = -1;
var th_temp: u8[..64];
var th_flux: u8[..64];
var th_mat = -1;
var th_source = 0.0;
var th_size = 0.0;
var th_energy_seen = false;
var p_offset = 0.0;
var p_offset_set = false;
var thermal_materials: i64[>..];

let STABILIZE_NAMES = ["none", "supg", "supg_pspg", "supg_pspg_lsic"];

fn reset() {
    codegen::reset();
    named.clear();
    fn_names.clear();
    fn_nodes.clear();
    bind_names.clear();
    bind_vals.clear();
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
    turbulent = false;
    spalart = false;
    sa_nt = [];
    sa_d = [];
    sa_d_name = [];
    sa_d_tok = -1;
    sa_lets = 0;
    sa_eqs = 0;
    v_start = [0.0, 0.0];
    nt_start = 0.0;
    layer_region = -1;
    layer_first = 0.0;
    layer_ratio = 0.0;
    layer_count = 0;
    turb_k = [];
    turb_w = [];
    period = [0.0, 0.0];
    period_high = -1;
    force = [0.0, 0.0];
    force_node = -1;
    mesh_size = 0.0;
    grade = 1.3;
    order = 2;
    stabilize = 1;
    stabilize_set = false;
    has_scales = false;
    solid_tok = -1;
    solid_name = [];
    solid_finite = false;
    solid_stress = [];
    solid_domain = -1;
    displacement = [];
    solid_material = -1;
    solid_size = 0.0;
    interface_node = -1;
    thermal = false;
    temp = [];
    temp_start = 0.0;
    flux = [];
    energy_seen = false;
    heat_source = 0.0;
    buoyancy = [0.0, 0.0];
    probe_on = false;
    th_tok = -1;
    th_name = [];
    th_domain = -1;
    th_temp = [];
    th_flux = [];
    th_mat = -1;
    th_source = 0.0;
    th_size = 0.0;
    th_energy_seen = false;
    p_offset = 0.0;
    p_offset_set = false;
    thermal_materials.clear();
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

fn named_of(n: i64, k: i64) -> i64 {
    guard kind(n) == syntax::N_NAME else { return -1; }
    let i = find_named(name(n));
    return if i >= 0 && named[i].kind == k { i } else { -1 };
}

fn property(word: u8[:]) -> i64 { return property_of(material_node, word); }

fn property_of(m: i64, word: u8[:]) -> i64 {
    guard m >= 0 else { return -1; }
    var p = first(m);
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

fn is_fill(n: i64) -> bool { return is_call(n, "fill") && syntax::child_count(n) == 2 && is_call(arg(n, 0), "spline"); }

recursive fn is_shape(n: i64) -> bool {
    if is_call(n, "rect") || is_call(n, "circle") || is_fill(n) { return true; }
    if kind(n) == syntax::N_NAME { return named_of(n, SHAPE) >= 0; }
    let op = if kind(n) == syntax::N_BINARY { syntax::text(tok(n)) } else { "" };
    return (op == "-" || op == "|" || op == "&") && is_shape(child(n, 0)) && is_shape(child(n, 1));
}

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
    } else if is_call(v, "boundary") || is_call(v, "point") || is_call(v, "line") || is_union(v) {
        g.kind = if is_call(v, "point") { POINT } else if is_call(v, "line") { LINE } else { BOUNDARY };
        named.push(g);
    } else {
        codegen::value_names.push(g.name);
        codegen::value_nodes.push(v);
    }
}

fn material(n: i64) {
    if property_of(n, "E") >= 0 {
        solid_materials(n);
        return;
    }
    if property_of(n, "mu") < 0 && property_of(n, "k") >= 0 && property_of(n, "rho") < 0 {
        thermal_materials.push(n);
        var p = first(n);
        while p >= 0 {
            if name(p) != "k" && name(p) != "cp" && name(p) != "rho" { unsupported(tok(p), str("the property '", name(p), "' is not in the flow generator"), "give the solid k"); }
            p = next(p);
        }
        return;
    }
    if material_node >= 0 {
        unsupported(tok(n), "a second material: the flow generator has one", "remove the second material");
        return;
    }
    material_node = n;
    set64(material_name, name(n));
    var p = first(n);
    while p >= 0 {
        if name(p) != "rho" && name(p) != "mu" && name(p) != "nu" && name(p) != "cp" && name(p) != "k" {
            unsupported(tok(p), str("the property '", name(p), "' is not in the flow generator"), "give the fluid rho and mu or nu, and cp and k for a temperature");
        }
        p = next(p);
    }
}

fn solid_materials(n: i64) {
    if solid_material >= 0 {
        unsupported(tok(n), "a second solid material: the flow generator has one", "remove the second material");
        return;
    }
    solid_material = n;
    var p = first(n);
    while p >= 0 {
        if name(p) != "E" && name(p) != "nu" { unsupported(tok(p), str("the property '", name(p), "' is not in the flow generator"), "give the solid E and nu"); }
        p = next(p);
    }
    if property_of(n, "nu") < 0 { unsupported(tok(n), "the solid needs nu", "add nu: Real = 0.3"); }
}

let SA_FIX = "write the Spalart-Allmaras model with the lets chi, fv1, fv2, nut, omg, s_t, r_sa, g_sa, fw, sigma and q_nt and the 3 eqs, as in examples/femahoy/naca0012_ladson.skn";

fn spalart_let_text(word: u8[:]) -> u8[>..] {
    let v = velocity;
    let nt = sa_nt;
    let d = sa_d;
    if word == "chi" { return str(nt, "/nu"); }
    if word == "fv1" { return "chi^3/(chi^3+cv1^3)"; }
    if word == "fv2" { return "1-chi/(1+chi*fv1)"; }
    if word == "nut" { return str(nt, "*fv1"); }
    if word == "omg" { return str("sqrt(2*ddot(skew(grad(", v, ")),skew(grad(", v, "))))"); }
    if word == "s_t" { return str("max(omg+fv2*", nt, "/(kap^2*", d, "^2),c_s*omg)"); }
    if word == "r_sa" { return str("min(", nt, "/(s_t*kap^2*", d, "^2),10)"); }
    if word == "g_sa" { return "r_sa+cw2*(r_sa^6-r_sa)"; }
    if word == "fw" { return "g_sa*exp(ln((1+cw3^6)/(g_sa^6+cw3^6))/6)"; }
    if word == "sigma" { return str("-", pressure, "*I+2*rho*(nu+nut)*sym(grad(", v, "))"); }
    if word == "q_nt" { return str("-(nu+", nt, ")/sig_sa*grad(", nt, ")"); }
    return "";
}

fn spalart_statement(c: i64) {
    let v = velocity;
    let nt = sa_nt;
    if kind(c) == syntax::N_LET {
        let want = spalart_let_text(name(c));
        if want.len == 0 || compact(child(c, 1)) != want {
            unsupported(tok(c), str("the let '", name(c), "' is not in the Spalart-Allmaras template"), SA_FIX);
            return;
        }
        if name(c) == "sigma" { set64(stress, name(c)); }
        sa_lets++;
        return;
    }
    let left = compact(child(c, 0));
    let right = compact(child(c, 1));
    sa_eqs++;
    if left == str("rho*grad(", v, ")*", v, "-div(sigma)") && right == "0" { return; }
    if left == str("div(", v, ")") && right == "0" { return; }
    if left == str("dot(", v, ",grad(", nt, "))+div(q_nt)-cb2/sig_sa*dot(grad(", nt, "),grad(", nt, "))") &&
       right == str("cb1*s_t*", nt, "-cw1*fw*(", nt, "/", sa_d, ")^2") {
        return;
    }
    unsupported(tok(c), "the eq is not in the Spalart-Allmaras template", SA_FIX);
}

fn is_walldist(n: i64) -> bool {
    var found = false;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNKNOWN {
            let t = child(c, 0);
            if kind(t) == syntax::N_TYPE_NAME && name(t) == "Length" { found = true; } else { return false; }
        }
        c = next(c);
    }
    return found;
}

fn walldist_item(n: i64) {
    if sa_d_tok >= 0 {
        unsupported(tok(n), "a second wall distance: the flow generator has one", "remove the second physics");
        return;
    }
    sa_d_tok = tok(n);
    set64(sa_d_name, name(n));
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_UNKNOWN { set64(sa_d, name(c)); }
        c = next(c);
    }
    c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_EQN {
            if compact(child(c, 0)) != str("norm(grad(", sa_d, "))") || compact(child(c, 1)) != "1" {
                unsupported(tok(c), "the eq of the wall distance is not the eikonal equation", str("write eq norm(grad(", sa_d, ")) = 1"));
            }
        } else if kind(c) != syntax::N_UNKNOWN && kind(c) != syntax::N_ON && kind(c) != syntax::N_MATASSIGN {
            unsupported(tok(c), "the physics of a wall distance has one eq", str("write eq norm(grad(", sa_d, ")) = 1"));
        }
        c = next(c);
    }
}

let PHYSICS_FIX = "write physics flow on dom: navier_stokes solves v: Vector<Velocity>, p: Pressure { material fluid }";

fn written_out(c: i64, eqs: i64&) {
    let k = kind(c);
    let v = velocity;
    let p = pressure;
    if spalart && (k == syntax::N_LET || k == syntax::N_EQN) {
        spalart_statement(c);
        return;
    }
    if k == syntax::N_LET && thermal && compact(child(c, 1)) == str("fourier(k,", temp, ")") {
        set64(flux, name(c));
        return;
    }
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
    if left == str("rho*grad(", v, ")*", v, "-div(", s, ")") && thermal {
        body_force_of(child(c, 1));
        return;
    }
    if thermal && flux.len > 0 && left == str("rho*cp*dot(", v, ",grad(", temp, "))+div(", flux, ")") {
        energy_seen = true;
        heat_source = if right == "0" { 0.0 } else { scaled_value(child(c, 1)) };
        return;
    }
    unsupported(tok(c), "the eq is not an equation of the flow generator",
                str("write eq rho*grad(", v, ")*", v, " - div(", s, ") = 0 and eq div(", v, ") = 0"));
}

fn body_force_of(f: i64) {
    probe_on = true;
    var q: Q[3] = [real(0.0); 3];
    for i in 3 {
        probe_value = i;
        q[i] = scaled(quantity(f));
    }
    probe_on = false;
    if !ok() { return; }
    for i in 3 {
        if q[i].n != 2 && !(q[i].n == 1 && q[i].v[0] == 0.0) {
            unsupported(tok(f), "the body force is not a vector", "write a vector, for example rho*(1 - beta*(T - T0))*(0, -g)");
            return;
        }
    }
    var f0: f64[2] = [0.0; 2];
    var f1: f64[2] = [0.0; 2];
    var f2: f64[2] = [0.0; 2];
    if q[0].n == 2 { f0 = q[0].v; }
    if q[1].n == 2 { f1 = q[1].v; }
    if q[2].n == 2 { f2 = q[2].v; }
    for j in 2 {
        let scale = max(abs(f0[j]), max(abs(f1[j]), abs(f2[j])));
        if abs(f2[j] - 2.0 * f1[j] + f0[j]) > 1.0e-9 * scale {
            unsupported(tok(f), str("the body force is not affine in ", temp), "write a force of the form f0 + f1*T, for example rho*(1 - beta*(T - T0))*(0, -g)");
            return;
        }
        force[j] = f0[j];
        buoyancy[j] = f1[j] - f0[j];
    }
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
            convective = name(c) != "stokes";
            turbulent = name(c) == "k_omega_sst";
            if name(c) != "navier_stokes" && name(c) != "stokes" && !turbulent {
                unsupported(tok(c), "the flow generator has the templates navier_stokes(f), stokes(f) and k_omega_sst(f)", PHYSICS_FIX);
            }
            if first(c) >= 0 { force_node = first(first(c)); }
        } else if k == syntax::N_UNKNOWN {
            let t = child(c, 0);
            if kind(t) == syntax::N_TYPE_NAME && name(t) == "Pressure" && pressure.len == 0 {
                set64(pressure, name(c));
            } else if kind(t) == syntax::N_TYPE_NAME && name(t) == "Temperature" && !thermal {
                thermal = true;
                set64(temp, name(c));
                if syntax::child_count(c) > 1 { temp_start = scaled_value(child(c, 1)); }
            } else if kind(t) == syntax::N_TYPE_NAME && name(t) == "KinViscosity" && !spalart {
                spalart = true;
                turbulent = true;
                set64(sa_nt, name(c));
                if syntax::child_count(c) > 1 { nt_start = scaled_value(child(c, 1)); }
            } else if kind(t) == syntax::N_TYPE_APPLY && velocity.len == 0 {
                set64(velocity, name(c));
                if syntax::child_count(c) > 1 {
                    let q = scaled(quantity(child(c, 1)));
                    if q.n != 2 { unsupported(tok(c), str(name(c), " starts at a vector"), "write v: Vector<Velocity> = U*e_d"); }
                    v_start = q.v;
                }
            } else if turbulent && turb_k.len == 0 {
                set64(turb_k, name(c));
            } else if turbulent && turb_w.len == 0 {
                set64(turb_w, name(c));
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
    if spalart && ok() && (sa_lets != 11 || sa_eqs != 3) {
        unsupported(tok(n), "the Spalart-Allmaras physics has the 11 lets and the 3 eqs of the template", SA_FIX);
    }
    if spalart { eqs = 2; }
    if !template && eqs != (if thermal { 3 } else { 2 }) { unsupported(tok(n), str("a written-out flow physics needs ", if thermal { "three" } else { "two" }, " eqs"), PHYSICS_FIX); }
    if thermal && ok() && !energy_seen {
        unsupported(tok(n), "the temperature has no energy equation", str("write let q: Vector<HeatFlux> = fourier(k, ", temp, ") and eq rho*cp*dot(", velocity, ", grad(", temp, ")) + div(q) = 0"));
    }
    if template && statements.len > 0 { unsupported(tok(statements[0]), "a physics with a template has no let or eq", "remove it"); }
    if domain < 0 && ok() { unsupported(tok(n), "the physics has no region", "write physics flow on dom"); }
}

let SOLID_FIX = "write physics flag on bar: linear_elastic solves u: Vector<Length> { material steel  assume plane_strain }";

fn is_thermal_solid(n: i64) -> bool {
    var temperature = false;
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_USES { return false; }
        if kind(c) == syntax::N_UNKNOWN {
            let t = child(c, 0);
            if !(kind(t) == syntax::N_TYPE_NAME && name(t) == "Temperature") { return false; }
            temperature = true;
        }
        c = next(c);
    }
    return temperature;
}

let THERMAL_SOLID_FIX = "write physics solid on wall solves Ts: Temperature = T0 { material pvc  let q: Vector<HeatFlux> = fourier(k, Ts)  eq div(q) = 0 }";

fn thermal_solid_item(n: i64) {
    if th_tok >= 0 {
        unsupported(tok(n), "a second solid with heat conduction: the flow generator has one", "remove the second physics");
        return;
    }
    th_tok = tok(n);
    set64(th_name, name(n));
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            th_domain = named_of(first(c), SHAPE);
            if th_domain < 0 { unsupported(tok(first(c)), "the region of the physics is not a shape that a let names", THERMAL_SOLID_FIX); }
        } else if k == syntax::N_UNKNOWN {
            if th_temp.len > 0 { unsupported(tok(c), "the solid solves one temperature", THERMAL_SOLID_FIX); }
            set64(th_temp, name(c));
        } else if k == syntax::N_MATASSIGN {
            for m in thermal_materials { if first(c) >= 0 && name(m) == name(first(c)) { th_mat = m; } }
        } else if k == syntax::N_LET && compact(child(c, 1)) == str("fourier(k,", th_temp, ")") {
            set64(th_flux, name(c));
        } else if k == syntax::N_EQN && th_flux.len > 0 && compact(child(c, 0)) == str("div(", th_flux, ")") {
            th_energy_seen = true;
            th_source = if compact(child(c, 1)) == "0" { 0.0 } else { scaled_value(child(c, 1)) };
        } else if k != syntax::N_ASSUME {
            unsupported(tok(c), "the statement is not in the solid with heat conduction", THERMAL_SOLID_FIX);
        }
        c = next(c);
    }
    if th_mat < 0 && ok() { unsupported(tok(n), "the solid needs a material with k", "add material pvc { k: Conductivity = 0.25 W/(m*K) }"); }
    if !th_energy_seen && ok() { unsupported(tok(n), "the solid has no energy equation", THERMAL_SOLID_FIX); }
}

fn assumes_finite(c: i64) -> bool {
    var a = first(c);
    while a >= 0 {
        if name(a) == "finite_strain" { return true; }
        a = next(a);
    }
    return false;
}

fn is_solid(n: i64) -> bool {
    var c = first(n);
    while c >= 0 {
        if kind(c) == syntax::N_USES && name(c) == "linear_elastic" { return true; }
        if kind(c) == syntax::N_ASSUME && assumes_finite(c) { return true; }
        c = next(c);
    }
    return false;
}

let SVK_FIX = "write assume finite_strain, isotropic, plane_strain, let S = hooke(E, nu, strain(u)) and eq div(defgrad(u)*S) = 0";
var solid_finite = false;
var solid_stress: u8[..64];

fn solid_finite_text(c: i64) -> bool {
    if compact(child(c, 1)) != str("hooke(E,nu,strain(", displacement, "))") { return false; }
    set64(solid_stress, name(c));
    return true;
}

fn solid_item(n: i64) {
    if solid_tok >= 0 {
        unsupported(tok(n), "a second solid: the flow generator has one", "remove the second solid");
        return;
    }
    solid_tok = tok(n);
    set64(solid_name, name(n));
    var plane = false;
    var hooke_seen = false;
    var balance_seen = false;
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            solid_domain = named_of(first(c), SHAPE);
            if solid_domain < 0 { unsupported(tok(first(c)), "the region of the solid is not a shape that a let names", SOLID_FIX); }
        } else if k == syntax::N_USES {
            if first(c) >= 0 { unsupported(tok(c), "the solid has no body force", "write linear_elastic"); }
        } else if k == syntax::N_UNKNOWN {
            set64(displacement, name(c));
        } else if k == syntax::N_ASSUME {
            var a = first(c);
            while a >= 0 {
                if name(a) == "plane_strain" { plane = true; }
                if name(a) == "finite_strain" { solid_finite = true; }
                if name(a) != "plane_strain" && name(a) != "small_strain" && name(a) != "isotropic" && name(a) != "finite_strain" {
                    unsupported(tok(a), str("the assumption '", name(a), "' is not in the flow generator"), "use plane_strain");
                }
                a = next(a);
            }
        } else if k == syntax::N_LET && solid_finite_text(c) {
            hooke_seen = true;
        } else if k == syntax::N_EQN && compact(child(c, 0)) == str("div(defgrad(", displacement, ")*", solid_stress, ")") && compact(child(c, 1)) == "0" {
            balance_seen = true;
        } else if k != syntax::N_MATASSIGN {
            unsupported(tok(c), "a physics with a template has no let or eq", SVK_FIX);
        }
        c = next(c);
    }
    if solid_finite && !(hooke_seen && balance_seen) {
        unsupported(tok(n), "the finite-strain solid is St. Venant-Kirchhoff", SVK_FIX);
    }
    if !plane { unsupported(tok(n), "the 2d solid of the flow generator is in plane strain", "add assume plane_strain"); }
}

fn interface_item(n: i64) {
    if interface_node >= 0 {
        unsupported(tok(n), "a second interface: the flow generator has one", "remove the second interface");
        return;
    }
    interface_node = n;
    var harmonic = false;
    var s = first(child(n, 2));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "mesh_motion" && kind(v) == syntax::N_NAME && name(v) == "harmonic" {
            harmonic = true;
        } else {
            unsupported(tok(s), "the interface of the flow generator has mesh_motion: harmonic only", "write interface fsi: wet_fluid, wet_solid { mesh_motion: harmonic }");
        }
        s = next(s);
    }
    if !harmonic { unsupported(tok(n), "the interface needs a mesh motion", "add mesh_motion: harmonic"); }
}

fn solid_method(n: i64) {
    if syntax::text(syntax::nodes[n].aux) != "cg" { unsupported(tok(n), "the solid has the cg method", "write method flag: cg { order: 2 }"); }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            if scalar_value(first(s)) != 2.0 { unsupported(tok(first(s)), "the solid has order 2", "write order: 2"); }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting of the solid is not in the flow generator", "use only mesh and order");
        }
        s = next(s);
    }
}

fn solid_mesh(n: i64) {
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            solid_size = scaled_value(v);
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "cell" && kind(v) == syntax::N_NAME && name(v) == "tri") {
            unsupported(tok(s), "the mesh of the solid has a size and tri cells", "write mesh ms on bar { size: 5 mm, cell: tri }");
        }
        s = next(s);
    }
    if solid_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 1 cm"); }
}

fn thermal_solid_mesh(n: i64) {
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            th_size = scaled_value(v);
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "cell" && kind(v) == syntax::N_NAME && name(v) == "tri") {
            unsupported(tok(s), "the mesh of the solid has a size and tri cells", "write mesh ms on wall { size: 1 mm, cell: tri }");
        }
        s = next(s);
    }
    if th_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 1 cm"); }
}

fn thermal_solid_method(n: i64) {
    if syntax::text(syntax::nodes[n].aux) != "cg" { unsupported(tok(n), "the solid with heat conduction has the cg method", "write method solid: cg { mesh: ms }"); }
    var s = first(child(n, 0));
    while s >= 0 {
        if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting of the solid is not in the flow generator", "use only mesh");
        }
        s = next(s);
    }
}

fn read_thermal(at: i64) {
    if th_tok >= 0 && !thermal { unsupported(th_tok, "a solid with heat conduction needs a flow with a temperature", "add Tf: Temperature to the flow"); }
    if !thermal { return; }
    if property("cp") < 0 || property("k") < 0 {
        unsupported(at, "a flow with a temperature needs a material with cp and k", "add cp: SpecificHeat = 1000 J/(kg*K) and k: Conductivity = 0.025 W/(m*K)");
        return;
    }
    if th_tok >= 0 && th_size <= 0.0 && ok() { unsupported(th_tok, "the solid has no mesh", "add mesh ms on wall { size: 1 mm, cell: tri }"); }
    if turbulent || solid_tok >= 0 || period_high >= 0 { unsupported(at, "a flow with a temperature is laminar, with no solid linear_elastic and not periodic", "remove them"); }
}

fn read_solid(at: i64) {
    if solid_tok < 0 {
        if interface_node >= 0 { unsupported(tok(interface_node), "an interface needs a solid physics linear_elastic", SOLID_FIX); }
        return;
    }
    if interface_node < 0 { unsupported(solid_tok, "the solid needs an interface with the flow", "write interface fsi: wet_fluid, wet_solid { mesh_motion: harmonic }"); }
    if solid_domain >= 0 && solid_size <= 0.0 { unsupported(solid_tok, "the solid has no mesh", "add mesh ms on bar { size: 5 mm, cell: tri }"); }
    if solid_material < 0 { unsupported(solid_tok, "the solid needs a material with E and nu", "add material steel { E: Pressure = 200 GPa, nu: Real = 0.3 }"); }
    if order != 2 { unsupported(at, "a flow with a solid has order 2", "write method flow: cg { order: { v: 2, p: 1 } }"); }
    if turbulent || period_high >= 0 { unsupported(solid_tok, "a flow with a solid is laminar and not periodic", "remove k_omega_sst and the symmetry"); }
}

fn assumptions(n: i64) {
    var a = first(n);
    while a >= 0 {
        let w = name(a);
        if w != "incompressible" && w != "newtonian" && w != "laminar" && !(turbulent && w == "turbulent") {
            unsupported(tok(a), str("the assumption '", w, "' is not in the flow generator"), "use incompressible, newtonian and laminar");
        }
        a = next(a);
    }
}

fn method(n: i64) {
    let kind_name = syntax::text(syntax::nodes[n].aux);
    if sa_d_tok >= 0 && name(n) == sa_d_name && kind_name == "cg" { return; }
    if spalart && name(n) == physics_name && kind_name == "fv" {
        order = 1;
        stabilize = 2;
        stabilize_set = true;
        return;
    }
    if name(n) != physics_name || kind_name != "cg" {
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

let REFINE_FIX = "refine a boundary on a side of a rect, boundary(box, y == 1 m), or on a circle, boundary(dom, hole)";

recursive fn refine_parts(shape: i64, cond: i64, h: f64) {
    if kind(cond) == syntax::N_BINARY && syntax::text(tok(cond)) == "||" {
        refine_parts(shape, child(cond, 0), h);
        refine_parts(shape, child(cond, 1), h);
        return;
    }
    let g = named_of(cond, SHAPE);
    if g >= 0 && is_call(named[g].value, "circle") {
        let c = named[g].value;
        let centre = scaled(quantity(arg(c, 0)));
        refines.push(Refine { circle: true, p: [centre.v[0], centre.v[1], scaled_value(arg(c, 1)), 0.0], h: h });
        return;
    }
    if g >= 0 && is_call(named[g].value, "rect") {
        edge_refines(named[g].value, h, refines);
        return;
    }
    if g >= 0 && is_fill(named[g].value) {
        fill_refines(named[g].value, h, refines);
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

fn fill_points(n: i64, xy: f64[>..]&) {
    let sp = arg(n, 0);
    var c = 1;
    while c < syntax::child_count(sp) {
        let q = scaled(quantity(arg(sp, c - 1)));
        if q.n != 2 { unsupported(tok(sp), "a point of a spline is a vector (x, y)", "write (0.5, 0.1) m"); return; }
        xy.push(q.v[0]);
        xy.push(q.v[1]);
        c++;
    }
}

fn fill_refines(n: i64, h: f64, into: Refine[>..]&) {
    var xy: f64[>..] = [];
    fill_points(n, xy);
    let m = xy.len / 2;
    for i in m {
        let j = (i + 1) % m;
        into.push(Refine { circle: false, p: [xy[2 * i], xy[2 * i + 1], xy[2 * j], xy[2 * j + 1]], h: h });
    }
}

recursive fn edge_refines(n: i64, h: f64, into: Refine[>..]&) {
    if kind(n) == syntax::N_NAME {
        edge_refines(named[named_of(n, SHAPE)].value, h, into);
    } else if is_call(n, "circle") {
        let c = scaled(quantity(arg(n, 0)));
        into.push(Refine { circle: true, p: [c.v[0], c.v[1], scaled_value(arg(n, 1)), 0.0], h: h });
    } else if is_call(n, "rect") {
        let a = scaled(quantity(arg(n, 0)));
        let b = scaled(quantity(arg(n, 1)));
        let x: f64[2] = [min(a.v[0], b.v[0]), max(a.v[0], b.v[0])];
        let y: f64[2] = [min(a.v[1], b.v[1]), max(a.v[1], b.v[1])];
        for k in 2 {
            into.push(Refine { circle: false, p: [x[k], y[0], x[k], y[1]], h: h });
            into.push(Refine { circle: false, p: [x[0], y[k], x[1], y[k]], h: h });
        }
    } else if is_fill(n) {
        fill_refines(n, h, into);
    } else {
        edge_refines(child(n, 0), h, into);
        edge_refines(child(n, 1), h, into);
    }
}

fn layers_item(s: i64) {
    let v = child(s, 0);
    let r = named_of(v, BOUNDARY);
    if r < 0 || !is_call(named[r].value, "boundary") || syntax::child_count(named[r].value) != 3 || named_of(arg(named[r].value, 1), SHAPE) < 0 {
        unsupported(tok(v), "the layers are not on a boundary(dom, shape) that a let names", "write layers wall: { first: 4 um, ratio: 1.2, count: 20 } with let wall: Curve = boundary(dom, hole)");
        return;
    }
    if layer_region >= 0 { unsupported(tok(s), "a second layers: the flow generator has one", "remove it"); return; }
    layer_region = r;
    var a = first(child(s, 1));
    while a >= 0 {
        let w = syntax::text(syntax::nodes[a].aux);
        if w == "first" { layer_first = scaled_value(first(a)); }
        else if w == "ratio" { layer_ratio = scalar_value(first(a)); }
        else if w == "count" { layer_count = scalar_value(first(a)) as! i64; }
        else { unsupported(tok(a), "layers has first, ratio and count", "write layers wall: { first: 4 um, ratio: 1.2, count: 20 }"); }
        a = next(a);
    }
    if layer_first <= 0.0 || layer_ratio < 1.0 || layer_count < 1 {
        unsupported(tok(s), "layers needs first > 0, ratio >= 1 and count >= 1", "write layers wall: { first: 4 um, ratio: 1.2, count: 20 }");
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
            } else if !is_call(named[r].value, "boundary") {
                unsupported(tok(v), "the refine is not on a boundary(shape, condition)", REFINE_FIX);
            } else {
                let b = named[r].value;
                refine_parts(arg(b, 0), arg(b, 1), scaled_value(child(s, 1)));
            }
        } else if kind(s) == syntax::N_LAYERS {
            layers_item(s);
        } else if kind(s) != syntax::N_SETTING {
            unsupported(tok(s), "crack_tip is not in the flow generator", "remove the line");
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

let BC_FIX_TEXT = "write bc wall: v = 0, bc lid: v = (1, 0) m/s, bc side: v.y = 0, bc corner: p = 0, or bc outlet: sigma*n = 0";

fn bc(n: i64) -> i64 {
    let region = child(n, 0);
    let target = child(n, 1);
    let v = child(n, 2);
    let r = if kind(region) == syntax::N_NAME { find_named(name(region)) } else { -1 };
    if r < 0 || named[r].kind == SHAPE || named[r].kind == LINE {
        unsupported(tok(region), "the region of the bc is not a boundary or a point that a let names", BC_FIX_TEXT);
        return -1;
    }
    var b = Bc { region: r, kind: BC_V, comp: -1, v: [0.0, 0.0], expr: -1 };
    let t = compact(target);
    if (spalart && t == sa_nt) || (sa_d_tok >= 0 && t == sa_d) {
        b.kind = if t == sa_d { BC_D } else { BC_NT };
        b.v[0] = scaled_value(v);
        if b.kind == BC_D && b.v[0] != 0.0 { unsupported(tok(v), "the wall distance is 0 on its wall", str("write bc ", named[r].name, ": ", sa_d, " = 0")); }
        use_region(r);
        bcs.push(b);
        return bcs.len - 1;
    }
    if thermal && is_temperature_target(t) { return temperature_bc(n, b, t); }
    if t == str(stress, "*n") && is_pressure_level(v) { return -1; }
    if solid_tok >= 0 && t == displacement { return solid_bc(n, b); }
    if mentions_xy(v) && (t == velocity || t == str(velocity, ".x") || t == str(velocity, ".y")) {
        if t == velocity && (kind(v) != syntax::N_VECTOR || syntax::child_count(v) != 2) { unsupported(tok(v), str(velocity, " is a vector"), "write v = (6*U*y*(H - y)/H^2, 0)"); }
        b.expr = v;
        b.comp = if t == velocity { -1 } else if t[t.len - 1] == 'x' { 0 } else { 1 };
        use_region(r);
        bcs.push(b);
        return bcs.len - 1;
    }
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

fn is_temperature_target(t: u8[:]) -> bool { return thermal_kind(t) != 0; }

fn thermal_kind(t: u8[:]) -> i64 {
    if t == temp || (th_tok >= 0 && t == th_temp) { return 1; }
    if t == str("dot(", flux, ",n)") || (th_tok >= 0 && t == str("dot(", th_flux, ",n)")) { return 2; }
    return 0;
}

fn temperature_bc(n: i64, b: Bc, t: u8[:]) -> i64 {
    let v = child(n, 2);
    let w = compact(v);
    if thermal_kind(w) != 0 && thermal_kind(w) == thermal_kind(t) && w != t {
        if th_tok < 0 { unsupported(tok(v), "an interface needs a solid with heat conduction", THERMAL_SOLID_FIX); }
        return -1;
    }
    var c = b;
    c.kind = if thermal_kind(t) == 1 { BC_T } else { BC_FLUX };
    c.v[0] = scaled(quantity(v)).v[0];
    use_region(b.region);
    bcs.push(c);
    return bcs.len - 1;
}

fn is_pressure_level(v: i64) -> bool {
    guard kind(v) == syntax::N_BINARY && syntax::text(tok(v)) == "*" && compact(child(v, 1)) == "n" else { return false; }
    let level = -scaled(quantity(child(v, 0))).v[0];
    if p_offset_set && abs(level - p_offset) > 1.0e-12 * max(abs(level), 1.0) {
        unsupported(tok(v), "the tractions -p*n of the bcs have different levels p", "give each traction bc the same pressure, or sigma*n = 0");
    }
    p_offset = level;
    p_offset_set = true;
    return true;
}

fn on_solid(g: i64) -> bool {
    let v = named[g].value;
    return solid_domain >= 0 && named[g].kind == BOUNDARY && is_call(v, "boundary") && named_of(arg(v, 0), SHAPE) == solid_domain;
}

fn on_thermal_solid(g: i64) -> bool {
    let v = named[g].value;
    return th_domain >= 0 && named[g].kind == BOUNDARY && is_call(v, "boundary") && named_of(arg(v, 0), SHAPE) == th_domain;
}

fn solid_bc(n: i64, b: Bc) -> i64 {
    let v = child(n, 2);
    let q = quantity(v);
    if q.v[0] != 0.0 || q.v[1] != 0.0 { unsupported(tok(v), "the flow generator holds the solid at 0", str("write bc root: ", displacement, " = 0")); }
    if !on_solid(b.region) { unsupported(tok(child(n, 0)), "the bc of the solid is not on a boundary of the solid", "write bc root: u = 0, with let root: Curve = boundary(bar, ...)"); }
    var c = b;
    c.kind = BC_U;
    use_region(b.region);
    bcs.push(c);
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
                   first_report: reports.len, report_count: 0, rtol: 1.0e-8, output: false, solid_output: false };
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
            s.solid_output = output_fields(c);
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

fn output_fields(n: i64) -> bool {
    var solid = false;
    var o = first(n);
    while o >= 0 {
        if solid_tok >= 0 && kind(o) == syntax::N_NAME && name(o) == displacement {
            solid = true;
        } else if kind(o) == syntax::N_EVERY {
            unsupported(tok(o), "a static step has no every", "remove it");
        } else if kind(o) == syntax::N_NAME && thermal && (name(o) == temp || (th_tok >= 0 && name(o) == th_temp)) {
        } else if kind(o) != syntax::N_NAME || (name(o) != velocity && name(o) != pressure && !(turbulent && (name(o) == turb_k || name(o) == turb_w)) && !(spalart && (name(o) == sa_nt || name(o) == sa_d))) {
            unsupported(tok(o), "the flow generator writes the velocity and the pressure, and k and w of k_omega_sst", str("write output ", velocity, ", ", pressure));
        }
        o = next(o);
    }
    return solid;
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

fn point_field(e: i64) -> i64, i64 {
    let f = compact(child(e, 0));
    if f == velocity { return 0, 1; }
    if f == str(velocity, ".x") { return 0, -1; }
    if f == str(velocity, ".y") { return 1, -1; }
    if f == pressure { return 2, -1; }
    if thermal && (f == temp || (th_tok >= 0 && f == th_temp)) { return 5, -1; }
    if solid_tok >= 0 && f == displacement { return 3, 4; }
    if solid_tok >= 0 && f == str(displacement, ".x") { return 3, -1; }
    if solid_tok >= 0 && f == str(displacement, ".y") { return 4, -1; }
    return -1, -1;
}

fn whole_vector(e: i64) -> bool {
    guard kind(e) == syntax::N_AT else { return false; }
    let what, other = point_field(e);
    return other >= 0;
}

fn at_point(p: i64) -> Q {
    let i = named_of(p, POINT);
    return quantity(if i >= 0 { arg(named[i].value, 0) } else { p });
}

fn point_ok(e: i64) -> bool {
    let before = syntax::errors.len;
    let p = at_point(child(e, 1));
    if p.n != 2 && syntax::errors.len == before {
        unsupported(tok(child(e, 1)), "the point of 'at' is not a vector (x, y)", "write at (x, y) m, or the name of a point((x, y) m)");
    }
    return syntax::errors.len == before;
}

var tf_a = 0.0;
var tf_b = 0.0;
var tf_m = 0.0;
var tf_x0 = 0.0;
var tf_y0 = 0.0;

fn traction_integral(e: i64) -> i64 {
    guard is_call(e, "integral") && syntax::child_count(e) == 3 && named_of(arg(e, 1), BOUNDARY) >= 0 else { return -1; }
    let f = compact(arg(e, 0));
    tf_a = 0.0;
    tf_b = 0.0;
    tf_m = 0.0;
    if f == str("(", stress, "*n).x") { tf_a = 1.0; return 0; }
    if f == str("(", stress, "*n).y") { tf_b = 1.0; return 0; }
    let w = arg(e, 0);
    if is_call(w, "dot") && syntax::child_count(w) == 3 && compact(arg(w, 0)) == str(stress, "*n") {
        let q = quantity(arg(w, 1));
        if q.n != 2 || !same_dim(q.d, [0, 0, 0]) {
            unsupported(tok(w), "the vector of dot(sigma*n, e) is two numbers", "write dot(sigma*n, e_l) with let e_l: Vector<Real> = (-sin(aoa), cos(aoa))");
            return -1;
        }
        tf_a = q.v[0];
        tf_b = q.v[1];
        return 0;
    }
    if kind(w) == syntax::N_BINARY && syntax::text(tok(w)) == "-" {
        let l = child(w, 0);
        let r = child(w, 1);
        let lx = if kind(l) == syntax::N_BINARY && syntax::text(tok(l)) == "*" { child(l, 0) } else { -1 };
        let ry = if kind(r) == syntax::N_BINARY && syntax::text(tok(r)) == "*" { child(r, 0) } else { -1 };
        if lx >= 0 && ry >= 0 && compact(l) != "" && compact(child(l, 1)) == str("(", stress, "*n).y") && compact(child(r, 1)) == str("(", stress, "*n).x") &&
           kind(lx) == syntax::N_BINARY && syntax::text(tok(lx)) == "-" && compact(child(lx, 0)) == "x" &&
           kind(ry) == syntax::N_BINARY && syntax::text(tok(ry)) == "-" && compact(child(ry, 0)) == "y" {
            tf_m = 1.0;
            tf_x0 = scaled_value(child(lx, 1));
            tf_y0 = scaled_value(child(ry, 1));
            return 1;
        }
    }
    return -1;
}

let REGION_FIX = "write mean(-grad(T).x, wall), integral(T*dot(v, n), outlet) or max(v.x, line_name), with let line_name: Curve = line((0.5, 0) m, (0.5, 1) m)";

fn field_text(e: i64, out: u8[>..]&) -> i64[3], bool {
    let f = compact(e);
    let zero: i64[3] = [0, 0, 0];
    if f == str(velocity, ".x") { out.append("s.vx"); return VELOCITY_D, true; }
    if f == str(velocity, ".y") { out.append("s.vy"); return VELOCITY_D, true; }
    if f == pressure {
        out.append("(s.p + ");
        lit(out, p_offset);
        out.push(')');
        return PRESSURE_D, true;
    }
    if thermal && (f == temp || (th_tok >= 0 && f == th_temp)) { out.append("s.t"); return zero, true; }
    if thermal && (f == str("grad(", temp, ").x") || (th_tok >= 0 && f == str("grad(", th_temp, ").x"))) { out.append("s.tx"); return [-2, 0, 0], true; }
    if thermal && (f == str("grad(", temp, ").y") || (th_tok >= 0 && f == str("grad(", th_temp, ").y"))) { out.append("s.ty"); return [-2, 0, 0], true; }
    if f == str("dot(", velocity, ",n)") { out.append("(s.vx * s.nx + s.vy * s.ny)"); return VELOCITY_D, true; }
    if thermal && f == str("dot(", flux, ",n)") {
        out.append("(-");
        lit(out, scaled_value(property("k")));
        out.append(" * (s.tx * s.nx + s.ty * s.ny))");
        return [0, 2, -6], true;
    }
    if th_tok >= 0 && f == str("dot(", th_flux, ",n)") {
        out.append("(-");
        lit(out, scaled_value(property_of(th_mat, "k")));
        out.append(" * (s.tx * s.nx + s.ty * s.ny))");
        return [0, 2, -6], true;
    }
    return zero, false;
}

fn try_field(out: u8[>..]&, e: i64) -> i64[3], bool {
    var text: u8[>..] = [];
    let d, field = field_text(e, text);
    if field { out.append(text); }
    return d, field;
}

recursive fn emit_integrand(out: u8[>..]&, e: i64) -> i64[3] {
    let d, field = try_field(out, e);
    if field { return d; }
    let k = kind(e);
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(e))) >= 0 {
        let op = syntax::text(tok(e));
        out.push('(');
        let a = emit_integrand(out, child(e, 0));
        format(out, " ", op, " ");
        let b = emit_integrand(out, child(e, 1));
        out.push(')');
        if op == "*" { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
        if op == "/" { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
        return a;
    }
    if k == syntax::N_UNARY && syntax::text(tok(e)) == "-" {
        out.append("(-");
        let a = emit_integrand(out, first(e));
        out.push(')');
        return a;
    }
    let q = scaled(quantity(e));
    out.push('(');
    lit(out, q.v[0]);
    out.push(')');
    return quantity(e).d;
}

fn region_ok(e: i64) -> bool {
    guard syntax::child_count(e) == 3 else {
        unsupported(tok(e), "mean, integral and max have a field and a region", REGION_FIX);
        return false;
    };
    let r = if kind(arg(e, 1)) == syntax::N_NAME { find_named(name(arg(e, 1))) } else { -1 };
    let wanted = if is_call(e, "max") { LINE } else { BOUNDARY };
    if r < 0 || named[r].kind != wanted {
        unsupported(tok(arg(e, 1)), str("the region of ", if is_call(e, "max") { "max" } else { "mean and integral" }, " is not a ", if is_call(e, "max") { "line" } else { "boundary" }, " that a let names"), REGION_FIX);
        return false;
    }
    if !thermal && !is_call(e, "max") {
        unsupported(tok(e), "mean and integral of a field are in a flow with a temperature", "use integral((sigma*n).x, wall)");
        return false;
    }
    var scratch: u8[>..] = [];
    emit_integrand(scratch, arg(e, 0));
    if wanted == BOUNDARY { use_region(r); }
    return ok();
}

recursive fn result_ok(e: i64) -> bool {
    let k = kind(e);
    if k == syntax::N_AT {
        let what, other = point_field(e);
        if what < 0 { unsupported(tok(e), str("'", source_of(child(e, 0)), "' at a point is not in the flow generator"), REPORT_FIX); }
        if other >= 0 { unsupported(tok(e), str("the vector '", source_of(child(e, 0)), "' at a point is only in a report, not in an expect or a + - * /"), REPORT_FIX); }
        return what >= 0 && other < 0 && point_ok(e);
    }
    if is_call(e, "integral") && traction_integral(e) < 0 && thermal { return region_ok(e); }
    if is_call(e, "mean") || (is_call(e, "max") && syntax::child_count(e) == 3) { return region_ok(e); }
    if is_call(e, "integral") {
        let form = traction_integral(e);
        if form >= 0 { use_region(named_of(arg(e, 1), BOUNDARY)); }
        if form < 0 && ok() { unsupported(tok(e), "the flow generator integrates (sigma*n).x, dot(sigma*n, e) or the moment about a point on a boundary", REPORT_FIX); }
        return form >= 0;
    }
    if k == syntax::N_BINARY && position(["+", "-", "*", "/"], syntax::text(tok(e))) >= 0 {
        return result_ok(child(e, 0)) && result_ok(child(e, 1));
    }
    if k == syntax::N_UNARY && syntax::text(tok(e)) == "-" { return result_ok(first(e)); }
    let q = quantity(e);
    return ok() && q.n >= 1;
}

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

fn row_count(t: i64) -> i64 {
    var k = 0;
    var c = first(t);
    while c >= 0 {
        if kind(c) == syntax::N_ROW { k++; }
        c = next(c);
    }
    return k;
}

fn read(prog: i64) -> bool {
    reset();
    var physics_nodes: i64[>..] = [];
    var mesh_nodes: i64[>..] = [];
    var methods: i64[>..] = [];
    var symmetry_node = -1;
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
        } else if k == syntax::N_FN {
            var w: u8[..64] = [];
            set64(w, name(c));
            fn_names.push(w);
            fn_nodes.push(c);
        } else if k == syntax::N_MATERIAL {
            material(c);
        } else if k == syntax::N_PHYSICS && is_thermal_solid(c) {
            thermal_solid_item(c);
        } else if k == syntax::N_PHYSICS && is_solid(c) {
            solid_item(c);
        } else if k == syntax::N_PHYSICS {
            physics_nodes.push(c);
        } else if k == syntax::N_MESH {
            mesh_nodes.push(c);
        } else if k == syntax::N_INTERFACE {
            interface_item(c);
        } else if k == syntax::N_METHOD {
            methods.push(c);
        } else if k == syntax::N_TABLE {
            tables.push(c);
        } else if k == syntax::N_SCALES {
            has_scales = true;
        } else if k == syntax::N_SYMMETRY && symmetry_node < 0 {
            symmetry_node = c;
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
    for n in physics_nodes {
        if is_walldist(n) { walldist_item(n); } else { physics_item(n); }
    }
    if physics_nodes.len == 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if !ok() { return false; }
    read_fluid(at);
    for m in mesh_nodes {
        let r = child(m, 0);
        if solid_domain >= 0 && kind(r) == syntax::N_NAME && name(r) == named[solid_domain].name { solid_mesh(m); }
        else if th_domain >= 0 && kind(r) == syntax::N_NAME && name(r) == named[th_domain].name { thermal_solid_mesh(m); }
        else { mesh(m); }
    }
    if mesh_size <= 0.0 && ok() { unsupported(at, "the model has no mesh", "add mesh m on dom { size: 1 cm }"); }
    for m in methods {
        if solid_tok >= 0 && name(m) == solid_name { solid_method(m); }
        else if th_tok >= 0 && name(m) == th_name { thermal_solid_method(m); }
        else { method(m); }
    }
    read_thermal(at);
    if !stabilize_set { stabilize = if convective { 1 } else { 0 }; }
    if order == 1 && stabilize < 2 { unsupported(at, "order 1 for the velocity needs pressure stabilization", "write stabilize: supg_pspg"); }
    if turbulent && (order != 1 || stabilize != 2) { unsupported(at, "k_omega_sst runs on P1 cells with supg_pspg", "write method flow: cg { order: 1, stabilize: supg_pspg }"); }
    if symmetry_node >= 0 { symmetry(symmetry_node); }
    read_solid(at);
    if !ok() { return false; }
    read_steps(prog, at);
    return ok();
}

let SYMMETRY_FIX = "write symmetry ends: bloch { period: (L, 0 m), low: inlet, high: outlet, k: (0, 0) m^-1 }, on a rect of width L";

fn symmetry(n: i64) {
    if name(child(n, 0)) != "bloch" {
        unsupported(tok(n), "the flow generator has the symmetry bloch with k = 0", SYMMETRY_FIX);
        return;
    }
    var s = first(syntax::nodes[child(n, 0)].next);
    while s >= 0 {
        let w = name(s);
        let v = first(s);
        if w == "period" {
            period = scaled(quantity(v)).v;
        } else if w == "k" {
            let k = quantity(v);
            if k.v[0] != 0.0 || k.v[1] != 0.0 { unsupported(tok(v), "the flow generator has bloch with k = 0 (periodic)", SYMMETRY_FIX); }
        } else if w == "high" {
            period_high = named_of(v, BOUNDARY);
            if period_high < 0 { unsupported(tok(v), "high is a boundary that a let names", SYMMETRY_FIX); }
        }
        s = next(s);
    }
    let rect = named[domain].value;
    let width = if is_call(rect, "rect") { abs(scaled(quantity(arg(rect, 1))).v[0] - scaled(quantity(arg(rect, 0))).v[0]) } else { -1.0 };
    if period[1] != 0.0 || abs(period[0] - width) > 1.0e-9 * width {
        unsupported(tok(n), "the period is the width of the rect of the physics, along x", SYMMETRY_FIX);
    }
    if period_high >= 0 { use_region(period_high); }
}

fn fluid_mu() -> f64 {
    if property("mu") >= 0 { return scaled_value(property("mu")); }
    return scaled_value(property("nu")) * scaled_value(property("rho"));
}

let SA_CONSTANTS = ["sig_sa", "kap", "cb1", "cb2", "cw1", "cw2", "cw3", "cv1", "c_s"];

fn constant(word: u8[:]) -> f64 {
    let v = codegen::find_value(word);
    if v < 0 {
        unsupported(tok(material_node), str("the Spalart-Allmaras model needs the constant '", word, "'"), str("write let ", word, ": Real = ..."));
        return 1.0;
    }
    return scalar_value(v);
}

fn read_fluid(at: i64) {
    if spalart { for i in 9 { constant(SA_CONSTANTS[i]); } }
    let r = property("rho");
    let m = if property("mu") >= 0 { property("mu") } else { property("nu") };
    if material_node < 0 || r < 0 || m < 0 {
        unsupported(at, "the flow needs a material with rho and mu or nu", "add material fluid { rho: Density = 1000 kg/m^3, mu: DynViscosity = 1e-3 Pa*s }");
        return;
    }
    if fluid_mu() <= 0.0 || scalar_value(r) <= 0.0 { unsupported(tok(material_node), "rho and mu are more than 0", "give rho > 0 and mu > 0"); }
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
    var scratch: u8[>..] = [];
    emit_selects(scratch);
}

fn merge(set: Bc[>..]&, b: Bc) {
    for i in 0..set.len {
        if set[i].region == b.region && set[i].kind == b.kind && set[i].comp == b.comp {
            set[i] = b;
            return;
        }
    }
    set.push(b);
}

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
    if is_fill(n) {
        var xy: f64[>..] = [];
        fill_points(n, xy);
        return mesher::fill(xy[..]);
    }
    let a = mesher_shape(child(n, 0));
    let b = mesher_shape(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { mesher::subtract(a, b) } else if op == "|" { mesher::union(a, b) } else { mesher::intersect(a, b) };
}

fn write_mesh(dir: u8[:]) -> bool {
    mesher::clear_shapes();
    mesher::clear_refines();
    mesher::set_grade(grade);
    for r in refines {
        if r.circle { mesher::refine_circle(r.p[0], r.p[1], r.p[2], r.h); } else { mesher::refine_segment(r.p[0], r.p[1], r.p[2], r.p[3], r.h); }
    }
    var shape = mesher_shape(named[domain].value);
    mesher::clear_interfaces();
    if solid_tok >= 0 || th_tok >= 0 {
        let sd = if solid_tok >= 0 { solid_domain } else { th_domain };
        let ssize = if solid_tok >= 0 { solid_size } else { th_size };
        let solid = mesher_shape(named[sd].value);
        var list: Refine[>..] = [];
        edge_refines(named[sd].value, ssize, list);
        for r in list { if !r.circle { mesher::refine_segment(r.p[0], r.p[1], r.p[2], r.p[3], r.h); } }
        shape = mesher::union(shape, solid);
        mesher::add_interface(solid);
    }
    var middle = 0.0;
    if period_high >= 0 {
        let r = named[domain].value;
        let a = scaled(quantity(arg(r, 0)));
        let b = scaled(quantity(arg(r, 1)));
        middle = 0.5 * (a.v[0] + b.v[0]);
        shape = mesher::rect(min(a.v[0], b.v[0]), a.v[1], middle, b.v[1]);
    }
    var meshed = false;
    if layer_region >= 0 {
        let wall = mesher_shape(named[named_of(arg(named[layer_region].value, 1), SHAPE)].value);
        meshed = mesher::mesh_2d_layers(shape, wall, mesh_size, layer_first, layer_ratio, layer_count);
    } else {
        meshed = mesher::mesh_2d(shape, mesh_size);
    }
    mesher::clear_interfaces();
    if !meshed {
        write_stderr(str("skein run: the mesher cannot mesh '", named[domain].name, "' with the size and the refines of the mesh\n"));
        return false;
    }
    guard solid_tok < 0 else { return write_coupled(dir); }
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
    if period_high >= 0 { mirror(xy, tri, middle); }
    if !write_masks(dir, xy) { return false; }
    cells = tri.len / 3;
    vertices = xy.len / 2;
    if !write_file(str(dir, "/flow.xy"), to_bytes(xy[..])) || !write_file(str(dir, "/flow.tri"), to_bytes(tri[..])) {
        write_stderr(str("skein run: cannot write the mesh in '", dir, "'\n"));
        return false;
    }
    return th_tok < 0 || write_parts(dir);
}

fn needs_mask(g: i64) -> bool {
    let b = named[g].value;
    if named[g].kind != BOUNDARY || !is_call(b, "boundary") || is_union(b) { return false; }
    return mask_shape(arg(b, 0)) >= 0 || (syntax::child_count(b) >= 3 && mask_shape(arg(b, 1)) >= 0);
}

fn mask_shape(n: i64) -> i64 {
    let s = named_of(n, SHAPE);
    if s < 0 || s == domain || s == solid_domain || s == th_domain { return -1; }
    let v = named[s].value;
    if is_call(v, "rect") || is_call(v, "circle") { return -1; }
    return s;
}

fn write_masks(dir: u8[:], xy: f64[:]) -> bool {
    for r in used_regions {
        if !needs_mask(r) { continue; }
        let b = named[r].value;
        let a = mask_shape(arg(b, 0));
        let c = if syntax::child_count(b) >= 3 { mask_shape(arg(b, 1)) } else { -1 };
        var mask: u8[>..] = [];
        let ha = if a >= 0 { mesher_shape(named[a].value) } else { -1 };
        let hc = if c >= 0 { mesher_shape(named[c].value) } else { -1 };
        let eps = 1.0e-6 * mesh_diameter();
        for n in xy.len / 2 {
            var on = true;
            if ha >= 0 { on = on && abs(mesher::phi(ha, xy[2 * n], xy[2 * n + 1], 0.0)) <= eps; }
            if hc >= 0 { on = on && abs(mesher::phi(hc, xy[2 * n], xy[2 * n + 1], 0.0)) <= eps; }
            mask.push(if on { 1 } else { 0 });
        }
        if !write_file(str(dir, "/", named[r].name, ".mask"), mask) {
            write_stderr(str("skein run: cannot write the mask of '", named[r].name, "' in '", dir, "'\n"));
            return false;
        }
    }
    return true;
}

fn mesh_diameter() -> f64 {
    let b = mesher::bounds(mesher_shape(named[domain].value));
    let w = b[3] - b[0];
    let h = b[4] - b[1];
    return sqrt(w * w + h * h);
}

fn write_parts(dir: u8[:]) -> bool {
    let solid = mesher_shape(named[th_domain].value);
    var part: u8[>..] = [];
    for c in mesher::cell_count() {
        var x = 0.0;
        var y = 0.0;
        for k in 3 {
            x += mesher::vertex(mesher::cell_node(c, k))[0] / 3.0;
            y += mesher::vertex(mesher::cell_node(c, k))[1] / 3.0;
        }
        part.push(if mesher::phi(solid, x, y, 0.0) < 0.0 { 1 } else { 0 });
    }
    if !write_file(str(dir, "/flow.part"), to_bytes(part[..])) {
        write_stderr(str("skein run: cannot write the mesh in '", dir, "'\n"));
        return false;
    }
    return true;
}

fn write_coupled(dir: u8[:]) -> bool {
    let solid = mesher_shape(named[solid_domain].value);
    var part: i64[>..] = [];
    for c in mesher::cell_count() {
        var x = 0.0;
        var y = 0.0;
        for k in 3 {
            x += mesher::vertex(mesher::cell_node(c, k))[0] / 3.0;
            y += mesher::vertex(mesher::cell_node(c, k))[1] / 3.0;
        }
        part.push(if mesher::phi(solid, x, y, 0.0) < 0.0 { 1 } else { 0 });
    }
    cells = 0;
    vertices = 0;
    for side in 2 {
        var xy: f64[>..] = [];
        var tri: i64[>..] = [];
        var index: i64[>..] = [];
        for n in mesher::point_count() { index.push(-1); }
        for c in mesher::cell_count() {
            if part[c] != side { continue; }
            for k in 3 {
                let n = mesher::cell_node(c, k);
                if index[n] < 0 {
                    index[n] = xy.len / 2;
                    xy.push(mesher::vertex(n)[0]);
                    xy.push(mesher::vertex(n)[1]);
                }
                tri.push(index[n]);
            }
        }
        let file = if side == 0 { "flow" } else { "solid" };
        if side == 0 {
            cells = tri.len / 3;
            vertices = xy.len / 2;
        }
        if tri.len == 0 || !write_file(str(dir, "/", file, ".xy"), to_bytes(xy[..])) || !write_file(str(dir, "/", file, ".tri"), to_bytes(tri[..])) {
            write_stderr(str("skein run: cannot write the mesh in '", dir, "'\n"));
            return false;
        }
    }
    return true;
}

fn mirror(xy: f64[>..]&, tri: i64[>..]&, middle: f64) {
    let n = xy.len / 2;
    let t = tri.len;
    let tol = 1.0e-9 * max(abs(middle), 1.0);
    var image: i64[>..] = [];
    for v in n {
        if abs(xy[2 * v] - middle) <= tol {
            image.push(v);
            continue;
        }
        image.push(xy.len / 2);
        xy.push(2.0 * middle - xy[2 * v]);
        xy.push(xy[2 * v + 1]);
    }
    for i in t { tri.push(image[tri[i]]); }
}

var cells = 0;
var vertices = 0;

fn dof_estimate() -> i64 {
    let area = shape_area(named[domain].value);
    var cell_count = area / (0.4330127 * mesh_size * mesh_size);
    for r in refines {
        let length = if r.circle { TAU * r.p[2] } else { sqrt((r.p[2] - r.p[0]) * (r.p[2] - r.p[0]) + (r.p[3] - r.p[1]) * (r.p[3] - r.p[1])) };
        let band = (mesh_size - r.h) / (grade - 1.0);
        cell_count += length * band / (0.4330127 * r.h * mesh_size);
    }
    if layer_region >= 0 {
        var h = mesh_size;
        for r in refines { h = min(h, r.h); }
        let perimeter = shape_perimeter(named[named_of(arg(named[layer_region].value, 1), SHAPE)].value);
        cell_count += 2.0 * layer_count * perimeter / h;
    }
    if th_tok >= 0 { cell_count += shape_area(named[th_domain].value) / (0.4330127 * th_size * th_size); }
    let vertex_count = 0.5 * cell_count;
    if thermal { return (4.0 * vertex_count + (if order == 2 { 3.0 * 1.5 * cell_count } else { 0.0 })) as! i64; }
    return (if turbulent { 5.0 * vertex_count } else if order == 1 { 3.0 * vertex_count } else { 3.0 * vertex_count + 2.0 * 1.5 * cell_count }) as! i64;
}

fn shape_perimeter(n: i64) -> f64 {
    var m = n;
    if kind(m) == syntax::N_NAME { m = named[named_of(m, SHAPE)].value; }
    if !is_fill(m) { return 0.0; }
    var xy: f64[>..] = [];
    fill_points(m, xy);
    var total = 0.0;
    let count = xy.len / 2;
    for i in count {
        let j = (i + 1) % count;
        total += sqrt((xy[2 * j] - xy[2 * i]) * (xy[2 * j] - xy[2 * i]) + (xy[2 * j + 1] - xy[2 * i + 1]) * (xy[2 * j + 1] - xy[2 * i + 1]));
    }
    return total;
}

recursive fn shape_area(n: i64) -> f64 {
    if kind(n) == syntax::N_NAME { return shape_area(named[named_of(n, SHAPE)].value); }
    if is_call(n, "rect") {
        let a = scaled(quantity(arg(n, 0)));
        let b = scaled(quantity(arg(n, 1)));
        return abs((b.v[0] - a.v[0]) * (b.v[1] - a.v[1]));
    }
    if is_call(n, "circle") { return PI * scaled_value(arg(n, 1)) * scaled_value(arg(n, 1)); }
    if is_fill(n) {
        var xy: f64[>..] = [];
        fill_points(n, xy);
        var twice = 0.0;
        let m = xy.len / 2;
        for i in m {
            let j = (i + 1) % m;
            twice += xy[2 * i] * xy[2 * j + 1] - xy[2 * j] * xy[2 * i + 1];
        }
        return abs(twice) / 2.0;
    }
    let a = shape_area(child(n, 0));
    let b = shape_area(child(n, 1));
    let op = syntax::text(tok(n));
    return if op == "-" { max(a - b, 0.0) } else if op == "|" { a + b } else { min(a, b) };
}

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

var vector_comp = -1;

let VALUE_NAMES = ["flow::VALUE_VX", "flow::VALUE_VY", "flow::VALUE_P", "flow::VALUE_UX", "flow::VALUE_UY", "flow::VALUE_T"];
let VELOCITY_D: i64[3] = [2, 0, -2];
let PRESSURE_D: i64[3] = [-2, 2, -4];
let LENGTH_D: i64[3] = [2, 0, 0];

recursive fn emit_result(out: u8[>..]&, e: i64) -> i64[3] {
    let k = kind(e);
    if k == syntax::N_AT {
        var what, other = point_field(e);
        if other >= 0 && vector_comp >= 0 { what += vector_comp; }
        let at = scaled(at_point(child(e, 1)));
        let d = if what == 2 { PRESSURE_D } else if what == 5 { [0, 0, 0] } else if what >= 3 { LENGTH_D } else { VELOCITY_D };
        format(out, "(at(", VALUE_NAMES[what], ", ");
        lit(out, at.v[0]);
        out.append(", ");
        lit(out, at.v[1]);
        out.append(", name)");
        if what == 2 && p_offset != 0.0 {
            out.append(" + ");
            lit(out, p_offset);
        }
        out.append(") * ");
        lit(out, reference(d));
        return d;
    }
    if is_call(e, "mean") || (is_call(e, "max") && syntax::child_count(e) == 3) || (is_call(e, "integral") && traction_integral(e) < 0) {
        return emit_region(out, e);
    }
    if is_call(e, "integral") {
        let r = named_of(arg(e, 1), BOUNDARY);
        let form = traction_integral(e);
        let d: i64[3] = if form == 1 { [2, 2, -4] } else { [0, 2, -4] };
        if form == 1 {
            out.append("wall_moment(select_");
            format(out, named[r].name, "(), ");
            lit(out, tf_x0);
            out.append(", ");
            lit(out, tf_y0);
            out.append(", ");
        } else {
            out.append("(");
            lit(out, tf_a);
            format(out, " * wall_force(select_", named[r].name, "(), 0, ");
            quoted(out, named[r].name);
            out.append(") + ");
            lit(out, tf_b);
            format(out, " * wall_force(select_", named[r].name, "(), 1, ");
            quoted(out, named[r].name);
            out.append("))");
        }
        if form == 1 {
            quoted(out, named[r].name);
            out.append(")");
        }
        out.append(" * ");
        lit(out, reference(d));
        return d;
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

fn emit_region(out: u8[>..]&, e: i64) -> i64[3] {
    let r = find_named(name(arg(e, 1)));
    var f: u8[>..] = [];
    var d = emit_integrand(f, arg(e, 0));
    if is_call(e, "max") {
        let a = scaled(quantity(arg(named[r].value, 0)));
        let b = scaled(quantity(arg(named[r].value, 1)));
        out.append("(flow::line_max(");
        lit(out, a.v[0]);
        out.append(", ");
        lit(out, a.v[1]);
        out.append(", ");
        lit(out, b.v[0]);
        out.append(", ");
        lit(out, b.v[1]);
        format(out, ") { s => ", f, " } * ");
    } else {
        format(out, "(flow::boundary_value(select_", named[r].name, "(), ", if is_call(e, "mean") { "true" } else { "false" }, ") { s => ", f, " } * ");
        if is_call(e, "integral") { d = [d[0] + 2, d[1], d[2]]; }
    }
    lit(out, reference(d));
    out.push(')');
    return d;
}

recursive fn emit_boundary(out: u8[>..]&, b: i64) {
    if is_union(b) {
        for k in 2 {
            out.append(if k == 0 { "(" } else { ") || (" });
            let c = child(b, k);
            let r = named_of(c, BOUNDARY);
            emit_boundary(out, if r >= 0 { named[r].value } else { c });
        }
        out.push(')');
        return;
    }
    let shape = arg(b, 0);
    let s = named_of(shape, SHAPE);
    if s != domain && (s < 0 || s != solid_domain) && (s < 0 || s != th_domain) {
        if s < 0 || !emit_shape_edge(out, named[s].value) {
            unsupported(tok(shape), "the shape of the boundary is not the domain, a rect or a circle", "write boundary(dom, ...)");
        }
        out.append(" && ");
    }
    if syntax::child_count(b) < 3 {
        out.append("true");
        return;
    }
    out.push('(');
    emit_cond(out, arg(b, 1));
    out.push(')');
}

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

    // The scaled moment of the force of the wall on the fluid about (x0, y0).
    fn wall_moment(found: bool, x0: f64, y0: f64, name: u8[:]) -> f64 {
        region(found, name);
        return flow::region_moment(x0, y0);
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
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

let COUPLED_PRELUDE = """
    fn solve_coupled(step: u8[:], rtol: f64) {
        let converged = flow::solve_fsi(rtol);
        format(solves, "  ", step, ",", flow::cell_count() + flow::solid_cells(), ",", flow::dof_count() + 2 * flow::solid_nodes(), ",",
               flow::newton_steps, ",", flow::linear_iterations, ",", flow::last_seconds, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        let r = flow::last_residual / max(flow::first_residual, 1.0e-300);
        check("solve", "residual", converged, str(r), str("<= ", rtol));
    }

    fn output_solid(path: u8[:], x_ref: f64) {
        if !flow::write_solid_vtu(path, x_ref) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    """;

var used_regions: i64[>..];

fn use_region(r: i64) { if position(used_regions, r) < 0 { used_regions.push(r); } }

fn emit_selects(out: u8[>..]&) {
    for r in used_regions {
        let g = named[r];
        if g.kind == BOUNDARY && needs_mask(r) {
            let b = g.value;
            format(out, "var mask_", g.name, ": u8[>..] = [];\n\nfn select_", g.name, "() -> bool {\n    return flow::select_masked(mask_", g.name, ") { x, y => ");
            if syntax::child_count(b) >= 3 && mask_shape(arg(b, 1)) < 0 { emit_cond(out, arg(b, 1)); } else { out.append("true"); }
            out.append(" };\n}\n\n");
            continue;
        }
        format(out, "fn select_", g.name, "() -> bool {\n    return ");
        if g.kind == POINT {
            let p = scaled(quantity(arg(g.value, 0)));
            out.append("flow::select_point(");
            lit(out, p.v[0]);
            out.append(", ");
            lit(out, p.v[1]);
            out.append(");\n}\n\n");
        } else {
            out.append(if on_solid(r) { "flow::select_solid_boundary() { x, y => " } else if on_thermal_solid(r) { "flow::select_solid_part() { x, y => " } else { "flow::select_boundary() { x, y => " });
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
    if b.kind == BC_U {
        out.append("    flow::fix_solid();\n");
        return;
    }
    if b.expr >= 0 {
        for j in 2 {
            if b.comp >= 0 && b.comp != j { continue; }
            format(out, "    flow::fix_velocity_at(", j, ") { x, y => ");
            emit_coordinate(out, if b.comp >= 0 { b.expr } else { child(b.expr, j) });
            out.append(" };\n");
        }
        return;
    }
    if b.kind == BC_T || b.kind == BC_FLUX {
        out.append(if b.kind == BC_T { "    flow::fix_temperature(" } else { "    flow::add_heat_flux(" });
        lit(out, b.v[0]);
        out.append(");\n");
        return;
    }
    if b.kind == BC_D {
        out.append("    flow::wall_distance_region();\n");
        return;
    }
    if b.kind == BC_NT {
        out.append("    flow::fix_spalart(");
        lit(out, b.v[0]);
        out.append(");\n");
        return;
    }
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
    if solid_tok >= 0 { out.append("    flow::clear_solid_bcs();\n"); }
    for b in active[s.first_bc..s.first_bc + s.bc_count] { emit_bc(out, b); }
    out.append(if solid_tok >= 0 { "    solve_coupled(" } else { "    solve(" });
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
    if s.solid_output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, "_", solid_name, ".vtu");
        out.append("    output_solid(");
        quoted(out, path);
        out.append(", ");
        lit(out, reference(LENGTH_D));
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
    out.append(";\nimport solvers from ");
    quoted(out, str(flow_path[..flow_path.len - 4], "solvers"));
    out.append(";\n");
    out.append(PRELUDE);
    if solid_tok >= 0 { out.append(COUPLED_PRELUDE); }
    emit_selects(out);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n");
    if spalart { out.append("    flow::set_spalart(true);\n"); } else if turbulent { out.append("    flow::set_turbulent(true);\n"); }
    if thermal { out.append("    flow::set_thermal(true);\n"); }
    format(out, "    // mesh on ", named[domain].name, ": P", order, "/P1 triangles, written by skein run\n    if !flow::read_mesh(");
    quoted(out, results_dir);
    format(out, ", ", order, ") {\n        write_stderr(\"skein run: cannot read the mesh in ");
    out.append(results_dir);
    out.append("\\n\");\n        exit(3);\n    }\n    tol = flow::tolerance();\n");
    for r in used_regions {
        if named[r].kind != BOUNDARY || !needs_mask(r) { continue; }
        format(out, "    if !read_file(");
        quoted(out, str(results_dir, "/", named[r].name, ".mask"));
        format(out, ", mask_", named[r].name, ") {\n        write_stderr(\"skein run: cannot read the mask of ", named[r].name, "\\n\");\n        exit(3);\n    }\n");
    }
    if solid_tok >= 0 {
        format(out, "    // mesh on ", named[solid_domain].name, ": P2 triangles that share the interface nodes\n    if !flow::read_solid(");
        quoted(out, results_dir);
        out.append(") {\n        write_stderr(\"skein run: cannot read the solid mesh in ");
        out.append(results_dir);
        format(out, "\\n\");\n        exit(3);\n    }\n    // physics ", solid_name, ", material ", name(solid_material), ": plane strain\n    flow::solid(");
        lit(out, scaled_value(property_of(solid_material, "E")));
        out.append(", ");
        lit(out, scalar_value(property_of(solid_material, "nu")));
        out.append(");\n");
        if solid_finite { out.append("    flow::solid_finite();\n"); }
    }
    format(out, "    // physics ", physics_name, ", material ", material_name, "\n    flow::fluid(");
    lit(out, scaled_value(property("rho")));
    out.append(", ");
    lit(out, fluid_mu());
    format(out, ", ", if convective { "true" } else { "false" }, ", flow::", ["NONE", "SUPG", "SUPG_PSPG", "SUPG_PSPG_LSIC"][stabilize], ");\n");
    if spalart {
        out.append("    flow::spalart_constants(");
        for i in 9 {
            if i > 0 { out.append(", "); }
            lit(out, constant(SA_CONSTANTS[i]));
        }
        out.append(");\n    flow::spalart_start(");
        lit(out, v_start[0]);
        out.append(", ");
        lit(out, v_start[1]);
        out.append(", ");
        lit(out, nt_start);
        out.append(");\n");
    }
    if force[0] != 0.0 || force[1] != 0.0 {
        out.append("    flow::body_force(");
        lit(out, force[0]);
        out.append(", ");
        lit(out, force[1]);
        out.append(");\n");
    }
    if thermal {
        format(out, "    // temperature ", temp, ": the energy equation of ", physics_name, "\n    flow::thermal_fluid(");
        lit(out, scaled_value(property("k")));
        out.append(", ");
        lit(out, scaled_value(property("rho")) * scaled_value(property("cp")));
        out.append(", ");
        lit(out, heat_source);
        out.append(");\n    flow::buoyant_force(");
        lit(out, buoyancy[0]);
        out.append(", ");
        lit(out, buoyancy[1]);
        out.append(");\n    flow::temperature_start(");
        lit(out, temp_start);
        out.append(");\n");
    }
    if th_tok >= 0 {
        format(out, "    // physics ", th_name, ", material ", name(th_mat), ": heat conduction in the solid cells\n    flow::thermal_solid(");
        lit(out, scaled_value(property_of(th_mat, "k")));
        out.append(", ");
        lit(out, th_source);
        out.append(");\n");
    }
    if period_high >= 0 {
        let g = named[period_high];
        format(out, "    // symmetry: bloch, periodic\n    region(select_", g.name, "(), ");
        quoted(out, g.name);
        out.append(");\n    if !flow::periodic(");
        lit(out, period[0]);
        out.append(", 0.0) {\n        write_stderr(\"skein run: the nodes of the ends of the period do not match\\n\");\n        exit(3);\n    }\n");
    }
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n    if failed { exit(1); }\n}\n");
}
