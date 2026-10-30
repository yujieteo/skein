// The code generator of the shell models: a Skein syntax tree to a Goose
// program that calls src/shells/.
//
// It takes a model whose physics uses the template reissner_mindlin (spec
// §4.3). The shell is a part of the side of a cylinder on the x axis,
// `boundary(c, y^2 + z^2 == R^2 && <cuts>)`, with cuts that keep one arc of
// θ (θ turns from +z to +y). The mesh is the structured MITC4 grid of
// mesh.goose, with `cell: quad` and order 1. A Curve is a side of the grid,
// for example `boundary(shell, x == 0 m)` or `boundary(shell, y == 0 m)`. A
// Point is the node at it. The section is the laminate of the template
// arguments `layup` and `mirror` and of the material. Each other construct
// gives an error `unsupported` at its token, with a fix.
//
// src/codegen/dispatch.goose has hooks of the shells lane for a model of
// the box generator. A shell model has 6 dofs at each node, a section in
// place of E and nu, and a direct solve (solve.goose), so
// src/codegen/run.goose sends the whole model here, as it does for
// src/elements/gen.goose. It uses the units and the values of
// src/codegen/codegen.goose.
import std;
import math;
namespace shells_gen;
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
fn op(n: i64) -> const u8[:] { syntax::text(syntax::nodes[n].tok) }
fn is_call(n: i64, callee: u8[:]) -> bool {
    kind(n) == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME && name(child(n, 0)) == callee
}
fn set64(into: u8[..64]&, s: u8[:]) { codegen::set64(into, s); }

// The value of argument i of call n.
fn arg(n: i64, i: i64) -> i64 { first(child(n, i + 1)) }

// --- the claim ---------------------------------------------------------------

// Whether this generator takes the model: a physics with a shell or a beam
// template.
fn claims(prog: i64) -> bool {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_PHYSICS {
            var u = first(c);
            while u >= 0 {
                if kind(u) == syntax::N_USES && (name(u) == "reissner_mindlin" || name(u) == "timoshenko") {
                    return true;
                }
                u = next(u);
            }
        }
        c = next(c);
    }
    return false;
}

// --- values --------------------------------------------------------------------

// A scalar (n 1), a vector (n 3) or a 3 × 3 tensor (n 9, row-major), in the
// unit system of the header.
struct Val {
    n: i64,
    v: f64[9],
}

fn scalar(x: f64) -> Val {
    var r = Val { n: 1, v: [0.0; 9] };
    r.v[0] = x;
    return r;
}

fn vector(a: f64[3]) -> Val {
    var r = Val { n: 3, v: [0.0; 9] };
    for j in 3 { r.v[j] = a[j]; }
    return r;
}

fn vec3(v: Val) -> f64[3] { return [v.v[0], v.v[1], v.v[2]]; }

// The names of the `let` and `param` items, of the lets of the physics and
// of the properties of its material, with their value nodes.
var value_names: u8[..64][>..];
var value_nodes: i64[>..];
var let_names: u8[..64][>..];
var let_nodes: i64[>..];
var let_types: i64[>..];
var prop_names: u8[..64][>..];
var prop_nodes: i64[>..];
var depth = 0;

// x, y and z are the point `px` while `at_point`. N and eps_m are the
// tensors `tn` and `te` while `in_physics`.
var at_point = false;
var px: f64[3] = [0.0, 0.0, 0.0];
var in_physics = false;
var tn: f64[9] = [0.0; 9];
var te: f64[9] = [0.0; 9];

fn find(names: u8[..64][:], word: u8[:]) -> i64 {
    for v, i in names {
        if v == word { return i; }
    }
    return -1;
}

let AXIS_NAMES = ["x", "y", "z"];

fn axis_of(word: u8[:]) -> i64 { return position(AXIS_NAMES, word); }

// The factor of an angle unit to rad, or -1 if `unit` is not deg or rad.
fn angle_unit(unit: i64) -> f64 {
    let u = first(unit);
    if kind(u) != syntax::N_UNIT_NAME { return -1.0; }
    if name(u) == "deg" { return PI / 180.0; }
    if name(u) == "rad" { return 1.0; }
    return -1.0;
}

fn unit_factor(unit: i64) -> f64 {
    let a = angle_unit(unit);
    return if a > 0.0 { a } else { codegen::to_system(unit) };
}

recursive fn named(n: i64, word: u8[:]) -> Val {
    if in_physics && (word == "N" || word == "eps_m") {
        var r = Val { n: 9, v: [0.0; 9] };
        for k in 9 { r.v[k] = if word == "N" { tn[k] } else { te[k] }; }
        return r;
    }
    var node = -1;
    if in_physics && find(let_names, word) >= 0 { node = let_nodes[find(let_names, word)]; }
    if node < 0 && in_physics && find(prop_names, word) >= 0 { node = prop_nodes[find(prop_names, word)]; }
    if node < 0 && find(value_names, word) >= 0 { node = value_nodes[find(value_names, word)]; }
    if node >= 0 {
        if depth > 100 {
            unsupported(tok(n), str("the value of '", word, "' depends on itself"), "give it a number");
            return scalar(0.0);
        }
        depth++;
        let r = value(node);
        depth--;
        return r;
    }
    if at_point && axis_of(word) >= 0 { return scalar(px[axis_of(word)]); }
    // n is the outward unit normal of the side of the cylinder.
    if at_point && word == "n" { return vector([0.0, px[1] / radius, px[2] / radius]); }
    if word == "pi" { return scalar(PI); }
    if axis_of(word) >= 0 {
        unsupported(tok(n), str("'", word, "' has no value here"), "use x, y and z only in a value at a point of the shell");
    } else {
        unsupported(tok(n), str("'", word, "' is not a value that this slice knows"),
                    "declare it before, with let or param and a number");
    }
    return scalar(0.0);
}

recursive fn value(n: i64) -> Val {
    let k = kind(n);
    if k == syntax::N_NUMBER && first(n) >= 0 && angle_unit(first(n)) > 0.0 {
        return scalar(codegen::number(tok(n)) * angle_unit(first(n)));
    }
    if k == syntax::N_NUMBER || k == syntax::N_PERCENT { return scalar(codegen::eval(n).v[0]); }
    if k == syntax::N_BOOL { return scalar(if op(n) == "true" { 1.0 } else { 0.0 }); }
    if k == syntax::N_NAME { return named(n, name(n)); }
    if k == syntax::N_VECTOR {
        var r = Val { n: 0, v: [0.0; 9] };
        var f = 1.0;
        var c = first(n);
        while c >= 0 {
            if kind(c) == syntax::N_UNIT {
                f = unit_factor(c);
            } else {
                let x = value(c);
                if x.n != 1 { unsupported(tok(c), "a component of a vector is a vector", "write a scalar"); }
                if r.n < 3 { r.v[r.n] = x.v[0]; }
                r.n++;
            }
            c = next(c);
        }
        if r.n != 3 { unsupported(tok(n), "a vector in this slice has 3 components", "write (x, y, z)"); }
        r.n = 3;
        for j in 3 { r.v[j] *= f; }
        return r;
    }
    if k == syntax::N_UNARY && op(n) == "-" {
        var r = value(first(n));
        for j in 9 { r.v[j] = -r.v[j]; }
        return r;
    }
    if k == syntax::N_MEMBER {
        let o = value(first(n));
        let a = axis_of(name(n));
        if o.n == 3 && a >= 0 { return scalar(o.v[a]); }
    }
    if k == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME {
        let f = name(child(n, 0));
        let count = syntax::child_count(n) - 1;
        if count == 1 && (f == "sin" || f == "cos" || f == "tan" || f == "sqrt" || f == "abs") {
            let x = value(arg(n, 0)).v[0];
            if f == "sin" { return scalar(sin(x)); }
            if f == "cos" { return scalar(cos(x)); }
            if f == "tan" { return scalar(tan(x)); }
            if f == "sqrt" { return scalar(sqrt(x)); }
            return scalar(abs(x));
        }
        if count == 2 && f == "dot" {
            let a = value(arg(n, 0));
            let b = value(arg(n, 1));
            if a.n == 3 && b.n == 3 { return scalar(a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2]); }
        }
    }
    if k == syntax::N_BINARY {
        let o = op(n);
        let a = value(child(n, 0));
        let b = value(child(n, 1));
        if (o == "+" || o == "-") && a.n == b.n {
            var r = a;
            for j in 9 { r.v[j] = if o == "+" { a.v[j] + b.v[j] } else { a.v[j] - b.v[j] }; }
            return r;
        }
        if o == "*" && a.n == 9 && b.n == 3 {
            var r = Val { n: 3, v: [0.0; 9] };
            for i in 3 { for j in 3 { r.v[i] += a.v[3 * i + j] * b.v[j]; } }
            return r;
        }
        if (o == "*" && (a.n == 1 || b.n == 1)) || (o == "/" && b.n == 1) {
            var r = if a.n == 1 && o == "*" { b } else { a };
            let s = if a.n == 1 && o == "*" { a.v[0] } else { b.v[0] };
            for j in 9 { r.v[j] = if o == "*" { r.v[j] * s } else { r.v[j] / s }; }
            return r;
        }
        if o == "^" && a.n == 1 && b.n == 1 { return scalar(pow(a.v[0], b.v[0])); }
    }
    unsupported(tok(n), str("the expression '", source_of(n), "' is not in this slice"),
                "write a number with a unit, a name, a vector (x, y, z), sin, cos, sqrt, abs, dot, or + - * / ^ of these");
    return scalar(0.0);
}

fn scalar_value(n: i64) -> f64 {
    let v = value(n);
    if v.n != 1 { unsupported(tok(n), "the value is not a scalar", "write a scalar"); }
    return v.v[0];
}

fn value_at(n: i64, p: f64[3]) -> Val {
    at_point = true;
    px = p;
    let v = value(n);
    at_point = false;
    return v;
}

// --- the model -------------------------------------------------------------------

// A region of the shell: sides of the grid (bit 0 is x = x0, 1 is x = x1,
// 2 is θ = t0, 3 is θ = t1) or the node at a point.
struct Region {
    name: u8[..64],
    sides: i64,     // 0 for a point
    at: f64[3],
}

// The bc kinds: a fixed dof, or a force per length `N*n` on a side.
let BC_FIX = 0;
let BC_EDGE = 1;

struct Bc {
    region: i64,
    kind: i64,
    dof: i64,       // 0 to 5, for BC_FIX
    load: f64[3],
}

// What a report or an expect gives: a component of the displacement or of
// the rotation at a point, a let of the physics at a point, or the reaction
// on a region.
let F_W = 0;
let F_THETA = 1;
let F_LET = 2;
let F_REACTION = 3;

struct Field {
    what: i64,
    comp: i64,      // 0 to 2, or -1 for the vector
    let_index: i64, // for F_LET
    region: i64,    // for F_REACTION
    at: f64[3],
}

struct Report {
    field: Field,
    label: u8[..256],
}

struct Step {
    name: u8[..64],
    first_own: i64,     // the bcs of the step: a range of `bcs`
    own_count: i64,
    first_report: i64,
    report_count: i64,
    output: bool,
}

struct Expect {
    name: u8[..64],
    field: Field,
    want: f64,
    within: f64,
    relative: bool,
    label: u8[..256],
}

var solid_names: u8[..64][>..];
var solid_nodes: i64[>..];
var table_names: u8[..64][>..];
var table_nodes: i64[>..];
var material_names: u8[..64][>..];
var material_nodes: i64[>..];
var regions: Region[>..];
var bcs: Bc[>..];
var global_bcs: i64[>..];
var reports: Report[>..];
var steps: Step[>..];
var expects: Expect[>..];

// The shell: the side of a cylinder of radius r about the x axis, from x0
// to x1 and from t0 to t1 in θ.
var shell_name: u8[..64];
var shell_tok = -1;
var radius = 0.0;
var x0 = 0.0;
var x1 = 0.0;
var t0 = 0.0;
var t1 = 0.0;
var mesh_size = 0.0;
var mesh_tok = -1;
var physics_tok = -1;
var physics_name: u8[..64];
var w_name: u8[..64];
var theta_name: u8[..64];
var material_name: u8[..64];
var args: i64[5] = [-1; 5];     // the nodes of layup, mirror, axis, top and q, or -1

let TEMPLATE_ARGS = ["layup", "mirror", "axis", "top", "q"];

fn reset() {
    codegen::reset();
    value_names.clear();
    value_nodes.clear();
    let_names.clear();
    let_nodes.clear();
    let_types.clear();
    prop_names.clear();
    prop_nodes.clear();
    depth = 0;
    at_point = false;
    in_physics = false;
    solid_names.clear();
    solid_nodes.clear();
    table_names.clear();
    table_nodes.clear();
    material_names.clear();
    material_nodes.clear();
    regions.clear();
    bcs.clear();
    global_bcs.clear();
    reports.clear();
    steps.clear();
    expects.clear();
    shell_name = [];
    shell_tok = -1;
    mesh_size = 0.0;
    mesh_tok = -1;
    physics_tok = -1;
    physics_name = [];
    w_name = [];
    theta_name = [];
    material_name = [];
    args = [-1; 5];
}

// The point of the shell at x and θ.
fn shell_point(x: f64, t: f64) -> f64[3] { return [x, radius * sin(t), radius * cos(t)]; }

let SHELL_FIX = "write the shell as boundary(c, y^2 + z^2 == R^2 && z >= 0 m), with c = cylinder((0, 0, 0) m, (1, 0, 0), R, L)";

// The comparisons of `cond`, a && b && ..., into `out`. False if a part is
// not a comparison.
recursive fn terms(cond: i64, out: i64[>..]&) -> bool {
    if kind(cond) != syntax::N_BINARY { return false; }
    let o = op(cond);
    if o == "&&" { return terms(child(cond, 0), out) && terms(child(cond, 1), out); }
    if o == "==" || o == ">=" || o == "<=" || o == ">" || o == "<" {
        out.push(cond);
        return true;
    }
    return false;
}

// The margin of comparison `term` at p, >= 0 where it holds, and the size
// of its sides. For == the margin is -|a - b|.
fn margin(term: i64, p: f64[3]) -> f64, f64 {
    let a = value_at(child(term, 0), p).v[0];
    let b = value_at(child(term, 1), p).v[0];
    let o = op(term);
    let size = abs(a) + abs(b);
    if o == ">=" || o == ">" { return a - b, size; }
    if o == "<=" || o == "<" { return b - a, size; }
    return -abs(a - b), size;
}

// Whether the cuts (the terms other than ==) hold at θ, at mid-length.
fn cuts_hold(cuts: i64[:], t: f64) -> bool {
    for c in cuts {
        let m, size = margin(c, shell_point(0.5 * (x0 + x1), t));
        if m < 0.0 { return false; }
    }
    return true;
}

// The θ between a (where the cuts do not hold) and b (where they hold) at
// which the cuts start to hold.
fn edge_of(cuts: i64[:], a: f64, b: f64) -> f64 {
    var lo = a;
    var hi = b;
    for k in 80 {
        let m = 0.5 * (lo + hi);
        if cuts_hold(cuts, m) { hi = m; } else { lo = m; }
    }
    return hi;
}

// `let roof: Surface = boundary(c, y^2 + z^2 == R^2 && z >= ...)`.
fn shell(n: i64, v: i64) {
    if shell_tok >= 0 {
        unsupported(tok(n), str("a second shell surface: this slice meshes ", shell_name, " only"),
                    "put all the shell on one surface");
        return;
    }
    guard is_call(v, "boundary") && syntax::child_count(v) == 3 && kind(arg(v, 0)) == syntax::N_NAME else {
        unsupported(tok(v), str("the shell '", name(n), "' is not the side of a cylinder"), SHELL_FIX);
        return;
    }
    let si = find(solid_names, name(arg(v, 0)));
    let c = if si >= 0 { solid_nodes[si] } else { -1 };
    guard c >= 0 && is_call(c, "cylinder") && syntax::child_count(c) == 5 else {
        unsupported(tok(arg(v, 0)), str("the shell '", name(n), "' is not on a cylinder that a let names"), SHELL_FIX);
        return;
    }
    let centre = vec3(value(arg(c, 0)));
    let axis = vec3(value(arg(c, 1)));
    radius = scalar_value(arg(c, 2));
    let length = scalar_value(arg(c, 3));
    let l = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if !ok() { return; }
    guard l > 0.0 && abs(axis[1]) <= 1.0e-12 * l && abs(axis[2]) <= 1.0e-12 * l && axis[0] > 0.0 &&
          centre[1] == 0.0 && centre[2] == 0.0 && radius > 0.0 && length > 0.0 else {
        unsupported(tok(c), "this slice has shells on a cylinder on the x axis only", SHELL_FIX);
        return;
    }
    x0 = centre[0];
    x1 = centre[0] + length;
    var parts: i64[>..] = [];
    guard terms(arg(v, 1), parts) else {
        unsupported(tok(arg(v, 1)), "the condition is not comparisons joined by &&", SHELL_FIX);
        return;
    }
    var cuts: i64[>..] = [];
    var side = false;
    for p in parts {
        if op(p) == "==" {
            // An == must hold on the side of the cylinder.
            var on = true;
            for t in [0.3, 1.7, -2.1] {
                for x in [x0, x1] {
                    let m, size = margin(p, shell_point(x, t));
                    if -m > 1.0e-9 * size { on = false; }
                }
            }
            if !on {
                unsupported(tok(p), "the condition is not on the side of the cylinder", SHELL_FIX);
                return;
            }
            side = true;
        } else {
            // A cut must not change along x.
            for t in [0.3, 1.7, -2.1, 2.9] {
                let a, sa = margin(p, shell_point(x0, t));
                let b, sb = margin(p, shell_point(x1, t));
                if (a >= 0.0) != (b >= 0.0) {
                    unsupported(tok(p), "this slice cuts the cylinder in θ only", "cut the length with the cylinder itself");
                    return;
                }
            }
            cuts.push(p);
        }
        if !ok() { return; }
    }
    if !side {
        unsupported(tok(arg(v, 1)), "the condition has no y^2 + z^2 == R^2", SHELL_FIX);
        return;
    }
    // The arc of θ where the cuts hold, from a scan of -π to π.
    let count = 1440;
    var holds: bool[>..] = [];
    var inside = 0;
    for k in count {
        holds.push(cuts_hold(cuts, -PI + 2.0 * PI * (k + 0.5) / count));
        if holds[k] { inside++; }
    }
    if inside == 0 || inside == count {
        unsupported(tok(arg(v, 1)), if inside == 0 { "the shell has no area" } else { "this slice has no closed tube" },
                    "cut the cylinder to an arc, for example with z >= 0 m and a plane of symmetry");
        return;
    }
    var starts = 0;
    var s = 0;
    for k in count {
        if holds[k] && !holds[(k + count - 1) % count] {
            starts++;
            s = k;
        }
    }
    if starts != 1 {
        unsupported(tok(arg(v, 1)), "the cuts keep more than one arc", "keep one arc of the cylinder");
        return;
    }
    let step = 2.0 * PI / count;
    let ts = -PI + step * (s + 0.5);
    t0 = edge_of(cuts, ts - step, ts);
    let last = ts + step * (inside - 1);
    t1 = edge_of(cuts, last + step, last);
    set64(shell_name, name(n));
    shell_tok = tok(n);
}

let REGION_FIX = "write boundary(shell, x == 0 m) for a side of the shell, or point((x, y, z) m)";

// `let end0: Curve = boundary(roof, x == 0 mm)`: the sides of the grid
// where the == holds, or `let tip: Point = point((x, y, z) mm)`.
fn region(n: i64, v: i64, curve: bool) {
    var r = Region { name: [], sides: 0, at: [0.0, 0.0, 0.0] };
    set64(r.name, name(n));
    if !curve {
        guard is_call(v, "point") && syntax::child_count(v) == 2 else {
            unsupported(tok(v), str("the point '", name(n), "' is not point((x, y, z))"), REGION_FIX);
            return;
        }
        r.at = vec3(value(arg(v, 0)));
        regions.push(r);
        return;
    }
    guard is_call(v, "boundary") && syntax::child_count(v) == 3 && kind(arg(v, 0)) == syntax::N_NAME &&
          shell_tok >= 0 && name(arg(v, 0)) == shell_name else {
        unsupported(tok(v), str("the curve '", name(n), "' is not a side of the shell"), REGION_FIX);
        return;
    }
    let cond = arg(v, 1);
    guard kind(cond) == syntax::N_BINARY && op(cond) == "==" else {
        unsupported(tok(cond), "the side needs one ==", REGION_FIX);
        return;
    }
    let xm = 0.5 * (x0 + x1);
    let tm = 0.5 * (t0 + t1);
    // The size of the sides at the middle of the shell, for the tolerance.
    let mm, scale = margin(cond, shell_point(xm, tm));
    for side in 4 {
        var on = true;
        for f in [0.0, 0.5, 1.0] {
            let x = if side == 0 { x0 } else if side == 1 { x1 } else { x0 + f * (x1 - x0) };
            let t = if side == 2 { t0 } else if side == 3 { t1 } else { t0 + f * (t1 - t0) };
            let m, size = margin(cond, shell_point(x, t));
            if -m > 1.0e-9 * (size + scale) { on = false; }
        }
        if on { r.sides |= 1 << side; }
    }
    if r.sides != 0 {
        regions.push(r);
        return;
    }
    if ok() { unsupported(tok(cond), str("the curve '", name(n), "' is not a side of the shell"), REGION_FIX); }
}

fn find_region(n: i64) -> i64 {
    if kind(n) != syntax::N_NAME { return -1; }
    for r, i in regions {
        if r.name == name(n) { return i; }
    }
    return -1;
}

fn value_item(n: i64) {
    let ty = child(n, 0);
    let v = child(n, 1);
    let t = if kind(ty) == syntax::N_TYPE_NAME { name(ty) } else { "" };
    var w: u8[..64] = [];
    set64(w, name(n));
    if t == "Solid" {
        solid_names.push(w);
        solid_nodes.push(v);
    } else if t == "Surface" {
        shell(n, v);
    } else if t == "Curve" || t == "Point" {
        region(n, v, t == "Curve");
    } else if t == "Volume" || t == "Region" {
        unsupported(tok(n), str("a ", t, " is not in this slice"), "use a Surface for the shell");
    } else {
        value_names.push(w);
        value_nodes.push(v);
    }
}

fn mesh(n: i64) {
    if mesh_tok >= 0 {
        unsupported(tok(n), "a second mesh: this slice has one", "remove the second mesh");
        return;
    }
    mesh_tok = tok(n);
    let r = child(n, 0);
    if kind(r) != syntax::N_NAME || name(r) != shell_name {
        unsupported(tok(r), str("the mesh is not on the shell ", shell_name), str("write mesh m on ", shell_name));
    }
    var s = first(child(n, 1));
    while s >= 0 {
        let v = first(s);
        if kind(s) == syntax::N_SETTING && name(s) == "size" {
            mesh_size = scalar_value(v);
        } else if kind(s) == syntax::N_SETTING && name(s) == "cell" {
            if kind(v) != syntax::N_NAME || name(v) != "quad" {
                unsupported(tok(v), "this slice meshes the shell with quad cells only", "write cell: quad");
            }
        } else {
            unsupported(tok(s), "the mesh setting is not in this slice", "use only size and cell: quad");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 50 mm"); }
}

let PHYSICS_FIX = "write physics shell on roof: reissner_mindlin(top: n) solves w: Vector<Length>, theta: Vector<Angle> { material ply }";

fn physics(n: i64) {
    if physics_tok >= 0 {
        unsupported(tok(n), "a second physics: this slice has one shell physics", "remove the second physics");
        return;
    }
    physics_tok = tok(n);
    set64(physics_name, name(n));
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            let r = first(c);
            if kind(r) != syntax::N_NAME || name(r) != shell_name {
                unsupported(tok(r), str("the physics is not on the shell ", shell_name), str("put it on ", shell_name));
            }
        } else if k == syntax::N_USES && name(c) == "reissner_mindlin" {
            var a = first(c);
            var i = 0;
            while a >= 0 {
                let at = if syntax::nodes[a].aux >= 0 { position(TEMPLATE_ARGS, syntax::text(syntax::nodes[a].aux)) } else { i };
                if at >= 0 && at < 5 { args[at] = first(a); }
                i++;
                a = next(a);
            }
        } else if k == syntax::N_USES {
            unsupported(tok(c), str("skein run has no generator for ", name(c), " yet"),
                        "run the beams with the shell lane tests (src/shells/test.goose)");
        } else if k == syntax::N_UNKNOWN {
            if w_name.len == 0 { set64(w_name, name(c)); } else if theta_name.len == 0 { set64(theta_name, name(c)); }
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            if syntax::child_count(c) != 1 || first(m) >= 0 || material_name.len > 0 {
                unsupported(tok(c), "this slice gives one material to the whole shell", PHYSICS_FIX);
            }
            set64(material_name, name(m));
        } else if k == syntax::N_LET {
            var w: u8[..64] = [];
            set64(w, name(c));
            let_names.push(w);
            let_nodes.push(child(c, 1));
            let_types.push(child(c, 0));
        } else {
            unsupported(tok(c), "the physics statement is not in this slice", "use only material and let");
        }
        c = next(c);
    }
    if args[3] < 0 { unsupported(tok(n), "reissner_mindlin needs top", "write top: n"); }
}

fn method(n: i64) {
    if name(n) != physics_name || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this slice has only the cg method of the physics", "write method shell: cg { order: 1 }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            if scalar_value(first(s)) != 1.0 {
                unsupported(tok(s), "the MITC4 shell is of order 1", "write order: 1");
            }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this slice", "use only mesh and order");
        }
        s = next(s);
    }
}

let BC_FIX_TEXT = "write bc end: w = 0, bc end: theta.x = 0, or bc end: N*n = (1, 0, 0) N/mm";

// The dof of `w`, `w.x`, `theta` or `theta.z`: the first dof and the
// count, or -1.
fn target_dofs(t: i64) -> i64, i64 {
    if kind(t) == syntax::N_NAME && name(t) == w_name { return 0, 3; }
    if kind(t) == syntax::N_NAME && name(t) == theta_name { return 3, 3; }
    if kind(t) == syntax::N_MEMBER && kind(first(t)) == syntax::N_NAME && axis_of(name(t)) >= 0 {
        if name(first(t)) == w_name { return axis_of(name(t)), 1; }
        if name(first(t)) == theta_name { return 3 + axis_of(name(t)), 1; }
    }
    return -1, 0;
}

// A `bc` item or statement: one Bc for each fixed dof, or one load.
fn bc(n: i64) {
    let r = find_region(child(n, 0));
    let target = child(n, 1);
    let v = child(n, 2);
    if r < 0 {
        unsupported(tok(child(n, 0)), "the region of the bc is not a Curve or a Point that a let names", REGION_FIX);
        return;
    }
    let d, count = target_dofs(target);
    if d >= 0 {
        let x = value(v);
        for j in 9 {
            if x.v[j] != 0.0 { unsupported(tok(v), "this slice fixes a value of 0 only", BC_FIX_TEXT); }
        }
        if x.n != 1 && (count == 1 || x.n != 3) { unsupported(tok(v), "the value has the wrong size", BC_FIX_TEXT); }
        for j in count {
            bcs.push(Bc { region: r, kind: BC_FIX, dof: d + j, load: [0.0, 0.0, 0.0] });
            if !in_step { global_bcs.push(bcs.len - 1); }
        }
        return;
    }
    if source_of(target) == "N*n" && regions[r].sides != 0 {
        let x = value(v);
        if x.n != 3 { unsupported(tok(v), "the force per length is a vector", BC_FIX_TEXT); }
        bcs.push(Bc { region: r, kind: BC_EDGE, dof: -1, load: vec3(x) });
        if !in_step { global_bcs.push(bcs.len - 1); }
        return;
    }
    unsupported(tok(target), str("the bc target '", source_of(target), "' is not in this slice"), BC_FIX_TEXT);
}

var in_step = false;

let REPORT_FIX = "write report w.z at (x, y, z) mm, report s at (x, y, z) mm with s a let of the physics, or report reaction(end)";

// The field of a report or an expect.
fn field(e: i64, one: bool) -> Field, bool {
    var f = Field { what: F_W, comp: -1, let_index: -1, region: -1, at: [0.0, 0.0, 0.0] };
    var left = e;
    if kind(e) == syntax::N_MEMBER && is_call(first(e), "reaction") {
        f.comp = axis_of(name(e));
        left = first(e);
    }
    if is_call(left, "reaction") && syntax::child_count(left) == 2 {
        f.what = F_REACTION;
        f.region = find_region(arg(left, 0));
        if f.region < 0 {
            unsupported(tok(left), "reaction needs a Curve or a Point that a let names", REPORT_FIX);
            return f, false;
        }
        return f, !one || f.comp >= 0;
    }
    guard kind(e) == syntax::N_AT else {
        unsupported(tok(e), str("the value '", source_of(e), "' is not in this slice"), REPORT_FIX);
        return f, false;
    }
    let v = child(e, 0);
    f.at = vec3(value(child(e, 1)));
    let d, count = target_dofs(v);
    if d >= 0 {
        f.what = if d < 3 { F_W } else { F_THETA };
        f.comp = if count == 1 { d % 3 } else { -1 };
        return f, !one || f.comp >= 0;
    }
    if kind(v) == syntax::N_NAME && find(let_names, name(v)) >= 0 {
        f.what = F_LET;
        f.let_index = find(let_names, name(v));
        return f, true;
    }
    unsupported(tok(v), str("'", source_of(v), "' at a point is not in this slice"), REPORT_FIX);
    return f, false;
}

fn step(n: i64) {
    let proc = child(n, 0);
    var s = Step { name: [], first_own: bcs.len, own_count: 0, first_report: reports.len, report_count: 0, output: false };
    set64(s.name, name(n));
    if name(proc) != "static" || first(proc) >= 0 {
        unsupported(tok(proc), "this slice has only static steps with no arguments", "write step load: static { ... }");
    }
    in_step = true;
    var c = next(proc);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            bc(c);
        } else if k == syntax::N_OUTPUT {
            s.output = true;
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                let f, good = field(first(r), false);
                if good {
                    var rep = Report { field: f, label: [] };
                    codegen::set256(rep.label, source_of(r));
                    reports.push(rep);
                }
                r = next(r);
            }
        } else {
            unsupported(tok(c), "the step statement is not in this slice", "use only bc, output and report");
        }
        c = next(c);
    }
    in_step = false;
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

fn expect(n: i64) {
    let f, good = field(child(n, 0), true);
    if !good {
        if ok() { unsupported(tok(child(n, 0)), "an expect compares one component", "write w.z, or reaction(end).x"); }
        return;
    }
    var x = Expect { name: [], field: f, want: scalar_value(child(n, 1)), within: scalar_value(child(n, 2)),
                     relative: kind(child(n, 2)) == syntax::N_PERCENT, label: [] };
    set64(x.name, name(n));
    codegen::set256(x.label, str(source_of(child(n, 1)), " within ", source_of(child(n, 2))));
    expects.push(x);
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// the errors.
fn read(prog: i64) -> bool {
    reset();
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            codegen::header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
            if !ok() { return false; }
        } else if k == syntax::N_TABLE {
            var w: u8[..64] = [];
            set64(w, name(c));
            table_names.push(w);
            table_nodes.push(c);
        } else if k == syntax::N_MATERIAL {
            var w: u8[..64] = [];
            set64(w, name(c));
            material_names.push(w);
            material_nodes.push(c);
        } else if k == syntax::N_MESH {
            mesh(c);
        } else if k == syntax::N_PHYSICS {
            physics(c);
        } else if k == syntax::N_METHOD {
            method(c);
        } else if k == syntax::N_BC {
            bc(c);
        } else if k == syntax::N_STEP {
            step(c);
        } else if k == syntax::N_EXPECT {
            expect(c);
        } else if k != syntax::N_SCALES {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this slice"),
                        "remove it: this slice runs one linear static shell");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if shell_tok < 0 { unsupported(at, "the model has no shell surface", SHELL_FIX); }
    if mesh_tok < 0 { unsupported(at, "the model has no mesh", "add mesh m on roof { size: 1000 mm, cell: quad }"); }
    if physics_tok < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { ... }"); }
    let m = find(material_names, material_name);
    if physics_tok >= 0 && m < 0 { unsupported(physics_tok, "the physics has no material that the model declares", PHYSICS_FIX); }
    if ok() { read_section(material_nodes[m]); }
    if ok() { read_loads(); }
    for r in reports { if ok() { check_let(r.field); } }
    for x in expects { if ok() { check_let(x.field); } }
    return ok();
}

// A let of the physics in a report or an expect must be linear in N and
// eps_m (let_coefficients).
fn check_let(f: Field) {
    if f.what != F_LET { return; }
    let c, linear = let_coefficients(f.let_index, f.at);
    if !linear {
        unsupported(physics_tok, str("the let ", let_names[f.let_index], " is not linear in N and eps_m"),
                    "write the let as dot(a, N*b)/h, with a and b vectors");
    }
}

// --- the section and the loads ----------------------------------------------------

var e1 = 0.0;
var e2 = 0.0;
var g12 = 0.0;
var nu12 = 0.0;
var t_ply = 0.0;
var angles: f64[>..] = [];     // the angle of each ply from the top face
var axis: f64[3] = [1.0, 0.0, 0.0];
// The load per area: a vector, or a pressure along the top normal.
var q: f64[3] = [0.0, 0.0, 0.0];
var q_normal = false;

let MATERIAL_FIX = "give the material E1, E2, G12, nu12, t_ply and h";

// The plies of `layup` (an angle, or a table column of angles) and the
// properties of the material: the plies repeat to fill h (h/2 with mirror),
// and mirror adds them again in the reverse sequence (spec §4.3).
fn read_section(m: i64) {
    in_physics = true;
    var p = first(m);
    while p >= 0 {
        var w: u8[..64] = [];
        set64(w, name(p));
        prop_names.push(w);
        prop_nodes.push(child(p, 1));
        p = next(p);
    }
    let props = ["E1", "E2", "G12", "nu12", "t_ply", "h"];
    var v: f64[6] = [0.0; 6];
    for word, i in props {
        let k = find(prop_names, word);
        if k < 0 {
            unsupported(tok(m), str("the material has no ", word), MATERIAL_FIX);
            in_physics = false;
            return;
        }
        v[i] = scalar_value(prop_nodes[k]);
    }
    e1 = v[0];
    e2 = v[1];
    g12 = v[2];
    nu12 = v[3];
    t_ply = v[4];
    let h = v[5];
    let mirror = args[1] >= 0 && scalar_value(args[1]) != 0.0;
    var base: f64[>..] = [];
    let layup = args[0];
    if layup < 0 {
        base.push(0.0);
    } else if kind(layup) == syntax::N_MEMBER && kind(first(layup)) == syntax::N_NAME &&
              find(table_names, name(first(layup))) >= 0 {
        table_column(table_nodes[find(table_names, name(first(layup)))], name(layup), tok(layup), base);
    } else {
        base.push(scalar_value(layup));
    }
    let fill = if mirror { 0.5 * h } else { h };
    let count = round(fill / max(t_ply, 1.0e-300)) as! i64;
    if !(e1 > 0.0 && e2 > 0.0 && g12 > 0.0 && t_ply > 0.0) || count < 1 || abs(count * t_ply - fill) > 1.0e-6 * fill {
        unsupported(tok(m), "h is not a whole number of plies of t_ply (2 t_ply with mirror)", MATERIAL_FIX);
    }
    angles.clear();
    if base.len > 0 {
        for k in max(count, 0) { angles.push(base[k % base.len]); }
    }
    if mirror {
        let c = angles.len;
        for k in c { angles.push(angles[c - 1 - k]); }
    }
    in_physics = false;
}

// The values of column `word` of table t, in the unit of the column.
fn table_column(t: i64, word: u8[:], at: i64, out: f64[>..]&) {
    var col = -1;
    var factor = 1.0;
    var i = 0;
    var c = first(t);
    while c >= 0 {
        if kind(c) == syntax::N_COLUMN {
            if name(c) == word {
                col = i;
                let u = syntax::child_of_kind(c, syntax::N_UNIT);
                if u >= 0 { factor = unit_factor(u); }
            }
            i++;
        } else if kind(c) == syntax::N_ROW && col >= 0 {
            out.push(scalar_value(child(c, col)) * factor);
        }
        c = next(c);
    }
    if col < 0 { unsupported(at, str("the table has no column ", word), "name a column of the table"); }
}

// The section axis, the top normal and the load per area, from the template
// arguments. Each is checked at points of the shell.
fn read_loads() {
    let xm = 0.5 * (x0 + x1);
    let samples = [shell_point(xm, t0 + 0.2 * (t1 - t0)), shell_point(x0 + 0.3 * (x1 - x0), 0.5 * (t0 + t1)),
                   shell_point(x1 - 0.1 * (x1 - x0), t1 - 0.3 * (t1 - t0))];
    if args[2] >= 0 { axis = vec3(value(args[2])); } else { axis = [1.0, 0.0, 0.0]; }
    let top = args[3];
    if kind(top) != syntax::N_NAME || name(top) != "n" {
        for p in samples {
            let v = value_at(top, p);
            if v.n != 3 || v.v[1] * p[1] + v.v[2] * p[2] <= 0.0 {
                unsupported(tok(top), "this slice has the top face outward", "write top: n, or a vector outward");
                return;
            }
        }
    }
    q = [0.0, 0.0, 0.0];
    q_normal = false;
    if args[4] < 0 { return; }
    var same = true;
    var radial = true;
    let a = vec3(value_at(args[4], samples[0]));
    let s = (a[1] * samples[0][1] + a[2] * samples[0][2]) / radius;
    for p in samples {
        let b = vec3(value_at(args[4], p));
        for j in 3 {
            let er = if j == 0 { 0.0 } else { p[j] / radius };
            if abs(b[j] - a[j]) > 1.0e-12 * (abs(a[0]) + abs(a[1]) + abs(a[2])) { same = false; }
            if abs(b[j] - s * er) > 1.0e-12 * abs(s) { radial = false; }
        }
    }
    if !ok() { return; }
    if same {
        q = a;
    } else if radial {
        q = [s, 0.0, 0.0];
        q_normal = true;
    } else {
        unsupported(tok(args[4]), "this slice has a load per area that is uniform or along the radius",
                    "write q: (0, 0, -1) kPa or q: p*n with n the outward normal");
    }
}

// The value of let i of the physics at p as c[0] + c[1..10] · N +
// c[10..19] · eps_m, from evaluations at the unit tensors. False if the
// value is not linear in N and eps_m.
fn let_coefficients(i: i64, p: f64[3]) -> f64[19], bool {
    var c: f64[19] = [0.0; 19];
    in_physics = true;
    tn = [0.0; 9];
    te = [0.0; 9];
    let v0 = value_at(let_nodes[i], p);
    c[0] = v0.v[0];
    for k in 18 {
        tn = [0.0; 9];
        te = [0.0; 9];
        if k < 9 { tn[k] = 1.0; } else { te[k - 9] = 1.0; }
        c[k + 1] = value_at(let_nodes[i], p).v[0] - c[0];
    }
    // A check with all the entries at once.
    var want = c[0];
    for k in 9 {
        tn[k] = 0.7 + 0.3 * k;
        te[k] = 1.3 - 0.2 * k;
        want += c[k + 1] * tn[k] + c[k + 10] * te[k];
    }
    let got = value_at(let_nodes[i], p).v[0];
    in_physics = false;
    var size = abs(c[0]);
    for k in 18 { size += abs(c[k + 1]); }
    return c, v0.n == 1 && abs(got - want) <= 1.0e-9 * size;
}

// --- the program -----------------------------------------------------------------

fn lit(out: u8[>..]&, x: f64) { codegen::lit(out, x); }

fn lits(out: u8[>..]&, v: f64[3]) {
    lit(out, v[0]);
    out.append(", ");
    lit(out, v[1]);
    out.append(", ");
    lit(out, v[2]);
}

fn quoted(out: u8[>..]&, s: u8[:]) { codegen::quoted(out, s); }

// Appends the unit of a declared type.
recursive fn type_unit(t: i64, out: u8[>..]&) {
    let sys = codegen::system;
    let k = kind(t);
    if k == syntax::N_TYPE_NAME {
        let w = name(t);
        if w == "Length" { out.append(sys.length); }
        if w == "Pressure" { out.append(sys.pressure); }
        if w == "Force" { out.append(sys.force); }
        if w == "Angle" { out.append("rad"); }
    } else if k == syntax::N_DIM_DIV || k == syntax::N_DIM_MUL {
        type_unit(child(t, 0), out);
        out.append(if k == syntax::N_DIM_DIV { "/" } else { "*" });
        type_unit(child(t, 1), out);
    }
}

fn field_unit(f: Field) -> u8[>..] {
    if f.what == F_W { return str(codegen::system.length); }
    if f.what == F_THETA { return str("rad"); }
    if f.what == F_REACTION { return str(codegen::system.force); }
    var u: u8[>..] = [];
    type_unit(let_types[f.let_index], u);
    return u;
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
    var grid = -1;

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

    // The nodes of the sides of the grid in the bits of `mask`, each once.
    fn sides(mask: i64) -> i64[>..] {
        var out: i64[>..] = [];
        for k in 4 {
            if mask & (1 << k) == 0 { continue; }
            for n in shells::grid_side(grid, k) { if position(out[..], n) < 0 { out.push(n); } }
        }
        return out;
    }

    fn fix_sides(mask: i64, dof: i64) {
        for n in sides(mask) { shells::fix(n, dof); }
    }

    fn edge_load(mask: i64, fx: f64, fy: f64, fz: f64) {
        for k in 4 {
            if mask & (1 << k) == 0 { continue; }
            let s = shells::grid_side(grid, k);
            for i in s.len - 1 { shells::edge_load(s[i], s[i + 1], fx, fy, fz); }
        }
    }

    // The node at (x, y, z), for a bc on a point. A point that is not a
    // node of the mesh is a fault of the model.
    fn node_at(x: f64, y: f64, z: f64, tol: f64, expr: u8[:]) -> i64 {
        let n, d = shells::nearest_node(x, y, z);
        if d > tol {
            write_stderr(str("skein run: the point of '", expr, "' is not a node of the mesh\n"));
            failed = true;
        }
        return n;
    }

    // The cell of the grid that holds (x, y, z), and the place of the point
    // in it. A point that is not on the shell is a fault of the model.
    fn locate(x: f64, y: f64, z: f64, tol: f64, expr: u8[:]) -> i64, i64, f64, f64 {
        let i, j, fx, ft, found = shells::grid_locate(grid, x, y, z, tol);
        if !found {
            write_stderr(str("skein run: the point of '", expr, "' is not on the shell\n"));
            failed = true;
        }
        return i, j, fx, ft;
    }

    // Dof d at (x, y, z), bilinear in the cell that holds the point.
    fn dof_at(x: f64, y: f64, z: f64, tol: f64, expr: u8[:], d: i64) -> f64 {
        let i, j, fx, ft = locate(x, y, z, tol, expr);
        return (1.0 - fx) * (1.0 - ft) * shells::displacement(shells::grid_node(grid, i, j), d)
            + fx * (1.0 - ft) * shells::displacement(shells::grid_node(grid, i + 1, j), d)
            + fx * ft * shells::displacement(shells::grid_node(grid, i + 1, j + 1), d)
            + (1.0 - fx) * ft * shells::displacement(shells::grid_node(grid, i, j + 1), d);
    }

    // N and eps_m at (x, y, z): the mean of the cells at a node, or the
    // value at the centre of the cell that holds the point.
    fn resultants_at(x: f64, y: f64, z: f64, tol: f64, expr: u8[:]) -> f64[9], f64[9] {
        let n, d = shells::nearest_node(x, y, z);
        if d <= tol { return shells::mean_resultants(n); }
        let i, j, fx, ft = locate(x, y, z, tol, expr);
        return shells::shell_resultants(shells::grid_cell(grid, i, j));
    }

    // The direct solve has no iterations: cg_iterations is 0.
    fn solve(step: u8[:]) {
        let converged = shells::solve();
        format(solves, "  ", step, ",", shells::shells.len, ",", shells::dof_count(), ",0,",
               shells::last_seconds, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of the sum of the sizes of the loads and the
    // reactions.
    fn equilibrium(unit: u8[:]) {
        let r = shells::reactions();
        var nodes: i64[>..] = [];
        for n in shells::node_count() { nodes.push(n); }
        let rx, ry, rz = shells::reaction_sum(r, nodes);
        let fx, fy, fz = shells::load_sum();
        let e = sqrt((rx + fx) * (rx + fx) + (ry + fy) * (ry + fy) + (rz + fz) * (rz + fz));
        let limit = 1.0e-6 * shells::load_scale(r);
        check("equilibrium", e <= limit, quantity(e, unit), quantity(limit, unit));
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let tol = if relative { within * abs(want) } else { within };
        check("expect", abs(got - want) <= tol, str(name, ": ", got, " ", unit), limit);
    }

    fn output(path: u8[:], w: u8[:], theta: u8[:]) {
        if !shells::write_vtu(path, w, theta) {
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

let DOF_NAMES = ["shells::UX", "shells::UY", "shells::UZ", "shells::RX", "shells::RY", "shells::RZ"];

fn emit_node_at(out: u8[>..]&, var_name: u8[:], at: f64[3], expr: u8[:]) {
    format(out, "        let ", var_name, " = node_at(");
    lits(out, at);
    out.append(", ");
    lit(out, 1.0e-3 * mesh_size);
    out.append(", ");
    quoted(out, expr);
    out.append(");\n");
}

fn emit_bc(out: u8[>..]&, b: Bc) {
    let r = regions[b.region];
    format(out, "    // bc ", r.name, "\n");
    if b.kind == BC_EDGE {
        format(out, "    edge_load(", r.sides, ", ");
        lits(out, b.load);
        out.append(");\n");
    } else if r.sides != 0 {
        format(out, "    fix_sides(", r.sides, ", ", DOF_NAMES[b.dof], ");\n");
    } else {
        out.append("    {\n");
        emit_node_at(out, "n", r.at, r.name);
        format(out, "        shells::fix(n, ", DOF_NAMES[b.dof], ");\n    }\n");
    }
}

// `let v = <the value of field f>;`, in a block.
fn emit_value(out: u8[>..]&, f: Field, label: u8[:], comp: i64) {
    if f.what == F_REACTION {
        let r = regions[f.region];
        out.append("        let r = shells::reactions();\n");
        if r.sides != 0 {
            format(out, "        let vx, vy, vz = shells::reaction_sum(r, sides(", r.sides, "));\n");
        } else {
            emit_node_at(out, "n", r.at, label);
            out.append("        let vx, vy, vz = shells::reaction_sum(r, [n]);\n");
        }
        return;
    }
    var at: u8[>..] = [];
    lits(at, f.at);
    at.append(", ");
    lit(at, 1.0e-3 * mesh_size);
    at.append(", ");
    quoted(at, label);
    if f.what == F_LET {
        let c, linear = let_coefficients(f.let_index, f.at);
        format(out, "        let tn, te = resultants_at(", at, ");\n        let v = ");
        lit(out, c[0]);
        for k in 18 {
            if c[k + 1] == 0.0 { continue; }
            out.append(" + ");
            lit(out, c[k + 1]);
            format(out, if k < 9 { " * tn[" } else { " * te[" }, k % 9, "]");
        }
        out.append(";\n");
        return;
    }
    let d0 = if f.what == F_W { 0 } else { 3 };
    for j in 3 {
        if comp >= 0 && j != comp { continue; }
        format(out, "        let v", AXIS_NAMES[j], " = dof_at(", at, ", ", DOF_NAMES[d0 + j], ");\n");
    }
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n");
    emit_value(out, r.field, r.label, r.field.comp);
    var v: u8[>..] = [];
    let unit = field_unit(r.field);
    if r.field.what == F_LET {
        format(v, "quantity(v, \"", unit, "\")");
    } else if r.field.comp >= 0 {
        format(v, "quantity(v", AXIS_NAMES[r.field.comp], ", \"", unit, "\")");
    } else {
        format(v, "vector(vx, vy, vz, \"", unit, "\")");
    }
    out.append("        report(");
    quoted(out, step);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", v, ");\n    }\n");
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    format(out, "// step ", s.name, ": static\n");
    format(out, "fn step_", s.name, "() {\n");
    out.append("    shells::clear_supports();\n    shells::clear_loads();\n");
    var own: i64[>..] = [];
    for i in global_bcs { own.push(i); }
    for i in s.first_own..s.first_own + s.own_count { own.push(i); }
    for i in own { if bcs[i].kind == BC_FIX { emit_bc(out, bcs[i]); } }
    if q_normal || q[0] != 0.0 || q[1] != 0.0 || q[2] != 0.0 {
        format(out, "    // physics ", physics_name, ": q\n    for c in shells::shells.len { shells::shell_load(c, ");
        lits(out, q);
        format(out, ", ", if q_normal { "true" } else { "false" }, "); }\n");
    }
    for i in own { if bcs[i].kind == BC_EDGE { emit_bc(out, bcs[i]); } }
    out.append("    solve(");
    quoted(out, s.name);
    format(out, ");\n    equilibrium(\"", codegen::system.force, "\");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s.name, r); }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(", ");
        quoted(out, w_name);
        out.append(", ");
        quoted(out, theta_name);
        out.append(");\n");
    }
    out.append("}\n\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n");
    emit_value(out, x.field, x.name, x.field.comp);
    let v = if x.field.what == F_LET { str("v") } else { str("v", AXIS_NAMES[x.field.comp]) };
    out.append("        expect(");
    quoted(out, x.name);
    format(out, ", ", v, ", ");
    lit(out, x.want);
    out.append(", ");
    lit(out, x.within);
    format(out, ", ", if x.relative { "true" } else { "false" }, ", \"", field_unit(x.field), "\", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}

// The number of cells along `length` for a cell size of at most `size`.
fn cells_for(length: f64, size: f64) -> i64 { return max(ceil(length / size - 1.0e-9) as! i64, 1); }

// Writes the Goose program of the model that `read` read. The program
// imports the shells from `shells_path` and writes its VTU files under
// `results_dir`.
fn generate(shells_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    let sys = codegen::system;
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", sys.length, ", force in ", sys.force,
           ", angle in rad.\n");
    out.append("import std;\nimport os;\nimport math;\nimport shells from ");
    quoted(out, shells_path);
    out.append(";\n");
    out.append(PRELUDE);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    let nx = cells_for(x1 - x0, mesh_size);
    let nt = cells_for(radius * (t1 - t0), mesh_size);
    out.append("fn main() {\n    shells::clear();\n");
    format(out, "    // material ", material_name, ": ", angles.len, " plies from the top face\n");
    out.append("    shells::laminate_clear();\n");
    for a in angles {
        out.append("    shells::laminate_ply(");
        lit(out, a);
        out.append(", ");
        lit(out, t_ply);
        out.append(", ");
        lit(out, e1);
        out.append(", ");
        lit(out, e2);
        out.append(", ");
        lit(out, g12);
        out.append(", ");
        lit(out, nu12);
        out.append(", ");
        lit(out, g12);
        out.append(", ");
        lit(out, g12);
        out.append(");\n");
    }
    out.append("    let section = shells::laminate_section();\n");
    format(out, "    // mesh on ", shell_name, ": ", nx, " x ", nt, " MITC4 cells of at most ", mesh_size, " ", sys.length, "\n");
    out.append("    grid = shells::mesh_cylinder(");
    lit(out, radius);
    out.append(", ");
    lit(out, x0);
    out.append(", ");
    lit(out, x1);
    out.append(", ");
    lit(out, t0);
    out.append(", ");
    lit(out, t1);
    format(out, ", ", nx, ", ", nt, ", section, ");
    lits(out, axis);
    out.append(");\n");
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { emit_expect(out, x); }
    out.append("    print_records();\n");
    out.append("    if failed { exit(1); }\n");
    out.append("}\n");
}
