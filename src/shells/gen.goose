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
// arguments `layup` and `mirror` and of the material.
//
// A physics with the template timoshenko is on a Curve of straight lines
// and of rings `boundary(shell, x == c)`, with pattern, pattern_polar,
// move, |, and & or - of a box. With a shell, each line is along x on the
// cylinder, the grid has a node line on each beam, and the tie
// `bc beams: w_b = w` and `theta_b = theta` makes the beams share the
// shell nodes. With no shell, each line is meshed alone. The 2-node beam
// is of order 1. Each other construct gives an error `unsupported` at its
// token, with a fix.
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
var prop_from = 0;      // the first property of the material in use
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
    if node < 0 && in_physics && find(prop_names[prop_from..], word) >= 0 { node = prop_nodes[prop_from + find(prop_names[prop_from..], word)]; }
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
    if at_point && word == "n" && radius > 0.0 { return vector([0.0, px[1] / radius, px[2] / radius]); }
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

// A region: sides of the grid (bit 0 is x = x0, 1 is x = x1, 2 is θ = t0,
// 3 is θ = t1), the node at a point, or pieces of a beam Curve. A side
// x == c is also a ring piece.
struct Region {
    name: u8[..64],
    sides: i64,     // 0 for a point
    at: f64[3],
    first_piece: i64,
    piece_count: i64,
}

// A piece of a beam Curve: the straight line from p0 to p1, or the ring of
// the shell at x = p0[0].
struct Piece {
    ring: bool,
    p0: f64[3],
    p1: f64[3],
}

// The bc kinds: a fixed dof, a force per length `N*n` on a side, or a
// force or a moment `F_b*n`, `M_b*n` on a point (dofs `dof` to `dof` + 2).
let BC_FIX = 0;
let BC_EDGE = 1;
let BC_FORCE = 2;

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
let F_NX = 4;

struct Field {
    what: i64,
    comp: i64,      // 0 to 2, or -1 for the vector
    let_index: i64, // for F_LET; F_NX is the axial force of the beam nearest to `at`
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
    step: i64,      // the step that checks it at its end, or -1 for an item
}

var solid_names: u8[..64][>..];
var solid_nodes: i64[>..];
var table_names: u8[..64][>..];
var table_nodes: i64[>..];
var material_names: u8[..64][>..];
var material_nodes: i64[>..];
var regions: Region[>..];
var pieces: Piece[>..];
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

// The beam physics: its region, unknowns, material and the nodes of its
// template arguments, and whether the tie of each unknown is there.
var beam_tok = -1;
var beam_name: u8[..64];
var beam_region = -1;
var wb_name: u8[..64];
var thb_name: u8[..64];
var beam_material: u8[..64];
var beam_args: i64[7] = [-1; 7];
var tie_w = false;
var tie_theta = false;

let BEAM_ARGS = ["A", "I_y", "I_z", "J", "y_axis", "offset", "q"];

fn reset() {
    codegen::reset();
    value_names.clear();
    value_nodes.clear();
    let_names.clear();
    let_nodes.clear();
    let_types.clear();
    prop_names.clear();
    prop_nodes.clear();
    prop_from = 0;
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
    pieces.clear();
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
    radius = 0.0;
    beam_tok = -1;
    beam_name = [];
    beam_region = -1;
    wb_name = [];
    thb_name = [];
    beam_material = [];
    beam_args = [-1; 7];
    tie_w = false;
    tie_theta = false;
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
    var r = Region { name: [], sides: 0, at: [0.0, 0.0, 0.0], first_piece: pieces.len, piece_count: 0 };
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
        if !curve_pieces(v, pieces) { return; }
        r.piece_count = pieces.len - r.first_piece;
        regions.push(r);
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
    // x == c is also the ring at c, for a frame.
    if kind(child(cond, 0)) == syntax::N_NAME && name(child(cond, 0)) == "x" {
        let c = scalar_value(child(cond, 1));
        let tol = 1.0e-9 * (x1 - x0);
        if c >= x0 - tol && c <= x1 + tol {
            pieces.push(Piece { ring: true, p0: [c, 0.0, 0.0], p1: [c, 0.0, 0.0] });
            r.piece_count = 1;
        }
    }
    if r.sides != 0 || r.piece_count > 0 {
        regions.push(r);
        return;
    }
    if ok() { unsupported(tok(cond), str("the curve '", name(n), "' is not a side of the shell"), REGION_FIX); }
}

let CURVE_FIX = "write the beam line as line(p0, p1), a pattern or pattern_polar of it, & or - a box, or boundary(shell, x == c) for a ring";

// The pieces of the beam Curve `v`, appended to `out`. False after an
// error. Each form changes the pieces that its operand appended, in place.
recursive fn curve_pieces(v: i64, out: Piece[>..]&) -> bool {
    let from = out.len;
    if is_call(v, "line") && syntax::child_count(v) == 3 {
        out.push(Piece { ring: false, p0: vec3(value(arg(v, 0))), p1: vec3(value(arg(v, 1))) });
        return ok();
    }
    if kind(v) == syntax::N_NAME && find_region(v) >= 0 && regions[find_region(v)].piece_count > 0 {
        let r = regions[find_region(v)];
        // `out` may be `pieces`: copy by index.
        for i in r.piece_count {
            let p = pieces[r.first_piece + i];
            out.push(p);
        }
        return true;
    }
    if (is_call(v, "pattern") || is_call(v, "pattern_polar")) && syntax::child_count(v) == 4 {
        if !curve_pieces(arg(v, 0), out) { return false; }
        let count = out.len - from;
        let d = vec3(value(arg(v, 1)));
        let n = scalar_value(arg(v, 2));
        if !ok() { return false; }
        for k in 1..max(round(n) as! i64, 1) {
            for i in count {
                let p = out[from + i];
                let q = if is_call(v, "pattern") { moved(p, [k * d[0], k * d[1], k * d[2]], tok(v)) } else { turned(p, d, 2.0 * PI * k / n, tok(v)) };
                if !(q.ring && has_ring(out[from..], q.p0[0])) { out.push(q); }
            }
        }
        return ok();
    }
    if is_call(v, "move") && syntax::child_count(v) == 3 {
        if !curve_pieces(arg(v, 0), out) { return false; }
        let d = vec3(value(arg(v, 1)));
        for i in from..out.len { out[i] = moved(out[i], d, tok(v)); }
        return ok();
    }
    if kind(v) == syntax::N_BINARY && op(v) == "|" {
        return curve_pieces(child(v, 0), out) && curve_pieces(child(v, 1), out);
    }
    if kind(v) == syntax::N_BINARY && (op(v) == "&" || op(v) == "-") {
        if !curve_pieces(child(v, 0), out) { return false; }
        let lo, hi, good = box_of(child(v, 1));
        if !good { return false; }
        let count = out.len - from;
        for i in count {
            if out[from + i].ring {
                unsupported(tok(v), "this slice cuts beam lines only, not a ring", "cut the shell, or move the ring");
                return false;
            }
            clip(out[from + i], lo, hi, op(v) == "&", out);
        }
        // The cut pieces replace the pieces of the operand.
        let kept = out.len - from - count;
        for i in kept { out[from + i] = out[from + count + i]; }
        for i in count { out.pop(); }
        return true;
    }
    unsupported(tok(v), str("the curve '", source_of(v), "' is not in this slice"), CURVE_FIX);
    return false;
}

fn has_ring(out: Piece[:], x: f64) -> bool {
    for p in out { if p.ring && abs(p.p0[0] - x) <= 1.0e-12 * (abs(x) + 1.0) { return true; } }
    return false;
}

fn moved(p: Piece, d: f64[3], at: i64) -> Piece {
    if p.ring && (d[1] != 0.0 || d[2] != 0.0) {
        unsupported(at, "a ring moves along x only", "move the ring by (dx, 0 m, 0 m)");
    }
    var q = p;
    for k in 3 {
        q.p0[k] += d[k];
        q.p1[k] += d[k];
    }
    return q;
}

// p turned by angle a about the axis d through the origin (Rodrigues).
fn turned(p: Piece, d: f64[3], a: f64, at: i64) -> Piece {
    let l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if p.ring {
        if abs(d[1]) > 1.0e-12 * l || abs(d[2]) > 1.0e-12 * l { unsupported(at, "a ring turns about the x axis only", "turn the ring about (1, 0, 0)"); }
        return p;
    }
    let k: f64[3] = [d[0] / l, d[1] / l, d[2] / l];
    return Piece { ring: false, p0: rodrigues(p.p0, k, a), p1: rodrigues(p.p1, k, a) };
}

fn rodrigues(v: f64[3], k: f64[3], a: f64) -> f64[3] {
    let kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    let c: f64[3] = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]];
    var r: f64[3] = [0.0, 0.0, 0.0];
    for j in 3 { r[j] = v[j] * cos(a) + c[j] * sin(a) + k[j] * kv * (1.0 - cos(a)); }
    return r;
}

// The corners of `box(p0, p1)`, or of a Solid let that is one.
recursive fn box_of(v: i64) -> f64[3], f64[3], bool {
    if kind(v) == syntax::N_NAME && find(solid_names, name(v)) >= 0 && depth < 100 {
        depth++;
        let lo, hi, good = box_of(solid_nodes[find(solid_names, name(v))]);
        depth--;
        return lo, hi, good;
    }
    guard is_call(v, "box") && syntax::child_count(v) == 3 else {
        unsupported(tok(v), str("this slice cuts a beam line with a box only, not '", source_of(v), "'"), CURVE_FIX);
        return [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], false;
    }
    let a = vec3(value(arg(v, 0)));
    let b = vec3(value(arg(v, 1)));
    return [min(a[0], b[0]), min(a[1], b[1]), min(a[2], b[2])], [max(a[0], b[0]), max(a[1], b[1]), max(a[2], b[2])], ok();
}

// The part of line p in the closed box lo..hi (keep), or out of it, into out.
fn clip(p: Piece, lo: f64[3], hi: f64[3], keep: bool, out: Piece[>..]&) {
    var ta = 0.0;
    var tb = 1.0;
    var size = 0.0;
    for k in 3 { size = max(size, hi[k] - lo[k]); }
    let tol = 1.0e-9 * size;
    for k in 3 {
        let d = p.p1[k] - p.p0[k];
        if abs(d) <= 1.0e-300 {
            if p.p0[k] < lo[k] - tol || p.p0[k] > hi[k] + tol { tb = -1.0; }
            continue;
        }
        let a = (lo[k] - tol - p.p0[k]) / d;
        let b = (hi[k] + tol - p.p0[k]) / d;
        ta = max(ta, min(a, b));
        tb = min(tb, max(a, b));
    }
    let length = sqrt(sq(p.p1[0] - p.p0[0]) + sq(p.p1[1] - p.p0[1]) + sq(p.p1[2] - p.p0[2]));
    let small = 1.0e-9 * max(length, size) / max(length, 1.0e-300);
    var parts: f64[4] = [0.0, 0.0, 0.0, 0.0];
    var count = 0;
    if keep && tb - ta > small {
        parts = [ta, tb, 0.0, 0.0];
        count = 1;
    } else if !keep && tb - ta <= small {
        parts = [0.0, 1.0, 0.0, 0.0];
        count = 1;
    } else if !keep {
        if ta > small { parts[2 * count] = 0.0; parts[2 * count + 1] = ta; count++; }
        if tb < 1.0 - small { parts[2 * count] = tb; parts[2 * count + 1] = 1.0; count++; }
    }
    for c in count { out.push(Piece { ring: false, p0: on_line(p, parts[2 * c]), p1: on_line(p, parts[2 * c + 1]) }); }
}

fn sq(x: f64) -> f64 { return x * x; }

fn on_line(p: Piece, t: f64) -> f64[3] {
    return [p.p0[0] + t * (p.p1[0] - p.p0[0]), p.p0[1] + t * (p.p1[1] - p.p0[1]), p.p0[2] + t * (p.p1[2] - p.p0[2])];
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

// Whether r is the shell, a beam Curve, or a | of these.
recursive fn mesh_region(r: i64) -> bool {
    if kind(r) == syntax::N_BINARY && op(r) == "|" { return mesh_region(child(r, 0)) && mesh_region(child(r, 1)); }
    if kind(r) != syntax::N_NAME { return false; }
    if shell_tok >= 0 && name(r) == shell_name { return true; }
    return find_region(r) >= 0 && regions[find_region(r)].piece_count > 0;
}

fn mesh(n: i64) {
    if mesh_tok >= 0 {
        unsupported(tok(n), "a second mesh: this slice has one", "remove the second mesh");
        return;
    }
    mesh_tok = tok(n);
    let r = child(n, 0);
    if !mesh_region(r) {
        unsupported(tok(r), "the mesh is not on the shell and the beams", "write mesh m on skin | beams, with the names of the physics regions");
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
    let u = syntax::child_of_kind(n, syntax::N_USES);
    if u >= 0 && name(u) == "timoshenko" {
        beam_physics(n);
        return;
    }
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
            unsupported(tok(c), str("skein run has no generator for ", name(c), " on a shell"),
                        "use reissner_mindlin for the shell and timoshenko for the beams");
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

let BEAM_FIX = "write physics beam on bars: timoshenko(A: A, I_y: I_y, I_z: I_z, J: J, y_axis: (0, 1, 0)) solves w: Vector<Length>, theta: Vector<Angle> { material sec }";

fn beam_physics(n: i64) {
    if beam_tok >= 0 {
        unsupported(tok(n), "a second beam physics: this slice has one", "put all the beams in one physics");
        return;
    }
    beam_tok = tok(n);
    set64(beam_name, name(n));
    var c = first(n);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_ON {
            beam_region = find_region(first(c));
            if beam_region < 0 || regions[beam_region].piece_count == 0 {
                unsupported(tok(first(c)), "the beam is not on a Curve of lines or rings that a let names", CURVE_FIX);
            }
        } else if k == syntax::N_USES {
            var a = first(c);
            var i = 0;
            while a >= 0 {
                let at = if syntax::nodes[a].aux >= 0 { position(BEAM_ARGS, syntax::text(syntax::nodes[a].aux)) } else { i };
                if at >= 0 && at < 7 { beam_args[at] = first(a); }
                i++;
                a = next(a);
            }
        } else if k == syntax::N_UNKNOWN {
            if wb_name.len == 0 { set64(wb_name, name(c)); } else if thb_name.len == 0 { set64(thb_name, name(c)); }
        } else if k == syntax::N_MATASSIGN {
            let m = first(c);
            if syntax::child_count(c) != 1 || first(m) >= 0 || beam_material.len > 0 {
                unsupported(tok(c), "this slice gives one material to all the beams", BEAM_FIX);
            }
            set64(beam_material, name(m));
        } else {
            unsupported(tok(c), "the beam physics statement is not in this slice", "use only material");
        }
        c = next(c);
    }
    for i in 5 {
        if beam_args[i] < 0 { unsupported(tok(n), str("timoshenko needs ", BEAM_ARGS[i]), BEAM_FIX); }
    }
}

fn method(n: i64) {
    if (name(n) != physics_name && name(n) != beam_name) || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this slice has only the cg method of the physics", "write method shell: cg { order: 1 }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            if scalar_value(first(s)) != 1.0 {
                unsupported(tok(s), "the MITC4 shell and the 2-node beam are of order 1", "write order: 1");
            }
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this slice", "use only mesh and order");
        }
        s = next(s);
    }
}

let BC_FIX_TEXT = "write bc end: w = 0, bc end: theta.x = 0, bc end: N*n = (1, 0, 0) N/mm, or bc tip: F_b*n = (0, 0, -1) kN";

// The dof of `w`, `w.x`, `theta` or `theta.z`: the first dof and the
// count, or -1.
// The beam unknowns share the dofs of the shell unknowns (the tie).
fn target_dofs(t: i64) -> i64, i64 {
    if kind(t) == syntax::N_NAME && is_w(name(t)) { return 0, 3; }
    if kind(t) == syntax::N_NAME && is_theta(name(t)) { return 3, 3; }
    if kind(t) == syntax::N_MEMBER && kind(first(t)) == syntax::N_NAME && axis_of(name(t)) >= 0 {
        if is_w(name(first(t))) { return axis_of(name(t)), 1; }
        if is_theta(name(first(t))) { return 3 + axis_of(name(t)), 1; }
    }
    return -1, 0;
}

fn is_w(word: u8[:]) -> bool { return word.len > 0 && (word == w_name || word == wb_name); }

fn is_theta(word: u8[:]) -> bool { return word.len > 0 && (word == theta_name || word == thb_name); }

// A `bc` item or statement: one Bc for each fixed dof, or one load.
fn bc(n: i64) {
    let r = find_region(child(n, 0));
    let target = child(n, 1);
    let v = child(n, 2);
    if r < 0 {
        unsupported(tok(child(n, 0)), "the region of the bc is not a Curve or a Point that a let names", REGION_FIX);
        return;
    }
    // The tie: bc beams: w_b = w and bc beams: theta_b = theta.
    if kind(target) == syntax::N_NAME && kind(v) == syntax::N_NAME && (name(target) == wb_name || name(target) == thb_name) &&
       (name(v) == w_name || name(v) == theta_name) {
        let pair_w = name(target) == wb_name && name(v) == w_name;
        let pair_theta = name(target) == thb_name && name(v) == theta_name;
        if !(pair_w || pair_theta) || r != beam_region || in_step {
            unsupported(tok(n), "the tie is bc beams: w_b = w and bc beams: theta_b = theta, at the top level, on the region of the beams",
                        "tie each beam unknown to the shell unknown of the same type");
            return;
        }
        if pair_w { tie_w = true; } else { tie_theta = true; }
        return;
    }
    let d, count = target_dofs(target);
    if d >= 0 && regions[r].sides == 0 && regions[r].piece_count > 0 {
        unsupported(tok(n), "this slice fixes a beam at a Point only", "write bc root: w = 0 with let root: Point = point(...)");
        return;
    }
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
    let lf = if kind(target) == syntax::N_BINARY && op(target) == "*" && kind(child(target, 1)) == syntax::N_NAME &&
                name(child(target, 1)) == "n" && kind(child(target, 0)) == syntax::N_NAME { name(child(target, 0)) } else { "" };
    if beam_tok >= 0 && (lf == "F_b" || lf == "M_b") && regions[r].sides == 0 && regions[r].piece_count == 0 {
        let x = value(v);
        if x.n != 3 { unsupported(tok(v), "the force or the moment at a point is a vector", BC_FIX_TEXT); }
        bcs.push(Bc { region: r, kind: BC_FORCE, dof: if lf == "F_b" { 0 } else { 3 }, load: vec3(x) });
        if !in_step { global_bcs.push(bcs.len - 1); }
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

let REPORT_FIX = "write report w.z at (x, y, z) mm, report s at (x, y, z) mm with s a let of the physics, report Nx at (x, y, z) mm, or report reaction(end)";

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
    if kind(v) == syntax::N_NAME && name(v) == "Nx" && beam_tok >= 0 {
        f.what = F_NX;
        return f, true;
    }
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
        } else if k == syntax::N_EXPECT {
            expect(c, steps.len);
        } else {
            unsupported(tok(c), "the step statement is not in this slice", "use only bc, output, report and expect");
        }
        c = next(c);
    }
    in_step = false;
    s.own_count = bcs.len - s.first_own;
    s.report_count = reports.len - s.first_report;
    steps.push(s);
}

// An `expect` item (step -1), or an `expect` of step `step`.
fn expect(n: i64, step: i64) {
    let f, good = field(child(n, 0), true);
    if !good {
        if ok() { unsupported(tok(child(n, 0)), "an expect compares one component", "write w.z, or reaction(end).x"); }
        return;
    }
    var x = Expect { name: [], field: f, want: scalar_value(child(n, 1)), within: scalar_value(child(n, 2)),
                     relative: kind(child(n, 2)) == syntax::N_PERCENT, label: [], step: step };
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
            expect(c, -1);
        } else if k != syntax::N_SCALES && k != syntax::N_USE && k != syntax::N_STUDY {
            unsupported(tok(c), str("the item '", syntax::kind_name(k), "' is not in this slice"),
                        "remove it: this slice runs one linear static shell");
        }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if shell_tok < 0 && beam_tok < 0 { unsupported(at, "the model has no shell surface", SHELL_FIX); }
    if mesh_tok < 0 { unsupported(at, "the model has no mesh", "add mesh m on roof { size: 1000 mm, cell: quad }"); }
    if physics_tok < 0 && beam_tok < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { ... }"); }
    let m = find(material_names, material_name);
    if physics_tok >= 0 && m < 0 { unsupported(physics_tok, "the physics has no material that the model declares", PHYSICS_FIX); }
    if ok() && physics_tok >= 0 { read_section(material_nodes[m]); }
    if ok() && physics_tok >= 0 { read_loads(); }
    if ok() && beam_tok >= 0 { read_beam(); }
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

// --- the beams -------------------------------------------------------------------

// E, nu, A, I_y, I_z and J of the beam section.
var b_section: f64[6] = [0.0; 6];
// y_axis, the offset and the load per length: a vector, or, with the
// *_radial flag, [0] times the outward normal of the shell.
var b_y: f64[3] = [0.0, 1.0, 0.0];
var b_y_radial = false;
var b_off: f64[3] = [0.0, 0.0, 0.0];
var b_off_radial = false;
var b_q: f64[3] = [0.0, 0.0, 0.0];

// The section and the vectors of the beams, from the beam material and the
// template arguments. The argument names are the properties of the beam
// material.
fn read_beam() {
    let m = find(material_names, beam_material);
    if m < 0 || beam_region < 0 {
        unsupported(beam_tok, "the beam physics has no region or no material that the model declares", BEAM_FIX);
        return;
    }
    let keep = prop_names.len;
    prop_from = keep;
    var p = first(material_nodes[m]);
    while p >= 0 {
        var w: u8[..64] = [];
        set64(w, name(p));
        prop_names.push(w);
        prop_nodes.push(child(p, 1));
        p = next(p);
    }
    in_physics = true;
    for word, i in ["E", "nu"] {
        let k = find(prop_names[keep..], word);
        if k < 0 { unsupported(tok(material_nodes[m]), str("the beam material has no ", word), "give the beam material E and nu"); }
        if k >= 0 { b_section[i] = scalar_value(prop_nodes[keep + k]); }
    }
    for i in 4 { b_section[2 + i] = scalar_value(beam_args[i]); }
    let samples = beam_samples();
    b_y, b_y_radial = beam_vector(beam_args[4], samples, "y_axis");
    b_off = [0.0, 0.0, 0.0];
    b_off_radial = false;
    if beam_args[5] >= 0 { b_off, b_off_radial = beam_vector(beam_args[5], samples, "offset"); }
    b_q = [0.0, 0.0, 0.0];
    if beam_args[6] >= 0 {
        let q, radial = beam_vector(beam_args[6], samples, "q");
        if radial { unsupported(tok(beam_args[6]), "this slice has a uniform load per length on the beams", "write q: (0, 0, -1) kN/m"); }
        b_q = q;
    }
    in_physics = false;
    while prop_names.len > keep {
        prop_names.pop();
        prop_nodes.pop();
    }
    prop_from = 0;
    if ok() { check_pieces(); }
}

fn beam_pieces() -> const Piece[:] {
    let r = regions[beam_region];
    return pieces[r.first_piece..r.first_piece + r.piece_count];
}

// Points on up to 8 pieces of the beams: the start and the middle of each.
fn beam_samples() -> f64[3][>..] {
    let ps = beam_pieces();
    var out: f64[3][>..] = [];
    for k in min(ps.len, 8) {
        let p = ps[k * ps.len / min(ps.len, 8)];
        if p.ring {
            out.push(shell_point(p.p0[0], t0 + 0.2 * (t1 - t0)));
            out.push(shell_point(p.p0[0], t0 + 0.7 * (t1 - t0)));
        } else {
            out.push(p.p0);
            out.push(on_line(p, 0.5));
        }
    }
    return out;
}

// The vector of template argument n on the beams: the same at each sample,
// or s times the outward normal of the shell ([s, 0, 0], true).
fn beam_vector(n: i64, samples: f64[3][:], what: u8[:]) -> f64[3], bool {
    let a = vec3(value_at(n, samples[0]));
    var same = true;
    var radial = radius > 0.0;
    let r0 = sqrt(sq(samples[0][1]) + sq(samples[0][2]));
    let s = if r0 > 0.0 { (a[1] * samples[0][1] + a[2] * samples[0][2]) / r0 } else { 0.0 };
    let size = abs(a[0]) + abs(a[1]) + abs(a[2]);
    for p in samples {
        let b = vec3(value_at(n, p));
        let rp = sqrt(sq(p[1]) + sq(p[2]));
        for j in 3 {
            if abs(b[j] - a[j]) > 1.0e-12 * size { same = false; }
            let er = if j == 0 || rp == 0.0 { 0.0 } else { p[j] / rp };
            if abs(b[j] - s * er) > 1.0e-12 * abs(s) { radial = false; }
        }
    }
    if same || !ok() { return a, false; }
    if radial && s != 0.0 { return [s, 0.0, 0.0], true; }
    unsupported(tok(n), str("this slice has a beam ", what, " that is uniform or along the normal of the shell"),
                "write a constant vector, or s*e_r with e_r the outward normal");
    return a, false;
}

// With a shell, each line is along x on the cylinder; with no shell, each
// piece is a line.
fn check_pieces() {
    for p in beam_pieces() {
        if physics_tok < 0 && p.ring {
            unsupported(beam_tok, "a ring of beams needs the shell", "use lines for the beams");
            return;
        }
        if physics_tok >= 0 && !on_shell(p) {
            unsupported(beam_tok, str("a beam line from (", p.p0[0], ", ", p.p0[1], ", ", p.p0[2], ") is not along x on the shell"),
                        "put each line on the cylinder at one θ, or cut the lines with & box(...)");
            return;
        }
        if !p.ring && on_line(p, 0.0) == on_line(p, 1.0) {
            unsupported(beam_tok, "a beam line has no length", "give each line two different points");
            return;
        }
    }
    if physics_tok >= 0 && (!tie_w || !tie_theta) {
        unsupported(beam_tok, "the beams on the shell need the tie", str("add bc ", regions[beam_region].name, ": ", wb_name, " = ", w_name,
                    " and bc ", regions[beam_region].name, ": ", thb_name, " = ", theta_name));
    }
    if physics_tok >= 0 && ok() { grid_counts(); }
}

fn on_shell(p: Piece) -> bool {
    let tol = 1.0e-6 * radius;
    let xt = 1.0e-9 * (x1 - x0);
    let in_x = min(p.p0[0], p.p1[0]) >= x0 - xt && max(p.p0[0], p.p1[0]) <= x1 + xt;
    if p.ring { return in_x; }
    if abs(p.p0[1] - p.p1[1]) > tol || abs(p.p0[2] - p.p1[2]) > tol { return false; }
    let t = theta_of(p.p0);
    return in_x && abs(sqrt(sq(p.p0[1]) + sq(p.p0[2])) - radius) <= tol && t >= t0 - 1.0e-9 && t <= t1 + 1.0e-9;
}

// θ of a point of the shell, from t0 up.
fn theta_of(q: f64[3]) -> f64 {
    var d = atan2(q[1], q[2]) - t0;
    d -= 2.0 * PI * floor(d / (2.0 * PI));
    if d > t1 - t0 && 2.0 * PI - d < 1.0e-9 { d -= 2.0 * PI; }
    return t0 + d;
}

// The cells of the grid along x and along θ: at least those of the mesh
// size, and a node line on each beam.
fn grid_counts() -> i64, i64 {
    var nx = cells_for(x1 - x0, mesh_size);
    var nt = cells_for(radius * (t1 - t0), mesh_size);
    if beam_region < 0 { return nx, nt; }
    let limit_x = 1000 * nx;
    let limit_t = 1000 * nt;
    while nx < limit_x && !fits(nx, true) { nx++; }
    while nt < limit_t && !fits(nt, false) { nt++; }
    if !fits(nx, true) || !fits(nt, false) {
        unsupported(beam_tok, "no grid near the mesh size has a node line on each beam", "put the beams at a common step of x and of θ");
    }
    return nx, nt;
}

fn fits(count: i64, along_x: bool) -> bool {
    for p in beam_pieces() {
        if along_x {
            for e in [p.p0[0], p.p1[0]] { if !on_step(count, (e - x0) / (x1 - x0)) { return false; } }
        } else if !p.ring && !on_step(count, (theta_of(p.p0) - t0) / (t1 - t0)) {
            return false;
        }
    }
    return true;
}

fn on_step(count: i64, f: f64) -> bool { return abs(f * count - round(f * count)) <= 1.0e-6; }

fn grid_index(count: i64, f: f64) -> i64 { return round(f * count) as! i64; }

fn line_cells(p: Piece) -> i64 { return cells_for(sqrt(sq(p.p1[0] - p.p0[0]) + sq(p.p1[1] - p.p0[1]) + sq(p.p1[2] - p.p0[2])), mesh_size); }

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
    if f.what == F_REACTION || f.what == F_NX { return str(codegen::system.force); }
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

    // Dof d of the node at (x, y, z), on a beam with no shell.
    fn node_dof(x: f64, y: f64, z: f64, tol: f64, expr: u8[:], d: i64) -> f64 {
        return shells::displacement(node_at(x, y, z, tol, expr), d);
    }

    // The axial force of the beam with its middle nearest to (x, y, z).
    fn axial_at(x: f64, y: f64, z: f64) -> f64 {
        var best = -1;
        var dist = 0.0;
        for b, i in shells::beams {
            let p = shells::node_point(b.n[0]);
            let q = shells::node_point(b.n[1]);
            let d = shells::norm([0.5 * (p[0] + q[0]) - x, 0.5 * (p[1] + q[1]) - y, 0.5 * (p[2] + q[2]) - z]);
            if best < 0 || d < dist {
                best = i;
                dist = d;
            }
        }
        if best < 0 {
            failed = true;
            return 0.0;
        }
        return shells::beam_axial_force(best);
    }

    // The outward unit normal of the cylinder at p.
    fn radial(p: f64[3]) -> f64[3] {
        let r = sqrt(p[1] * p[1] + p[2] * p[2]);
        return [0.0, p[1] / r, p[2] / r];
    }

    // A beam from node a to node b. With `y_radial`, y is the normal at the
    // middle; with `o_radial`, the offset at each node is o[0] times the
    // normal there.
    fn beam(a: i64, b: i64, section: i64, y: f64[3], y_radial: bool, o: f64[3], o_radial: bool) {
        let pa = shells::node_point(a);
        let pb = shells::node_point(b);
        var ya = y;
        if y_radial { ya = radial([0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]); }
        var oa = o;
        var ob = o;
        if o_radial {
            let na = radial(pa);
            let nb = radial(pb);
            oa = [o[0] * na[0], o[0] * na[1], o[0] * na[2]];
            ob = [o[0] * nb[0], o[0] * nb[1], o[0] * nb[2]];
        }
        shells::add_beam(a, b, section, ya, oa, ob);
    }

    // Beams on the grid nodes from (i0, j0) to (i1, j1), a line of x or of θ.
    fn grid_beams(i0: i64, j0: i64, i1: i64, j1: i64, section: i64, y: f64[3], y_radial: bool, o: f64[3], o_radial: bool) {
        let along_x = j0 == j1;
        for k in (if along_x { i1 - i0 } else { j1 - j0 }) {
            let a = if along_x { shells::grid_node(grid, i0 + k, j0) } else { shells::grid_node(grid, i0, j0 + k) };
            let b = if along_x { shells::grid_node(grid, i0 + k + 1, j0) } else { shells::grid_node(grid, i0, j0 + k + 1) };
            beam(a, b, section, y, y_radial, o, o_radial);
        }
    }

    // The node at (x, y, z) if one is within tol, else a new node.
    fn node_or_new(p: f64[3], tol: f64) -> i64 {
        if shells::node_count() > 0 {
            let n, d = shells::nearest_node(p[0], p[1], p[2]);
            if d <= tol { return n; }
        }
        return shells::add_node(p[0], p[1], p[2]);
    }

    // `count` beams on the line from p to q. An end at a node shares it.
    fn line_beams(p: f64[3], q: f64[3], count: i64, tol: f64, section: i64, y: f64[3], o: f64[3]) {
        var a = node_or_new(p, tol);
        for k in count {
            let f = 1.0 * (k + 1) / count;
            let at: f64[3] = [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])];
            let b = if k == count - 1 { node_or_new(q, tol) } else { shells::add_node(at[0], at[1], at[2]) };
            shells::add_beam(a, b, section, y, o, o);
            a = b;
        }
    }

    // The load per length (qx, qy, qz) on each beam, half at each node.
    fn beam_load(qx: f64, qy: f64, qz: f64) {
        for b in shells::beams {
            let p = shells::node_point(b.n[0]);
            let q = shells::node_point(b.n[1]);
            let half = 0.5 * shells::norm([q[0] - p[0], q[1] - p[1], q[2] - p[2]]);
            for n in b.n {
                shells::force(n, shells::UX, half * qx);
                shells::force(n, shells::UY, half * qy);
                shells::force(n, shells::UZ, half * qz);
            }
        }
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
        format(solves, "  ", step, ",", shells::shells.len + shells::beams.len, ",", shells::dof_count(), ",0,",
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
        s.append(solvers::memory_records());
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
    } else if b.kind == BC_FORCE {
        out.append("    {\n");
        emit_node_at(out, "n", r.at, r.name);
        for j in 3 {
            if b.load[j] == 0.0 { continue; }
            format(out, "        shells::force(n, ", DOF_NAMES[b.dof + j], ", ");
            lit(out, b.load[j]);
            out.append(");\n");
        }
        out.append("    }\n");
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
    if f.what == F_NX {
        out.append("        let v = axial_at(");
        lits(out, f.at);
        out.append(");\n");
        return;
    }
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
        format(out, "        let v", AXIS_NAMES[j], " = ", if physics_tok >= 0 { "dof_at(" } else { "node_dof(" }, at, ", ",
               DOF_NAMES[d0 + j], ");\n");
    }
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n");
    emit_value(out, r.field, r.label, r.field.comp);
    var v: u8[>..] = [];
    let unit = field_unit(r.field);
    if r.field.what == F_LET || r.field.what == F_NX {
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
    if b_q[0] != 0.0 || b_q[1] != 0.0 || b_q[2] != 0.0 {
        format(out, "    // physics ", beam_name, ": q\n    beam_load(");
        lits(out, b_q);
        out.append(");\n");
    }
    for i in own { if bcs[i].kind != BC_FIX { emit_bc(out, bcs[i]); } }
    out.append("    solve(");
    quoted(out, s.name);
    format(out, ");\n    equilibrium(\"", codegen::system.force, "\");\n");
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s.name, r); }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output {
        var path: u8[>..] = [];
        format(path, results_dir, "/", s.name, ".vtu");
        out.append("    output(");
        quoted(out, path);
        out.append(", ");
        quoted(out, if physics_tok >= 0 { w_name[..] } else { wb_name[..] });
        out.append(", ");
        quoted(out, if physics_tok >= 0 { theta_name[..] } else { thb_name[..] });
        out.append(");\n");
    }
    out.append("}\n\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n");
    emit_value(out, x.field, x.name, x.field.comp);
    let v = if x.field.what == F_LET || x.field.what == F_NX { str("v") } else { str("v", AXIS_NAMES[x.field.comp]) };
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

// The dofs of the mesh of the model that `read` read: 6 at each node.
fn dof_estimate() -> i64 {
    if physics_tok >= 0 {
        let nx, nt = grid_counts();
        return 6 * (nx + 1) * (nt + 1);
    }
    var nodes = 0;
    for p in beam_pieces() { nodes += line_cells(p) + 1; }
    return 6 * nodes;
}

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
    out.append(";\nimport solvers from ");
    quoted(out, str(shells_path, "/../solvers"));
    out.append(";   // the shells import it\n");
    out.append(PRELUDE);
    for si in 0..steps.len { emit_step(out, si, results_dir); }
    out.append("fn main() {\n    shells::clear();\n");
    if physics_tok >= 0 { emit_shell(out); }
    if beam_tok >= 0 { emit_beams(out); }
    for s in steps { format(out, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(out, x); } }
    out.append("    print_records();\n");
    out.append("    if failed { exit(1); }\n");
    out.append("}\n");
}

// The laminate and the grid of the shell.
fn emit_shell(out: u8[>..]&) {
    let sys = codegen::system;
    let nx, nt = grid_counts();
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
}

// The beam section and the beams of each piece: on the grid lines of the
// shell, or on lines of their own.
fn emit_beams(out: u8[>..]&) {
    format(out, "    // material ", beam_material, ": the beam section\n    let beam_section = shells::beam_section(");
    for v, i in b_section {
        if i > 0 { out.append(", "); }
        lit(out, v);
    }
    out.append(");\n");
    var tail: u8[>..] = [];
    tail.append("beam_section, [");
    lits(tail, b_y);
    format(tail, "], ", if b_y_radial { "true" } else { "false" }, ", [");
    lits(tail, b_off);
    format(tail, "], ", if b_off_radial { "true" } else { "false" }, ");\n");
    let ps = beam_pieces();
    format(out, "    // physics ", beam_name, " on ", regions[beam_region].name, ": ", ps.len, " lines and rings\n");
    if physics_tok < 0 {
        for p in ps {
            out.append("    line_beams([");
            lits(out, p.p0);
            out.append("], [");
            lits(out, p.p1);
            format(out, "], ", line_cells(p), ", ");
            lit(out, 1.0e-3 * mesh_size);
            out.append(", beam_section, [");
            lits(out, b_y);
            out.append("], [");
            lits(out, b_off);
            out.append("]);\n");
        }
        return;
    }
    let nx, nt = grid_counts();
    for p in ps {
        let ia = grid_index(nx, (min(p.p0[0], p.p1[0]) - x0) / (x1 - x0));
        if p.ring {
            format(out, "    grid_beams(", ia, ", 0, ", ia, ", ", nt, ", ", tail);
        } else {
            let ib = grid_index(nx, (max(p.p0[0], p.p1[0]) - x0) / (x1 - x0));
            let j = grid_index(nt, (theta_of(p.p0) - t0) / (t1 - t0));
            format(out, "    grid_beams(", ia, ", ", j, ", ", ib, ", ", j, ", ", tail);
        }
    }
}
