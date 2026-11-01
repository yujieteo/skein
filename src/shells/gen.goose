// The code generator of the shell models: a Skein syntax tree to a Goose
// program that calls src/shells/.
//
// It takes a model whose physics uses the template reissner_mindlin (spec
// §4.3), the template timoshenko, or both. The first Surface is the side of
// a cylinder on the x axis, `boundary(c, y^2 + z^2 == R^2 && <cuts>)`, with
// cuts that keep one arc of θ (θ turns from +z to +y). The physics of the
// shell is on it or on a later Surface: a part of it with & and - of
// solids, for example a cutout (gen_geo.goose). The mesh is a grid of x and
// θ lines with a line on each beam, finer near a `refine` curve, with the
// nodes moved onto the cuts and onto the boundaries of the materials. It
// has MITC4 cells (order 1) or MITC9 cells (order 2). A Curve of the shell
// is the edges of the mesh on it. A Point is the node at it. Each material
// of the shell is a laminate of the template arguments `layup` and
// `mirror` and of the material.
//
// The beams are on a Curve of straight lines and of rings
// `boundary(shell, x == c)`, with pattern, pattern_polar, move, |, and & or
// - of a solid. With a shell, each line is along x on the cylinder, and
// the tie `bc beams: w_b = w` and `theta_b = theta` makes the beams share
// the shell nodes. With no shell, each line is meshed alone, and a modal
// or a buckling step can follow. Each other construct gives an error
// `unsupported` at its token, with a fix.
//
// gen_expr.goose writes the reports, the expects and the outputs as code,
// and gen_emit.goose writes the program. A submodel takes the state of
// another shell model through a results value (spec §15.1).
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
import .gen_geo;
import .gen_expr;
import .gen_emit;

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

// The names of the `let` and `param` items and of the lets of the physics,
// with their value nodes and their types.
var value_names: u8[..64][>..];
var value_nodes: i64[>..];
var value_types: i64[>..];
var let_names: u8[..64][>..];
var let_nodes: i64[>..];
var let_types: i64[>..];
var depth = 0;

// x, y and z are the point `px` while `at_point`. N and eps_m are the
// tensors `tn` and `te` while `in_physics`. n is `edge_n`, the outward
// normal of an edge on the surface, while `edge_mode`.
var at_point = false;
var edge_mode = false;
var edge_n: f64[3] = [0.0, 0.0, 0.0];
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
    let b = bound(word);
    if b >= 0 { return bind_vals[b]; }
    var node = -1;
    if in_physics && find(let_names, word) >= 0 { node = let_nodes[find(let_names, word)]; }
    if node < 0 && in_physics && cur_mat >= 0 { node = prop_node(cur_mat, word); }
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
    if at_point && word == "n" && edge_mode { return vector(edge_n); }
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
    if user_fn(n) >= 0 {
        let f = user_fn(n);
        let body = bind_call(n, f);
        if body < 0 { return scalar(0.0); }
        let r = value(body);
        unbind(f);
        return r;
    }
    if k == syntax::N_CALL && kind(child(n, 0)) == syntax::N_NAME {
        let f = name(child(n, 0));
        let count = syntax::child_count(n) - 1;
        if count == 1 && (f == "sin" || f == "cos" || f == "sqrt" || f == "abs") {
            let x = value(arg(n, 0)).v[0];
            if f == "sin" { return scalar(sin(x)); }
            if f == "cos" { return scalar(cos(x)); }
            if f == "sqrt" { return scalar(sqrt(x)); }
            return scalar(abs(x));
        }
        if count == 1 && f == "copy" { return scalar(copy_index(arg(n, 0))); }
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

// The copies of the patterns that copy(r) has met: the pattern node and
// the shape of its copy 0.
var copy_nodes: i64[>..] = [];
var copy_shapes: i64[>..] = [];

// copy(r) at the point px (spec §8.1): the index of the copy of the
// pattern r, a Solid, whose φ is the lowest.
fn copy_index(r: i64) -> f64 {
    var v = r;
    if kind(v) == syntax::N_NAME {
        let si = find(solid_names, name(v));
        if si >= 0 { v = solid_nodes[si]; }
    }
    let polar = is_call(v, "pattern_polar");
    if !at_point || !(polar || is_call(v, "pattern")) || syntax::child_count(v) != 4 {
        unsupported(tok(r), "this slice has copy(r) of a Solid of a pattern, in the value of a bc", "write let windows: Solid = pattern(window_0, (f_pitch, 0 m, 0 m), 8)");
        return 0.0;
    }
    var k = position(copy_nodes[..], v);
    if k < 0 {
        k = copy_nodes.len;
        copy_nodes.push(v);
        copy_shapes.push(shape_of(arg(v, 0)));
    }
    let n = round(scalar_value(arg(v, 2))) as! i64;
    let d = vec3(value(arg(v, 1)));
    let l = sqrt(dot3(d, d));
    var best = 0;
    var low = 1.0e300;
    for c in n {
        let q = if polar { rodrigues(px, [d[0] / l, d[1] / l, d[2] / l], -2.0 * PI * c / n) }
                else { [px[0] - c * d[0], px[1] - c * d[1], px[2] - c * d[2]] };
        let f = phi(copy_shapes[k], q);
        if f < low {
            low = f;
            best = c;
        }
    }
    return 1.0 * best;
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

// A region: a curve of the shell (its edges in the mesh, gen_geo.goose), the
// node at a point, or the pieces of a beam Curve. boundary(s, x == c) is
// also a ring piece, for a frame.
struct Region {
    name: u8[..64],
    node: i64,      // the value of a curve of the shell, or -1
    point: bool,
    at: f64[3],
    first_piece: i64,
    piece_count: i64,
}

// A piece of a beam Curve: the straight line from p0 to p1, or the arc of
// the ring of the shell at x = p0[0] from θ = ta to tb.
struct Piece {
    ring: bool,
    p0: f64[3],
    p1: f64[3],
    ta: f64,
    tb: f64,
}

// The bc kinds: fixed dofs, a force per length `N*n` on the edges of a
// curve, or a force or a moment `F_b*n`, `M_b*n` at a point.
let BC_FIX = 0;
let BC_EDGE = 1;
let BC_FORCE = 2;

struct Bc {
    region: i64,
    kind: i64,
    dof: i64,       // the first dof, 0 to 5
    count: i64,     // the dofs from dof: 1 or 3
    value: i64,     // the value: -1 for 0, a field of the results value, or the load
}

let STEP_STATIC = 0;
let STEP_MODAL = 1;
let STEP_BUCKLING = 2;

struct Step {
    name: u8[..64],
    kind: i64,
    count: i64,         // the modes of a modal or a buckling step
    first_own: i64,     // the bcs of the step: a range of `bcs`
    own_count: i64,
    first_report: i64,
    report_count: i64,
    first_output: i64,  // the values of its output: a range of `outputs`
    output_count: i64,
}

struct Report {
    node: i64,
    label: u8[..256],
}

struct Expect {
    name: u8[..64],
    node: i64,
    want: i64,
    within: i64,
    relative: bool,
    label: u8[..256],
    step: i64,      // the step that checks it at its end, or -1 for an item
}

var solid_names: u8[..64][>..];
var solid_nodes: i64[>..];
var solid_shapes: i64[>..];     // the shape of each Solid, or -1 before its first use
var shape_nodes: i64[>..];      // the nodes that node_shape made a shape for
var shape_ids: i64[>..];        // the shape of each of shape_nodes
var surface_names: u8[..64][>..];
var surface_shapes: i64[>..];
var table_names: u8[..64][>..];
var table_nodes: i64[>..];
var material_names: u8[..64][>..];
var material_nodes: i64[>..];
var regions: Region[>..];
var pieces: Piece[>..];
var bcs: Bc[>..];
var global_bcs: i64[>..];
var reports: Report[>..];
var outputs: i64[>..];
var steps: Step[>..];
var expects: Expect[>..];

// The cylinder of the shell: radius r about the x axis, from x0 to x1 and
// from t0 to t1 in θ. The meshed surface is mesh_g, a part of it.
var shell_name: u8[..64];
var shell_tok = -1;
var radius = 0.0;
var x0 = 0.0;
var x1 = 0.0;
var t0 = 0.0;
var t1 = 0.0;
var mesh_g = -1;
var mesh_size = 0.0;
var mesh_tok = -1;
var physics_tok = -1;
var physics_name: u8[..64];
var w_name: u8[..64];
var theta_name: u8[..64];
var args: i64[5] = [-1; 5];     // the nodes of layup, mirror, axis, top and q, or -1

let TEMPLATE_ARGS = ["layup", "mirror", "axis", "top", "q"];

// The materials of the shell and of the beams: an index of material_nodes,
// and the shape (shell) or the region (beams) where each one is, or -1 for
// the rest.
var shell_mats: i64[>..];
var shell_mat_g: i64[>..];
var beam_mats: i64[>..];
var beam_mat_r: i64[>..];
var cur_mat = -1;       // the material whose properties `named` gives

// The beam physics: its region, unknowns and the nodes of its template
// arguments, and whether the tie of each unknown is there.
var beam_tok = -1;
var beam_name: u8[..64];
var beam_region = -1;
var wb_name: u8[..64];
var thb_name: u8[..64];
var beam_args: i64[6] = [-1; 6];
var tie_w = false;
var tie_theta = false;

let BEAM_ARGS = ["A", "I_y", "I_z", "J", "y_axis", "offset"];

// The results value of a submodel (spec §15.1): its name, the model file,
// the model name in the file, and the step.
var source_name: u8[..64];
var source_path: u8[..256];
var source_model: u8[..64];
var source_step: u8[..64];

fn reset() {
    codegen::reset();
    value_names.clear();
    value_nodes.clear();
    value_types.clear();
    let_names.clear();
    let_nodes.clear();
    let_types.clear();
    cur_mat = -1;
    copy_nodes.clear();
    copy_shapes.clear();
    depth = 0;
    at_point = false;
    in_physics = false;
    edge_mode = false;
    solid_names.clear();
    solid_nodes.clear();
    solid_shapes.clear();
    shape_nodes.clear();
    shape_ids.clear();
    surface_names.clear();
    surface_shapes.clear();
    shapes.clear();
    patch_shape = -1;
    fn_names.clear();
    fn_nodes.clear();
    bind_names.clear();
    bind_vals.clear();
    refine_nodes.clear();
    refine_sizes.clear();
    refine_x.clear();
    refine_t.clear();
    refine_of.clear();
    order = 1;
    xs.clear();
    ts.clear();
    edges.clear();
    cell_mat.clear();
    table_names.clear();
    table_nodes.clear();
    material_names.clear();
    material_nodes.clear();
    regions.clear();
    pieces.clear();
    bcs.clear();
    global_bcs.clear();
    reports.clear();
    outputs.clear();
    steps.clear();
    expects.clear();
    fields.clear();
    field_ranks.clear();
    shell_name = [];
    shell_tok = -1;
    mesh_g = -1;
    mesh_size = 0.0;
    mesh_tok = -1;
    physics_tok = -1;
    physics_name = [];
    w_name = [];
    theta_name = [];
    args = [-1; 5];
    radius = 0.0;
    shell_mats.clear();
    shell_mat_g.clear();
    beam_mats.clear();
    beam_mat_r.clear();
    lams.clear();
    angles.clear();
    bmats.clear();
    beam_tok = -1;
    beam_name = [];
    beam_region = -1;
    wb_name = [];
    thb_name = [];
    beam_args = [-1; 6];
    tie_w = false;
    tie_theta = false;
    source_name = [];
    source_path = [];
    source_model = [];
    source_step = [];
}

// The value node of property `word` of material m (an index of
// material_nodes), or -1.
fn prop_node(m: i64, word: u8[:]) -> i64 {
    var p = first(material_nodes[m]);
    while p >= 0 {
        if name(p) == word { return child(p, 1); }
        p = next(p);
    }
    return -1;
}

// The declared type of property `word` of a material of the shell, or -1.
fn prop_type(word: u8[:]) -> i64 {
    for m in shell_mats {
        var p = first(material_nodes[m]);
        while p >= 0 {
            if name(p) == word { return child(p, 0); }
            p = next(p);
        }
    }
    return -1;
}

// Property `word` of shell material k.
fn prop_value(k: i64, word: u8[:], at: i64) -> Val {
    let p = prop_node(shell_mats[k], word);
    if p < 0 {
        unsupported(at, str("the material ", material_names[shell_mats[k]], " has no ", word), "give each material of the shell the property");
        return scalar(0.0);
    }
    let keep = cur_mat;
    cur_mat = shell_mats[k];
    in_physics = true;
    let v = value(p);
    in_physics = false;
    cur_mat = keep;
    return v;
}

fn is_source(n: i64) -> bool { return source_name.len > 0 && kind(n) == syntax::N_NAME && name(n) == source_name; }

// Whether n is `global.f`, with global the results value.
fn source_member(n: i64) -> bool { return kind(n) == syntax::N_MEMBER && is_source(first(n)); }

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

// `let end0: Curve = boundary(roof, x == 0 mm)`, a curve of the shell
// (gen_geo.goose), `let tip: Point = point((x, y, z) mm)`, or the lines of
// beams.
fn region(n: i64, v: i64, curve: bool) {
    var r = Region { name: [], node: -1, point: !curve, at: [0.0, 0.0, 0.0], first_piece: pieces.len, piece_count: 0 };
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
    if shell_curve(v) {
        r.node = v;
        // boundary(s, x == c) is also the ring at c, for a frame.
        let cond = if is_call(v, "boundary") { arg(v, 1) } else { -1 };
        if cond >= 0 && kind(cond) == syntax::N_BINARY && op(cond) == "==" && kind(child(cond, 0)) == syntax::N_NAME &&
           name(child(cond, 0)) == "x" {
            let c = scalar_value(child(cond, 1));
            pieces.push(Piece { ring: true, p0: [c, 0.0, 0.0], p1: [c, 0.0, 0.0], ta: t0, tb: t1 });
            r.piece_count = 1;
        }
        regions.push(r);
        return;
    }
    if !curve_pieces(v, pieces) { return; }
    r.piece_count = pieces.len - r.first_piece;
    regions.push(r);
}

let CURVE_FIX = "write the beam line as line(p0, p1), a pattern or pattern_polar of it, & or - a box, or boundary(shell, x == c) for a ring";

// The pieces of the beam Curve `v`, appended to `out`. False after an
// error. Each form changes the pieces that its operand appended, in place.
recursive fn curve_pieces(v: i64, out: Piece[>..]&) -> bool {
    let from = out.len;
    if is_call(v, "line") && syntax::child_count(v) == 3 {
        out.push(Piece { ring: false, p0: vec3(value(arg(v, 0))), p1: vec3(value(arg(v, 1))), ta: 0.0, tb: 0.0 });
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
        let g = shape_of(child(v, 1));
        if g < 0 { return false; }
        let count = out.len - from;
        for i in count { cut_piece(out[from + i], g, op(v) == "&", out); }
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
        // A turn by a about +x takes θ to θ - a.
        var q = p;
        let s = if d[0] > 0.0 { a } else { -a };
        q.ta -= s;
        q.tb -= s;
        return q;
    }
    let k: f64[3] = [d[0] / l, d[1] / l, d[2] / l];
    return Piece { ring: false, p0: rodrigues(p.p0, k, a), p1: rodrigues(p.p1, k, a), ta: 0.0, tb: 0.0 };
}

fn rodrigues(v: f64[3], k: f64[3], a: f64) -> f64[3] {
    let kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    let c: f64[3] = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]];
    var r: f64[3] = [0.0, 0.0, 0.0];
    for j in 3 { r[j] = v[j] * cos(a) + c[j] * sin(a) + k[j] * kv * (1.0 - cos(a)); }
    return r;
}

// The parts of piece p inside solid g (keep_in), or out of it, into out.
fn cut_piece(p: Piece, g: i64, keep_in: bool, out: Piece[>..]&) {
    let n = 400;
    var start = -1.0;
    var prev = 0.0;
    var prev_keep = false;
    for k in n + 1 {
        let s = 1.0 * k / n;
        let keep = (phi(g, piece_point(p, s)) < 0.0) == keep_in;
        if k == 0 {
            if keep { start = 0.0; }
        } else if keep != prev_keep {
            let c = piece_cross(p, g, prev, s);
            if keep { start = c; } else { push_part(p, start, c, out); }
        }
        prev = s;
        prev_keep = keep;
    }
    if prev_keep { push_part(p, start, 1.0, out); }
}

fn piece_point(p: Piece, s: f64) -> f64[3] {
    if p.ring { return shell_point(p.p0[0], p.ta + s * (p.tb - p.ta)); }
    return on_line(p, s);
}

// The place between a and b where φ of g changes sign along piece p.
fn piece_cross(p: Piece, g: i64, a: f64, b: f64) -> f64 {
    let inside = phi(g, piece_point(p, a)) < 0.0;
    var lo = a;
    var hi = b;
    for k in 100 {
        let m = 0.5 * (lo + hi);
        if (phi(g, piece_point(p, m)) < 0.0) == inside { lo = m; } else { hi = m; }
    }
    return 0.5 * (lo + hi);
}

fn push_part(p: Piece, a: f64, b: f64, out: Piece[>..]&) {
    if b - a <= 1.0e-9 { return; }
    if p.ring {
        out.push(Piece { ring: true, p0: p.p0, p1: p.p1, ta: p.ta + a * (p.tb - p.ta), tb: p.ta + b * (p.tb - p.ta) });
    } else {
        out.push(Piece { ring: false, p0: on_line(p, a), p1: on_line(p, b), ta: 0.0, tb: 0.0 });
    }
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
        solid_shapes.push(-1);
    } else if t == "Surface" && shell_tok < 0 {
        shell(n, v);
        if ok() { patch_shape = add_shape(S_PATCH, -1, -1, [0.0; 8]); }
    } else if t == "Surface" {
        let g = shape_of(v);
        if g >= 0 {
            surface_names.push(w);
            surface_shapes.push(g);
        }
    } else if t == "Curve" || t == "Point" {
        region(n, v, t == "Curve");
    } else if t == "Results" {
        results_item(n, v);
    } else if t == "Volume" || t == "Region" {
        unsupported(tok(n), str("a ", t, " is not in this slice"), "use a Surface for the shell");
    } else {
        value_names.push(w);
        value_nodes.push(v);
        value_types.push(ty);
    }
}

let RESULTS_FIX = "write let global: Results = results(\"model.skn\", step: \"load\"), with a shell model in the same unit system";

// `let global: Results = results(file, step: s)`. The file is relative to
// the model file, and its model is a shell model in the same unit system.
fn results_item(n: i64, v: i64) {
    guard source_name.len == 0 && is_call(v, "results") && syntax::child_count(v) == 3 && kind(arg(v, 0)) == syntax::N_STRING &&
          kind(arg(v, 1)) == syntax::N_STRING else {
        unsupported(tok(n), "this slice has one results value of a file and a step", RESULTS_FIX);
        return;
    };
    let file = arg(v, 0);
    let model_path = syntax::files[syntax::tokens[tok(file)].file].path;
    var cut = 0;
    for c, i in model_path { if c == '/' { cut = i + 1; } }
    let quoted_file = name(file);
    codegen::set256(source_path, str(model_path[..cut], quoted_file[1..quoted_file.len - 1]));
    let s = name(arg(v, 1));
    set64(source_step, s[1..s.len - 1]);
    var text: u8[>..] = [];
    guard read_file(source_path, text) else {
        unsupported(tok(file), str("cannot read '", source_path, "'"), "give the path of a model, relative to this file");
        return;
    };
    let f = syntax::add_file(source_path, text);
    syntax::lex(f);
    let prog = syntax::parse(f);
    var c = first(prog);
    var system = "";
    while c >= 0 {
        if kind(c) == syntax::N_HEADER {
            set64(source_model, name(c));
            system = syntax::text(tok(child(c, 1)));
        }
        c = next(c);
    }
    var same = false;
    for sys, i in codegen::SYSTEM_NAMES { if sys == system && codegen::SYSTEMS[i].l == codegen::system.l && codegen::SYSTEMS[i].m == codegen::system.m && codegen::SYSTEMS[i].t == codegen::system.t { same = true; } }
    if source_model.len == 0 || !same {
        unsupported(tok(v), str("the model of '", source_path, "' has no header in the unit system of this model"), RESULTS_FIX);
        return;
    }
    set64(source_name, name(n));
}

// Whether r is the meshed surface, a beam Curve, or a | of these.
recursive fn mesh_region(r: i64) -> bool {
    if kind(r) == syntax::N_BINARY && op(r) == "|" { return mesh_region(child(r, 0)) && mesh_region(child(r, 1)); }
    if kind(r) != syntax::N_NAME { return false; }
    if shell_tok >= 0 && (name(r) == shell_name || find(surface_names, name(r)) >= 0) { return true; }
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
        } else if kind(s) == syntax::N_REFINE && find_region(v) >= 0 && regions[find_region(v)].node >= 0 {
            refine_nodes.push(regions[find_region(v)].node);
            refine_sizes.push(scalar_value(child(s, 1)));
        } else {
            unsupported(tok(s), "the mesh setting is not in this slice",
                        "use only size, cell: quad, and refine c: 5 mm with c a Curve of the shell");
        }
        s = next(s);
    }
    if mesh_size <= 0.0 { unsupported(tok(n), "the mesh needs a size > 0", "add size: 50 mm"); }
    for h in refine_sizes { if h <= 0.0 || h > mesh_size { unsupported(tok(n), "a refine size is not in 0 < size <= the mesh size", "give a smaller size"); } }
}

// The shape of Surface name n (the shell, or a later Surface), or -1.
fn surface_shape(n: i64) -> i64 {
    if kind(n) != syntax::N_NAME || shell_tok < 0 { return -1; }
    if name(n) == shell_name { return patch_shape; }
    let i = find(surface_names, name(n));
    return if i >= 0 { surface_shapes[i] } else { -1 };
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
            mesh_g = surface_shape(first(c));
            if mesh_g < 0 { unsupported(tok(first(c)), "the physics is not on a Surface of the cylinder", PHYSICS_FIX); }
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
            var m = first(c);
            while m >= 0 {
                let i = find(material_names, name(m));
                let g = if first(m) >= 0 { surface_shape(first(m)) } else { -1 };
                if i < 0 || (first(m) >= 0 && g < 0) {
                    unsupported(tok(m), "the material is not one that the model declares, on a Surface of the shell", PHYSICS_FIX);
                }
                // The first match wins, so the material of the rest comes last.
                if first(m) < 0 && next(m) >= 0 {
                    unsupported(tok(m), "a material with no on takes the rest, so it comes last", "write material a on pad, b");
                }
                shell_mats.push(i);
                shell_mat_g.push(g);
                m = next(m);
            }
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
    if shell_mats.len == 0 { unsupported(tok(n), "the shell has no material", PHYSICS_FIX); }
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
                if at >= 0 && at < 6 {
                    beam_args[at] = first(a);
                } else {
                    unsupported(tok(a), "timoshenko takes A, I_y, I_z, J, y_axis and offset only", BEAM_FIX);
                }
                i++;
                a = next(a);
            }
        } else if k == syntax::N_UNKNOWN {
            if wb_name.len == 0 { set64(wb_name, name(c)); } else if thb_name.len == 0 { set64(thb_name, name(c)); }
        } else if k == syntax::N_MATASSIGN {
            var m = first(c);
            while m >= 0 {
                let i = find(material_names, name(m));
                let r = if first(m) >= 0 { find_region(first(m)) } else { -1 };
                if i < 0 || (first(m) >= 0 && (r < 0 || regions[r].piece_count == 0)) {
                    unsupported(tok(m), "the material is not one that the model declares, on a beam Curve", BEAM_FIX);
                }
                if first(m) < 0 && next(m) >= 0 {
                    unsupported(tok(m), "a material with no on takes the rest, so it comes last", "write material s on stringers, f");
                }
                beam_mats.push(i);
                beam_mat_r.push(r);
                m = next(m);
            }
        } else {
            unsupported(tok(c), "the beam physics statement is not in this slice", "use only material");
        }
        c = next(c);
    }
    for i in 5 {
        if beam_args[i] < 0 { unsupported(tok(n), str("timoshenko needs ", BEAM_ARGS[i]), BEAM_FIX); }
    }
    if beam_mats.len == 0 { unsupported(tok(n), "the beams have no material", BEAM_FIX); }
}

var method_order = -1;

fn method(n: i64) {
    if (name(n) != physics_name && name(n) != beam_name) || syntax::text(syntax::nodes[n].aux) != "cg" {
        unsupported(tok(n), "this slice has only the cg method of the physics", "write method shell: cg { order: 1 }");
        return;
    }
    var s = first(child(n, 0));
    while s >= 0 {
        if kind(s) == syntax::N_SETTING && name(s) == "order" {
            let o = scalar_value(first(s));
            if (o != 1.0 && o != 2.0) || (method_order > 0 && o != method_order) {
                unsupported(tok(s), "the shell and the beams have one order, 1 (MITC4) or 2 (MITC9)", "write order: 1 or order: 2 in each method");
            }
            method_order = o as! i64;
            order = method_order;
        } else if !(kind(s) == syntax::N_SETTING && name(s) == "mesh") {
            unsupported(tok(s), "the method setting is not in this slice", "use only mesh and order");
        }
        s = next(s);
    }
}

let BC_FIX_TEXT = "write bc end: w = 0, bc end: theta.x = 0, bc cut: w = global.w, bc end: N*n = (1, 0, 0) N/mm, or bc tip: F_b*n = (0, 0, -1) kN";

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

var in_step = false;

fn add_bc(b: Bc) {
    bcs.push(b);
    if !in_step { global_bcs.push(bcs.len - 1); }
}

// A `bc` item or statement.
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
    let reg = regions[r];
    let d, count = target_dofs(target);
    if d >= 0 {
        if !reg.point && (reg.node < 0 || physics_tok < 0) {
            unsupported(tok(n), "this slice fixes a Curve of the shell, or a Point", "write bc root: w = 0 with let root: Point = point(...)");
            return;
        }
        if source_member(v) {
            if physics_tok < 0 || reg.point { unsupported(tok(v), "a value of the results value needs a Curve of the shell", BC_FIX_TEXT); }
            if !(if d < 3 { is_w(name(v)) } else { is_theta(name(v)) }) {
                unsupported(tok(v), "the results value has the type of the unknown: w = global.w, theta = global.theta", BC_FIX_TEXT);
            }
            add_bc(Bc { region: r, kind: BC_FIX, dof: d, count: count, value: v });
            return;
        }
        let x = value(v);
        for j in 9 {
            if x.v[j] != 0.0 { unsupported(tok(v), "this slice fixes a value of 0, or a value of a results value", BC_FIX_TEXT); }
        }
        if x.n != 1 && (count == 1 || x.n != 3) { unsupported(tok(v), "the value has the wrong size", BC_FIX_TEXT); }
        add_bc(Bc { region: r, kind: BC_FIX, dof: d, count: count, value: -1 });
        return;
    }
    let lf = if kind(target) == syntax::N_BINARY && op(target) == "*" && kind(child(target, 1)) == syntax::N_NAME &&
                name(child(target, 1)) == "n" && kind(child(target, 0)) == syntax::N_NAME { name(child(target, 0)) } else { "" };
    if beam_tok >= 0 && (lf == "F_b" || lf == "M_b") && reg.point {
        if value(v).n != 3 { unsupported(tok(v), "the force or the moment at a point is a vector", BC_FIX_TEXT); }
        add_bc(Bc { region: r, kind: BC_FORCE, dof: if lf == "F_b" { 0 } else { 3 }, count: 3, value: v });
        return;
    }
    if source_of(target) == "N*n" && reg.node >= 0 && physics_tok >= 0 {
        add_bc(Bc { region: r, kind: BC_EDGE, dof: 0, count: 3, value: v });
        return;
    }
    unsupported(tok(target), str("the bc target '", source_of(target), "' is not in this slice"), BC_FIX_TEXT);
}

fn step(n: i64) {
    let proc = child(n, 0);
    var s = Step { name: [], kind: STEP_STATIC, count: 0, first_own: bcs.len, own_count: 0, first_report: reports.len, report_count: 0,
                   first_output: outputs.len, output_count: 0 };
    set64(s.name, name(n));
    let p = name(proc);
    if p == "modal" || p == "buckling" {
        s.kind = if p == "modal" { STEP_MODAL } else { STEP_BUCKLING };
        let a = first(proc);
        if a < 0 || syntax::nodes[a].aux < 0 || syntax::text(syntax::nodes[a].aux) != "count" || next(a) >= 0 {
            unsupported(tok(proc), str("a ", p, " step of the beams takes count only"), str("write step s: ", p, "(count: 4) { ... }"));
        } else {
            s.count = round(scalar_value(first(a))) as! i64;
        }
        if physics_tok >= 0 || beam_tok < 0 {
            unsupported(tok(proc), str("this slice has ", p, " steps of the beams only, with no shell"), "remove the shell, or use a static step");
        }
        if s.count < 1 || s.count > 50 { unsupported(tok(proc), "the count is not 1 to 50", "write count: 4"); }
        // The geometric stiffness comes from the axial force of the solve before.
        if s.kind == STEP_BUCKLING && (steps.len == 0 || steps[steps.len - 1].kind != STEP_STATIC) {
            unsupported(tok(proc), "a buckling step needs a static step just before it", "add step preload: static { ... } before it");
        }
    } else if p != "static" || first(proc) >= 0 {
        unsupported(tok(proc), "this slice has static steps with no arguments, and modal and buckling steps", "write step load: static { ... }");
    }
    in_step = true;
    var c = next(proc);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_BC {
            bc(c);
            if s.kind != STEP_STATIC { unsupported(tok(c), "a modal or a buckling step has no bc of its own", "put the bc in the static step before"); }
        } else if k == syntax::N_OUTPUT {
            var o = first(c);
            while o >= 0 {
                if kind(o) != syntax::N_EVERY { outputs.push(o); }
                o = next(o);
            }
        } else if k == syntax::N_REPORT {
            var r = first(c);
            while r >= 0 {
                var rep = Report { node: first(r), label: [] };
                codegen::set256(rep.label, if syntax::nodes[r].aux >= 0 { syntax::text(syntax::nodes[r].aux) } else { source_of(first(r)) });
                reports.push(rep);
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
    s.output_count = outputs.len - s.first_output;
    steps.push(s);
}

// An `expect` item (step -1), or an `expect` of step `step`.
fn expect(n: i64, step: i64) {
    var x = Expect { name: [], node: child(n, 0), want: child(n, 1), within: child(n, 2),
                     relative: kind(child(n, 2)) == syntax::N_PERCENT, label: [], step: step };
    set64(x.name, name(n));
    codegen::set256(x.label, str(source_of(child(n, 1)), " within ", source_of(child(n, 2))));
    expects.push(x);
}

// Reads the program `prog`. False on an error, and `syntax::errors` holds
// the errors.
fn read(prog: i64) -> bool {
    reset();
    method_order = -1;
    var c = first(prog);
    while c >= 0 {
        let k = kind(c);
        if k == syntax::N_HEADER {
            codegen::header(c);
        } else if k == syntax::N_LET || k == syntax::N_PARAM {
            value_item(c);
            if !ok() { return false; }
        } else if k == syntax::N_FN {
            var w: u8[..64] = [];
            set64(w, name(c));
            fn_names.push(w);
            fn_nodes.push(c);
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
                        "remove it: this slice runs linear static shells and beams, and modes and buckling of beams");
        }
        if !ok() { return false; }
        c = next(c);
    }
    let at = if first(prog) >= 0 { tok(first(prog)) } else { tok(prog) };
    if shell_tok < 0 && beam_tok < 0 { unsupported(at, "the model has no shell surface", SHELL_FIX); }
    if mesh_tok < 0 { unsupported(at, "the model has no mesh", "add mesh m on roof { size: 1000 mm, cell: quad }"); }
    if physics_tok < 0 && beam_tok < 0 { unsupported(at, "the model has no physics", PHYSICS_FIX); }
    if steps.len == 0 { unsupported(at, "the model has no step", "add step load: static { ... }"); }
    if refine_nodes.len > 0 && physics_tok < 0 { unsupported(mesh_tok, "refine needs the shell", "remove refine"); }
    if ok() && beam_tok >= 0 { check_pieces(); }
    if ok() && physics_tok >= 0 { build_grid(); }
    if ok() && physics_tok >= 0 { read_sections(); }
    if ok() && physics_tok >= 0 { read_loads(); }
    if ok() && beam_tok >= 0 { read_beams(); }
    if ok() { check_values(); }
    return ok();
}

// The lines, the cells and the edges of the mesh of the shell, with a line
// on each beam.
fn build_grid() {
    var ring_x: f64[>..] = [];
    var line_t: f64[>..] = [];
    if beam_region >= 0 {
        for p in beam_pieces() {
            if p.ring { ring_x.push(p.p0[0]); } else { line_t.push(theta_of(p.p0)); }
        }
    }
    make_lines(mesh_g, ring_x, line_t);
    if lat_nx() * lat_nt() > 4000000 {
        unsupported(mesh_tok, str("the mesh has ", (lat_nx() + 1) * (lat_nt() + 1), " nodes, more than this slice solves"), "make the size larger");
        return;
    }
    make_cells(mesh_g, shell_mat_g);
    var kept_cells = 0;
    for m in cell_mat {
        if m == -2 {
            unsupported(physics_tok, "a part of the shell has no material", "give the rest of the shell a material with no on");
            return;
        }
        if m >= 0 { kept_cells++; }
    }
    if kept_cells == 0 { unsupported(physics_tok, "the meshed surface has no cell", "check the surface of the physics"); }
    make_edges();
}

// Each report, expect and output value can run: a dry emission.
fn check_values() {
    var scratch: u8[>..] = [];
    fields.clear();
    field_ranks.clear();
    if beam_tok >= 0 { make_runs(); }
    for s in steps {
        eigen_step = s.kind != STEP_STATIC;
        for r in reports[s.first_report..s.first_report + s.report_count] {
            let rank, d = emit_global_expr(r.node, scratch);
            if rank == 9 { unsupported(tok(r.node), "a report of a tensor is not in this slice", "report a component"); }
        }
    }
    for x in expects {
        // A top level expect reads the results of the last step.
        eigen_step = steps[if x.step >= 0 { x.step } else { steps.len - 1 }].kind != STEP_STATIC;
        for v in [x.node, x.want, x.within] {
            let rank, d = emit_global_expr(v, scratch);
            if rank != 1 { unsupported(tok(v), "an expect compares scalars", "write a component, such as (reaction(end)).x"); }
        }
    }
    eigen_step = false;
    for o in outputs {
        if is_output_word(o) { continue; }
        if physics_tok < 0 { unsupported(tok(o), "an output of a value needs the shell", "write output w"); }
        let k = field_index(o);
    }
    // The bcs, so that a value of a bc that cannot run stops here.
    for b in bcs { emit_bc(scratch, b); }
    for k in fields.len { emit_field_fn(scratch, k); }
    fields.clear();
    field_ranks.clear();
}

fn is_output_word(o: i64) -> bool {
    return kind(o) == syntax::N_NAME && (is_w(name(o)) || is_theta(name(o)) || (name(o) == "Nx" && beam_tok >= 0));
}

// --- the section and the loads ----------------------------------------------------

// The laminate of a shell material: the ply properties, and its plies at
// angles[first..first + count] from the top face.
struct Lam {
    e1: f64,
    e2: f64,
    g12: f64,
    nu12: f64,
    t_ply: f64,
    first: i64,
    count: i64,
}

var lams: Lam[>..];
var angles: f64[>..];
var axis: f64[3] = [1.0, 0.0, 0.0];
// The load per area: a vector, or a pressure along the top normal.
var q: f64[3] = [0.0, 0.0, 0.0];
var q_normal = false;

let MATERIAL_FIX = "give the material E1, E2, G12, nu12, t_ply and h";
let LAM_PROPS = ["E1", "E2", "G12", "nu12", "t_ply", "h"];

fn read_sections() {
    for k in shell_mats.len { read_section(k); }
}

// The plies of `layup` (an angle, or a table column of angles) and the
// properties of shell material k: the plies repeat to fill h (h/2 with
// mirror), and mirror adds them again in the reverse sequence (spec §4.3).
fn read_section(k: i64) {
    let m = material_nodes[shell_mats[k]];
    var v: f64[6] = [0.0; 6];
    for i in 6 { v[i] = prop_value(k, LAM_PROPS[i], tok(m)).v[0]; }
    if !ok() { return; }
    var lam = Lam { e1: v[0], e2: v[1], g12: v[2], nu12: v[3], t_ply: v[4], first: angles.len, count: 0 };
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
    let count = round(fill / max(lam.t_ply, 1.0e-300)) as! i64;
    if !(lam.e1 > 0.0 && lam.e2 > 0.0 && lam.g12 > 0.0 && lam.t_ply > 0.0) || count < 1 || abs(count * lam.t_ply - fill) > 1.0e-6 * fill {
        unsupported(tok(m), "h is not a whole number of plies of t_ply (2 t_ply with mirror)", MATERIAL_FIX);
        return;
    }
    if base.len > 0 {
        for j in count { angles.push(base[j % base.len]); }
    }
    if mirror {
        let c = angles.len - lam.first;
        for j in c { angles.push(angles[lam.first + c - 1 - j]); }
    }
    lam.count = angles.len - lam.first;
    lams.push(lam);
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

// --- the beams -------------------------------------------------------------------

// A beam material: E, nu, A, I_y, I_z and J, the density (0 if none), and
// y_axis and the offset: a vector, or, with the *_radial flag, [0] times
// the outward normal of the shell.
struct BeamMat {
    section: f64[6],
    rho: f64,
    y: f64[3],
    y_radial: bool,
    off: f64[3],
    off_radial: bool,
}

var bmats: BeamMat[>..];

let BEAM_PROPS = ["E", "nu"];

// The section and the vectors of each beam material, from the material
// and the template arguments. The argument values can use the properties
// of the material.
fn read_beams() {
    let samples = beam_samples();
    for k in beam_mats.len {
        cur_mat = beam_mats[k];
        in_physics = true;
        var b = BeamMat { section: [0.0; 6], rho: 0.0, y: [0.0, 1.0, 0.0], y_radial: false, off: [0.0, 0.0, 0.0], off_radial: false };
        for i in 2 {
            let p = prop_node(cur_mat, BEAM_PROPS[i]);
            if p < 0 { unsupported(tok(material_nodes[cur_mat]), str("the beam material has no ", BEAM_PROPS[i]), "give the beam material E and nu"); }
            if p >= 0 { b.section[i] = scalar_value(p); }
        }
        for i in 4 { b.section[2 + i] = scalar_value(beam_args[i]); }
        let rho = prop_node(cur_mat, "rho");
        if rho >= 0 { b.rho = scalar_value(rho); }
        b.y, b.y_radial = beam_vector(beam_args[4], samples, "y_axis");
        if beam_args[5] >= 0 { b.off, b.off_radial = beam_vector(beam_args[5], samples, "offset"); }
        in_physics = false;
        cur_mat = -1;
        bmats.push(b);
    }
    for s in steps {
        if s.kind != STEP_MODAL { continue; }
        for b in bmats { if b.rho <= 0.0 { unsupported(beam_tok, "a modal step needs the density rho of each beam material", "add rho: Density = 7850 kg/m^3"); } }
    }
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
        out.push(piece_point(p, 0.0));
        out.push(piece_point(p, 0.5));
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
    if beam_region < 0 { return; }
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

// Whether piece p is on a piece of region r: its middle is on one.
fn piece_in(r: i64, p: Piece) -> bool {
    if r < 0 { return true; }
    let reg = regions[r];
    let m = piece_point(p, 0.5);
    let tol = 1.0e-6 * max(radius, 1.0e-300) + 1.0e-9 * mesh_size;
    for q in pieces[reg.first_piece..reg.first_piece + reg.piece_count] {
        if q.ring != p.ring { continue; }
        if q.ring {
            let t = theta_of(m);
            if abs(q.p0[0] - p.p0[0]) <= tol && t >= min(q.ta, q.tb) - 1.0e-9 && t <= max(q.ta, q.tb) + 1.0e-9 { return true; }
            continue;
        }
        let d = [q.p1[0] - q.p0[0], q.p1[1] - q.p0[1], q.p1[2] - q.p0[2]];
        let l2 = sq(d[0]) + sq(d[1]) + sq(d[2]);
        let s = ((m[0] - q.p0[0]) * d[0] + (m[1] - q.p0[1]) * d[1] + (m[2] - q.p0[2]) * d[2]) / l2;
        let c = on_line(q, s);
        if s >= -1.0e-9 && s <= 1.0 + 1.0e-9 && sqrt(sq(c[0] - m[0]) + sq(c[1] - m[1]) + sq(c[2] - m[2])) <= tol { return true; }
    }
    return false;
}

// The beam material of piece p: the first whose region holds it.
fn piece_mat(p: Piece) -> i64 {
    for r, k in beam_mat_r { if piece_in(r, p) { return k; } }
    return -1;
}

fn line_cells(p: Piece) -> i64 { return order * cells_for(sqrt(sq(p.p1[0] - p.p0[0]) + sq(p.p1[1] - p.p0[1]) + sq(p.p1[2] - p.p0[2])), mesh_size); }

// The number of cells along `length` for a cell size of at most `size`.
fn cells_for(length: f64, size: f64) -> i64 { return max(ceil(length / size - 1.0e-9) as! i64, 1); }

// The dofs of the mesh of the model that `read` read: 6 at each node.
fn dof_estimate() -> i64 {
    if physics_tok >= 0 { return 6 * (lat_nx() + 1) * (lat_nt() + 1); }
    var nodes = 0;
    for p in beam_pieces() { nodes += line_cells(p) + 1; }
    return 6 * nodes;
}
