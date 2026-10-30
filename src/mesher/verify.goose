// `skein verify model.skn` (docs/command_line.md §1): `check`, then the
// mesh of each `mesh` block and the pre-solve check mesh_quality. The
// other pre-solve checks are not in this version.
//
// The reader evaluates the geometry of the model: numbers with a length
// or an angle unit, vectors, + - * / and the names of `let` items, the
// shapes box, cylinder, rect, circle, polygon, extrude and revolve, and
// the operators |, & and -. Each length is in the length unit of the
// header (m for SI, mm for the mm systems), and each angle is in radians.
import std;
import os;
import math;
namespace mesher;
import syntax from "../syntax";
import fmt.fmt;
import check from "../check";
import .geometry;
import .mesh;
import .quality;
import .tri;
import .tet;

let VERIFY_USAGE = "usage: skein verify <file.skn>\n";

// A value: a geometry node (n = 0), a scalar (n = 1), or a vector of 2 or
// 3 components (n = 2 or 3).
struct Value {
    n: i64,
    v: f64[3],
    g: i64,
}

var metre = 1.0;            // one metre in the length unit of the header
var length_name: u8[..4] = [];
var let_names: u8[..64][>..] = [];
var let_values: i64[>..] = [];
var eval_depth = 0;
var polygon_scratch: f64[>..<] = [];

fn reader_fault(t: i64, msg: u8[:], fix: u8[:]) {
    if syntax::errors.len == 0 { syntax::error_at(t, "verify_unsupported", msg, fix); }
}

fn node_text(n: i64) -> const u8[:] { syntax::text(syntax::nodes[n].tok) }

fn number_value(x: f64) -> Value { Value { n: 1, v: [x, 0.0, 0.0], g: -1 } }

// The factor of a unit: a length to the header unit, an angle to radians.
fn unit_factor(u: i64) -> f64 {
    let n = syntax::nodes[u].first;
    if syntax::nodes[n].kind == syntax::N_UNIT_NAME {
        let w = node_text(n);
        if w == "rad" { return 1.0; }
        if w == "deg" { return PI / 180.0; }
        if w == "m" { return metre; }
        let prefixes = "numck";
        let factors: f64[5] = [1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2, 1.0e3];
        if w.len == 2 && w[1] == 'm' {
            let p = position(prefixes, w[0]);
            if p >= 0 { return factors[p] * metre; }
        }
    }
    reader_fault(syntax::nodes[u].tok, str("the unit '", node_text(n), "' is not a length or an angle"),
                 "write a length unit such as mm or m, or an angle unit, deg or rad");
    return 1.0;
}

fn let_node(word: u8[:]) -> i64 {
    for w, i in let_names {
        if w == word { return let_values[i]; }
    }
    return -1;
}

// The argument i of the call c (0 is the first argument).
fn argument(c: i64, i: i64) -> i64 { syntax::nodes[syntax::child(c, i + 1)].first }

recursive fn eval_value(n: i64) -> Value {
    let k = syntax::nodes[n].kind;
    let t = syntax::nodes[n].tok;
    if k == syntax::N_NUMBER {
        let x, good = parse_flt(syntax::text(t));
        if !good { reader_fault(t, "this number cannot be read", "write the number as digits"); }
        let u = syntax::nodes[n].first;
        return number_value(if u >= 0 { x * unit_factor(u) } else { x });
    }
    if k == syntax::N_NAME {
        if syntax::text(t) == "pi" { return number_value(PI); }
        let v = let_node(syntax::text(t));
        guard v >= 0 && eval_depth < 100 else {
            reader_fault(t, str("'", syntax::text(t), "' is not a let with a geometry or a number"),
                         "declare it before with let");
            return number_value(0.0);
        }
        eval_depth++;
        let r = eval_value(v);
        eval_depth--;
        return r;
    }
    if k == syntax::N_VECTOR {
        var r = Value { n: 0, v: [0.0, 0.0, 0.0], g: -1 };
        var f = 1.0;
        var c = syntax::nodes[n].first;
        while c >= 0 {
            if syntax::nodes[c].kind == syntax::N_UNIT {
                f = unit_factor(c);
            } else {
                let x = eval_value(c);
                if x.n != 1 || r.n >= 3 { reader_fault(syntax::nodes[c].tok, "a vector has 2 or 3 numbers", "write (x, y) or (x, y, z)"); }
                if r.n < 3 { r.v[r.n] = x.v[0]; }
                r.n++;
            }
            c = syntax::nodes[c].next;
        }
        for j in 3 { r.v[j] *= f; }
        return r;
    }
    if k == syntax::N_UNARY && syntax::text(t) == "-" {
        var r = eval_value(syntax::nodes[n].first);
        if r.n == 0 { reader_fault(t, "a shape has no negative", "remove the -"); }
        for j in 3 { r.v[j] = -r.v[j]; }
        return r;
    }
    if k == syntax::N_BINARY {
        let op = syntax::text(t);
        let a = eval_value(syntax::child(n, 0));
        let b = eval_value(syntax::child(n, 1));
        if a.n == 0 || b.n == 0 {
            guard a.n == 0 && b.n == 0 else {
                reader_fault(t, "a shape and a number", "combine two shapes with |, & or -");
                return a;
            }
            if op == "|" { return Value { n: 0, v: a.v, g: union(a.g, b.g) }; }
            if op == "&" { return Value { n: 0, v: a.v, g: intersect(a.g, b.g) }; }
            if op == "-" { return Value { n: 0, v: a.v, g: subtract(a.g, b.g) }; }
            reader_fault(t, str("'", op, "' of two shapes"), "combine two shapes with |, & or -");
            return a;
        }
        if op == "+" || op == "-" {
            if a.n != b.n { reader_fault(t, "a number and a vector are added", "add two numbers or two vectors"); }
            var r = a;
            for j in 3 { r.v[j] = if op == "+" { a.v[j] + b.v[j] } else { a.v[j] - b.v[j] }; }
            return r;
        }
        if op == "*" && (a.n == 1 || b.n == 1) {
            var r = if a.n == 1 { b } else { a };
            let s = if a.n == 1 { a.v[0] } else { b.v[0] };
            for j in 3 { r.v[j] *= s; }
            return r;
        }
        if op == "/" && b.n == 1 {
            var r = a;
            for j in 3 { r.v[j] /= b.v[0]; }
            return r;
        }
        if op == "^" && a.n == 1 && b.n == 1 { return number_value(pow(a.v[0], b.v[0])); }
    }
    if k == syntax::N_CALL { return shape_call(n); }
    reader_fault(t, "this expression is not a shape or a number that the mesher reads",
                 "write a number with a unit, a vector, a let name, + - * / ^, a shape, or |, & and - of shapes");
    return number_value(0.0);
}

recursive fn shape_call(n: i64) -> Value {
    let t = syntax::nodes[n].tok;
    let name = node_text(syntax::child(n, 0));
    let count = syntax::child_count(n) - 1;
    fn arg(i: i64, want: i64) -> Value {
        let v = eval_value(argument(n, i));
        if v.n != want {
            let what: u8[>..] = if want == 0 { "a shape" } else if want == 1 { "a number" } else { str("a vector of ", want) };
            reader_fault(syntax::nodes[argument(n, i)].tok, str("argument ", i + 1, " of ", name, " is not ", what),
                         "see the shapes of spec §8.1");
        }
        return v;
    }
    fn shape(g: i64) -> Value { Value { n: 0, v: [0.0, 0.0, 0.0], g: g } }
    if name == "box" && count == 2 {
        let a = arg(0, 3);
        let b = arg(1, 3);
        return shape(box(a.v[0], a.v[1], a.v[2], b.v[0], b.v[1], b.v[2]));
    }
    if name == "cylinder" && count == 4 {
        let c = arg(0, 3);
        let a = arg(1, 3);
        return shape(cylinder(c.v[0], c.v[1], c.v[2], a.v[0], a.v[1], a.v[2], arg(2, 1).v[0], arg(3, 1).v[0]));
    }
    if name == "rect" && count == 2 {
        let a = arg(0, 2);
        let b = arg(1, 2);
        return shape(rect(a.v[0], a.v[1], b.v[0], b.v[1]));
    }
    if name == "circle" && count == 2 {
        let c = arg(0, 2);
        return shape(circle(c.v[0], c.v[1], arg(1, 1).v[0]));
    }
    if name == "polygon" && count >= 3 {
        // The points go on a global stack: a local array cannot live
        // across the recursive call of arg.
        let mark = polygon_scratch.len;
        for i in count {
            let p = arg(i, 2);
            polygon_scratch.push(p.v[0]);
            polygon_scratch.push(p.v[1]);
        }
        let g = polygon(polygon_scratch[mark..]);
        polygon_scratch.resize(mark);
        return shape(g);
    }
    if name == "extrude" && count == 2 {
        let s = arg(0, 0);
        return shape(extrude(s.g, arg(1, 1).v[0]));
    }
    if name == "revolve" && count == 2 {
        let s = arg(0, 0);
        return shape(revolve(s.g, arg(1, 1).v[0]));
    }
    reader_fault(t, str("the shape '", name, "' with ", count, " arguments is not in the mesher"),
                 "use box, cylinder, rect, circle, polygon, extrude or revolve, with |, & and -");
    return number_value(0.0);
}

// The record of one mesh, for the table of the report.
struct MeshRecord {
    name: u8[..64],
    cell: u8[..4],
    size: f64,
    q: Quality,
}

// Meshes block n. False on an error in the block.
fn mesh_block(n: i64, dim: i64, records: MeshRecord[>..]&) -> bool {
    let region = syntax::child(n, 0);
    guard syntax::nodes[region].kind == syntax::N_NAME else {
        reader_fault(syntax::nodes[region].tok, "the region of the mesh is not a let name", "name a let shape");
        return false;
    }
    let g = eval_value(region);
    guard syntax::errors.len == 0 else { return false; }
    guard g.n == 0 else {
        reader_fault(syntax::nodes[region].tok, "the region of the mesh is not a shape", "name a let shape");
        return false;
    }
    var size = 0.0;
    // A copy, not a slice of a string: mesh_3d and mesh_2d reallocate the
    // byte arrays of the runtime mesh.
    var cell: u8[..3] = [];
    cell.append(if dim == 2 { "tri" } else { "tet" });
    var s = syntax::nodes[syntax::child(n, 1)].first;
    while s >= 0 {
        let value = syntax::nodes[s].first;
        if syntax::nodes[s].kind != syntax::N_SETTING {
            reader_fault(syntax::nodes[s].tok, "refine and layers are not in the mesher yet: it makes a uniform mesh",
                         "remove the line, and give a size that resolves the smallest feature");
        } else if node_text(s) == "size" {
            size = eval_value(value).v[0];
        } else if node_text(s) == "cell" {
            let c = node_text(value);
            if c != cell {
                reader_fault(syntax::nodes[value].tok, str("the mesher makes ", cell, " cells in a ", dim, "d model"),
                             str("write cell: ", cell));
            }
        } else {
            reader_fault(syntax::nodes[s].tok, str("the mesh setting '", node_text(s), "' is not in the mesher yet"),
                         "use size and cell");
        }
        s = syntax::nodes[s].next;
    }
    if size <= 0.0 { reader_fault(syntax::nodes[n].tok, "the mesh needs a size > 0", "add size: 5 mm"); }
    guard syntax::errors.len == 0 else { return false; }
    let b = bounds(g.g);
    let three = b[5] > b[2];
    guard three == (dim == 3) else {
        reader_fault(syntax::nodes[region].tok, str("the region is not a ", dim, "d shape"),
                     if dim == 3 { "mesh a Solid" } else { "mesh a 2D Surface" });
        return false;
    }
    let done = if dim == 3 { mesh_3d(g.g, size) } else { mesh_2d(g.g, size) };
    if !done { reader_fault(syntax::nodes[n].tok, "the mesher did not finish in its budget of points", "give a larger size"); }
    var r = MeshRecord { name: [], cell: [], size: size, q: mesh_quality() };
    r.name.append(node_text(n));
    r.cell.append(cell);
    records.push(r);
    return done;
}

fn verify(args: u8[:][:]) -> i64 {
    guard args.len == 3 && args[2].len > 0 && args[2][0] != '-' else {
        write_stderr(str("skein verify: give one model file\n", VERIFY_USAGE));
        return 2;
    }
    let status = check::command(args);
    if status != 0 { return status; }
    let path = args[2];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein verify: cannot read '", path, "'\n"));
        return 1;
    }
    let prog = fmt::parse_text(path, text);
    guard prog >= 0 else {
        fmt::write_errors(path);
        return 1;
    }
    clear_shapes();
    let_names.clear();
    let_values.clear();
    var dim = 3;
    metre = 1.0;
    length_name.clear();
    length_name.append("m");
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        let k = syntax::nodes[item].kind;
        if k == syntax::N_HEADER {
            let space = syntax::child_of_kind(item, syntax::N_SPACE);
            if space >= 0 && node_text(space) == "2d" { dim = 2; }
            let system = syntax::child_of_kind(item, syntax::N_UNITSYS);
            if system >= 0 && starts_with(node_text(system), "mm_") {
                metre = 1000.0;
                length_name.clear();
                length_name.append("mm");
            }
        } else if k == syntax::N_LET {
            var w: u8[..64] = [];
            w.append(node_text(item)[..min(node_text(item).len, 64)]);
            let_names.push(w);
            let_values.push(syntax::child(item, 1));
        }
        item = syntax::nodes[item].next;
    }
    var records: MeshRecord[>..] = [];
    item = syntax::nodes[prog].first;
    while item >= 0 {
        if syntax::nodes[item].kind == syntax::N_MESH && !mesh_block(item, dim, records) { break; }
        item = syntax::nodes[item].next;
    }
    if syntax::errors.len > 0 {
        fmt::write_errors(path);
        return 1;
    }
    let unit = length_name;
    var out: u8[>..] = [];
    format(out, "mesh[", records.len, "]{name,cell,size,cells,nodes,min_angle,max_angle,max_aspect,",
           if dim == 2 { "area" } else { "volume" }, "}:\n");
    for r in records {
        format(out, "  ", r.name, ",", r.cell, ",\"", r.size, " ", unit, "\",", r.q.cells, ",", r.q.nodes, ",\"");
        format_flt(out, r.q.min_angle, 2);
        out.append(" deg\",\"");
        format_flt(out, r.q.max_angle, 2);
        out.append(" deg\",");
        format_flt(out, r.q.max_aspect, 3);
        format(out, ",\"", r.q.measure, " ", unit, "^", dim, "\"\n");
    }
    var failed = false;
    format(out, "verify[", records.len, "]{stage,id,result,value,limit}:\n");
    for r in records {
        let pass = quality_passes(r.q);
        if !pass { failed = true; }
        format(out, "  pre,mesh_quality,", if pass { "pass" } else { "fail" }, ",\"");
        format_flt(out, r.q.min_angle, 2);
        if dim == 2 {
            out.append(" deg\",\"10 deg\"\n");
        } else {
            out.append(" deg, ");
            format_flt(out, r.q.max_angle, 2);
            out.append(" deg\",\"5 deg, 175 deg\"\n");
        }
    }
    write_stdout(out);
    return if failed { 1 } else { 0 };
}
