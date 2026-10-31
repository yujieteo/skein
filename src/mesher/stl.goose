// `skein export model.skn --geometry name --format stl` (docs/command_line.md
// §7.1, spec §22.5): `check`, then the mesh of one `let` shape, and its
// surface as ASCII STL on standard output. A Solid gives the boundary
// faces of its tetrahedra, with the normals out of the solid. A 2D
// Surface gives its triangles, in the plane z = 0. The coordinates are in
// the length unit of the header, as for `skein verify`.
import std;
import os;
import math;
import dictionary;
namespace mesher;
import syntax from "../syntax";
import check from "../check";
import .geometry;
import .mesh;
import .quality;
import .tri;
import .tet;
import .verify;

let EXPORT_USAGE = "usage: skein export <file.skn> --geometry <name> [--format stl] [--size <length>]\n";

// A length option such as "2 mm" or "0.002m", in the length unit of the
// header. 0 when it is not a number > 0 with a length unit.
fn length_option(s: u8[:]) -> f64 {
    var i = 0;
    while i < s.len && (syntax::is_digit(s[i]) || s[i] == '.' || s[i] == 'e' || s[i] == '+' || s[i] == '-') { i++; }
    let x, good = parse_flt(s[..i]);
    let u = trim(s[i..]);
    guard good && x > 0.0 else { return 0.0; }
    if u == "m" { return x * metre; }
    let prefixes = "numck";
    let factors: f64[5] = [1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2, 1.0e3];
    guard u.len == 2 && u[1] == 'm' && position(prefixes, u[0]) >= 0 else { return 0.0; }
    return x * factors[position(prefixes, u[0])] * metre;
}

// The `size` of mesh block n, or 0 when it has none.
fn block_size(n: i64) -> f64 {
    var s = syntax::nodes[syntax::child(n, 1)].first;
    while s >= 0 {
        if syntax::nodes[s].kind == syntax::N_SETTING && node_text(s) == "size" {
            return eval_value(syntax::nodes[s].first).v[0];
        }
        s = syntax::nodes[s].next;
    }
    return 0.0;
}

// The default size: the size of the first mesh block on `name`, else of
// the first mesh block. 0 when the model has no mesh block.
fn model_size(prog: i64, name: u8[:]) -> f64 {
    var first = -1;
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        if syntax::nodes[item].kind == syntax::N_MESH {
            if first < 0 { first = item; }
            if node_text(syntax::child(item, 0)) == name { return block_size(item); }
        }
        item = syntax::nodes[item].next;
    }
    return if first >= 0 { block_size(first) } else { 0.0 };
}

// Writes one facet a, b, c, with the unit normal of the right-hand rule.
fn facet(out: u8[>..]&, a: f64[3], b: f64[3], c: f64[3]) {
    let n = cross3(sub3(b, a), sub3(c, a));
    let l = norm3(n);
    let s = if l > 0.0 { 1.0 / l } else { 0.0 };
    format(out, "  facet normal ", n[0] * s, " ", n[1] * s, " ", n[2] * s, "\n    outer loop\n");
    for p in [a, b, c] { format(out, "      vertex ", p[0], " ", p[1], " ", p[2], "\n"); }
    out.append("    endloop\n  endfacet\n");
}

// The faces of the tetrahedra that are in one cell only, each with its
// nodes in the order that turns its normal out of the cell.
fn boundary_facets(out: u8[>..]&) {
    let n = point_count();
    fn key(c: i64, skip: i64) -> i64 {
        var f: i64[..3] = [];
        for i in 4 {
            if i != skip { f.push(cell_node(c, i)); }
        }
        sort(f);
        return (f[0] * n + f[1]) * n + f[2];
    }
    var faces = dictionary<i64, i64> {};
    for c in cell_count() {
        for skip in 4 { faces.update(key(c, skip), 0) { it += 1; }; }
    }
    for c in cell_count() {
        for skip in 4 {
            let count = faces.get(key(c, skip));
            guard count != null && count == 1 else { continue; };
            var f: f64[3][..3] = [];
            for i in 4 {
                if i != skip { f.push(vertex(cell_node(c, i))); }
            }
            let inside = sub3(vertex(cell_node(c, skip)), f[0]);
            if dot3(cross3(sub3(f[1], f[0]), sub3(f[2], f[0])), inside) > 0.0 {
                facet(out, f[0], f[2], f[1]);
            } else {
                facet(out, f[0], f[1], f[2]);
            }
        }
    }
}

// `skein export`, with the STL in `out`. The result is the exit status
// (command_line.md §10). A fault of the model goes to standard error.
fn export_stl(args: u8[:][:], out: u8[>..]&) -> i64 {
    // The index in args of the file, of the geometry and of the size.
    var at = -1;
    var geo = -1;
    var sized = -1;
    var i = 2;
    while i < args.len {
        let a = args[i];
        if (a == "--geometry" || a == "--format" || a == "--size") && i + 1 < args.len {
            let v = args[i + 1];
            if a == "--geometry" { geo = i + 1; }
            if a == "--size" { sized = i + 1; }
            guard a != "--format" || v == "stl" else {
                write_stderr(str("skein export: the format '", v, "' is not in this version: write --format stl\n"));
                return 2;
            }
            i += 2;
            continue;
        }
        guard a.len > 0 && a[0] != '-' && at < 0 else {
            write_stderr(EXPORT_USAGE);
            return 2;
        }
        at = i;
        i++;
    }
    guard at >= 0 && geo >= 0 else {
        write_stderr(str("skein export: give one model file and --geometry\n", EXPORT_USAGE));
        return 2;
    }
    let path = args[at];
    let name = args[geo];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein export: cannot read '", path, "'\n"));
        return 2;
    }
    check::check_text(path, text);
    guard check::error_count() == 0 else { return model_fault(false); }
    let dim = read_geometry(check::model_prog);
    let g = let_node(name);
    guard g >= 0 else {
        write_stderr(str("skein export: '", name, "' is not a let of the model\n"));
        return 2;
    }
    let size = if sized >= 0 { length_option(args[sized]) } else { model_size(check::model_prog, name) };
    guard size > 0.0 else {
        let msg: u8[>..] = if sized >= 0 { str("skein export: the size '", args[sized], "' is not a length > 0, such as 2 mm\n") }
                           else { "skein export: the model has no mesh block: give --size, such as --size 2 mm\n" };
        write_stderr(msg);
        return 2;
    }
    let v = eval_value(g);
    if syntax::errors.len == 0 && v.n != 0 {
        reader_fault(syntax::nodes[g].tok, str("'", name, "' is not a shape"), "export a let with a Solid or a 2D Surface");
    }
    guard syntax::errors.len == 0 else { return model_fault(true); }
    let b = bounds(v.g);
    guard (b[5] > b[2]) == (dim == 3) else {
        reader_fault(syntax::nodes[g].tok, str("'", name, "' is not a ", dim, "d shape"),
                     if dim == 3 { "export a Solid" } else { "export a 2D Surface" });
        return model_fault(true);
    }
    let done = if dim == 3 { mesh_3d(v.g, size) } else { mesh_2d(v.g, size) };
    guard done else {
        reader_fault(syntax::nodes[g].tok, "the mesher did not finish in its budget of points", "give a larger --size");
        return model_fault(true);
    }
    format(out, "solid ", name, "\n");
    if dim == 3 {
        boundary_facets(out);
    } else {
        for c in cell_count() { facet(out, vertex(cell_node(c, 0)), vertex(cell_node(c, 1)), vertex(cell_node(c, 2))); }
    }
    format(out, "endsolid ", name, "\n");
    return 0;
}

// Writes the diagnostics of the model to standard error, after the faults
// of the reader when `take` is true. Standard output is for the STL only.
// The result is the exit status 1.
fn model_fault(take: bool) -> i64 {
    if take { check::take_syntax_errors(); }
    check::sort_diags();
    var err: u8[>..] = [];
    check::write_toon(err);
    write_stderr(err);
    return 1;
}

fn export_command(args: u8[:][:]) -> i64 {
    var out: u8[>..] = [];
    let status = export_stl(args, out);
    if status == 0 { write_stdout(out); }
    return status;
}
