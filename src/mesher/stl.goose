// `skein export model.skn --geometry name --format stl` (docs/command_line.md
// §7.1, spec §22): `check`, then the mesh of one `let` shape, and its
// surface as ASCII STL on standard output, or as STEP AP242 with `--format step`. A Solid gives the boundary
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

let EXPORT_USAGE = "usage: skein export <file.skn> --geometry <name> [--format stl|step] [--size <length>]\n";

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

// Adds one triangle a, b, c to `tris`, three points for each triangle.
fn facet(tris: f64[3][>..]&, a: f64[3], b: f64[3], c: f64[3]) {
    tris.push(a);
    tris.push(b);
    tris.push(c);
}

// Writes the triangle that starts at tris[3 * t] as one STL facet, with
// the unit normal of the right-hand rule.
fn write_facet(out: u8[>..]&, tris: f64[3][:], t: i64) {
    let a = tris[3 * t];
    let b = tris[3 * t + 1];
    let c = tris[3 * t + 2];
    let n = cross3(sub3(b, a), sub3(c, a));
    let l = norm3(n);
    let s = if l > 0.0 { 1.0 / l } else { 0.0 };
    format(out, "  facet normal ", n[0] * s, " ", n[1] * s, " ", n[2] * s, "\n    outer loop\n");
    for p in [a, b, c] { format(out, "      vertex ", p[0], " ", p[1], " ", p[2], "\n"); }
    out.append("    endloop\n  endfacet\n");
}

// The faces of the tetrahedra that are in one cell only, each with its
// nodes in the order that turns its normal out of the cell.
fn boundary_facets(out: f64[3][>..]&) {
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

// A real number of a STEP file: it has a point, and an exponent is E.
fn step_real(out: u8[>..]&, x: f64) {
    for ch in str(x) { out.push(if ch == 'e' { 'E' } else { ch }); }
}

// The STEP AP242 file of the triangles: one tessellated shell with one
// triangulated face, in the length unit of the header. A 2D Surface and a
// Solid use the same entities. Entity numbers: #1 to #25 are the product
// and its units, #30 to #33 the geometry.
fn write_step(out: u8[>..]&, name: u8[:], tris: f64[3][:]) {
    let unit = if metre == 1000.0 { ".MILLI." } else { "$" };
    format(out, "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('skein tessellated export'),'2;1');\n",
           "FILE_NAME('", name, ".stp','',(''),(''),'skein','skein','');\n",
           "FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nENDSEC;\nDATA;\n",
           "#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n",
           "#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2011,#1);\n",
           "#3=MECHANICAL_CONTEXT('',#1,'mechanical');\n",
           "#4=PRODUCT('", name, "','", name, "','',(#3));\n",
           "#5=PRODUCT_DEFINITION_FORMATION('','',#4);\n",
           "#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n",
           "#7=PRODUCT_DEFINITION('design','',#5,#6);\n",
           "#8=PRODUCT_DEFINITION_SHAPE('','',#7);\n",
           "#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(", unit, ",.METRE.));\n",
           "#11=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n",
           "#12=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n",
           "#13=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-9),#10,'distance_accuracy_value','');\n",
           "#14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))",
           "GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''));\n",
           "#15=CARTESIAN_POINT('',(0.,0.,0.));\n#16=DIRECTION('',(0.,0.,1.));\n#17=DIRECTION('',(1.,0.,0.));\n",
           "#18=AXIS2_PLACEMENT_3D('',#15,#16,#17);\n");
    format(out, "#30=COORDINATES_LIST('',", tris.len, ",(");
    for i in tris.len {
        if i > 0 { out.append(","); }
        out.append("(");
        for j in 3 {
            if j > 0 { out.append(","); }
            step_real(out, tris[i][j]);
        }
        out.append(")");
    }
    format(out, "));\n#31=TRIANGULATED_FACE('',#30,0,(),$,(),(");
    for t in tris.len / 3 {
        if t > 0 { out.append(","); }
        format(out, "(", 3 * t + 1, ",", 3 * t + 2, ",", 3 * t + 3, ")");
    }
    out.append("));\n#32=TESSELLATED_SHELL('',(#31),$);\n");
    format(out, "#33=TESSELLATED_SHAPE_REPRESENTATION('", name, "',(#32,#18),#14);\n",
           "#34=SHAPE_DEFINITION_REPRESENTATION(#8,#33);\nENDSEC;\nEND-ISO-10303-21;\n");
}

// `skein export`, with the STL in `out`. The result is the exit status
// (command_line.md §10). A fault of the model goes to standard error.
fn export_stl(args: u8[:][:], out: u8[>..]&) -> i64 {
    // The index in args of the file, of the geometry and of the size.
    var at = -1;
    var geo = -1;
    var sized = -1;
    var step = false;
    var i = 2;
    while i < args.len {
        let a = args[i];
        if (a == "--geometry" || a == "--format" || a == "--size") && i + 1 < args.len {
            let v = args[i + 1];
            if a == "--geometry" { geo = i + 1; }
            if a == "--size" { sized = i + 1; }
            guard a != "--format" || v == "stl" || v == "step" else {
                write_stderr(str("skein export: the format '", v, "' is not known: write --format stl or --format step\n"));
                return 2;
            }
            if a == "--format" { step = v == "step"; }
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
        let msg: u8[>..] = if sized >= 0 { str("skein export: the size '", args[sized], "' is not a length > 0, such as 2mm\n") }
                           else { "skein export: the model has no mesh block: give --size, such as --size 2mm\n" };
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
    var tris: f64[3][>..] = [];
    if dim == 3 {
        boundary_facets(tris);
    } else {
        for c in cell_count() { facet(tris, vertex(cell_node(c, 0)), vertex(cell_node(c, 1)), vertex(cell_node(c, 2))); }
    }
    if step {
        write_step(out, name, tris[..]);
        return 0;
    }
    format(out, "solid ", name, "\n");
    for t in tris.len / 3 { write_facet(out, tris[..], t); }
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
