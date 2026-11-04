import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import mesher from "../mesher";
import .gen;

fn emit_pml(out: u8[>..]&) {
    guard pml_g >= 0 else { return; }
    format(out, "    elements::pml(", pml_g, ", ");
    lits(out, pml_lo);
    out.append(", ");
    lits(out, pml_hi);
    out.append(", ");
    lit(out, pml_h);
    out.append(");\n");
}

fn wave_region_terms(w: Law) -> f64[4] {
    let base = material_node;
    material_node = w.item;
    let c = if wave_c_node >= 0 { cval(wave_c_node) } else { [0.0, 0.0] };
    let src = if wave_s_node >= 0 && !mentions_xyz(wave_s_node) { cval(wave_s_node) } else { [0.0, 0.0] };
    material_node = base;
    return [c[0], c[1], src[0], src[1]];
}

fn emit_wave_laws(out: u8[>..]&) {
    for w in region_laws {
        let t = wave_region_terms(w);
        format(out, "    cells(elements::wave_in(", w.g, ", 1.0");
        for v in t {
            out.append(", ");
            lit(out, v);
        }
        out.append("), ");
        quoted(out, w.label);
        out.append(");\n");
    }
    guard wave_s_node >= 0 && mentions_xyz(wave_s_node) else { return; }
    let base = material_node;
    format(out, "    elements::wave_source_in(-1) { p => ");
    emit_complex(out, wave_s_node, true);
    out.append(" };\n");
    for w in region_laws {
        material_node = w.item;
        format(out, "    elements::wave_source_in(", w.g, ") { p => ");
        emit_complex(out, wave_s_node, true);
        out.append(" };\n");
    }
    material_node = base;
}

recursive fn emit_complex(out: u8[>..]&, n: i64, spatial: bool) {
    let k = kind(n);
    let op = if k == syntax::N_BINARY || k == syntax::N_UNARY { syntax::text(tok(n)) } else { "" };
    if !(if spatial { mentions_xyz(n) } else { !ci_free(n) }) {
        let c = cval(n);
        out.append("[");
        lit(out, c[0]);
        out.append(", ");
        lit(out, c[1]);
        out.append("]");
        return;
    }
    if spatial && k == syntax::N_NAME && axis_of(name(n)) >= 0 {
        format(out, "[p[", axis_of(name(n)), "], 0.0]");
        return;
    }
    if !spatial && k == syntax::N_NAME && name(n) == unknown_name {
        out.append("[s[0], s[1]]");
        return;
    }
    if k == syntax::N_NAME && codegen::depth < 100 {
        let v = if codegen::find_value(name(n)) >= 0 { codegen::find_value(name(n)) } else { property(name(n)) };
        if v >= 0 {
            codegen::depth++;
            emit_complex(out, v, spatial);
            codegen::depth--;
            return;
        }
    }
    if k == syntax::N_UNARY && op == "-" {
        out.append("wsub([0.0, 0.0], ");
        emit_complex(out, first(n), spatial);
        out.append(")");
        return;
    }
    if !spatial && is_call(n, "dot") && syntax::child_count(n) == 3 {
        var g = -1;
        var m = -1;
        if is_call(arg(n, 0), "grad") && source_of(arg(arg(n, 0), 0)) == unknown_name && source_of(arg(n, 1)) == "n" { g = arg(n, 0); m = arg(n, 1); }
        if is_call(arg(n, 1), "grad") && source_of(arg(arg(n, 1), 0)) == unknown_name && source_of(arg(n, 0)) == "n" { g = arg(n, 1); m = arg(n, 0); }
        if g >= 0 && m >= 0 {
            if dim == 3 { unsupported(tok(n), "dot with the normal n on a curve in 3d: the curve has no normal", "integrate over a line or a circle in a 2d model"); }
            out.append("[s[2] * n[0] + s[3] * n[1] + s[4] * n[2], s[5] * n[0] + s[6] * n[1] + s[7] * n[2]]");
            return;
        }
    }
    if k == syntax::N_CALL && syntax::child_count(n) == 2 && kind(child(n, 0)) == syntax::N_NAME {
        let f = name(child(n, 0));
        if spatial && position(["sin", "cos", "exp"], f) >= 0 && !mentions_im(arg(n, 0)) {
            out.append("[");
            emit_space(out, n, -1);
            out.append(", 0.0]");
            return;
        }
        if !spatial && f == "conj" {
            out.append("wconj(");
            emit_complex(out, arg(n, 0), spatial);
            out.append(")");
            return;
        }
        if !spatial {
            for w, i in PARTS {
                if f == w {
                    format(out, "[wave_pick(", i, ", ");
                    emit_complex(out, arg(n, 0), spatial);
                    out.append("), 0.0]");
                    return;
                }
            }
        }
    }
    guard op == "+" || op == "-" || op == "*" || op == "/" else {
        if spatial {
            unsupported(tok(n), str("the expression '", source_of(n), "' is not in this generator"),
                        "write + - * / of numbers, im, lets, properties, x, y, z, and sin, cos or exp of a real function of x, y and z");
        } else {
            unsupported(tok(n), str("the integrand '", source_of(n), "' is not in this generator"),
                        str("write the unknown ", unknown_name, ", dot(grad(", unknown_name, "), n), re, imag, abs, arg and conj of them, joined by + - * / and constants"));
        }
        out.append("[0.0, 0.0]");
        return;
    };
    out.append(if op == "+" { "wadd(" } else if op == "-" { "wsub(" } else if op == "*" { "wmul(" } else { "wdiv(" });
    emit_complex(out, child(n, 0), spatial);
    out.append(", ");
    emit_complex(out, child(n, 1), spatial);
    out.append(")");
}

recursive fn mentions_im(n: i64) -> bool {
    if kind(n) == syntax::N_NAME && name(n) == "im" { return true; }
    if kind(n) == syntax::N_NAME && codegen::find_value(name(n)) >= 0 && codegen::depth < 100 {
        codegen::depth++;
        let r = mentions_im(codegen::find_value(name(n)));
        codegen::depth--;
        return r;
    }
    var c = first(n);
    while c >= 0 {
        if mentions_im(c) { return true; }
        c = next(c);
    }
    return false;
}

struct Seg {
    kind: i64,
    a: f64[3],
    b: f64[3],
    r: f64,
}

fn curve_value(n: i64) -> i64 {
    guard kind(n) == syntax::N_NAME else { return -1; }
    let v = codegen::find_value(name(n));
    if v >= 0 && is_call(v, "line") && syntax::child_count(v) == 3 { return v; }
    for l in surface_lets {
        if name(l) != name(n) { continue; }
        let b = child(l, 1);
        if is_call(b, "boundary") && syntax::child_count(b) == 2 && (is_call(arg(b, 0), "circle") || is_call(arg(b, 0), "rect")) { return b; }
    }
    return -1;
}

fn is_curve(n: i64) -> bool { return curve_value(n) >= 0; }

fn curve_segments(n: i64) -> Seg[>..] {
    var segs: Seg[>..] = [];
    let v = curve_value(n);
    guard v >= 0 else { return segs; }
    if is_call(v, "line") {
        segs.push(Seg { kind: 0, a: value(arg(v, 0)).v, b: value(arg(v, 1)).v, r: 0.0 });
        return segs;
    }
    let shape = arg(v, 0);
    if is_call(shape, "circle") {
        segs.push(Seg { kind: 1, a: value(arg(shape, 0)).v, b: [0.0; 3], r: scalar_value(arg(shape, 1)) });
        return segs;
    }
    let lo = value(arg(shape, 0)).v;
    let hi = value(arg(shape, 1)).v;
    let c: f64[3][4] = [[lo[0], lo[1], 0.0], [hi[0], lo[1], 0.0], [hi[0], hi[1], 0.0], [lo[0], hi[1], 0.0]];
    for i in 4 { segs.push(Seg { kind: 0, a: c[i], b: c[(i + 1) % 4], r: 0.0 }); }
    return segs;
}

recursive fn ci_free(e: i64) -> bool {
    if kind(e) == syntax::N_NAME && (name(e) == unknown_name || name(e) == "n") { return false; }
    var c = first(e);
    while c >= 0 {
        if !ci_free(c) { return false; }
        c = next(c);
    }
    return true;
}


fn emit_curve_integral(pre: u8[>..]&, n: i64, z: i64) {
    let segs = curve_segments(arg(z, 1));
    format(pre, "        var wz", n, ": f64[2] = [0.0, 0.0];\n        var wf", n, " = true;\n");
    for sg in segs {
        var length = sg.r * 2.0 * PI;
        if sg.kind == 0 {
            length = 0.0;
            for t in 3 { length += (sg.b[t] - sg.a[t]) * (sg.b[t] - sg.a[t]); }
            length = sqrt(length);
        }
        let count = max(16, ceil(2.0 * length / mesh_size) as! i64);
        format(pre, "        {\n            let ci, ok = elements::curve_integral(", sg.kind, ", [");
        lits(pre, sg.a);
        pre.append("], [");
        lits(pre, sg.b);
        pre.append("], ");
        lit(pre, sg.r);
        format(pre, ", ", count, ", false) { p, n, s => ");
        emit_complex(pre, arg(z, 0), false);
        format(pre, " };\n            wz", n, " = wadd(wz", n, ", ci);\n            if !ok { wf", n, " = false; }\n        }\n");
    }
}
