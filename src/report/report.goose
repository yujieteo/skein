import std;
import os;
import base64;
namespace report;
import syntax from "../syntax";
import check from "../check";

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

let USAGE = "usage: skein report <file.skn> [--study <name>] [--no-run] [-o <file.html>]\n";

fn shell(command: u8[:]) -> i64 {
    var c: u8[>..] = [];
    c.append(command);
    c.push(0);
    let status = system(c[0]);
    if status < 0 || (status & 0x7f) != 0 { return 3; }
    return (status >> 8) & 0xff;
}

fn quoted(s: u8[:]) -> u8[>..] {
    var out: u8[>..] = "'";
    for c in s {
        if c == '\'' { out.append("'\\''"); } else { out.push(c); }
    }
    out.push('\'');
    return out;
}

fn make_dirs(path: u8[:]) -> bool {
    for c, i in path {
        if c == '/' && i > 0 { make_dir(path[..i]); }
    }
    return make_dir(path);
}

fn json_string(out: u8[>..]&, s: u8[:]) {
    out.push('"');
    for c, i in s {
        if c == '"' || c == '\\' {
            out.push('\\');
            out.push(c);
        } else if c == '\n' {
            out.append("\\n");
        } else if c < 32 {
            format(out, "\\u00", if c < 16 { "0" } else { "1" }, "0123456789abcdef"[c % 16..c % 16 + 1]);
        } else if c == '/' && i > 0 && s[i - 1] == '<' {
            out.append("\\/");
        } else {
            out.push(c);
        }
    }
    out.push('"');
}

fn json_field(out: u8[>..]&, key: u8[:], value: u8[:]) {
    json_string(out, key);
    out.push(':');
    json_string(out, value);
    out.push(',');
}

fn corner_type(t: i64, nodes: i64) -> i64, i64 {
    match t {
        21 => { return 3, 2; },
        22 => { return 5, 3; },
        23, 28 => { return 9, 4; },
        24 => { return 10, 4; },
        25, 29 => { return 12, 8; },
        26 => { return 13, 6; },
        _ => { return t, nodes; },
    }
}

fn attribute(tag: u8[:], key: u8[:]) -> const u8[:] {
    let at = find(tag, str(" ", key, "=\""));
    if at < 0 { return tag[..0]; }
    let from = at + key.len + 3;
    let end = find(tag[from..], "\"");
    return if end < 0 { tag[..0] } else { tag[from..from + end] };
}

fn numbers(s: u8[:], out: f64[>..]&) {
    var i = 0;
    while i < s.len {
        while i < s.len && is_space(s[i]) { i++; }
        let from = i;
        while i < s.len && !is_space(s[i]) { i++; }
        if i > from {
            let x, good = parse_flt(s[from..i]);
            out.push(if good { x } else { 0.0 });
        }
    }
}

struct Array {
    section: u8[..16],
    name: u8[..64],
    comps: i64,
    lo: i64,
    hi: i64,
}

var values: f64[>..];

fn arrays(text: u8[:]) -> Array[>..] {
    var out: Array[>..] = [];
    values.clear();
    var section: u8[..16] = [];
    var i = 0;
    while i < text.len {
        let lt = find(text[i..], "<");
        if lt < 0 { break; }
        i += lt;
        let gt = find(text[i..], ">");
        if gt < 0 { break; }
        let tag = text[i..i + gt + 1];
        i += gt + 1;
        for s in ["Points", "Cells", "PointData", "CellData"] {
            if starts_with(tag, str("<", s, ">")) || starts_with(tag, str("<", s, " ")) {
                section.clear();
                section.append(s);
            }
        }
        if !starts_with(tag, "<DataArray") { continue; }
        let end = find(text[i..], "</DataArray>");
        if end < 0 { break; }
        var a = Array { section: section, name: [], comps: 1, lo: values.len, hi: 0 };
        numbers(text[i..i + end], values);
        a.hi = values.len;
        let name = attribute(tag, "Name");
        a.name.append(name[..min(name.len, 64)]);
        let comps, good = parse_int(attribute(tag, "NumberOfComponents"));
        if good && comps > 0 { a.comps = comps; }
        out.push(a);
        i += end;
    }
    return out;
}

fn array_at(xs: Array[:], section: u8[:], name: u8[:]) -> i64 {
    for a, k in xs {
        if a.section[..] == section && (name.len == 0 || a.name[..] == name) { return k; }
    }
    return -1;
}

fn append_block(data: u8[>..]&, bytes: u8[:]) {
    let n = bytes.len;
    for k in 4 { data.push(((n >> (8 * k)) & 0xff) as u8); }
    data.append(bytes);
}

fn floats(xs: f64[:]) -> f32[>..] {
    var out: f32[>..] = [];
    for x in xs { out.push(x as f32); }
    return out;
}

fn array_tag(out: u8[>..]&, kind: u8[:], name: u8[:], comps: i64, offset: i64) {
    format(out, "<DataArray type=\"", kind, "\"");
    if name.len > 0 { format(out, " Name=\"", name, "\""); }
    format(out, " NumberOfComponents=\"", comps, "\" format=\"appended\" offset=\"", offset, "\"/>\n");
}

let INT16_MAX = 32767.0;

fn shorts(xs: f64[:]) -> i16[>..], f64 {
    var m = 0.0;
    for x in xs { m = max(m, abs(x)); }
    let scale = if m > 0.0 { m / INT16_MAX } else { 1.0 };
    var out: i16[>..] = [];
    for x in xs {
        let half = if x >= 0.0 { 0.5 } else { -0.5 };
        out.push((x / scale + half) as! i16);
    }
    return out, scale;
}

fn short_array(out: u8[>..]&, data: u8[>..]&, name: u8[:], comps: i64, xs: f64[:]) {
    let n, scale = shorts(xs);
    format(out, "<DataArray type=\"Int16\" Name=\"", name, "\" NumberOfComponents=\"", comps,
           "\" Scale=\"", scale, "\" format=\"appended\" offset=\"", data.len, "\"/>\n");
    append_block(data, bytes_of(n));
}

fn vtu_open(out: u8[>..]&, points: i64, cells: i64) {
    format(out, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
           "<UnstructuredGrid>\n<Piece NumberOfPoints=\"", points, "\" NumberOfCells=\"", cells, "\">\n");
}

fn vtu_close(out: u8[>..]&, data: u8[:]) {
    out.append("</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_");
    out.append(data);
    out.append("\n</AppendedData>\n</VTKFile>\n");
}

fn small_vtu(text: u8[:], mesh: u8[>..]&, fields: u8[>..]&) -> bool {
    let xs = arrays(text);
    let pa = array_at(xs, "Points", "");
    let ca = array_at(xs, "Cells", "connectivity");
    let oa = array_at(xs, "Cells", "offsets");
    let ta = array_at(xs, "Cells", "types");
    guard pa >= 0 && ca >= 0 && oa >= 0 && ta >= 0 else { return false; };
    let conn = values[xs[ca].lo..xs[ca].hi];
    let offsets = values[xs[oa].lo..xs[oa].hi];
    let pts = values[xs[pa].lo..xs[pa].hi];
    let cells = offsets.len;
    var keep: i64[>..] = [];
    for n in pts.len / 3 { keep.push(-1); }
    var old: i64[>..] = [];
    var nodes: i32[>..] = [];
    var ends: i32[>..] = [];
    var types: u8[>..] = [];
    var from = 0;
    for c in cells {
        let to = offsets[c] as i64;
        let t, corners = corner_type(values[xs[ta].lo + c] as i64, to - from);
        for k in min(corners, to - from) {
            let n = conn[from + k] as i64;
            guard n >= 0 && n < keep.len else { return false; };
            if keep[n] < 0 {
                keep[n] = old.len;
                old.push(n);
            }
            nodes.push(keep[n] as i32);
        }
        ends.push(nodes.len as i32);
        types.push(t as u8);
        from = to;
    }
    var data: u8[>..] = [];
    vtu_open(mesh, old.len, cells);
    var points: f64[>..] = [];
    for n in old { for k in 3 { points.push(pts[3 * n + k]); } }
    mesh.append("<Points>\n");
    array_tag(mesh, "Float32", "", 3, data.len);
    append_block(data, bytes_of(floats(points)));
    mesh.append("</Points>\n<Cells>\n");
    array_tag(mesh, "Int32", "connectivity", 1, data.len);
    append_block(data, bytes_of(nodes));
    array_tag(mesh, "Int32", "offsets", 1, data.len);
    append_block(data, bytes_of(ends));
    array_tag(mesh, "UInt8", "types", 1, data.len);
    append_block(data, types);
    mesh.append("</Cells>\n");
    vtu_close(mesh, data);
    small_fields(xs, old, cells, keep.len, fields);
    return true;
}

fn small_fields(xs: Array[:], old: i64[:], cells: i64, nodes: i64, out: u8[>..]&) {
    var data: u8[>..] = [];
    vtu_open(out, old.len, cells);
    out.append("<PointData>\n");
    for a in xs {
        guard a.section[..] == "PointData" && a.hi - a.lo == a.comps * nodes else { continue; };
        var v: f64[>..] = [];
        for n in old { for k in a.comps { v.push(values[a.lo + a.comps * n + k]); } }
        short_array(out, data, a.name, a.comps, v);
    }
    out.append("</PointData>\n<CellData>\n");
    for a in xs {
        guard a.section[..] == "CellData" && a.hi - a.lo == a.comps * cells else { continue; };
        short_array(out, data, a.name, a.comps, values[a.lo..a.hi]);
    }
    out.append("</CellData>\n");
    vtu_close(out, data);
}

fn result_files(toon: u8[:]) -> const u8[:][>..] {
    var out: const u8[:][>..] = [];
    var inside = false;
    for line in split(toon, '\n') {
        if starts_with(line, "results[") {
            inside = true;
        } else if inside && starts_with(line, "  ") {
            out.push(trim(line));
        } else {
            inside = false;
        }
    }
    return out;
}

fn names(prog: i64) -> u8[>..], u8[>..] {
    var model: u8[>..] = [];
    var study: u8[>..] = [];
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        let k = syntax::nodes[c].kind;
        let word = syntax::text(syntax::nodes[c].tok);
        if k == syntax::N_HEADER { model.append(word); }
        if k == syntax::N_STUDY && study.len == 0 { study.append(word); }
        c = syntax::nodes[c].next;
    }
    return model, study;
}

fn text_of(path: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    read_file(path, out);
    return out;
}

var json: u8[>..];
var blobs: u8[>..];
var files = 0;
var mesh_texts: u8[>..];
var mesh_ends: i64[>..];
var mesh_blobs: i64[>..];
var cases = 0;

fn add_result(work: u8[:], vtu: u8[:]) -> i64, i64 {
    var mesh: u8[>..] = [];
    var fields: u8[>..] = [];
    guard small_vtu(text_of(vtu), mesh, fields) else {
        write_stderr(str("skein: cannot read the VTU file '", vtu, "'\n"));
        return -1, -1;
    };
    var m = -1;
    for end, k in mesh_ends {
        let start = if k > 0 { mesh_ends[k - 1] } else { 0 };
        if mesh_texts[start..end] == mesh { m = mesh_blobs[k]; }
    }
    if m < 0 {
        m = add_blob(work, mesh);
        mesh_texts.append(mesh);
        mesh_ends.push(mesh_texts.len);
        mesh_blobs.push(m);
    }
    return m, add_blob(work, fields);
}

fn add_blob(work: u8[:], bytes: u8[:]) -> i64 {
    guard write_blob(work, bytes) else { return -1; };
    let k = files;
    let path = str(work, "/", k, ".vtu");
    format(blobs, "<script type=\"text/plain\" id=\"vtu-", k, "\">");
    base64::encode(blobs, text_of(str(path, ".gz")));
    blobs.append("</script>\n");
    files++;
    return k;
}

fn write_blob(work: u8[:], bytes: u8[:]) -> bool {
    let path = str(work, "/", files, ".vtu");
    guard write_file(path, bytes) && shell(str("gzip -9 -n -f ", quoted(path), " -c > ", quoted(str(path, ".gz")))) == 0 else {
        write_stderr(str("skein: gzip cannot compress '", path, "'\n"));
        return false;
    };
    return true;
}

fn add_case(work: u8[:], root: u8[:], id: i64, toon: u8[:]) -> bool {
    if cases > 0 { json.push(','); }
    cases++;
    format(json, "{\"id\":", id, ",");
    json_field(json, "run", toon);
    json.append("\"results\":[");
    for vtu, i in result_files(toon) {
        let m, k = add_result(work, str(root, vtu));
        if m < 0 || k < 0 { return false; }
        format(json, if i > 0 { "," } else { "" }, "{\"mesh\":", m, ",\"blob\":", k, ",");
        json_field(json, "path", vtu);
        json.pop();
        json.push('}');
    }
    json.append("]}");
    return true;
}

fn case_count(toon: u8[:]) -> i64 {
    let at = find(toon, "cases[");
    guard at >= 0 else { return 0; };
    let end = find(toon[at..], "]");
    let n, good = parse_int(toon[at + 6..at + max(end, 6)]);
    return if good { n } else { 0 };
}

fn last_run(model: u8[:], study: u8[:], toon: u8[:]) -> i64 {
    let file = str(".skein/studies/", model, "/", study);
    let last = text_of(str(file, ".toon"));
    guard study.len > 0 && last.len > 0 && write_file(toon, last) else {
        write_stderr(str("skein report: --no-run needs the table of a study run in .skein/studies/", model, "\n"));
        return 2;
    };
    guard text_of(str(file, ".hash")) == syntax::source_hash() else {
        write_stderr("skein report: the results are of an older source; run skein report without --no-run\n");
        return 2;
    };
    return if find(last, ",fail,") >= 0 { 1 } else { 0 };
}

fn add_checks(path: u8[:], work: u8[:]) {
    var me: u8[>..] = [];
    arg(0, me);
    for key in ["check", "verify"] {
        let file = str(work, "/", key, ".toon");
        shell(str(quoted(me), " ", key, " ", quoted(path), " > ", quoted(file), " 2>&1"));
        json_field(json, key, text_of(file));
    }
}

fn run_model(path: u8[:], model: u8[:], study: u8[:], work: u8[:], no_run: bool) -> i64 {
    var me: u8[>..] = [];
    arg(0, me);
    let skein = quoted(me);
    let toon = str(work, "/run.toon");
    let study_arg = if study.len > 0 { str(" --study ", quoted(study)) } else { str("") };
    var status = 0;
    if no_run {
        status = last_run(model, study, toon);
    } else {
        status = shell(str(skein, " run ", quoted(path), study_arg, " > ", quoted(toon), " 2> ", quoted(str(work, "/run.err"))));
        write_stderr(text_of(str(work, "/run.err")));
    }
    if status > 1 { return status; }
    add_checks(path, work);
    json_field(json, "study", study);
    let table: u8[>..] = if study.len > 0 { text_of(toon) } else { "" };
    json_field(json, "table", table);
    json.append("\"cases\":[");
    if study.len == 0 {
        if !add_case(work, "", 1, text_of(toon)) { return 3; }
    } else {
        for c in 1..case_count(text_of(toon)) + 1 {
            let case_toon = str(".skein/studies/", model, "/", study, "/", c, "/run.toon");
            if !add_case(work, "", c, text_of(case_toon)) { return 3; }
        }
    }
    json.append("]}");
    return status;
}

fn page_solve(path: u8[:], model: u8[:], work: u8[:]) -> u8[>..], bool {
    var me: u8[>..] = [];
    arg(0, me);
    let out = str(work, "/build.out");
    let status = shell(str(quoted(me), " build ", quoted(path), " --target wasm > ", quoted(out), " 2> ", quoted(str(work, "/build.err"))));
    guard status == 0 else {
        write_stderr(text_of(str(work, "/build.err")));
        return "", false;
    };
    return text_of(str(".skein/build/", model, "/web/solve.html")), find(text_of(out), "solve_in_page: yes") >= 0;
}

fn put(doc: u8[>..]&, mark: u8[:], with: u8[:]) {
    let at = find(doc, mark);
    guard at >= 0 else { return; };
    var out: u8[>..] = [];
    out.append(doc[..at]);
    out.append(with);
    out.append(doc[at + mark.len..]);
    doc.clear();
    doc.append(out);
}

fn token_classes(out: u8[>..]&) {
    out.append("[");
    var first = true;
    var prev = -1;
    for i in 0..syntax::tokens.len {
        let t = syntax::tokens[i];
        let kind = t.kind;
        var class = "";
        if t.file == 0 {
            if kind >= syntax::KW_ALLOW { class = "k"; }
            else if kind == syntax::TK_NUMBER || kind == syntax::TK_PERCENT { class = "n"; }
            else if kind == syntax::TK_STRING { class = "s"; }
            else if kind == syntax::TK_IDENT {
                let word = syntax::text(i);
                if (prev == syntax::TK_NUMBER || prev == syntax::TK_PERCENT) && syntax::is_unit_name(word) { class = "u"; }
                else if (prev == syntax::TK_COLON || prev == syntax::TK_LANGLE) && word[0] >= 'A' && word[0] <= 'Z' { class = "t"; }
            }
        }
        if kind != syntax::TK_COMMENT && kind != syntax::TK_DOC_COMMENT { prev = kind; }
        if class.len == 0 { continue; }
        var chars = 0;
        for k in 0..t.len {
            if (syntax::source[t.start + k] & 0xC0) != 0x80 { chars++; }
        }
        format(out, if first { "" } else { "," }, "[", t.line, ",", t.col - 1, ",", chars, ",\"", class, "\"]");
        first = false;
    }
    out.append("]");
}

fn page(model: u8[:], solve: u8[:]) -> u8[>..] {
    let dir = str(resource_dir(), "../src/report/");
    var html = text_of(str(dir, "page.html"));
    guard html.len > 0 else {
        write_stderr(str("skein: cannot read '", dir, "page.html'\n"));
        return html;
    };
    put(html, "<!--title-->", model);
    var css = text_of(str(dir, "kowalski.css"));
    css.append(text_of(str(dir, "viewer.css")));
    var inter: u8[>..] = [];
    base64::encode(inter, text_of(str(dir, "fonts/inter.woff2")));
    put(css, "/fonts/inter.woff2", str("data:font/woff2;base64,", inter));
    var plex: u8[>..] = [];
    base64::encode(plex, text_of(str(dir, "fonts/ibm-plex-mono.woff2")));
    put(css, "/fonts/ibm-plex-mono.woff2", str("data:font/woff2;base64,", plex));
    put(html, "/*css*/", css);
    put(html, "/*js*/", text_of(str(dir, "viewer.js")));
    put(html, "<!--data-->", str("<script type=\"application/json\" id=\"skein-data\">", json, "</script>\n", blobs, solve));
    return html;
}

fn web_page(path: u8[:], text: u8[:], model: u8[:], dir: u8[:], toon: u8[:], solve: u8[:]) -> u8[>..] {
    syntax::parse_model(path, text);
    var classes: u8[>..] = [];
    token_classes(classes);
    json.push('{');
    json_field(json, "model", model);
    json_field(json, "file", path);
    json_field(json, "source", text);
    format(json, "\"tokens\":", classes, ",");
    add_checks(path, dir);
    for key in ["study", "table"] { json_field(json, key, ""); }
    json.append("\"cases\":[");
    guard add_case(dir, str(dir, "/"), 1, toon) else { return ""; };
    json.append("]}");
    return page(model, solve);
}

fn run<F>(args: const u8[:][:]) -> i64 {
    var path: u8[>..] = [];
    var out_path: u8[>..] = [];
    var study: u8[>..] = [];
    var no_run = false;
    var i = 0;
    while i < args.len {
        let more = i + 1 < args.len;
        if args[i] == "-o" && more {
            out_path.append(args[i + 1]);
            i += 2;
        } else if args[i] == "--no-run" {
            no_run = true;
            i++;
        } else if args[i] == "--study" && more {
            study.append(args[i + 1]);
            i += 2;
        } else if args[i].len > 0 && args[i][0] != '-' && path.len == 0 {
            path.append(args[i]);
            i++;
        } else {
            write_stderr(str("skein report: give one model file\n", USAGE));
            return 2;
        }
    }
    guard path.len > 0 else {
        write_stderr(str("skein report: give one model file\n", USAGE));
        return 2;
    };
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein report: cannot read '", path, "'\n"));
        return 1;
    };
    check::check_text(path, text);
    if check::error_count() > 0 { return check::finish(false); }
    let prog = syntax::parse_model(path, text);
    let model, first_study = names(prog);
    let takes = F(prog);
    var classes: u8[>..] = [];
    token_classes(classes);
    if study.len == 0 { study.append(first_study); }
    if out_path.len == 0 { out_path.append(str(model, ".html")); }
    let work = str(".skein/report/", model);
    guard make_dirs(work) else {
        write_stderr("skein report: cannot make the directory .skein/report\n");
        return 1;
    };
    json.push('{');
    json_field(json, "model", model);
    json_field(json, "file", path);
    json_field(json, "source", text);
    format(json, "\"tokens\":", classes, ",");
    let status = run_model(path, model, study, work, no_run);
    if status > 1 { return status; }
    var solve: u8[>..] = [];
    var solve_in_page = false;
    if takes {
        let page, in_page = page_solve(path, model, work);
        solve.append(page);
        solve_in_page = in_page;
    }
    let html = page(model, solve);
    guard html.len > 0 else { return 3; };
    guard write_file(out_path, html) else {
        write_stderr(str("skein report: cannot write '", out_path, "'\n"));
        return 1;
    };
    print("report: ", out_path);
    print("bytes: ", html.len);
    print("cases: ", cases);
    print("solve_in_page: ", if solve_in_page { "yes" } else { "no" });
    return status;
}
