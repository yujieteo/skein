// `skein check` and `skein lint` (docs/command_line.md §2, §4, §6): parses
// a model, applies the rules of spec §19 and prints the diagnostics and the
// free boundary report and the default methods in TOON or JSON.
import std;
import os;
import dictionary;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .expr;
import .items;
import .physics;
import .steps;
import .templates;
import .lint;

// Checks the model `text` from the file `path`. The result is in `diags`
// and `free`.
fn check_text(path: u8[:], text: u8[:]) {
    syntax::reset();
    reset_state();
    objectives.clear();
    high_derivs.clear();
    globals.clear();
    std_prog = -1;
    let f = syntax::add_file(path, text);
    syntax::lex(f);
    let prog = syntax::parse(f);
    model_prog = prog;
    syntax::read_uses(prog);
    take_syntax_errors();
    // The other rules need a tree with no syntax errors (spec §19).
    if diags.len > 0 { return; }
    let s = syntax::add_file(STD_PATH, STD_SOURCE);
    syntax::lex(s);
    std_prog = syntax::parse(s);
    // The header first: it sets the space and the units of the built-in names.
    let h = syntax::child_of_kind(prog, syntax::N_HEADER);
    if h >= 0 { check_header(h); }
    add_builtins();
    note_declarations(prog);
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        check_item(item);
        item = syntax::nodes[item].next;
    }
    scales_end();
    objectives_end();
    method_missing();
    crack_tip_order();
    free_boundary();
    sort_diags();
}

// Adds each error of `syntax::errors` to `diags`.
fn take_syntax_errors() {
    for e in syntax::errors {
        let k = syntax::tokens[e.tok];
        var d = Diag { id: e.id, file: [], fnum: k.file, line: k.line, col: k.col, msg: e.msg, fix: e.fix };
        let file = syntax::files[k.file].path;
        d.file.append(file[..min(file.len, 256)]);
        diags.push(d);
    }
}

// `skein lint`: checks the model `text`, then applies the lint rules when
// it has no error. The result is in `diags`.
fn lint_text(path: u8[:], text: u8[:]) {
    check_text(path, text);
    if diags.len > 0 { return; }
    lint_rules();
    sort_diags();
}

fn check_item(n: i64) {
    match syntax::nodes[n].kind {
        syntax::N_PARAM => { check_value_decl(n, S_PARAM); },
        syntax::N_LET => { check_value_decl(n, S_LET); },
        syntax::N_FN => { check_fn(n); },
        syntax::N_TABLE => { check_table(n); },
        syntax::N_MATERIAL => { check_material(n); },
        syntax::N_MESH => { check_discretization(n, S_MESH); },
        syntax::N_POINTS => { check_discretization(n, S_POINTS); },
        syntax::N_PHYSICS => { check_physics(n); },
        syntax::N_METHOD => { check_method(n); },
        syntax::N_BC => { check_bc(n); },
        syntax::N_CONTACT, syntax::N_INTERFACE => { check_pair(n); },
        syntax::N_SYMMETRY => { check_symmetry(n); },
        syntax::N_STEP => { check_step(n); },
        syntax::N_EXPECT => { check_expect(n); },
        syntax::N_ALLOW => { check_allow(n); },
        syntax::N_SCALES => { check_scales(n); },
        syntax::N_STUDY => { check_study(n); },
        syntax::N_MANUFACTURE => { check_manufacture(n); },
        syntax::N_TOLERANCE => { check_tolerance(n); },
        syntax::N_TEMPLATE => {
            error(syntax::nodes[n].tok, "duplicate_name", "only the standard library declares templates",
                  "write the physics in full");
        },
        _ => {},
    }
}

// Sorts the diagnostics by file, line and column (command_line.md §6.1),
// stably. The model file comes first, then each library file in `use` order.
fn after(a: Diag, b: Diag) -> bool {
    if a.fnum != b.fnum { return a.fnum > b.fnum; }
    if a.line != b.line { return a.line > b.line; }
    return a.col > b.col;
}

fn sort_diags() {
    for i in 1..diags.len {
        var j = i;
        while j > 0 && after(diags[j - 1], diags[j]) {
            let tmp = diags[j - 1];
            diags[j - 1] = diags[j];
            diags[j] = tmp;
            j--;
        }
    }
}

// --- output (docs/command_line.md §6) -------------------------------------------------------

let HEX = "0123456789abcdef";

// `s` in double quotes, with `"`, `\` and the control bytes escaped (TOON and JSON).
fn quoted(out: u8[>..]&, s: u8[:]) {
    out.push('"');
    for c in s {
        if c == '"' || c == '\\' {
            out.push('\\');
            out.push(c);
        } else if c == '\n' {
            out.append("\\n");
        } else if c == '\t' {
            out.append("\\t");
        } else if c == '\r' {
            out.append("\\r");
        } else if c < 0x20 {
            out.append("\\u00");
            out.push(HEX[c / 16]);
            out.push(HEX[c % 16]);
        } else {
            out.push(c);
        }
    }
    out.push('"');
}

// A TOON value: quoted when it holds a comma, a quote or a space at an end.
fn toon_value(out: u8[>..]&, s: u8[:]) {
    var plain = s.len > 0 && s[0] != ' ' && s[s.len - 1] != ' ';
    for c in s {
        if c == ',' || c == '"' || c == '\\' || c == ':' { plain = false; }
    }
    if plain { out.append(s); } else { quoted(out, s); }
}

// The number of diagnostics that are errors, not lint warnings.
fn error_count() -> i64 {
    var n = 0;
    for d in diags {
        if !d.lint { n++; }
    }
    return n;
}

fn write_toon(out: u8[>..]&) {
    if nondim { out.append("units: nondim\n"); }
    format(out, "diagnostics[", diags.len, "]{class,id,file,line,col,msg,fix}:\n");
    for d in diags {
        out.append(if d.lint { "  lint," } else { "  error," });
        out.append(d.id);
        out.push(',');
        toon_value(out, d.file);
        format(out, ",", d.line, ",", d.col, ",");
        quoted(out, d.msg);
        out.push(',');
        quoted(out, d.fix);
        out.push('\n');
    }
    let errors = error_count();
    format(out, "summary: errors=", errors, " lints=", diags.len - errors, "\n");
    if errors > 0 { return; }
    format(out, "free_boundary[", free.len, "]{physics,field,region,condition}:\n");
    for f in free {
        out.append("  ");
        toon_value(out, f.physics);
        out.push(',');
        toon_value(out, f.field);
        out.push(',');
        quoted(out, f.region);
        out.push(',');
        quoted(out, f.condition);
        out.push('\n');
    }
    format(out, "method[", chosen.len, "]{physics,kind,mesh,order,stabilize}:\n");
    for c in chosen {
        out.append("  ");
        toon_value(out, c.physics);
        out.append(",cg,");
        toon_value(out, c.mesh);
        out.push(',');
        toon_value(out, c.order);
        out.push(',');
        out.append(c.stabilize);
        out.push('\n');
    }
}

fn write_json(out: u8[>..]&) {
    out.append("{");
    if nondim { out.append("\"units\":\"nondim\","); }
    out.append("\"diagnostics\":[");
    for d, i in diags {
        if i > 0 { out.push(','); }
        out.append(if d.lint { "{\"class\":\"lint\",\"id\":" } else { "{\"class\":\"error\",\"id\":" });
        quoted(out, d.id);
        out.append(",\"file\":");
        quoted(out, d.file);
        format(out, ",\"line\":", d.line, ",\"col\":", d.col, ",\"msg\":");
        quoted(out, d.msg);
        out.append(",\"fix\":");
        quoted(out, d.fix);
        out.push('}');
    }
    let errors = error_count();
    format(out, "],\"summary\":{\"errors\":", errors, ",\"lints\":", diags.len - errors, "}");
    if errors == 0 {
        out.append(",\"free_boundary\":[");
        for f, i in free {
            if i > 0 { out.push(','); }
            out.append("{\"physics\":");
            quoted(out, f.physics);
            out.append(",\"field\":");
            quoted(out, f.field);
            out.append(",\"region\":");
            quoted(out, f.region);
            out.append(",\"condition\":");
            quoted(out, f.condition);
            out.push('}');
        }
        out.append("],\"method\":[");
        for c, i in chosen {
            if i > 0 { out.push(','); }
            out.append("{\"physics\":");
            quoted(out, c.physics);
            out.append(",\"kind\":\"cg\",\"mesh\":");
            quoted(out, c.mesh);
            out.append(",\"order\":");
            quoted(out, c.order);
            out.append(",\"stabilize\":");
            quoted(out, c.stabilize);
            out.push('}');
        }
        out.push(']');
    }
    out.append("}\n");
}

// `skein check`, and `skein lint` with no mesh. `args[1]` is the command.
// The result is the exit status (command_line.md §10): a lint warning
// alone gives 0.
fn command(args: u8[:][:]) -> i64 {
    let status, json, mesh = start(args);
    if status >= 0 { return status; }
    return finish(json);
}

// Reads the options of `args` and the model, then checks it, or lints it
// when `args[1]` is "lint". The result is an exit status, or -1 when the
// command continues, then true for `--format json` and for `--mesh`. Only
// `lint` has `--mesh`.
fn start(args: u8[:][:]) -> i64, bool, bool {
    let lint = args[1] == "lint";
    let usage = str("usage: skein ", args[1], if lint { " [--mesh]" } else { "" }, " [--format toon|json] model.skn\n");
    var json = false;
    var mesh = false;
    var at = -1;
    var i = 2;
    while i < args.len {
        let a = args[i];
        if a == "--format" && i + 1 < args.len {
            if args[i + 1] == "json" {
                json = true;
            } else if args[i + 1] != "toon" {
                write_stderr(usage);
                return 2, json, mesh;
            }
            i += 2;
            continue;
        }
        if a == "--mesh" && lint {
            mesh = true;
            i++;
            continue;
        }
        if (a.len > 0 && a[0] == '-') || at >= 0 {
            write_stderr(usage);
            return 2, json, mesh;
        }
        at = i;
        i++;
    }
    guard at >= 0 else {
        write_stderr(usage);
        return 2, json, mesh;
    };
    let path = args[at];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein: cannot read '", path, "'\n"));
        return 2, json, mesh;
    };
    if lint { lint_text(path, text); } else { check_text(path, text); }
    return -1, json, mesh;
}

// Prints the diagnostics, in JSON when `json` is true. The result is the
// exit status.
fn finish(json: bool) -> i64 {
    var out: u8[>..] = [];
    if json { write_json(out); } else { write_toon(out); }
    write_stdout(out);
    return if error_count() > 0 { 1 } else { 0 };
}
