// Tests of the syntax tree, the lexer and the parser. `./build.sh test` runs
// this program from the root of the repository, and a non-zero exit status
// is a failure.
import std;
import os;
import .tree;
import .lexer;
import .parser;

var failures = 0;

fn expect_eq(what: u8[:], got: u8[:], want: u8[:]) {
    if got == want { return; }
    failures++;
    print("FAIL ", what, "\n  got:  ", got, "\n  want: ", want);
}

fn expect_int(what: u8[:], got: i64, want: i64) {
    if got == want { return; }
    failures++;
    print("FAIL ", what, ": got ", got, ", want ", want);
}

// Adds a token of kind `kind` at bytes `start` to `start + len` of the last file.
fn tok(kind: i64, start: i64, len: i64) -> i64 {
    syntax::tokens.push(syntax::Token { kind: kind, start: start, len: len, line: 1, col: start + 1, file: 0 });
    return syntax::tokens.len - 1;
}

fn test_keywords() {
    expect_int("keyword let", syntax::keyword_kind("let"), syntax::KW_LET);
    expect_int("keyword within", syntax::keyword_kind("within"), syntax::KW_WITHIN);
    expect_int("keyword adapt", syntax::keyword_kind("adapt"), syntax::KW_ADAPT);
    expect_int("not a keyword", syntax::keyword_kind("Length"), syntax::TK_IDENT);
    expect_int("keyword count", syntax::KEYWORDS.len, syntax::KW_WITHIN - syntax::KW_ADAPT + 1);
}

fn test_tree() {
    syntax::reset();
    // let L: Length = 200 mm
    let f = syntax::add_file("t.skn", "let L: Length = 200 mm");
    expect_int("file index", f, 0);
    let t_let = tok(syntax::KW_LET, 0, 3);
    let t_l = tok(syntax::TK_IDENT, 4, 1);
    let t_len = tok(syntax::TK_IDENT, 7, 6);
    let t_num = tok(syntax::TK_NUMBER, 16, 3);
    let t_mm = tok(syntax::TK_IDENT, 20, 2);
    expect_eq("token text", syntax::text(t_len), "Length");

    let n = syntax::node(syntax::N_LET, t_l);
    syntax::add(n, syntax::node(syntax::N_TYPE_NAME, t_len));
    let num = syntax::node(syntax::N_NUMBER, t_num);
    let unit = syntax::node(syntax::N_UNIT, t_mm);
    syntax::add(unit, syntax::node(syntax::N_UNIT_NAME, t_mm));
    syntax::add(num, unit);
    syntax::add(n, num);

    expect_int("child count", syntax::child_count(n), 2);
    expect_int("child 1", syntax::child(n, 1), num);
    expect_int("no child 2", syntax::child(n, 2), -1);
    expect_int("child of kind", syntax::child_of_kind(num, syntax::N_UNIT), unit);
    expect_int("no child of kind", syntax::child_of_kind(n, syntax::N_UNIT), -1);
    expect_int("unused token", t_let, 0);

    var out: u8[>..] = [];
    syntax::dump(n, out);
    expect_eq("dump", out, "(let \"L\" (type \"Length\") (number \"200\" (unit \"mm\" (unit_name \"mm\"))))");

    syntax::reset();
    expect_int("reset", syntax::nodes.len + syntax::tokens.len + syntax::source.len, 0);
}

fn test_kind_names() {
    expect_eq("kind name let", syntax::kind_name(syntax::N_LET), "let");
    expect_eq("kind name at", syntax::kind_name(syntax::N_AT), "at");
    expect_eq("kind name unit_pow", syntax::kind_name(syntax::N_UNIT_POW), "unit_pow");
    expect_eq("kind name block", syntax::kind_name(syntax::N_BLOCK), "block");
    expect_eq("kind name type_apply", syntax::kind_name(syntax::N_TYPE_APPLY), "type_apply");
    expect_eq("kind name number", syntax::kind_name(syntax::N_NUMBER), "number");
    expect_eq("kind name out of range", syntax::kind_name(1000), "?");
}

// Parses `src` as one file. It gives the dump of each item on its own line,
// then each error as `id line:col`.
fn parse_text(src: u8[:]) -> u8[>..] {
    syntax::reset();
    let f = syntax::add_file("t.skn", src);
    syntax::lex(f);
    let prog = syntax::parse(f);
    var out: u8[>..] = [];
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        syntax::dump(c, out);
        out.push('\n');
        c = syntax::nodes[c].next;
    }
    for e in syntax::errors {
        let t = syntax::tokens[e.tok];
        out.append(str(e.id, " ", t.line, ":", t.col, "\n"));
        if e.fix.len == 0 { out.append("no fix\n"); }
    }
    return out;
}

let HEADER = "model m: 3d, SI\n";

// `src` after the header gives the item `want`.
fn expect_item(src: u8[:], want: u8[:]) {
    let got = parse_text(str(HEADER, src));
    expect_eq(src, got, str("(header \"m\" (space \"3d\") (unitsys \"SI\"))\n", want, "\n"));
}

fn test_parse() {
    // A unit expression ends where a unit name cannot follow (spec §6.3).
    expect_item("let w: Real = 2 m/s*(1 - y/R)",
        "(let \"w\" (type \"Real\") (binary \"*\" (number \"2\" (unit \"m\" (unit_div \"/\" (unit_name \"m\") (unit_name \"s\")))) (binary \"-\" (number \"1\") (binary \"/\" (name \"y\") (name \"R\")))))");
    expect_item("let k: Real = 1.2 W/(m*K)",
        "(let \"k\" (type \"Real\") (number \"1.2\" (unit \"W\" (unit_div \"/\" (unit_name \"W\") (unit_mul \"*\" (unit_name \"m\") (unit_name \"K\"))))))");
    expect_item("let k: Force/Length^2 = 1 N/mm^2",
        "(let \"k\" (dim_div \"/\" (type \"Force\") (dim_pow \"^\" 2 (type \"Length\"))) (number \"1\" (unit \"N\" (unit_div \"/\" (unit_name \"N\") (unit_pow \"^\" 2 (unit_name \"mm\"))))))");
    // An exponent of a type or a unit is a multiple of 0.5 (femahoy gap G14).
    expect_item("let k: Pressure*Length^-0.5 = 1 MPa*mm^1.5",
        "(let \"k\" (dim_mul \"*\" (type \"Pressure\") (dim_pow \"^\" -0.5 (type \"Length\"))) (number \"1\" (unit \"MPa\" (unit_mul \"*\" (unit_name \"MPa\") (unit_pow \"^\" 1.5 (unit_name \"mm\"))))))");
    expect_item("let z: Vector<Complex<Velocity>> = 0",
        "(let \"z\" (type_apply \"Vector\" (type_apply \"Complex\" (type \"Velocity\"))) (number \"0\"))");
    expect_item("let ok: Bool = !a && b || c != d",
        "(let \"ok\" (type \"Bool\") (binary \"||\" (binary \"&&\" (unary \"!\" (name \"a\")) (name \"b\")) (binary \"!=\" (name \"c\") (name \"d\"))))");
    // A line that ends with an operator or an open bracket continues (spec §2.2).
    expect_item("let a: Real = b +\n    c // note\n",
        "(let \"a\" (type \"Real\") (binary \"+\" (name \"b\") (name \"c\")))");
    expect_item("let a: Real = f(\n    b,\n    c,\n)",
        "(let \"a\" (type \"Real\") (call \"(\" (name \"f\") (arg \"b\" (name \"b\")) (arg \"c\" (name \"c\"))))");
    // A keyword can be the name of an argument.
    expect_item("let r: Real = results(\"a.skn\", mesh: m)",
        "(let \"r\" (type \"Real\") (call \"(\" (name \"results\") (arg \"\"a.skn\"\" (string \"\"a.skn\"\")) (arg \"mesh\" aux=\"mesh\" (name \"m\"))))");
    expect_item("table t { x: Length in mm, y: Real }:\n    1, 2\n    3, 4\nlet b: Real = 1",
        "(table \"t\" (column \"x\" (type \"Length\") (unit \"mm\" (unit_name \"mm\"))) (column \"y\" (type \"Real\")) (row \"1\" (number \"1\") (number \"2\")) (row \"3\" (number \"3\") (number \"4\")))\n(let \"b\" (type \"Real\") (number \"1\"))");
    expect_item("study s: sample(n: 4) {\n    D: [1, 2] mm\n    H: 1 m..2 m\n}",
        "(study \"s\" (proc \"sample\" (arg \"n\" aux=\"n\" (number \"4\"))) (block \"{\" (setting \"D\" (list \"[\" (number \"1\") (number \"2\") (unit \"mm\" (unit_name \"mm\")))) (setting \"H\" (range \"..\" (number \"1\" (unit \"m\" (unit_name \"m\"))) (number \"2\" (unit \"m\" (unit_name \"m\")))))))");
    // A comment line in a table is no row.
    expect_item("table t { x: Real }:\n    1\n    // Ghia et al.\n// note\n    2\n",
        "(table \"t\" (column \"x\" (type \"Real\")) (row \"1\" (number \"1\")) (row \"2\" (number \"2\")))");
    // A procedure argument can be a list (femahoy gap G9).
    expect_item("step f: harmonic(freq: [50] Hz) {\n    output u\n}",
        "(step \"f\" (proc \"harmonic\" (arg \"freq\" aux=\"freq\" (list \"[\" (number \"50\") (unit \"Hz\" (unit_name \"Hz\"))))) (output \"output\" (name \"u\")))");
    // A positional procedure argument (spec §22, `minimize`).
    expect_item("study o: minimize(drag, method: bobyqa) {\n    D: 1 m..2 m\n}",
        "(study \"o\" (proc \"minimize\" (arg \"drag\" (name \"drag\")) (arg \"method\" aux=\"method\" (name \"bobyqa\"))) (block \"{\" (setting \"D\" (range \"..\" (number \"1\" (unit \"m\" (unit_name \"m\"))) (number \"2\" (unit \"m\" (unit_name \"m\")))))))");
    // `-` before a `degC` literal is part of the literal (femahoy gap G29).
    expect_item("let T: Temperature = -5 degC",
        "(let \"T\" (type \"Temperature\") (number \"5\" aux=\"-\" (unit \"degC\" (unit_name \"degC\"))))");
    expect_item("let T: Temperature = -(5 degC)",
        "(let \"T\" (type \"Temperature\") (unary \"-\" (number \"5\" (unit \"degC\" (unit_name \"degC\")))))");
    expect_item("let T: Temperature = - 5 degC",
        "(let \"T\" (type \"Temperature\") (unary \"-\" (number \"5\" (unit \"degC\" (unit_name \"degC\")))))");
    expect_item("let T: Temperature = -5 K",
        "(let \"T\" (type \"Temperature\") (unary \"-\" (number \"5\" (unit \"K\" (unit_name \"K\")))))");
    // An exponent can be any expression (femahoy gap G30).
    expect_item("let a: Real = j^(-al/3)",
        "(let \"a\" (type \"Real\") (binary \"^\" (name \"j\") (binary \"/\" (unary \"-\" (name \"al\")) (number \"3\"))))");
    expect_item("let a: Real = x^1.3",
        "(let \"a\" (type \"Real\") (binary \"^\" (name \"x\") (number \"1.3\")))");
    expect_item("expect e: u.x at (0 mm, a) == 3*S within 5%",
        "(expect \"e\" (at \"at\" (member \"x\" (name \"u\")) (vector \"(\" (number \"0\" (unit \"mm\" (unit_name \"mm\"))) (name \"a\"))) (binary \"*\" (number \"3\") (name \"S\")) (percent \"5%\"))");
}

// `src` gives the errors `want`, each as `id line:col`.
fn expect_errors(src: u8[:], want: u8[:]) {
    let got = parse_text(src);
    var errs: u8[>..] = [];
    for line in split(got, '\n') {
        if line.len > 0 && line[0] != '(' { errs.append(line); errs.push('\n'); }
    }
    expect_eq(src, errs, want);
}

fn test_errors() {
    expect_errors("model m: 3d, SI\nlet E Pressure = 1 Pa\nlet b: Real = 2\n", "syntax 2:7\n");
    expect_errors("let a: Real = 1\nmodel m: 3d, SI\n", "header_order 1:1\nheader_order 2:1\n");
    expect_errors("model m: 3d, SI\nmodel n: 3d, SI\n", "header_order 2:1\n");
    expect_errors("", "header_order 1:1\n");
    // The lexer reports these, and the parser adds no second error.
    expect_errors("model m: 3d, SI\nlet s: String = \"abc\n", "syntax 2:17\n");
    expect_errors("model m: 3d, SI\nlet a: Real = 1 # 2\n", "syntax 2:17\n");
    expect_errors("model m: 3d, SI\nlet a: Length^0.25 = 1 m\n", "syntax 2:20\n");
    expect_errors("model m: 3d, SI\nlet a: Bool = a < b < c\n", "syntax 2:21\n");
    // A list is only a study entry or a procedure argument.
    expect_errors("model m: 3d, SI\nlet a: Real = [1, 2]\n", "syntax 2:15\n");
    expect_errors("model m: 3d, SI\nlet a: Real = f([1, 2])\n", "syntax 2:17\n");
    expect_errors("model m: 3d, SI\nmesh m on x { size: [1, 2] mm }\n", "syntax 2:21\n");
    expect_errors("model m: 3d, SI\nmesh m on x { size: 1 mm..2 mm }\n", "syntax 2:25\n");
    // An open bracket does not hide the next item.
    expect_errors("model m: 3d, SI\nlet a: Real = (1, 2\nlet b: Real = 3\nlet c: Real =\n", "syntax 2:20\nsyntax 5:1\n");
    expect_errors("model m: 3d, SI\nmaterial s {\n    E: Pressure = 1 Pa\nlet b: Real = c d e\n", "syntax 4:1\nsyntax 4:17\n");
    expect_errors("model m: 3d, SI\nmesh m on x {\n    size: 1 mm\nlet b: Real = c d e\n", "syntax 4:1\nsyntax 4:17\n");
    expect_errors("model m: 3d, SI\nlet r: R = {\n    a: 1\nlet b: Real = c d e\n", "syntax 4:1\nsyntax 4:17\n");
    expect_errors("model m: 3d, SI\nphysics p solves u: Length {\n    eq div(u = 0\n}\nlet z: Real = 1\n", "syntax 3:14\n");
}

// Each `.skn` file in `dir` parses with no error.
fn test_corpus(dir: u8[:]) {
    var names: u8[>..] = [];
    guard list_dir(dir, names) else { return; }
    for name in split(names, '\n') {
        guard name.len > 4 && name[name.len - 4..] == ".skn" else { continue; }
        let path = str(dir, "/", name);
        var src: u8[>..] = [];
        guard read_file(path, src) else {
            failures++;
            print("FAIL ", path, ": cannot read the file");
            continue;
        }
        syntax::reset();
        let f = syntax::add_file(path, src);
        syntax::lex(f);
        syntax::parse(f);
        for e in syntax::errors {
            let t = syntax::tokens[e.tok];
            failures++;
            print("FAIL ", path, ":", t.line, ":", t.col, ": ", e.id, ": ", e.msg);
        }
    }
}

fn main() {
    test_keywords();
    test_tree();
    test_kind_names();
    test_parse();
    test_errors();
    test_corpus("samples");
    test_corpus("examples");
    test_corpus("examples/femahoy");
    if failures > 0 {
        print(failures, " failed");
        exit(1);
    }
    print("syntax: ok");
}
