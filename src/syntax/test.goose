// Tests of the syntax tree tables. `./build.sh test` runs this program, and
// a non-zero exit status is a failure.
import std;
import os;
import .tree;

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

fn main() {
    test_keywords();
    test_tree();
    test_kind_names();
    if failures > 0 {
        print(failures, " failed");
        exit(1);
    }
    print("syntax: ok");
}
