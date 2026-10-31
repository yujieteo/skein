// Tests of `skein fmt`. `./build.sh test` runs this program from the root of
// the repository, and a non-zero exit status is a failure.
import std;
import os;
import syntax from "../syntax";
import .fmt;

var failures = 0;

fn fail(what: u8[:], got: u8[:], want: u8[:]) {
    failures++;
    print("FAIL ", what, "\n--- got:\n", got, "--- want:\n", want);
}

// Formats `text` and gives the result. A result that is not OK is a failure.
fn formatted(what: u8[:], text: u8[:]) -> u8[>..] {
    var result: u8[>..] = [];
    let r = fmt::format_text("t.skn", text, result);
    if r != fmt::OK {
        failures++;
        print("FAIL ", what, ": format_text gives ", r);
    }
    return result;
}

// `fmt` changes `text` to `want`, and changes `want` not at all.
fn expect_fmt(what: u8[:], text: u8[:], want: u8[:]) {
    let once = formatted(what, text);
    if once != want { fail(what, once, want); return; }
    let twice = formatted(what, once);
    if twice != once { fail(str(what, ", again"), twice, once); }
}

let HEAD = "model t: 3d, SI\n\n";

fn test_spaces() {
    expect_fmt("spaces",
        str(HEAD, "let  a:Real=2 * 3+ -x/y ^2\nlet b: Vector<Length> = ( 0,0, 1 )mm\n"),
        str(HEAD, "let a: Real = 2*3 + -x/y^2\nlet b: Vector<Length> = (0, 0, 1) mm\n"));
    expect_fmt("units and signs",
        str(HEAD, "let k: Real = 3 kg/m^3*f(a)-  4 mm^-1\nlet t: Temperature = -5 degC\n"),
        str(HEAD, "let k: Real = 3 kg/m^3*f(a) - 4 mm^-1\nlet t: Temperature = -5 degC\n"));
    expect_fmt("range and block",
        str(HEAD, "param Re: Real = 20 in 1 .. 200\n\nmethod m: cg {mesh: m,order: {v: 1}}\n"),
        str(HEAD, "param Re: Real = 20 in 1..200\n\nmethod m: cg { mesh: m, order: { v: 1 } }\n"));
    expect_fmt("table",
        str(HEAD, "table g {y: Length in m, u: Real}:\n  1.0,2\n\t0.5 ,  -3\n\nlet a: Real = 1\n"),
        str(HEAD, "table g{y: Length in m, u: Real}:\n    1.0, 2\n    0.5, -3\n\nlet a: Real = 1\n"));
    expect_fmt("comment between table rows",
        str(HEAD, "table g{y: Length in m, u: Real}:\n 1.0, 2\n // mid\n// low\n 0.5, -3\n// after\nlet a: Real = 1\n"),
        str(HEAD, "table g{y: Length in m, u: Real}:\n    1.0, 2\n    // mid\n    // low\n    0.5, -3\n\n// after\nlet a: Real = 1\n"));
}

fn test_blank_lines() {
    expect_fmt("let runs keep one blank line",
        str(HEAD, "let a: Real = 1\n\n\n\nlet b: Real = 2\nlet c: Real = 3\n"),
        str(HEAD, "let a: Real = 1\n\nlet b: Real = 2\nlet c: Real = 3\n"));
    expect_fmt("one blank line between other items",
        str(HEAD, "let a: Real = 1\nexpect e: a == 1 within 1%\nexpect f: a == 1 within 1%\n"),
        str(HEAD, "let a: Real = 1\n\nexpect e: a == 1 within 1%\n\nexpect f: a == 1 within 1%\n"));
    expect_fmt("bodies",
        str(HEAD, "step s: static {\n\n  bc a: u = 0\n\n  bc b: u = 0\n\n\n  output u\n\n}\n"),
        str(HEAD, "step s: static {\n    bc a: u = 0\n    bc b: u = 0\n\n    output u\n}\n"));
    expect_fmt("the header and its comment",
        "// About t.\n\n\nmodel t: 3d, SI\nlet a: Real = 1\n",
        "// About t.\n\nmodel t: 3d, SI\n\nlet a: Real = 1\n");
}

fn test_comments() {
    expect_fmt("trailing comments",
        str(HEAD, "let a: Real = 1 // one\nmaterial s {   // steel\n    E: Pressure = 1 GPa  //E\n}\n"),
        str(HEAD, "let a: Real = 1  // one\n\nmaterial s {  // steel\n    E: Pressure = 1 GPa  //E\n}\n"));
    expect_fmt("a comment inside a statement keeps the lines",
        str(HEAD, "let a: Real = f(1,  // one\n  2)\n"),
        str(HEAD, "let a: Real = f(1,  // one\n    2)\n"));
    expect_fmt("a comment before a statement in a body",
        str(HEAD, "material s {\n// the modulus\n    E: Pressure = 1 GPa\n\n  // end\n}\n"),
        str(HEAD, "material s {\n    // the modulus\n    E: Pressure = 1 GPa\n\n    // end\n}\n"));
}

fn test_long_lines() {
    // A statement on many lines that fits joins into one line.
    expect_fmt("join",
        str(HEAD, "let a: Real = f(\n    1,\n    2,\n)\nlet b: Real = 1 +\n    2\n"),
        str(HEAD, "let a: Real = f(1, 2)\nlet b: Real = 1 + 2\n"));
    let terms = "alpha_one*beta_one + alpha_two*beta_two - alpha_three*beta_three + alpha_four*beta_four*gamma";
    expect_fmt("break after an operator",
        str(HEAD, "let a: Real = ", terms, "\n"),
        str(HEAD, "let a: Real = alpha_one*beta_one +\n    alpha_two*beta_two -\n",
            "    alpha_three*beta_three +\n    alpha_four*beta_four*gamma\n"));
    let items = "first_argument_value, second_argument_value, third_argument_value, fourth_argument_value";
    expect_fmt("one item on each line",
        str(HEAD, "let a: Real = f(", items, ")\n"),
        str(HEAD, "let a: Real = f(\n    first_argument_value,\n    second_argument_value,\n",
            "    third_argument_value,\n    fourth_argument_value,\n)\n"));
    expect_fmt("the sides of an eq",
        str(HEAD, "physics p solves u: Real {\n    eq alpha_one*beta_one + alpha_two*beta_two - alpha_three*beta_three",
            " = alpha_four*beta_four*gamma_four\n}\n"),
        str(HEAD, "physics p solves u: Real {\n    eq alpha_one*beta_one + alpha_two*beta_two - alpha_three*beta_three",
            " =\n        alpha_four*beta_four*gamma_four\n}\n"));
    expect_fmt("a long block goes on many lines",
        str(HEAD, "step s: static {\n    solve { nonlinear: newton, linear: gmres, precond: ilu0, ",
            "rtol: 1e-8, atol: 1e-10, max_iterations: 50 }\n}\n"),
        str(HEAD, "step s: static {\n    solve {\n        nonlinear: newton\n        linear: gmres\n",
            "        precond: ilu0\n        rtol: 1e-8\n        atol: 1e-10\n",
            "        max_iterations: 50\n    }\n}\n"));
    // A line breaks at an operator inside a list that ends with an added comma.
    let calls = "g(alpha_one, alpha_one) + g(alpha_one, alpha_one) - g(alpha_one, alpha_one) - g(alpha_one, alpha_one)";
    expect_fmt("an operator after an added comma",
        str(HEAD, "expect e: f(", calls, ") == 1 within 1\n"),
        str(HEAD, "expect e: f(\n    g(alpha_one, alpha_one) +\n        g(alpha_one, alpha_one) -\n",
            "        g(alpha_one, alpha_one) -\n        g(alpha_one, alpha_one),\n) == 1 within 1\n"));
}

fn test_errors() {
    var result: u8[>..] = [];
    let r = fmt::format_text("t.skn", str(HEAD, "let a: = 1\n"), result);
    if r != fmt::SYNTAX_ERROR { failures++; print("FAIL syntax error: format_text gives ", r); }
}

// The names of the files in directory `dir` that end with `.skn`.
fn models(dir: u8[:], into: u8[>..]&) -> const u8[:][>..] {
    var names: const u8[:][>..] = [];
    if !list_dir(dir, into) { return names; }
    for name in split(into, '\n') {
        if ends_with(name, ".skn") { names.push(name); }
    }
    return names;
}

// Each sample is in the standard layout, and a change of one value changes
// only its line (implementation.md §10).
fn test_samples() {
    var list: u8[>..] = [];
    let names = models("samples", list);
    if names.len < 6 { failures++; print("FAIL samples: found ", names.len); }
    for name in names {
        let path = str("samples/", name);
        var text: u8[>..] = [];
        read_file(path, text);
        let once = formatted(path, text);
        if once != text { fail(path, once, text); continue; }
        // The first value that is a number.
        var at = -1;
        for i in 0..text.len - 4 {
            if text[i..i + 3] == " = " && text[i + 3] >= '0' && text[i + 3] <= '9' { at = i; break; }
        }
        if at < 0 { failures++; print("FAIL ", path, ": no number to change"); continue; }
        var edited: u8[>..] = [];
        edited.append(text[..at + 3]);
        edited.append("9");
        edited.append(text[at + 3..]);
        let after = formatted(str(path, ", edited"), edited);
        if after != edited { fail(str(path, ", edited"), after, edited); }
    }
}

// `fmt` is idempotent on each femahoy model that parses.
fn test_femahoy() {
    var list: u8[>..] = [];
    let names = models("examples/femahoy", list);
    if names.len == 0 { failures++; print("FAIL femahoy: no models"); }
    for name in names {
        let path = str("examples/femahoy/", name);
        var text: u8[>..] = [];
        read_file(path, text);
        var once: u8[>..] = [];
        let r = fmt::format_text(path, text, once);
        if r == fmt::SYNTAX_ERROR { continue; }
        if r != fmt::OK { failures++; print("FAIL ", path, ": format_text gives ", r); continue; }
        var twice: u8[>..] = [];
        fmt::format_text(path, once, twice);
        if twice != once { fail(str(path, ", again"), twice, once); }
    }
}

// `skein fmt --check` and `skein fmt` on a file.
fn test_command() {
    let path = "build/fmt_test.skn";
    let messy = str(HEAD, "let a:Real=1\n");
    write_file(path, messy);
    let check: const u8[:][2] = ["--check", path];
    let write: const u8[:][1] = [path];
    if fmt::run(check) != 1 { failures++; print("FAIL --check on a messy file"); }
    var text: u8[>..] = [];
    read_file(path, text);
    if text != messy { fail("--check changes no file", text, messy); }
    if fmt::run(write) != 0 { failures++; print("FAIL fmt on a file"); }
    text.clear();
    read_file(path, text);
    let want = str(HEAD, "let a: Real = 1\n");
    if text != want { fail("fmt rewrites the file", text, want); }
    if fmt::run(check) != 0 { failures++; print("FAIL --check on a formatted file"); }
    delete_file(path);
}

fn main() {
    test_spaces();
    test_blank_lines();
    test_comments();
    test_long_lines();
    test_errors();
    test_samples();
    test_femahoy();
    test_command();
    if failures > 0 {
        print("fmt: ", failures, " failures");
        exit(1);
    }
    print("fmt: ok");
}
