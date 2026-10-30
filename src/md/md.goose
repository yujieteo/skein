// `skein md` (docs/command_line.md §11). It writes a model as one Markdown
// document that a person reads: the header, then each item in the order of
// the file, with the line of each item, so that a reader can find each part
// in the source.
//
// The math is TeX inside `$...$` and `$$...$$`, as GitHub and the common
// Markdown readers render it. The document holds the model only, never
// results. It holds no date, so the same model always gives the same text.
import std;
import os;
namespace md;
import syntax from "../syntax";
import check from "../check";

// The document, as `write_model` makes it.
var out: u8[>..];

fn node_kind(n: i64) -> i64 { syntax::nodes[n].kind }
fn node_text(n: i64) -> const u8[:] { syntax::text(syntax::nodes[n].tok) }
fn node_line(n: i64) -> i64 { syntax::tokens[syntax::nodes[n].tok].line }
fn first(n: i64) -> i64 { syntax::nodes[n].first }
fn next(n: i64) -> i64 { syntax::nodes[n].next }

// --- text -----------------------------------------------------------------------------

// Appends `s` as Markdown text: each punctuation character that Markdown
// gives a meaning to has a backslash.
fn text(s: u8[:]) {
    for c in s {
        if contains("\\`*_[]<>#|$~", c) { out.push('\\'); }
        out.push(c);
    }
}

// The math holds no backslash before a punctuation character, such as `\,`
// or `\_`. Reason: GitHub applies the Markdown escapes inside `$...$`, so
// `\_` arrives as `_`, and the other readers do not.

// Appends the identifier `s` inside math, with each `_` as a comma, as in a
// subscript `cell,x`.
fn ident(s: u8[:]) {
    for c in s { out.push(if c == '_' { ',' } else { c }); }
}

// Appends `s` inside \text. The characters that TeX gives a meaning to
// become spaces, and `_` is a line under a space, outside the \text.
fn plain(s: u8[:]) {
    for c in s {
        if c == '_' {
            out.append("}\\underline{\\ }\\text{");
        } else {
            out.push(if contains("\\{}$%#&~^", c) { ' ' } else { c });
        }
    }
}

// Appends `s` as inline code. Skein names and types hold no backquote.
fn code(s: u8[:]) {
    out.push('`');
    out.append(s);
    out.push('`');
}

fn code_node(n: i64) { code(node_text(n)); }

// The source line of `n`, in parentheses.
fn line_note(n: i64) { format(out, " (line ", node_line(n), ")"); }

// --- types ----------------------------------------------------------------------------

// Appends the type `n` in Skein notation, for example `Vector<Length>`.
recursive fn type_str(n: i64, s: u8[>..]&) {
    match node_kind(n) {
        syntax::N_TYPE_APPLY => {
            s.append(node_text(n));
            s.push('<');
            type_str(first(n), s);
            s.push('>');
        },
        syntax::N_DIM_MUL, syntax::N_DIM_DIV => {
            type_str(first(n), s);
            s.append(node_text(n));
            // The products are left to right, so a product on the right
            // had parentheses.
            let r = next(first(n));
            let wrap = node_kind(r) == syntax::N_DIM_MUL || node_kind(r) == syntax::N_DIM_DIV;
            if wrap { s.push('('); }
            type_str(r, s);
            if wrap { s.push(')'); }
        },
        syntax::N_DIM_POW => {
            type_str(first(n), s);
            format(s, "^", syntax::nodes[n].aux);
        },
        _ => { s.append(node_text(n)); },
    }
}

// The type `n` as inline code.
fn type_code(n: i64) {
    var s: u8[>..] = [];
    type_str(n, s);
    code(s);
}

// True for the type of a geometry value, which the document gives in words.
fn is_geometry_type(n: i64) -> bool {
    if node_kind(n) != syntax::N_TYPE_NAME { return false; }
    let s = node_text(n);
    return s == "Solid" || s == "Surface" || s == "Curve" || s == "Point";
}

// --- math -----------------------------------------------------------------------------

let GREEK = [
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
    "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi",
    "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Upsilon", "Phi", "Psi",
    "Omega",
];

fn all_digits(s: u8[:]) -> bool {
    for c in s {
        if c < '0' || c > '9' { return false; }
    }
    return s.len > 0;
}

// A name without `_`: a Greek letter, one letter, a letter and digits, or an
// italic word.
fn math_word(s: u8[:]) {
    if contains(GREEK, s) {
        out.push('\\');
        out.append(s);
    } else if s.len == 1 {
        out.append(s);
    } else if all_digits(s[1..]) {
        out.push(s[0]);
        format(out, "_{", s[1..], "}");
    } else {
        out.append("\\mathit{");
        ident(s);
        out.push('}');
    }
}

// A name: `h_air` is h with the subscript "air", and `sigma` is σ.
fn math_name(s: u8[:]) {
    let i = find(s, "_");
    if i <= 0 || i == s.len - 1 {
        math_word(s);
        return;
    }
    math_word(s[..i]);
    out.append("_{\\mathrm{");
    ident(s[i + 1..]);
    out.append("}}");
}

// A number. `1e7` is 1 × 10^7.
fn math_number(s: u8[:]) {
    var e = find(s, "e");
    if e < 0 { e = find(s, "E"); }
    if e < 0 {
        out.append(s);
        return;
    }
    let exp = if e + 1 < s.len && s[e + 1] == '+' { s[e + 2..] } else { s[e + 1..] };
    if s[..e] != "1" { format(out, s[..e], "\\times "); }
    format(out, "10^{", exp, "}");
}

// The units that take an SI prefix (syntax/units.goose).
fn micro(s: u8[:]) -> bool {
    s.len > 1 && s[0] == 'u' && contains(syntax::PREFIXED, s[1..])
}

recursive fn math_unit(n: i64) {
    match node_kind(n) {
        syntax::N_UNIT => { math_unit(first(n)); },
        syntax::N_UNIT_NAME => {
            let s = node_text(n);
            if s == "degC" {
                out.append("{}^{\\circ}\\mathrm{C}");
            } else if s == "deg" {
                out.append("{}^{\\circ}");
            } else if micro(s) {
                format(out, "\\mu\\mathrm{", s[1..], "}");
            } else {
                format(out, "\\mathrm{", s, "}");
            }
        },
        syntax::N_UNIT_MUL, syntax::N_UNIT_DIV => {
            math_unit(first(n));
            out.append(if node_kind(n) == syntax::N_UNIT_MUL { "\\cdot " } else { "/" });
            let r = next(first(n));
            let wrap = node_kind(r) == syntax::N_UNIT_MUL || node_kind(r) == syntax::N_UNIT_DIV;
            if wrap { out.push('('); }
            math_unit(r);
            if wrap { out.push(')'); }
        },
        syntax::N_UNIT_POW => {
            let b = first(n);
            let wrap = node_kind(b) != syntax::N_UNIT_NAME;
            if wrap { out.push('('); }
            math_unit(b);
            if wrap { out.push(')'); }
            format(out, "^{", syntax::nodes[n].aux, "}");
        },
        _ => {},
    }
}

// The unit of a literal, after a thin space, if it has one.
fn unit_after(n: i64) {
    let u = syntax::child_of_kind(n, syntax::N_UNIT);
    if u < 0 { return; }
    out.append("\\ ");
    math_unit(u);
}

// The precedence of an operator of kind `k` (spec §7.1), times 10.
fn op_prec(k: i64) -> i64 {
    match k {
        syntax::TK_OROR, syntax::TK_BAR => 20,
        syntax::TK_ANDAND, syntax::TK_AMP => 30,
        syntax::TK_EQ, syntax::TK_NE, syntax::TK_LT, syntax::TK_LE, syntax::TK_GT, syntax::TK_GE => 40,
        syntax::TK_PLUS, syntax::TK_MINUS => 50,
        syntax::TK_STAR, syntax::TK_SLASH => 60,
        syntax::TK_CARET => 70,
        _ => 0,
    }
}

// How tightly the expression `n` binds. A `-x` binds between `+` and `*`.
fn prec(n: i64) -> i64 {
    match node_kind(n) {
        syntax::N_BINARY => op_prec(syntax::tokens[syntax::nodes[n].tok].kind),
        syntax::N_UNARY => 55,
        _ => 90,
    }
}

// `n` in parentheses when it binds less tightly than `least`.
recursive fn math_at(n: i64, least: i64) {
    if prec(n) >= least {
        math(n);
        return;
    }
    out.append("(");
    math(n);
    out.append(")");
}

let OPS = ["||", "\\lor", "|", "\\cup", "&&", "\\land", "&", "\\cap", "==", "=", "!=", "\\neq",
           "<=", "\\leq", ">=", "\\geq", "<", "<", ">", ">", "+", "+", "-", "-"];

recursive fn math_binary(n: i64) {
    let l = first(n);
    let r = next(l);
    let op = node_text(n);
    let p = prec(n);
    if op == "/" {
        out.append("\\frac{");
        math(l);
        out.append("}{");
        math(r);
        out.push('}');
        return;
    }
    if op == "^" {
        out.push('{');
        math_at(l, 90);
        out.append("}^{");
        math(r);
        out.push('}');
        return;
    }
    if op == "*" {
        math_at(l, 55);
        let number = node_kind(r) == syntax::N_NUMBER || node_kind(r) == syntax::N_PERCENT;
        out.append(if number { " \\cdot " } else { " " });
        math_at(r, 60);
        return;
    }
    // A comparison does not chain, and `a - (b - c)` keeps its parentheses.
    math_at(l, if p == 40 { 41 } else { p });
    for o, i in OPS {
        if i % 2 == 0 && o == op { format(out, " ", OPS[i + 1], " "); }
    }
    math_at(r, p + 1);
}

// The arguments of call `n`, from `a`, separated by commas.
recursive fn math_args(a: i64) {
    var c = a;
    var count = 0;
    while c >= 0 {
        if count > 0 { out.append(", "); }
        math(c);
        count++;
        c = next(c);
    }
}

let FUNCTIONS = ["exp", "sin", "cos", "tan", "tanh", "ln", "max", "min", "det"];

recursive fn math_call(n: i64) {
    let callee = first(n);
    let args = next(callee);
    // The value of the first argument.
    let a = if args >= 0 { first(args) } else { -1 };
    let name = if node_kind(callee) == syntax::N_NAME { node_text(callee) } else { "" };
    let one = args >= 0 && next(args) < 0;
    let two = args >= 0 && next(args) >= 0 && next(next(args)) < 0;
    if one && (name == "div" || name == "grad" || name == "curl" || name == "lap") {
        out.append(if name == "div" { "\\nabla\\cdot " }
                   else if name == "grad" { "\\nabla " }
                   else if name == "curl" { "\\nabla\\times " }
                   else { "\\nabla^{2} " });
        math_at(a, 90);
        return;
    }
    if one && name == "ddt" {
        if node_kind(a) == syntax::N_NAME {
            out.append("\\frac{\\partial ");
            math(a);
            out.append("}{\\partial t}");
        } else {
            out.append("\\frac{\\partial}{\\partial t}(");
            math(a);
            out.append(")");
        }
        return;
    }
    if two && (name == "dot" || name == "cross") {
        math_at(a, 90);
        out.append(if name == "dot" { " \\cdot " } else { " \\times " });
        math_at(first(next(args)), 90);
        return;
    }
    if one && name == "sqrt" {
        out.append("\\sqrt{");
        math(a);
        out.push('}');
        return;
    }
    if one && (name == "abs" || name == "norm") {
        out.append(if name == "abs" { "\\lvert " } else { "\\lVert " });
        math(a);
        out.append(if name == "abs" { "\\rvert " } else { "\\rVert " });
        return;
    }
    if contains(FUNCTIONS, name) {
        format(out, "\\", name);
    } else if name.len > 0 {
        // `von_mises` is von with the subscript "mises", as for a name.
        let i = find(name, "_");
        out.append("\\operatorname{");
        ident(if i > 0 { name[..i] } else { name });
        out.push('}');
        if i > 0 {
            out.append("_{\\mathrm{");
            ident(name[i + 1..]);
            out.append("}}");
        }
    } else {
        math(callee);
    }
    out.append("(");
    math_args(args);
    out.append(")");
}

// Appends the expression `n` in TeX math. Parentheses come from the
// precedence, as the syntax tree has none.
recursive fn math(n: i64) {
    let k = syntax::nodes[n];
    match k.kind {
        syntax::N_NUMBER => {
            if k.aux >= 0 { out.push('-'); }
            math_number(node_text(n));
            unit_after(n);
        },
        syntax::N_PERCENT => {
            let s = node_text(n);
            // `%` starts a comment in TeX, so a percent inside an
            // expression is a word. `inline` writes a whole percent value
            // with `%` after the math.
            math_number(s[..s.len - 1]);
            out.append("\\ \\text{percent}");
        },
        syntax::N_STRING => {
            out.append("\\text{");
            plain(node_text(n));
            out.append("}");
        },
        syntax::N_BOOL => { format(out, "\\mathrm{", node_text(n), "}"); },
        syntax::N_NAME => { math_name(node_text(n)); },
        syntax::N_VECTOR => {
            out.append("(");
            var c = k.first;
            var count = 0;
            while c >= 0 {
                if node_kind(c) != syntax::N_UNIT {
                    if count > 0 { out.append(", "); }
                    math(c);
                    count++;
                }
                c = next(c);
            }
            out.append(")");
            unit_after(n);
        },
        syntax::N_LIST => {
            out.append("[");
            var c = k.first;
            var count = 0;
            while c >= 0 {
                if node_kind(c) != syntax::N_UNIT {
                    if count > 0 { out.append(", "); }
                    math(c);
                    count++;
                }
                c = next(c);
            }
            out.append("]");
            unit_after(n);
        },
        syntax::N_RECORD => {
            out.append("\\lbrace ");
            math_args(k.first);
            out.append("\\rbrace ");
        },
        syntax::N_RANGE => {
            math(k.first);
            out.append(" \\ldots ");
            math(next(k.first));
        },
        syntax::N_BINARY => { math_binary(n); },
        syntax::N_UNARY => {
            out.append(if node_text(n) == "!" { "\\neg " } else { "-" });
            math_at(k.first, 60);
        },
        syntax::N_CALL => { math_call(n); },
        syntax::N_ARG => {
            if k.aux >= 0 {
                math_name(syntax::text(k.aux));
                out.append(" = ");
            }
            math(k.first);
        },
        syntax::N_MEMBER => {
            out.push('{');
            math_at(k.first, 90);
            out.append("}_{");
            let m = node_text(n);
            var axes = true;
            for c in m { axes = axes && contains("xyz", c); }
            if axes { out.append(m); } else { format(out, "\\mathrm{", m, "}"); }
            out.push('}');
        },
        syntax::N_AT => {
            math_at(k.first, 90);
            out.append("\\ \\text{at}\\ ");
            math(next(k.first));
        },
        _ => {
            out.append("\\text{");
            plain(node_text(n));
            out.push('}');
        },
    }
}

// The expression `n` in inline math.
fn inline(n: i64) {
    if node_kind(n) == syntax::N_PERCENT {
        let s = node_text(n);
        out.push('$');
        math_number(s[..s.len - 1]);
        out.append("$ %");
        return;
    }
    out.push('$');
    math(n);
    out.push('$');
}

// The children of `n` from `c` in inline math, separated by commas.
fn inline_list(c0: i64) {
    var c = c0;
    var count = 0;
    while c >= 0 {
        if count > 0 { out.append(", "); }
        inline(c);
        count++;
        c = next(c);
    }
}

// `left = right` as a displayed equation.
fn equation(left: i64, right: i64) {
    out.append("$$\n");
    math(left);
    out.append(" = ");
    math(right);
    out.append("\n$$\n\n");
}

// --- geometry in words ------------------------------------------------------------------

// The value of argument `i` of call `n`, or -1.
fn arg(n: i64, i: i64) -> i64 {
    var c = next(first(n));
    var k = 0;
    while c >= 0 && k < i { c = next(c); k++; }
    if c < 0 { return -1; }
    return first(c);
}

fn arg_count(n: i64) -> i64 { syntax::child_count(n) - 1 }

// Appends the geometry expression `n` in words, with the values in math.
recursive fn words(n: i64) {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_BINARY {
        let op = node_text(n);
        let l = k.first;
        let r = next(l);
        if op == "-" || op == "|" || op == "&" {
            words(l);
            out.append(if op == "-" { " minus " } else if op == "|" { " joined to " } else { " intersected with " });
            words(r);
            return;
        }
    }
    if k.kind == syntax::N_CALL && node_kind(k.first) == syntax::N_NAME {
        let name = node_text(k.first);
        let count = arg_count(n);
        if (name == "box" || name == "rect") && count == 2 {
            out.append(if name == "box" { "the box from " } else { "the rectangle from " });
            inline(arg(n, 0));
            out.append(" to ");
            inline(arg(n, 1));
            return;
        }
        if name == "cylinder" && count == 4 {
            out.append("the cylinder with its base centre at ");
            inline(arg(n, 0));
            out.append(", axis ");
            inline(arg(n, 1));
            out.append(", radius ");
            inline(arg(n, 2));
            out.append(" and length ");
            inline(arg(n, 3));
            return;
        }
        if (name == "sphere" || name == "circle") && count == 2 {
            format(out, "the ", name, " with its centre at ");
            inline(arg(n, 0));
            out.append(" and radius ");
            inline(arg(n, 1));
            return;
        }
        if name == "boundary" && count == 2 {
            out.append("the part of the boundary of ");
            words(arg(n, 0));
            let where = arg(n, 1);
            if node_kind(where) == syntax::N_NAME {
                out.append(" that touches ");
                words(where);
            } else {
                out.append(" where ");
                inline(where);
            }
            return;
        }
        if name == "boundary" && count == 1 {
            out.append("the boundary of ");
            words(arg(n, 0));
            return;
        }
    }
    if k.kind == syntax::N_NAME {
        code_node(n);
        return;
    }
    inline(n);
}

// --- blocks -------------------------------------------------------------------------------

// The entries of block `b` as a list: `size: 5 mm`, `refine hole: 1 mm`.
fn block_list(b: i64) {
    if b < 0 || first(b) < 0 { return; }
    var c = first(b);
    while c >= 0 {
        out.append("- ");
        let k = node_kind(c);
        if k == syntax::N_SETTING {
            code_node(c);
            out.append(": ");
            setting_value(first(c));
        } else {
            out.append(if k == syntax::N_REFINE { "refine " } else { "layers " });
            words(first(c));
            out.append(": ");
            inline(next(first(c)));
        }
        out.push('\n');
        c = next(c);
    }
    out.push('\n');
}

// A setting value: a name is a word, as in `cell: tet`, and a block is a
// list of settings on one line.
recursive fn setting_value(v: i64) {
    if node_kind(v) == syntax::N_NAME {
        code_node(v);
    } else if node_kind(v) == syntax::N_RECORD {
        var a = first(v);
        var count = 0;
        while a >= 0 {
            if count > 0 { out.append(", "); }
            code(syntax::text(syntax::nodes[a].aux));
            out.append(": ");
            setting_value(first(a));
            count++;
            a = next(a);
        }
    } else {
        inline(v);
    }
}

// A procedure with its arguments, for example `static` or `transient` with
// $t_{\mathrm{end}} = 1\,\mathrm{s}$.
fn procedure(p: i64) {
    code_node(p);
    if first(p) < 0 { return; }
    out.append(" with ");
    inline_list(first(p));
}

// --- items --------------------------------------------------------------------------------

// The groups of items. A new section starts where the group changes.
let G_NONE = 0;
let G_VALUES = 1;
let G_GEOMETRY = 2;
let G_MATERIALS = 3;
let G_MESH = 4;
let G_PHYSICS = 5;
let G_BC = 6;
let G_STEPS = 7;
let G_CHECKS = 8;
let G_OTHER = 9;

let GROUP_TITLES = ["", "Parameters and constants", "Geometry and regions", "Materials",
                    "Discretization", "Physics", "Boundary conditions", "Steps", "Checks",
                    "Other items"];

fn group(n: i64) -> i64 {
    match node_kind(n) {
        syntax::N_PARAM, syntax::N_FN, syntax::N_TABLE => G_VALUES,
        syntax::N_LET => if is_geometry_type(first(n)) { G_GEOMETRY } else { G_VALUES },
        syntax::N_MATERIAL => G_MATERIALS,
        syntax::N_MESH, syntax::N_POINTS, syntax::N_METHOD => G_MESH,
        syntax::N_PHYSICS => G_PHYSICS,
        syntax::N_BC, syntax::N_CONTACT, syntax::N_INTERFACE, syntax::N_SYMMETRY => G_BC,
        syntax::N_STEP => G_STEPS,
        syntax::N_EXPECT, syntax::N_ALLOW, syntax::N_TOLERANCE => G_CHECKS,
        syntax::N_SCALES, syntax::N_HEADER => G_NONE,
        _ => G_OTHER,
    }
}

// True when `n` is a `param` or a `let` of a value: a row of the table of values.
fn value_row(n: i64) -> bool {
    let k = node_kind(n);
    return k == syntax::N_PARAM || (k == syntax::N_LET && !is_geometry_type(first(n)));
}

// The table of a run of `param` and `let` items from `n`. It gives the item
// after the run.
fn value_table(n0: i64) -> i64 {
    out.append("| Name | Kind | Type | Value | Line |\n|---|---|---|---|---|\n");
    var n = n0;
    while n >= 0 && value_row(n) {
        out.append("| ");
        code_node(n);
        out.append(if node_kind(n) == syntax::N_PARAM { " | param | " } else { " | let | " });
        type_code(first(n));
        out.append(" | ");
        let v = next(first(n));
        inline(v);
        var c = next(v);
        while c >= 0 {
            if node_kind(c) == syntax::N_TOL {
                out.append(", tolerance ");
                inline(first(c));
            } else if node_kind(c) == syntax::N_RANGE {
                out.append(", in ");
                inline(c);
            }
            c = next(c);
        }
        format(out, " | ", node_line(n), " |\n");
        n = next(n);
    }
    out.push('\n');
    return n;
}

// The list of a run of geometry `let` items from `n`. It gives the item
// after the run.
fn geometry_list(n0: i64) -> i64 {
    var n = n0;
    while n >= 0 && node_kind(n) == syntax::N_LET && is_geometry_type(first(n)) {
        out.append("- ");
        code_node(n);
        line_note(n);
        format(out, " is a ", node_text(first(n)), ": ");
        words(next(first(n)));
        out.append(".\n");
        n = next(n);
    }
    out.push('\n');
    return n;
}

fn material_item(n: i64) {
    out.append("### Material ");
    code_node(n);
    line_note(n);
    out.append("\n\n| Property | Type | Value | Line |\n|---|---|---|---|\n");
    var p = first(n);
    while p >= 0 {
        out.append("| ");
        code_node(p);
        out.append(" | ");
        type_code(first(p));
        out.append(" | ");
        inline(next(first(p)));
        format(out, " | ", node_line(p), " |\n");
        p = next(p);
    }
    out.push('\n');
}

fn discretization_item(n: i64) {
    let k = node_kind(n);
    if k == syntax::N_METHOD {
        out.append("### Method of ");
        code_node(n);
        line_note(n);
        out.append("\n\nThe method is ");
        code(syntax::text(syntax::nodes[n].aux));
        out.append(".\n\n");
        block_list(first(n));
        return;
    }
    out.append(if k == syntax::N_MESH { "### Mesh " } else { "### Points " });
    code_node(n);
    line_note(n);
    out.append("\n\nOn ");
    words(first(n));
    out.append(".\n\n");
    block_list(next(first(n)));
}

// The note after a statement of a physics: its line in the model. A
// statement of a template has none.
fn statement_note(c: i64, model: bool) {
    if model { line_note(c); }
}

// The statements of a physics (`model`) or a template, from `c`.
fn physics_statements(c0: i64, model: bool) {
    var c = c0;
    while c >= 0 {
        let k = node_kind(c);
        if k == syntax::N_MATASSIGN {
            out.append("Material: ");
            var m = first(c);
            var count = 0;
            while m >= 0 {
                if count > 0 { out.append(", "); }
                code_node(m);
                if first(m) >= 0 {
                    out.append(" on ");
                    words(first(m));
                }
                count++;
                m = next(m);
            }
            statement_note(c, model);
            out.append(".\n\n");
        } else if k == syntax::N_ASSUME {
            out.append("Assumptions");
            statement_note(c, model);
            out.append(":\n\n");
            var a = first(c);
            while a >= 0 {
                out.append("- ");
                code_node(a);
                if first(a) >= 0 {
                    out.append(" with ");
                    inline_list(first(a));
                }
                out.push('\n');
                a = next(a);
            }
            out.push('\n');
        } else if k == syntax::N_LET || k == syntax::N_STATE {
            out.append(if k == syntax::N_LET { "Let " } else { "State " });
            out.push('$');
            math_name(node_text(c));
            out.append("$, of type ");
            type_code(first(c));
            if k == syntax::N_STATE { out.append(", with the initial value"); }
            statement_note(c, model);
            out.append(":\n\n");
            out.append("$$\n");
            math_name(node_text(c));
            out.append(" = ");
            math(next(first(c)));
            out.append("\n$$\n\n");
        } else if k == syntax::N_EQN {
            out.append("Equation");
            statement_note(c, model);
            out.append(":\n\n");
            equation(first(c), next(first(c)));
        } else if k == syntax::N_PLASTIC {
            out.append("Plasticity");
            statement_note(c, model);
            out.append(":\n\n");
            block_list(first(c));
        }
        c = next(c);
    }
}

// The unknowns of a physics, from `c`. It gives the node after them.
fn unknowns(c0: i64) -> i64 {
    var c = c0;
    out.append("Unknowns: ");
    var count = 0;
    while c >= 0 && node_kind(c) == syntax::N_UNKNOWN {
        if count > 0 { out.append(", "); }
        out.push('$');
        math_name(node_text(c));
        out.append("$ of type ");
        type_code(first(c));
        if next(first(c)) >= 0 {
            out.append(", initial value ");
            inline(next(first(c)));
        }
        count++;
        c = next(c);
    }
    out.append(".\n\n");
    return c;
}

fn find_template(name: u8[:]) -> i64 {
    if check::std_prog < 0 { return -1; }
    var c = first(check::std_prog);
    while c >= 0 {
        if node_kind(c) == syntax::N_TEMPLATE && node_text(c) == name { return c; }
        c = next(c);
    }
    return -1;
}

// The value of template parameter `p` (the `i`-th) in the uses `u`: the
// argument by name or by position, or the default.
fn template_arg(u: i64, p: i64, i: i64) -> i64 {
    var a = first(u);
    var k = 0;
    while a >= 0 {
        let named = syntax::nodes[a].aux;
        if named >= 0 && syntax::text(named) == node_text(p) { return first(a); }
        if named < 0 && k == i { return first(a); }
        k++;
        a = next(a);
    }
    return next(first(p));
}

fn physics_item(n: i64) {
    out.append("### Physics ");
    code_node(n);
    line_note(n);
    out.append("\n\n");
    var c = first(n);
    if c >= 0 && node_kind(c) == syntax::N_ON {
        out.append("On ");
        words(first(c));
        out.append(".\n\n");
        c = next(c);
    }
    var uses = -1;
    if c >= 0 && node_kind(c) == syntax::N_USES {
        uses = c;
        c = next(c);
    }
    c = unknowns(c);
    physics_statements(c, true);
    guard uses >= 0 else { return; }
    let t = find_template(node_text(uses));
    out.append("The physics uses the template ");
    code_node(uses);
    out.append(" of the standard library");
    guard t >= 0 else {
        out.append(".\n\n");
        return;
    };
    var p = first(t);
    var i = 0;
    while p >= 0 && node_kind(p) == syntax::N_TPARAM {
        out.append(if i == 0 { ", with " } else { ", " });
        let v = template_arg(uses, p, i);
        out.push('$');
        math_name(node_text(p));
        out.append(" = ");
        if v >= 0 { math(v); }
        out.push('$');
        i++;
        p = next(p);
    }
    out.append(". Its statements are:\n\n");
    while p >= 0 && node_kind(p) == syntax::N_UNKNOWN { p = next(p); }
    physics_statements(p, false);
}

// The table of a run of `bc` statements from `n`. It gives the node after
// the run.
fn bc_table(n0: i64) -> i64 {
    out.append("| Region | Condition | Line |\n|---|---|---|\n");
    var n = n0;
    while n >= 0 && node_kind(n) == syntax::N_BC {
        let region = first(n);
        let target = next(region);
        out.append("| ");
        words(region);
        out.append(" | $");
        math(target);
        out.append(" = ");
        math(next(target));
        format(out, "$ | ", node_line(n), " |\n");
        n = next(n);
    }
    out.push('\n');
    return n;
}

fn pair_item(n: i64) {
    out.append(if node_kind(n) == syntax::N_CONTACT { "### Contact " } else { "### Interface " });
    code_node(n);
    line_note(n);
    out.append("\n\nBetween ");
    words(first(n));
    out.append(" and ");
    words(next(first(n)));
    out.append(".\n\n");
    block_list(next(next(first(n))));
}

// An item with a procedure and a block: `symmetry`, `study`, `manufacture`.
fn procedure_item(n: i64, title: u8[:]) {
    format(out, "### ", title, " ");
    var c = first(n);
    if node_kind(n) == syntax::N_MANUFACTURE {
        out.append("of ");
        words(c);
        c = next(c);
    } else {
        code_node(n);
    }
    line_note(n);
    out.append("\n\nProcedure: ");
    procedure(c);
    out.append(".\n\n");
    block_list(next(c));
}

fn step_item(n: i64) {
    out.append("### Step ");
    code_node(n);
    line_note(n);
    out.append("\n\nProcedure: ");
    procedure(first(n));
    out.append(".\n\n");
    var c = next(first(n));
    while c >= 0 {
        let k = node_kind(c);
        if k == syntax::N_BC {
            c = bc_table(c);
            continue;
        }
        if k == syntax::N_OUTPUT {
            out.append("Output");
            line_note(c);
            out.append(": ");
            var v = first(c);
            var count = 0;
            while v >= 0 && node_kind(v) != syntax::N_EVERY {
                if count > 0 { out.append(", "); }
                inline(v);
                count++;
                v = next(v);
            }
            if v >= 0 {
                out.append(", every ");
                inline(first(v));
            }
            out.append(".\n\n");
        } else if k == syntax::N_REPORT {
            out.append("Report");
            line_note(c);
            out.append(": ");
            var r = first(c);
            var count = 0;
            while r >= 0 {
                if count > 0 { out.append(", "); }
                if syntax::nodes[r].aux >= 0 {
                    code(syntax::text(syntax::nodes[r].aux));
                    out.append(" = ");
                }
                inline(first(r));
                count++;
                r = next(r);
            }
            out.append(".\n\n");
        } else if k == syntax::N_RELEASE {
            out.append("Release ");
            words(first(c));
            line_note(c);
            out.append(".\n\n");
        } else if k == syntax::N_BIRTH || k == syntax::N_DEATH {
            out.append(if k == syntax::N_BIRTH { "Birth of " } else { "Death of " });
            var r = first(c);
            var count = 0;
            while r >= 0 {
                if node_kind(r) == syntax::N_WHERE {
                    out.append(" where ");
                    inline(first(r));
                } else {
                    if count > 0 { out.append(", "); }
                    words(r);
                    count++;
                }
                r = next(r);
            }
            line_note(c);
            out.append(".\n\n");
        } else if k == syntax::N_SOLVE || k == syntax::N_ADAPT {
            out.append(if k == syntax::N_SOLVE { "Solver settings" } else { "Adaptivity" });
            line_note(c);
            out.append(":\n\n");
            block_list(first(c));
        }
        c = next(c);
    }
}

// The table of a run of `expect` items from `n`. It gives the item after
// the run.
fn expect_table(n0: i64) -> i64 {
    out.append("| Name | Quantity | Expected value | Limit | Line |\n|---|---|---|---|---|\n");
    var n = n0;
    while n >= 0 && node_kind(n) == syntax::N_EXPECT {
        let l = first(n);
        let r = next(l);
        out.append("| ");
        code_node(n);
        out.append(" | ");
        inline(l);
        out.append(" | ");
        inline(r);
        out.append(" | ");
        inline(next(r));
        format(out, " | ", node_line(n), " |\n");
        n = next(n);
    }
    out.push('\n');
    return n;
}

fn table_item(n: i64) {
    out.append("### Table ");
    code_node(n);
    line_note(n);
    out.append("\n\n|");
    var c = first(n);
    var columns = 0;
    while c >= 0 && node_kind(c) == syntax::N_COLUMN {
        out.push(' ');
        code_node(c);
        out.append(" (");
        type_code(first(c));
        let u = next(first(c));
        if u >= 0 {
            out.append(", $");
            math_unit(u);
            out.push('$');
        }
        out.append(") |");
        columns++;
        c = next(c);
    }
    out.append("\n|");
    for i in columns { out.append("---|"); }
    out.push('\n');
    while c >= 0 {
        out.push('|');
        var v = first(c);
        while v >= 0 {
            out.push(' ');
            inline(v);
            out.append(" |");
            v = next(v);
        }
        out.push('\n');
        c = next(c);
    }
    out.push('\n');
}

fn other_item(n: i64) {
    let k = node_kind(n);
    if k == syntax::N_ALLOW {
        out.append("### Allowed check ");
        code_node(n);
        line_note(n);
        out.append("\n\n");
        let s = syntax::text(syntax::nodes[n].aux);
        text(s[1..s.len - 1]);
        out.append("\n\n");
    } else if k == syntax::N_TOLERANCE {
        out.append("### Tolerance ");
        code_node(n);
        line_note(n);
        out.append("\n\nThe value ");
        let v = first(n);
        let lo = next(v);
        inline(v);
        out.append(" is in ");
        inline(lo);
        out.append(" to ");
        inline(next(lo));
        out.append(".\n\n");
        block_list(next(next(lo)));
    } else if k == syntax::N_FN {
        out.append("### Function ");
        code_node(n);
        line_note(n);
        out.append("\n\n$$\n");
        math_name(node_text(n));
        out.push('(');
        var p = first(n);
        var count = 0;
        while node_kind(p) == syntax::N_FNPARAM {
            if count > 0 { out.append(", "); }
            math_name(node_text(p));
            count++;
            p = next(p);
        }
        out.append(") = ");
        math(next(p));
        out.append("\n$$\n\n");
    } else if k == syntax::N_TABLE {
        table_item(n);
    } else if k == syntax::N_SYMMETRY || k == syntax::N_STUDY || k == syntax::N_MANUFACTURE {
        procedure_item(n, if k == syntax::N_SYMMETRY { "Symmetry" } else if k == syntax::N_STUDY { "Study" } else { "Manufacture" });
    } else if k == syntax::N_CONTACT || k == syntax::N_INTERFACE {
        pair_item(n);
    } else if k == syntax::N_USE {
        out.append("### Use");
        line_note(n);
        out.append("\n\n");
        code_node(n);
        out.append("\n\n");
    }
}

// --- the document -------------------------------------------------------------------------

// The base units of the named unit systems (spec §6.5).
let SYSTEMS = [
    "SI", "m", "kg", "s", "N", "Pa", "K",
    "mm_N_MPa_s", "mm", "t", "s", "N", "MPa", "K",
    "mm_kg_ms", "mm", "kg", "ms", "kN", "GPa", "K",
];
let QUANTITIES = ["length", "mass", "time", "force", "stress", "temperature"];

fn space_words(s: u8[:]) -> const u8[:] {
    if s == "1d" { return "one-dimensional"; }
    if s == "2d" { return "two-dimensional"; }
    if s == "3d" { return "three-dimensional"; }
    if s == "axi" { return "axisymmetric"; }
    return s;
}

fn header(h: i64) {
    let sp = syntax::child_of_kind(h, syntax::N_SPACE);
    let us = syntax::child_of_kind(h, syntax::N_UNITSYS);
    out.append("## Model\n\nThe model ");
    code_node(h);
    line_note(h);
    format(out, " is ", space_words(node_text(sp)), ". ");
    let name = node_text(us);
    if name == "nondim" {
        out.append("It has no units: each quantity is a `Real`.\n\n");
        return;
    }
    if name == "{" {
        out.append("Its unit system gives the base units:");
        var f = first(us);
        var count = 0;
        while f >= 0 {
            out.append(if count > 0 { ", " } else { " " });
            text(node_text(f));
            out.append(" in $");
            math_unit(first(f));
            out.push('$');
            count++;
            f = next(f);
        }
        out.append(".\n\n");
        return;
    }
    out.append("Its unit system is ");
    code(name);
    for s, i in SYSTEMS {
        if i % 7 != 0 || s != name { continue; }
        out.append(":");
        for q, j in QUANTITIES {
            out.append(if j > 0 { ", " } else { " " });
            format(out, q, " in ", SYSTEMS[i + 1 + j]);
        }
    }
    out.append(". The unit system sets the units of the output. Each value in the source carries its own unit.\n\n");
}

fn scales(n: i64) {
    out.append("The reference values of the scales line");
    line_note(n);
    out.append(" are ");
    inline_list(first(n));
    out.append(". The solver divides each quantity by the product of powers of these values that has its dimension.\n\n");
}

// The `//` comments at the top of the file, before the header, as the
// introduction. The lines of one comment block make one paragraph.
fn introduction() {
    var started = false;
    var last_line = 0;
    for t, i in syntax::tokens {
        if t.file != 0 { break; }
        if t.kind == syntax::TK_NL { continue; }
        if t.kind != syntax::TK_COMMENT && t.kind != syntax::TK_DOC_COMMENT { break; }
        var s = syntax::text(i);
        s = trim(if starts_with(s, "///") { s[3..] } else { s[2..] });
        if started { out.append(if t.line > last_line + 1 { "\n\n" } else { " " }); }
        if s.len > 0 && (s[0] == '-' || s[0] == '+') { out.push('\\'); }
        text(s);
        started = true;
        last_line = t.line;
    }
    if started { out.append("\n\n"); }
}

// Writes the model `prog` of `path` to `out`. The checker has checked it.
fn write_model(path: u8[:], prog: i64) {
    out.clear();
    let h = syntax::child_of_kind(prog, syntax::N_HEADER);
    out.append("# The model ");
    code_node(h);
    out.append("\n\nThis document is the model in the file ");
    code(path);
    out.append(", written by `skein md`. The items are in the order of the file, and each item gives its line.\n\n");
    introduction();
    var current = G_NONE;
    var n = first(prog);
    while n >= 0 {
        let k = node_kind(n);
        if k == syntax::N_HEADER {
            header(n);
            n = next(n);
            continue;
        }
        let g = group(n);
        if g != G_NONE && g != current {
            format(out, "## ", GROUP_TITLES[g], "\n\n");
            current = g;
        }
        if value_row(n) {
            n = value_table(n);
            continue;
        }
        match k {
            syntax::N_LET => { n = geometry_list(n); },
            syntax::N_BC => { n = bc_table(n); },
            syntax::N_EXPECT => { n = expect_table(n); },
            _ => {
                match k {
                    syntax::N_SCALES => { scales(n); },
                    syntax::N_MATERIAL => { material_item(n); },
                    syntax::N_MESH, syntax::N_POINTS, syntax::N_METHOD => { discretization_item(n); },
                    syntax::N_PHYSICS => { physics_item(n); },
                    syntax::N_STEP => { step_item(n); },
                    _ => { other_item(n); },
                }
                n = next(n);
            },
        }
    }
    // One new line at the end, not a blank line.
    while out.len > 1 && out[out.len - 1] == '\n' && out[out.len - 2] == '\n' { out.resize(out.len - 1); }
}

// The first program in the syntax tables: the model. The standard library
// comes after it.
fn model_program() -> i64 {
    for node, i in syntax::nodes {
        if node.kind == syntax::N_PROGRAM { return i; }
    }
    return -1;
}

// Checks the model `text` of `path`, and makes its document in `result`.
// It gives false, with the diagnostics in `check::diags`, when the model has
// errors.
fn md_text(path: u8[:], text: u8[:], result: u8[>..]&) -> bool {
    check::check_text(path, text);
    if check::diags.len > 0 { return false; }
    write_model(path, model_program());
    result.clear();
    result.append(out);
    return true;
}

let USAGE = "usage: skein md model.skn > model.md\n";

// `skein md`. `args` are the words after `md`. It writes the document to
// standard output, and gives the exit status (command_line.md §10).
fn run(args: const u8[:][:]) -> i64 {
    guard args.len == 1 && args[0].len > 0 && args[0][0] != '-' else {
        write_stderr(str("skein md: give one model file\n", USAGE));
        return 2;
    }
    let path = args[0];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein md: cannot read '", path, "'\n"));
        return 1;
    };
    var result: u8[>..] = [];
    guard md_text(path, text, result) else {
        var report: u8[>..] = [];
        check::write_toon(report);
        write_stderr(report);
        return 1;
    };
    write_stdout(result);
    return 0;
}
