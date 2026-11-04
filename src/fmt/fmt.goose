import std;
import os;
namespace fmt;
import syntax from "../syntax";

let WIDTH = 100;
let INDENT = 4;
let ADDED_COMMA = -2;
let NONE = -1;
let COMMENT_ENTRY = -3;

let EQ_SIDES = 1;

var prec: i64[>..];
var listopen: bool[>..];
var blockopen: bool[>..];
var tight: bool[>..];
var multi: bool[>..];
var parens: i64[>..];
var partner: i64[>..];

var pool: i64[>..];
var out: u8[>..];
var last_tok = -1;

fn kind(t: i64) -> i64 {
    if t == ADDED_COMMA { return syntax::TK_COMMA; }
    return syntax::tokens[t].kind;
}

fn line(t: i64) -> i64 { syntax::tokens[t].line }

fn precedence(t: i64) -> i64 { if t == ADDED_COMMA { 0 } else { prec[t] } }

fn opens_list(t: i64) -> bool { t != ADDED_COMMA && listopen[t] }

fn is_comment(t: i64) -> bool {
    let k = kind(t);
    return k == syntax::TK_COMMENT || k == syntax::TK_DOC_COMMENT;
}

fn is_open(k: i64) -> bool {
    k == syntax::TK_LPAREN || k == syntax::TK_LBRACKET || k == syntax::TK_LBRACE || k == syntax::TK_LANGLE
}

fn is_close(k: i64) -> bool {
    k == syntax::TK_RPAREN || k == syntax::TK_RBRACKET || k == syntax::TK_RBRACE || k == syntax::TK_RANGLE
}

fn token_text(t: i64) -> const u8[:] {
    if t == ADDED_COMMA { return ","; }
    let s = syntax::text(t);
    if is_comment(t) { return trim_end(s); }
    return s;
}

fn code_points(s: u8[:]) -> i64 {
    var n = 0;
    for c in s {
        if (c & 0xC0) != 0x80 { n++; }
    }
    return n;
}

fn next_code(t: i64) -> i64 {
    var k = t + 1;
    while k < syntax::tokens.len && is_comment(k) { k++; }
    return k;
}

fn next_row(t: i64) -> bool {
    var k = t + 1;
    while k < syntax::tokens.len && (is_comment(k) || kind(k) == syntax::TK_NL) { k++; }
    return kind(k) == syntax::TK_INDENT;
}

fn op_prec(k: i64) -> i64 {
    match k {
        syntax::TK_OROR, syntax::TK_BAR => 2,
        syntax::TK_ANDAND, syntax::TK_AMP => 3,
        syntax::TK_EQ, syntax::TK_NE, syntax::TK_LT, syntax::TK_LE, syntax::TK_GT, syntax::TK_GE => 4,
        syntax::TK_PLUS, syntax::TK_MINUS => ADDITIVE,
        syntax::TK_STAR, syntax::TK_SLASH => 6,
        syntax::TK_CARET => 7,
        _ => 0,
    }
}

fn annotate() {
    let n = syntax::tokens.len;
    prec.clear();
    listopen.clear();
    blockopen.clear();
    tight.clear();
    multi.clear();
    parens.clear();
    partner.clear();
    var stack: i64[>..] = [];
    var open = 0;
    for t in 0..n {
        let k = kind(t);
        let binary = k == syntax::TK_MINUS || k == syntax::TK_STAR || k == syntax::TK_SLASH || k == syntax::TK_CARET;
        prec.push(if binary { 0 } else { op_prec(k) });
        listopen.push(false);
        blockopen.push(false);
        tight.push(false);
        multi.push(false);
        partner.push(-1);
        if is_close(k) {
            if k != syntax::TK_RBRACE && open > 0 { open--; }
            if stack.len > 0 {
                let o = stack[stack.len - 1];
                stack.resize(stack.len - 1);
                partner[o] = t;
                partner[t] = o;
            }
        }
        parens.push(open);
        if is_open(k) {
            stack.push(t);
            if k != syntax::TK_LBRACE { open++; }
        }
    }
    for node in syntax::nodes {
        let t = node.tok;
        match node.kind {
            syntax::N_BINARY => { prec[t] = op_prec(kind(t)); },
            syntax::N_CALL, syntax::N_VECTOR, syntax::N_LIST => { listopen[t] = true; },
            syntax::N_PROC, syntax::N_USES, syntax::N_FN, syntax::N_TEMPLATE => {
                let p = next_code(t);
                if p < n && kind(p) == syntax::TK_LPAREN { listopen[p] = true; }
            },
            syntax::N_BLOCK, syntax::N_RECORD => { blockopen[t] = true; },
            syntax::N_TABLE => {
                let b = next_code(t);
                if b < n && kind(b) == syntax::TK_LBRACE { tight[b] = true; }
            },
            _ => {},
        }
    }
    for t in 0..n {
        if kind(t) != syntax::KW_EQ { continue; }
        var level = 0;
        var k = t + 1;
        while k < n && kind(k) != syntax::TK_NL && kind(k) != syntax::TK_EOF {
            let kk = kind(k);
            if is_open(kk) { level++; }
            if is_close(kk) { level--; }
            if kk == syntax::TK_ASSIGN && level == 0 { prec[k] = EQ_SIDES; break; }
            k++;
        }
    }
    for t in 0..n {
        if kind(t) != syntax::TK_LBRACE || partner[t] < 0 { continue; }
        let c = partner[t];
        if t + 1 < c && (is_comment(t + 1) || line(t + 1) > line(t)) { multi[t] = true; continue; }
        for k in t + 1..c {
            if kind(k) == syntax::TK_NL { multi[t] = true; break; }
        }
    }
}

fn space(a: i64, b: i64) -> bool {
    let ka = kind(a);
    let kb = kind(b);
    if kb == syntax::TK_COMMA || kb == syntax::TK_COLON || kb == syntax::TK_DOT { return false; }
    if kb == syntax::TK_RPAREN || kb == syntax::TK_RBRACKET || kb == syntax::TK_RANGLE || kb == syntax::TK_LANGLE {
        return false;
    }
    if ka == syntax::TK_LPAREN || ka == syntax::TK_LBRACKET || ka == syntax::TK_LANGLE || ka == syntax::TK_DOT {
        return false;
    }
    if ka == syntax::TK_DOTDOT || kb == syntax::TK_DOTDOT { return false; }
    if ka == syntax::TK_LBRACE { return kb != syntax::TK_RBRACE && !tight[a]; }
    if kb == syntax::TK_RBRACE { return partner[b] < 0 || !tight[partner[b]]; }
    if kb == syntax::TK_LBRACE { return !tight[b]; }
    if ka == syntax::TK_STAR || ka == syntax::TK_SLASH || ka == syntax::TK_CARET { return false; }
    if kb == syntax::TK_STAR || kb == syntax::TK_SLASH || kb == syntax::TK_CARET { return false; }
    if ka == syntax::TK_MINUS { return prec[a] > 0; }
    if ka == syntax::TK_BANG { return false; }
    if kb == syntax::TK_LPAREN {
        return !(ka == syntax::TK_IDENT || ka == syntax::TK_RPAREN || ka == syntax::TK_RBRACKET);
    }
    return true;
}

fn width(s: i64, n: i64) -> i64 {
    var w = 0;
    for i in s..s + n {
        let t = pool[i];
        if i > s && space(pool[i - 1], t) { w++; }
        w += code_points(token_text(t));
    }
    return w;
}

fn write_indent(ind: i64) {
    for i in 0..ind { out.push(' '); }
}

fn write_line(s: i64, n: i64, ind: i64) {
    write_indent(ind);
    for i in s..s + n {
        let t = pool[i];
        if i > s && space(pool[i - 1], t) { out.push(' '); }
        out.append(token_text(t));
    }
    out.push('\n');
}

fn trailing_comment(c: i64) {
    if out.len > 0 && out[out.len - 1] == '\n' { out.resize(out.len - 1); }
    out.append("  ");
    out.append(token_text(c));
    out.push('\n');
    last_tok = c;
}

fn depths(s: i64, n: i64, d: i64[>..]&) -> i64 {
    d.clear();
    var cur = 0;
    var low = 0;
    for i in 0..n {
        let k = kind(pool[s + i]);
        if is_close(k) { cur--; }
        d.push(cur);
        low = min(low, cur);
        if is_open(k) { cur++; }
    }
    return low;
}

fn close_at(s: i64, n: i64, i: i64) -> i64 {
    let c = partner[pool[s + i]];
    for j in i + 1..n {
        if pool[s + j] == c { return j; }
    }
    return -1;
}

var cut_pool: i64[>..];

let WRITE = 0;
let CUT = 1;
let EXPLODE = 2;

let ADDITIVE = 5;

fn add_cuts(s: i64, n: i64, d: i64[:], level: i64, p: i64) {
    for i in 0..n - 1 {
        if d[i] != level { continue; }
        let t = pool[s + i];
        if (p == 0 && kind(t) == syntax::TK_COMMA) || (p > 0 && precedence(t) == p) { cut_pool.push(i); }
    }
}

fn low_prec(s: i64, n: i64, d: i64[:], level: i64, most: i64) -> i64 {
    var low = 0;
    for i in 0..n - 1 {
        let p = precedence(pool[s + i]);
        if d[i] == level && p > 0 && p <= most && (low == 0 || p < low) { low = p; }
    }
    return low;
}

fn widest_list(s: i64, n: i64, ind: i64, d: i64[:], top: i64, fit: bool) -> i64, i64 {
    var best = -1;
    var best_close = -1;
    var best_width = -1;
    let all = width(s, n);
    for i in 0..n {
        if d[i] != top || !opens_list(pool[s + i]) { continue; }
        let j = close_at(s, n, i);
        if j < 0 { continue; }
        let inner = width(s + i + 1, j - i - 1);
        if fit {
            let short = ind + width(s, i + 1) <= WIDTH && ind + width(s + j, n - j) <= WIDTH;
            if !short || inner < all - inner { continue; }
        }
        if inner > best_width {
            best = i;
            best_close = j;
            best_width = inner;
        }
    }
    return best, best_close;
}

fn plan(s: i64, n: i64, ind: i64) -> i64, i64, i64 {
    if ind + width(s, n) <= WIDTH { return WRITE, 0, 0; }
    var d: i64[>..] = [];
    let top = depths(s, n, d);
    var deepest = top;
    for x in d { deepest = max(deepest, x); }
    let c0 = cut_pool.len;

    add_cuts(s, n, d, top, 0);
    if cut_pool.len > c0 { return CUT, c0, cut_pool.len - c0; }

    let o, c = widest_list(s, n, ind, d, top, true);
    if o >= 0 { return EXPLODE, o, c; }

    let low = low_prec(s, n, d, top, ADDITIVE);
    if low > 0 {
        add_cuts(s, n, d, top, low);
        return CUT, c0, cut_pool.len - c0;
    }

    let o2, c2 = widest_list(s, n, ind, d, top, false);
    if o2 >= 0 { return EXPLODE, o2, c2; }

    for level in top + 1..deepest + 1 {
        let p = low_prec(s, n, d, level, ADDITIVE);
        if p > 0 {
            add_cuts(s, n, d, level, p);
            return CUT, c0, cut_pool.len - c0;
        }
    }

    let high = low_prec(s, n, d, top, 100);
    if high > 0 {
        add_cuts(s, n, d, top, high);
        return CUT, c0, cut_pool.len - c0;
    }

    return WRITE, 0, 0;
}

recursive fn lay(s: i64, n: i64, ind: i64) {
    if n == 0 { return; }
    let how, a, b = plan(s, n, ind);
    if how == CUT {
        cut(s, n, ind, a, b);
    } else if how == EXPLODE {
        explode(s, n, ind, a, b);
    } else {
        write_line(s, n, ind);
    }
}

recursive fn cut(s: i64, n: i64, ind: i64, cs: i64, cn: i64) {
    var from = 0;
    for k in 0..cn {
        let c = cut_pool[cs + k];
        lay(s + from, c + 1 - from, if k == 0 { ind } else { ind + INDENT });
        from = c + 1;
    }
    lay(s + from, n - from, ind + INDENT);
}

recursive fn explode(s: i64, n: i64, ind: i64, o: i64, c: i64) {
    lay(s, o + 1, ind);
    var level = 0;
    var from = o + 1;
    for i in o + 1..c + 1 {
        let k = kind(pool[s + i]);
        let end = i == c || (level == 0 && k == syntax::TK_COMMA);
        if is_open(k) { level++; }
        if is_close(k) { level--; }
        if !end { continue; }
        if i > from {
            let item = pool.len;
            for j in from..i { pool.push(pool[s + j]); }
            pool.push(ADDED_COMMA);
            lay(item, i - from + 1, ind + INDENT);
        }
        from = i + 1;
    }
    lay(s + c, n - c, ind);
}

fn write_as_written(a: i64, b: i64, ind: i64, add_comma: bool) {
    var prev = -1;
    var last_code = -1;
    for t in a..b {
        if !is_comment(t) { last_code = t; }
    }
    for t in a..b {
        if prev >= 0 && line(t) == line(prev) {
            if is_comment(t) {
                out.append("  ");
            } else if space(prev, t) {
                out.push(' ');
            }
        } else {
            if prev >= 0 { out.push('\n'); }
            write_indent(if prev < 0 || is_close(kind(t)) { ind } else { ind + INDENT });
        }
        out.append(token_text(t));
        if t == last_code && add_comma { out.push(','); }
        prev = t;
    }
    out.push('\n');
    last_tok = b - 1;
}

fn ends_statement(t: i64, level: i64, sep_comma: bool) -> bool {
    let k = kind(t);
    if k == syntax::TK_EOF { return true; }
    if level != 0 { return false; }
    return k == syntax::TK_NL || k == syntax::TK_RBRACE || (k == syntax::TK_COMMA && sep_comma);
}

fn dropped_comma(t: i64) -> bool {
    let c = next_code(t);
    if c >= syntax::tokens.len { return false; }
    let k = kind(c);
    if k == syntax::TK_RBRACE { return partner[c] >= 0 && !multi[partner[c]]; }
    return k == syntax::TK_RPAREN || k == syntax::TK_RBRACKET;
}

recursive fn statement(i0: i64, ind: i64, sep_comma: bool, add_comma: bool) -> i64 {
    let first = kind(i0);
    let one_line = first == syntax::KW_MODEL || first == syntax::KW_TABLE;
    var i = i0;
    var level = 0;
    loop {
        let seg = pool.len;
        let seg_first = i;
        let seg_level = level;
        var as_written = false;
        var comment = -1;
        var last_code = -1;
        var brace = false;
        loop {
            if is_comment(i) {
                if last_code >= 0 && line(i) == line(last_code) {
                    comment = i;
                } else if ends_statement(next_code(i), level, sep_comma) {
                    break;
                } else {
                    as_written = true;
                }
                i++;
                continue;
            }
            if ends_statement(i, level, sep_comma) { break; }
            if comment >= 0 { as_written = true; }
            let k = kind(i);
            if k == syntax::TK_COMMA && dropped_comma(i) { i++; continue; }
            pool.push(i);
            last_code = i;
            if is_open(k) { level++; }
            if is_close(k) { level--; }
            i++;
            if k == syntax::TK_LBRACE && multi[i - 1] { brace = true; break; }
        }
        let n = pool.len - seg;
        let last = !brace;
        if !as_written && !one_line && ind + width(seg, n) > WIDTH {
            var found = -1;
            for p in seg..seg + n {
                let t = pool[p];
                if blockopen[t] && !multi[t] { found = t; break; }
            }
            if found >= 0 {
                multi[found] = true;
                pool.resize(seg);
                i = seg_first;
                level = seg_level;
                continue;
            }
        }
        if as_written {
            write_as_written(seg_first, i, ind, last && add_comma);
        } else {
            if last && add_comma { pool.push(ADDED_COMMA); }
            let m = pool.len - seg;
            if one_line { write_line(seg, m, ind); } else { lay(seg, m, ind); }
            last_tok = if last_code >= 0 { last_code } else { last_tok };
            if comment >= 0 { trailing_comment(comment); }
        }
        if last { break; }
        let b = i - 1;
        if i < syntax::tokens.len && is_comment(i) && line(i) == line(b) {
            trailing_comment(i);
            i++;
        }
        i = body(i, ind + INDENT, blockopen[b], blockopen[b] && parens[b] > 0, false);
    }
    if kind(i) == syntax::TK_NL || (sep_comma && level == 0 && kind(i) == syntax::TK_COMMA) {
        i++;
        if i < syntax::tokens.len && is_comment(i) && line(i) == line(i - 1) {
            trailing_comment(i);
            i++;
        }
    }
    if first == syntax::KW_TABLE {
        while kind(i) == syntax::TK_INDENT || kind(i) == syntax::TK_NL || (is_comment(i) && next_row(i)) {
            if kind(i) == syntax::TK_NL { i++; continue; }
            if is_comment(i) {
                write_indent(INDENT);
                out.append(token_text(i));
                out.push('\n');
                last_tok = i;
                i++;
                continue;
            }
            i++;
            let seg = pool.len;
            var comment = -1;
            while kind(i) != syntax::TK_NL && kind(i) != syntax::TK_EOF {
                if is_comment(i) { comment = i; } else { pool.push(i); last_tok = i; }
                i++;
            }
            write_line(seg, pool.len - seg, INDENT);
            if comment >= 0 { trailing_comment(comment); }
        }
    }
    return i;
}

fn is_let(k: i64) -> bool { k == syntax::KW_LET || k == syntax::KW_PARAM }

fn blank_between(prev: i64, next: i64, author: bool, top: bool) -> bool {
    if prev == NONE { return false; }
    if prev == COMMENT_ENTRY || next == COMMENT_ENTRY { return author; }
    if prev == syntax::KW_BC && next == syntax::KW_BC { return false; }
    if !top { return author; }
    if is_let(prev) && is_let(next) { return author; }
    return true;
}

recursive fn body(i0: i64, ind: i64, sep_comma: bool, add_comma: bool, top: bool) -> i64 {
    var i = i0;
    var prev = NONE;
    loop {
        while kind(i) == syntax::TK_NL || (sep_comma && kind(i) == syntax::TK_COMMA) { i++; }
        let k = kind(i);
        if k == syntax::TK_EOF || k == syntax::TK_RBRACE { break; }
        let author = last_tok >= 0 && line(i) - line(last_tok) >= 2;
        if is_comment(i) {
            var j = i;
            while j + 1 < syntax::tokens.len && is_comment(j + 1) && line(j + 1) == line(j) + 1 { j++; }
            let next = next_code(j);
            let own = line(next) > line(j) + 1 || kind(next) == syntax::TK_RBRACE || kind(next) == syntax::TK_EOF;
            let entry = if own { COMMENT_ENTRY } else { kind(next) };
            if blank_between(prev, entry, author, top) { out.push('\n'); }
            for c in i..j + 1 {
                write_indent(ind);
                out.append(token_text(c));
                out.push('\n');
            }
            last_tok = j;
            i = j + 1;
            if own { prev = COMMENT_ENTRY; continue; }
            i = statement(i, ind, sep_comma, add_comma);
            prev = entry;
            continue;
        }
        if blank_between(prev, k, author, top) { out.push('\n'); }
        i = statement(i, ind, sep_comma, add_comma);
        prev = k;
    }
    return i;
}

fn shape(prog: i64, into: u8[>..]&) {
    syntax::dump(prog, into);
    for t, i in syntax::tokens {
        if is_comment(i) {
            into.push('\n');
            into.append(token_text(i));
        }
    }
}

fn parse_text(path: u8[:], text: u8[:]) -> i64 {
    syntax::reset();
    let f = syntax::add_file(path, text);
    syntax::lex(f);
    let prog = syntax::parse(f);
    if syntax::errors.len > 0 { return -1; }
    return prog;
}

let OK = 0;
let SYNTAX_ERROR = 1;
let INTERNAL_FAULT = 3;

fn format_text(path: u8[:], text: u8[:], result: u8[>..]&) -> i64 {
    let prog = parse_text(path, text);
    if prog < 0 { return SYNTAX_ERROR; }
    var before: u8[>..] = [];
    shape(prog, before);
    annotate();
    pool.clear();
    cut_pool.clear();
    out.clear();
    last_tok = -1;
    body(0, 0, false, false, true);
    result.clear();
    result.append(out);
    let again = parse_text(path, result);
    if again < 0 { return INTERNAL_FAULT; }
    var after: u8[>..] = [];
    shape(again, after);
    if after != before { return INTERNAL_FAULT; }
    return OK;
}

fn toon_value(into: u8[>..]&, s: u8[:], quote: bool) {
    var needs = quote || s.len == 0;
    for c in s {
        if c == ',' || c == '"' || c == ':' || c == '\\' || c == ' ' { needs = true; }
    }
    if !needs { into.append(s); return; }
    into.push('"');
    for c in s {
        if c == '"' || c == '\\' { into.push('\\'); }
        into.push(c);
    }
    into.push('"');
}

fn write_errors() {
    var s: u8[>..] = [];
    format(s, "diagnostics[", syntax::errors.len, "]{class,id,file,line,col,msg,fix}:\n");
    for e in syntax::errors {
        let t = syntax::tokens[e.tok];
        s.append("  error,");
        s.append(e.id);
        s.push(',');
        toon_value(s, syntax::files[t.file].path, false);
        format(s, ",", t.line, ",", t.col, ",");
        toon_value(s, e.msg, true);
        s.push(',');
        toon_value(s, e.fix, true);
        s.push('\n');
    }
    format(s, "summary: errors=", syntax::errors.len, " lints=0\n");
    write_stdout(s);
}

let USAGE = "usage: skein fmt [--check] [files]\n";

fn run(args: const u8[:][:]) -> i64 {
    var check = false;
    var paths: const u8[:][>..] = [];
    for a in args {
        if a == "--check" {
            check = true;
        } else if a.len > 1 && a[0] == '-' {
            write_stderr(str("skein fmt: unknown option '", a, "'\n", USAGE));
            return 2;
        } else {
            paths.push(a);
        }
    }
    if paths.len == 0 { paths.push("-"); }
    var status = 0;
    for path in paths {
        var text: u8[>..] = [];
        if path == "-" {
            read_stdin(text);
        } else if !read_file(path, text) {
            write_stderr(str("skein fmt: cannot read '", path, "'\n"));
            status = max(status, 1);
            continue;
        }
        var result: u8[>..] = [];
        let r = format_text(path, text, result);
        if r == SYNTAX_ERROR {
            write_errors();
            status = max(status, 1);
            continue;
        }
        if r == INTERNAL_FAULT {
            write_stderr(str("skein fmt: internal fault: the layout of '", path,
                             "' changes its meaning. Please report it.\n"));
            status = 3;
            continue;
        }
        let changed = result != text;
        if check {
            if changed {
                write_stdout(str(path, "\n"));
                status = max(status, 1);
            }
        } else if path == "-" {
            write_stdout(result);
        } else if changed && !write_file(path, result) {
            write_stderr(str("skein fmt: cannot write '", path, "'\n"));
            status = max(status, 1);
        }
    }
    return status;
}
