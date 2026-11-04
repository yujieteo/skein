import std;
namespace codegen;
import syntax from "../syntax";

let MAX_FN_LINES = 60;

struct Stmt { lo: i64, hi: i64 }

fn depth_change(line: u8[:]) -> i64 {
    var d = 0;
    var quoted = false;
    var i = 0;
    while i < line.len {
        let c = line[i];
        if quoted {
            if c == '\\' { i++; } else if c == '"' { quoted = false; }
        } else if c == '"' {
            quoted = true;
        } else if c == '/' && i + 1 < line.len && line[i + 1] == '/' {
            return d;
        } else if c == '{' || c == '(' || c == '[' {
            d++;
        } else if c == '}' || c == ')' || c == ']' {
            d--;
        }
        i++;
    }
    return d;
}

fn statements(body: const u8[:][:], out: Stmt[>..]&) -> bool {
    var lo = 0;
    var depth = 0;
    for line, i in body {
        let t = trim(line);
        if has_word(t, "return") || has_word(t, "guard") { return false; }
        if depth == 0 && starts_with(t, "//") { continue; }
        depth += depth_change(line);
        if depth == 0 {
            out.push(Stmt { lo: lo, hi: i + 1 });
            lo = i + 1;
        }
    }
    return depth == 0 && lo == body.len;
}

fn has_word(s: u8[:], w: u8[:]) -> bool {
    var at = 0;
    while at < s.len {
        let i = find(s[at..], w);
        if i < 0 { return false; }
        let lo = at + i;
        let hi = lo + w.len;
        let before = lo > 0 && (syntax::is_letter(s[lo - 1]) || syntax::is_digit(s[lo - 1]));
        let after = hi < s.len && (syntax::is_letter(s[hi]) || syntax::is_digit(s[hi]));
        if !before && !after { return true; }
        at = hi;
    }
    return false;
}

fn used_below(body: const u8[:][:], s: Stmt, from: i64) -> bool {
    var t = trim(body[s.hi - 1]);
    for i in s.lo..s.hi {
        let u = trim(body[i]);
        if !starts_with(u, "//") { t = u; break; }
    }
    if !starts_with(t, "let ") && !starts_with(t, "var ") { return false; }
    var n = 4;
    while n < t.len && t[n] != '=' && t[n] != ':' { n++; }
    for name in split(t[4..n], ',') {
        for line in body[from..] {
            if has_word(line, trim(name)) { return true; }
        }
    }
    return false;
}

fn cuts(body: const u8[:][:], stmts: Stmt[:]) -> i64[>..] {
    var free: bool[>..] = [];
    for k in 0..stmts.len {
        var ok = true;
        for j in 0..k {
            if used_below(body, stmts[j], stmts[k].lo) { ok = false; break; }
        }
        free.push(ok);
    }
    var starts: i64[>..] = [0];
    var last = 0;
    for k in 0..stmts.len {
        if stmts[k].hi - stmts[k].lo + 2 > MAX_FN_LINES { return []; }
        if k > 0 && free[k] { last = k; }
        if stmts[k].hi - stmts[starts[starts.len - 1]].lo + 2 <= MAX_FN_LINES { continue; }
        if last == starts[starts.len - 1] || stmts[k].hi - stmts[last].lo + 2 > MAX_FN_LINES { return []; }
        starts.push(last);
    }
    return starts;
}

fn split_long_functions(program: u8[>..]&) {
    var text: u8[>..] = [];
    text.append(program);
    program.clear();
    let lines = split(text, '\n');
    var i = 0;
    while i < lines.len {
        var end = i;
        if starts_with(lines[i], "fn ") && ends_with(lines[i], "() {") {
            while end < lines.len && lines[end] != "}" { end++; }
        }
        if end > i && end < lines.len && end - i + 1 > MAX_FN_LINES && split_function(lines[i..end + 1], program) {
            i = end + 1;
            continue;
        }
        program.append(lines[i]);
        if i + 1 < lines.len { program.push('\n'); }
        i++;
    }
}

fn split_function(f: const u8[:][:], out: u8[>..]&) -> bool {
    let name = f[0]["fn ".len..f[0].len - "() {".len];
    let body = f[1..f.len - 1];
    var stmts: Stmt[>..] = [];
    if !statements(body, stmts) { return false; }
    let starts = cuts(body, stmts);
    if starts.len < 2 || starts.len + 2 > MAX_FN_LINES { return false; }
    for p, k in starts {
        let hi = if k + 1 < starts.len { stmts[starts[k + 1]].lo } else { body.len };
        format(out, "fn part", k + 1, "_", name, "() {\n");
        for line in body[stmts[p].lo..hi] { format(out, line, "\n"); }
        out.append("}\n");
    }
    format(out, "fn ", name, "() {\n");
    for k in 0..starts.len { format(out, "    part", k + 1, "_", name, "();\n"); }
    out.append("}\n");
    return true;
}
