import std;
namespace codegen;

struct Item { lo: i64, hi: i64, keep: bool }

fn is_word_byte(c: u8) -> bool {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
}

fn has_name(s: u8[:], w: u8[:]) -> bool {
    var i = 0;
    var prev = 0;
    var prev_end = 0;
    while i < s.len {
        if is_word_byte(s[i]) {
            var j = i;
            while j < s.len && is_word_byte(s[j]) { j++; }
            if j - i == w.len && s[i..j] == w {
                var k = j;
                while k < s.len && s[k] == ' ' { k++; }
                let key = k < s.len && s[k] == ':' && !(k + 1 < s.len && s[k + 1] == ':');
                let after_bc = prev_end - prev == 2 && s[prev..prev_end] == "bc";
                if !key || after_bc { return true; }
            }
            prev = i;
            prev_end = j;
            i = j;
        } else if s[i] == '/' && i + 1 < s.len && s[i + 1] == '/' {
            while i < s.len && s[i] != '\n' { i++; }
        } else {
            i++;
        }
    }
    return false;
}

fn first_word(s: u8[:], from: i64) -> u8[:], i64 {
    var i = from;
    while i < s.len && s[i] == ' ' { i++; }
    var j = i;
    while j < s.len && is_word_byte(s[j]) { j++; }
    return s[i..j], j;
}

fn word_after(s: u8[:], w: u8[:]) -> u8[:] {
    var i = 0;
    while i < s.len {
        if is_word_byte(s[i]) {
            var j = i;
            while j < s.len && is_word_byte(s[j]) { j++; }
            if s[i..j] == w {
                var f = j;
                while f < s.len && (s[f] == ' ' || s[f] == ':') { f++; }
                let next, _ = first_word(s, f);
                return next;
            }
            i = j;
        } else {
            i++;
        }
    }
    return s[0..0];
}

fn items_of(text: u8[:]) -> Item[>..] {
    var out: Item[>..] = [];
    var starts: i64[>..] = [];
    var pending = -1;
    var at = 0;
    while at < text.len {
        var e = at;
        while e < text.len && text[e] != '\n' { e++; }
        let line = text[at..e];
        let top = line.len > 0 && line[0] != ' ' && line[0] != '\t' && line[0] != '}' && line[0] != ')' && line[0] != '\r';
        if top && line.len > 1 && line[0] == '/' && line[1] == '/' {
            if pending < 0 { pending = at; }
        } else if top {
            starts.push(if pending >= 0 { pending } else { at });
            pending = -1;
        }
        at = e + 1;
    }
    for s, i in starts {
        let hi = if i + 1 < starts.len { starts[i + 1] } else { text.len };
        out.push(Item { lo: s, hi: hi, keep: true });
    }
    return out;
}

fn item_head(text: u8[:], it: Item) -> u8[:], u8[:] {
    var at = it.lo;
    while at < it.hi && text[at] == '/' {
        while at < it.hi && text[at] != '\n' { at++; }
        at++;
    }
    let kw, e = first_word(text[..it.hi], at);
    let name, _ = first_word(text[..it.hi], e);
    return kw, name;
}

fn physics_names(text: u8[:], it: Item, names: u8[:][>..]&) {
    let body = text[it.lo..it.hi];
    var s = 0;
    while s + 6 <= body.len && body[s..s + 6] != "solves" { s++; }
    var i = s + 6;
    var depth = 0;
    while i < body.len && body[i] != '{' {
        if body[i] == '<' { depth++; } else if body[i] == '>' { depth--; }
        if is_word_byte(body[i]) && depth == 0 && (i == 0 || !is_word_byte(body[i - 1])) {
            var j = i;
            while j < body.len && is_word_byte(body[j]) { j++; }
            var k = j;
            while k < body.len && body[k] == ' ' { k++; }
            if k < body.len && body[k] == ':' { names.push(body[i..j]); }
            i = j;
        } else {
            i++;
        }
    }
    var at = 0;
    while at < body.len {
        var e = at;
        while e < body.len && body[e] != '\n' { e++; }
        let line = body[at..e];
        let w, after = first_word(line, 0);
        if line.len > 0 && line[0] == ' ' && w == "let" {
            let n, _ = first_word(line, after);
            names.push(n);
        }
        at = e + 1;
    }
}

fn in_names(names: u8[:][:], w: u8[:]) -> bool {
    for n in names { if n == w { return true; } }
    return false;
}

fn mentions_any(s: u8[:], names: u8[:][:]) -> bool {
    for n in names { if n.len > 0 && has_name(s, n) { return true; } }
    return false;
}

fn step_line(line: u8[:], gone: u8[:][:], out: u8[>..]&) {
    let w, after = first_word(line, 0);
    if w == "bc" {
        var c = after;
        while c < line.len && line[c] != ':' { c++; }
        if mentions_any(line[c..], gone) { return; }
    } else if w == "output" {
        var kept: u8[:][>..] = [];
        var i = after;
        while i < line.len {
            var j = i;
            while j < line.len && line[j] != ',' { j++; }
            var lo = i;
            var hi = j;
            while lo < hi && line[lo] == ' ' { lo++; }
            while hi > lo && (line[hi - 1] == ' ' || line[hi - 1] == '\r') { hi--; }
            if hi > lo && !in_names(gone, line[lo..hi]) { kept.push(line[lo..hi]); }
            i = j + 1;
        }
        if kept.len == 0 { return; }
        out.append(line[..after]);
        out.push(' ');
        for k, n in kept {
            if n > 0 { out.append(", "); }
            out.append(k);
        }
        out.push('\n');
        return;
    }
    out.append(line);
    out.push('\n');
}

fn statement_end(text: u8[:], at: i64, hi: i64) -> i64 {
    var indent = 0;
    while at + indent < hi && text[at + indent] == ' ' { indent++; }
    var e = at;
    while e < hi && text[e] != '\n' { e++; }
    while indent > 0 && e + 1 < hi {
        var c = e + 1;
        while c < hi && text[c] == ' ' { c++; }
        guard c < hi && text[c] != '\n' && (c - (e + 1) > indent || text[c] == ')' || text[c] == ']') else { break; };
        e = c;
        while e < hi && text[e] != '\n' { e++; }
    }
    return e;
}

fn item_text(text: u8[:], its: Item[:], phys: i64[:], it: Item, gone: u8[:][:]) -> u8[>..] {
    var out: u8[>..] = [];
    let kw, _ = item_head(text, it);
    if kw != "step" {
        out.append(text[it.lo..it.hi]);
        return out;
    }
    var at = it.lo;
    while at < it.hi {
        let e = statement_end(text, at, it.hi);
        let w, _ = first_word(text[at..e], 0);
        if w != "expect" {
            step_line(text[at..e], gone, out);
        } else if !mentions_any(text[at..e], gone) && owner_count(text, its, phys, text[at..e]) == 1 {
            out.append(text[at..e]);
            out.push('\n');
        }
        at = e + 1;
    }
    return out;
}

fn owner_count(text: u8[:], its: Item[:], phys: i64[:], s: u8[:]) -> i64 {
    var owners = 0;
    for p in phys {
        var names: u8[:][>..] = [];
        physics_names(text, its[p], names);
        if mentions_any(s, names) { owners++; }
    }
    return owners;
}

fn lane_text(text: u8[:], its: Item[>..]&, phys: i64[:], k: i64) -> u8[>..] {
    var kept: u8[:][>..] = [];
    var dropped: u8[:][>..] = [];
    var gone: u8[:][>..] = [];
    var drop_methods: u8[:][>..] = [];
    var drop_meshes: u8[:][>..] = [];
    var keep_mesh: u8[:] = text[0..0];
    for p, j in phys {
        let _, name = item_head(text, its[p]);
        for it in its {
            let kw, mname = item_head(text, it);
            if kw == "method" && mname == name {
                let mesh = word_after(text[it.lo..it.hi], "mesh");
                if j == k { keep_mesh = mesh; } else { drop_meshes.push(mesh); }
                if j != k { drop_methods.push(mname); }
            }
        }
        if j == k { physics_names(text, its[p], kept); } else { physics_names(text, its[p], dropped); }
    }
    for n in dropped { if !in_names(kept, n) { gone.push(n); } }
    for it in its.len {
        let kw, name = item_head(text, its[it]);
        var drop = false;
        for p, j in phys { if p == it && j != k { drop = true; } }
        if kw == "method" && in_names(drop_methods, name) { drop = true; }
        if kw == "mesh" && in_names(drop_meshes, name) && name != keep_mesh { drop = true; }
        if kw == "expect" && (mentions_any(text[its[it].lo..its[it].hi], gone) ||
                              owner_count(text, its, phys, text[its[it].lo..its[it].hi]) != 1) { drop = true; }
        its[it].keep = !drop;
    }
    var changed = true;
    while changed {
        changed = false;
        for it in its.len {
            guard its[it].keep else { continue; };
            let kw, name = item_head(text, its[it]);
            let head = text[its[it].lo..its[it].hi];
            var line_end = 0;
            while line_end < head.len && head[line_end] != '\n' { line_end++; }
            let geometry = kw == "let" && (has_name(head[..line_end], "Solid") || has_name(head[..line_end], "Surface") ||
                                           has_name(head[..line_end], "Curve") || has_name(head[..line_end], "Point"));
            guard geometry || kw == "material" || kw == "table" else { continue; };
            var used = false;
            for o in its.len {
                if o != it && its[o].keep && has_name(item_text(text, its, phys, its[o], gone), name) { used = true; }
            }
            if !used {
                its[it].keep = false;
                changed = true;
            }
        }
    }
    var out: u8[>..] = [];
    for it in its {
        if it.keep { out.append(item_text(text, its, phys, it, gone)); }
    }
    return out;
}

fn lane_texts(text: u8[:], buf: u8[>..]&, ends: i64[>..]&) -> u8[:] {
    var its = items_of(text);
    var phys: i64[>..] = [];
    for it, i in its {
        let kw, _ = item_head(text, it);
        if kw == "physics" { phys.push(i); }
    }
    guard phys.len >= 2 else { return text[0..0]; };
    var bad: u8[:] = text[0..0];
    for it in its {
        let kw, name = item_head(text, it);
        if kw == "expect" && bad.len == 0 && owner_count(text, its, phys, text[it.lo..it.hi]) != 1 { bad = name; }
        guard kw == "step" else { continue; };
        var at = it.lo;
        while at < it.hi && bad.len == 0 {
            let e = statement_end(text, at, it.hi);
            let w, after = first_word(text[at..e], 0);
            if w == "expect" && owner_count(text, its, phys, text[at..e]) != 1 {
                let n, _ = first_word(text[at..e], after);
                bad = n;
            }
            at = e + 1;
        }
    }
    for k in phys.len {
        buf.append(lane_text(text, its, phys, k));
        ends.push(buf.len);
    }
    return bad;
}

fn header_count(line: u8[:]) -> i64 {
    var i = 0;
    while i < line.len && line[i] != '[' { i++; }
    var n = 0;
    i++;
    while i < line.len && line[i] >= '0' && line[i] <= '9' {
        n = n * 10 + (line[i] - '0');
        i++;
    }
    return n;
}

fn line_key(line: u8[:]) -> u8[:] {
    var i = 0;
    while i < line.len && line[i] != '[' && line[i] != ':' { i++; }
    return line[..i];
}

fn header_tail(line: u8[:]) -> u8[:] {
    var i = 0;
    while i < line.len && line[i] != ']' { i++; }
    return line[i..];
}

fn number_after_colon(line: u8[:]) -> i64 {
    var i = 0;
    while i < line.len && line[i] != ':' { i++; }
    i++;
    while i < line.len && line[i] == ' ' { i++; }
    var n = 0;
    while i < line.len && line[i] >= '0' && line[i] <= '9' {
        n = n * 10 + (line[i] - '0');
        i++;
    }
    return n;
}

fn merge_toon(outs: u8[:][:]) -> u8[>..] {
    var keys: u8[:][>..] = [];
    var heads: u8[:][>..] = [];
    var rows: u8[:][>..] = [];
    var row_key: i64[>..] = [];
    var counts: i64[>..] = [];
    var checks = 0;
    var estimate = 0;
    var peak = 0;
    var scalars: u8[:][>..] = [];
    for out in outs {
        var cur = -1;
        var at = 0;
        while at < out.len {
            var e = at;
            while e < out.len && out[e] != '\n' { e++; }
            let line = out[at..e];
            at = e + 1;
            if line.len == 0 { continue; }
            if line[0] == ' ' {
                if cur >= 0 {
                    rows.push(line);
                    row_key.push(cur);
                }
                continue;
            }
            cur = -1;
            let key = line_key(line);
            if key == "memory_checks" { checks += number_after_colon(line); continue; }
            if key == "estimate_mb" { estimate = max(estimate, number_after_colon(line)); continue; }
            if key == "peak_mb" { peak = max(peak, number_after_colon(line)); continue; }
            var is_table = false;
            for c in line { if c == '[' { is_table = true; } }
            if !is_table {
                scalars.push(line);
                continue;
            }
            for kk, j in keys {
                if kk == key && header_tail(heads[j]) == header_tail(line) { cur = j; }
            }
            if cur < 0 {
                keys.push(key);
                heads.push(line);
                counts.push(0);
                cur = keys.len - 1;
            }
            counts[cur] += header_count(line);
        }
    }
    var merged: u8[>..] = [];
    for key, j in keys {
        let head = heads[j];
        format(merged, key, "[", counts[j], header_tail(head), "\n");
        for r, ri in rows { if row_key[ri] == j { format(merged, r, "\n"); } }
    }
    for s in scalars { format(merged, s, "\n"); }
    format(merged, "memory_checks: ", checks, "\n");
    format(merged, "estimate_mb: ", estimate, "\n");
    format(merged, "peak_mb: ", peak, "\n");
    return merged;
}
