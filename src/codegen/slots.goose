// The slots of the program of the page (web.goose): the number literals
// of a Goose program, and the values of the literals that change with
// the params.
import std;
namespace codegen;

// A number literal of a program: its bytes lo..hi, and whether it is an
// f64 (with a point or an exponent).
struct Literal {
    lo: i64,
    hi: i64,
    float: bool,
}

fn is_ident_byte(c: u8) -> bool { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'; }
fn is_digit(c: u8) -> bool { return c >= '0' && c <= '9'; }

// The number literals of the Goose program `text`, out of comments,
// strings and characters. A `-` just before a digit is part of the literal
// when no name or closing bracket is in front of it: the generator writes
// a negative value so, and a subtraction with spaces.
fn literals(text: u8[:]) -> Literal[>..] {
    var out: Literal[>..] = [];
    var i = 0;
    while i < text.len {
        let c = text[i];
        if c == '/' && i + 1 < text.len && text[i + 1] == '/' {
            while i < text.len && text[i] != '\n' { i++; }
        } else if c == '"' || c == '\'' {
            i++;
            while i < text.len && text[i] != c {
                if text[i] == '\\' { i++; }
                i++;
            }
            i++;
        } else if is_ident_byte(c) && !is_digit(c) {
            while i < text.len && is_ident_byte(text[i]) { i++; }
        } else if is_digit(c) || (c == '-' && i + 1 < text.len && is_digit(text[i + 1]) &&
                                  !(i > 0 && (is_ident_byte(text[i - 1]) || text[i - 1] == ')' || text[i - 1] == ']'))) {
            var l = Literal { lo: i, hi: i + 1, float: false };
            while l.hi < text.len {
                let d = text[l.hi];
                let sign = (d == '-' || d == '+') && (text[l.hi - 1] == 'e' || text[l.hi - 1] == 'E');
                guard is_ident_byte(d) || d == '.' || sign else { break; };
                if d == '.' || d == 'e' || d == 'E' { l.float = true; }
                l.hi++;
            }
            out.push(l);
            i = l.hi;
        } else {
            i++;
        }
    }
    return out;
}

// Whether programs `a` and `b`, with the literals `la` and `lb`, differ
// only in the literals.
fn same_text(a: u8[:], la: Literal[:], b: u8[:], lb: Literal[:]) -> bool {
    guard la.len == lb.len else { return false; };
    var at_a = 0;
    var at_b = 0;
    for l, k in la {
        guard a[at_a..l.lo] == b[at_b..lb[k].lo] else { return false; };
        at_a = l.hi;
        at_b = lb[k].hi;
    }
    return a[at_a..] == b[at_b..];
}

// The values of the slots `slots` (indices into the literals of `base`)
// in the program `text`, into `f` and `i` in the sequence of the slots.
// False, with a message on stderr, if `text` differs from `base` in
// another place than a slot.
fn slot_values(base: u8[:], slots: i64[:], text: u8[:], f: f64[>..]&, i: i64[>..]&) -> bool {
    let lb = literals(base);
    let lt = literals(text);
    guard same_text(base, lb, text, lt) else {
        write_stderr("skein web: these values change the structure of the program: run the model with skein run --set\n");
        return false;
    };
    f.clear();
    i.clear();
    var s = 0;
    for l, k in lt {
        let word = text[l.lo..l.hi];
        if s < slots.len && slots[s] == k {
            s++;
            if lb[k].float {
                let v, good = parse_flt(word);
                guard good else { return false; };
                f.push(v);
            } else {
                let v, good = parse_int(word);
                guard good else { return false; };
                i.push(v);
            }
        } else if word != base[lb[k].lo..lb[k].hi] {
            write_stderr(str("skein web: these values change the literal ", base[lb[k].lo..lb[k].hi],
                             " that the module holds fixed: run the model with skein run --set\n"));
            return false;
        }
    }
    return true;
}

// The program of the page: `text` with each literal of `slots` replaced
// by its slot.
fn slot_program(text: u8[:], slots: i64[:]) -> u8[>..] {
    let ls = literals(text);
    var out: u8[>..] = [];
    var at = 0;
    var nf = 0;
    var ni = 0;
    for k in slots {
        let l = ls[k];
        out.append(text[at..l.lo]);
        if l.float {
            format(out, "web_f[", nf, "]");
            nf++;
        } else {
            format(out, "web_i[", ni, "]");
            ni++;
        }
        at = l.hi;
    }
    out.append(text[at..]);
    return out;
}

// The slots of the text `s`, indices with a comma between them.
fn parse_slots(s: u8[:]) -> i64[>..] {
    var out: i64[>..] = [];
    var start = 0;
    for i in 0..s.len + 1 {
        if i == s.len || s[i] == ',' {
            let v, good = parse_int(s[start..i]);
            if good { out.push(v); }
            start = i + 1;
        }
    }
    return out;
}
