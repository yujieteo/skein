import std;
namespace syntax;
import .tree;

fn is_letter(c: u8) -> bool { (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' }
fn is_digit(c: u8) -> bool { c >= '0' && c <= '9' }

fn is_type_constructor(word: u8[:]) -> bool {
    word == "Vector" || word == "Tensor2" || word == "Tensor4" || word == "Complex"
}

fn continues(kind: i64) -> bool {
    match kind {
        TK_COMMA, TK_COLON, TK_ASSIGN, TK_LPAREN, TK_LBRACKET, TK_LBRACE, TK_LANGLE => true,
        TK_EQ, TK_NE, TK_LT, TK_LE, TK_GT, TK_GE, TK_PLUS, TK_MINUS, TK_STAR, TK_SLASH => true,
        TK_CARET, TK_OROR, TK_ANDAND, TK_BAR, TK_AMP, TK_ARROW, TK_DOTDOT => true,
        _ => false,
    }
}

fn starts_item_line(text: u8[:], i: i64) -> bool {
    var k = i;
    while k < text.len && (is_letter(text[k]) || is_digit(text[k])) { k++; }
    if k == i || !starts_item(keyword_kind(text[i..k])) { return false; }
    while k < text.len && (text[k] == ' ' || text[k] == '\t') { k++; }
    return k >= text.len || text[k] != ':';
}

fn punct(text: u8[:], i: i64) -> i64, i64 {
    let c = text[i];
    let d = if i + 1 < text.len { text[i + 1] } else { 0 as u8 };
    if c == '.' && d == '.' { return TK_DOTDOT, 2; }
    if c == '-' && d == '>' { return TK_ARROW, 2; }
    if c == '=' && d == '=' { return TK_EQ, 2; }
    if c == '!' && d == '=' { return TK_NE, 2; }
    if c == '<' && d == '=' { return TK_LE, 2; }
    if c == '>' && d == '=' { return TK_GE, 2; }
    if c == '|' && d == '|' { return TK_OROR, 2; }
    if c == '&' && d == '&' { return TK_ANDAND, 2; }
    let k = match c {
        '(' => TK_LPAREN, ')' => TK_RPAREN, '[' => TK_LBRACKET, ']' => TK_RBRACKET,
        '{' => TK_LBRACE, '}' => TK_RBRACE, ',' => TK_COMMA, ':' => TK_COLON,
        '.' => TK_DOT, '=' => TK_ASSIGN, '<' => TK_LT, '>' => TK_GT, '+' => TK_PLUS,
        '-' => TK_MINUS, '*' => TK_STAR, '/' => TK_SLASH, '^' => TK_CARET, '!' => TK_BANG,
        '|' => TK_BAR, '&' => TK_AMP,
        _ => TK_ERROR,
    };
    return k, 1;
}

fn lex(f: i64) {
    let base = files[f].start;
    let text = source[base..base + files[f].len];
    var i = 0;
    var line = 1;
    var line_start = 0;
    var depth = 0;
    var angle = 0;
    var last = TK_NL;
    var stmt_first = -1;
    var brace = 0;
    var in_table = false;
    var at_line_start = true;

    fn col_of(at: i64) -> i64 {
        var c = 1;
        for k in line_start..at {
            if (text[k] & 0xC0) != 0x80 { c++; }
        }
        return c;
    }

    fn emit(kind: i64, at: i64, len: i64) {
        tokens.push(Token { kind: kind, start: base + at, len: len, line: line, col: col_of(at), file: f });
        if kind == TK_COMMENT || kind == TK_DOC_COMMENT { return; }
        if last == TK_NL || last == TK_INDENT { stmt_first = kind; }
        last = kind;
    }

    fn newline(at: i64) {
        if starts_item_line(text, at + 1) {
            if last != TK_NL && last != TK_INDENT { emit(TK_NL, at, 0); }
            depth = 0;
            angle = 0;
            brace = 0;
            in_table = false;
        } else if in_table {
            if last != TK_NL { emit(TK_NL, at, 0); }
        } else if depth == 0 && last != TK_NL && last != TK_INDENT {
            if last == TK_COLON && stmt_first == KW_TABLE && brace == 0 {
                emit(TK_NL, at, 0);
                in_table = true;
            } else if !continues(last) {
                emit(TK_NL, at, 0);
            }
        }
        line++;
        line_start = at + 1;
        at_line_start = true;
    }

    while i < text.len {
        let c = text[i];
        if at_line_start {
            at_line_start = false;
            if in_table {
                var k = i;
                while k < text.len && (text[k] == ' ' || text[k] == '\t' || text[k] == '\r') { k++; }
                let blank = k >= text.len || text[k] == '\n' ||
                    (text[k] == '/' && k + 1 < text.len && text[k + 1] == '/');
                if !blank {
                    if k > i { emit(TK_INDENT, k, 0); i = k; continue; }
                    in_table = false;
                }
            }
        }
        if c == '\n' { newline(i); i++; continue; }
        if c == ' ' || c == '\t' || c == '\r' { i++; continue; }
        let start = i;
        if c == '/' && i + 1 < text.len && text[i + 1] == '/' {
            while i < text.len && text[i] != '\n' { i++; }
            var end = i;
            while end > start && text[end - 1] == '\r' { end--; }
            let doc = end - start >= 3 && text[start + 2] == '/';
            emit(if doc { TK_DOC_COMMENT } else { TK_COMMENT }, start, end - start);
            continue;
        }
        if is_letter(c) {
            while i < text.len && (is_letter(text[i]) || is_digit(text[i])) { i++; }
            let word = text[start..i];
            emit(keyword_kind(word), start, i - start);
            continue;
        }
        if is_digit(c) {
            while i < text.len && is_digit(text[i]) { i++; }
            if i - start == 1 && c >= '1' && c <= '3' && i < text.len && text[i] == 'd' &&
                    (i + 1 >= text.len || !(is_letter(text[i + 1]) || is_digit(text[i + 1]))) {
                i++;
                emit(TK_SPACE, start, 2);
                continue;
            }
            if i + 1 < text.len && text[i] == '.' && is_digit(text[i + 1]) {
                i++;
                while i < text.len && is_digit(text[i]) { i++; }
            }
            if i < text.len && (text[i] == 'e' || text[i] == 'E') {
                var k = i + 1;
                if k < text.len && (text[k] == '+' || text[k] == '-') { k++; }
                if k < text.len && is_digit(text[k]) {
                    i = k;
                    while i < text.len && is_digit(text[i]) { i++; }
                }
            }
            if i < text.len && text[i] == '%' {
                i++;
                emit(TK_PERCENT, start, i - start);
            } else {
                emit(TK_NUMBER, start, i - start);
            }
            continue;
        }
        if c == '"' {
            i++;
            var closed = false;
            while i < text.len && text[i] != '\n' {
                if text[i] == '\\' && i + 1 < text.len && (text[i + 1] == '"' || text[i + 1] == '\\') {
                    i += 2;
                    continue;
                }
                if text[i] == '"' { i++; closed = true; break; }
                i++;
            }
            if closed {
                emit(TK_STRING, start, i - start);
            } else {
                emit(TK_ERROR, start, i - start);
                error_at(tokens.len - 1, "syntax", "the string has no closing quote",
                         "add '\"' at the end of the string");
            }
            continue;
        }
        var kind, len = punct(text, i);
        if kind == TK_LT && last == TK_IDENT && is_type_constructor(syntax::text(tokens.len - 1)) {
            kind = TK_LANGLE;
        } else if kind == TK_GT && angle > 0 {
            kind = TK_RANGLE;
        }
        if kind == TK_ERROR {
            i++;
            while i < text.len && (text[i] & 0xC0) == 0x80 { i++; }
            emit(TK_ERROR, start, i - start);
            error_at(tokens.len - 1, "syntax", str("the character '", text[start..i], "' is not valid here"),
                     "remove it, or put it in a comment or a string");
            continue;
        }
        if kind == TK_LPAREN || kind == TK_LBRACKET || kind == TK_LANGLE { depth++; }
        if (kind == TK_RPAREN || kind == TK_RBRACKET || kind == TK_RANGLE) && depth > 0 { depth--; }
        if kind == TK_LANGLE { angle++; }
        if kind == TK_RANGLE { angle--; }
        if kind == TK_LBRACE { brace++; }
        if kind == TK_RBRACE && brace > 0 { brace--; }
        emit(kind, start, len);
        i += len;
    }
    if last != TK_NL { emit(TK_NL, text.len, 0); }
    emit(TK_EOF, text.len, 0);
}
