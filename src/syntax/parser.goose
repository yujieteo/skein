// The parser (spec §3). It reads the tokens of one file and writes the
// syntax tree into `nodes`. It is recursive descent, with one function for
// each production of the grammar.
//
// After a syntax error, the parser records the error and skips to the next
// line that starts at column 1 with a keyword of an item, and continues
// there. Thus one run reports the errors of each item.
import std;
namespace syntax;
import .tree;
import .units;

// The tokens of the file, without the comments, as indices into `tokens`.
var ptoks: i64[>..];
// The position in `ptoks`.
var pos = 0;

fn tok() -> i64 { ptoks[pos] }
fn kind() -> i64 { tokens[ptoks[pos]].kind }
fn kind_at(k: i64) -> i64 { tokens[ptoks[min(pos + k, ptoks.len - 1)]].kind }
fn at(k: i64) -> bool { kind() == k }

// The current token, and moves past it.
fn advance() -> i64 {
    let t = ptoks[pos];
    if pos < ptoks.len - 1 { pos++; }
    return t;
}

fn accept(k: i64) -> bool {
    if kind() != k { return false; }
    pos++;
    return true;
}

// Words for a token in a message.
fn describe(t: i64) -> u8[>..] {
    let k = tokens[t].kind;
    if k == TK_NL { return str("the end of the line"); }
    if k == TK_EOF { return str("the end of the file"); }
    if k == TK_INDENT { return str("an indented line"); }
    return str("'", text(t), "'");
}

// Records a syntax error at the current token, and leaves the item.
// The lexer has already reported an error at a TK_ERROR token.
fn fail(msg: u8[:], fix: u8[:]) -> i64 {
    if kind() != TK_ERROR { error_at(tok(), "syntax", msg, fix); }
    return -1 from parse_item;
}

fn fail_expected(what: u8[:], fix: u8[:]) -> i64 {
    return fail(str("expected ", what, ", found ", describe(tok())), fix);
}

// The current token, which must be of kind `k`, and moves past it.
fn expect(k: i64, what: u8[:]) -> i64 {
    if kind() != k { fail_expected(what, str("insert ", what)); }
    return advance();
}

fn expect_ident(what: u8[:]) -> i64 {
    if kind() != TK_IDENT {
        if kind() >= KW_ADAPT {
            fail(str("expected ", what, ", found the keyword '", text(tok()), "'"),
                 str("use a name that is not a keyword (spec §2.5)"));
        }
        fail_expected(what, str("insert ", what));
    }
    return advance();
}

// True when the cursor is at a name and `:`. The name of a setting or of an
// argument may be a keyword, as in `mesh: m` and `results(f, step: "load")`.
fn at_label() -> bool {
    (kind() == TK_IDENT || kind() >= KW_ADAPT) && kind_at(1) == TK_COLON
}

// The name of a setting, an argument or a field, and its `:`.
fn label(what: u8[:]) -> i64 {
    if !at_label() { expect_ident(what); expect(TK_COLON, "':'"); }
    let t = advance();
    advance();
    return t;
}

// The end of a statement.
fn end_statement() {
    if at(TK_EOF) { return; }
    if kind() != TK_NL {
        fail(str("expected the end of the statement, found ", describe(tok())),
             "put the rest on a new line, or remove it");
    }
    advance();
}

fn skip_newlines() { while at(TK_NL) { advance(); } }

// A node of kind `k` at token `t` with the children `a` and `b`.
fn node2(k: i64, t: i64, a: i64, b: i64) -> i64 {
    let n = node(k, t);
    add(n, a);
    add(n, b);
    return n;
}

fn node1(k: i64, t: i64, a: i64) -> i64 {
    let n = node(k, t);
    add(n, a);
    return n;
}

// --- units (spec §6.3) ---------------------------------------------------------

// True when the tokens at the cursor continue a unit expression after `*` or
// `/`: a unit name, or `(` and a unit name.
fn unit_follows() -> bool {
    if at(TK_IDENT) { return is_unit_name(text(tok())); }
    return at(TK_LPAREN) && kind_at(1) == TK_IDENT && is_unit_name(text(ptoks[pos + 1]));
}

// The exponent after `^`: an optional `-` and digits.
fn exponent() -> i64 {
    let neg = accept(TK_MINUS);
    let t = expect(TK_NUMBER, "an integer exponent");
    let v, ok = parse_int(text(t));
    if !ok { fail(str("the exponent '", text(t), "' is not an integer"), "write a whole number"); }
    return if neg { -v } else { v };
}

recursive fn unit_factor() -> i64 {
    if at(TK_LPAREN) {
        advance();
        let u = unit_product(true);
        expect(TK_RPAREN, "')'");
        return u;
    }
    let t = expect_ident("a unit name");
    var u = node(N_UNIT_NAME, t);
    if at(TK_CARET) {
        let c = advance();
        let e = exponent();
        u = node1(N_UNIT_POW, c, u);
        nodes[u].aux = e;
    }
    return u;
}

// A unit expression. Inside parentheses each `*` and `/` continues it. After
// a literal, it continues only before a unit name (spec §6.3).
recursive fn unit_product(inside: bool) -> i64 {
    var u = unit_factor();
    while at(TK_STAR) || at(TK_SLASH) {
        if !inside {
            let save = pos;
            pos++;
            let more = unit_follows();
            pos = save;
            if !more { break; }
        }
        let op = advance();
        let r = unit_factor();
        u = node2(if tokens[op].kind == TK_STAR { N_UNIT_MUL } else { N_UNIT_DIV }, op, u, r);
    }
    return u;
}

// An N_UNIT node for the unit expression at the cursor.
fn unit() -> i64 {
    let t = tok();
    return node1(N_UNIT, t, unit_product(false));
}

// Adds the unit after a literal to `n`, if there is one (spec §2.6).
fn literal_unit(n: i64) {
    if at(TK_IDENT) { add(n, unit()); }
}

// --- types (spec §5) ---------------------------------------------------------------

recursive fn dim_term() -> i64 {
    if at(TK_LPAREN) {
        advance();
        let d = dim_product();
        expect(TK_RPAREN, "')'");
        return d;
    }
    let t = expect_ident("a type");
    var d = node(N_TYPE_NAME, t);
    if at(TK_CARET) {
        let c = advance();
        let e = exponent();
        d = node1(N_DIM_POW, c, d);
        nodes[d].aux = e;
    }
    return d;
}

recursive fn dim_product() -> i64 {
    var d = dim_term();
    while at(TK_STAR) || at(TK_SLASH) {
        let op = advance();
        let r = dim_term();
        d = node2(if tokens[op].kind == TK_STAR { N_DIM_MUL } else { N_DIM_DIV }, op, d, r);
    }
    return d;
}

recursive fn parse_type() -> i64 {
    if at(TK_IDENT) && kind_at(1) == TK_LANGLE {
        let t = advance();
        advance();
        let inner = parse_type();
        expect(TK_RANGLE, "'>'");
        return node1(N_TYPE_APPLY, t, inner);
    }
    return dim_product();
}

// --- expressions (spec §3, §7) ---------------------------------------------------

// An argument list after `(`, up to and with the `)`. Each argument is added
// to `n` as an N_ARG. With `lists`, an argument can be a list (a procedure).
recursive fn args(n: i64, lists: bool) {
    while !at(TK_RPAREN) {
        let start = tok();
        let a = node(N_ARG, start);
        if at_label() { nodes[a].aux = label("a name"); }
        add(a, if lists && at(TK_LBRACKET) { list() } else { expr() });
        add(n, a);
        if !accept(TK_COMMA) { break; }
    }
    expect(TK_RPAREN, "')' or ','");
}

// The items of a vector, a list or a record, after the open bracket, up to
// and with the close bracket of kind `close`.
recursive fn items(n: i64, close: i64, what: u8[:]) {
    while !at(close) {
        add(n, expr());
        if !accept(TK_COMMA) { break; }
    }
    expect(close, what);
}

// A list `[a, b]` and its unit. It is a study entry or a procedure
// argument (spec §3, `list`).
recursive fn list() -> i64 {
    let l = node(N_LIST, advance());
    items(l, TK_RBRACKET, "']' or ','");
    literal_unit(l);
    return l;
}

recursive fn primary() -> i64 {
    let t = tok();
    match kind() {
        TK_NUMBER => {
            advance();
            let n = node(N_NUMBER, t);
            literal_unit(n);
            return n;
        },
        TK_PERCENT => { advance(); return node(N_PERCENT, t); },
        TK_STRING => { advance(); return node(N_STRING, t); },
        KW_TRUE, KW_FALSE => { advance(); return node(N_BOOL, t); },
        TK_IDENT => { advance(); return node(N_NAME, t); },
        TK_LPAREN => {
            advance();
            let first = expr();
            if !at(TK_COMMA) {
                expect(TK_RPAREN, "')'");
                return first;
            }
            let v = node(N_VECTOR, t);
            add(v, first);
            advance();
            items(v, TK_RPAREN, "')' or ','");
            literal_unit(v);
            return v;
        },
        TK_LBRACE => {
            advance();
            let r = node(N_RECORD, t);
            accept(TK_NL);
            while !at(TK_RBRACE) {
                let name = label("a field name");
                let a = node(N_ARG, name);
                nodes[a].aux = name;
                add(a, expr());
                add(r, a);
                if !accept(TK_COMMA) && !accept(TK_NL) { break; }
                accept(TK_NL);
            }
            expect(TK_RBRACE, "'}'");
            return r;
        },
        _ => {
            return fail_expected("a value", "write a number, a name or an expression here");
        },
    }
}

recursive fn postfix() -> i64 {
    var e = primary();
    loop {
        if at(TK_DOT) {
            advance();
            let name = expect_ident("a member name after '.'");
            e = node1(N_MEMBER, name, e);
        } else if at(TK_LPAREN) {
            let t = advance();
            let c = node1(N_CALL, t, e);
            args(c, false);
            e = c;
        } else if at(KW_AT) {
            let t = advance();
            e = node2(N_AT, t, e, primary());
        } else {
            break;
        }
    }
    return e;
}

recursive fn power() -> i64 {
    let b = postfix();
    if !at(TK_CARET) { return b; }
    let t = advance();
    return node2(N_BINARY, t, b, unary());
}

recursive fn unary() -> i64 {
    if at(TK_MINUS) || at(TK_BANG) {
        let t = advance();
        let e = unary();
        // `-5 degC` is one literal, -5 °C, and not -(5 degC).
        if tokens[t].kind == TK_MINUS && nodes[e].kind == N_NUMBER && nodes[e].aux < 0 {
            let u = child_of_kind(e, N_UNIT);
            if u >= 0 && nodes[nodes[u].first].kind == N_UNIT_NAME && text(nodes[u].tok) == "degC" {
                nodes[e].aux = t;
                return e;
            }
        }
        return node1(N_UNARY, t, e);
    }
    return power();
}

recursive fn mul_expr() -> i64 {
    var e = unary();
    while at(TK_STAR) || at(TK_SLASH) {
        let t = advance();
        e = node2(N_BINARY, t, e, unary());
    }
    return e;
}

recursive fn add_expr() -> i64 {
    var e = mul_expr();
    while at(TK_PLUS) || at(TK_MINUS) {
        let t = advance();
        e = node2(N_BINARY, t, e, mul_expr());
    }
    return e;
}

fn is_comparison(k: i64) -> bool {
    k == TK_EQ || k == TK_NE || k == TK_LT || k == TK_LE || k == TK_GT || k == TK_GE
}

recursive fn cmp_expr() -> i64 {
    let e = add_expr();
    if !is_comparison(kind()) { return e; }
    let t = advance();
    let c = node2(N_BINARY, t, e, add_expr());
    if is_comparison(kind()) {
        fail("a comparison does not chain (spec §7.1)", "join the two comparisons with '&&'");
    }
    return c;
}

recursive fn and_expr() -> i64 {
    var e = cmp_expr();
    while at(TK_ANDAND) || at(TK_AMP) {
        let t = advance();
        e = node2(N_BINARY, t, e, cmp_expr());
    }
    return e;
}

recursive fn expr() -> i64 {
    var e = and_expr();
    while at(TK_OROR) || at(TK_BAR) {
        let t = advance();
        e = node2(N_BINARY, t, e, and_expr());
    }
    return e;
}

// --- blocks and procedures ---------------------------------------------------------

// The value of a setting: an expression, a range or a study list (spec §3,
// `study`).
fn setting_value() -> i64 {
    if at(TK_LBRACKET) { return list(); }
    let v = expr();
    if !at(TK_DOTDOT) { return v; }
    let t = advance();
    return node2(N_RANGE, t, v, expr());
}

// A block: `{`, entries on separate lines or with commas, `}`.
fn parse_block() -> i64 {
    let b = node(N_BLOCK, expect(TK_LBRACE, "'{'"));
    accept(TK_NL);
    while !at(TK_RBRACE) {
        if at(KW_REFINE) || at(KW_LAYERS) {
            let t = advance();
            let region = expr();
            expect(TK_COLON, "':'");
            add(b, node2(if tokens[t].kind == KW_REFINE { N_REFINE } else { N_LAYERS }, t, region, expr()));
        } else {
            let name = label("a setting name");
            add(b, node1(N_SETTING, name, setting_value()));
        }
        if !accept(TK_COMMA) && !at(TK_NL) { break; }
        accept(TK_NL);
    }
    expect(TK_RBRACE, "'}'");
    return b;
}

// A procedure: a name and an optional argument list.
fn procedure() -> i64 {
    let p = node(N_PROC, expect_ident("a procedure name"));
    if accept(TK_LPAREN) { args(p, true); }
    return p;
}

// `name: type`, the start of many declarations. It gives the node of kind
// `k` with the type as its child.
fn typed_name(k: i64, what: u8[:]) -> i64 {
    let n = node(k, expect_ident(what));
    expect(TK_COLON, "':'");
    add(n, parse_type());
    return n;
}

// --- items (spec §3) -----------------------------------------------------------

fn header() -> i64 {
    advance();
    let h = node(N_HEADER, expect_ident("the model name"));
    expect(TK_COLON, "':'");
    if at(TK_SPACE) || (at(TK_IDENT) && text(tok()) == "axi") {
        add(h, node(N_SPACE, advance()));
    } else {
        fail_expected("the space '1d', '2d', '3d' or 'axi'", "insert the space, for example '3d'");
    }
    expect(TK_COMMA, "','");
    let us = node(N_UNITSYS, tok());
    if accept(TK_LBRACE) {
        loop {
            let f = node(N_UNITFIELD, expect_ident("a unit field name"));
            expect(TK_COLON, "':'");
            add(f, unit());
            add(us, f);
            if !accept(TK_COMMA) { break; }
        }
        expect(TK_RBRACE, "'}'");
    } else {
        expect_ident("a unit system, for example 'SI'");
    }
    add(h, us);
    end_statement();
    return h;
}

fn param_decl() -> i64 {
    advance();
    let n = typed_name(N_PARAM, "a parameter name");
    expect(TK_ASSIGN, "'='");
    add(n, expr());
    if at(KW_TOL) {
        let t = advance();
        add(n, node1(N_TOL, t, expr()));
    }
    if accept(KW_IN) {
        let lo = expr();
        let t = expect(TK_DOTDOT, "'..'");
        add(n, node2(N_RANGE, t, lo, expr()));
    }
    end_statement();
    return n;
}

fn let_decl() -> i64 {
    advance();
    let n = typed_name(N_LET, "a name");
    expect(TK_ASSIGN, "'='");
    add(n, expr());
    end_statement();
    return n;
}

fn fn_decl() -> i64 {
    advance();
    let n = node(N_FN, expect_ident("a function name"));
    expect(TK_LPAREN, "'('");
    loop {
        add(n, typed_name(N_FNPARAM, "a parameter name"));
        if !accept(TK_COMMA) || at(TK_RPAREN) { break; }
    }
    expect(TK_RPAREN, "')' or ','");
    expect(TK_ARROW, "'->'");
    add(n, parse_type());
    expect(TK_ASSIGN, "'='");
    add(n, expr());
    end_statement();
    return n;
}

fn table_decl() -> i64 {
    advance();
    let n = node(N_TABLE, expect_ident("a table name"));
    expect(TK_LBRACE, "'{'");
    loop {
        let c = typed_name(N_COLUMN, "a column name");
        if accept(KW_IN) { add(c, unit()); }
        add(n, c);
        if !accept(TK_COMMA) { break; }
    }
    expect(TK_RBRACE, "'}'");
    expect(TK_COLON, "':'");
    end_statement();
    while at(TK_INDENT) {
        advance();
        let r = node(N_ROW, tok());
        loop {
            add(r, expr());
            if !accept(TK_COMMA) { break; }
        }
        add(n, r);
        end_statement();
    }
    return n;
}

fn material_decl() -> i64 {
    advance();
    let n = node(N_MATERIAL, expect_ident("a material name"));
    expect(TK_LBRACE, "'{'");
    accept(TK_NL);
    while !at(TK_RBRACE) && !at(TK_EOF) {
        let p = typed_name(N_PROP, "a property name");
        expect(TK_ASSIGN, "'='");
        add(p, expr());
        add(n, p);
        end_statement();
    }
    expect(TK_RBRACE, "'}'");
    end_statement();
    return n;
}

// `mesh` or `points`.
fn discretization(k: i64) -> i64 {
    advance();
    let n = node(k, expect_ident("a name"));
    expect(KW_ON, "'on'");
    add(n, expr());
    add(n, parse_block());
    end_statement();
    return n;
}

fn unknowns(n: i64) {
    loop {
        let u = typed_name(N_UNKNOWN, "an unknown");
        if accept(TK_ASSIGN) { add(u, expr()); }
        add(n, u);
        if !accept(TK_COMMA) { break; }
    }
}

// One statement of a physics block or a template (spec §3, `physstmt`).
fn physics_statement() -> i64 {
    let t = tok();
    match kind() {
        KW_MATERIAL => {
            advance();
            let m = node(N_MATASSIGN, t);
            loop {
                let one = node(N_MATON, expect_ident("a material name"));
                if accept(KW_ON) { add(one, expr()); }
                add(m, one);
                if !accept(TK_COMMA) { break; }
            }
            return m;
        },
        KW_ASSUME => {
            advance();
            let a = node(N_ASSUME, t);
            loop {
                let one = node(N_ASSUMPTION, expect_ident("an assumption"));
                if accept(TK_LPAREN) {
                    loop {
                        add(one, expr());
                        if !accept(TK_COMMA) { break; }
                    }
                    expect(TK_RPAREN, "')' or ','");
                }
                add(a, one);
                if !accept(TK_COMMA) { break; }
            }
            return a;
        },
        KW_LET => {
            advance();
            let n = typed_name(N_LET, "a name");
            expect(TK_ASSIGN, "'='");
            add(n, expr());
            return n;
        },
        KW_EQ => {
            advance();
            let l = expr();
            expect(TK_ASSIGN, "'='");
            return node2(N_EQN, t, l, expr());
        },
        KW_STATE => {
            advance();
            let n = typed_name(N_STATE, "a state variable");
            expect(TK_ASSIGN, "'='");
            add(n, expr());
            return n;
        },
        KW_PLASTIC => {
            advance();
            return node1(N_PLASTIC, t, parse_block());
        },
        _ => {
            return fail_expected("a statement of a physics block",
                                 "start the line with material, assume, let, eq, state or plastic");
        },
    }
}

// `{`, the statements of a physics block, `}`.
fn physics_body(n: i64) {
    expect(TK_LBRACE, "'{'");
    accept(TK_NL);
    while !at(TK_RBRACE) && !at(TK_EOF) {
        add(n, physics_statement());
        end_statement();
    }
    expect(TK_RBRACE, "'}'");
}

fn physics_decl() -> i64 {
    advance();
    let n = node(N_PHYSICS, expect_ident("a physics name"));
    if at(KW_ON) {
        let t = advance();
        add(n, node1(N_ON, t, expr()));
    }
    if accept(TK_COLON) {
        let u = node(N_USES, expect_ident("a template name"));
        if accept(TK_LPAREN) { args(u, false); }
        add(n, u);
    }
    expect(KW_SOLVES, "'solves'");
    unknowns(n);
    if at(TK_LBRACE) { physics_body(n); }
    end_statement();
    return n;
}

fn template_decl() -> i64 {
    advance();
    expect(KW_PHYSICS, "'physics'");
    let n = node(N_TEMPLATE, expect_ident("a template name"));
    expect(TK_LPAREN, "'('");
    while !at(TK_RPAREN) {
        let p = typed_name(N_TPARAM, "a parameter name");
        if accept(TK_ASSIGN) { add(p, expr()); }
        add(n, p);
        if !accept(TK_COMMA) { break; }
    }
    expect(TK_RPAREN, "')' or ','");
    expect(KW_SOLVES, "'solves'");
    unknowns(n);
    physics_body(n);
    end_statement();
    return n;
}

fn method_decl() -> i64 {
    advance();
    let n = node(N_METHOD, expect_ident("the name of a physics"));
    expect(TK_COLON, "':'");
    nodes[n].aux = expect_ident("a method kind, for example 'cg'");
    add(n, parse_block());
    end_statement();
    return n;
}

// `bc region: target = value`, without the end of the statement.
fn bc_statement() -> i64 {
    let n = node(N_BC, advance());
    add(n, expr());
    expect(TK_COLON, "':'");
    add(n, expr());
    expect(TK_ASSIGN, "'='");
    add(n, expr());
    return n;
}

// `contact` or `interface`.
fn pair_decl(k: i64) -> i64 {
    advance();
    let n = node(k, expect_ident("a name"));
    expect(TK_COLON, "':'");
    add(n, expr());
    expect(TK_COMMA, "','");
    add(n, expr());
    add(n, parse_block());
    end_statement();
    return n;
}

// `symmetry` or `study`: a name, a procedure and a block.
fn procedure_decl(k: i64) -> i64 {
    advance();
    let n = node(k, expect_ident("a name"));
    expect(TK_COLON, "':'");
    add(n, procedure());
    add(n, parse_block());
    end_statement();
    return n;
}

// `birth` or `death`.
fn active_set(k: i64) -> i64 {
    let n = node(k, advance());
    loop {
        add(n, expr());
        if !accept(TK_COMMA) { break; }
    }
    if at(KW_WHERE) {
        let t = advance();
        add(n, node1(N_WHERE, t, expr()));
    }
    return n;
}

// One statement of a step (spec §3, `stepstmt`).
fn step_statement() -> i64 {
    let t = tok();
    match kind() {
        KW_BC => { return bc_statement(); },
        KW_RELEASE => { advance(); return node1(N_RELEASE, t, expr()); },
        KW_BIRTH => { return active_set(N_BIRTH); },
        KW_DEATH => { return active_set(N_DEATH); },
        KW_SOLVE => { advance(); return node1(N_SOLVE, t, parse_block()); },
        KW_ADAPT => { advance(); return node1(N_ADAPT, t, parse_block()); },
        KW_OUTPUT => {
            advance();
            let n = node(N_OUTPUT, t);
            loop {
                add(n, expr());
                if !accept(TK_COMMA) { break; }
            }
            if at(KW_EVERY) {
                let e = advance();
                add(n, node1(N_EVERY, e, expr()));
            }
            return n;
        },
        KW_REPORT => {
            advance();
            let n = node(N_REPORT, t);
            loop {
                let r = node(N_REPORTITEM, tok());
                if at(TK_IDENT) && kind_at(1) == TK_ASSIGN {
                    nodes[r].aux = advance();
                    advance();
                }
                add(r, expr());
                add(n, r);
                if !accept(TK_COMMA) { break; }
            }
            return n;
        },
        _ => {
            return fail_expected("a statement of a step",
                                 "start the line with bc, release, birth, death, solve, adapt, output or report");
        },
    }
}

fn step_decl() -> i64 {
    advance();
    let n = node(N_STEP, expect_ident("a step name"));
    expect(TK_COLON, "':'");
    add(n, procedure());
    expect(TK_LBRACE, "'{'");
    accept(TK_NL);
    while !at(TK_RBRACE) && !at(TK_EOF) {
        add(n, step_statement());
        end_statement();
    }
    expect(TK_RBRACE, "'}'");
    end_statement();
    return n;
}

fn expect_decl() -> i64 {
    advance();
    let n = node(N_EXPECT, expect_ident("a name"));
    expect(TK_COLON, "':'");
    add(n, add_expr());
    expect(TK_EQ, "'=='");
    add(n, add_expr());
    expect(KW_WITHIN, "'within'");
    add(n, expr());
    end_statement();
    return n;
}

fn allow_decl() -> i64 {
    advance();
    let n = node(N_ALLOW, expect_ident("the id of a check"));
    expect(TK_COLON, "':'");
    nodes[n].aux = expect(TK_STRING, "the reason, as a string");
    end_statement();
    return n;
}

fn scales_decl() -> i64 {
    let n = node(N_SCALES, advance());
    loop {
        add(n, expr());
        if !accept(TK_COMMA) { break; }
    }
    end_statement();
    return n;
}

fn manufacture_decl() -> i64 {
    let n = node(N_MANUFACTURE, advance());
    add(n, expr());
    expect(TK_COLON, "':'");
    add(n, procedure());
    add(n, parse_block());
    end_statement();
    return n;
}

fn tolerance_decl() -> i64 {
    advance();
    let n = node(N_TOLERANCE, expect_ident("a name"));
    expect(TK_COLON, "':'");
    add(n, expr());
    expect(KW_IN, "'in'");
    add(n, expr());
    expect(TK_DOTDOT, "'..'");
    add(n, expr());
    if at(TK_LBRACE) { add(n, parse_block()); }
    end_statement();
    return n;
}

// One item. It gives -1 after a syntax error.
fn parse_item(first: bool) -> i64 {
    if first != at(KW_MODEL) {
        if first {
            error_at(tok(), "header_order", "the header line is missing or not first",
                     "start the file with a header, for example 'model name: 3d, SI'");
        } else {
            error_at(tok(), "header_order", "the file has a second header line", "remove this header line");
        }
        if first { return 0; }
        return -1;
    }
    match kind() {
        KW_MODEL => header(),
        KW_USE => {
            let n = node(N_USE, advance());
            nodes[n].tok = expect(TK_STRING, "the path of the file, as a string");
            end_statement();
            n
        },
        KW_PARAM => param_decl(),
        KW_LET => let_decl(),
        KW_FN => fn_decl(),
        KW_TABLE => table_decl(),
        KW_MATERIAL => material_decl(),
        KW_MESH => discretization(N_MESH),
        KW_POINTS => discretization(N_POINTS),
        KW_PHYSICS => physics_decl(),
        KW_TEMPLATE => template_decl(),
        KW_METHOD => method_decl(),
        KW_BC => {
            let n = bc_statement();
            end_statement();
            n
        },
        KW_CONTACT => pair_decl(N_CONTACT),
        KW_INTERFACE => pair_decl(N_INTERFACE),
        KW_SYMMETRY => procedure_decl(N_SYMMETRY),
        KW_STEP => step_decl(),
        KW_EXPECT => expect_decl(),
        KW_ALLOW => allow_decl(),
        KW_SCALES => scales_decl(),
        KW_STUDY => procedure_decl(N_STUDY),
        KW_MANUFACTURE => manufacture_decl(),
        KW_TOLERANCE => tolerance_decl(),
        _ => fail_expected("an item", "start the line with a keyword such as let, material or step"),
    }
}

// Moves to the next line that starts at column 1 with the keyword of an
// item, or to the end of the file.
fn recover() {
    loop {
        if at(TK_EOF) { return; }
        let t = advance();
        if tokens[t].kind == TK_NL && starts_item(kind()) && tokens[tok()].col == 1 { return; }
    }
}

// Parses file `f`, which `lex` has read. It gives the N_PROGRAM node.
fn parse(f: i64) -> i64 {
    ptoks.clear();
    for t, i in tokens {
        if t.file == f && t.kind != TK_COMMENT && t.kind != TK_DOC_COMMENT { ptoks.push(i); }
    }
    pos = 0;
    let prog = node(N_PROGRAM, tok());
    skip_newlines();
    if at(TK_EOF) {
        error_at(tok(), "header_order", "the file is empty", "start the file with a header, for example 'model name: 3d, SI'");
    }
    var first = true;
    while !at(TK_EOF) {
        let start = pos;
        let n = parse_item(first);
        if n > 0 {
            add(prog, n);
        } else if n < 0 || pos == start {
            recover();
        }
        first = false;
        skip_newlines();
    }
    return prog;
}
