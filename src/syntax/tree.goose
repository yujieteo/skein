import std;
namespace syntax;

struct Token {
    kind: i64,
    start: i64,
    len: i64,
    line: i64,
    col: i64,
    file: i64,
}

let TK_EOF = 0;
let TK_NL = 1;
let TK_INDENT = 2;
let TK_COMMENT = 3;
let TK_DOC_COMMENT = 4;
let TK_IDENT = 5;
let TK_NUMBER = 6;
let TK_PERCENT = 7;
let TK_STRING = 8;
let TK_SPACE = 9;
let TK_ERROR = 10;

let TK_LPAREN = 20;
let TK_RPAREN = 21;
let TK_LBRACKET = 22;
let TK_RBRACKET = 23;
let TK_LBRACE = 24;
let TK_RBRACE = 25;
let TK_LANGLE = 26;
let TK_RANGLE = 27;

let TK_COMMA = 30;
let TK_COLON = 31;
let TK_DOT = 32;
let TK_DOTDOT = 33;
let TK_ARROW = 34;
let TK_ASSIGN = 35;
let TK_EQ = 36;
let TK_NE = 37;
let TK_LT = 38;
let TK_LE = 39;
let TK_GT = 40;
let TK_GE = 41;
let TK_PLUS = 42;
let TK_MINUS = 43;
let TK_STAR = 44;
let TK_SLASH = 45;
let TK_CARET = 46;
let TK_BANG = 47;
let TK_OROR = 48;
let TK_ANDAND = 49;
let TK_BAR = 50;
let TK_AMP = 51;

let KW_ALLOW = 100;
let KW_ASSUME = 101;
let KW_AT = 102;
let KW_BC = 103;
let KW_BIRTH = 104;
let KW_CONTACT = 105;
let KW_CRACK_TIP = 106;
let KW_DEATH = 107;
let KW_EQ = 108;
let KW_EVERY = 109;
let KW_EXPECT = 110;
let KW_FALSE = 111;
let KW_FN = 112;
let KW_IN = 113;
let KW_INTERFACE = 114;
let KW_LAYERS = 115;
let KW_LET = 116;
let KW_MATERIAL = 117;
let KW_MESH = 118;
let KW_METHOD = 119;
let KW_MODEL = 120;
let KW_ON = 121;
let KW_OUTPUT = 122;
let KW_PARAM = 123;
let KW_PHYSICS = 124;
let KW_PLASTIC = 125;
let KW_POINTS = 126;
let KW_REFINE = 127;
let KW_RELEASE = 128;
let KW_REPORT = 129;
let KW_SCALES = 130;
let KW_SOLVE = 131;
let KW_SOLVES = 132;
let KW_STATE = 133;
let KW_STEP = 134;
let KW_STUDY = 135;
let KW_SYMMETRY = 136;
let KW_TABLE = 137;
let KW_TEMPLATE = 138;
let KW_TRUE = 139;
let KW_USE = 140;
let KW_WHERE = 141;
let KW_WITHIN = 142;

let KEYWORDS = [
    "allow", "assume", "at", "bc", "birth", "contact", "crack_tip", "death", "eq", "every",
    "expect", "false", "fn", "in", "interface", "layers", "let", "material", "mesh", "method",
    "model", "on", "output", "param", "physics", "plastic", "points", "refine", "release", "report",
    "scales", "solve", "solves", "state", "step", "study", "symmetry", "table", "template", "true",
    "use", "where", "within",
];

fn keyword_kind(word: u8[:]) -> i64 {
    for k, i in KEYWORDS {
        if k == word { return KW_ALLOW + i; }
    }
    return TK_IDENT;
}

fn starts_item(kind: i64) -> bool {
    match kind {
        KW_MODEL, KW_USE, KW_PARAM, KW_LET, KW_FN, KW_TABLE, KW_MATERIAL, KW_MESH => true,
        KW_POINTS, KW_PHYSICS, KW_TEMPLATE, KW_METHOD, KW_BC, KW_CONTACT, KW_INTERFACE => true,
        KW_SYMMETRY, KW_STEP, KW_EXPECT, KW_ALLOW, KW_SCALES, KW_STUDY => true,
        _ => false,
    }
}

struct Node {
    kind: i64,
    tok: i64,
    aux: i64,
    first: i64,
    last: i64,
    next: i64,
}

let N_PROGRAM = 1;
let N_HEADER = 2;
let N_SPACE = 3;
let N_UNITSYS = 4;
let N_UNITFIELD = 5;
let N_USE = 6;
let N_PARAM = 7;
let N_LET = 9;
let N_FN = 10;
let N_FNPARAM = 11;
let N_TABLE = 12;
let N_COLUMN = 13;
let N_ROW = 14;
let N_MATERIAL = 15;
let N_PROP = 16;
let N_MESH = 17;
let N_POINTS = 18;
let N_PHYSICS = 19;
let N_ON = 20;
let N_USES = 21;
let N_UNKNOWN = 22;
let N_TEMPLATE = 23;
let N_TPARAM = 24;
let N_METHOD = 25;
let N_BC = 26;
let N_CONTACT = 27;
let N_INTERFACE = 28;
let N_SYMMETRY = 29;
let N_STEP = 30;
let N_EXPECT = 31;
let N_ALLOW = 32;
let N_SCALES = 33;
let N_STUDY = 34;

let N_MATASSIGN = 40;
let N_MATON = 41;
let N_ASSUME = 42;
let N_ASSUMPTION = 43;
let N_EQN = 44;
let N_STATE = 45;
let N_PLASTIC = 46;

let N_RELEASE = 50;
let N_BIRTH = 51;
let N_DEATH = 52;
let N_WHERE = 53;
let N_SOLVE = 54;
let N_OUTPUT = 56;
let N_EVERY = 57;
let N_REPORT = 58;
let N_REPORTITEM = 59;

let N_BLOCK = 60;
let N_SETTING = 61;
let N_REFINE = 62;
let N_LAYERS = 63;
let N_PROC = 64;
let N_LIST = 65;
let N_RANGE = 66;
let N_CRACK_TIP = 67;

let N_TYPE_NAME = 70;
let N_TYPE_APPLY = 71;
let N_DIM_MUL = 72;
let N_DIM_DIV = 73;
let N_DIM_POW = 74;

let N_UNIT = 80;
let N_UNIT_NAME = 81;
let N_UNIT_MUL = 82;
let N_UNIT_DIV = 83;
let N_UNIT_POW = 84;

let N_NUMBER = 90;
let N_PERCENT = 91;
let N_STRING = 92;
let N_BOOL = 93;
let N_NAME = 94;
let N_VECTOR = 95;
let N_RECORD = 96;
let N_BINARY = 97;
let N_UNARY = 98;
let N_CALL = 99;
let N_ARG = 100;
let N_MEMBER = 101;
let N_AT = 102;

let NODE_NAMES = [
    "", "program", "header", "space", "unitsys", "unitfield", "use", "param", "",
    "let", "fn", "fnparam", "table", "column", "row", "material", "prop", "mesh",
    "points", "physics", "on", "uses", "unknown", "template", "tparam", "method", "bc",
    "contact", "interface", "symmetry", "step", "expect", "allow", "scales", "study",
    "", "", "", "", "", "matassign", "maton", "assume", "assumption",
    "eq", "state", "plastic", "", "", "", "release", "birth", "death", "where", "solve",
    "", "output", "every", "report", "reportitem", "block", "setting", "refine",
    "layers", "proc", "list", "range", "crack_tip", "", "", "type", "type_apply", "dim_mul",
    "dim_div", "dim_pow", "", "", "", "", "", "unit", "unit_name", "unit_mul", "unit_div",
    "unit_pow", "", "", "", "", "", "number", "percent", "string", "bool", "name",
    "vector", "record", "binary", "unary", "call", "arg", "member", "at",
];

fn kind_name(kind: i64) -> const u8[:] {
    if kind < 0 || kind >= NODE_NAMES.len { return "?"; }
    return NODE_NAMES[kind];
}

fn format_halves(out: u8[>..]&, h: i64) {
    if h % 2 == 0 { format(out, h / 2); return; }
    if h < 0 { out.push('-'); }
    format(out, abs(h) / 2, ".5");
}

struct File {
    start: i64,
    len: i64,
    path: u8[..1024],
}

var source: u8[>..];
var files: File[>..];
var tokens: Token[>..];
var nodes: Node[>..];

struct Error {
    id: u8[..32],
    tok: i64,
    msg: u8[..512],
    fix: u8[..512],
}

var errors: Error[>..];

fn reset() {
    source.clear();
    files.clear();
    tokens.clear();
    nodes.clear();
    errors.clear();
}

fn add_file(path: u8[:], text: u8[:]) -> i64 {
    var f = File { start: source.len, len: text.len, path: [] };
    f.path.append(path[..min(path.len, 1024)]);
    files.push(f);
    source.append(text);
    return files.len - 1;
}

fn source_hash() -> u8[>..] {
    return str(hash(source));
}

fn error_at(tok: i64, id: u8[:], msg: u8[:], fix: u8[:]) {
    guard fix.len > 0 else { abort("error_at: an error has a fix"); }
    var e = Error { id: [], tok: tok, msg: [], fix: [] };
    e.id.append(id[..min(id.len, 32)]);
    e.msg.append(msg[..min(msg.len, 512)]);
    e.fix.append(fix[..min(fix.len, 512)]);
    errors.push(e);
}

fn text(t: i64) -> const u8[:] {
    let k = tokens[t];
    return source[k.start..k.start + k.len];
}

fn node(kind: i64, tok: i64) -> i64 {
    nodes.push(Node { kind: kind, tok: tok, aux: -1, first: -1, last: -1, next: -1 });
    return nodes.len - 1;
}

fn add(parent: i64, child: i64) {
    let p = nodes[parent].last;
    if p < 0 { nodes[parent].first = child; } else { nodes[p].next = child; }
    nodes[parent].last = child;
}

fn child_count(n: i64) -> i64 {
    var count = 0;
    var c = nodes[n].first;
    while c >= 0 { count++; c = nodes[c].next; }
    return count;
}

fn child(n: i64, i: i64) -> i64 {
    var c = nodes[n].first;
    var k = 0;
    while c >= 0 && k < i { c = nodes[c].next; k++; }
    return c;
}

fn child_of_kind(n: i64, kind: i64) -> i64 {
    var c = nodes[n].first;
    while c >= 0 && nodes[c].kind != kind { c = nodes[c].next; }
    return c;
}

recursive fn dump(n: i64, out: u8[>..]&) {
    let k = nodes[n];
    out.push('(');
    out.append(kind_name(k.kind));
    if k.tok >= 0 && k.kind != N_PROGRAM {
        out.append(" \"");
        out.append(text(k.tok));
        out.push('"');
    }
    if k.aux >= 0 && k.kind != N_DIM_POW && k.kind != N_UNIT_POW {
        out.append(" aux=\"");
        out.append(text(k.aux));
        out.push('"');
    }
    if k.kind == N_DIM_POW || k.kind == N_UNIT_POW { out.push(' '); format_halves(out, k.aux); }
    var c = k.first;
    while c >= 0 {
        out.push(' ');
        dump(c, out);
        c = nodes[c].next;
    }
    out.push(')');
}
