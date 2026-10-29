// The syntax tree of a Skein program (spec §2, §3).
//
// The lexer writes the tokens and the parser writes the nodes into the
// global arrays below. All links are indices into these arrays, and -1 is
// "none". A node has an ordered list of children: `first` is the first
// child, `next` is the next sibling and `last` is the last child. The list
// of node kinds below gives the children of each kind in their order.
//
// The source text of all files is in `source`, one file after the other, so
// a token offset is unique over all files. `files` gives the range of each
// file. A checker can parse a second file (for `results`, spec §15.1) into
// the same arrays.
import std;
namespace syntax;

// --- tokens ----------------------------------------------------------------

struct Token {
    kind: i64,     // a TK_ or KW_ constant
    start: i64,    // offset of the first byte in `source`
    len: i64,      // length in bytes
    line: i64,     // 1-based line in its file
    col: i64,      // 1-based column in Unicode code points
    file: i64,     // index into `files`
}

// End of the file.
let TK_EOF = 0;
// The end of a statement (spec §2.2).
let TK_NL = 1;
// The start of a row in the body of a table (spec §3, `tabledecl`).
let TK_INDENT = 2;
// A `//` comment, up to the end of the line. The parser skips it.
let TK_COMMENT = 3;
// A `///` documentation comment. The parser skips it.
let TK_DOC_COMMENT = 4;
let TK_IDENT = 5;
// A number with no sign (spec §2.4).
let TK_NUMBER = 6;
// A number and `%`, for example `5%`.
let TK_PERCENT = 7;
// A string with its quotes.
let TK_STRING = 8;
// `1d`, `2d` or `3d`. A token only in the header (spec §3).
let TK_SPACE = 9;
// A byte that starts no token.
let TK_ERROR = 10;

// Brackets. `<` and `>` around the argument of a tensor or number-system
// type (`Vector<Length>`) are TK_LANGLE and TK_RANGLE. Elsewhere they are
// the comparisons TK_LT and TK_GT.
let TK_LPAREN = 20;
let TK_RPAREN = 21;
let TK_LBRACKET = 22;
let TK_RBRACKET = 23;
let TK_LBRACE = 24;
let TK_RBRACE = 25;
let TK_LANGLE = 26;
let TK_RANGLE = 27;

// Punctuation and operators.
let TK_COMMA = 30;
let TK_COLON = 31;
let TK_DOT = 32;
let TK_DOTDOT = 33;     // ..
let TK_ARROW = 34;      // ->
let TK_ASSIGN = 35;     // =
let TK_EQ = 36;         // ==
let TK_NE = 37;         // !=
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
let TK_OROR = 48;       // ||
let TK_ANDAND = 49;     // &&
let TK_BAR = 50;        // | (union)
let TK_AMP = 51;        // & (intersection)

// Keywords (spec §2.5), in the order of that list.
let KW_ADAPT = 100;
let KW_ALLOW = 101;
let KW_ASSUME = 102;
let KW_AT = 103;
let KW_BC = 104;
let KW_BIRTH = 105;
let KW_CONTACT = 106;
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
let KW_MANUFACTURE = 117;
let KW_MATERIAL = 118;
let KW_MESH = 119;
let KW_METHOD = 120;
let KW_MODEL = 121;
let KW_ON = 122;
let KW_OUTPUT = 123;
let KW_PARAM = 124;
let KW_PHYSICS = 125;
let KW_PLASTIC = 126;
let KW_POINTS = 127;
let KW_REFINE = 128;
let KW_RELEASE = 129;
let KW_REPORT = 130;
let KW_SCALES = 131;
let KW_SOLVE = 132;
let KW_SOLVES = 133;
let KW_STATE = 134;
let KW_STEP = 135;
let KW_STUDY = 136;
let KW_SYMMETRY = 137;
let KW_TABLE = 138;
let KW_TEMPLATE = 139;
let KW_TOL = 140;
let KW_TOLERANCE = 141;
let KW_TRUE = 142;
let KW_USE = 143;
let KW_WHERE = 144;
let KW_WITHIN = 145;

let KEYWORDS = [
    "adapt", "allow", "assume", "at", "bc", "birth", "contact", "death", "eq", "every",
    "expect", "false", "fn", "in", "interface", "layers", "let", "manufacture",
    "material", "mesh", "method", "model", "on", "output", "param", "physics", "plastic",
    "points", "refine", "release", "report", "scales", "solve", "solves", "state", "step",
    "study", "symmetry", "table", "template", "tol", "tolerance", "true", "use", "where",
    "within",
];

// The keyword kind of `word`, or TK_IDENT.
fn keyword_kind(word: u8[:]) -> i64 {
    for k, i in KEYWORDS {
        if k == word { return KW_ADAPT + i; }
    }
    return TK_IDENT;
}

// --- nodes -------------------------------------------------------------------

struct Node {
    kind: i64,     // an N_ constant
    tok: i64,      // the main token: the table below says which
    aux: i64,      // a second token or a value: the table below says which, else -1
    first: i64,    // first child
    last: i64,     // last child
    next: i64,     // next sibling
}

// The kinds. "tok" and "aux" give the meaning of those fields. "[x]" is a
// child that can be absent, "x*" is zero or more children. A child is an
// expression where no kind is given.

// Items (spec §3, `item`).
let N_PROGRAM = 1;       // N_HEADER, then the items in source order
let N_HEADER = 2;        // tok: model name. N_SPACE, N_UNITSYS
let N_SPACE = 3;         // tok: TK_SPACE, or the ident `axi`
let N_UNITSYS = 4;       // tok: the system name, `nondim` or `{`. N_UNITFIELD* for `{`
let N_UNITFIELD = 5;     // tok: field name. N_UNIT
let N_USE = 6;           // tok: the string
let N_PARAM = 7;         // tok: name. type, value, [N_TOL], [N_RANGE]
let N_TOL = 8;           // tok: `tol`. the tolerance
let N_LET = 9;           // tok: name. type, value
let N_FN = 10;           // tok: name. N_FNPARAM*, result type, body
let N_FNPARAM = 11;      // tok: name. type
let N_TABLE = 12;        // tok: name. N_COLUMN*, N_ROW*
let N_COLUMN = 13;       // tok: name. type, [N_UNIT]
let N_ROW = 14;          // tok: first token of the row. the values
let N_MATERIAL = 15;     // tok: name. N_PROP*
let N_PROP = 16;         // tok: name. type, value
let N_MESH = 17;         // tok: name. region, N_BLOCK
let N_POINTS = 18;       // tok: name. region, N_BLOCK
let N_PHYSICS = 19;      // tok: name. [N_ON], [N_USES], N_UNKNOWN*, statements
let N_ON = 20;           // tok: `on`. the region
let N_USES = 21;         // tok: template name. N_ARG*
let N_UNKNOWN = 22;      // tok: name. type, [initial value]
let N_TEMPLATE = 23;     // tok: name. N_TPARAM*, N_UNKNOWN*, statements
let N_TPARAM = 24;       // tok: name. type, [default]
let N_METHOD = 25;       // tok: physics name, aux: method kind. N_BLOCK
let N_BC = 26;           // tok: `bc`. region, target, value
let N_CONTACT = 27;      // tok: name. region, region, N_BLOCK
let N_INTERFACE = 28;    // tok: name. region, region, N_BLOCK
let N_SYMMETRY = 29;     // tok: name. N_PROC, N_BLOCK
let N_STEP = 30;         // tok: name. N_PROC, step statements
let N_EXPECT = 31;       // tok: name. left, right, within
let N_ALLOW = 32;        // tok: check id, aux: the reason string
let N_SCALES = 33;       // tok: `scales`. the scales
let N_STUDY = 34;        // tok: name. N_PROC, N_BLOCK
let N_MANUFACTURE = 35;  // tok: `manufacture`. region, N_PROC, N_BLOCK
let N_TOLERANCE = 36;    // tok: name. value, low, high, [N_BLOCK]

// Physics statements (spec §3, `physstmt`). A physics also holds N_LET.
let N_MATASSIGN = 40;    // tok: `material`. N_MATON*
let N_MATON = 41;        // tok: material name. [region]
let N_ASSUME = 42;       // tok: `assume`. N_ASSUMPTION*
let N_ASSUMPTION = 43;   // tok: name. arguments
let N_EQN = 44;          // tok: `eq`. left, right
let N_STATE = 45;        // tok: name. type, initial value
let N_PLASTIC = 46;      // tok: `plastic`. N_BLOCK

// Step statements (spec §3, `stepstmt`). A step also holds N_BC.
let N_RELEASE = 50;      // tok: `release`. region
let N_BIRTH = 51;        // tok: `birth`. regions, [N_WHERE]
let N_DEATH = 52;        // tok: `death`. regions, [N_WHERE]
let N_WHERE = 53;        // tok: `where`. condition
let N_SOLVE = 54;        // tok: `solve`. N_BLOCK
let N_ADAPT = 55;        // tok: `adapt`. N_BLOCK
let N_OUTPUT = 56;       // tok: `output`. values, [N_EVERY]
let N_EVERY = 57;        // tok: `every`. count
let N_REPORT = 58;       // tok: `report`. N_REPORTITEM*
let N_REPORTITEM = 59;   // tok: first token, aux: name or -1. value

// Blocks and procedures (spec §3, `block`, `procedure`).
let N_BLOCK = 60;        // tok: `{`. N_SETTING*, N_REFINE*, N_LAYERS* in source order
let N_SETTING = 61;      // tok: name. value: an expression, N_LIST or N_RANGE
let N_REFINE = 62;       // tok: `refine`. region, value
let N_LAYERS = 63;       // tok: `layers`. region, value
let N_PROC = 64;         // tok: procedure name. N_ARG*
let N_LIST = 65;         // tok: `[`. items, [N_UNIT]
let N_RANGE = 66;        // tok: `..`. low, high

// Types (spec §3, `type`).
let N_TYPE_NAME = 70;    // tok: name, for example `Real`, `Solid` or `Length`
let N_TYPE_APPLY = 71;   // tok: `Vector`, `Tensor2`, `Tensor4`, `Complex`, `Quaternion` or `Dual`. type
let N_DIM_MUL = 72;      // tok: `*`. left, right
let N_DIM_DIV = 73;      // tok: `/`. left, right
let N_DIM_POW = 74;      // tok: `^`, aux: the exponent. base

// Units (spec §6.3). N_UNIT wraps the unit expression of a literal.
let N_UNIT = 80;         // tok: first token. N_UNIT_NAME, N_UNIT_MUL, N_UNIT_DIV or N_UNIT_POW
let N_UNIT_NAME = 81;    // tok: name, for example `mm` or `GPa`
let N_UNIT_MUL = 82;     // tok: `*`. left, right
let N_UNIT_DIV = 83;     // tok: `/`. left, right
let N_UNIT_POW = 84;     // tok: `^`, aux: the exponent. base

// Expressions (spec §3, `expr`). Parentheses make no node.
let N_NUMBER = 90;       // tok: TK_NUMBER. [N_UNIT]
let N_PERCENT = 91;      // tok: TK_PERCENT
let N_STRING = 92;       // tok: TK_STRING
let N_BOOL = 93;         // tok: `true` or `false`
let N_NAME = 94;         // tok: the identifier
let N_VECTOR = 95;       // tok: `(`. components, [N_UNIT]
let N_RECORD = 96;       // tok: `{`. N_ARG*, each with a name
let N_BINARY = 97;       // tok: the operator. left, right
let N_UNARY = 98;        // tok: `-` or `!`. operand
let N_CALL = 99;         // tok: `(`. callee, N_ARG*
let N_ARG = 100;         // tok: first token, aux: name or -1. value
let N_MEMBER = 101;      // tok: the member name. object
let N_AT = 102;          // tok: `at`. field, point

let NODE_NAMES = [
    "", "program", "header", "space", "unitsys", "unitfield", "use", "param", "tol",
    "let", "fn", "fnparam", "table", "column", "row", "material", "prop", "mesh",
    "points", "physics", "on", "uses", "unknown", "template", "tparam", "method", "bc",
    "contact", "interface", "symmetry", "step", "expect", "allow", "scales", "study",
    "manufacture", "tolerance", "", "", "", "matassign", "maton", "assume", "assumption",
    "eq", "state", "plastic", "", "", "", "release", "birth", "death", "where", "solve",
    "adapt", "output", "every", "report", "reportitem", "block", "setting", "refine",
    "layers", "proc", "list", "range", "", "", "", "type", "type_apply", "dim_mul",
    "dim_div", "dim_pow", "", "", "", "", "", "unit", "unit_name", "unit_mul", "unit_div",
    "unit_pow", "", "", "", "", "", "number", "percent", "string", "bool", "name",
    "vector", "record", "binary", "unary", "call", "arg", "member", "at",
];

// The name of a node kind, as `dump` writes it.
fn kind_name(kind: i64) -> const u8[:] {
    if kind < 0 || kind >= NODE_NAMES.len { return "?"; }
    return NODE_NAMES[kind];
}

// --- the tables ----------------------------------------------------------------

struct File {
    start: i64,    // offset of the first byte in `source`
    len: i64,
    path: u8[..1024],
}

var source: u8[>..];
var files: File[>..];
var tokens: Token[>..];
var nodes: Node[>..];

// Empties all the tables, for a new run.
fn reset() {
    source.clear();
    files.clear();
    tokens.clear();
    nodes.clear();
}

// Appends the text of a file to `source`, and gives its index in `files`.
fn add_file(path: u8[:], text: u8[:]) -> i64 {
    var f = File { start: source.len, len: text.len, path: [] };
    f.path.append(path[..min(path.len, 1024)]);
    files.push(f);
    source.append(text);
    return files.len - 1;
}

// The text of token `t`.
fn text(t: i64) -> const u8[:] {
    let k = tokens[t];
    return source[k.start..k.start + k.len];
}

// A new node with no children. It gives the index of the node.
fn node(kind: i64, tok: i64) -> i64 {
    nodes.push(Node { kind: kind, tok: tok, aux: -1, first: -1, last: -1, next: -1 });
    return nodes.len - 1;
}

// Appends `child` to the children of `parent`.
fn add(parent: i64, child: i64) {
    let p = nodes[parent].last;
    if p < 0 { nodes[parent].first = child; } else { nodes[p].next = child; }
    nodes[parent].last = child;
}

// The number of children of node `n`.
fn child_count(n: i64) -> i64 {
    var count = 0;
    var c = nodes[n].first;
    while c >= 0 { count++; c = nodes[c].next; }
    return count;
}

// Child `i` of node `n` (0-based), or -1.
fn child(n: i64, i: i64) -> i64 {
    var c = nodes[n].first;
    var k = 0;
    while c >= 0 && k < i { c = nodes[c].next; k++; }
    return c;
}

// The first child of `n` with kind `kind`, or -1.
fn child_of_kind(n: i64, kind: i64) -> i64 {
    var c = nodes[n].first;
    while c >= 0 && nodes[c].kind != kind { c = nodes[c].next; }
    return c;
}

// Writes node `n` and its children as one S-expression, for example
// `(let "L" (type "Length") (number "200" (unit "mm" (unit_name "mm"))))`.
// The tests of the parser compare these.
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
    if k.kind == N_DIM_POW || k.kind == N_UNIT_POW { format(out, " ", k.aux); }
    var c = k.first;
    while c >= 0 {
        out.push(' ');
        dump(c, out);
        c = nodes[c].next;
    }
    out.push(')');
}
