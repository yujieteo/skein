// The state of one run of the checker: the diagnostics, the symbols and
// the facts about the model that the rules need.
import std;
import dictionary;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;

// --- diagnostics (docs/checks.md §1, docs/command_line.md §6.1) -----------------

struct Diag {
    id: u8[..32],
    file: u8[..256],
    line: i64,
    col: i64,
    msg: u8[..512],
    fix: u8[..512],
}

var diags: Diag[>..];

// When it is a token, each error goes to it. The checker sets it while it
// checks the source of a template for a physics of the model.
var report_tok = -1;

// Records the error `id` at token `tok`. `fix` is never empty.
fn error(tok: i64, id: u8[:], msg: u8[:], fix: u8[:]) {
    let at = if report_tok >= 0 { report_tok } else { tok };
    let k = syntax::tokens[max(at, 0)];
    for d in diags {
        if d.line == k.line && d.col == k.col && d.id == id { return; }
    }
    var d = Diag { id: [], file: [], line: k.line, col: k.col, msg: [], fix: [] };
    let path = syntax::files[k.file].path;
    d.id.append(id[..min(id.len, 32)]);
    d.file.append(path[..min(path.len, 256)]);
    d.msg.append(msg[..min(msg.len, 512)]);
    let f = if fix.len > 0 { fix } else { "see docs/checks.md" };
    d.fix.append(f[..min(f.len, 512)]);
    diags.push(d);
}

// A name, kept inline: a global must not hold a slice of the source text,
// because syntax::reset clears it.
type Name = u8[..64];

fn key(s: u8[:]) -> Name {
    var k: Name = [];
    k.append(s[..min(s.len, 64)]);
    return k;
}

// --- symbols ----------------------------------------------------------------------

// The kinds of a symbol.
let S_BUILTIN = 0;    // a built-in constant or function (Appendix A)
let S_LET = 1;
let S_PARAM = 2;
let S_FN = 3;
let S_TABLE = 4;
let S_MATERIAL = 5;
let S_MESH = 6;
let S_POINTS = 7;
let S_PHYSICS = 8;
let S_UNKNOWN = 9;
let S_STATE = 10;
let S_SYMMETRY = 11;
let S_STEP = 12;
let S_EXPECT = 13;
let S_TEMPLATE = 14;
let S_PLET = 15;      // a `let` of a physics block or of a template

struct Sym {
    kind: i64,
    tok: i64,         // the token of the name, or -1 for a built-in
    node: i64,        // the declaration, or -1
    ty: Ty,
    has_val: bool,    // `val` is the constant value in SI
    val: f64,
    phys: i64,        // the index in `physics` of an unknown, state or physics let, else -1
}

var syms: Sym[>..];
// The names that are visible now, and their symbols.
var visible = dictionary<Name, i64> {};
// Each name that the model file declares, and the token of its declaration.
var declared = dictionary<Name, i64> {};

// A name of a function parameter, a template parameter or a template unknown.
struct Local {
    name: Name,
    ty: Ty,
}

var locals: Local[>..<];

// A material property.
struct Prop {
    mat: i64,      // the symbol of the material
    tok: i64,
    ty: Ty,
}

var props: Prop[>..];

// A physics block.
struct Phys {
    sym: i64,
    node: i64,
    region: i64,      // the region expression, or -1
    template: i64,    // the N_TEMPLATE node of the standard library, or -1
    first: i64,       // the symbol of the first unknown
    count: i64,       // the number of unknowns
    has_method: bool,
    plastic: bool,
    incompressible: bool,
}

var physics: Phys[>..];
// The materials of each physics: pairs of a physics index and a material symbol.
struct PhysMat { phys: i64, mat: i64 }
var phys_mats: PhysMat[>..];
// The fluxes: `div(F)` in equation `k` of a physics makes `F*n` a natural
// condition for unknown `k` (spec §12.1).
struct Flux { phys: i64, name: Name, unknown: i64 }
var fluxes: Flux[>..];
// The pairs (physics, property name) that already have a `material_prop_*` error.
struct PropSeen { phys: i64, name: Name }
var prop_seen: PropSeen[>..];

// The unknowns of another model, for `results` (spec §15.1).
struct ResultField { res: i64, name: Name, ty: Ty }
var result_fields: ResultField[>..];
var results_count = 0;

// The boundary conditions, for `bc_conflict` and the free boundary report.
struct Bc {
    step: i64,        // the step index, or -1 at the top level
    region: i64,      // the region expression
    unknown: i64,     // the symbol of the unknown
    comp: Name,       // the component, or ""
    natural: bool,
    node: i64,
}

var bcs: Bc[>..];

// The facts about the model.
var space = 3;            // the space dimension: 1, 2 or 3. `axi` gives 2.
var axi = false;
var nondim = false;
var has_scales = false;
var thickness = false;    // 2D with a thickness (plane stress or strain), or `axi`
var model_file = 0;
var header_tok = 0;
var cur_phys = -1;        // the physics that the checker is in, or -1
var cur_step = -1;        // the step that the checker is in, or -1
var cur_proc: u8[..32];   // the procedure of that step
var last_proc: u8[..32];  // the procedure of the last step of the file
var step_count = 0;
var setting_mode = false; // an undeclared name is a word (a setting value)
var in_expect = false;
var last_reports: u8[..64][>..];   // the named reports of the last step
var scale_dims: i64[7][>..];
var used_dims: i64[7][>..];

// The default natural conditions, for the check report (command_line.md §6.2).
struct Free {
    physics: u8[..64],
    field: u8[..64],
    region: u8[..512],
    condition: u8[..256],
}

var free: Free[>..];

// Empties the state for a new file.
fn reset_state() {
    diags.clear();
    syms.clear();
    visible.clear();
    declared.clear();
    locals.clear();
    props.clear();
    physics.clear();
    phys_mats.clear();
    fluxes.clear();
    prop_seen.clear();
    result_fields.clear();
    results_count = 0;
    bcs.clear();
    space = 3;
    axi = false;
    nondim = false;
    set_nondim(false);
    has_scales = false;
    thickness = false;
    cur_phys = -1;
    cur_step = -1;
    cur_proc.clear();
    last_proc.clear();
    step_count = 0;
    setting_mode = false;
    in_expect = false;
    report_tok = -1;
    last_reports.clear();
    scale_dims.clear();
    used_dims.clear();
    free.clear();
}

fn add_sym(kind: i64, tok: i64, node: i64, ty: Ty) -> i64 {
    syms.push(Sym { kind: kind, tok: tok, node: node, ty: ty, has_val: false, val: 0.0, phys: cur_phys });
    return syms.len - 1;
}

// The text of the source from the first to the last token of node `n`,
// with its brackets balanced.
fn node_text(n: i64) -> const u8[:] {
    var a, b = span(n);
    let src = syntax::source;
    loop {
        var open = 0;
        var close = 0;
        for i in a..b {
            if src[i] == '(' { open++; }
            if src[i] == ')' { close++; }
        }
        if open > close && b < src.len && src[b] != '\n' {
            b++;
        } else if close > open && a > 0 && src[a - 1] != '\n' {
            a--;
        } else {
            break;
        }
    }
    return src[a..b];
}

recursive fn span(n: i64) -> i64, i64 {
    let k = syntax::nodes[n];
    var a = syntax::tokens[k.tok].start;
    var b = a + syntax::tokens[k.tok].len;
    var c = k.first;
    while c >= 0 {
        let ca, cb = span(c);
        a = min(a, ca);
        b = max(b, cb);
        c = syntax::nodes[c].next;
    }
    return a, b;
}
