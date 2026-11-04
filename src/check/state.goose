import std;
import dictionary;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;

struct Diag {
    id: u8[..32],
    file: u8[..256],
    fnum: i64,
    line: i64,
    col: i64,
    msg: u8[..512],
    fix: u8[..512],
    lint: bool = false,
}

var diags: Diag[>..];

var report_tok = -1;

fn error(tok: i64, id: u8[:], msg: u8[:], fix: u8[:]) { add_diag(tok, id, msg, fix, false); }

fn add_diag(tok: i64, id: u8[:], msg: u8[:], fix: u8[:], lint: bool) {
    guard id.len <= 32 else { abort(str("add_diag: the id '", id, "' is longer than the 32 bytes of Diag.id")); }
    let at = if report_tok >= 0 { report_tok } else { tok };
    let k = syntax::tokens[max(at, 0)];
    for d in diags {
        if d.fnum == k.file && d.line == k.line && d.col == k.col && d.id == id { return; }
    }
    var d = Diag { id: [], file: [], fnum: k.file, line: k.line, col: k.col, msg: [], fix: [], lint: lint };
    let path = syntax::files[k.file].path;
    d.id.append(id[..min(id.len, 32)]);
    d.file.append(path[..min(path.len, 256)]);
    d.msg.append(msg[..min(msg.len, 512)]);
    let f = if fix.len > 0 { fix } else { "see docs/checks.md" };
    d.fix.append(f[..min(f.len, 512)]);
    diags.push(d);
}

type Name = u8[..64];

fn key(s: u8[:]) -> Name {
    var k: Name = [];
    k.append(s[..min(s.len, 64)]);
    return k;
}

let S_BUILTIN = 0;
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

struct Sym {
    kind: i64,
    tok: i64,
    node: i64,
    ty: Ty,
    has_val: bool,
    val: f64,
    phys: i64,
}

var syms: Sym[>..];
var visible = dictionary<Name, i64> {};
var declared = dictionary<Name, i64> {};

struct Local {
    name: Name,
    ty: Ty,
}

var locals: Local[>..<];

struct Prop {
    mat: i64,
    tok: i64,
    ty: Ty,
}

var props: Prop[>..];

struct Phys {
    sym: i64,
    node: i64,
    region: i64,
    template: i64,
    first: i64,
    count: i64,
    has_method: bool,
    plastic: bool,
    incompressible: bool,
    convective: bool,
}

var physics: Phys[>..];
struct PhysMat { phys: i64, mat: i64 }
var phys_mats: PhysMat[>..];
struct Flux { phys: i64, node: i64, unknown: i64, curl: bool, curve: bool }
var fluxes: Flux[>..];
struct PropSeen { phys: i64, name: Name }
var prop_seen: PropSeen[>..];

struct ResultField { res: i64, name: Name, ty: Ty }
var result_fields: ResultField[>..];
var results_count = 0;

struct Bc {
    step: i64,
    region: i64,
    unknown: i64,
    comp: Name,
    natural: bool,
    node: i64,
}

var bcs: Bc[>..];
var rigid_points: u8[..64][>..];

var space = 3;
var axi = false;
var nondim = false;
var has_scales = false;
var thickness = false;
var model_prog = -1;
var header_tok = 0;
var cur_phys = -1;
var cur_step = -1;
var cur_proc: u8[..32];
var last_proc: u8[..32];
var cur_count = -1;
var eig_proc: u8[..32];
var eig_count = -1;
var step_count = 0;
var setting_mode = false;
var in_expect = false;
var in_report = false;
var in_bc = false;
var in_natural = false;
var last_reports: u8[..64][>..];
var scale_dims: i64[7][>..];
var used_dims: i64[7][>..];

struct Free {
    physics: u8[..64],
    field: u8[..64],
    region: u8[..512],
    condition: u8[..256],
}

var free: Free[>..];

struct Chosen {
    physics: u8[..64],
    mesh: u8[..64],
    order: u8[..256],
    stabilize: u8[..16],
}

var chosen: Chosen[>..];

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
    rigid_points.clear();
    space = 3;
    axi = false;
    nondim = false;
    set_nondim(false);
    has_scales = false;
    thickness = false;
    model_prog = -1;
    cur_phys = -1;
    cur_step = -1;
    cur_proc.clear();
    last_proc.clear();
    cur_count = -1;
    eig_proc.clear();
    eig_count = -1;
    step_count = 0;
    setting_mode = false;
    in_expect = false;
    in_report = false;
    in_bc = false;
    in_natural = false;
    report_tok = -1;
    last_reports.clear();
    scale_dims.clear();
    used_dims.clear();
    free.clear();
    chosen.clear();
}

fn add_sym(kind: i64, tok: i64, node: i64, ty: Ty) -> i64 {
    syms.push(Sym { kind: kind, tok: tok, node: node, ty: ty, has_val: false, val: 0.0, phys: cur_phys });
    return syms.len - 1;
}

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

fn norm_text(n: i64) -> u8[>..] {
    var out: u8[>..] = [];
    for c in node_text(n) {
        if c != ' ' && c != '\t' && c != '\r' && c != '\n' { out.push(c); }
    }
    return out;
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
