// The types of Skein values (spec §5) and the rules that join them (spec
// §5.7, §6.1, §7.2).
import std;
namespace check;
import .dims;

// The kinds of a type.
let K_ANY = 0;       // the type of a value after an error. It matches every type.
let K_NUM = 1;       // a quantity: `Int`, `Real`, `Angle` or a dimension, of a rank
let K_BOOL = 2;
let K_STR = 3;
let K_GEOM = 4;      // a geometry value: `geo` is its topological dimension
let K_MESH = 5;
let K_POINTS = 6;
let K_RESULTS = 7;   // `aux` is the index into `results`
let K_WORD = 8;      // a word in a setting, for example `tet` or `newton`
let K_RECORD = 9;    // `{ a: 1, b: 2 }`
let K_LIST = 10;     // `[a, b]`
let K_DECL = 11;     // a name that is not a value: a material, a physics, a step

// The scalar kind of a dimensionless quantity (spec §5.1).
let B_REAL = 0;
let B_INT = 1;
let B_ANGLE = 2;

// The number systems (spec §5.1).
let NS_REAL = 0;
let NS_COMPLEX = 1;
let NS_QUATERNION = 2;
let NS_DUAL = 3;
let NS_NAMES = ["", "Complex", "Quaternion", "Dual"];

struct Ty {
    kind: i64,
    dim: i64[7],
    rank: i64,        // 0, 1, 2 or 4 (spec §5.2)
    num: i64,         // an NS_ constant
    base: i64,        // a B_ constant
    geo: i64,         // the topological dimension of a K_GEOM
    aux: i64,         // the results index of a K_RESULTS, else -1
    zero: bool,       // the literal 0 with no unit (spec §5.7, exception 1)
    list: bool,       // a table column, or an expression of one (spec §5.6)
    varying: bool,    // depends on x, y, z, t or a field (spec §7.6)
    percent: bool,    // a percent literal
}

fn any_ty() -> Ty {
    return Ty { kind: K_ANY, dim: NO_DIM, rank: 0, num: NS_REAL, base: B_REAL, geo: 0, aux: -1,
                zero: false, list: false, varying: false, percent: false };
}

fn kind_ty(kind: i64) -> Ty {
    var t = any_ty();
    t.kind = kind;
    return t;
}

fn num_ty(dim: i64[7], rank: i64) -> Ty {
    var t = kind_ty(K_NUM);
    t.dim = dim;
    t.rank = rank;
    return t;
}

fn base_ty(base: i64) -> Ty {
    var t = kind_ty(K_NUM);
    t.base = base;
    return t;
}

fn geom_ty(geo: i64) -> Ty {
    var t = kind_ty(K_GEOM);
    t.geo = geo;
    return t;
}

fn is_num(t: Ty) -> bool { t.kind == K_NUM }
fn is_any(t: Ty) -> bool { t.kind == K_ANY }
fn is_scalar(t: Ty) -> bool { t.kind == K_NUM && t.rank == 0 }
fn is_int(t: Ty) -> bool { t.kind == K_NUM && t.rank == 0 && t.base == B_INT && is_dimless(t.dim) && t.num == NS_REAL }
fn is_dimless_scalar(t: Ty) -> bool { t.kind == K_NUM && t.rank == 0 && is_dimless(t.dim) }

let GEO_NAMES = ["Point", "Curve", "Surface", "Solid"];

// Writes the type `t` as a Skein type, for example `Vector<Pressure>`.
fn format_ty(out: u8[>..]&, t: Ty) {
    match t.kind {
        K_ANY => out.append("an unknown type"),
        K_BOOL => out.append("Bool"),
        K_STR => out.append("Str"),
        K_GEOM => out.append(GEO_NAMES[clamp(t.geo, 0, 3)]),
        K_MESH => out.append("Mesh"),
        K_POINTS => out.append("Points"),
        K_RESULTS => out.append("Results"),
        K_WORD => out.append("a word"),
        K_RECORD => out.append("a record"),
        K_LIST => out.append("a list"),
        K_DECL => out.append("a declaration that is not a value"),
        _ => {
            if t.rank == 1 { out.append("Vector<"); }
            if t.rank == 2 { out.append("Tensor2<"); }
            if t.rank == 4 { out.append("Tensor4<"); }
            if t.num != NS_REAL { format(out, NS_NAMES[t.num], "<"); }
            if !is_dimless(t.dim) {
                format_dim(out, t.dim);
            } else if t.base == B_INT {
                out.append("Int");
            } else if t.base == B_ANGLE {
                out.append("Angle");
            } else {
                out.append("Real");
            }
            if t.num != NS_REAL { out.push('>'); }
            if t.rank != 0 { out.push('>'); }
        },
    }
}

fn ty_str(t: Ty) -> u8[>..] {
    var s: u8[>..] = [];
    format_ty(s, t);
    return s;
}

fn dim_str(d: i64[7]) -> u8[>..] {
    var s: u8[>..] = [];
    format_dim(s, d);
    return s;
}

// The result of joining `a` and `b` in a value: the varying and list flags
// of both.
fn joined(t: Ty, a: Ty, b: Ty) -> Ty {
    var r = t;
    r.varying = a.varying || b.varying;
    r.list = a.list || b.list;
    return r;
}

fn rank_name(rank: i64) -> const u8[:] {
    match rank {
        0 => "a scalar",
        1 => "a Vector",
        2 => "a Tensor2",
        _ => "a Tensor4",
    }
}
