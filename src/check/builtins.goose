// The types of calls of the built-in functions (Appendix A, spec §7.3, §8,
// §14.5, §15.1) and of the phase 1 gaps G6, G7 and G12.
import std;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .expr;
import .results;

// The arguments of a call. Only the first `MAX_ARGS` keep their type: the
// geometry functions with more points check the others as they go.
let MAX_ARGS = 8;

struct Args {
    count: i64,
    ty: Ty[8],
    node: i64[8],    // the value
    name: i64[8],    // the token of the argument name, or -1
}

fn args_varying(a: Args) -> bool {
    var v = false;
    for i in min(a.count, MAX_ARGS) { v = v || a.ty[i].varying; }
    return v;
}

fn arity(tok: i64, name: u8[:], a: Args, lo: i64, hi: i64) -> bool {
    if a.count >= lo && a.count <= hi { return true; }
    let want = if lo == hi { str(lo) } else { str(lo, " to ", hi) };
    error(tok, "type_mismatch", str(name, " takes ", want, " arguments, not ", a.count),
          str("give ", name, " ", want, " arguments (Appendix A)"));
    return false;
}

// A quantity of rank `rank`, or an error at argument `i`.
fn want_rank(tok: i64, name: u8[:], a: Args, i: i64, rank: i64) -> bool {
    let t = a.ty[i];
    if is_any(t) || (is_num(t) && (t.rank == rank || t.zero)) { return true; }
    error(first_tok(a.node[i]), "type_mismatch",
          str("argument ", i + 1, " of ", name, " is ", ty_str(t), ", not ", rank_name(rank)),
          str("give ", name, " ", rank_name(rank)));
    return false;
}

fn want_geom(tok: i64, name: u8[:], a: Args, i: i64) -> bool {
    let t = a.ty[i];
    if is_any(t) || t.kind == K_GEOM { return true; }
    error(first_tok(a.node[i]), "type_mismatch",
          str("argument ", i + 1, " of ", name, " is ", ty_str(t), ", not a geometry value"),
          "give a Solid, Surface, Curve or Point");
    return false;
}

// A point: a vector of `Length` (or 0).
fn want_point(name: u8[:], t: Ty, node: i64, i: i64) {
    if is_any(t) || t.zero { return; }
    if is_num(t) && t.rank == 1 && t.dim == LENGTH { return; }
    if t.kind == K_GEOM && t.geo == 0 { return; }
    error(first_tok(node), "type_mismatch",
          str("argument ", i + 1, " of ", name, " is ", ty_str(t), ", not a point"),
          "give a vector of Length, for example (0, 0, 0) mm");
}

fn want_length(name: u8[:], t: Ty, node: i64, i: i64) {
    if is_any(t) || t.zero { return; }
    if is_num(t) && t.rank == 0 && t.dim == LENGTH { return; }
    error(first_tok(node), "unit_mismatch",
          str("argument ", i + 1, " of ", name, " is ", ty_str(t), ", not a Length"),
          "give a Length, for example 5 mm");
}

fn want_dimless(tok: i64, name: u8[:], a: Args, i: i64) {
    let t = a.ty[i];
    if is_any(t) || (is_num(t) && is_dimless(t.dim)) { return; }
    error(first_tok(a.node[i]), "unit_transcendental",
          str(name, " takes a Real or an Angle, not ", ty_str(t)),
          "divide the argument by a reference value of its dimension");
}

// Checks one point argument after the first MAX_ARGS of a geometry function.
fn extra_point(name: u8[:], t: Ty, node: i64, i: i64) {
    if name == "polygon" || name == "spline" { want_point(name, t, node, i); }
}

fn space_2d(tok: i64, name: u8[:]) {
    if space != 2 {
        error(tok, "geometry_space", str(name, " is a 2D function, and the model is ", space, "D"),
              "use box, cylinder or sphere in a 3d model");
    }
}

fn space_3d(tok: i64, name: u8[:]) {
    if space != 3 {
        error(tok, "geometry_space", str(name, " is a 3D function, and the model is ", space, "D"),
              "use rect, circle or polygon in a 2d model");
    }
}

// The extra exponent of Length in an integral over a region of a 2D model
// with a thickness (spec §11.1, §14.5).
fn depth() -> i64 { if space == 2 && thickness { 1 } else { 0 } }

fn same(t: Ty, varying: bool) -> Ty {
    var r = t;
    r.zero = false;
    r.percent = false;
    r.varying = varying;
    return r;
}

fn real_scalar(a: Args) -> Ty {
    var r = base_ty(B_REAL);
    r.varying = args_varying(a);
    for i in min(a.count, MAX_ARGS) { r.list = r.list || a.ty[i].list; r.table = rows(r.table, a.ty[i].table); }
    return r;
}

// The type of the call of the built-in `name` at token `tok`.
fn builtin_ty(tok: i64, name: u8[:], a: Args) -> Ty {
    let v = args_varying(a);
    let t0 = if a.count > 0 { a.ty[0] } else { any_ty() };
    let t1 = if a.count > 1 { a.ty[1] } else { any_ty() };
    var listed = false;
    var table = -1;
    for i in min(a.count, MAX_ARGS) { listed = listed || a.ty[i].list; table = rows(table, a.ty[i].table); }
    var r = builtin_ty_inner(tok, name, a, t0, t1, v);
    // laminate sums over the plies of its layup column (spec §4.3), and sum
    // over the rows of a table (§5.6).
    if is_num(r) && name != "laminate" && name != "sum" {
        r.list = r.list || listed;
        r.table = rows(r.table, table);
    }
    return r;
}

fn builtin_ty_inner(tok: i64, name: u8[:], a: Args, t0: Ty, t1: Ty, v: bool) -> Ty {
    // Tensor algebra.
    if name == "sym" || name == "skew" || name == "dev" || name == "transpose" {
        guard arity(tok, name, a, 1, 1) && want_rank(tok, name, a, 0, 2) else { return any_ty(); }
        return same(t0, v);
    }
    if name == "inv" {
        guard arity(tok, name, a, 1, 1) && want_rank(tok, name, a, 0, 2) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, v);
        r.dim = NO_DIM - t0.dim;
        return r;
    }
    if name == "tr" || name == "von_mises" || name == "tresca" || name == "det" {
        guard arity(tok, name, a, 1, 1) && want_rank(tok, name, a, 0, 2) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, v);
        r.rank = 0;
        if name == "det" { r.dim = dim_pow(t0.dim, if axi { 3 } else { space }); }
        return r;
    }
    if name == "principal" {
        guard arity(tok, name, a, 1, 1) && want_rank(tok, name, a, 0, 2) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, v);
        r.rank = 1;
        return r;
    }
    if name == "dot" || name == "ddot" || name == "outer" || name == "cross" {
        guard arity(tok, name, a, 2, 2) && !is_any(t0) && !is_any(t1) else { return any_ty(); }
        if !is_num(t0) || !is_num(t1) || t0.rank == 0 || t1.rank == 0 {
            error(tok, "rank_mismatch", str(name, " takes two tensors, not ", ty_str(t0), " and ", ty_str(t1)),
                  "use '*' for a product with a scalar");
            return any_ty();
        }
        var rank = 0;
        if name == "dot" { rank = t0.rank + t1.rank - 2; }
        if name == "ddot" { rank = t0.rank + t1.rank - 4; }
        if name == "outer" { rank = t0.rank + t1.rank; }
        if name == "cross" {
            if t0.rank != 1 || t1.rank != 1 {
                error(tok, "rank_mismatch", "cross takes two vectors", "give cross two vectors");
                return any_ty();
            }
            rank = if space == 3 { 1 } else { 0 };
        }
        if rank < 0 || rank == 3 || rank > 4 {
            error(tok, "rank_mismatch", str(name, " of ", ty_str(t0), " and ", ty_str(t1), " is not defined"),
                  "check the ranks of the two arguments");
            return any_ty();
        }
        var r = num_ty(t0.dim + t1.dim, rank);
        r.num = max(t0.num, t1.num);
        r.varying = v;
        return r;
    }
    if name == "norm" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, v);
        r.rank = 0;
        return r;
    }
    if name == "strain" || name == "defgrad" {
        guard arity(tok, name, a, 1, 1) && want_rank(tok, name, a, 0, 1) && !is_any(t0) else { return any_ty(); }
        var r = num_ty(if name == "strain" { t0.dim - LENGTH } else { NO_DIM }, 2);
        r.num = t0.num;
        r.varying = v;
        return r;
    }
    if name == "hooke" {
        guard arity(tok, name, a, 3, 3) && !is_any(t0) && !is_any(a.ty[2]) else { return any_ty(); }
        want_rank(tok, name, a, 2, 2);
        want_dimless(tok, name, a, 1);
        var r = num_ty(t0.dim + a.ty[2].dim, 2);
        r.num = max(t0.num, a.ty[2].num);
        r.varying = v;
        return r;
    }
    if name == "elasticity" {
        guard arity(tok, name, a, 2, 2) && !is_any(t0) else { return any_ty(); }
        var r = num_ty(t0.dim, 4);
        r.varying = v;
        return r;
    }
    if name == "interp" {
        guard arity(tok, name, a, 3, 3) && !is_any(t1) else { return any_ty(); }
        var r = same(t1, v);
        r.list = false;
        r.table = -1;
        r.varying = a.ty[2].varying;
        let x = a.ty[0];
        if is_num(x) && is_num(a.ty[2]) && !a.ty[2].zero && x.dim != a.ty[2].dim {
            error(first_tok(a.node[2]), "unit_mismatch",
                  str("interp reads x of dimension ", dim_str(a.ty[2].dim), " in a column of ", dim_str(x.dim)),
                  "give x the dimension of the first column");
        }
        return r;
    }
    if name == "sum" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        guard t0.list else {
            error(tok, "sum_no_column", "sum takes an expression of a table column", "use a column, for example sum(t.col)");
            return any_ty();
        }
        guard t0.table != MIXED else {
            error(tok, "sum_mixed_tables", "sum takes the columns of one table", "sum over each table on its own");
            return any_ty();
        }
        var r = same(t0, v);
        r.list = false;
        r.table = -1;
        return r;
    }
    if name == "fourier" {
        guard arity(tok, name, a, 2, 2) && !is_any(t0) && !is_any(t1) else { return any_ty(); }
        var r = num_ty(t0.dim + t1.dim - LENGTH, t1.rank + 1);
        r.varying = v;
        return r;
    }
    // Scalar functions (spec §6.1).
    if name == "abs" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, v);
        r.num = NS_REAL;
        return r;
    }
    if name == "sqrt" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        for e in t0.dim {
            if e % 2 != 0 {
                error(tok, "unit_fractional", str("sqrt of ", dim_str(t0.dim), " has a fractional exponent"),
                      "the exponents of the argument must be even");
                return any_ty();
            }
        }
        var r = same(t0, v);
        for i in 7 { r.dim[i] = t0.dim[i] / 2; }
        r.base = B_REAL;
        return r;
    }
    if name == "exp" || name == "ln" || name == "sin" || name == "cos" || name == "tan" || name == "tanh" {
        guard arity(tok, name, a, 1, 1) else { return any_ty(); }
        want_dimless(tok, name, a, 0);
        var r = real_scalar(a);
        r.num = t0.num;
        return r;
    }
    if name == "atan" {
        guard arity(tok, name, a, 1, 1) else { return any_ty(); }
        want_dimless(tok, name, a, 0);
        var r = real_scalar(a);
        r.base = B_ANGLE;
        return r;
    }
    if name == "atan2" {
        guard arity(tok, name, a, 2, 2) else { return any_ty(); }
        sum_ty(tok, t0, t1, ",");
        var r = real_scalar(a);
        r.base = B_ANGLE;
        return r;
    }
    if name == "min" || name == "max" {
        guard arity(tok, name, a, 2, MAX_ARGS) else { return any_ty(); }
        if a.count == 2 && t1.kind == K_GEOM {
            // max(f, r): the largest value of `f` on region `r` (spec §14.5).
            return same(t0, false);
        }
        var r = t0;
        for i in 1..a.count { r = sum_ty(first_tok(a.node[i]), r, a.ty[i], ","); }
        return same(r, v);
    }
    if name == "diff" {
        guard arity(tok, name, a, 2, 2) && !is_any(t0) && !is_any(t1) else { return any_ty(); }
        var r = num_ty(t0.dim - t1.dim, t1.rank);
        r.varying = v;
        return r;
    }
    // Differential operators (spec §7.3).
    if name == "grad" || name == "div" || name == "curl" || name == "lap" || name == "ddt" || name == "ddt2" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) else { return any_ty(); }
        guard is_num(t0) else {
            error(tok, "type_mismatch", str(name, " takes a quantity, not ", ty_str(t0)), "give a field");
            return any_ty();
        }
        var r = same(t0, true);
        if name == "grad" {
            if t0.rank >= 2 {
                error(tok, "rank_mismatch", "grad of a Tensor2 is not a type of Skein", "take grad of a scalar or a vector");
                return any_ty();
            }
            r.rank = t0.rank + 1;
            r.dim = t0.dim - LENGTH;
        }
        if name == "div" {
            if t0.rank == 0 {
                error(tok, "rank_mismatch", "div of a scalar", "take div of a vector or a Tensor2");
                return any_ty();
            }
            r.rank = if t0.rank == 4 { 2 } else { t0.rank - 1 };
            r.dim = t0.dim - LENGTH;
        }
        if name == "curl" {
            if t0.rank != 1 {
                error(tok, "rank_mismatch", "curl takes a vector", "take curl of a vector");
                return any_ty();
            }
            r.dim = t0.dim - LENGTH;
        }
        if name == "lap" { r.dim = t0.dim - dim_pow(LENGTH, 2); }
        if name == "ddt" { r.dim = t0.dim - TIME; }
        if name == "ddt2" { r.dim = t0.dim - dim_pow(TIME, 2); }
        return r;
    }
    // Shells and beams (spec §4.3).
    if name == "along" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        var r = same(t0, true);
        r.dim = t0.dim - LENGTH;
        return r;
    }
    if name == "tangent" {
        guard arity(tok, name, a, 2, 2) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        if t0.rank != 1 && t0.rank != 2 {
            error(tok, "rank_mismatch", str("tangent takes a vector or a Tensor2, not ", ty_str(t0)),
                  "give tangent a vector or a Tensor2 field");
            return any_ty();
        }
        want_rank(tok, name, a, 1, 1);
        want_dimless(tok, name, a, 1);
        return same(t0, v);
    }
    if name == "ply_stiffness" {
        guard arity(tok, name, a, 4, 4) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        for i in 1..3 { sum_ty(first_tok(a.node[i]), t0, a.ty[i], ","); }
        want_dimless(tok, name, a, 3);
        var r = num_ty(t0.dim, 4);
        r.varying = v;
        return r;
    }
    if name == "laminate" {
        guard arity(tok, name, a, 8, 8) && !is_any(t0) && want_rank(tok, name, a, 0, 4) else { return any_ty(); }
        let k, known = const_value(a.node[7]);
        guard known && (k == 0.0 || k == 1.0 || k == 2.0) else {
            error(first_tok(a.node[7]), "type_mismatch", "the moment of laminate is not 0, 1 or 2",
                  "give 0 for A, 1 for B or 2 for D");
            return any_ty();
        }
        want_length(name, a.ty[3], a.node[3], 3);
        want_length(name, a.ty[4], a.node[4], 4);
        var r = num_ty(t0.dim + dim_pow(LENGTH, k as i64 + 1), 4);
        r.varying = v;
        return r;
    }
    if name == "argmax" {
        // argmax(f, r): the point of the largest value of `f` on region `r` (spec §14.5).
        guard arity(tok, name, a, 2, 2) && want_geom(tok, name, a, 1) && !is_any(t0) else { return any_ty(); }
        if is_num(t0) && t0.rank != 0 {
            error(tok, "rank_mismatch", str("argmax takes a scalar, not ", ty_str(t0)), "take argmax of a scalar, for example norm(u)");
        }
        return num_ty(LENGTH, 1);
    }
    // Complex parts (gap G7).
    if name == "re" || name == "imag" || name == "arg" || name == "conj" {
        guard arity(tok, name, a, 1, 1) && !is_any(t0) && is_num(t0) else { return any_ty(); }
        var r = same(t0, v);
        if name != "conj" { r.num = NS_REAL; }
        if name == "arg" {
            r.dim = NO_DIM;
            r.base = B_ANGLE;
        }
        return r;
    }
    // Reports (spec §14.5) and eigenvalues (gap G6).
    if name == "reaction" {
        guard arity(tok, name, a, 1, 1) && want_geom(tok, name, a, 0) else { return any_ty(); }
        var d = FORCE;
        if space == 2 && !thickness { d = d - LENGTH; }
        return num_ty(d, 1);
    }
    if name == "integral" || name == "mean" {
        guard arity(tok, name, a, 2, 2) && want_geom(tok, name, a, 1) && !is_any(t0) else { return any_ty(); }
        var r = same(t0, false);
        if name == "integral" && t1.kind == K_GEOM { r.dim = t0.dim + dim_pow(LENGTH, t1.geo + depth()); }
        return r;
    }
    if name == "eigenvalue" {
        guard arity(tok, name, a, 1, 1) else { return any_ty(); }
        if !is_any(t0) && !is_int(t0) {
            error(first_tok(a.node[0]), "type_mismatch", str("the mode number is ", ty_str(t0), ", not Int"),
                  "give the number of the mode, for example eigenvalue(1)");
        }
        let proc = if in_expect { last_proc } else { cur_proc };
        if proc == "modal" { return num_ty(dim_pow(TIME, -1), 0); }
        if proc == "buckling" { return base_ty(B_REAL); }
        error(tok, "type_mismatch", "eigenvalue needs a modal or a buckling step",
              "use eigenvalue in a step with modal or buckling, or in an expect after one");
        return any_ty();
    }
    if name == "results" {
        guard arity(tok, name, a, 1, 2) else { return any_ty(); }
        return results_ty(tok, a.node[0], if a.count > 1 { a.node[1] } else { -1 });
    }
    // Geometry (spec §8).
    if name == "boundary" {
        guard arity(tok, name, a, 1, 2) && want_geom(tok, name, a, 0) && !is_any(t0) else { return any_ty(); }
        if a.count == 2 && !is_any(t1) && t1.kind != K_BOOL && t1.kind != K_GEOM {
            error(first_tok(a.node[1]), "type_mismatch",
                  str("the selection of boundary is ", ty_str(t1), ", not a Bool or a geometry value"),
                  "select with a condition such as x == 0 mm, or with a geometry value");
        }
        if t0.geo == 0 {
            error(tok, "geometry_dim", "a Point has no boundary", "take the boundary of a Curve, Surface or Solid");
            return any_ty();
        }
        return geom_ty(t0.geo - 1);
    }
    if name == "box" || name == "rect" {
        if name == "box" { space_3d(tok, name); } else { space_2d(tok, name); }
        guard arity(tok, name, a, 2, 2) else { return any_ty(); }
        want_point(name, t0, a.node[0], 0);
        want_point(name, t1, a.node[1], 1);
        return geom_ty(if name == "box" { 3 } else { 2 });
    }
    if name == "cylinder" || name == "cone" || name == "torus" {
        space_3d(tok, name);
        guard arity(tok, name, a, 4, 5) else { return any_ty(); }
        want_point(name, t0, a.node[0], 0);
        for i in 2..a.count { want_length(name, a.ty[i], a.node[i], i); }
        return geom_ty(3);
    }
    if name == "sphere" || name == "circle" {
        if name == "sphere" { space_3d(tok, name); } else { space_2d(tok, name); }
        guard arity(tok, name, a, 2, 2) else { return any_ty(); }
        want_point(name, t0, a.node[0], 0);
        want_length(name, t1, a.node[1], 1);
        return geom_ty(if name == "sphere" { 3 } else { 2 });
    }
    if name == "polygon" {
        space_2d(tok, name);
        guard arity(tok, name, a, 3, 100000) else { return any_ty(); }
        for i in min(a.count, MAX_ARGS) { want_point(name, a.ty[i], a.node[i], i); }
        return geom_ty(2);
    }
    if name == "extrude" || name == "revolve" || name == "sweep" || name == "loft" {
        space_3d(tok, name);
        guard arity(tok, name, a, 2, 3) else { return any_ty(); }
        return geom_ty(3);
    }
    if name == "line" || name == "arc" || name == "spline" {
        guard arity(tok, name, a, 2, 100000) else { return any_ty(); }
        if name != "arc" {
            for i in min(a.count, MAX_ARGS) { want_point(name, a.ty[i], a.node[i], i); }
        }
        return geom_ty(1);
    }
    if name == "point" {
        guard arity(tok, name, a, 1, 1) else { return any_ty(); }
        want_point(name, t0, a.node[0], 0);
        return geom_ty(0);
    }
    if name == "shell" || name == "pattern" || name == "pattern_polar" || name == "fillet" ||
       name == "chamfer" || name == "move" || name == "rotate" || name == "scale" {
        guard arity(tok, name, a, 2, 3) && want_geom(tok, name, a, 0) else { return any_ty(); }
        // A pattern keeps its copies: max(f, r) on it gives one value for
        // each copy (spec §14.5).
        var r = t0;
        if name == "pattern" || name == "pattern_polar" { r.list = true; }
        return r;
    }
    return any_ty();
}
