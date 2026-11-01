// The rules of steps, expectations, `allow`, scales, studies and the other
// items (spec §12.4 to §22).
import std;
import dictionary;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .expr;
import .typeexpr;
import .items;
import .physics;

// --- procedures (spec §14.2) ----------------------------------------------------------------

// The arguments of each procedure. The first one is also the positional one.
let STATIC_ARGS = ["ramp"];
let TRANSIENT_ARGS = ["end", "dt", "scheme"];
let EXPLICIT_SCHEMES = ["ssprk3", "central_difference"];
let HARMONIC_ARGS = ["freq", "damping"];
let MODAL_ARGS = ["count", "shift"];
let BUCKLING_ARGS = ["count"];
let PROCEDURES = ["static", "transient", "harmonic", "modal", "buckling"];

fn check_procedure(pr: i64) {
    let name = tok_text(syntax::nodes[pr].tok);
    guard contains_name(PROCEDURES, name) else {
        error(syntax::nodes[pr].tok, "procedure_args", str("'", name, "' is not a procedure"),
              str("use one of: ", join(PROCEDURES)));
        return;
    }
    let valid = if name == "static" { STATIC_ARGS[..] } else if name == "transient" { TRANSIENT_ARGS[..] }
                else if name == "harmonic" { HARMONIC_ARGS[..] } else if name == "modal" { MODAL_ARGS[..] } else { BUCKLING_ARGS[..] };
    var have_end = false;
    var have_dt = false;
    var have_freq = false;
    var stable = -1;
    var explicit = false;
    var a = syntax::nodes[pr].first;
    var i = 0;
    while a >= 0 {
        let k = syntax::nodes[a];
        let arg = key(if k.aux >= 0 { tok_text(k.aux) } else if i == 0 { valid[0] } else { "" });
        if !contains_name(valid, arg) {
            let msg: u8[>..] = if arg.len > 0 { str("'", arg, "' is not an argument of ", name) } else { str("only the first argument of ", name, " can have no name") };
            error(k.tok, "procedure_args", msg, str("the arguments of ", name, ": ", join(valid)));
            setting_ty(k.first);
        } else {
            proc_arg(name, arg, k.first);
        }
        if arg == "end" { have_end = true; }
        if arg == "dt" { have_dt = true; }
        if arg == "dt" && is_stable(k.first) { stable = k.first; }
        if arg == "scheme" && syntax::nodes[k.first].kind == syntax::N_NAME
           && contains_name(EXPLICIT_SCHEMES, tok_text(syntax::nodes[k.first].tok)) { explicit = true; }
        if arg == "freq" { have_freq = true; }
        i++;
        a = k.next;
    }
    if name == "transient" && (!have_end || !have_dt) {
        error(syntax::nodes[pr].tok, "procedure_args", str("transient needs ", if have_end { "dt" } else { "end" }),
              "write transient(end: <time>, dt: <time>)");
    }
    if stable >= 0 && !explicit {
        error(first_tok(stable), "procedure_args", "dt: stable needs an explicit scheme",
              "scheme: central_difference or ssprk3, or give dt a time");
    }
    if name == "harmonic" && !have_freq {
        error(syntax::nodes[pr].tok, "procedure_args", "harmonic needs freq", "write harmonic(freq: [<frequency>])");
    }
}

// `dt: stable`: the solver chooses the stable increment of an explicit
// scheme (§14.2).
fn is_stable(v: i64) -> bool {
    return syntax::nodes[v].kind == syntax::N_NAME && tok_text(syntax::nodes[v].tok) == "stable";
}

fn proc_arg(name: u8[:], arg: u8[:], v: i64) {
    if arg == "dt" && is_stable(v) { return; }
    if arg == "freq" {
        guard syntax::nodes[v].kind == syntax::N_LIST else {
            error(first_tok(v), "procedure_args", "freq is a list of frequencies", str("write freq: [", node_text(v), "]"));
            return;
        };
        let unit = syntax::child_of_kind(v, syntax::N_UNIT);
        var c = syntax::nodes[v].first;
        while c >= 0 {
            if syntax::nodes[c].kind != syntax::N_UNIT {
                var t = expr_ty(c);
                if unit >= 0 && is_bare(c) { t = bare_in_unit(t, unit); }
                want_constant(c, t, "a frequency");
                assign(first_tok(c), num_ty(dim_pow(TIME, -1), 0), t, "a frequency");
            }
            c = syntax::nodes[c].next;
        }
        return;
    }
    let t = setting_ty(v);
    want_constant(v, t, str("the argument ", arg));
    if arg == "end" || arg == "dt" { assign(first_tok(v), num_ty(TIME, 0), t, str("the argument ", arg)); }
    if arg == "ramp" || arg == "count" { assign(first_tok(v), base_ty(B_INT), t, str("the argument ", arg)); }
    if arg == "damping" { assign(first_tok(v), base_ty(B_REAL), t, str("the argument ", arg)); }
}

// --- steps (spec §14) ----------------------------------------------------------------------

fn check_step(n: i64) {
    let k = syntax::nodes[n];
    declare(S_STEP, k.tok, n, kind_ty(K_DECL));
    cur_step = step_count;
    step_count++;
    let pr = k.first;
    check_procedure(pr);
    cur_proc.clear();
    let pn = tok_text(syntax::nodes[pr].tok);
    cur_proc.append(pn[..min(pn.len, 32)]);
    cur_count = proc_count(pr);
    last_reports.clear();
    var c = syntax::nodes[pr].next;
    while c >= 0 {
        let s = syntax::nodes[c];
        match s.kind {
            syntax::N_BC => { check_bc(c); },
            syntax::N_RELEASE => { region_ty(s.first); },
            syntax::N_BIRTH, syntax::N_DEATH => {
                var r = s.first;
                while r >= 0 {
                    if syntax::nodes[r].kind == syntax::N_WHERE {
                        let t = expr_ty(syntax::nodes[r].first);
                        if !is_any(t) && t.kind != K_BOOL {
                            error(syntax::nodes[r].tok, "type_mismatch", str("the condition is ", ty_str(t)), "write a comparison");
                        }
                    } else {
                        region_ty(r);
                    }
                    r = syntax::nodes[r].next;
                }
            },
            syntax::N_SOLVE => { check_values(s.first); },
            syntax::N_EXPECT => { check_expect(c, 0); },
            syntax::N_OUTPUT => {
                var o = s.first;
                while o >= 0 {
                    if syntax::nodes[o].kind == syntax::N_EVERY {
                        let t = expr_ty(syntax::nodes[o].first);
                        assign(first_tok(syntax::nodes[o].first), base_ty(B_INT), t, "the count of every");
                    } else {
                        expr_ty(o);
                    }
                    o = syntax::nodes[o].next;
                }
            },
            syntax::N_REPORT => {
                var r = s.first;
                while r >= 0 {
                    in_report = true;
                    let t = expr_ty(syntax::nodes[r].first);
                    in_report = false;
                    if syntax::nodes[r].aux >= 0 && is_scalar(t) {
                        var nm: u8[..64] = [];
                        let rn = tok_text(syntax::nodes[r].aux);
                        nm.append(rn[..min(rn.len, 64)]);
                        last_reports.push(nm);
                    }
                    r = syntax::nodes[r].next;
                }
            },
            _ => {},
        }
        c = s.next;
    }
    last_proc = cur_proc;
    if cur_proc == "modal" || cur_proc == "buckling" {
        eig_proc = cur_proc;
        eig_count = cur_count;
    }
    cur_proc.clear();
    cur_count = -1;
    cur_step = -1;
}

// The literal `count` of a `modal` or `buckling` step, its first argument,
// or -1 when no literal gives it (spec §14.2).
fn proc_count(pr: i64) -> i64 {
    var a = syntax::nodes[pr].first;
    var i = 0;
    while a >= 0 {
        let k = syntax::nodes[a];
        if (k.aux >= 0 && tok_text(k.aux) == "count") || (k.aux < 0 && i == 0) {
            let v, ok = literal_value(k.first);
            return if ok { v as i64 } else { -1 };
        }
        i++;
        a = k.next;
    }
    return -1;
}

// The value of a literal number, which no param or --set can change.
fn literal_value(n: i64) -> f64, bool {
    guard syntax::nodes[n].kind == syntax::N_NUMBER else { return 0.0, false; }
    return const_value(n);
}

// --- expect (spec §15.3) -------------------------------------------------------------------

// An `expect` item, an `expect` of a step, or an `expect` of a study whose
// lists hold one value for each of its `cases` cases. `cases` is 0 where a
// list cannot stand, and -1 when the checker cannot count the cases.
fn check_expect(n: i64, cases: i64) {
    let k = syntax::nodes[n];
    in_expect = true;
    let lhs = k.first;
    let rhs = syntax::nodes[lhs].next;
    let within = syntax::nodes[rhs].next;
    let l = expr_ty(lhs);
    let rs = expect_values(rhs, cases);
    let ws = expect_values(within, cases);
    in_expect = false;
    for r in rs {
        guard is_any(l) || is_any(r) || (l.kind == K_NUM && r.kind == K_NUM) else {
            error(first_tok(lhs), "expect_form", str("an expect compares ", ty_str(l), " with ", ty_str(r)),
                  "compare two quantities, or a field with a value at a point with 'at'");
            break;
        };
        let joined_ty = sum_ty(first_tok(rhs), l, r, "==");
        for w in ws {
            if !w.percent && !is_any(joined_ty) && is_num(w) && !w.zero && w.dim != joined_ty.dim {
                error(first_tok(within), "unit_mismatch",
                      str("within is ", dim_str(w.dim), ", but the sides are ", dim_str(joined_ty.dim)),
                      str("give within in a unit of ", dim_str(joined_ty.dim), ", or as a percent"));
                break;
            }
        }
    }
    declare(S_EXPECT, k.tok, n, kind_ty(K_DECL));
}

// The types of the expected value or the `within` `v` of an expect: one for
// each value of a list, else one.
fn expect_values(v: i64, cases: i64) -> Ty[>..] {
    var tys: Ty[>..] = [];
    guard syntax::nodes[v].kind == syntax::N_LIST else {
        tys.push(expr_ty(v));
        return tys;
    };
    let unit = syntax::child_of_kind(v, syntax::N_UNIT);
    var c = syntax::nodes[v].first;
    while c >= 0 {
        if syntax::nodes[c].kind != syntax::N_UNIT {
            var t = expr_ty(c);
            if unit >= 0 && is_bare(c) { t = bare_in_unit(t, unit); }
            tys.push(t);
        }
        c = syntax::nodes[c].next;
    }
    if cases == 0 {
        error(syntax::nodes[v].tok, "expect_cases", "a list of values needs a grid or a zip study",
              "give one value, or put the expect in a grid or a zip study");
    } else if cases > 0 && tys.len != cases {
        error(syntax::nodes[v].tok, "expect_cases", str("the list has ", tys.len, " values, but the study has ", cases, " cases"),
              "give one value for each case, in the order of the cases");
    }
    return tys;
}

// --- allow (spec §19) ----------------------------------------------------------------------

// The ids that `allow` can turn off: lint rules and verification checks
// (docs/checks.md §3, §4).
let ALLOWABLE = [
    "unused_let", "unused_material", "unused_field", "load_no_cells", "bending_lock", "volumetric_lock",
    "dt_unstable", "stabilization_missing", "mesh_feature_small", "rigid_body", "pressure_level",
    "mesh_quality", "points_quality", "material_admissible", "residual", "equilibrium", "energy_balance", "conservation", "submodel_cut", "expect",
    "convergence",
];

fn check_allow(n: i64) {
    let k = syntax::nodes[n];
    let id = tok_text(k.tok);
    if !contains_name(ALLOWABLE, id) {
        let near = nearest_of(id, ALLOWABLE);
        let fix: u8[>..] = if near.len > 0 { str("did you mean '", near, "'?") } else { "see docs/checks.md §3 and §4 for the ids" };
        error(k.tok, "allow_reason", str("'", id, "' is not a lint rule or a verification check"), fix);
    }
    if k.aux >= 0 {
        let reason = tok_text(k.aux);
        var blank = true;
        for ch in reason[1..reason.len - 1] {
            if ch != ' ' { blank = false; }
        }
        if blank { error(k.aux, "allow_reason", "the reason of allow is empty", "give the reason in words"); }
    }
}

// --- scales (spec §21.2) -------------------------------------------------------------------

// The rank of the dimensions in `rows`, by Gaussian elimination.
fn dim_rank(rows: i64[7][:]) -> i64 {
    var m: f64[7][>..] = [];
    for r in rows {
        var row: f64[7] = [0.0; 7];
        for j in 7 { row[j] = r[j]; }
        m.push(row);
    }
    var rank = 0;
    for col in 7 {
        var piv = -1;
        for i in rank..m.len {
            if abs(m[i][col]) > 1e-9 { piv = i; break; }
        }
        if piv < 0 { continue; }
        let tmp = m[piv];
        m[piv] = m[rank];
        m[rank] = tmp;
        for i in 0..m.len {
            if i == rank { continue; }
            let f = m[i][col] / m[rank][col];
            for j in 7 { m[i][j] -= f * m[rank][j]; }
        }
        rank++;
    }
    return rank;
}

fn check_scales(n: i64) {
    if has_scales {
        error(syntax::nodes[n].tok, "duplicate_name", "the model has two scales lines", "keep one scales line");
    }
    has_scales = true;
    var c = syntax::nodes[n].first;
    while c >= 0 {
        let t = expr_ty(c);
        if is_num(t) {
            want_constant(c, t, "a scale");
            scale_dims.push(t.dim);
            if is_dimless(t.dim) || dim_rank(scale_dims) < scale_dims.len {
                error(first_tok(c), "scales_dependent",
                      str("the scale '", node_text(c), "' (", dim_str(t.dim), ") depends on the scales before it"),
                      str("remove '", node_text(c), "'"));
                scale_dims.pop();
            }
        } else if !is_any(t) {
            error(first_tok(c), "type_mismatch", str("a scale is ", ty_str(t)), "give a quantity");
        }
        c = syntax::nodes[c].next;
    }
}

// `scales_missing` and `scales_incomplete`, after all the items.
fn scales_end() {
    if nondim || used_dims.len == 0 { return; }
    if !has_scales {
        error(header_tok, "scales_missing", "the model has dimensions but no scales line",
              "add a scales line after the lets, for example 'scales L, E' with a length and a stress of the model");
        return;
    }
    let base = scale_dims.len;
    for d in used_dims {
        scale_dims.push(d);
        let r = dim_rank(scale_dims);
        scale_dims.pop();
        if r > base {
            var missing: u8[>..] = [];
            for i in 7 {
                var unit = NO_DIM;
                unit[i] = 1;
                scale_dims.push(unit);
                if dim_rank(scale_dims) > base {
                    if missing.len > 0 { missing.append(", "); }
                    missing.append(BASE_DIMS[i]);
                }
                scale_dims.pop();
            }
            error(header_tok, "scales_incomplete",
                  str("the scales do not make the dimension ", dim_str(d), " of the model"),
                  str("add a scale with the base dimension ", missing));
            return;
        }
    }
}

// --- study (spec §20.2) --------------------------------------------------------------------

let STUDIES = ["grid", "zip", "sample", "minimize", "topology"];
let TOPOLOGY_ARGS = ["objective", "region", "volume", "penalty", "filter", "max_iterations", "threshold"];
var objectives: i64[>..];

fn check_study(n: i64) {
    let k = syntax::nodes[n];
    let pr = k.first;
    let b = syntax::nodes[pr].next;
    let pname = tok_text(syntax::nodes[pr].tok);
    if !contains_name(STUDIES, pname) {
        error(syntax::nodes[pr].tok, "procedure_args", str("'", pname, "' is not a study procedure"),
              str("use one of: ", join(STUDIES)));
    }
    let optimizes = pname == "minimize" || pname == "topology";
    var a = syntax::nodes[pr].first;
    while a >= 0 {
        let v = syntax::nodes[a].first;
        let word = if syntax::nodes[a].aux >= 0 { tok_text(syntax::nodes[a].aux) } else { "objective" };
        if optimizes && word == "objective" && syntax::nodes[v].kind == syntax::N_NAME {
            objectives.push(syntax::nodes[v].tok);
        } else if pname == "minimize" && (word == "xtol" || word == "step") {
            fraction_arg(word, v);
        } else if pname == "topology" {
            topology_arg(syntax::nodes[a].aux, word, v);
        } else {
            setting_ty(v);
        }
        a = syntax::nodes[a].next;
    }
    if pname == "topology" { topology_needs(pr); }
    var c = syntax::nodes[b].first;
    while c >= 0 {
        let s = syntax::nodes[c];
        if s.kind == syntax::N_SETTING && pname == "topology" {
            error(s.tok, "procedure_args", "a topology study varies the density of each cell, not a param",
                  "remove the entry, or make a second study");
        } else if s.kind == syntax::N_SETTING {
            study_entry(s.tok, s.first);
        }
        if s.kind == syntax::N_EXPECT { check_expect(c, study_case_count(pname, b)); }
        c = s.next;
    }
}

// The argument `xtol` or `step` of minimize: a fraction of each range, a
// Real above 0 and at most 1 (femahoy gap G25).
fn fraction_arg(name: u8[:], v: i64) {
    let t = setting_ty(v);
    let x, known = const_value(v);
    if !(is_any(t) || (is_dimless_scalar(t) && (!known || (x > 0.0 && x <= 1.0)))) {
        error(first_tok(v), "procedure_args", str(name, " of minimize is a fraction of each range, above 0 and at most 1"),
              str("write ", name, ": ", if name == "xtol" { "0.001" } else { "0.1" }));
    }
}

// An argument of topology (femahoy gap G1): the region of the design, the
// volume fraction above 0 and at most 1, the SIMP power of 1 or more, the
// filter radius, a Length above 0, and the largest number of iterations.
fn topology_arg(at: i64, word: u8[:], v: i64) {
    let where = if at >= 0 { at } else { first_tok(v) };
    if !contains_name(TOPOLOGY_ARGS, word) {
        error(where, "procedure_args", str("'", word, "' is not an argument of topology"), str("use: ", join(TOPOLOGY_ARGS)));
        return;
    }
    if word == "objective" {
        error(first_tok(v), "procedure_args", "the objective of topology names a report", "write objective: compliance");
        return;
    }
    if word == "region" {
        let t = region_ty(v);
        if !is_any(t) && t.geo != 2 {
            error(first_tok(v), "procedure_args", "the region of topology is a Surface", "name the Surface of the mesh");
        }
        return;
    }
    let t = setting_ty(v);
    let x, known = const_value(v);
    let length, _ = named_dim("Length");
    var good = is_any(t) || !known;
    if word == "volume" { good = good || is_dimless_scalar(t) && x > 0.0 && x <= 1.0; }
    if word == "penalty" { good = good || is_dimless_scalar(t) && x >= 1.0; }
    if word == "filter" { good = good || t.kind == K_NUM && t.rank == 0 && t.dim == length && x > 0.0; }
    if word == "max_iterations" { good = good || is_int(t) && x >= 1.0; }
    if word == "threshold" { good = good || is_dimless_scalar(t) && x > 0.0 && x < 1.0; }
    if !good {
        let fix = if word == "volume" { "write volume: 0.5" } else if word == "penalty" { "write penalty: 3" }
                  else if word == "filter" { "write filter: 2.4 mm" }
                  else if word == "threshold" { "write threshold: 0.5" } else { "write max_iterations: 200" };
        error(first_tok(v), "procedure_args", str("the ", word, " of topology is not valid"), fix);
    }
}

// The arguments that topology needs: the objective, the region, the
// volume and the filter.
fn topology_needs(pr: i64) {
    for word in ["objective", "region", "volume", "filter"] {
        var found = false;
        var a = syntax::nodes[pr].first;
        while a >= 0 {
            let aux = syntax::nodes[a].aux;
            if (aux >= 0 && tok_text(aux) == word) || (aux < 0 && word == "objective") { found = true; }
            a = syntax::nodes[a].next;
        }
        if !found {
            error(syntax::nodes[pr].tok, "procedure_args", str("topology needs the argument ", word),
                  "write topology(objective: compliance, region: beam, volume: 0.5, filter: 2.4 mm)");
        }
    }
}

// The number of cases of a `grid` or a `zip` study with the block `b`, -1
// when a range makes it unknown, and 0 for the other procedures.
fn study_case_count(pname: u8[:], b: i64) -> i64 {
    if pname != "grid" && pname != "zip" { return 0; }
    var cases = 1;
    var c = syntax::nodes[b].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_SETTING {
            let v = syntax::nodes[c].first;
            if syntax::nodes[v].kind == syntax::N_RANGE { return -1; }
            var count = 1;
            if syntax::nodes[v].kind == syntax::N_LIST {
                count = 0;
                var i = syntax::nodes[v].first;
                while i >= 0 {
                    if syntax::nodes[i].kind != syntax::N_UNIT { count++; }
                    i = syntax::nodes[i].next;
                }
            }
            // A zip takes value i of each list, so its lists have one length.
            if pname == "zip" { return count; }
            cases *= count;
        }
        c = syntax::nodes[c].next;
    }
    return cases;
}

fn study_entry(tok: i64, v: i64) {
    let name = tok_text(tok);
    let s = visible.get(key(name));
    guard s && syms[s].kind == S_PARAM else {
        if s {
            error(tok, "study_param", str("'", name, "' is not a param"), str("change 'let ", name, "' to 'param ", name, "'"));
        } else {
            unknown_name(tok);
        }
        return;
    };
    let decl = syms[s].ty;
    let range = param_range_node(s);
    if syntax::nodes[v].kind == syntax::N_RANGE {
        let lo = syntax::nodes[v].first;
        let hi = syntax::nodes[lo].next;
        study_value(lo, -1, decl, range, name);
        study_value(hi, -1, decl, range, name);
        return;
    }
    guard syntax::nodes[v].kind == syntax::N_LIST else {
        study_value(v, -1, decl, range, name);
        return;
    };
    let unit = syntax::child_of_kind(v, syntax::N_UNIT);
    var c = syntax::nodes[v].first;
    while c >= 0 {
        if syntax::nodes[c].kind != syntax::N_UNIT { study_value(c, unit, decl, range, name); }
        c = syntax::nodes[c].next;
    }
}

fn study_value(c: i64, unit: i64, decl: Ty, range: i64, name: u8[:]) {
    var t = expr_ty(c);
    var v, ok = const_value(c);
    if unit >= 0 && is_bare(c) {
        t = bare_in_unit(t, unit);
        let ut, scale = unit_ty(unit);
        v *= scale;
    }
    assign(first_tok(c), decl, t, str("a value of '", name, "'"));
    if ok && range >= 0 { in_range(first_tok(c), v, range, name); }
}

fn objectives_end() {
    for tok in objectives {
        var found = false;
        for r in last_reports {
            if r == tok_text(tok) { found = true; }
        }
        if !found {
            var names: u8[>..] = [];
            for r, i in last_reports {
                if i > 0 { names.append(", "); }
                names.append(r);
            }
            let fix: u8[>..] = if names.len > 0 { str("use one of: ", names) } else { str("add 'report ", tok_text(tok), " = …' to the last step") };
            error(tok, "study_objective", str("'", tok_text(tok), "' is not a named scalar report of the last step"), fix);
        }
    }
}

// --- other items (spec §12.4, §12.6, §10.7) -----------------------------------------------

let GROUPS = ["bloch"];

fn check_pair(n: i64) {
    let a = syntax::nodes[n].first;
    let b = syntax::nodes[a].next;
    region_ty(a);
    region_ty(b);
    check_values(syntax::nodes[b].next);
}

fn check_symmetry(n: i64) {
    let k = syntax::nodes[n];
    let pr = k.first;
    let g = tok_text(syntax::nodes[pr].tok);
    if !contains_name(GROUPS, g) {
        error(syntax::nodes[pr].tok, "unknown_name", str("'", g, "' is not a symmetry group"),
              str("use one of: ", join(GROUPS)));
    }
    check_values(syntax::nodes[pr].next);
    declare(S_SYMMETRY, k.tok, n, kind_ty(K_DECL));
}
