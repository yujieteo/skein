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
let TRANSIENT_ARGS = ["end", "dt", "scheme", "adapt"];
let HARMONIC_ARGS = ["freq"];
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
        if arg == "freq" { have_freq = true; }
        i++;
        a = k.next;
    }
    if name == "transient" && (!have_end || !have_dt) {
        error(syntax::nodes[pr].tok, "procedure_args", str("transient needs ", if have_end { "dt" } else { "end" }),
              "write transient(end: <time>, dt: <time>)");
    }
    if name == "harmonic" && !have_freq {
        error(syntax::nodes[pr].tok, "procedure_args", "harmonic needs freq", "write harmonic(freq: [<frequency>])");
    }
}

fn proc_arg(name: u8[:], arg: u8[:], v: i64) {
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
            syntax::N_SOLVE, syntax::N_ADAPT => { check_values(s.first); },
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
    cur_proc.clear();
    cur_step = -1;
}

// --- expect (spec §15.3) -------------------------------------------------------------------

fn check_expect(n: i64) {
    let k = syntax::nodes[n];
    in_expect = true;
    let lhs = k.first;
    let rhs = syntax::nodes[lhs].next;
    let within = syntax::nodes[rhs].next;
    let l = expr_ty(lhs);
    let r = expr_ty(rhs);
    let w = expr_ty(within);
    in_expect = false;
    if !is_any(l) && !is_any(r) && (l.kind != K_NUM || r.kind != K_NUM) {
        error(first_tok(lhs), "expect_form", str("an expect compares ", ty_str(l), " with ", ty_str(r)),
              "compare two quantities, or a field with a value at a point with 'at'");
    } else {
        let joined_ty = sum_ty(first_tok(rhs), l, r, "==");
        if !w.percent && !is_any(joined_ty) && is_num(w) && !w.zero && w.dim != joined_ty.dim {
            error(first_tok(within), "unit_mismatch",
                  str("within is ", dim_str(w.dim), ", but the sides are ", dim_str(joined_ty.dim)),
                  str("give within in a unit of ", dim_str(joined_ty.dim), ", or as a percent"));
        }
    }
    declare(S_EXPECT, k.tok, n, kind_ty(K_DECL));
}

// --- allow (spec §19) ----------------------------------------------------------------------

// The ids that `allow` can turn off: lint rules and verification checks
// (docs/checks.md §3, §4).
let ALLOWABLE = [
    "unused_let", "unused_material", "unused_field", "load_no_cells", "bending_lock", "volumetric_lock",
    "dt_unstable", "stabilization_missing", "mesh_feature_small", "rigid_body", "pressure_level",
    "mesh_quality", "points_quality", "material_admissible", "dfm_wall", "dfm_draft", "dfm_undercut",
    "dfm_inner_radius", "dfm_hole", "dfm_overhang", "dfm_bend", "dfm_stack", "residual",
    "newton_converged", "equilibrium", "energy_balance", "conservation", "submodel_cut", "expect",
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

let STUDIES = ["grid", "zip", "sample", "minimize"];
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
    var a = syntax::nodes[pr].first;
    while a >= 0 {
        let v = syntax::nodes[a].first;
        if pname == "minimize" && syntax::nodes[a].aux < 0 && syntax::nodes[v].kind == syntax::N_NAME {
            objectives.push(syntax::nodes[v].tok);
        } else {
            setting_ty(v);
        }
        a = syntax::nodes[a].next;
    }
    var c = syntax::nodes[b].first;
    while c >= 0 {
        let s = syntax::nodes[c];
        if s.kind == syntax::N_SETTING { study_entry(s.tok, s.first); }
        c = s.next;
    }
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

// --- other items (spec §12.4, §12.6, §10.7, §22) -------------------------------------------

let GROUPS = ["cyclic", "dihedral", "bloch"];
let PROCESSES = ["cnc_mill", "injection_mould", "casting", "sheet_metal", "additive", "turning"];
let CNC_SETTINGS = ["axes", "tool_radius", "max_depth_ratio", "approach"];
let MOULD_SETTINGS = ["pull", "draft", "min_wall", "max_wall"];
let CASTING_SETTINGS = ["pull", "draft", "min_wall"];
let SHEET_SETTINGS = ["thickness", "min_bend_radius", "min_flange"];
let ADDITIVE_SETTINGS = ["build", "max_overhang", "min_wall", "min_hole"];
let TURNING_SETTINGS = ["axis", "tool_radius"];

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

fn check_manufacture(n: i64) {
    let region = syntax::nodes[n].first;
    region_ty(region);
    let pr = syntax::nodes[region].next;
    let p = tok_text(syntax::nodes[pr].tok);
    let b = syntax::nodes[pr].next;
    guard contains_name(PROCESSES, p) else {
        error(syntax::nodes[pr].tok, "process_setting", str("'", p, "' is not a process"), str("use one of: ", join(PROCESSES)));
        check_values(b);
        return;
    };
    let valid = if p == "cnc_mill" { CNC_SETTINGS[..] } else if p == "injection_mould" { MOULD_SETTINGS[..] }
                else if p == "casting" { CASTING_SETTINGS[..] } else if p == "sheet_metal" { SHEET_SETTINGS[..] }
                else if p == "additive" { ADDITIVE_SETTINGS[..] } else { TURNING_SETTINGS[..] };
    check_block(b, valid, "process_setting", str("the process ", p));
}

fn check_tolerance(n: i64) {
    let v = syntax::nodes[n].first;
    let lo = syntax::nodes[v].next;
    let hi = syntax::nodes[lo].next;
    let t = expr_ty(v);
    want_constant(v, t, "the value of a tolerance");
    sum_ty(first_tok(lo), t, expr_ty(lo), "the value and its low limit");
    sum_ty(first_tok(hi), t, expr_ty(hi), "the value and its high limit");
    check_block(syntax::nodes[hi].next, ["method", "samples"], "type_mismatch", "a tolerance");
}
