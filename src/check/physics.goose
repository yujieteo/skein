// The rules of physics blocks, templates, assumptions, conditions and
// methods (spec §4.2, §10 to §13).
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

// The node of the template `name` in the standard library, or -1.
var std_prog = -1;

fn find_template(name: u8[:]) -> i64 {
    if std_prog < 0 { return -1; }
    var c = syntax::nodes[std_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_TEMPLATE && tok_text(syntax::nodes[c].tok) == name { return c; }
        c = syntax::nodes[c].next;
    }
    return -1;
}

// --- assumptions (spec §11.1) --------------------------------------------------------------

let ASSUMPTIONS = [
    "small_strain", "finite_strain", "isotropic", "orthotropic", "plane_stress", "plane_strain",
    "incompressible", "compressible", "newtonian", "laminar", "turbulent", "boussinesq", "isothermal",
];
let CONFLICTS = [
    "small_strain", "finite_strain", "isotropic", "orthotropic", "plane_stress", "plane_strain",
    "incompressible", "compressible", "laminar", "turbulent", "boussinesq", "compressible",
];

// The assumptions of the physics that the checker is in.
var assumed: Name[>..];

fn check_assume(n: i64) {
    var c = syntax::nodes[n].first;
    while c >= 0 {
        let a = syntax::nodes[c];
        let name = tok_text(a.tok);
        var arg = a.first;
        while arg >= 0 {
            setting_ty(arg);
            arg = syntax::nodes[arg].next;
        }
        if !contains_name(ASSUMPTIONS, name) {
            let near = nearest_of(name, ASSUMPTIONS);
            let fix: u8[>..] = if near.len > 0 { str("did you mean '", near, "'?") } else { str("use one of: ", join(ASSUMPTIONS)) };
            error(a.tok, "assumption_unknown", str("'", name, "' is not an assumption (spec §11.1)"), fix);
        } else {
            for i in 0..CONFLICTS.len / 2 {
                let x = CONFLICTS[2 * i];
                let y = CONFLICTS[2 * i + 1];
                let other = if name == x { y } else if name == y { x } else { "" };
                if other.len == 0 { continue; }
                for b in assumed {
                    if b == other {
                        error(a.tok, "assumption_conflict", str("'", name, "' conflicts with '", other, "'"),
                              str("remove '", name, "' or '", other, "'"));
                    }
                }
            }
            if (name == "plane_stress" || name == "plane_strain") && (space != 2 || axi) {
                error(a.tok, "assumption_space", str("'", name, "' needs a 2d model"),
                      "write the model as 2d, or remove the assumption");
            }
            // Only an explicit thickness makes integrals per body, not per depth.
            if (name == "plane_stress" || name == "plane_strain") && a.first >= 0 { thickness = true; }
            if name == "incompressible" { physics[cur_phys].incompressible = true; }
            assumed.push(key(name));
        }
        c = a.next;
    }
}

fn nearest_of(name: u8[:], names: const u8[:][:]) -> u8[>..] {
    var best: u8[>..] = [];
    var best_d = max(2, name.len / 3) + 1;
    for k in names {
        let d = edit_distance(name, k);
        if d < best_d { best_d = d; best.clear(); best.append(k); }
    }
    return best;
}

// --- physics (spec §10) ----------------------------------------------------------------------

// Records `div(F)` in the equation of unknown `k` (spec §12.1).
recursive fn note_fluxes(n: i64, k: i64) {
    let nd = syntax::nodes[n];
    if nd.kind == syntax::N_CALL {
        let callee = nd.first;
        let arg = syntax::nodes[callee].next;
        if syntax::nodes[callee].kind == syntax::N_NAME && tok_text(syntax::nodes[callee].tok) == "div" && arg >= 0 {
            fluxes.push(Flux { phys: cur_phys, node: syntax::nodes[arg].first, unknown: k });
        }
    }
    var c = nd.first;
    while c >= 0 {
        note_fluxes(c, k);
        c = syntax::nodes[c].next;
    }
}

// Checks the statement `s` of a physics or a template. `eqs` counts the
// equations.
fn check_phys_stmt(s: i64, eqs: i64&) {
    let k = syntax::nodes[s];
    match k.kind {
        syntax::N_ASSUME => { check_assume(s); },
        syntax::N_LET => {
            let sym = check_value_decl(s, S_LET);
            syms[sym].phys = cur_phys;
        },
        syntax::N_STATE => {
            let decl = type_of(k.first, true);
            let value = syntax::nodes[k.first].next;
            assign(first_tok(value), decl, expr_ty(value), str("'", tok_text(k.tok), "'"));
            use_dim(decl);
            var t = decl;
            t.varying = true;
            declare(S_STATE, k.tok, s, t);
        },
        syntax::N_EQN => {
            let l = expr_ty(k.first);
            let r = expr_ty(syntax::nodes[k.first].next);
            sum_ty(k.tok, l, r, "=");
            let p = physics[cur_phys];
            if eqs < p.count { note_fluxes(s, p.first + eqs); }
            note_high_derivs(s, 0);
            if convection(k.first, 0) || convection(syntax::nodes[k.first].next, 0) { physics[cur_phys].convective = true; }
            eqs++;
        },
        syntax::N_PLASTIC => {
            physics[cur_phys].plastic = true;
            let b = k.first;
            check_block(b, ["yield", "strain", "hardening"], "type_mismatch", "plastic");
            let y = setting(b, "yield");
            if y >= 0 {
                let t = expr_ty(y);
                let stress = num_ty(FORCE - dim_pow(LENGTH, 2), 0);
                sum_ty(first_tok(y), t, stress, "the yield function and a stress");
            }
        },
        _ => {},
    }
}

fn check_physics(n: i64) {
    let k = syntax::nodes[n];
    let psym = declare(S_PHYSICS, k.tok, n, kind_ty(K_DECL));
    physics.push(Phys { sym: psym, node: n, region: -1, template: -1, first: syms.len, count: 0,
                        has_method: false, plastic: false, incompressible: false, convective: false });
    cur_phys = physics.len - 1;
    assumed.clear();
    var uses = -1;
    var rty = any_ty();
    var c = k.first;
    while c >= 0 {
        let ck = syntax::nodes[c];
        if ck.kind == syntax::N_ON {
            physics[cur_phys].region = ck.first;
            rty = region_ty(ck.first);
        }
        if ck.kind == syntax::N_USES { uses = c; }
        c = ck.next;
    }
    // The unknowns.
    physics[cur_phys].first = syms.len;
    c = k.first;
    while c >= 0 {
        let u = syntax::nodes[c];
        if u.kind == syntax::N_UNKNOWN {
            var t = type_of(u.first, true);
            use_dim(t);
            let init = syntax::nodes[u.first].next;
            if init >= 0 { assign(first_tok(init), t, expr_ty(init), str("the initial value of '", tok_text(u.tok), "'")); }
            t.varying = true;
            declare(S_UNKNOWN, u.tok, c, t);
            physics[cur_phys].count++;
        }
        c = u.next;
    }
    // The materials first: a `let` may use their properties.
    var whole = 0;
    var parts = 0;
    c = k.first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_MATASSIGN {
            var m = syntax::nodes[c].first;
            while m >= 0 {
                if syntax::nodes[m].first >= 0 { parts++; } else { whole++; }
                if whole > 1 || (whole == 1 && parts > 0) {
                    error(syntax::nodes[m].tok, "material_cover", "two materials cover one point of the region",
                          "give each material its part with 'on', for example 'material a on r1, b on r2'");
                }
                let mt = syntax::nodes[m].tok;
                let ms = visible.get(key(tok_text(mt)));
                if !ms {
                    unknown_name(mt);
                } else if syms[ms].kind != S_MATERIAL {
                    error(mt, "type_mismatch", str("'", tok_text(mt), "' is not a material"), "name a material");
                } else {
                    phys_mats.push(PhysMat { phys: cur_phys, mat: ms });
                }
                if syntax::nodes[m].first >= 0 { region_ty(syntax::nodes[m].first); }
                m = syntax::nodes[m].next;
            }
        }
        c = syntax::nodes[c].next;
    }
    var eqs = 0;
    if uses >= 0 {
        check_uses(uses, eqs);
        // A shell is on a Surface and a beam is on a Curve (spec §4.3).
        let ut = syntax::nodes[uses].tok;
        let want = if tok_text(ut) == "reissner_mindlin" { 2 } else if tok_text(ut) == "timoshenko" { 1 } else { -1 };
        if want >= 0 && rty.kind == K_GEOM && rty.geo != want {
            error(ut, "geometry_dim", str("'", tok_text(ut), "' needs a ", GEO_NAMES[want], ", but the region is a ", GEO_NAMES[rty.geo]),
                  str("put the physics on a ", GEO_NAMES[want]));
        }
    }
    c = k.first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_EQN && uses >= 0 {
            error(syntax::nodes[c].tok, "template_eq", "a physics that uses a template cannot add an eq",
                  "write the physics in full, from the source of the template (skein help std)");
        } else {
            check_phys_stmt(c, eqs);
        }
        c = syntax::nodes[c].next;
    }
    let p = physics[cur_phys];
    if p.incompressible && uses < 0 && !has_continuity(n) {
        error(k.tok, "assumption_invalid", "incompressible needs an eq div(v) = 0 or a pressure unknown",
              "add the eq div(v) = 0 with a pressure unknown, or remove incompressible");
    }
    if eqs != p.count && !(uses >= 0 && physics[cur_phys].template < 0) {
        if eqs < p.count {
            error(k.tok, "eq_count", str("the physics solves ", p.count, " unknowns but has ", eqs, " equations"),
                  str("add the eq of '", tok_text(syms[p.first + eqs].tok), "'"));
        } else {
            error(k.tok, "eq_count", str("the physics solves ", p.count, " unknowns but has ", eqs, " equations"),
                  "remove an eq, or add its unknown to solves");
        }
    }
    cur_phys = -1;
}

// A physics that uses a template: its unknowns, arguments and statements (spec §4.2).
fn check_uses(uses: i64, eqs: i64&) {
    let ut = syntax::nodes[uses].tok;
    let tm = find_template(tok_text(ut));
    guard tm >= 0 else {
        error(ut, "unknown_name", str("'", tok_text(ut), "' is not a template of the standard library"),
              "use linear_elastic or navier_stokes, or write the physics in full");
        return;
    }
    physics[cur_phys].template = tm;
    let p = physics[cur_phys];
    let mark = locals.len;
    // The parameters, with the arguments of the model.
    var params = 0;
    var c = syntax::nodes[tm].first;
    while c >= 0 && syntax::nodes[c].kind == syntax::N_TPARAM {
        let pt = type_of(syntax::nodes[c].first, false);
        locals.push(Local { name: key(tok_text(syntax::nodes[c].tok)), ty: pt });
        params++;
        c = syntax::nodes[c].next;
    }
    var a = syntax::nodes[uses].first;
    var i = 0;
    while a >= 0 {
        let at = syntax::nodes[a];
        var slot = if at.aux < 0 { i } else { -1 };
        if at.aux >= 0 {
            for j in params {
                if locals[mark + j].name == tok_text(at.aux) { slot = j; }
            }
        }
        let t = expr_ty(at.first);
        if slot < 0 || slot >= params {
            error(at.tok, "type_mismatch", str("the template '", tok_text(ut), "' has no parameter '", if at.aux >= 0 { tok_text(at.aux) } else { "" }, "'"),
                  "use a parameter of the template (skein help std)");
        } else {
            assign(first_tok(at.first), locals[mark + slot].ty, t, str("the parameter '", locals[mark + slot].name, "'"));
        }
        i++;
        a = at.next;
    }
    // The unknowns: the model names them, with the types of the template.
    var j = 0;
    while c >= 0 && syntax::nodes[c].kind == syntax::N_UNKNOWN {
        let tt = type_of(syntax::nodes[c].first, false);
        if j < p.count {
            let ms = syms[p.first + j];
            let mt = type_of(syntax::nodes[ms.node].first, false);
            if mt.kind != tt.kind || mt.dim != tt.dim || mt.rank != tt.rank || mt.num != tt.num {
                error(ms.tok, "type_mismatch",
                      str("the unknown '", tok_text(ms.tok), "' is ", ty_str(mt), ", but the template solves ", ty_str(tt)),
                      str("declare it as ", ty_str(tt)));
            }
            var lt = tt;
            lt.varying = true;
            locals.push(Local { name: key(tok_text(syntax::nodes[c].tok)), ty: lt });
        }
        j++;
        c = syntax::nodes[c].next;
    }
    if j != p.count {
        error(ut, "type_mismatch", str("the template '", tok_text(ut), "' solves ", j, " unknowns, not ", p.count),
              "name one unknown for each unknown of the template");
    }
    // The statements, with each error at the template name.
    report_tok = ut;
    while c >= 0 {
        check_phys_stmt(c, eqs);
        c = syntax::nodes[c].next;
    }
    report_tok = -1;
    locals.resize(mark);
}

// True when physics `n` has an eq `div(v) = …` of an unknown and a
// pressure unknown (spec §11.1, `incompressible`).
fn has_continuity(n: i64) -> bool {
    var div_eq = false;
    var c = syntax::nodes[n].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        if k.kind == syntax::N_EQN && syntax::nodes[k.first].kind == syntax::N_CALL {
            let callee = syntax::nodes[k.first].first;
            let arg = syntax::nodes[callee].next;
            if tok_text(syntax::nodes[callee].tok) == "div" && arg >= 0 && unknown_of(syntax::nodes[arg].first) >= 0 {
                div_eq = true;
            }
        }
        c = k.next;
    }
    let p = physics[cur_phys];
    var pressure = false;
    for j in p.count {
        let t = syms[p.first + j].ty;
        if is_num(t) && t.rank == 0 && t.dim == FORCE - dim_pow(LENGTH, 2) { pressure = true; }
    }
    return div_eq && pressure;
}

// The unknowns under two nested derivatives, for `order_too_low`.
struct HighDeriv { phys: i64, sym: i64, tok: i64 }
var high_derivs: HighDeriv[>..];

let DERIVS = ["grad", "div", "curl", "lap"];

// Walks `n`. `depth` counts the derivatives around it: `lap` counts 2.
recursive fn note_high_derivs(n: i64, depth: i64) {
    let k = syntax::nodes[n];
    var d = depth;
    if k.kind == syntax::N_CALL && syntax::nodes[k.first].kind == syntax::N_NAME {
        let name = tok_text(syntax::nodes[k.first].tok);
        if contains_name(DERIVS, name) { d += if name == "lap" { 2 } else { 1 }; }
    }
    if k.kind == syntax::N_NAME && depth > 2 {
        let s = unknown_of(n);
        if s >= 0 { high_derivs.push(HighDeriv { phys: cur_phys, sym: s, tok: k.tok }); }
    }
    var c = k.first;
    while c >= 0 {
        note_high_derivs(c, d);
        c = syntax::nodes[c].next;
    }
}

// `order_too_low`: a `cg` order under 2 for an unknown that the equations
// differentiate more than twice (spec §13).
fn order_too_low(p: i64, b: i64) {
    let order = setting(b, "order");
    guard order >= 0 else { return; };
    for h in high_derivs {
        if h.phys != p { continue; }
        var v = 2.0;
        if syntax::nodes[order].kind == syntax::N_RECORD {
            var c = syntax::nodes[order].first;
            while c >= 0 {
                let ov, ok = const_value(syntax::nodes[c].first);
                if ok && syntax::nodes[c].aux >= 0 && tok_text(syntax::nodes[c].aux) == tok_text(syms[h.sym].tok) { v = ov; }
                c = syntax::nodes[c].next;
            }
        } else {
            let ov, ok = const_value(order);
            if ok { v = ov; }
        }
        if v < 2.0 {
            error(first_tok(order), "order_too_low",
                  str("the equations differentiate '", tok_text(syms[h.sym].tok), "' more than twice, but its order is ", v),
                  str("give '", tok_text(syms[h.sym].tok), "' the order 2 or more"));
            return;
        }
    }
}

// --- boundary conditions (spec §12) ------------------------------------------------------

// The unknown symbol of name node `n`, or -1.
fn unknown_of(n: i64) -> i64 {
    if syntax::nodes[n].kind != syntax::N_NAME { return -1; }
    let s = visible.get(key(tok_text(syntax::nodes[n].tok)));
    if s && syms[s].kind == S_UNKNOWN { return s; }
    return -1;
}

// The unknown of the flux `f` (an expression), or -1.
fn flux_unknown(f: i64) -> i64 {
    let text = norm_text(f);
    for fl in fluxes {
        if norm_text(fl.node) == text { return fl.unknown; }
    }
    return -1;
}

// For `F*n` or `dot(F, n)`, the node of F, else -1.
fn flux_of(n: i64) -> i64 {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_BINARY && syntax::tokens[k.tok].kind == syntax::TK_STAR {
        let r = syntax::nodes[k.first].next;
        if syntax::nodes[r].kind == syntax::N_NAME && tok_text(syntax::nodes[r].tok) == "n" { return k.first; }
    }
    if k.kind == syntax::N_CALL {
        let callee = k.first;
        let a1 = syntax::nodes[callee].next;
        if a1 >= 0 && syntax::nodes[callee].kind == syntax::N_NAME && tok_text(syntax::nodes[callee].tok) == "dot" {
            let a2 = syntax::nodes[a1].next;
            if a2 >= 0 && syntax::nodes[syntax::nodes[a2].first].kind == syntax::N_NAME
                && tok_text(syntax::nodes[syntax::nodes[a2].first].tok) == "n" {
                return syntax::nodes[a1].first;
            }
        }
    }
    return -1;
}

// For a call `name(a)` with one argument, the node of `a`, else -1.
fn one_arg(n: i64, name: u8[:]) -> i64 {
    let k = syntax::nodes[n];
    guard k.kind == syntax::N_CALL && call_name(k.first) == name else { return -1; };
    let a = syntax::nodes[k.first].next;
    guard a >= 0 && syntax::nodes[a].next < 0 else { return -1; };
    return syntax::nodes[a].first;
}

// A rigid face, `bc r: u = rigid(p)` (spec §12.1, gap G17): `u` is a whole
// displacement, and `p` is a `Point`.
fn check_rigid(rhs: i64, l: Ty, comp: u8[:], natural: bool) {
    let p = one_arg(rhs, "rigid");
    let pt = if p >= 0 { expr_ty(p) } else { any_ty() };
    let point = p >= 0 && pt.kind == K_GEOM && pt.geo == 0;
    if natural || comp.len > 0 || !point || !is_num(l) || l.rank != 1 || l.dim != LENGTH {
        error(first_tok(rhs), "rigid_form", "rigid(p) needs a whole displacement on the left and a Point p",
              "write bc face: u = rigid(p), where u is a Vector<Length> and p is a Point");
        return;
    }
    var nm: u8[..64] = [];
    let text = norm_text(p);
    nm.append(text[..min(text.len, 64)]);
    rigid_points.push(nm);
}

fn is_rigid_point(region: i64) -> bool {
    let text = norm_text(region);
    for p in rigid_points {
        if p == text[..min(text.len, 64)] { return true; }
    }
    return false;
}

fn check_bc(n: i64) {
    let k = syntax::nodes[n];
    let region = k.first;
    var lhs = syntax::nodes[region].next;
    let rhs = syntax::nodes[lhs].next;
    let rt = region_ty(region);
    // `moment(F*n) = M` on the point of a rigid face is the moment M about it.
    let moment = one_arg(lhs, "moment") >= 0;
    if moment { lhs = one_arg(lhs, "moment"); }
    var unknown = unknown_of(lhs);
    var comp: Name = [];
    var natural = false;
    if unknown < 0 && !moment && syntax::nodes[lhs].kind == syntax::N_MEMBER {
        unknown = unknown_of(syntax::nodes[lhs].first);
        if unknown >= 0 { comp = key(tok_text(syntax::nodes[lhs].tok)); }
    }
    if unknown < 0 || moment {
        unknown = -1;
        let f = flux_of(lhs);
        if f >= 0 {
            unknown = flux_unknown(f);
            natural = unknown >= 0;
        }
    }
    guard unknown >= 0 else {
        in_bc = true;
        expr_ty(lhs);
        expr_ty(rhs);
        in_bc = false;
        error(first_tok(lhs), "bc_form", str("the left side '", node_text(lhs), "' is not an unknown or a flux"),
              "write an unknown (u = …), a component (u.x = …), or a flux F*n where div(F) is in an eq");
        return;
    };
    cur_phys = syms[unknown].phys;
    let l = expr_ty(lhs);
    let point = rt.kind == K_GEOM && rt.geo == 0;
    if moment {
        comp = key("moment");
        guard point && is_rigid_point(region) else {
            expr_ty(rhs);
            cur_phys = -1;
            error(first_tok(syntax::nodes[region].next), "rigid_form",
                  str("moment needs a rigid(", node_text(region), ") before it"),
                  str("write bc face: u = rigid(", node_text(region), ") before this bc, on the region of the point"));
            return;
        };
    }
    if one_arg(rhs, "rigid") >= 0 {
        check_rigid(rhs, l, comp, natural);
    } else {
        let r = expr_ty(rhs);
        if natural {
            var lt = l;
            if point && is_num(lt) {
                // Gap G16: on a Point, `F*n = P` is the force P.
                lt.dim = lt.dim + dim_pow(LENGTH, space - 1 + depth());
                // Gap G17: the moment is a scalar in 2D, and a vector in 3D.
                if moment {
                    lt.dim = lt.dim + LENGTH;
                    lt.rank = if space == 3 { 1 } else { 0 };
                }
            }
            sum_ty(k.tok, lt, r, "=");
        } else {
            assign(first_tok(rhs), l, r, str("'", node_text(lhs), "'"));
        }
    }
    cur_phys = -1;
    if rt.kind == K_GEOM && rt.geo > 0 && rt.geo >= space {
        error(first_tok(region), "bc_off_boundary",
              str("the region '", node_text(region), "' is a ", GEO_NAMES[rt.geo], ", not a part of a boundary"),
              str("give a part of the boundary, for example boundary(", node_text(region), ")"));
    }
    for b in bcs {
        if b.step != cur_step || b.unknown != unknown || !same_region(b.region, region) { continue; }
        if b.comp.len > 0 && comp.len > 0 && b.comp != comp { continue; }
        // A moment is on the turn of the point, not on its displacement.
        if (b.comp == "moment") != (comp == "moment") { continue; }
        let other = syntax::nodes[b.node];
        let orhs = syntax::nodes[syntax::nodes[other.first].next].next;
        // On the point of a rigid(p), a force with 0 in the component of an essential condition does not fight it.
        if point && is_rigid_point(region) && b.natural != natural {
            let ec = if natural { b.comp } else { comp };
            let v, ok = rhs_value(if natural { rhs } else { orhs }, ec);
            if ec.len > 0 && ok && v == 0.0 { continue; }
        }
        if b.natural != natural || (!natural && !same_value(orhs, b.comp, rhs, comp)) {
            error(k.tok, "bc_conflict",
                  str("this bc and the bc on line ", syntax::tokens[other.tok].line, " fight on '", node_text(region), "'"),
                  "remove one of the two lines");
        }
    }
    bcs.push(Bc { step: cur_step, region: region, unknown: unknown, comp: comp, natural: natural, node: n });
}

// The components of a vector, in the order of its literal (spec §5.2).
let COMPONENTS = ["x", "y", "z"];

fn same_region(a: i64, b: i64) -> bool { norm_text(a) == norm_text(b) }

// True when the right side `a` of a bc on component `ca` and the right
// side `b` of a bc on component `cb` ("" for all components) give one
// value on the components that they share.
fn same_value(a: i64, ca: u8[:], b: i64, cb: u8[:]) -> bool {
    let va, oka = rhs_value(a, cb);
    let vb, okb = rhs_value(b, ca);
    if oka && okb { return va == vb; }
    return ca != cb || norm_text(a) == norm_text(b);
}

// The constant value in SI of component `comp` of the right side `rhs`
// ("" for all of it), and true when it is a constant.
fn rhs_value(rhs: i64, comp: u8[:]) -> f64, bool {
    guard comp.len > 0 && syntax::nodes[rhs].kind == syntax::N_VECTOR else { return const_value(rhs); };
    let unit = syntax::child_of_kind(rhs, syntax::N_UNIT);
    var i = 0;
    var c = syntax::nodes[rhs].first;
    while c >= 0 && i < COMPONENTS.len {
        if syntax::nodes[c].kind != syntax::N_UNIT {
            if COMPONENTS[i] == comp {
                var v, ok = const_value(c);
                if unit >= 0 && is_bare(c) {
                    var affine = 0;
                    let d, s = unit_value(unit, affine);
                    v *= s;
                }
                return v, ok;
            }
            i++;
        }
        c = syntax::nodes[c].next;
    }
    return 0.0, false;
}

// --- the free boundary report (command_line.md §6.2) -------------------------------------

// True when a bc on unknown `u` before bc `i` has the region of bc `i`.
fn seen_before(u: i64, i: i64) -> bool {
    for j in i {
        if bcs[j].unknown == u && same_region(bcs[j].region, bcs[i].region) { return true; }
    }
    return false;
}

// True when a bc on unknown `u` and component `comp` ("" for all) has the region `region`.
fn has_bc(u: i64, region: i64, comp: u8[:]) -> bool {
    for b in bcs {
        if b.unknown == u && b.comp == comp && same_region(b.region, region) { return true; }
    }
    return false;
}

fn add_free(psym: i64, field: u8[:], region: u8[:], condition: u8[:]) {
    var f = Free { physics: [], field: [], region: [], condition: [] };
    let name = tok_text(syms[psym].tok);
    f.physics.append(name[..min(name.len, 64)]);
    f.field.append(field[..min(field.len, 64)]);
    f.region.append(region[..min(region.len, 512)]);
    f.condition.append(condition[..min(condition.len, 256)]);
    free.push(f);
}

// The text of expression `e`, in parentheses when it is not a name or a call.
fn grouped(e: i64) -> u8[>..] {
    let k = syntax::nodes[e].kind;
    if k == syntax::N_NAME || k == syntax::N_CALL { return node_text(e); }
    return str("(", node_text(e), ")");
}

// A region with a bc on some components of an unknown gets the natural
// condition on each other component (spec §12.3).
fn free_boundary() {
    for p in physics {
        if p.region < 0 { continue; }
        for j in p.count {
            let u = p.first + j;
            var flux = -1;
            for fl in fluxes {
                if fl.unknown == u { flux = fl.node; break; }
            }
            if flux < 0 { continue; }
            let name = tok_text(syms[u].tok);
            let ftext = grouped(flux);
            var region: u8[>..] = str("boundary(", node_text(p.region), ")");
            for b, i in bcs {
                if b.unknown != u || seen_before(u, i) { continue; }
                region.append(" - ");
                region.append(grouped(b.region));
            }
            add_free(p.sym, name, region, str(ftext, "*n = 0"));
            for b, i in bcs {
                if b.unknown != u || b.comp.len == 0 || seen_before(u, i) || has_bc(u, b.region, "") { continue; }
                for c in COMPONENTS[..space] {
                    if !has_bc(u, b.region, c) {
                        add_free(p.sym, str(name, ".", c), node_text(b.region), str("(", ftext, "*n).", c, " = 0"));
                    }
                }
            }
        }
    }
}

// --- methods (spec §13) --------------------------------------------------------------------

let METHOD_KINDS = ["cg", "dg", "fv", "rkpm", "gfd", "ode"];
let CG_SETTINGS = ["mesh", "order", "family", "stabilize", "quadrature", "couple"];
let DG_SETTINGS = ["mesh", "order", "flux", "penalty", "shock_capturing", "couple"];
let FV_SETTINGS = ["mesh", "gradient", "flux", "limiter", "pv_coupling", "couple"];
let RKPM_SETTINGS = ["points", "basis", "support", "integrate", "bc", "couple"];
let GFD_SETTINGS = ["points", "stencil", "weight", "couple"];
let ODE_SETTINGS = ["couple"];
let STABLE = ["supg_pspg", "vms", "supg_pspg_lsic", "pspg"];
let FAMILIES = ["lagrange", "serendipity", "incompatible", "nedelec", "raviart_thomas"];
let QUADRATURES = ["full", "reduced"];

fn check_method(n: i64) {
    let k = syntax::nodes[n];
    let name = tok_text(k.tok);
    let kind = if k.aux >= 0 { tok_text(k.aux) } else { "" };
    let b = k.first;
    var p = -1;
    let s = visible.get(key(name));
    if !s {
        unknown_name(k.tok);
    } else if syms[s].kind != S_PHYSICS {
        error(k.tok, "type_mismatch", str("'", name, "' is not a physics"), "name a physics");
    } else {
        for q, i in physics {
            if q.sym == s { p = i; }
        }
    }
    if p >= 0 {
        if physics[p].has_method {
            error(k.tok, "duplicate_name", str("the physics '", name, "' has two methods"), "remove one method");
        }
        physics[p].has_method = true;
    }
    guard contains_name(METHOD_KINDS, kind) else {
        error(k.aux, "method_setting", str("'", kind, "' is not a method kind"), str("use one of: ", join(METHOD_KINDS)));
        check_values(b);
        return;
    }
    let valid = if kind == "cg" { CG_SETTINGS[..] } else if kind == "dg" { DG_SETTINGS[..] } else if kind == "fv" { FV_SETTINGS[..] }
                else if kind == "rkpm" { RKPM_SETTINGS[..] } else if kind == "gfd" { GFD_SETTINGS[..] } else { ODE_SETTINGS[..] };
    check_block(b, valid, "method_setting", str("the method ", kind));
    let mesh = setting(b, "mesh");
    if mesh >= 0 && contains_name(valid, "mesh") { want_kind(mesh, K_MESH, "a mesh"); }
    let pts = setting(b, "points");
    if pts >= 0 && contains_name(valid, "points") { want_kind(pts, K_POINTS, "a point cloud"); }
    if p >= 0 && kind == "gfd" && physics[p].plastic {
        error(k.aux, "method_physics", "gfd cannot collocate a physics with a plastic block",
              "use cg or rkpm for this physics");
    }
    if p >= 0 && kind == "cg" && physics[p].incompressible { inf_sup(p, b, k.aux); }
    if p >= 0 && kind == "cg" { order_too_low(p, b); }
    if kind == "cg" { family_cell(b); }
}

// The word of setting `name` in block `b`, or `default` when the block does
// not have it. A value that is not in `words` is `method_setting`, and it
// gives "".
fn setting_word(b: i64, name: u8[:], words: const u8[:][:], default: u8[:]) -> u8[>..] {
    let v = setting(b, name);
    guard v >= 0 else { return str(default); };
    let w = node_text(v);
    guard syntax::nodes[v].kind == syntax::N_NAME && contains_name(words, w) else {
        error(first_tok(v), "method_setting", str("'", w, "' is not a value of ", name), str("use one of: ", join(words)));
        return [];
    };
    return w;
}

// `family_cell`: a serendipity cell needs order 2, an incompatible cell
// needs order 1 and full quadrature, and both need quad or hex cells
// (spec §13).
fn family_cell(b: i64) {
    let family = setting_word(b, "family", FAMILIES, "lagrange");
    let quadrature = setting_word(b, "quadrature", QUADRATURES, "full");
    guard family == "serendipity" || family == "incompatible" else { return; };
    let at = first_tok(setting(b, "family"));
    let want = if family == "serendipity" { 2.0 } else { 1.0 };
    let order = setting(b, "order");
    var orders: i64[>..] = [];
    if order >= 0 && syntax::nodes[order].kind == syntax::N_RECORD {
        var c = syntax::nodes[order].first;
        while c >= 0 {
            orders.push(syntax::nodes[c].first);
            c = syntax::nodes[c].next;
        }
    } else if order >= 0 {
        orders.push(order);
    }
    var bad = order < 0 && want != 2.0;
    for o in orders {
        let v, ok = const_value(o);
        if ok && v != want { bad = true; }
    }
    if bad {
        error(at, "family_cell", str("the family ", family, " needs order ", want as i64),
              str("write order: ", want as i64, ", or remove the family"));
    }
    if family == "incompatible" && quadrature == "reduced" {
        error(at, "family_cell", "the incompatible modes need full quadrature", "remove quadrature: reduced");
    }
    let ms = setting(b, "mesh");
    guard ms >= 0 && syntax::nodes[ms].kind == syntax::N_NAME else { return; };
    let s = visible.get(key(tok_text(syntax::nodes[ms].tok)));
    guard s && syms[s].kind == S_MESH else { return; };
    let cell = setting(syntax::child_of_kind(syms[s].node, syntax::N_BLOCK), "cell");
    let name = if cell >= 0 { node_text(cell) } else { "" };
    if name != "quad" && name != "hex" {
        error(at, "family_cell", str("the family ", family, " needs quad or hex cells, but the mesh '",
                                     node_text(ms), "' has ", if cell >= 0 { name } else { "the default cells" }),
              str("add cell: ", if space == 3 { "hex" } else { "quad" }, " to the mesh"));
    }
}

fn want_kind(n: i64, kind: i64, what: u8[:]) {
    let t = setting_ty(n);
    if !is_any(t) && t.kind != kind {
        error(first_tok(n), "type_mismatch", str("'", node_text(n), "' is not ", what), str("name ", what));
    }
}

// `inf_sup`: equal orders of velocity and pressure in `cg` with no
// pressure stabilization.
fn inf_sup(p: i64, b: i64, tok: i64) {
    let order = setting(b, "order");
    guard order >= 0 else { return; };
    let stab = setting(b, "stabilize");
    if stab >= 0 && contains_name(STABLE, node_text(stab)) { return; }
    let ph = physics[p];
    var vel = -1;
    var pre = -1;
    for j in ph.count {
        let t = syms[ph.first + j].ty;
        if is_num(t) && t.rank == 1 && t.dim == LENGTH - TIME { vel = ph.first + j; }
        if is_num(t) && t.rank == 0 && t.dim == FORCE - dim_pow(LENGTH, 2) { pre = ph.first + j; }
    }
    guard vel >= 0 && pre >= 0 else { return; };
    var ov = 0.0;
    var op = 0.0;
    if syntax::nodes[order].kind == syntax::N_RECORD {
        var c = syntax::nodes[order].first;
        while c >= 0 {
            let v, ok = const_value(syntax::nodes[c].first);
            if syntax::nodes[c].aux >= 0 && ok {
                let f = tok_text(syntax::nodes[c].aux);
                if f == tok_text(syms[vel].tok) { ov = v; }
                if f == tok_text(syms[pre].tok) { op = v; }
            }
            c = syntax::nodes[c].next;
        }
    } else {
        let v, ok = const_value(order);
        if ok { ov = v; op = v; }
    }
    if ov > 0.0 && ov == op {
        error(first_tok(order), "inf_sup", "equal orders of velocity and pressure with no pressure stabilization",
              "add stabilize: supg_pspg, or give the velocity a higher order: { v: 2, p: 1 }");
    }
}

// `method_missing`: a physics with a region and no method needs exactly
// one mesh on its region (spec §13). The checker has no geometry, so it
// counts the meshes on the same region, or all meshes when no mesh names it.
// A physics that takes the defaults goes to `chosen` for the check report.
fn method_missing() {
    for p, i in physics {
        if p.has_method || p.region < 0 { continue; }
        var all = 0;
        var same = 0;
        var any_mesh = -1;
        var same_mesh = -1;
        for s in syms {
            if s.kind != S_MESH { continue; }
            all++;
            any_mesh = s.tok;
            if same_region(syntax::nodes[s.node].first, p.region) {
                same++;
                same_mesh = s.tok;
            }
        }
        let name = tok_text(syms[p.sym].tok);
        if same == 1 || (same == 0 && all == 1) {
            choose(i, if same == 1 { same_mesh } else { any_mesh });
            continue;
        }
        let msg: u8[>..] = if all == 0 { str("the physics '", name, "' has no method and the model has no mesh") }
                           else { str("the physics '", name, "' has no method, and ", all, " meshes can hold its region") };
        error(syms[p.sym].tok, "method_missing", msg, str("add method ", name, ": cg { mesh: <name> }"));
    }
}

// Adds the default method of physics `p` on the mesh named by `mesh_tok` to
// `chosen` (spec §13).
fn choose(p: i64, mesh_tok: i64) {
    let ph = physics[p];
    var c = Chosen { physics: [], mesh: [], order: [], stabilize: [] };
    let name = tok_text(syms[ph.sym].tok);
    c.physics.append(name[..min(name.len, 64)]);
    c.mesh.append(tok_text(mesh_tok)[..min(tok_text(mesh_tok).len, 64)]);
    let first = if ph.count > 0 { default_order(ph.first) } else { 2 };
    var same = true;
    for j in ph.count {
        if default_order(ph.first + j) != first { same = false; }
    }
    var order: u8[>..] = [];
    if same {
        format(order, first);
    } else {
        order.append("{ ");
        for j in ph.count {
            if j > 0 { order.append(", "); }
            format(order, tok_text(syms[ph.first + j].tok), ": ", default_order(ph.first + j));
        }
        order.append(" }");
    }
    c.order.append(order[..min(order.len, 256)]);
    c.stabilize.append(if ph.convective { "supg" } else { "none" });
    chosen.push(c);
}

// The default `cg` order of unknown `u` (spec §13): 2, except 1 for a
// pressure that pairs with a velocity under `incompressible` (Taylor-Hood).
fn default_order(u: i64) -> i64 {
    let p = syms[u].phys;
    guard p >= 0 && physics[p].incompressible && is_pressure(syms[u].ty) else { return 2; };
    let ph = physics[p];
    for j in ph.count {
        if is_velocity(syms[ph.first + j].ty) { return 1; }
    }
    return 2;
}

fn is_velocity(t: Ty) -> bool { is_num(t) && t.rank == 1 && t.dim == LENGTH - TIME }

fn is_pressure(t: Ty) -> bool { is_num(t) && t.rank == 0 && t.dim == FORCE - dim_pow(LENGTH, 2) }

// True when the expression `n` has a convection term: `grad(a)*w`,
// `w*grad(a)`, `dot(w, grad(a))` or `dot(grad(a), w)`, with `w` a name of
// a vector of velocities and each side perhaps a product with other
// factors (spec §13, the default `stabilize`). A name of a physics `let`
// gives the convection terms of its value, to `depth` 16.
recursive fn convection(n: i64, depth: i64) -> bool {
    guard n >= 0 && depth <= 16 else { return false; };
    let k = syntax::nodes[n];
    if k.kind == syntax::N_NAME {
        let s = visible.get(key(tok_text(k.tok)));
        return s && syms[s].kind == S_LET && syms[s].phys == cur_phys && syms[s].node >= 0
               && convection(syntax::nodes[syms[s].node].last, depth + 1);
    }
    var a = -1;
    var b = -1;
    if k.kind == syntax::N_BINARY && syntax::tokens[k.tok].kind == syntax::TK_STAR {
        a = k.first;
        b = syntax::nodes[a].next;
    } else if k.kind == syntax::N_CALL && tok_text(syntax::nodes[k.first].tok) == "dot" {
        let x = syntax::nodes[k.first].next;
        if x >= 0 && syntax::nodes[x].next >= 0 {
            a = syntax::nodes[x].first;
            b = syntax::nodes[syntax::nodes[x].next].first;
        }
    }
    if a >= 0 && b >= 0 && ((factor(a, true) && factor(b, false)) || (factor(b, true) && factor(a, false))) { return true; }
    var c = k.first;
    while c >= 0 {
        if convection(c, depth) { return true; }
        c = syntax::nodes[c].next;
    }
    return false;
}

// True when the product `n` has a factor `grad(…)` (`grad` true) or a
// velocity name (`grad` false).
recursive fn factor(n: i64, grad: bool) -> bool {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_BINARY && syntax::tokens[k.tok].kind == syntax::TK_STAR {
        return factor(k.first, grad) || factor(syntax::nodes[k.first].next, grad);
    }
    return if grad { is_grad(n) } else { velocity_name(n) };
}

fn is_grad(n: i64) -> bool {
    let k = syntax::nodes[n];
    return k.kind == syntax::N_CALL && syntax::nodes[k.first].kind == syntax::N_NAME && tok_text(syntax::nodes[k.first].tok) == "grad";
}

// True when `n` is a name of a vector of velocities: a local of a template,
// or a symbol of the model.
fn velocity_name(n: i64) -> bool {
    guard syntax::nodes[n].kind == syntax::N_NAME else { return false; };
    let name = tok_text(syntax::nodes[n].tok);
    for i in 0..locals.len {
        let l = locals[locals.len - 1 - i];
        if l.name == name { return is_velocity(l.ty); }
    }
    let s = visible.get(key(name));
    return s && is_velocity(syms[s].ty);
}
