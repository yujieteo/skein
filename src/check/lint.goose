// The lint rules of `skein lint` (docs/checks.md §3) that need no mesh.
// Each one is a warning on a model with no error. `allow <id>: "reason"` in
// the model turns a rule off. The rules that need the mesh, `bending_lock`
// and `mesh_feature_small`, are in src/mesher/lint.goose and run only with
// `--mesh`.
//
// `dt_unstable` estimates the stable limit of solids only, and
// `stabilization_missing` needs a velocity in the `scales` line: flow
// physics is phase 3 (docs/roadmap.md).
import std;
import dictionary;
import math;
namespace check;
import syntax from "../syntax";
import .dims;
import .types;
import .state;
import .expr;
import .items;

let EXPLICIT_SCHEMES = ["ssprk3", "central_difference"];

// The names that a part of the model uses.
var used_names = dictionary<Name, i64> {};
// The regions with no active cells in the current step, by their text.
var dead = dictionary<Name, i64> {};

fn lint_rules() {
    unused_let();
    unused_material();
    unused_field();
    load_no_cells();
    volumetric_lock();
    dt_unstable();
    stabilization_missing();
}

// Records the lint warning `id` at token `tok`, unless an `allow` turns it off.
fn warn(tok: i64, id: u8[:], msg: u8[:], fix: u8[:]) {
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_ALLOW && tok_text(syntax::nodes[c].tok) == id { return; }
        c = syntax::nodes[c].next;
    }
    add_diag(tok, id, msg, fix, true);
}

// Adds each name of the subtree of `n` to `used_names`, except in the
// subtree `skip`. With `targets` false, the left side of a `bc` is not a use.
recursive fn note_names(n: i64, skip: i64, targets: bool) {
    if n == skip { return; }
    let k = syntax::nodes[n];
    if k.kind == syntax::N_NAME { used_names.insert(key(tok_text(k.tok)), 1); }
    var c = k.first;
    var i = 0;
    while c >= 0 {
        if targets || k.kind != syntax::N_BC || i != 1 { note_names(c, skip, targets); }
        c = syntax::nodes[c].next;
        i++;
    }
}

// `unused_let`: a `let` of the model, at the top level or in a physics,
// that no name uses. A `let` of a `use` file is not a part of the model.
fn unused_let() {
    used_names.clear();
    note_names(model_prog, -1, true);
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        if k.kind == syntax::N_LET { unused_let_at(c); }
        if k.kind == syntax::N_PHYSICS {
            var s = k.first;
            while s >= 0 {
                if syntax::nodes[s].kind == syntax::N_LET { unused_let_at(s); }
                s = syntax::nodes[s].next;
            }
        }
        c = k.next;
    }
}

fn unused_let_at(n: i64) {
    let t = syntax::nodes[n].tok;
    if used_names.contains(key(tok_text(t))) { return; }
    warn(t, "unused_let", str("the let '", tok_text(t), "' has no use"), "remove it");
}

// `unused_material`: a material of the model that no physics assigns.
fn unused_material() {
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        var assigned = false;
        for pm in phys_mats {
            if syms[pm.mat].node == c { assigned = true; }
        }
        if k.kind == syntax::N_MATERIAL && !assigned {
            warn(k.tok, "unused_material", str("material ", tok_text(k.tok), " is not assigned"),
                 "assign it in a physics block or remove it");
        }
        c = k.next;
    }
}

// `unused_field`: no part of the model outside a physics uses a name of it
// (an unknown, a state or a `let`), other than as the left side of a `bc`.
// Then the solve of its unknowns gives nothing.
fn unused_field() {
    for p, i in physics {
        used_names.clear();
        note_names(model_prog, p.node, false);
        var used = false;
        for s in syms {
            if s.phys == i && s.tok >= 0 && used_names.contains(key(tok_text(s.tok))) { used = true; }
        }
        if used { continue; }
        let name = tok_text(syms[p.sym].tok);
        for j in p.count {
            let u = tok_text(syms[p.first + j].tok);
            warn(syms[p.first + j].tok, "unused_field",
                 str("the unknown '", u, "' of physics '", name, "' is solved but never output, reported or used"),
                 str("output it, for example 'output ", u, "' in a step, or remove the physics"));
        }
    }
}

// --- load_no_cells -----------------------------------------------------------------------

// `load_no_cells`: a `bc` whose region has no active cells in a step. A
// region is dead after a `death` and before its `birth`. A `birth` or a
// `death` with `where` kills only a part, so it is not counted.
fn load_no_cells() {
    dead.clear();
    // A region whose first event is a birth is dead before it.
    var seen = dictionary<Name, i64> {};
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_STEP {
            each_event(c) { r, birth =>
                let k = key(norm_text(r));
                if !seen.contains(k) {
                    seen.insert(k, 1);
                    if birth { dead.insert(k, 1); }
                }
            }
        }
        c = syntax::nodes[c].next;
    }
    var warned = dictionary<i64, i64> {};
    c = syntax::nodes[model_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_STEP {
            each_event(c) { r, birth =>
                if birth { dead.remove(key(norm_text(r))); } else { dead.insert(key(norm_text(r)), 1); }
            }
            var s = syntax::nodes[c].first;
            while s >= 0 {
                if syntax::nodes[s].kind == syntax::N_BC { dead_bc(s, c, warned); }
                s = syntax::nodes[s].next;
            }
            var t = syntax::nodes[model_prog].first;
            while t >= 0 {
                if syntax::nodes[t].kind == syntax::N_BC { dead_bc(t, c, warned); }
                t = syntax::nodes[t].next;
            }
        }
        c = syntax::nodes[c].next;
    }
}

// Calls `f(region, birth)` for each region of each `birth` and `death` of
// step `step` that has no `where`.
fn each_event<F>(step: i64) {
    var s = syntax::nodes[step].first;
    while s >= 0 {
        let k = syntax::nodes[s];
        if (k.kind == syntax::N_BIRTH || k.kind == syntax::N_DEATH) && syntax::child_of_kind(s, syntax::N_WHERE) < 0 {
            var r = k.first;
            while r >= 0 {
                F(r, k.kind == syntax::N_BIRTH);
                r = syntax::nodes[r].next;
            }
        }
        s = k.next;
    }
}

fn dead_bc(bc: i64, step: i64, warned: dictionary<i64, i64>&) {
    let region = syntax::nodes[bc].first;
    if warned.contains(bc) || !region_dead(region, 0) { return; }
    warned.insert(bc, 1);
    let sname = tok_text(syntax::nodes[step].tok);
    warn(first_tok(region), "load_no_cells",
         str("the region '", node_text(region), "' has no active cells in step '", sname, "', so the bc does nothing there"),
         "move the bc to the step of the birth of its region, or to a later step");
}

// True when region expression `n` has no active cells: a dead region, a
// `let` of one, the boundary of one, a union of dead regions, or an
// intersection or a difference with a dead region on the left.
recursive fn region_dead(n: i64, depth: i64) -> bool {
    if depth > 16 { return false; }
    if dead.contains(key(norm_text(n))) { return true; }
    let k = syntax::nodes[n];
    if k.kind == syntax::N_NAME {
        let s = visible.get(key(tok_text(k.tok)));
        if s && syms[s].kind == S_LET && syms[s].node >= 0 {
            return region_dead(syntax::nodes[syms[s].node].last, depth + 1);
        }
        return false;
    }
    if k.kind == syntax::N_CALL && call_name(k.first) == "boundary" {
        let a = syntax::nodes[k.first].next;
        return a >= 0 && region_dead(syntax::nodes[a].first, depth + 1);
    }
    if k.kind == syntax::N_BINARY {
        let op = tok_text(k.tok);
        let l = region_dead(k.first, depth + 1);
        if op == "|" { return l && region_dead(syntax::nodes[k.first].next, depth + 1); }
        if op == "&" { return l || region_dead(syntax::nodes[k.first].next, depth + 1); }
        if op == "-" { return l; }
    }
    return false;
}

// --- rules on the numbers of the model ---------------------------------------------------

// The method item of physics `p`, or -1.
fn method_of(p: i64) -> i64 {
    let name = tok_text(syms[physics[p].sym].tok);
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_METHOD && tok_text(syntax::nodes[c].tok) == name { return c; }
        c = syntax::nodes[c].next;
    }
    return -1;
}

// The displacement unknown (a vector of lengths) of physics `p`, or -1.
fn displacement_of(p: i64) -> i64 {
    let ph = physics[p];
    for j in ph.count {
        let t = syms[ph.first + j].ty;
        if is_num(t) && t.rank == 1 && t.dim == LENGTH { return ph.first + j; }
    }
    return -1;
}

// The cell order of unknown `u` in the method `m` (spec §13: default 2).
fn order_of(m: i64, u: i64) -> f64 {
    if m < 0 { return default_order(u); }
    let order = setting(syntax::nodes[m].first, "order");
    if order < 0 { return 2.0; }
    if syntax::nodes[order].kind == syntax::N_RECORD {
        var c = syntax::nodes[order].first;
        while c >= 0 {
            let a = syntax::nodes[c];
            if a.aux >= 0 && tok_text(a.aux) == tok_text(syms[u].tok) {
                let v, ok = const_value(a.first);
                return if ok { v } else { 2.0 };
            }
            c = a.next;
        }
        return 2.0;
    }
    let v, ok = const_value(order);
    return if ok { v } else { 2.0 };
}

// The constant value of property `name` of material symbol `mat`, and true.
fn prop_value(mat: i64, name: u8[:]) -> f64, bool {
    var c = syntax::nodes[syms[mat].node].first;
    while c >= 0 {
        let k = syntax::nodes[c];
        if k.kind == syntax::N_PROP && tok_text(k.tok) == name { return const_value(k.last); }
        c = k.next;
    }
    return 0.0, false;
}

// The mesh item of physics `p`: the mesh that its method names, or else the
// one mesh that a cg, dg or fv method with no `mesh` takes (spec §13), or -1.
fn mesh_of(p: i64) -> i64 {
    var mesh = -1;
    let m = method_of(p);
    let ms = if m >= 0 { setting(syntax::nodes[m].first, "mesh") } else { -1 };
    if ms >= 0 && syntax::nodes[ms].kind == syntax::N_NAME {
        let s = visible.get(key(tok_text(syntax::nodes[ms].tok)));
        if s && syms[s].kind == S_MESH { mesh = syms[s].node; }
    } else if ms < 0 && (m < 0 || syntax::nodes[m].aux >= 0 && contains_name(["cg", "dg", "fv"], tok_text(syntax::nodes[m].aux))) {
        var all = 0;
        for s in syms {
            if s.kind != S_MESH { continue; }
            all++;
            if all == 1 || (physics[p].region >= 0 && norm_text(syntax::nodes[s.node].first) == norm_text(physics[p].region)) {
                mesh = s.node;
            }
        }
    }
    return mesh;
}

// The smallest constant cell size of the mesh of physics `p`: the `size`
// and each `refine` and `crack_tip` of its mesh. False when it is unknown.
fn cell_size(p: i64) -> f64, bool {
    let mesh = mesh_of(p);
    guard mesh >= 0 else { return 0.0, false; };
    let b = syntax::child_of_kind(mesh, syntax::N_BLOCK);
    let size = setting(b, "size");
    guard size >= 0 else { return 0.0, false; };
    var h, ok = const_value(size);
    var c = syntax::nodes[b].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_REFINE || syntax::nodes[c].kind == syntax::N_CRACK_TIP {
            let v, vok = const_value(syntax::nodes[c].last);
            if vok && v < h { h = v; }
        }
        c = syntax::nodes[c].next;
    }
    return h, ok && h > 0.0;
}

// `crack_tip_form` (spec §9.1): a crack tip moves the mid-side nodes, so
// the displacement of each physics on the mesh needs order 2.
fn crack_tip_order() {
    for p, i in physics {
        let u = displacement_of(i);
        let mesh = mesh_of(i);
        guard u >= 0 && mesh >= 0 else { continue; };
        let order = order_of(method_of(i), u);
        var c = syntax::nodes[syntax::child_of_kind(mesh, syntax::N_BLOCK)].first;
        while c >= 0 {
            if syntax::nodes[c].kind == syntax::N_CRACK_TIP && order != 2.0 {
                error(syntax::nodes[c].tok, "crack_tip_form",
                      str("crack_tip moves the mid-side nodes, but ", tok_text(syms[u].tok), " has order ", order as i64),
                      str("write order: 2 in method ", tok_text(syms[p.sym].tok)));
            }
            c = syntax::nodes[c].next;
        }
    }
}

// A value rounded to 3 significant digits, for a message.
fn approx(v: f64) -> f64 {
    if v <= 0.0 { return v; }
    let scale = pow(10.0, floor(log10(v)) - 2.0);
    return round(v / scale) * scale;
}

// `volumetric_lock`: an order-1 displacement in `cg` or `dg`, with no
// pressure unknown, and a material with `nu` above 0.49.
fn volumetric_lock() {
    for p, i in physics {
        let u = displacement_of(i);
        let m = method_of(i);
        guard u >= 0 && m >= 0 && syntax::nodes[m].aux >= 0 else { continue; };
        let kind = tok_text(syntax::nodes[m].aux);
        guard (kind == "cg" || kind == "dg") && order_of(m, u) == 1.0 else { continue; };
        var pressure = false;
        for j in p.count {
            let t = syms[p.first + j].ty;
            if is_num(t) && t.rank == 0 && t.dim == FORCE - dim_pow(LENGTH, 2) { pressure = true; }
        }
        if pressure { continue; }
        for pm in phys_mats {
            if pm.phys != i { continue; }
            let nu, ok = prop_value(pm.mat, "nu");
            if ok && nu > 0.49 {
                let order = setting(syntax::nodes[m].first, "order");
                warn(first_tok(order), "volumetric_lock",
                     str("the displacement '", tok_text(syms[u].tok), "' has order 1 and material ", tok_text(syms[pm.mat].tok),
                         " has nu = ", nu, ", with no pressure unknown"),
                     "order: 2, or a mixed form with a pressure unknown");
                break;
            }
        }
    }
}

// `dt_unstable`: an explicit scheme with a `dt` above the estimated stable
// limit of a solid, h/(q*c), for the smallest cell size h, the cell order q
// and the largest P-wave speed c = sqrt(E*(1 - nu)/(rho*(1 + nu)*(1 - 2*nu))).
fn dt_unstable() {
    var c = syntax::nodes[model_prog].first;
    while c >= 0 {
        if syntax::nodes[c].kind == syntax::N_STEP { dt_unstable_step(c); }
        c = syntax::nodes[c].next;
    }
}

fn dt_unstable_step(step: i64) {
    let pr = syntax::nodes[step].first;
    guard tok_text(syntax::nodes[pr].tok) == "transient" else { return; };
    var dt = -1;
    var a = syntax::nodes[pr].first;
    var explicit = false;
    while a >= 0 {
        let k = syntax::nodes[a];
        let name = if k.aux >= 0 { tok_text(k.aux) } else { "" };
        if name == "dt" { dt = k.first; }
        if name == "scheme" && syntax::nodes[k.first].kind == syntax::N_NAME
           && contains_name(EXPLICIT_SCHEMES, tok_text(syntax::nodes[k.first].tok)) { explicit = true; }
        a = k.next;
    }
    guard explicit && dt >= 0 else { return; };
    let dtv, dok = const_value(dt);
    guard dok else { return; };
    for p, i in physics {
        let u = displacement_of(i);
        guard u >= 0 else { continue; };
        let h, hok = cell_size(i);
        guard hok else { continue; };
        var speed = 0.0;
        for pm in phys_mats {
            if pm.phys != i { continue; }
            let e, eok = prop_value(pm.mat, "E");
            let rho, rok = prop_value(pm.mat, "rho");
            var nu, nok = prop_value(pm.mat, "nu");
            if !nok || nu < 0.0 || nu >= 0.5 { nu = 0.0; }
            if eok && rok && rho > 0.0 && e > 0.0 {
                speed = max(speed, sqrt(e * (1.0 - nu) / (rho * (1.0 + nu) * (1.0 - 2.0 * nu))));
            }
        }
        guard speed > 0.0 else { continue; };
        let limit = h / (order_of(method_of(i), u) * speed);
        if dtv > limit {
            warn(first_tok(dt), "dt_unstable",
                 str("dt = ", approx(dtv), " s is above the estimated stable limit ", approx(limit), " s of physics '",
                     tok_text(syms[p.sym].tok), "' (cell size ", approx(h), " m, wave speed ", approx(speed), " m/s)"),
                 str("dt: ", approx(limit), " s or less"));
            return;
        }
    }
}

// `stabilization_missing`: `stabilize: none` in a method whose cell
// Reynolds number rho*U*h/mu or cell Peclet number rho*c_p*U*h/k is above
// 1. U is the velocity of the `scales` line, and h the smallest cell size.
fn stabilization_missing() {
    let speed = scale_velocity();
    guard speed > 0.0 else { return; };
    for p, i in physics {
        let m = method_of(i);
        guard m >= 0 else { continue; };
        let stab = setting(syntax::nodes[m].first, "stabilize");
        guard stab >= 0 && node_text(stab) == "none" else { continue; };
        let h, hok = cell_size(i);
        guard hok else { continue; };
        var cell = 0.0;
        for pm in phys_mats {
            if pm.phys != i { continue; }
            let rho, rok = prop_value(pm.mat, "rho");
            let mu, mok = prop_value(pm.mat, "mu");
            let k, kok = prop_value(pm.mat, "k");
            var cp, cok = prop_value(pm.mat, "c_p");
            if !cok { cp, cok = prop_value(pm.mat, "cp"); }
            if rok && mok && mu > 0.0 { cell = max(cell, rho * speed * h / mu); }
            if rok && kok && cok && k > 0.0 { cell = max(cell, rho * cp * speed * h / k); }
        }
        if cell > 1.0 {
            warn(first_tok(stab), "stabilization_missing",
                 str("the cell Reynolds or Peclet number of physics '", tok_text(syms[p.sym].tok), "' is about ", approx(cell),
                     " (cell size ", approx(h), " m, velocity ", approx(speed), " m/s)"),
                 "stabilize: supg, or a finer mesh");
        }
    }
}

// The constant velocity of the `scales` line in m/s, or 0.
fn scale_velocity() -> f64 {
    let sc = syntax::child_of_kind(model_prog, syntax::N_SCALES);
    guard sc >= 0 else { return 0.0; };
    var c = syntax::nodes[sc].first;
    while c >= 0 {
        let t = expr_ty(c);
        if is_num(t) && t.rank == 0 && t.dim == LENGTH - TIME {
            let v, ok = const_value(c);
            if ok { return abs(v); }
        }
        c = syntax::nodes[c].next;
    }
    return 0.0;
}
