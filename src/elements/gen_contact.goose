// The part of the simplex generator for contact (spec §12.4) and the rigid
// face (spec §12.1), on src/elements/contact.goose.
//
// `contact c: master, slave { friction: mu, enforce: penalty, stiffness:
// kn, stick_stiffness: kt }` names two boundaries `boundary(A, ...)` and
// `boundary(B, ...)` of two bodies A and B, and the body of the physics is
// their union `A | B`. The mesher meshes each body alone, so the bodies
// share no node (elements::mesh_bodies), and a boundary of a body holds
// only its nodes. With no stick_stiffness, kt is kn.
//
// `bc face: u = rigid(p)` makes face a rigid face, and `bc p: moment(sigma*n)
// = M` puts the moment M on it. `release r` in a step removes the bcs on r
// from that step on.
//
// A static step with a contact or a rigid face is Newton's method on the
// linear law (elements::solve_contact), or on the law of the nonlinear
// solid; `static(ramp: n)` takes n increments.
import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import mesher from "../mesher";
import .gen;

let BC_RIGID = 6;   // `bc face: u = rigid(p)`: node is the boundary of p
let BC_MOMENT = 7;  // `bc p: moment(sigma*n) = M`: v holds M

let CONTACT_FIX = "write contact c: block_face, pad_face { friction: 0.3, enforce: penalty, stiffness: 1e7 MPa/mm }, with block_face = boundary(block, y == 0 mm) and the body block | pad";

struct Contact {
    name: u8[..64],
    master: i64,    // the boundaries in surfs
    slave: i64,
    kn: f64,
    kt: f64,
    mu: f64,
}

// `release r` in step `step`.
struct Release {
    step: i64,
    surf: i64,
}

var contact_items: i64[>..];
var contacts: Contact[>..];
var body_names: u8[..64][>..];  // the bodies of the contacts
var body_shapes: i64[>..];      // their shapes of the mesher
var releases: Release[>..];
var has_rigid = false;

fn contact_reset() {
    contact_items.clear();
    contacts.clear();
    body_names.clear();
    body_shapes.clear();
    releases.clear();
    has_rigid = false;
}

// Whether a static step solves with Newton's method for the ties.
fn has_ties() -> bool { return contact_items.len > 0 || has_rigid; }

fn find_body(word: u8[:]) -> i64 {
    for b, i in body_names { if b == word { return i; } }
    return -1;
}

// The let of the boundary `word`, or -1.
fn surface_let(word: u8[:]) -> i64 {
    for s in surface_lets { if name(s) == word { return s; } }
    return -1;
}

// The bodies of the contacts, from the first argument of the boundary of
// each region, before the boundaries: the body of the physics is their
// union. It also finds the bcs with rigid(p).
fn contact_bodies(prog: i64) {
    find_rigid(prog);
    guard contact_items.len > 0 else { return; }
    for c in contact_items {
        for i in 2 {
            let r = child(c, i);
            let l = if kind(r) == syntax::N_NAME { surface_let(name(r)) } else { -1 };
            let v = if l >= 0 { child(l, 1) } else { -1 };
            let b = if v >= 0 && is_call(v, "boundary") && syntax::child_count(v) == 3 { arg(v, 0) } else { -1 };
            if b < 0 || kind(b) != syntax::N_NAME || name(b) == body_name {
                unsupported(tok(r), "a region of a contact is boundary(A, ...) of one body A of the physics", CONTACT_FIX);
                continue;
            }
            if find_body(name(b)) >= 0 { continue; }
            let g = mesher::eval_value(b);
            guard ok() && g.n == 0 else { return; }
            var word: u8[..64] = [];
            set64(word, name(b));
            body_names.push(word);
            body_shapes.push(g.g);
        }
    }
    let n = find_geometry(body_name);
    var parts: i64[>..] = [];
    if n >= 0 { union_parts(n, parts); }
    var same = parts.len == body_names.len;
    for p in parts { if find_body(name(p)) < 0 { same = false; } }
    if !same || geo != GEO_SHAPE {
        unsupported(tok(contact_items[0]), str("the body ", body_name, " of the physics is not the union of the bodies of the contacts"),
                    "write let parts: Surface = block | pad, with physics solid on parts");
    }
}

// The names of a | b | c into parts.
recursive fn union_parts(n: i64, parts: i64[>..]&) {
    if kind(n) == syntax::N_BINARY && syntax::text(tok(n)) == "|" {
        union_parts(child(n, 0), parts);
        union_parts(child(n, 1), parts);
    } else {
        parts.push(n);
    }
}

// Sets has_rigid if a bc of the program or of a step has u = rigid(p).
fn find_rigid(prog: i64) {
    var c = first(prog);
    while c >= 0 {
        if kind(c) == syntax::N_BC && is_call(child(c, 2), "rigid") { has_rigid = true; }
        if kind(c) == syntax::N_STEP {
            var s = next(child(c, 0));
            while s >= 0 {
                if kind(s) == syntax::N_BC && is_call(child(s, 2), "rigid") { has_rigid = true; }
                s = next(s);
            }
        }
        c = next(c);
    }
}

// A contact item, after the boundaries.
fn contact_item(n: i64) {
    var c = Contact { name: [], master: -1, slave: -1, kn: 0.0, kt: -1.0, mu: 0.0 };
    set64(c.name, name(n));
    c.master = find_surf(name(child(n, 0)));
    c.slave = find_surf(name(child(n, 1)));
    var penalty = false;
    var s = first(child(n, 2));
    while s >= 0 {
        let w = name(s);
        let v = first(s);
        if w == "friction" {
            c.mu = scalar_value(v);
        } else if w == "enforce" {
            penalty = kind(v) == syntax::N_NAME && name(v) == "penalty";
            if !penalty { unsupported(tok(v), "this generator enforces contact by penalty", "write enforce: penalty"); }
        } else if w == "stiffness" {
            c.kn = scalar_value(v);
        } else if w == "stick_stiffness" {
            c.kt = scalar_value(v);
        } else {
            unsupported(tok(s), str("the contact setting '", w, "' is not in this generator"), CONTACT_FIX);
        }
        s = next(s);
    }
    if c.kt < 0.0 { c.kt = c.kn; }
    if !penalty || c.kn <= 0.0 || c.mu < 0.0 {
        unsupported(tok(n), "a contact needs enforce: penalty, stiffness > 0 and friction >= 0", CONTACT_FIX);
    }
    if c.master < 0 || c.slave < 0 || surfs[c.master].body < 0 || surfs[c.slave].body < 0 || surfs[c.master].body == surfs[c.slave].body {
        unsupported(tok(n), "the master and the slave are boundaries of two bodies", CONTACT_FIX);
    }
    contacts.push(c);
}

// `bc face: u = rigid(p)` and `bc p: moment(sigma*n) = M` into b. False if
// the bc is neither.
fn rigid_bc(target: i64, v: i64, b: Bc&) -> bool {
    guard physics == ELASTIC else { return false; }
    let fix = "write bc top: u = rigid(c) and bc c: moment(sigma*n) = M, with let c: Point";
    if source_of(target) == unknown_name && is_call(v, "rigid") {
        let p = if syntax::child_count(v) == 2 && kind(arg(v, 0)) == syntax::N_NAME { find_surf(name(arg(v, 0))) } else { -1 };
        if p < 0 || surfs[p].points.len != 3 { unsupported(tok(v), "rigid(p) needs a Point p that a let names", fix); }
        b.kind = BC_RIGID;
        b.node = p;
        return true;
    }
    guard is_call(target, "moment") && syntax::child_count(target) == 2 && source_of(arg(target, 0)) == traction_target() else { return false; }
    let m = value(v);
    if dim == 2 && m.n == 1 {
        b.v = [0.0, 0.0, m.v[0]];
    } else if dim == 3 && m.n == 3 {
        b.v = m.v;
    } else {
        unsupported(tok(v), "the moment is a scalar about z in 2d and a vector in 3d", fix);
    }
    b.kind = BC_MOMENT;
    return true;
}

// `release r` of the step that step() reads.
fn release(n: i64) {
    let s = region_surf(first(n));
    if s < 0 {
        unsupported(tok(n), "release needs a boundary that a let names", "write release top");
        return;
    }
    releases.push(Release { step: steps.len, surf: s });
}

// Removes the bcs that step i releases from set.
fn release_bcs(i: i64, set: Bc[>..]&) {
    for r in releases {
        if r.step != i { continue; }
        var kept: Bc[>..] = [];
        for b in set { if b.surf != r.surf { kept.push(b); } }
        set.clear();
        set.append(kept);
    }
}

// --- the program ---------------------------------------------------------------

// The mesh of the bodies of the contacts, each alone.
fn emit_bodies(out: u8[>..]&) {
    format(out, ", by the mesher, each body alone\n    if !elements::mesh_bodies(", dim, ", ", order, ", [");
    for g, i in body_shapes {
        if i > 0 { out.append(", "); }
        format(out, g);
    }
    format(out, "], ", body_shape, ", ");
    lit(out, mesh_size);
    out.append(") {\n        write_stderr(\"skein run: the mesher did not finish in its budget of points: give a larger size\\n\");\n");
    out.append("        exit(1);\n    }\n");
}

// The pairs of each contact.
fn emit_contacts(out: u8[>..]&) {
    for c in contacts {
        format(out, "    // contact ", c.name, ": the master ", surfs[c.master].name, ", the slave ", surfs[c.slave].name, "\n");
        emit_region(out, c.master);
        format(out, "    region(elements::contact_master(), \"", surfs[c.master].name, "\");\n");
        emit_region(out, c.slave);
        out.append("    cells(elements::contact_slave(");
        lit(out, c.kn);
        out.append(", ");
        lit(out, c.kt);
        out.append(", ");
        lit(out, c.mu);
        format(out, ") > 0, \"contact ", c.name, "\");\n");
    }
}

// The calls of a BC_RIGID or a BC_MOMENT after its region.
fn emit_rigid(out: u8[>..]&, b: Bc) {
    if b.kind == BC_RIGID {
        format(out, "    region(elements::rigid_face(), \"", surfs[b.surf].name, "\");\n");
        return;
    }
    out.append("    elements::rigid_moment(");
    lits(out, b.v);
    out.append(");\n");
}
