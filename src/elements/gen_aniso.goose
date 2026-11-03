// The written-out linear elastic law of a physics: `let sigma` is any
// expression that is linear in `strain(u)`, for example the orthotropic
// stiffness written out with `outer`, a thermal strain tensor, or the
// engineering constants of a ply turned by its angle. The generator
// evaluates the expression at the 6 unit strains and at the zero strain, so
// that the law becomes one 6 x 6 stiffness and one initial stress for each
// material (src/elements/assemble.goose, set_aniso). The same code serves
// each physics that has an elastic flux, and each region with its own material.
import std;
import math;
namespace elements_gen;
import syntax from "../syntax";
import codegen.codegen;
import .gen;
import .gen_solid;

let ANISO_FIX = "write the flux as a sum of + - * / of numbers, material properties, lets, strain(u), I, outer, tr, transpose, sym, dev, dot, ddot and hooke, and the members .xx .xy of a tensor";

// A scalar (rank 0), a vector (rank 1, 3 components) or a tensor (rank 2,
// 3 x 3 in rows) in the unit system of the header.
struct Ten {
    rank: i64,
    v: f64[9],
}

var aniso_node = -1;        // the value of the flux let when the law is written out, or -1
var probe_eps: f64[9] = [0.0; 9];   // strain(u) of the evaluation now

fn aniso_reset() {
    aniso_node = -1;
    probe_eps = [0.0; 9];
}

fn t_scalar(x: f64) -> Ten { return Ten { rank: 0, v: [x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] }; }


fn t_eye() -> Ten { return Ten { rank: 2, v: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0] }; }

fn t_scale(a: Ten, f: f64) -> Ten {
    var r = a;
    for i in 9 { r.v[i] = a.v[i] * f; }
    return r;
}

fn t_trans(a: Ten) -> Ten {
    var r = a;
    for i in 3 {
        for j in 3 { r.v[3 * i + j] = a.v[3 * j + i]; }
    }
    return r;
}

fn t_trace(a: Ten) -> f64 { return a.v[0] + a.v[4] + a.v[8]; }

// The physics let `word`, or -1.
fn aniso_let(word: u8[:]) -> i64 {
    for l in solid_lets { if name(l) == word { return l; } }
    return -1;
}

fn member_index(w: u8[:], rank: i64) -> i64 {
    if rank == 1 && w == "x" { return 0; }
    if rank == 1 && w == "y" { return 1; }
    if rank == 1 && w == "z" { return 2; }
    guard rank == 2 && w.len == 2 else { return -1; }
    var idx = 0;
    for k in 2 {
        // x, y and z are the bytes 120, 121 and 122.
        let a = w[k] as! i64 - 120;
        if a < 0 || a > 2 { return -1; }
        idx = 3 * idx + a;
    }
    return idx;
}

recursive fn tval(n: i64, depth: i64) -> Ten {
    let k = kind(n);
    if depth > 60 {
        unsupported(tok(n), "the law nests lets more than 60 deep", ANISO_FIX);
        return t_scalar(0.0);
    }
    if k == syntax::N_NAME {
        let w = name(n);
        if w == "I" { return t_eye(); }
        let l = aniso_let(w);
        if l >= 0 { return tval(child(l, 1), depth + 1); }
        let p = property(w);
        let g = codegen::find_value(w);
        if g >= 0 { return tval(g, depth + 1); }
        if p >= 0 { return tval(p, depth + 1); }
        return t_scalar(scalar_value(n));
    }
    if k == syntax::N_VECTOR {
        let x = value(n);
        var r = Ten { rank: 1, v: [0.0; 9] };
        for j in 3 { r.v[j] = x.v[j]; }
        return r;
    }
    if k == syntax::N_NUMBER || k == syntax::N_PERCENT { return t_scalar(scalar_value(n)); }
    if k == syntax::N_UNARY && syntax::text(tok(n)) == "-" { return t_scale(tval(first(n), depth + 1), -1.0); }
    if k == syntax::N_MEMBER {
        let base = first(n);
        if kind(base) == syntax::N_NAME && position(codegen::table_names, name(base)) >= 0 { return t_scalar(scalar_value(n)); }
        let b = tval(base, depth + 1);
        let i = member_index(name(n), b.rank);
        if i < 0 {
            unsupported(tok(n), str("the member '", name(n), "' of a value of rank ", b.rank), ANISO_FIX);
            return t_scalar(0.0);
        }
        return t_scalar(b.v[i]);
    }
    if k == syntax::N_BINARY { return tbinary(n, depth); }
    if k == syntax::N_CALL { return tcall(n, depth); }
    unsupported(tok(n), "this expression is not in the written-out elastic law of this generator", ANISO_FIX);
    return t_scalar(0.0);
}

recursive fn tbinary(n: i64, depth: i64) -> Ten {
    let op = syntax::text(tok(n));
    let a = tval(child(n, 0), depth + 1);
    let b = tval(child(n, 1), depth + 1);
    if op == "+" || op == "-" {
        if a.rank != b.rank {
            unsupported(tok(n), "a sum of values of two ranks", ANISO_FIX);
            return a;
        }
        var r = a;
        for i in 9 { r.v[i] = if op == "+" { a.v[i] + b.v[i] } else { a.v[i] - b.v[i] }; }
        return r;
    }
    if op == "*" {
        if a.rank == 0 { return t_scale(b, a.v[0]); }
        if b.rank == 0 { return t_scale(a, b.v[0]); }
        if a.rank == 2 && b.rank == 2 {
            var r = Ten { rank: 2, v: [0.0; 9] };
            for i in 3 {
                for j in 3 {
                    for m in 3 { r.v[3 * i + j] += a.v[3 * i + m] * b.v[3 * m + j]; }
                }
            }
            return r;
        }
        if a.rank == 2 && b.rank == 1 {
            var r = Ten { rank: 1, v: [0.0; 9] };
            for i in 3 {
                for m in 3 { r.v[i] += a.v[3 * i + m] * b.v[m]; }
            }
            return r;
        }
        unsupported(tok(n), "a product of these two ranks", ANISO_FIX);
        return a;
    }
    if op == "/" {
        if b.rank != 0 || b.v[0] == 0.0 {
            unsupported(tok(n), "a division by a tensor or by 0", ANISO_FIX);
            return a;
        }
        return t_scale(a, 1.0 / b.v[0]);
    }
    if op == "^" {
        guard a.rank == 0 && b.rank == 0 else {
            unsupported(tok(n), "a power of a tensor", ANISO_FIX);
            return a;
        };
        let e = b.v[0];
        if e == floor(e) && abs(e) <= 16.0 {
            var r = 1.0;
            for i in abs(e) as! i64 { r *= a.v[0]; }
            return t_scalar(if e < 0.0 { 1.0 / r } else { r });
        }
        guard a.v[0] > 0.0 else {
            unsupported(tok(n), "a power of a value that is not positive", ANISO_FIX);
            return a;
        };
        return t_scalar(exp(e * log(a.v[0])));
    }
    unsupported(tok(n), str("the operator '", op, "' in the elastic law"), ANISO_FIX);
    return a;
}

recursive fn tcall(n: i64, depth: i64) -> Ten {
    let count = syntax::child_count(n) - 1;
    let f = if kind(child(n, 0)) == syntax::N_NAME { name(child(n, 0)) } else { "" };
    if f == "strain" && count == 1 && source_of(arg(n, 0)) == unknown_name { return Ten { rank: 2, v: probe_eps }; }
    if f == "outer" && count == 2 {
        let a = tval(arg(n, 0), depth + 1);
        let b = tval(arg(n, 1), depth + 1);
        if a.rank == 1 && b.rank == 1 {
            var r = Ten { rank: 2, v: [0.0; 9] };
            for i in 3 {
                for j in 3 { r.v[3 * i + j] = a.v[i] * b.v[j]; }
            }
            return r;
        }
    }
    if (f == "tr" || f == "transpose" || f == "sym" || f == "dev") && count == 1 {
        let a = tval(arg(n, 0), depth + 1);
        if a.rank == 2 {
            if f == "tr" { return t_scalar(t_trace(a)); }
            if f == "transpose" { return t_trans(a); }
            let t = t_trans(a);
            var r = a;
            for i in 9 { r.v[i] = if f == "sym" { 0.5 * (a.v[i] + t.v[i]) } else { a.v[i] }; }
            if f == "dev" {
                let m = t_trace(a) / 3.0;
                r.v[0] -= m;
                r.v[4] -= m;
                r.v[8] -= m;
            }
            return r;
        }
    }
    if (f == "dot" || f == "ddot") && count == 2 {
        let a = tval(arg(n, 0), depth + 1);
        let b = tval(arg(n, 1), depth + 1);
        if f == "dot" && a.rank == 1 && b.rank == 1 { return t_scalar(a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2]); }
        if f == "ddot" && a.rank == 2 && b.rank == 2 {
            var s = 0.0;
            for i in 9 { s += a.v[i] * b.v[i]; }
            return t_scalar(s);
        }
        if f == "dot" && a.rank == 2 && b.rank == 1 {
            var r = Ten { rank: 1, v: [0.0; 9] };
            for i in 3 {
                for m in 3 { r.v[i] += a.v[3 * i + m] * b.v[m]; }
            }
            return r;
        }
    }
    if f == "hooke" && count == 3 {
        let e = tval(arg(n, 0), depth + 1);
        let nu = tval(arg(n, 1), depth + 1);
        let t = tval(arg(n, 2), depth + 1);
        if e.rank == 0 && nu.rank == 0 && t.rank == 2 {
            let mu = e.v[0] / (2.0 * (1.0 + nu.v[0]));
            let lam = e.v[0] * nu.v[0] / ((1.0 + nu.v[0]) * (1.0 - 2.0 * nu.v[0]));
            var r = t_scale(t, 2.0 * mu);
            r.v[0] += lam * t_trace(t);
            r.v[4] += lam * t_trace(t);
            r.v[8] += lam * t_trace(t);
            return r;
        }
    }
    unsupported(tok(n), str("the call '", f, "' in the elastic law of this generator, or its arguments"), ANISO_FIX);
    return t_scalar(0.0);
}

// The stress of the strain probe_eps by the written-out law, in Voigt order.
fn aniso_stress_now() -> f64[6] {
    let s = tval(aniso_node, 0);
    var r: f64[6] = [0.0; 6];
    if s.rank != 2 {
        unsupported(tok(aniso_node), "the flux is not a tensor of rank 2", ANISO_FIX);
        return r;
    }
    let vi: i64[6] = [0, 1, 2, 0, 1, 0];
    let vj: i64[6] = [0, 1, 2, 1, 2, 2];
    for v in 6 { r[v] = s.v[3 * vi[v] + vj[v]]; }
    return r;
}

// The stiffness c (the stress for each unit strain xx, yy, zz, xy, yz, xz, with
// the tensor shear eps_xy = eps_yx = 1, in rows) and the stress s0 of the zero
// strain of the law of the material now. A law that is not linear in strain(u)
// is not in this generator.
fn aniso_law(c: f64[36]&, s0: f64[6]&) {
    let vi: i64[6] = [0, 1, 2, 0, 1, 0];
    let vj: i64[6] = [0, 1, 2, 1, 2, 2];
    probe_eps = [0.0; 9];
    s0 = aniso_stress_now();
    var big = 0.0;
    for w in 6 {
        probe_eps = [0.0; 9];
        probe_eps[3 * vi[w] + vj[w]] = 1.0;
        probe_eps[3 * vj[w] + vi[w]] = 1.0;
        let s = aniso_stress_now();
        for v in 6 {
            c[6 * v + w] = s[v] - s0[v];
            big = max(big, abs(c[6 * v + w]));
        }
    }
    // A second strain checks that the law is linear: the stress of the sum of the strains.
    probe_eps = [1.0, 0.25, 0.0, 0.25, 0.5, 0.0, 0.0, 0.0, -0.5];
    let s = aniso_stress_now();
    let e: f64[6] = [1.0, 0.5, -0.5, 0.25, 0.0, 0.0];
    for v in 6 {
        var want = s0[v];
        for w in 6 { want += c[6 * v + w] * e[w]; }
        if abs(want - s[v]) > 1.0e-9 * (big + abs(want)) {
            unsupported(tok(aniso_node), "the flux is not linear in strain(u)", ANISO_FIX);
            return;
        }
    }
    probe_eps = [0.0; 9];
    if big <= 0.0 { unsupported(tok(aniso_node), "the flux does not depend on strain(u)", ANISO_FIX); }
    if dim == 2 && plane == PLANE_STRESS { condense_plane_stress(c, s0); }
}

// Plane stress: the stress zz is 0, so the strain zz leaves the law, which the
// 2d cells use in the components xx, yy and xy.
fn condense_plane_stress(c: f64[36]&, s0: f64[6]&) {
    let czz = c[14];
    if czz <= 0.0 || c[6 * 2 + 4] != 0.0 || c[6 * 2 + 5] != 0.0 {
        unsupported(tok(aniso_node), "a plane stress law needs C_zz > 0 and no coupling of the strain zz with yz and xz", ANISO_FIX);
        return;
    }
    let idx: i64[3] = [0, 1, 3];
    var d: f64[36] = c;
    for a in 3 {
        for b in 3 { d[6 * idx[a] + idx[b]] = c[6 * idx[a] + idx[b]] - c[6 * idx[a] + 2] * c[6 * 2 + idx[b]] / czz; }
    }
    var t: f64[6] = s0;
    for a in 3 { t[idx[a]] = s0[idx[a]] - c[6 * idx[a] + 2] * s0[2] / czz; }
    t[2] = 0.0;
    for a in 6 { d[6 * 2 + a] = 0.0; d[6 * a + 2] = 0.0; }
    c = d;
    s0 = t;
}

// The parts of this generator that read lame and shear, and so not the written-out law.
fn aniso_limits() {
    if contacts.len > 0 { unsupported(tok(aniso_node), "a contact of a written-out elastic law is not in this generator", "use hooke(E, nu, strain(u))"); }
    if rkpm_support > 0.0 { unsupported(tok(aniso_node), "the method rkpm of a written-out elastic law is not in this generator", "use method cg"); }
    if cut_source >= 0 { unsupported(tok(aniso_node), "a submodel source of a written-out elastic law is not in this generator", "use hooke(E, nu, strain(u))"); }
    if solid_law != 0 { unsupported(tok(aniso_node), "a nonlinear solid with a written-out elastic law is not in this generator", "use hooke(E, nu, strain(u) - eps_p)"); }
}
