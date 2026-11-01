// The near null space of each physics for the AMG of amg.goose: the
// vectors that the matrix of the physics takes near 0. Each solve with
// CG and AMG (the elements package, the box runtime) builds its levels
// with these modes, so that a physics gives them in one place.
//
// - A scalar field (heat): the constant.
// - Elasticity: the rigid-body modes, the translations and then the
//   rotations about z, x and y (3 in 2D, 6 in 3D).
// - Elasticity of a thin 3D body: also 3 plate bending modes about the
//   thinnest axis of the bounding box. Reason: on interposer_warpage
//   (139449 dofs, layers of 10 um to 786 um under cells of 0.5 mm) the
//   rigid-body modes need 1425 CG iterations and the bending modes 698.
//
// A mode is 0 on the fixed dofs.
import std;
import math;
namespace solvers;

let ELASTIC_MODES = 0;
let SCALAR_MODES = 1;

// A body is thin when its thinnest side is less than THIN times its
// middle side.
let THIN = 0.1;

// The thinnest axis of the bounding box of `coords` (x, y and z of each
// node), if a body in dim dimensions is thin there, or -1. The bounding box
// in lo and hi.
fn thin_axis(coords: f64[:], dim: i64, lo: f64[3]&, hi: f64[3]&) -> i64 {
    guard coords.len >= 3 && coords.len % 3 == 0 else { abort("thin_axis: 3 coordinates for each node"); }
    for t in 3 {
        lo[t] = coords[t];
        hi[t] = coords[t];
    }
    for n in coords.len / 3 {
        for t in 3 {
            lo[t] = min(lo[t], coords[3 * n + t]);
            hi[t] = max(hi[t], coords[3 * n + t]);
        }
    }
    var a = 0;
    for t in 3 { if hi[t] - lo[t] < hi[a] - lo[a] { a = t; } }
    let middle = min(hi[(a + 1) % 3] - lo[(a + 1) % 3], hi[(a + 2) % 3] - lo[(a + 2) % 3]);
    return if dim == 3 && hi[a] - lo[a] < THIN * middle { a } else { -1 };
}

// The number of modes of near_null_space: the vectors of one value for
// each dof that it gives.
fn mode_count(kind: i64, coords: f64[:], dim: i64, ncomp: i64) -> i64 {
    guard kind == ELASTIC_MODES && ncomp == dim else { return ncomp; }
    var lo: f64[3] = [0.0; 3];
    var hi: f64[3] = [0.0; 3];
    let rigid = if dim == 3 { 6 } else { 3 };
    return if thin_axis(coords, dim, lo, hi) >= 0 { rigid + 3 } else { rigid };
}

// The modes of `kind` for the dofs of `fixed`, in blocks of ncomp dofs,
// one block for each node of `coords` (x, y and z of each node). It gives
// the modes, nb values for each dof by rows, and nb.
fn near_null_space(kind: i64, coords: f64[:], dim: i64, ncomp: i64, fixed: u8[:]) -> f64[>..], i64 {
    guard ncomp >= 1 && fixed.len % ncomp == 0 && coords.len == 3 * (fixed.len / ncomp) && (dim == 2 || dim == 3) else { abort("near_null_space: ncomp dofs and 3 coordinates for each node, in 2D or 3D"); }
    var modes: f64[>..] = [];
    guard kind == ELASTIC_MODES && ncomp == dim else {
        for d in fixed.len {
            for c in ncomp { modes.push(if d % ncomp == c && fixed[d] == 0 { 1.0 } else { 0.0 }); }
        }
        return modes, ncomp;
    }
    // The modes are about the centre of the bounding box, so that their
    // values are of the size of the body.
    var lo: f64[3] = [0.0; 3];
    var hi: f64[3] = [0.0; 3];
    let a = thin_axis(coords, dim, lo, hi);
    let bending = a >= 0;
    // The other two axes: the plane of a thin body.
    let b = (a + 1) % 3;
    let c = (a + 2) % 3;
    let rigid = if dim == 3 { 6 } else { 3 };
    let nb = if bending { rigid + 3 } else { rigid };
    for d in fixed.len {
        let n = d / ncomp;
        let k = d % ncomp;
        var p: f64[3] = [0.0; 3];
        for t in 3 { p[t] = coords[3 * n + t] - 0.5 * (lo[t] + hi[t]); }
        var m: f64[9] = [0.0; 9];
        m[k] = 1.0;
        if k == 0 { m[dim] = -p[1]; }
        if k == 1 { m[dim] = p[0]; }
        if dim == 3 && k == 1 { m[4] = -p[2]; }
        if dim == 3 && k == 2 { m[4] = p[1]; }
        if dim == 3 && k == 0 { m[5] = p[2]; }
        if dim == 3 && k == 2 { m[5] = -p[0]; }
        // The bending modes: the curvatures of the plate along b, along c,
        // and the twist.
        if bending && k == a {
            m[6] = p[b] * p[b] / 2.0;
            m[7] = p[c] * p[c] / 2.0;
            m[8] = p[b] * p[c];
        }
        if bending && k == b {
            m[6] = -p[b] * p[a];
            m[8] = -p[c] * p[a];
        }
        if bending && k == c {
            m[7] = -p[c] * p[a];
            m[8] = -p[b] * p[a];
        }
        for j in nb { modes.push(if fixed[d] != 0 { 0.0 } else { m[j] }); }
    }
    return modes, nb;
}
