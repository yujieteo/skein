// The 3D tetrahedron mesher (spec §9.2): isosurface stuffing after Labelle
// and Shewchuk (2007), on a uniform body-centred cubic (BCC) lattice.
//
// 1. The lattice covers the bounding box of the shape, with one cube more
//    on each side. The spacing along each axis divides the side of the box
//    into a whole number of cubes no larger than `size`, so the planes of
//    an axis-aligned box fall on lattice points. The lattice points are
//    the cube corners and the cube centres. Its tetrahedra each have 2
//    centres of neighbour cubes and an edge of their common face.
// 2. φ at each lattice point gives its sign: inside (-), on the surface
//    (0) or outside (+). Each lattice edge from - to + has a cut point on
//    the surface, found by bisection of φ.
// 3. Warping: a cut point nearer to an end of its edge than a fraction
//    alpha of the edge moves that lattice point onto the surface, which
//    makes the point 0. alpha is 0.24999 on the long edges (the axis
//    edges) and 0.41189 on the short edges (corner to centre), the values
//    of the paper.
// 4. Each lattice tetrahedron gives the part of it inside the surface: a
//    tetrahedron, a pyramid or a prism. A pyramid or a prism is cut into
//    tetrahedra along the diagonal of each quadrilateral face from its
//    point of smallest index (Dompierre et al. 1999), so two cells that
//    share the face cut it the same way, and the mesh is conforming.
//
// The paper's own cuts of the pyramids and the prisms give dihedral angles
// from 10.7 to 164.8 degrees. These cuts do not have that bound, so
// mesh_quality measures each mesh. The sharp edges of the shape are cut
// off (chamfered) by about one lattice cell.
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .quality;

let ALPHA_LONG = 0.24999;
let ALPHA_SHORT = 0.41189;

var lattice_cubes: i64[3] = [0, 0, 0];
var lattice_origin: f64[3] = [0.0, 0.0, 0.0];   // corner (0, 0, 0)
var lattice_step: f64[3] = [0.0, 0.0, 0.0];

var lx: f64[>..] = [];          // x, y and z of each lattice point
var lphi: f64[>..] = [];        // φ, or 0 after warping
var lnode: i64[>..] = [];       // the mesh node of each lattice point, or -1
var cuts = dictionary<i64, i64> {};

fn corner_count() -> i64 {
    return (lattice_cubes[0] + 1) * (lattice_cubes[1] + 1) * (lattice_cubes[2] + 1);
}

fn corner_id(i: i64, j: i64, k: i64) -> i64 {
    return i + (lattice_cubes[0] + 1) * (j + (lattice_cubes[1] + 1) * k);
}

fn centre_id(i: i64, j: i64, k: i64) -> i64 {
    return corner_count() + i + lattice_cubes[0] * (j + lattice_cubes[1] * k);
}

fn lattice_point_count() -> i64 {
    return corner_count() + lattice_cubes[0] * lattice_cubes[1] * lattice_cubes[2];
}

// The lattice point id, and its integer position in half steps.
fn half_position(p: i64) -> i64, i64, i64 {
    let n = corner_count();
    if p < n {
        let a = lattice_cubes[0] + 1;
        let b = lattice_cubes[1] + 1;
        return 2 * (p % a), 2 * (p / a % b), 2 * (p / (a * b));
    }
    let q = p - n;
    return 2 * (q % lattice_cubes[0]) + 1, 2 * (q / lattice_cubes[0] % lattice_cubes[1]) + 1,
           2 * (q / (lattice_cubes[0] * lattice_cubes[1])) + 1;
}

// The lattice point at half-step position (i, j, k), or -1.
fn point_at_half(i: i64, j: i64, k: i64) -> i64 {
    let ci = i % 2;
    guard ci == j % 2 && ci == k % 2 else { return -1; }
    guard i >= 0 && j >= 0 && k >= 0 else { return -1; }
    if ci == 0 {
        guard i / 2 <= lattice_cubes[0] && j / 2 <= lattice_cubes[1] && k / 2 <= lattice_cubes[2] else { return -1; }
        return corner_id(i / 2, j / 2, k / 2);
    }
    guard i / 2 < lattice_cubes[0] && j / 2 < lattice_cubes[1] && k / 2 < lattice_cubes[2] else { return -1; }
    return centre_id(i / 2, j / 2, k / 2);
}

// The 14 lattice neighbours, in half steps: 6 along the axes (2 half
// steps, the long edges) and 8 along the diagonals (the short edges).
let NEIGHBOURS: i64[3][14] = [
    [2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2],
    [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
    [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1],
];

fn sign_of(v: f64) -> i64 {
    if v < 0.0 { return -1; }
    if v > 0.0 { return 1; }
    return 0;
}

// The point on the surface between lattice points a (inside) and b
// (outside), as the fraction of the way from a to b.
fn cut_fraction(g: i64, a: i64, b: i64) -> f64 {
    var lo = 0.0;
    var hi = 1.0;
    for k in 60 {
        let t = 0.5 * (lo + hi);
        let x = lx[3 * a] + t * (lx[3 * b] - lx[3 * a]);
        let y = lx[3 * a + 1] + t * (lx[3 * b + 1] - lx[3 * a + 1]);
        let z = lx[3 * a + 2] + t * (lx[3 * b + 2] - lx[3 * a + 2]);
        if phi(g, x, y, z) < 0.0 { lo = t; } else { hi = t; }
    }
    return 0.5 * (lo + hi);
}

// The mesh node of lattice point p, made when first used.
fn lattice_node(p: i64) -> i64 {
    if lnode[p] < 0 { lnode[p] = add_point(lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]); }
    return lnode[p];
}

// The mesh node at the cut point of the edge from inside point a to
// outside point b.
fn cut_node(g: i64, a: i64, b: i64) -> i64 {
    let n = lattice_point_count();
    let key = min(a, b) * n + max(a, b);
    let have = cuts.get(key);
    if have { return have; }
    let t = cut_fraction(g, a, b);
    let node = add_point(lx[3 * a] + t * (lx[3 * b] - lx[3 * a]),
                         lx[3 * a + 1] + t * (lx[3 * b + 1] - lx[3 * a + 1]),
                         lx[3 * a + 2] + t * (lx[3 * b + 2] - lx[3 * a + 2]));
    cuts.insert(key, node);
    return node;
}

fn add_tet(a: i64, b: i64, c: i64, d: i64) {
    let pa = vertex(a);
    let u = sub3(vertex(b), pa);
    let v = sub3(vertex(c), pa);
    let w = sub3(vertex(d), pa);
    if dot3(u, cross3(v, w)) < 0.0 {
        cells.push(a); cells.push(c); cells.push(b); cells.push(d);
    } else {
        cells.push(a); cells.push(b); cells.push(c); cells.push(d);
    }
}

// A pyramid with apex z over the quadrilateral q0 q1 q2 q3, in cyclic
// sequence.
fn add_pyramid(z: i64, q0: i64, q1: i64, q2: i64, q3: i64) {
    let m = min(min(q0, q1), min(q2, q3));
    if m == q0 || m == q2 {
        add_tet(z, q0, q1, q2);
        add_tet(z, q0, q2, q3);
    } else {
        add_tet(z, q1, q2, q3);
        add_tet(z, q1, q3, q0);
    }
}

// A prism with triangles a0 a1 a2 and b0 b1 b2, and side edges a_i b_i.
fn add_prism(a0: i64, a1: i64, a2: i64, b0: i64, b1: i64, b2: i64) {
    var a: i64[3] = [a0, a1, a2];
    var b: i64[3] = [b0, b1, b2];
    // Put the point of smallest index at a[m].
    let low = min(min(min(a0, a1), min(a2, b0)), min(b1, b2));
    if low == b0 || low == b1 || low == b2 { swap(a, b); }
    var m = 0;
    for i in 3 {
        if a[i] == low { m = i; }
    }
    let r: i64[3] = [a[m], a[(m + 1) % 3], a[(m + 2) % 3]];
    let s: i64[3] = [b[m], b[(m + 1) % 3], b[(m + 2) % 3]];
    add_tet(r[0], s[0], s[1], s[2]);
    // The quadrilateral r1 r2 s2 s1: its diagonal from its smallest point.
    if min(r[1], s[2]) < min(r[2], s[1]) {
        add_tet(r[0], r[1], r[2], s[2]);
        add_tet(r[0], r[1], s[2], s[1]);
    } else {
        add_tet(r[0], r[1], r[2], s[1]);
        add_tet(r[0], r[2], s[2], s[1]);
    }
}

// The part inside the surface of the lattice tetrahedron p0 p1 p2 p3.
fn stuff(g: i64, p: i64[4]) {
    var minus: i64[..4] = [];
    var zero: i64[..4] = [];
    var plus: i64[..4] = [];
    for v in p {
        let s = sign_of(lphi[v]);
        if s < 0 { minus.push(v); } else if s == 0 { zero.push(v); } else { plus.push(v); }
    }
    if minus.len == 0 {
        // Four points on the surface: keep the cell if its centre is inside.
        guard zero.len == 4 else { return; }
        var c: f64[3] = [0.0, 0.0, 0.0];
        for v in p {
            for i in 3 { c[i] += 0.25 * lx[3 * v + i]; }
        }
        guard phi(g, c[0], c[1], c[2]) < 0.0 else { return; }
    }
    if plus.len == 0 {
        add_tet(lattice_node(p[0]), lattice_node(p[1]), lattice_node(p[2]), lattice_node(p[3]));
        return;
    }
    if minus.len == 3 {
        add_prism(lattice_node(minus[0]), lattice_node(minus[1]), lattice_node(minus[2]),
                  cut_node(g, minus[0], plus[0]), cut_node(g, minus[1], plus[0]), cut_node(g, minus[2], plus[0]));
    } else if minus.len == 2 && plus.len == 2 {
        add_prism(lattice_node(minus[0]), cut_node(g, minus[0], plus[0]), cut_node(g, minus[0], plus[1]),
                  lattice_node(minus[1]), cut_node(g, minus[1], plus[0]), cut_node(g, minus[1], plus[1]));
    } else if minus.len == 2 {
        // One point on the surface and one outside.
        add_pyramid(lattice_node(zero[0]), lattice_node(minus[0]), lattice_node(minus[1]),
                    cut_node(g, minus[1], plus[0]), cut_node(g, minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 1 {
        add_tet(lattice_node(minus[0]), lattice_node(zero[0]), lattice_node(zero[1]),
                cut_node(g, minus[0], plus[0]));
    } else if minus.len == 1 && plus.len == 2 {
        add_tet(lattice_node(minus[0]), lattice_node(zero[0]),
                cut_node(g, minus[0], plus[0]), cut_node(g, minus[0], plus[1]));
    } else if minus.len == 1 {
        add_tet(lattice_node(minus[0]), cut_node(g, minus[0], plus[0]),
                cut_node(g, minus[0], plus[1]), cut_node(g, minus[0], plus[2]));
    }
}

// Meshes the solid g with tetrahedra of edge length about `size`.
fn mesh_3d(g: i64, size: f64) -> bool {
    clear_mesh(3);
    let b = bounds(g);
    for a in 3 {
        let side = b[3 + a] - b[a];
        guard side > 0.0 else { return false; }
        let n = max(ceil(side / size - 1.0e-9) as i64, 1);
        lattice_step[a] = side / n;
        lattice_cubes[a] = n + 2;
        lattice_origin[a] = b[a] - lattice_step[a];
    }
    let points = lattice_point_count();
    lx.clear();
    lphi.clear();
    lnode.clear();
    cuts.clear();
    let h = min(min(lattice_step[0], lattice_step[1]), lattice_step[2]);
    for p in points {
        let i, j, k = half_position(p);
        let x = lattice_origin[0] + 0.5 * i * lattice_step[0];
        let y = lattice_origin[1] + 0.5 * j * lattice_step[1];
        let z = lattice_origin[2] + 0.5 * k * lattice_step[2];
        lx.push(x);
        lx.push(y);
        lx.push(z);
        let f = phi(g, x, y, z);
        lphi.push(if abs(f) <= 1.0e-9 * h { 0.0 } else { f });
        lnode.push(-1);
    }
    // Warping: move each lattice point that a cut point violates onto the
    // nearest such cut point.
    for p in points {
        guard lphi[p] != 0.0 else { continue; }
        let i, j, k = half_position(p);
        var best = 2.0;
        var best_x: f64[3] = [0.0, 0.0, 0.0];
        for d, e in NEIGHBOURS {
            let q = point_at_half(i + d[0], j + d[1], k + d[2]);
            guard q >= 0 && lphi[q] != 0.0 && (lphi[q] < 0.0) != (lphi[p] < 0.0) else { continue; }
            let inside = if lphi[p] < 0.0 { p } else { q };
            let outside = if lphi[p] < 0.0 { q } else { p };
            let t = cut_fraction(g, inside, outside);
            let from_p = if inside == p { t } else { 1.0 - t };
            let alpha = if e < 6 { ALPHA_LONG } else { ALPHA_SHORT };
            if from_p < alpha && from_p < best {
                best = from_p;
                for c in 3 { best_x[c] = lx[3 * p + c] + from_p * (lx[3 * q + c] - lx[3 * p + c]); }
            }
        }
        if best < 2.0 {
            for c in 3 { lx[3 * p + c] = best_x[c]; }
            lphi[p] = 0.0;
        }
    }
    // The lattice tetrahedra: for each pair of cubes with a common face,
    // the 2 centres and each edge of the face.
    let n = lattice_cubes;
    for k in n[2] {
        for j in n[1] {
            for i in n[0] {
                let c = centre_id(i, j, k);
                for axis in 3 {
                    var o: i64[3] = [i, j, k];
                    o[axis] += 1;
                    guard o[axis] < n[axis] else { continue; }
                    let c2 = centre_id(o[0], o[1], o[2]);
                    // The corners of the common face, in cyclic sequence.
                    let u = (axis + 1) % 3;
                    let v = (axis + 2) % 3;
                    var f: i64[4] = [0, 0, 0, 0];
                    let du: i64[4] = [0, 1, 1, 0];
                    let dv: i64[4] = [0, 0, 1, 1];
                    for q in 4 {
                        var at: i64[3] = o;
                        at[u] += du[q];
                        at[v] += dv[q];
                        f[q] = corner_id(at[0], at[1], at[2]);
                    }
                    for q in 4 { stuff(g, [c, c2, f[q], f[(q + 1) % 4]]); }
                }
            }
        }
    }
    return true;
}
