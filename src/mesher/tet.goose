// The 3D tetrahedron mesher (spec §9.2): isosurface stuffing after Labelle
// and Shewchuk (2007), on a uniform body-centred cubic (BCC) lattice.
//
// 1. The lattice covers the bounding box of the shape, with one cube more
//    on each side. Each axis-aligned plane of a box of the CSG tree in the
//    bounds is a plane of cube corners, and between two neighbour planes
//    the lattice has the fewest equal steps that are not longer than
//    `size` (grid_lines of grid.goose). Thus an edge of a box, where 2 such
//    planes meet, is on edges of the mesh. The lattice points are the cube
//    corners and the cube centres. Its tetrahedra each have 2 centres of
//    neighbour cubes and an edge of their common face. refine_cells
//    (refine.goose) then bisects them to the refine targets.
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
// 5. improve (improve.goose) repairs the cells with γ < 0.35.
//
// The paper's own cuts of the pyramids and the prisms give dihedral angles
// from 10.7 to 164.8 degrees. These cuts do not have that bound, so
// mesh_quality measures each mesh. A curved or an oblique sharp edge of
// the shape is cut off (chamfered) by about one third of a lattice cell. So a prism along z with a
// cylinder or an extrude in it is not stuffed: sweep_tets (sweep.goose)
// sweeps the triangles of its section, and keeps its sharp edges, and an
// octant of a spherical shell is swept radially (shell_tets).
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .quality;
import .refine;
import .grid;
import .improve;
import .feature;

let ALPHA_LONG = 0.24999;
let ALPHA_SHORT = 0.41189;

var lattice_cubes: i64[3] = [0, 0, 0];
// The planes of the cube corners along each axis: those of axis a are
// lattice_lines[line_start[a]..line_start[a + 1]], lattice_cubes[a] + 1.
var lattice_lines: f64[>..] = [];
var line_start: i64[4] = [0, 0, 0, 0];

var lx: f64[>..] = [];          // x, y and z of each lattice point
var lphi: f64[>..] = [];        // φ, or 0 after warping
var lnode: i64[>..] = [];       // the mesh node of each lattice point, or -1
var on_surface: bool[>..] = [];  // the mesh node is on the surface
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

// The integer position of lattice point p in half cubes along each axis.
fn half_position(p: i64) -> i64[3] {
    let n = corner_count();
    if p < n {
        let a = lattice_cubes[0] + 1;
        let b = lattice_cubes[1] + 1;
        return [2 * (p % a), 2 * (p / a % b), 2 * (p / (a * b))];
    }
    let q = p - n;
    return [2 * (q % lattice_cubes[0]) + 1, 2 * (q / lattice_cubes[0] % lattice_cubes[1]) + 1,
            2 * (q / (lattice_cubes[0] * lattice_cubes[1])) + 1];
}

// Appends the lo and hi corners of each box that the unions, the
// intersects and the subtracts of shape g reach, 6 numbers for each box.
fn box_leaves(g: i64, out: f64[>..]&) {
    var todo: i64[>..] = [g];
    while todo.len > 0 {
        let s = shapes[todo[todo.len - 1]];
        todo.pop();
        if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT {
            todo.push(s.b);
            todo.push(s.a);
        } else if s.kind == BOX {
            for i in 6 { out.push(s.p[i]); }
        }
    }
}

fn sign_of(v: f64) -> i64 {
    if v < 0.0 { return -1; }
    if v > 0.0 { return 1; }
    return 0;
}

// The point on the surface between lattice points a (inside) and b
// (outside), as the fraction of the way from a to b.
fn cut_fraction(g: i64, a: i64, b: i64) -> f64 { return cut_fraction_steps(g, a, b, 60); }

fn cut_fraction_steps(g: i64, a: i64, b: i64, steps: i64) -> f64 {
    var lo = 0.0;
    var hi = 1.0;
    for k in steps {
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
    if lnode[p] < 0 {
        lnode[p] = add_point(lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        on_surface.push(lphi[p] == 0.0);
    }
    return lnode[p];
}

// The mesh node at the cut point of the edge from inside point a to
// outside point b.
fn cut_node(g: i64, a: i64, b: i64) -> i64 {
    let key = edge_key(a, b);
    let have = cuts.get(key);
    if have { return have; }
    let t = cut_fraction(g, a, b);
    let node = add_point(lx[3 * a] + t * (lx[3 * b] - lx[3 * a]),
                         lx[3 * a + 1] + t * (lx[3 * b + 1] - lx[3 * a + 1]),
                         lx[3 * a + 2] + t * (lx[3 * b + 2] - lx[3 * a + 2]));
    cuts.insert(key, node);
    on_surface.push(true);
    return node;
}

fn add_tet(a: i64, b: i64, c: i64, d: i64) {
    let pa = vertex(a);
    let u = sub3(vertex(b), pa);
    let v = sub3(vertex(c), pa);
    let w = sub3(vertex(d), pa);
    if dot3(u, cross3(v, w)) < 0.0 {
        add_simplex([a, c, b, d]);
    } else {
        add_simplex([a, b, c, d]);
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

// The background tetrahedra, 4 lattice points each: the lattice cells,
// after refine_cells.
var btets: i64[>..] = [];

// The edges of the background tetrahedra (ends ea and eb), whether each is
// the longest edge of a tetrahedron (a long edge of the paper), and the
// edges at each point: edges_at[edge_start[p]..edge_start[p + 1]].
var edge_index = dictionary<i64, i64> {};
var ea: i64[>..] = [];
var eb: i64[>..] = [];
var elong: bool[>..] = [];
var edge_start: i64[>..] = [];
var edges_at: i64[>..] = [];

fn find_edges() {
    edge_index.clear();
    ea.clear();
    eb.clear();
    elong.clear();
    for t in btets.len / 4 {
        let v = btets[4 * t..4 * t + 4];
        let a, b, longest = longest_edge(lx, v);
        for i in 4 {
            for j in i + 1..4 {
                let e = edge_index.get_or_insert(edge_key(v[i], v[j]), ea.len);
                if e == ea.len {
                    ea.push(v[i]);
                    eb.push(v[j]);
                    elong.push(false);
                }
                if gap(lx, v[i], v[j]) >= longest * (1.0 - 1.0e-9) { elong[e] = true; }
            }
        }
    }
    let points = lx.len / 3;
    edge_start.clear();
    for p in points + 1 { edge_start.push(0); }
    for e in ea.len {
        edge_start[ea[e] + 1]++;
        edge_start[eb[e] + 1]++;
    }
    for p in points { edge_start[p + 1] += edge_start[p]; }
    edges_at.clear();
    for e in 2 * ea.len { edges_at.push(0); }
    var fill: i64[>..] = [];
    fill.append(edge_start[..points]);
    for e in ea.len {
        edges_at[fill[ea[e]]] = e;
        fill[ea[e]]++;
        edges_at[fill[eb[e]]] = e;
        fill[eb[e]]++;
    }
}

// Sharp edges (feature.goose): a curved or an oblique sharp edge of the
// shape is not on the edges of the lattice, so the stuffing chamfers it by
// about one third of a lattice cell. protect_edges puts it on a chain of
// lattice edges instead. For each curve it takes the lattice points near the
// curve, with each lattice edge between 2 of them that goes about as far
// along the curve as it is long, and finds the chain of such edges from one
// end of the curve to the other (or round a closed curve) whose points are
// the nearest to the curve. Then it moves each point of the chain onto the
// curve and makes its φ 0. The stuffing keeps each edge between 2 points
// with φ = 0, so the mesh has the edges of the chain, and its boundary
// triangles on the 2 sides of it are the cuts of the 2 surfaces. The
// lattice cells keep their orientation and a γ above MIN_GAMMA, or the
// chain is not used and the next, thinner tube is tried.
let TUBES: f64[4] = [0.3, 0.45, 0.6, 0.8];   // the radius of the tube about a curve, in lattice edges
let MIN_GAMMA = 0.03;
let STEP_LOW = 0.25;                    // the least and the largest advance along the curve, in edge lengths
let STEP_HIGH = 1.35;
let MAX_TRIES = 8;                      // chains for a curve, each without the points of the cells of the last
let FLAT_GAMMA = 0.01;                   // a stuffed cell with a smaller γ shows a chain that made a flat cell
let RETRY_SCALES = [1.0, 0.97, 1.03, 0.94, 1.06];
let RELAX_DEPTH = 1;                   // the free lattice points this many edges from the chain move with it
let RELAX_ROUNDS = 6;
let RELAX = 0.6;
let SUPPORT = 1.05;                     // the margin of the cut of warp
let END_REACH = 0.6;                    // a chain starts and ends this near the ends of an open curve

var keep_edges = true;          // false: mesh_3d does not protect the sharp edges
var curves_found = 0;           // the sharp edges of the last mesh_3d, and those not on a chain
var curves_missed = 0;
var chain_a: i64[>..] = [];             // the lattice edges of the chains, from lattice point a to b
var chain_b: i64[>..] = [];
var chain_owner: i64[>..] = []; // the curve whose chain has the point, or -1

// The mean length of the lattice edges at point p.
fn node_edge_mean(p: i64) -> f64 {
    var sum = 0.0;
    var count = 0.0;
    for e in edges_at[edge_start[p]..edge_start[p + 1]] {
        sum += gap(lx[..], ea[e], eb[e]);
        count += 1.0;
    }
    return sum / count;
}

fn lattice_at(p: i64) -> f64[3] { return [lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]]; }

// Whether x is on 2 planes of the lattice, so on a lattice edge.
fn on_lattice_edge(x: f64[3], d: f64) -> bool {
    var planes = 0;
    for a in 3 {
        for i in line_start[a]..line_start[a + 1] {
            if abs(x[a] - lattice_lines[i]) <= 1.0e-9 * d {
                planes++;
                break;
            }
        }
    }
    return planes >= 2;
}

// The position and the unit normal of the point of the surface on the
// lattice edge a b (a inside, b outside), or on the lattice point a = b.
// Made when first used, in `slots` (6 numbers each).
fn surface_slot(g: i64, key: i64, a: i64, b: i64, d: f64, slot_of: dictionary<i64, i64>&, slots: f64[>..]&) -> i64 {
    let have = slot_of.get(key);
    if have { return have; }
    var x = lattice_at(a);
    if a != b {
        let t = cut_fraction_steps(g, a, b, 26);   // 1.5e-8 of the edge: enough for a normal
        x = axpy3(x, t, sub3(lattice_at(b), x));
    }
    // Forward differences: 4 evaluations of φ, not 6, as this runs for each surface point.
    let eps = FD * d;
    let f0 = phi(g, x[0], x[1], x[2]);
    let n: f64[3] = [(phi(g, x[0] + eps, x[1], x[2]) - f0) / eps, (phi(g, x[0], x[1] + eps, x[2]) - f0) / eps,
                     (phi(g, x[0], x[1], x[2] + eps) - f0) / eps];
    let norm = norm3(n);
    let slot = slots.len / 6;
    slots.append(x);
    slots.append(if norm > 0.0 { [n[0] / norm, n[1] / norm, n[2] / norm] } else { [0.0, 0.0, 0.0] });
    slot_of.insert(key, slot);
    return slot;
}

// Finds the sharp edges of g in the lattice. A lattice tetrahedron with 2
// points of the surface (a lattice point with φ = 0, or the cut point of
// a lattice edge) whose normals differ by more than PAIR_ANGLE may have an edge
// between them. The edges on 2 planes of the lattice have no error, and
// are not curves.
fn find_curves(g: i64, size: f64) {
    clear_curves();
    let d = shape_scale(g);
    var slot_of = dictionary<i64, i64> {};
    var slots: f64[>..] = [];
    for t in btets.len / 4 {
        let v = btets[4 * t..4 * t + 4];
        var found: i64[10] = [0; 10];
        var n = 0;
        for i in 4 {
            if lphi[v[i]] == 0.0 {
                found[n] = surface_slot(g, -(v[i] + 1), v[i], v[i], d, slot_of, slots);
                n++;
            }
        }
        for i in 4 {
            for j in i + 1..4 {
                guard lphi[v[i]] * lphi[v[j]] < 0.0 else { continue; }
                let inside = if lphi[v[i]] < 0.0 { v[i] } else { v[j] };
                let outside = if lphi[v[i]] < 0.0 { v[j] } else { v[i] };
                found[n] = surface_slot(g, edge_key(v[i], v[j]), inside, outside, d, slot_of, slots);
                n++;
            }
        }
        var best = PAIR_COS;
        var bi = -1;
        var bj = -1;
        for i in n {
            for j in i + 1..n {
                let c = slots[6 * found[i] + 3] * slots[6 * found[j] + 3] + slots[6 * found[i] + 4] * slots[6 * found[j] + 4] +
                        slots[6 * found[i] + 5] * slots[6 * found[j] + 5];
                if c < best && c > -1.5 {
                    best = c;
                    bi = i;
                    bj = j;
                }
            }
        }
        guard bi >= 0 else { continue; }
        let u: f64[3] = [slots[6 * found[bi]], slots[6 * found[bi] + 1], slots[6 * found[bi] + 2]];
        let nu: f64[3] = [slots[6 * found[bi] + 3], slots[6 * found[bi] + 4], slots[6 * found[bi] + 5]];
        let w: f64[3] = [slots[6 * found[bj]], slots[6 * found[bj] + 1], slots[6 * found[bj] + 2]];
        let nw: f64[3] = [slots[6 * found[bj] + 3], slots[6 * found[bj] + 4], slots[6 * found[bj] + 5]];
        guard abs(dot3(nu, nu) - 1.0) < 1.0e-6 && abs(dot3(nw, nw) - 1.0) < 1.0e-6 else { continue; }
        let mid = axpy3(u, 0.5, sub3(w, u));
        let local = wanted_size(mid[0], mid[1], mid[2], 0.0, size);
        guard !near_curve(mid, local) else { continue; }
        let r = find_ridge(g, u, nu, w, nw, d);
        guard r.found && !on_lattice_edge(r.x, d) && !near_curve(r.x, 0.75 * local / CURVE_DIV) else { continue; }
        trace_curve(g, r, size, d);
    }
}

// The cumulative length along curve c, one number for each of its segments
// and a first 0. A closed curve has a last segment to its first sample.
fn curve_arc(c: i64, arc: f64[>..]&) {
    arc.clear();
    arc.push(0.0);
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = if curve_closed[c] { n } else { n - 1 };
    for i in segs {
        arc.push(arc[i] + norm3(sub3(sample_at(first + (i + 1) % n), sample_at(first + i))));
    }
}

// The distance from x to curve c, and the length along the curve to the
// nearest point.
fn curve_foot(c: i64, arc: f64[:], x: f64[3]) -> f64, f64 {
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = arc.len - 1;
    var best = 1.0e300;
    var best_s = 0.0;
    for i in segs {
        let a = sample_at(first + i);
        let ab = sub3(sample_at(first + (i + 1) % n), a);
        let len2 = dot3(ab, ab);
        let lam = if len2 > 0.0 { clamp(dot3(sub3(x, a), ab) / len2, 0.0, 1.0) } else { 0.0 };
        let dist = norm3(sub3(x, axpy3(a, lam, ab)));
        if dist < best {
            best = dist;
            best_s = arc[i] + lam * sqrt(len2);
        }
    }
    return best, best_s;
}

// The point of curve c at the length s along it, found again on the
// sharp edge between 2 points on the sides of it. Its true point, or the
// point on the segment if the edge is not found there.
fn curve_exact(g: i64, c: i64, arc: f64[:], s: f64, d: f64) -> f64[3] {
    let first = curve_start[c];
    let n = curve_start[c + 1] - first;
    let segs = arc.len - 1;
    var i = 0;
    while i < segs - 1 && arc[i + 1] < s { i++; }
    let a = sample_at(first + i);
    let b = sample_at(first + (i + 1) % n);
    let len = arc[i + 1] - arc[i];
    let lam = if len > 0.0 { clamp((s - arc[i]) / len, 0.0, 1.0) } else { 0.0 };
    let y = axpy3(a, lam, sub3(b, a));
    guard len > 0.0 else { return y; }
    if !curve_closed[c] && (s <= 0.0 || s >= arc[segs]) { return y; }
    let na = sample_normal(first + i);
    let nb = sample_normal(first + (i + 1) % n);
    let nav = unit3([na[0] + nb[0], na[1] + nb[1], na[2] + nb[2]]);
    let t = unit3(sub3(b, a));
    let w = unit3(cross3(t, nav));
    let delta = 0.2 * len;
    let u, nu, ok_u = settle(g, axpy3(y, -delta, w), d);
    let v, nv, ok_v = settle(g, axpy3(y, delta, w), d);
    guard ok_u && ok_v else { return y; }
    let r = ridge_between(g, u, nu, v, nv, d);
    guard r.found && norm3(sub3(r.x, y)) <= 2.0 * delta else { return y; }
    return r.x;
}

// The cost of the node with index q, at the distance d from the curve.
fn node_cost(d: f64, l: f64) -> f64 { return (d / l) * (d / l); }

// The cost of an edge of length le that advances by `advance` along the curve.
fn step_cost(advance: f64, le: f64, l: f64) -> f64 { return (advance - le) * (advance - le) / (l * l); }

// The candidates: the lattice points within the tube of radius `tube` lattice
// edges about curve c, found from the lattice point nearest to its first
// sample through the lattice edges. For each: the lattice point, the length
// along the curve of its nearest point, its distance from the curve and the
// mean length of its edges.
fn find_candidates(c: i64, arc: f64[:], tube: f64, node: i64[>..]&, at_s: f64[>..]&, at_d: f64[>..]&, at_l: f64[>..]&) {
    node.clear();
    at_s.clear();
    at_d.clear();
    at_l.clear();
    let start = sample_at(curve_start[c]);
    var near = 0;
    var near_d = 1.0e300;
    for p in lx.len / 3 {
        let dist = norm3(sub3(lattice_at(p), start));
        if dist < near_d {
            near_d = dist;
            near = p;
        }
    }
    var seen: bool[>..] = [];
    for p in lx.len / 3 { seen.push(false); }
    var queue: i64[>..] = [near];
    seen[near] = true;
    var head = 0;
    while head < queue.len {
        let p = queue[head];
        head++;
        let l = node_edge_mean(p);
        let dist, s = curve_foot(c, arc, lattice_at(p));
        guard dist <= 1.2 * l else { continue; }
        if dist <= tube * l && chain_owner[p] == -1 {
            node.push(p);
            at_s.push(s);
            at_d.push(dist);
            at_l.push(l);
        }
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            if !seen[q] {
                seen[q] = true;
                queue.push(q);
            }
        }
    }
}

// The cells at each lattice point: tets_of[tet_start[p]..tet_start[p + 1]].
var tet_start: i64[>..] = [];
var tets_of: i64[>..] = [];

fn find_tets_of() {
    let points = lx.len / 3;
    tet_start.clear();
    for p in points + 1 { tet_start.push(0); }
    for t in btets.len / 4 {
        for i in 4 { tet_start[btets[4 * t + i] + 1]++; }
    }
    for p in points { tet_start[p + 1] += tet_start[p]; }
    tets_of.clear();
    for i in tet_start[points] { tets_of.push(0); }
    var fill: i64[>..] = [];
    fill.append(tet_start[..points]);
    for t in btets.len / 4 {
        for i in 4 {
            let p = btets[4 * t + i];
            tets_of[fill[p]] = t;
            fill[p]++;
        }
    }
}

// Whether lattice point m is inside and warp does not move it onto the
// surface: no edge of m to an outside point has its cut nearer to m than the
// fraction alpha of the edge. The cut is found from φ at the 2 ends as if
// φ were linear on the edge (with a margin of SUPPORT).
fn stays_inside(m: i64) -> bool {
    guard lphi[m] < 0.0 else { return false; }
    for e in edges_at[edge_start[m]..edge_start[m + 1]] {
        let q = if ea[e] == m { eb[e] } else { ea[e] };
        guard lphi[q] > 0.0 else { continue; }
        let t = -lphi[m] / (lphi[q] - lphi[m]);
        if t < SUPPORT * (if elong[e] { ALPHA_LONG } else { ALPHA_SHORT }) { return false; }
    }
    return true;
}

// Whether the lattice edge a b, with φ = 0 at both ends, is an edge of a
// cell of the mesh. The stuffing makes a cell on it if a lattice cell
// of the edge has a point inside that warp does not move onto the surface.
fn edge_supported(a: i64, b: i64) -> bool {
    for t in tets_of[tet_start[a]..tet_start[a + 1]] {
        var has_b = false;
        for i in 4 {
            if btets[4 * t + i] == b { has_b = true; }
        }
        guard has_b else { continue; }
        for i in 4 {
            let m = btets[4 * t + i];
            if m != a && m != b && stays_inside(m) { return true; }
        }
    }
    return false;
}

// The state of the search of best_chain from the start p0 (-1 for an open
// curve): the position of each candidate along the curve from the start, the
// candidates in that sequence, and the cost and the previous candidate.
var dp_pos: f64[>..] = [];
var dp_order: i64[>..] = [];
var dp_cost: f64[>..] = [];
var dp_prev: i64[>..] = [];

fn chain_start(count: i64, closed: bool, p0: i64, total: f64, at_s: f64[:], at_d: f64[:], at_l: f64[:]) {
    // Positions from the start, in sequence.
    dp_pos.clear();
    dp_order.clear();
    dp_cost.clear();
    dp_prev.clear();
    var keys: i64[>..] = [];
    for k in count {
        let origin = if closed { at_s[p0] } else { 0.0 };
        var x = at_s[k] - origin;
        if x < 0.0 { x += total; }
        dp_pos.push(x);
        dp_cost.push(1.0e300);
        dp_prev.push(-1);
        keys.push(trunc(x / total * 1.0e9) as i64 * 4194304 + k);
    }
    sort(keys);
    for k in count { dp_order.push(keys[k] % 4194304); }
    if closed {
        dp_cost[p0] = node_cost(at_d[p0], at_l[p0]);
    } else {
        for k in count {
            if at_s[k] <= END_REACH * at_l[k] { dp_cost[k] = node_cost(at_d[k], at_l[k]); }
        }
    }
}

// The least cost of a chain to each candidate, over the lattice edges in the
// sequence along the curve.
fn chain_forward(count: i64, closed: bool, p0: i64, node: i64[:], index_of: dictionary<i64, i64>&, at_d: f64[:], at_l: f64[:]) {
    for oi in count {
        let q = dp_order[oi];
        for e in edges_at[edge_start[node[q]]..edge_start[node[q] + 1]] {
            let pn = if ea[e] == node[q] { eb[e] } else { ea[e] };
            let have = index_of.get(pn);
            guard have else { continue; }
            let pk = have;
            guard dp_cost[pk] < 1.0e299 && dp_pos[pk] < dp_pos[q] else { continue; }
            let le = gap(lx[..], node[q], pn);
            let advance = dp_pos[q] - dp_pos[pk];
            guard advance >= STEP_LOW * le && advance <= STEP_HIGH * le else { continue; }
            if closed && q == p0 { continue; }
            guard edge_supported(node[pk], node[q]) else { continue; }
            let l = 0.5 * (at_l[q] + at_l[pk]);
            let sum = dp_cost[pk] + step_cost(advance, le, l) + node_cost(at_d[q], at_l[q]);
            if sum < dp_cost[q] {
                dp_cost[q] = sum;
                dp_prev[q] = pk;
            }
        }
    }
}

// The best chain of lattice edges along curve c from the candidates, as the
// index in the candidates of each point of it, in sequence. False if there
// is none. `total` is the length of the curve.
fn best_chain(c: i64, total: f64, node: i64[:], at_s: f64[:], at_d: f64[:], at_l: f64[:], path: i64[>..]&) -> bool {
    path.clear();
    let count = node.len;
    guard count >= 2 else { return false; }
    let closed = curve_closed[c];
    var index_of = dictionary<i64, i64> {};
    for k in count { index_of.insert(node[k], k); }
    var starts: i64[>..] = [];
    if closed {
        for k in count {
            if at_s[k] < 2.5 * at_l[k] || at_s[k] > total - 0.5 * at_l[k] { starts.push(k); }
        }
        guard starts.len > 0 else { return false; }
    } else {
        starts.push(-1);
    }
    var best_total = 1.0e300;
    var best_prev: i64[>..] = [];
    var best_last = -1;
    for p0 in starts {
        chain_start(count, closed, p0, total, at_s, at_d, at_l);
        chain_forward(count, closed, p0, node, index_of, at_d, at_l);
        // The end of the chain.
        for q in count {
            guard dp_cost[q] < 1.0e299 else { continue; }
            var sum = dp_cost[q];
            if closed {
                guard q != p0 else { continue; }
                var closes = false;
                for e in edges_at[edge_start[node[q]]..edge_start[node[q] + 1]] {
                    if ea[e] == node[p0] || eb[e] == node[p0] { closes = true; }
                }
                let le = gap(lx[..], node[q], node[p0]);
                let advance = total - dp_pos[q];
                guard closes && advance >= STEP_LOW * le && advance <= STEP_HIGH * le && edge_supported(node[q], node[p0]) else { continue; }
                sum += step_cost(advance, le, 0.5 * (at_l[q] + at_l[p0]));
            } else {
                guard total - at_s[q] <= END_REACH * at_l[q] else { continue; }
            }
            if sum < best_total {
                best_total = sum;
                best_last = q;
                best_prev.clear();
                best_prev.append(dp_prev);
            }
        }
    }
    guard best_last >= 0 else { return false; }
    var rev: i64[>..] = [];
    var k = best_last;
    while k >= 0 {
        rev.push(k);
        k = best_prev[k];
    }
    for i in rev.len { path.push(rev[rev.len - 1 - i]); }
    return true;
}

// Whether the lattice point p stays where it is when the chain moves: it is
// on the surface, on a plane of a box of the shape (`planes`: the axis and
// the coordinate of each), or on the chain of another curve.
fn pinned(p: i64, planes: f64[:]) -> bool {
    if lphi[p] == 0.0 || chain_owner[p] >= 0 { return true; }
    for k in planes.len / 2 {
        if abs(lx[3 * p + planes[2 * k] as i64] - planes[2 * k + 1]) <= 1.0e-12 { return true; }
    }
    return false;
}

// The free lattice points near the chain (RELAX_DEPTH edges): not on a plane
// of `planes` and not on the chain.
fn near_chain(node: i64[:], path: i64[:], planes: f64[:], free: i64[>..]&) {
    var queue: i64[>..] = [];
    var depth: i64[>..] = [];
    var seen = dictionary<i64, i64> {};
    for k in path.len {
        queue.push(node[path[k]]);
        depth.push(0);
        seen.insert(node[path[k]], 0);
    }
    var head = 0;
    while head < queue.len {
        let p = queue[head];
        let dp = depth[head];
        head++;
        guard dp < RELAX_DEPTH else { continue; }
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            guard !seen.get(q) else { continue; }
            seen.insert(q, dp + 1);
            queue.push(q);
            depth.push(dp + 1);
            if !pinned(q, planes) { free.push(q); }
        }
    }
}

// Moves each free point part of the way to the mean of its neighbours, and
// keeps it on its side of the surface: a point that would cross it does not
// move. `free_old` and `free_phi` are the positions and φ before.
fn relax_free(g: i64, free: i64[:], free_old: f64[:], free_phi: f64[:]) {
    for round in RELAX_ROUNDS {
        for p in free {
            var mean: f64[3] = [0.0, 0.0, 0.0];
            var count = 0.0;
            for e in edges_at[edge_start[p]..edge_start[p + 1]] {
                let q = if ea[e] == p { eb[e] } else { ea[e] };
                for a in 3 { mean[a] += lx[3 * q + a]; }
                count += 1.0;
            }
            for a in 3 { lx[3 * p + a] += RELAX * (mean[a] / count - lx[3 * p + a]); }
        }
    }
    for k in free.len {
        let p = free[k];
        let f = phi(g, lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        if f == 0.0 || (f < 0.0) != (free_phi[k] < 0.0) {
            for a in 3 { lx[3 * p + a] = free_old[3 * k + a]; }
        } else {
            lphi[p] = f;
        }
    }
}

// Whether the lattice cells `cells` (with the orientation `sign`) are all
// positive with γ of at least MIN_GAMMA. The point of the chain in the worst
// cell that moved the most is not a candidate for the next try.
fn cells_valid(cells: i64[:], sign: f64[:], marked: dictionary<i64, i64>&, target: f64[:], old: f64[:]) -> bool {
    var ok = true;
    for i in cells.len {
        let v = btets[4 * cells[i]..4 * cells[i] + 4];
        let q = if sign[i] < 0.0 { [v[1], v[0], v[2], v[3]] } else { [v[0], v[1], v[2], v[3]] };
        guard gamma_of(lattice_at(q[0]), lattice_at(q[1]), lattice_at(q[2]), lattice_at(q[3])) < MIN_GAMMA else { continue; }
        ok = false;
        var far = -1;
        var far_move = -1.0;
        for n in v {
            let k = marked.get(n);
            if k {
                let move = norm3(sub3([target[3 * k], target[3 * k + 1], target[3 * k + 2]],
                                      [old[3 * k], old[3 * k + 1], old[3 * k + 2]]));
                if move > far_move {
                    far_move = move;
                    far = n;
                }
            }
        }
        if far >= 0 { chain_owner[far] = -2; }
    }
    return ok;
}

// Moves the lattice points of the chain onto the curve, and the free
// lattice points near them part of the way after them, by smoothing, so that
// the cells next to the chain keep their shape. False, with everything as
// before, if a lattice cell then has an inverted orientation or a γ below
// MIN_GAMMA.
fn snap_chain(g: i64, c: i64, arc: f64[:], d: f64, node: i64[:], at_s: f64[:], path: i64[:], planes: f64[:]) -> bool {
    let total = arc[arc.len - 1];
    var target: f64[>..] = [];
    for k in path.len {
        var s = at_s[path[k]];
        if !curve_closed[c] {
            if k == 0 { s = 0.0; }
            if k == path.len - 1 { s = total; }
        }
        target.append(curve_exact(g, c, arc, s, d));
    }
    var marked = dictionary<i64, i64> {};
    for k in path.len { marked.insert(node[path[k]], k); }
    var free: i64[>..] = [];
    near_chain(node, path, planes, free);
    // The cells that move, and their orientation.
    var cells: i64[>..] = [];
    var sign: f64[>..] = [];
    var moves = dictionary<i64, i64> {};
    for p in free { moves.insert(p, 0); }
    for t in btets.len / 4 {
        var touched = false;
        for i in 4 {
            if marked.get(btets[4 * t + i]) || moves.get(btets[4 * t + i]) { touched = true; }
        }
        guard touched else { continue; }
        cells.push(t);
        sign.push(if lattice_volume(t) < 0.0 { -1.0 } else { 1.0 });
    }
    var old: f64[>..] = [];
    var old_phi: f64[>..] = [];
    for k in path.len {
        let p = node[path[k]];
        old.append(lattice_at(p));
        old_phi.push(lphi[p]);
        for a in 3 { lx[3 * p + a] = target[3 * k + a]; }
        lphi[p] = 0.0;
    }
    var free_old: f64[>..] = [];
    var free_phi: f64[>..] = [];
    for p in free {
        free_old.append(lattice_at(p));
        free_phi.push(lphi[p]);
    }
    relax_free(g, free[..], free_old[..], free_phi[..]);
    let ok = cells_valid(cells[..], sign[..], marked, target[..], old[..]);
    if !ok {
        for k in path.len {
            let p = node[path[k]];
            for a in 3 { lx[3 * p + a] = old[3 * k + a]; }
            lphi[p] = old_phi[k];
        }
        for k in free.len {
            for a in 3 { lx[3 * free[k] + a] = free_old[3 * k + a]; }
            lphi[free[k]] = free_phi[k];
        }
    }
    return ok;
}

fn lattice_volume(t: i64) -> f64 {
    let v = btets[4 * t..4 * t + 4];
    let a = lattice_at(v[0]);
    return dot3(sub3(lattice_at(v[1]), a), cross3(sub3(lattice_at(v[2]), a), sub3(lattice_at(v[3]), a)));
}

// Protects the sharp edges of g. It returns the number of curves that it
// did not put on a chain of lattice edges.
fn protect_edges(g: i64, size: f64) -> i64 {
    curves_found = 0;
    curves_missed = 0;
    chain_a.clear();
    chain_b.clear();
    find_curves(g, size);
    chain_owner.clear();
    for p in lx.len / 3 { chain_owner.push(-1); }
    guard curve_count() > 0 else { return 0; }
    find_tets_of();
    let d = shape_scale(g);
    var missed = 0;
    var arc: f64[>..] = [];
    var node: i64[>..] = [];
    var at_s: f64[>..] = [];
    var at_d: f64[>..] = [];
    var at_l: f64[>..] = [];
    var path: i64[>..] = [];
    var boxes: f64[>..] = [];
    box_leaves(g, boxes);
    var planes: f64[>..] = [];
    for k in boxes.len / 6 {
        for a in 3 {
            planes.append([a, boxes[6 * k + a], a, boxes[6 * k + 3 + a]]);
        }
    }
    for c in curve_count() {
        curve_arc(c, arc);
        var done = false;
        for tube in TUBES {
            for attempt in MAX_TRIES {
                find_candidates(c, arc[..], tube, node, at_s, at_d, at_l);
                guard best_chain(c, arc[arc.len - 1], node[..], at_s[..], at_d[..], at_l[..], path) else { break; }
                if snap_chain(g, c, arc[..], d, node[..], at_s[..], path[..], planes[..]) {
                    for k in path { chain_owner[node[k]] = c; }
                    for k in path.len - 1 {
                        chain_a.push(node[path[k]]);
                        chain_b.push(node[path[k + 1]]);
                    }
                    if curve_closed[c] {
                        chain_a.push(node[path[path.len - 1]]);
                        chain_b.push(node[path[0]]);
                    }
                    done = true;
                    break;
                }
            }
            if done { break; }
        }
        // The points that this curve could not use are free for the next.
        for p in chain_owner.len {
            if chain_owner[p] == -2 { chain_owner[p] = -1; }
        }
        if !done { missed++; }
    }
    curves_found = curve_count();
    curves_missed = missed;
    return missed;
}

// Warping: moves each lattice point that a cut point violates onto the
// nearest such cut point.
fn warp(g: i64) {
    for p in lx.len / 3 {
        guard lphi[p] != 0.0 else { continue; }
        var best = 2.0;
        var best_x: f64[3] = [0.0, 0.0, 0.0];
        for e in edges_at[edge_start[p]..edge_start[p + 1]] {
            let q = if ea[e] == p { eb[e] } else { ea[e] };
            guard lphi[q] != 0.0 && (lphi[q] < 0.0) != (lphi[p] < 0.0) else { continue; }
            let inside = if lphi[p] < 0.0 { p } else { q };
            let outside = if lphi[p] < 0.0 { q } else { p };
            let t = cut_fraction(g, inside, outside);
            let from_p = if inside == p { t } else { 1.0 - t };
            let alpha = if elong[e] { ALPHA_LONG } else { ALPHA_SHORT };
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
}

// The lattice tetrahedra: for each pair of cubes with a common face, the 2
// centres and each edge of the face.
fn lattice_tets() {
    btets.clear();
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
                    for q in 4 { btets.append([c, c2, f[q], f[(q + 1) % 4]]); }
                }
            }
        }
    }
}

// The lattice points of shape g of bounds b: the planes of grid_lines on
// each axis, with one cube more on each side. False for a flat shape.
fn build_lattice(g: i64, b: f64[6], size: f64) -> bool {
    var boxes: f64[>..] = [];
    box_leaves(g, boxes);
    lattice_lines.clear();
    var line: f64[>..] = [];
    for a in 3 {
        guard b[3 + a] > b[a] else { return false; }
        grid_lines(boxes[..], a, b[a], b[3 + a], size, line);
        line_start[a] = lattice_lines.len;
        lattice_lines.push(2.0 * line[0] - line[1]);
        lattice_lines.append(line);
        lattice_lines.push(2.0 * line[line.len - 1] - line[line.len - 2]);
        lattice_cubes[a] = line.len + 1;
    }
    line_start[3] = lattice_lines.len;
    lx.clear();
    lphi.clear();
    lnode.clear();
    cuts.clear();
    on_surface.clear();
    for p in lattice_point_count() {
        let at = half_position(p);
        for a in 3 {
            let lo = lattice_lines[line_start[a] + at[a] / 2];
            lx.push(if at[a] % 2 == 0 { lo } else { 0.5 * (lo + lattice_lines[line_start[a] + at[a] / 2 + 1]) });
        }
    }
    return true;
}

// Meshes the solid g with tetrahedra of edge length about `size`, and
// smaller at the targets of refine.goose. False if the refinement does not
// finish.
fn has_bad_cell() -> bool {
    for c in cell_count() {
        let a = vertex(cell_node(c, 0));
        let b = vertex(cell_node(c, 1));
        let d = vertex(cell_node(c, 2));
        let e = vertex(cell_node(c, 3));
        if gamma_of(a, b, d, e) < FLAT_GAMMA { return true; }
    }
    return false;
}

// The lattice of shape g of bounds b for the size `scale`, refined, with φ
// at each point, and the sharp lattice edges. False if it cannot be built.
fn lattice_pass(g: i64, b: f64[6], scale: f64) -> bool {
    clear_mesh(3);
    guard build_lattice(g, b, scale) else { return false; }
    lattice_tets();
    guard refine_cells(lx, btets, 4, scale) else { return false; }
    var h = scale;
    for a in 3 {
        for i in lattice_cubes[a] { h = min(h, lattice_lines[line_start[a] + i + 1] - lattice_lines[line_start[a] + i]); }
    }
    for t in targets { h = min(h, t.size); }
    lphi.clear();
    lnode.clear();
    for p in lx.len / 3 {
        let f = phi(g, lx[3 * p], lx[3 * p + 1], lx[3 * p + 2]);
        lphi.push(if abs(f) <= 1.0e-9 * h { 0.0 } else { f });
        lnode.push(-1);
    }
    find_edges();
    return true;
}

// Whether each lattice edge of the chains is an edge of the mesh.
fn chains_in_mesh() -> bool {
    let n = point_count();
    var want = dictionary<i64, i64> {};
    for k in chain_a.len {
        let a = lnode[chain_a[k]];
        let b = lnode[chain_b[k]];
        guard a >= 0 && b >= 0 else { return false; }
        want.insert(min(a, b) * n + max(a, b), 0);
    }
    for c in cell_count() {
        for i in 4 {
            for j in i + 1..4 {
                let a = cell_node(c, i);
                let b = cell_node(c, j);
                let key = min(a, b) * n + max(a, b);
                if want.get(key) { want.update(key, 0) { it = 1; }; }
            }
        }
    }
    var missing = 0;
    want.each() { key, v => if v == 0 { missing++; } };
    return missing == 0;
}

fn mesh_3d(g: i64, size: f64) -> bool {
    clear_mesh(3);
    if shell_tets(g, size) { return true; }
    let b = bounds(g);
    if has_round(g) && targets_span_z(b[2], b[5]) {
        let s = section(g, b[2], b[5]);
        if s >= 0 { return sweep_tets(s, b[2], b[5], size); }
    }
    // A prism along x or y: turn the shapes so that its axis is z, sweep,
    // and turn the shapes and the nodes back.
    for shift in [1, 2] {
        guard has_round(g) && turnable() else { break; }
        let count = shapes.len;
        turn_shapes(count, shift);
        let t = bounds(g);
        let s = if targets_span_z(t[2], t[5]) { section(g, t[2], t[5]) } else { -1 };
        let done = s >= 0 && sweep_tets(s, t[2], t[5], size);
        turn_shapes(count, 3 - shift);
        if s >= 0 {
            turn_nodes(3 - shift);
            return done;
        }
    }
    // A curve that no chain of lattice edges fits at this size may fit when
    // the lattice is a little finer or coarser, so the lattice is built again
    // with a size that is up to a few percent different. The last pass does
    // not protect the curves: it is the mesh without them, which is also the
    // result when a chain makes a cell with γ below FLAT_GAMMA, leaves a
    // chain edge out of the mesh, or gives a mesh that mesh_quality fails.
    let passes = RETRY_SCALES.len + 1;
    for pass in passes {
        let protect = keep_edges && pass < RETRY_SCALES.len;
        let scale = if protect { size * RETRY_SCALES[pass] } else { size };
        guard lattice_pass(g, b, scale) else { return false; }
        if protect {
            let missed = protect_edges(g, scale);
            if missed > 0 && (pass + 1 < RETRY_SCALES.len || missed == curves_found) { continue; }
        }
        warp(g);
        for t in btets.len / 4 { stuff(g, [btets[4 * t], btets[4 * t + 1], btets[4 * t + 2], btets[4 * t + 3]]); }
        if protect && curves_found > 0 && (has_bad_cell() || !chains_in_mesh()) {
            curves_missed = curves_found;
            continue;
        }
        if improve_cells { improve(on_surface[..]); }
        if protect && curves_found > 0 && !quality_passes(mesh_quality()) {
            curves_missed = curves_found;
            continue;
        }
        break;
    }
    return true;
}
