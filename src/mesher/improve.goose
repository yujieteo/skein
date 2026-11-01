// The improvement of the bad tetrahedra of the stuffed mesh (spec §9.2),
// after the schedule of Marot and Remacle (2020): first the cheap
// operations, smoothing of a vertex and then edge removal and the 2-3
// flip, on each cell with γ < BAD_GAMMA, again until no operation helps.
//
// γ is 2 √6 times the inradius over the longest edge: 1 for a regular
// tetrahedron, 1 / max_aspect of quality.goose. Each operation changes a
// small cavity of cells, and the mesh keeps the change only if the worst γ
// of the new cells is more than that of the old cells. Thus the worst cell
// of the mesh never gets worse. A node on the surface does not move, and
// no operation changes a face on the surface, so the boundary and its
// sharp edges stay as the stuffing made them.
//
// 1. Smoothing moves a free vertex of the bad cell to the mean of its
//    neighbours, or to the apex of a regular tetrahedron on the face
//    opposite to it in its worst cell, or part of the way to one of
//    them, if that makes its worst cell better.
// 2. Edge removal (Klingner and Shewchuk 2007) replaces the m cells about
//    an inside edge a b by the 2 (m - 2) cells of a triangulation of the
//    ring of the m other points, joined to a and to b. A dynamic program
//    over the triangulations of the ring finds the one with the best
//    worst cell. m = 3 is the 3-2 flip.
// 3. The 2-3 flip replaces the 2 cells on an inside face p q r by the 3
//    cells about the edge between their other points d and e.
import std;
import math;
namespace mesher;
import .mesh;
import .quality;
import runtime from "../runtime";

let BAD_GAMMA = 0.35;
let MAX_PASSES = 20;
let MAX_RING = 7;
let GAIN = 1.0e-6;
let CLIMB_START = 0.1;
let CLIMB_LEVELS = 4;
let CLIMB_TRIES = 4;     // the least rise of the worst γ that counts as better

var improve_cells = true;   // false: mesh_3d keeps the stuffed cells

var ix: f64[>..] = [];      // x, y and z of each node
var it: i64[>..] = [];      // the 4 nodes of each cell, positive
var alive: bool[>..] = [];
var fixed: bool[>..] = [];  // the node is on the surface
var locked: bool[>..] = []; // an operation of this pass changed the cells of the node
// The cells at each node of a bad cell: star[star_start[s]..star_start[s + 1]]
// with s = star_of[node], or star_of[node] = -1.
var star_of: i64[>..] = [];
var star_start: i64[>..] = [];
var star: i64[>..] = [];

fn at(n: i64) -> f64[3] { return [ix[3 * n], ix[3 * n + 1], ix[3 * n + 2]]; }

// γ of the tetrahedron a b c d, or -1 if it is not positive.
fn gamma_of(a: f64[3], b: f64[3], c: f64[3], d: f64[3]) -> f64 {
    let u = sub3(b, a);
    let v = sub3(c, a);
    let w = sub3(d, a);
    let volume = dot3(u, cross3(v, w)) / 6.0;
    var longest = max(max(norm3(u), norm3(v)), norm3(w));
    longest = max(longest, max(max(norm3(sub3(c, b)), norm3(sub3(d, b))), norm3(sub3(d, c))));
    guard volume > 1.0e-12 * longest * longest * longest else { return -1.0; }
    let faces = 0.5 * (norm3(cross3(u, v)) + norm3(cross3(v, w)) + norm3(cross3(w, u)) +
                       norm3(cross3(sub3(c, b), sub3(d, b))));
    return 2.0 * sqrt(6.0) * 3.0 * volume / (faces * longest);
}

fn gamma4(a: i64, b: i64, c: i64, d: i64) -> f64 { return gamma_of(at(a), at(b), at(c), at(d)); }

fn cell_gamma(t: i64) -> f64 { return gamma4(it[4 * t], it[4 * t + 1], it[4 * t + 2], it[4 * t + 3]); }

fn has(t: i64, n: i64) -> bool {
    return it[4 * t] == n || it[4 * t + 1] == n || it[4 * t + 2] == n || it[4 * t + 3] == n;
}

fn star_cells(n: i64) -> i64[:] {
    let s = star_of[n];
    guard s >= 0 else { abort("star_cells: a node of a bad cell"); }
    return star[star_start[s]..star_start[s + 1]];
}

fn add_cell(a: i64, b: i64, c: i64, d: i64) {
    it.append([a, b, c, d]);
    alive.push(true);
}

// Improves the tetrahedra of the mesh. `surface` marks the nodes on the
// surface. It returns the number of operations.
fn improve(surface: bool[:]) -> i64 {
    guard mesh_dim == 3 && surface.len == point_count() else { abort("improve: a 3D mesh and a mark for each node"); }
    ix.clear();
    ix.append(runtime::coords);
    it.clear();
    alive.clear();
    for c in cell_count() {
        guard runtime::cell_type[c] as i64 == runtime::TET4 else { return 0; }
        add_cell(cell_node(c, 0), cell_node(c, 1), cell_node(c, 2), cell_node(c, 3));
    }
    fixed.clear();
    fixed.append(surface);
    var done = 0;
    for pass in MAX_PASSES {
        let changes = improve_pass();
        done += changes;
        if changes == 0 { break; }
    }
    guard done > 0 else { return 0; }
    runtime::mesh_clear();
    for n in ix.len / 3 { add_point(ix[3 * n], ix[3 * n + 1], ix[3 * n + 2]); }
    for t in alive.len {
        if alive[t] { add_simplex(it[4 * t..4 * t + 4]); }
    }
    return done;
}

// One pass over the bad cells. It returns the number of operations.
fn improve_pass() -> i64 {
    var bad: i64[>..] = [];
    for t in alive.len {
        if alive[t] && cell_gamma(t) < BAD_GAMMA { bad.push(t); }
    }
    guard bad.len > 0 else { return 0; }
    find_stars(bad[..]);
    locked.clear();
    for n in ix.len / 3 { locked.push(false); }
    var changes = 0;
    for t in bad {
        guard alive[t] && cell_gamma(t) < BAD_GAMMA else { continue; }
        let v = it[4 * t..4 * t + 4];
        guard !locked[v[0]] && !locked[v[1]] && !locked[v[2]] && !locked[v[3]] else { continue; }
        if improve_cell(t) { changes++; }
    }
    return changes;
}

// The cells at each node of the cells `bad`.
fn find_stars(bad: i64[:]) {
    star_of.clear();
    for n in ix.len / 3 { star_of.push(-1); }
    var count = 0;
    for t in bad {
        for n in it[4 * t..4 * t + 4] {
            if star_of[n] < 0 {
                star_of[n] = count;
                count++;
            }
        }
    }
    star_start.clear();
    for s in count + 1 { star_start.push(0); }
    for t in alive.len {
        guard alive[t] else { continue; }
        for n in it[4 * t..4 * t + 4] {
            if star_of[n] >= 0 { star_start[star_of[n] + 1]++; }
        }
    }
    for s in count { star_start[s + 1] += star_start[s]; }
    star.clear();
    for i in star_start[count] { star.push(0); }
    var fill: i64[>..] = [];
    fill.append(star_start[..count]);
    for t in alive.len {
        guard alive[t] else { continue; }
        for n in it[4 * t..4 * t + 4] {
            let s = star_of[n];
            if s >= 0 {
                star[fill[s]] = t;
                fill[s]++;
            }
        }
    }
}

// Tries the operations on the bad cell t, in the sequence of the schedule.
fn improve_cell(t: i64) -> bool {
    let v: i64[4] = [it[4 * t], it[4 * t + 1], it[4 * t + 2], it[4 * t + 3]];
    for n in v {
        if !fixed[n] && smooth(n) { return true; }
    }
    for i in 4 {
        for j in i + 1..4 {
            if remove_edge(v[i], v[j]) { return true; }
        }
    }
    for i in 4 {
        if flip23(t, i) { return true; }
    }
    return false;
}

// The worst γ of the cells at n with n at x.
fn worst_at(n: i64, x: f64[3]) -> f64 {
    var worst = 1.0;
    for t in star_cells(n) {
        guard alive[t] else { continue; }
        var p: f64[3][4] = default<f64[3][4]>();
        for i in 4 {
            let m = it[4 * t + i];
            p[i] = if m == n { x } else { at(m) };
        }
        worst = min(worst, gamma_of(p[0], p[1], p[2], p[3]));
    }
    return worst;
}

// Moves the free node n to a better place, if there is one.
fn smooth(n: i64) -> bool {
    let x = at(n);
    let old = worst_at(n, x);
    // The mean of the neighbours, and the apex of the regular tetrahedron
    // on the face opposite to n in its worst cell.
    var mean: f64[3] = [0.0, 0.0, 0.0];
    var count = 0.0;
    var apex = x;
    var worst = 2.0;
    for t in star_cells(n) {
        guard alive[t] else { continue; }
        var face: i64[..3] = [];
        for i in 4 {
            let m = it[4 * t + i];
            if m != n {
                face.push(m);
                for c in 3 { mean[c] += ix[3 * m + c]; }
                count += 1.0;
            }
        }
        let g = cell_gamma(t);
        if g < worst {
            worst = g;
            apex = regular_apex(at(face[0]), at(face[1]), at(face[2]), x);
        }
    }
    for c in 3 { mean[c] /= count; }
    var best = old;
    var best_x = x;
    for target in [mean, apex] {
        for step in [1.0, 0.5, 0.25] {
            let y: f64[3] = [x[0] + step * (target[0] - x[0]), x[1] + step * (target[1] - x[1]), x[2] + step * (target[2] - x[2])];
            let q = worst_at(n, y);
            if q > best {
                best = q;
                best_x = y;
            }
        }
    }
    best_x, best = climb(n, best_x, best, edge_mean(n));
    guard best > old + GAIN else { return false; }
    for c in 3 { ix[3 * n + c] = best_x[c]; }
    locked[n] = true;
    return true;
}

// The mean length of the edges of the cells at n.
fn edge_mean(n: i64) -> f64 {
    var sum = 0.0;
    var count = 0.0;
    for t in star_cells(n) {
        guard alive[t] else { continue; }
        for i in 4 {
            let m = it[4 * t + i];
            if m != n {
                sum += norm3(sub3(at(m), at(n)));
                count += 1.0;
            }
        }
    }
    return sum / count;
}

// A pattern search from x, where the worst γ of the cells at n is q: steps
// along the axes of `CLIMB_START` times the edge length h, halved
// CLIMB_LEVELS times. It returns the best point and its worst γ.
fn climb(n: i64, x: f64[3], q: f64, h: f64) -> f64[3], f64 {
    var best_x = x;
    var best = q;
    var step = CLIMB_START * h;
    for level in CLIMB_LEVELS {
        var moved = true;
        for tries in CLIMB_TRIES {
            guard moved else { break; }
            moved = false;
            for d in 6 {
                var y = best_x;
                y[d / 2] += if d % 2 == 0 { step } else { -step };
                let r = worst_at(n, y);
                if r > best + GAIN {
                    best = r;
                    best_x = y;
                    moved = true;
                }
            }
        }
        step *= 0.5;
    }
    return best_x, best;
}

// The apex of the regular tetrahedron on the triangle a b c, on the side
// of x.
fn regular_apex(a: f64[3], b: f64[3], c: f64[3], x: f64[3]) -> f64[3] {
    var normal = cross3(sub3(b, a), sub3(c, a));
    let l = norm3(normal);
    guard l > 0.0 else { return x; }
    let centre: f64[3] = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0];
    if dot3(normal, sub3(x, centre)) < 0.0 { normal = [-normal[0], -normal[1], -normal[2]]; }
    let edge = (norm3(sub3(b, a)) + norm3(sub3(c, b)) + norm3(sub3(a, c))) / 3.0;
    let h = sqrt(2.0 / 3.0) * edge / l;
    return [centre[0] + h * normal[0], centre[1] + h * normal[1], centre[2] + h * normal[2]];
}

fn orient(a: i64, b: i64, c: i64, d: i64) -> f64 {
    let p = at(a);
    return dot3(sub3(at(b), p), cross3(sub3(at(c), p), sub3(at(d), p)));
}

// The ring of the inside edge a b: the other 2 points of each of its cells
// in a cycle, so that each a b ring[i] ring[i + 1] is positive, and the
// cells. False for an edge on the surface, or a ring of more than MAX_RING.
fn edge_ring(a: i64, b: i64, ring: i64[..8]&, cells: i64[..8]&) -> bool {
    var from: i64[..8] = [];
    var to: i64[..8] = [];
    for t in star_cells(a) {
        guard alive[t] && has(t, b) else { continue; }
        guard cells.len < MAX_RING else { return false; }
        cells.push(t);
        var o: i64[..2] = [];
        for i in 4 {
            let m = it[4 * t + i];
            if m != a && m != b { o.push(m); }
        }
        let flip = orient(a, b, o[0], o[1]) < 0.0;
        from.push(if flip { o[1] } else { o[0] });
        to.push(if flip { o[0] } else { o[1] });
    }
    guard cells.len >= 3 else { return false; }
    ring.push(from[0]);
    var next = to[0];
    for k in cells.len - 1 {
        guard next != ring[0] else { return false; }
        ring.push(next);
        var found = -1;
        for e in from.len {
            if from[e] == next { found = e; }
        }
        guard found >= 0 else { return false; }
        next = to[found];
    }
    return next == ring[0];
}

// Removes the inside edge a b, if a triangulation of its ring is better.
fn remove_edge(a: i64, b: i64) -> bool {
    var ring: i64[..8] = [];
    var cells: i64[..8] = [];
    guard edge_ring(a, b, ring, cells) else { return false; }
    var old = 1.0;
    for t in cells { old = min(old, cell_gamma(t)); }
    let m = ring.len;
    // best[i][j]: the best worst γ of the triangulations of ring i..j, and
    // split[i][j] the third point of its triangle on the side i j.
    var best: f64[8][8] = default<f64[8][8]>();
    var split: i64[8][8] = default<i64[8][8]>();
    for d in 2..m {
        for i in m - d {
            let j = i + d;
            best[i][j] = -2.0;
            for k in i + 1..j {
                let left = if k - i > 1 { best[i][k] } else { 1.0 };
                let right = if j - k > 1 { best[k][j] } else { 1.0 };
                let q = min(min(left, right), min(gamma4(a, ring[i], ring[k], ring[j]), gamma4(b, ring[i], ring[j], ring[k])));
                if q > best[i][j] {
                    best[i][j] = q;
                    split[i][j] = k;
                }
            }
        }
    }
    guard best[0][m - 1] > old + GAIN else { return false; }
    for t in cells { alive[t] = false; }
    var todo: i64[..16] = [0, m - 1];
    while todo.len > 0 {
        let j = todo[todo.len - 1];
        let i = todo[todo.len - 2];
        todo.pop();
        todo.pop();
        let k = split[i][j];
        add_cell(a, ring[i], ring[k], ring[j]);
        add_cell(b, ring[i], ring[j], ring[k]);
        if k - i > 1 { todo.append([i, k]); }
        if j - k > 1 { todo.append([k, j]); }
    }
    locked[a] = true;
    locked[b] = true;
    for n in ring { locked[n] = true; }
    return true;
}

// Flips the face of cell t opposite to its point `skip` into 3 cells about
// the edge between the other points of the 2 cells on it, if that is
// better.
fn flip23(t: i64, skip: i64) -> bool {
    let d = it[4 * t + skip];
    var f: i64[..3] = [];
    for i in 4 {
        if i != skip { f.push(it[4 * t + i]); }
    }
    if orient(f[0], f[1], f[2], d) < 0.0 {
        let first = f[0];
        f[0] = f[1];
        f[1] = first;
    }
    var other = -1;
    for u in star_cells(f[0]) {
        if u != t && alive[u] && has(u, f[1]) && has(u, f[2]) { other = u; }
    }
    guard other >= 0 else { return false; }
    var e = -1;
    for i in 4 {
        let m = it[4 * other + i];
        if m != f[0] && m != f[1] && m != f[2] { e = m; }
    }
    guard !locked[e] else { return false; }
    let old = min(cell_gamma(t), cell_gamma(other));
    var q = 1.0;
    for i in 3 { q = min(q, gamma4(f[i], f[(i + 1) % 3], e, d)); }
    guard q > old + GAIN else { return false; }
    alive[t] = false;
    alive[other] = false;
    for i in 3 { add_cell(f[i], f[(i + 1) % 3], e, d); }
    locked[d] = true;
    locked[e] = true;
    for n in f { locked[n] = true; }
    return true;
}
