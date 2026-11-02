// The local refinement of spec §9.1 (`refine r: size`): longest-edge
// bisection after Rivara (1984), on the background lattice of mesh_3d
// before it cuts the shape, and the wanted size of the triangles of
// mesh_2d.
//
// 1. Each target is a shape g with a size. A surface target wants that
//    size at |φ_g| = 0, and a region target in all of φ_g <= 0. Away from
//    a target, the wanted size grows by grade - 1 for each unit of
//    distance (the mesh setting grade, spec §9.1), up to the size of the
//    mesh.
// 2. Each round marks the cells whose longest edge is more than SLACK
//    times the least wanted size in the cell, and each cell with a point
//    at the midpoint of one of its edges. The least wanted size is that at
//    the distance `reach` nearer to the target than the centre, with reach
//    the distance from the centre to the farthest point of the cell. Each
//    marked cell is cut in two at the midpoint of its longest edge. Equal
//    lengths go to the edge of the smaller key, so two cells that share an
//    edge cut it the same way.
// 3. The rounds stop when no cell is marked. Then no edge holds a point at
//    its midpoint, and the cells are conforming.
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;

let SLACK = 1.5;
let MAX_ROUNDS = 60;
let MAX_POINTS = 4000000;

struct Target {
    g: i64,
    size: f64,
    surface: bool,
}

var targets: Target[>..] = [];
var grade = 1.3;       // the mesh setting grade (spec §9.1)
var grade_given = false;
var section_z = 0.0;   // the z at which the 2D mesher measures a target
var sweeping = false;  // true while the 2D mesher meshes the section of a 3D sweep

fn clear_refine() {
    targets.clear();
    grade = 1.3;
    grade_given = false;
}

fn set_grade(g: f64) {
    guard g > 1.0 else { abort("set_grade: g > 1"); }
    grade = g;
    grade_given = true;
}

// Wants the cells of edge `size` at the surface of g (surface) or in g.
fn add_refine(g: i64, size: f64, surface: bool) {
    guard g >= 0 && g < shapes.len && size > 0.0 else { abort("add_refine: a shape and a size > 0"); }
    targets.push(Target { g: g, size: size, surface: surface });
}

// The wanted size at a point at most `reach` from (x, y, z). With no
// grade given, the slope is 0.3 in 2D and 0.4 in 3D (spec §9.1).
fn wanted_size(x: f64, y: f64, z: f64, reach: f64, size: f64, solid: bool) -> f64 {
    let slope = if grade_given || !solid { grade - 1.0 } else { 0.4 };
    var h = size;
    for t in targets {
        let f = phi(t.g, x, y, z);
        let d = if t.surface { abs(f) } else { max(f, 0.0) };
        h = min(h, t.size + slope * max(d - reach, 0.0));
    }
    return h;
}

fn edge_key(a: i64, b: i64) -> i64 { return min(a, b) * 4294967296 + max(a, b); }

fn gap(xyz: f64[:], a: i64, b: i64) -> f64 {
    var s = 0.0;
    for t in 3 { s += (xyz[3 * a + t] - xyz[3 * b + t]) * (xyz[3 * a + t] - xyz[3 * b + t]); }
    return sqrt(s);
}

// The local indices of the longest edge of the cell `v`, and its length.
fn longest_edge(xyz: f64[:], v: i64[:]) -> i64, i64, f64 {
    var ia = 0;
    var ib = 1;
    var best = -1.0;
    for i in v.len {
        for j in i + 1..v.len {
            let l = gap(xyz, v[i], v[j]);
            if l > best || l == best && edge_key(v[i], v[j]) < edge_key(v[ia], v[ib]) {
                best = l;
                ia = i;
                ib = j;
            }
        }
    }
    return ia, ib, best;
}

// Whether the cell `v` has a point at the midpoint of an edge.
fn hanging(mids: dictionary<i64, i64>&, v: i64[:]) -> bool {
    for i in v.len {
        for j in i + 1..v.len {
            if mids.get(edge_key(v[i], v[j])) { return true; }
        }
    }
    return false;
}

// Refines the cells `cells` of k points each on the points `xyz` (x, y and
// z of each) to the targets, with the size `size` away from them. It adds
// the midpoints to `xyz`. False if the points would pass MAX_POINTS or the
// cells do not settle in MAX_ROUNDS rounds.
fn refine_cells(xyz: f64[>..]&, cells: i64[>..]&, k: i64, size: f64) -> bool {
    guard targets.len > 0 else { return true; }
    guard k >= 3 && k <= 4 && cells.len % k == 0 else { abort("refine_cells: triangles or tetrahedra"); }
    var mids = dictionary<i64, i64> {};
    var next: i64[>..] = [];
    for round in MAX_ROUNDS {
        next.clear();
        var cut = 0;
        for c in cells.len / k {
            let v = cells[k * c..k * (c + 1)];
            let ia, ib, l = longest_edge(xyz, v);
            var centre: f64[3] = [0.0, 0.0, 0.0];
            for n in v {
                for t in 3 { centre[t] += xyz[3 * n + t] / k; }
            }
            // Each point of the cell is at most `reach` from its centre.
            var reach = 0.0;
            for n in v {
                var s = 0.0;
                for t in 3 { s += (xyz[3 * n + t] - centre[t]) * (xyz[3 * n + t] - centre[t]); }
                reach = max(reach, sqrt(s));
            }
            guard l > SLACK * wanted_size(centre[0], centre[1], centre[2], reach, size, true) || hanging(mids, v) else {
                next.append(v);
                continue;
            }
            let key = edge_key(v[ia], v[ib]);
            var m = xyz.len / 3;
            let have = mids.get(key);
            if have {
                m = have;
            } else {
                for t in 3 { xyz.push(0.5 * (xyz[3 * v[ia] + t] + xyz[3 * v[ib] + t])); }
                mids.insert(key, m);
            }
            // Each child keeps the orientation: a point moves to a point
            // on one of its own edges.
            var first: i64[4] = [0, 0, 0, 0];
            var second: i64[4] = [0, 0, 0, 0];
            for i in k {
                first[i] = if i == ia { m } else { v[i] };
                second[i] = if i == ib { m } else { v[i] };
            }
            next.append(first[..k]);
            next.append(second[..k]);
            cut++;
        }
        cells.clear();
        cells.append(next);
        if xyz.len / 3 > MAX_POINTS { return false; }
        if cut == 0 { return true; }
    }
    return false;
}
