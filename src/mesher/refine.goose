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
    axis: i64,
    planes: f64[2],
    plane_count: i64,
    round: bool,
    ring: f64[3],
}

var targets: Target[>..] = [];
let DEFAULT_GRADE = 1.3;
let DEFAULT_SLOPE = 0.4;
var grade = DEFAULT_GRADE;
var grade_given = false;
var section_z = 0.0;

fn clear_refine() {
    targets.clear();
    grade = DEFAULT_GRADE;
    grade_given = false;
}

fn set_grade(g: f64) {
    guard g > 1.0 else { abort("set_grade: g > 1"); }
    grade = g;
    grade_given = true;
}

fn add_refine(g: i64, size: f64, surface: bool) {
    guard g >= 0 && g < shapes.len && size > 0.0 else { abort("add_refine: a shape and a size > 0"); }
    targets.push(Target { g: g, size: size, surface: surface, axis: -1, planes: [0.0; 2], plane_count: 0, round: false, ring: [0.0; 3] });
}

fn add_refine_clip(g: i64, size: f64, surface: bool, axis: i64, planes: f64[2], count: i64, round: bool, ring: f64[3]) {
    guard g >= 0 && g < shapes.len && size > 0.0 && axis >= 0 && axis < 3 && count >= 1 && count <= 2 && (!round || ring[2] > 0.0) else { abort("add_refine_clip: a shape, a size > 0, an axis, 1 or 2 planes and a ring radius > 0"); }
    targets.push(Target { g: g, size: size, surface: surface, axis: axis, planes: planes, plane_count: count, round: round, ring: ring });
}

fn wanted_size(x: f64, y: f64, z: f64, reach: f64, size: f64) -> f64 {
    let slope = if grade_given { grade - 1.0 } else { DEFAULT_SLOPE };
    var h = size;
    for t in targets {
        let f = phi(t.g, x, y, z);
        var d = if t.surface { abs(f) } else { max(f, 0.0) };
        if t.axis >= 0 {
            let c: f64[3] = [x, y, z];
            var near = 1.0e300;
            for i in t.plane_count { near = min(near, abs(c[t.axis] - t.planes[i])); }
            var sq = d * d + near * near;
            if t.round {
                let ring = sqrt((x - t.ring[0]) * (x - t.ring[0]) + (y - t.ring[1]) * (y - t.ring[1])) - t.ring[2];
                sq += ring * ring;
            }
            d = sqrt(sq);
        }
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

fn hanging(mids: dictionary<i64, i64>&, v: i64[:]) -> bool {
    for i in v.len {
        for j in i + 1..v.len {
            if mids.get(edge_key(v[i], v[j])) { return true; }
        }
    }
    return false;
}

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
            var reach = 0.0;
            for n in v {
                var s = 0.0;
                for t in 3 { s += (xyz[3 * n + t] - centre[t]) * (xyz[3 * n + t] - centre[t]); }
                reach = max(reach, sqrt(s));
            }
            guard l > SLACK * wanted_size(centre[0], centre[1], centre[2], reach, size) || hanging(mids, v) else {
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
