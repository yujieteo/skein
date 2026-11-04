import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .grid;
import .quad;
import runtime from "../runtime";

fn mesh_cells(g: i64, dim: i64, size: f64, planes: i64[:]) -> bool {
    if mesh_grid(g, dim, size, planes) { return true; }
    if mesh_sector(g, dim, size) { return true; }
    if dim == 2 { return mesh_blossom(g, size); }
    let b = bounds(g);
    let s = section(g, b[2], b[5]);
    if s >= 0 { return mesh_sweep(s, b[2], b[5], size); }
    mesh_mapped(g, size, planes);
    return true;
}

fn round_shape(g: i64, dim: i64) -> bool {
    let s = shapes[g];
    if dim == 2 { return s.kind == CIRCLE; }
    return s.kind == CYLINDER && abs(s.p[3]) < 1.0e-12 && abs(s.p[4]) < 1.0e-12;
}

fn sector_parts(g: i64) -> i64, i64, i64 {
    let s = shapes[g];
    if s.kind == INTERSECT {
        if shapes[s.a].kind == SUBTRACT { return shapes[s.a].a, shapes[s.a].b, s.b; }
        if shapes[s.b].kind == SUBTRACT { return shapes[s.b].a, shapes[s.b].b, s.a; }
    }
    if s.kind == SUBTRACT && shapes[s.a].kind == INTERSECT {
        let x = shapes[s.a];
        if shapes[x.b].kind == BOX || shapes[x.b].kind == RECT { return x.a, s.b, x.b; }
        return x.b, s.b, x.a;
    }
    return -1, -1, -1;
}

fn mesh_sector(g: i64, dim: i64, size: f64) -> bool {
    let o, i, c = sector_parts(g);
    guard o >= 0 && round_shape(o, dim) && round_shape(i, dim) else { return false; };
    guard shapes[c].kind == (if dim == 2 { RECT } else { BOX }) else { return false; };
    let po = shapes[o].p;
    let pn = shapes[i].p;
    let r0 = if dim == 2 { pn[2] } else { pn[6] };
    let r1 = if dim == 2 { po[2] } else { po[6] };
    let cb = bounds(c);
    let tol = 1.0e-9 * r1;
    guard abs(po[0] - pn[0]) < tol && abs(po[1] - pn[1]) < tol && r0 < r1 else { return false; };
    guard abs(cb[0] - po[0]) < tol && abs(cb[1] - po[1]) < tol && cb[3] >= po[0] + r1 - tol && cb[4] >= po[1] + r1 - tol else {
        return false;
    };
    var z0 = 0.0;
    var z1 = 0.0;
    if dim == 3 {
        let bo = bounds(o);
        let bi = bounds(i);
        z0 = max(bo[2], cb[2]);
        z1 = min(bo[5], cb[5]);
        guard z1 > z0 && bi[2] <= z0 + tol && bi[5] >= z1 - tol else { return false; };
    }
    let nr = max(1, ceil((r1 - r0) / size - 1.0e-9) as i64);
    let nt = max(1, ceil(0.5 * PI * r1 / size - 1.0e-9) as i64);
    let nz = if dim == 3 { max(1, ceil((z1 - z0) / size - 1.0e-9) as i64) } else { 0 };
    clear_mesh(dim);
    for k in nz + 1 {
        for j in nt + 1 {
            for m in nr + 1 {
                let r = r0 + (r1 - r0) * m / nr;
                let t = 0.5 * PI * j / nt;
                add_point(po[0] + r * cos(t), po[1] + r * sin(t), if dim == 3 { z0 + (z1 - z0) * k / nz } else { 0.0 });
            }
        }
    }
    fn node(m: i64, j: i64, k: i64) -> i64 { return m + (nr + 1) * (j + (nt + 1) * k); }
    for k in max(nz, 1) {
        for j in nt {
            for m in nr {
                var nodes: i64[..8] = [node(m, j, k), node(m + 1, j, k), node(m + 1, j + 1, k), node(m, j + 1, k)];
                if dim == 3 {
                    for n in [node(m, j, k + 1), node(m + 1, j, k + 1), node(m + 1, j + 1, k + 1), node(m, j + 1, k + 1)] {
                        nodes.push(n);
                    }
                }
                add_grid_cell(nodes[..]);
            }
        }
    }
    return true;
}

fn own_section(g: i64, z0: f64, z1: f64) -> i64 {
    let s = shapes[g];
    let p = s.p;
    let b = bounds(g);
    let tol = 1.0e-9 * max(1.0, abs(z1 - z0));
    guard b[2] <= z0 + tol && b[5] >= z1 - tol else { return -1; };
    if s.kind == BOX { return rect(p[0], p[1], p[3], p[4]); }
    if s.kind == CYLINDER { return if round_shape(g, 3) { circle(p[0], p[1], p[6]) } else { -1 }; }
    if s.kind == EXTRUDE { return s.a; }
    return if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT { -2 } else { -1 };
}

fn section(g: i64, z0: f64, z1: f64) -> i64 {
    var todo: i64[>..] = [g];
    var done: i64[>..] = [];
    while todo.len > 0 {
        let n = todo[todo.len - 1];
        todo.pop();
        if n < 0 {
            let c = done[done.len - 1];
            done.pop();
            let a = done[done.len - 1];
            done.pop();
            done.push(if a >= 0 && c >= 0 { add_shape(shapes[-1 - n].kind, a, c, default<f64[8]>()) } else { -1 });
            continue;
        }
        let own = own_section(n, z0, z1);
        if own != -2 {
            done.push(own);
            continue;
        }
        todo.push(-1 - n);
        todo.push(shapes[n].b);
        todo.push(shapes[n].a);
    }
    return done[0];
}

fn mesh_sweep(s: i64, z0: f64, z1: f64, size: f64) -> bool {
    let done = mesh_blossom(s, 2.0 * size);
    all_quads(s);
    var xy: f64[>..] = [];
    for v in point_count() {
        let p = point(v);
        xy.push(p[0]);
        xy.push(p[1]);
    }
    var quads: i64[>..] = [];
    for c in cell_count() {
        for i in 4 { quads.push(cell_node(c, i)); }
    }
    let n = xy.len / 2;
    let nz = max(1, ceil((z1 - z0) / size - 1.0e-9) as i64);
    clear_mesh(3);
    for k in nz + 1 {
        for v in n { add_point(xy[2 * v], xy[2 * v + 1], z0 + (z1 - z0) * k / nz); }
    }
    for k in nz {
        for q in quads.len / 4 {
            var nodes: i64[..8] = [];
            for i in 4 { nodes.push(quads[4 * q + i] + n * k); }
            for i in 4 { nodes.push(quads[4 * q + i] + n * (k + 1)); }
            add_grid_cell(nodes[..]);
        }
    }
    return done;
}

fn has_round(g: i64) -> bool {
    var todo: i64[>..] = [g];
    while todo.len > 0 {
        let s = shapes[todo[todo.len - 1]];
        todo.pop();
        if s.kind == CYLINDER || s.kind == EXTRUDE { return true; }
        if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT {
            todo.push(s.b);
            todo.push(s.a);
        }
    }
    return false;
}

fn targets_span_z(z0: f64, z1: f64) -> bool {
    for t in targets {
        if t.axis >= 0 { return false; }
        let tb = bounds(t.g);
        if tb[2] > z0 || tb[5] < z1 { return false; }
    }
    return true;
}

fn sweep_tets(s: i64, z0: f64, z1: f64, size: f64) -> bool {
    section_z = 0.5 * (z0 + z1);
    let done = mesh_2d(s, size);
    section_z = 0.0;
    var xy: f64[>..] = [];
    for v in point_count() {
        let p = point(v);
        xy.push(p[0]);
        xy.push(p[1]);
    }
    var tris: i64[>..] = [];
    for c in cell_count() {
        for i in 3 { tris.push(cell_node(c, i)); }
    }
    let n = xy.len / 2;
    let nz = max(1, ceil((z1 - z0) / size - 1.0e-9) as i64);
    clear_mesh(3);
    for k in nz + 1 {
        for v in n { add_point(xy[2 * v], xy[2 * v + 1], z0 + (z1 - z0) * k / nz); }
    }
    for k in nz {
        for t in tris.len / 3 {
            let a = n * k;
            let b = n * (k + 1);
            add_prism(tris[3 * t] + a, tris[3 * t + 1] + a, tris[3 * t + 2] + a, tris[3 * t] + b, tris[3 * t + 1] + b, tris[3 * t + 2] + b);
        }
    }
    return done;
}

fn turnable() -> bool {
    for s in shapes {
        if s.kind != BOX && s.kind != CYLINDER && s.kind != UNION && s.kind != INTERSECT && s.kind != SUBTRACT { return false; }
    }
    return true;
}

fn turn_shapes(count: i64, shift: i64) {
    for i in count {
        guard shapes[i].kind == BOX || shapes[i].kind == CYLINDER else { continue; }
        let p = shapes[i].p;
        var q: f64[8] = p;
        for c in 3 {
            q[c] = p[(c + shift) % 3];
            q[3 + c] = p[3 + (c + shift) % 3];
        }
        shapes[i].p = q;
    }
}

fn turn_nodes(shift: i64) {
    for n in point_count() {
        let x = point(n);
        for c in 3 { runtime::coords[3 * n + c] = x[(c + shift) % 3]; }
    }
}

fn shell_parts(g: i64) -> f64[3], f64, f64 {
    let o, i, c = sector_parts(g);
    guard o >= 0 && shapes[o].kind == SPHERE && shapes[i].kind == SPHERE && shapes[c].kind == BOX else { return [0.0; 3], 0.0, 0.0; };
    let po = shapes[o].p;
    let pn = shapes[i].p;
    let cb = bounds(c);
    let tol = 1.0e-9 * po[3];
    for a in 3 {
        guard abs(po[a] - pn[a]) < tol && abs(cb[a] - po[a]) < tol && cb[3 + a] >= po[a] + po[3] - tol else { return [0.0; 3], 0.0, 0.0; };
    }
    guard pn[3] < po[3] else { return [0.0; 3], 0.0, 0.0; };
    return [po[0], po[1], po[2]], pn[3], po[3];
}

fn shell_direction(k: i64, s: f64, t: f64) -> f64[3] {
    var e: f64[3] = [0.0; 3];
    e[k] = 1.0;
    var a = e;
    a[(k + 1) % 3] = 1.0;
    var b = e;
    b[(k + 2) % 3] = 1.0;
    let c: f64[3] = [1.0, 1.0, 1.0];
    let q = e * ((1.0 - s) * (1.0 - t)) + a * (s * (1.0 - t)) + c * (s * t) + b * ((1.0 - s) * t);
    return q * (1.0 / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]));
}

fn shell_layers(centre: f64[3], r0: f64, r1: f64, size: f64, radii: f64[>..]&, cells: i64&) -> bool {
    var n = 1;
    radii.clear();
    radii.push(r0);
    n = max(1, ceil(0.25 * PI * r1 / size - 1.0e-9) as! i64);
    if targets.len == 0 {
        let nr = max(1, ceil((r1 - r0) / size - 1.0e-9) as! i64);
        radii.clear();
        for l in nr + 1 { radii.push(r0 + (r1 - r0) * l / nr); }
    } else {
        let inv = 1.0 / sqrt(3.0);
        let dirs: f64[12] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, inv, inv, inv];
        fn wanted(r: f64) -> f64 {
            var h = size;
            for k in 4 { h = min(h, wanted_size(shell_centre[0] + r * dirs[3 * k], shell_centre[1] + r * dirs[3 * k + 1], shell_centre[2] + r * dirs[3 * k + 2], 0.0, size)); }
            return h;
        }
        shell_centre = centre;
        n = 1;
        var r = r0;
        while r < r1 - 1.0e-9 * r1 && radii.len < SHELL_LAYERS {
            n = max(n, ceil(0.25 * PI * r / wanted(r) - 1.0e-9) as! i64);
            r += max(wanted(r + 0.5 * wanted(r)), 1.0e-6 * r1);
            radii.push(min(r, r1));
        }
        guard radii[radii.len - 1] >= r1 - 1.0e-9 * r1 else { return false; }
        n = max(n, ceil(0.25 * PI * r1 / wanted(r1) - 1.0e-9) as! i64);
        if radii.len >= 3 && radii[radii.len - 1] - radii[radii.len - 2] < 0.5 * (radii[radii.len - 2] - radii[radii.len - 3]) {
            radii.pop();
            radii[radii.len - 1] = r1;
        }
        guard 3 * n * n * (radii.len - 1) * 6 <= SHELL_CELLS else { return false; }
    }
    cells = n;
    return true;
}

var shell_centre: f64[3] = [0.0; 3];
let SHELL_LAYERS = 400;
let SHELL_CELLS = 4000000;

fn shell_tets(g: i64, size: f64) -> bool {
    let centre, r0, r1 = shell_parts(g);
    guard r1 > 0.0 else { return false; };
    var radii: f64[>..] = [r0];
    var n = 1;
    guard shell_layers(centre, r0, r1, size, radii, n) else { return false; };
    let nr = radii.len - 1;
    var dirs: f64[>..] = [];
    var at: i64[>..] = [];
    for k in 3 {
        for j in n + 1 {
            for i in n + 1 {
                let d = shell_direction(k, i * 1.0 / n, j * 1.0 / n);
                var found = -1;
                for v in dirs.len / 3 {
                    let w: f64[3] = [dirs[3 * v], dirs[3 * v + 1], dirs[3 * v + 2]];
                    if abs(w[0] - d[0]) + abs(w[1] - d[1]) + abs(w[2] - d[2]) < 1.0e-12 { found = v; }
                }
                if found < 0 {
                    found = dirs.len / 3;
                    dirs.append(d);
                }
                at.push(found);
            }
        }
    }
    let m = dirs.len / 3;
    clear_mesh(3);
    for l in nr + 1 {
        let r = radii[l];
        for v in m { add_point(centre[0] + r * dirs[3 * v], centre[1] + r * dirs[3 * v + 1], centre[2] + r * dirs[3 * v + 2]); }
    }
    fn node(k: i64, i: i64, j: i64) -> i64 { return at[k * (n + 1) * (n + 1) + j * (n + 1) + i]; }
    for l in nr {
        for k in 3 {
            for j in n {
                for i in n {
                    let q: i64[4] = [node(k, i, j), node(k, i + 1, j), node(k, i + 1, j + 1), node(k, i, j + 1)];
                    add_prism(q[0] + m * l, q[1] + m * l, q[2] + m * l, q[0] + m * (l + 1), q[1] + m * (l + 1), q[2] + m * (l + 1));
                    add_prism(q[0] + m * l, q[2] + m * l, q[3] + m * l, q[0] + m * (l + 1), q[2] + m * (l + 1), q[3] + m * (l + 1));
                }
            }
        }
    }
    return true;
}

fn mesh_mapped(g: i64, size: f64, planes: i64[:]) {
    var boxes: f64[>..] = [];
    for s in planes {
        var extra: f64[>..] = [];
        if box_union(s, extra) { for p in extra { boxes.push(p); } }
    }
    fill_lattice(g, 3, size, boxes[..], []);
    let n = point_count();
    let faces: i64[4][6] = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
    var count = dictionary<i64, i64> {};
    fn key(c: i64, f: i64) -> i64 {
        var v: i64[..4] = [];
        for i in 4 { v.push(cell_node(c, faces[f][i])); }
        sort(v);
        return (v[0] * n + v[1]) * n + v[2];
    }
    for c in cell_count() {
        for f in 6 { count.update(key(c, f), 0) { it += 1; }; }
    }
    var outer: bool[>..] = [];
    for v in n { outer.push(false); }
    for c in cell_count() {
        for f in 6 {
            let k = count.get(key(c, f));
            guard k && k == 1 else { continue; };
            for i in 4 { outer[cell_node(c, faces[f][i])] = true; }
        }
    }
    for v in n {
        guard outer[v] else { continue; };
        let x = project(g, point(v));
        for c in 3 { runtime::coords[3 * v + c] = x[c]; }
    }
}
