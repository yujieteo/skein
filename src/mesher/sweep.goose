// The `quad` and `hex` cells of spec §9.2 on a region that is not a union
// of rects or boxes:
//
// - a quarter of a ring, (circle - circle) & rect, or of a tube,
//   (cylinder - cylinder) & box, with the corner of the rect or the box on
//   the axis: a grid in r, θ (and z);
// - in 2D, another shape: Blossom-Quad (quad.goose);
// - in 3D, a prism along z (each primitive is a box, a cylinder along z or
//   an extrude, and it holds the z range of the region): the quads of
//   Blossom-Quad on the section at twice the size, each cut into 4 quads
//   or a triangle into 3, swept through layers of equal height;
// - in 3D, another shape: the grid of the planes of the boxes in the
//   bounds of the region, the cells whose centre is inside, and each node
//   on the boundary of the cells moved onto the boundary of the region
//   (Schneiders 1996).
import std;
import math;
import dictionary;
namespace mesher;
import .geometry;
import .mesh;
import .grid;
import .quad;
import runtime from "../runtime";

// Meshes g with quad (2D) or hex (3D) cells of about `size`, with the
// planes of each shape of `planes`. False if mesh_2d does not finish.
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

// A circle in 2D, or a cylinder along +z in 3D.
fn round_shape(g: i64, dim: i64) -> bool {
    let s = shapes[g];
    if dim == 2 { return s.kind == CIRCLE; }
    return s.kind == CYLINDER && abs(s.p[3]) < 1.0e-12 && abs(s.p[4]) < 1.0e-12;
}

// The outer shape, the inner shape and the box (or rect) of a quarter of a
// ring or a tube, in either form (o - i) & b or (o & b) - i, or -1.
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

// Meshes a quarter of a ring or a tube. False if g is not one.
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

// The 2D section of node g alone, or -2 if it is a union, intersect or
// subtract: its section is the same operation on the sections of its parts.
fn own_section(g: i64, z0: f64, z1: f64) -> i64 {
    let s = shapes[g];
    let p = s.p;
    let b = bounds(g);
    let tol = 1.0e-9 * max(1.0, abs(z1 - z0));
    guard b[2] <= z0 + tol && b[5] >= z1 - tol else { return -1; };
    if s.kind == BOX { return rect(p[0], p[1], p[3], p[4]); }
    if s.kind == CYLINDER { return if round_shape(g, 3) { circle(p[0], p[1], p[6]) } else { -1 }; }
    if s.kind == EXTRUDE { return s.a; }
    // A shape that is subtracted may hold more than [z0, z1], not less.
    return if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT { -2 } else { -1 };
}

// The 2D section of g between z0 and z1, if g is a prism along z there:
// each box, cylinder along z and extrude of g holds [z0, z1]. Else -1. A
// stack walks the tree of g, with no recursion (rule 1).
fn section(g: i64, z0: f64, z1: f64) -> i64 {
    var todo: i64[>..] = [g];     // the nodes; -1 - n after the parts of node n
    var done: i64[>..] = [];      // the sections, in the sequence of the nodes
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

// Sweeps the quads of the 2D section s from z0 to z1. False as mesh_2d.
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

// Meshes g with the grid of the planes of its bounds and of the boxes of
// `planes`, and moves each node on the boundary of the cells onto the
// boundary of g.
fn mesh_mapped(g: i64, size: f64, planes: i64[:]) {
    var boxes: f64[>..] = [];
    for s in planes {
        var extra: f64[>..] = [];
        if box_union(s, extra) { for p in extra { boxes.push(p); } }
    }
    fill_lattice(g, 3, size, boxes[..], []);
    // The faces in one cell. A face is keyed by its 3 lowest nodes.
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
