// The structured grid of spec §9.2 on a union of axis-aligned boxes (3D)
// or rects (2D): hex8 or quad4 cells. Each face of each box of the region,
// and of each box of `planes` inside the bounds of the region, is a grid
// plane. Between two neighbour planes, the grid has the fewest equal cells
// that are not longer than the size. A cell of the lattice is in the mesh
// if its centre is in a box of the region. Thus the cells fill the region
// exactly, and no cell crosses a face of a box.
import std;
import math;
namespace mesher;
import .geometry;
import .mesh;

// Appends the lo and hi corners of each box of shape g to `out`, 6 numbers
// for each box (z = 0 for a rect). False if g is not a union of boxes or of
// rects.
recursive fn box_union(g: i64, out: f64[>..]&) -> bool {
    let s = shapes[g];
    if s.kind == UNION { return box_union(s.a, out) && box_union(s.b, out); }
    if s.kind == BOX {
        for i in 6 { out.push(s.p[i]); }
        return true;
    }
    if s.kind == RECT {
        for i in [0, 1] { out.push(s.p[i]); }
        out.push(0.0);
        for i in [2, 3] { out.push(s.p[i]); }
        out.push(0.0);
        return true;
    }
    return false;
}

// The grid lines along axis a, from lo to hi: the planes of `boxes` in
// (lo, hi), and the equal cells of at most `size` between two of them.
fn grid_lines(boxes: f64[:], a: i64, lo: f64, hi: f64, size: f64, out: f64[>..]&) {
    var cuts: f64[>..] = [lo, hi];
    for b in boxes.len / 6 {
        for c in [boxes[6 * b + a], boxes[6 * b + 3 + a]] {
            if c > lo && c < hi { cuts.push(c); }
        }
    }
    sort(cuts);
    let tol = 1.0e-9 * (hi - lo);
    out.clear();
    out.push(lo);
    var last = lo;
    for c in cuts {
        guard c - last > tol else { continue; };
        let n = max(1, ceil((c - last) / size - 1.0e-9) as i64);
        for k in 1..n + 1 { out.push(last + (c - last) * k / n); }
        last = c;
    }
}

// Meshes shape g with a grid of cells of at most `size`, with the planes
// of g and of each shape of `planes` that is a union of boxes. False if g
// is not a union of boxes or of rects.
fn mesh_grid(g: i64, dim: i64, size: f64, planes: i64[:]) -> bool {
    var region: f64[>..] = [];
    guard box_union(g, region) else { return false; };
    var boxes: f64[>..] = [];
    for p in region { boxes.push(p); }
    for s in planes {
        var extra: f64[>..] = [];
        if box_union(s, extra) { for p in extra { boxes.push(p); } }
    }
    fill_lattice(g, dim, size, boxes[..], region[..]);
    return true;
}

// The grid lines of the planes of `boxes` in the bounds of g, and a cell
// of the lattice where its centre is in a box of `region`, or where φ of g
// is negative at its centre if `region` is empty.
fn fill_lattice(g: i64, dim: i64, size: f64, boxes: f64[:], region: f64[:]) {
    let b = bounds(g);
    var lx: f64[>..] = [];
    var ly: f64[>..] = [];
    var lz: f64[>..] = [0.0];
    grid_lines(boxes, 0, b[0], b[3], size, lx);
    grid_lines(boxes, 1, b[1], b[4], size, ly);
    if dim == 3 { grid_lines(boxes, 2, b[2], b[5], size, lz); }
    let nx = lx.len;
    let ny = ly.len;
    let nz = lz.len;
    clear_mesh(dim);
    var ids: i64[>..] = [];
    for i in nx * ny * nz { ids.push(-1); }
    fn node(i: i64, j: i64, k: i64) -> i64 {
        let at = i + nx * (j + ny * k);
        if ids[at] < 0 { ids[at] = add_point(lx[i], ly[j], lz[k]); }
        return ids[at];
    }
    fn inside(x: f64[3]) -> bool {
        if region.len == 0 { return phi(g, x[0], x[1], x[2]) < 0.0; }
        for r in region.len / 6 {
            var all = true;
            for a in dim {
                if x[a] < region[6 * r + a] || x[a] > region[6 * r + 3 + a] { all = false; }
            }
            if all { return true; }
        }
        return false;
    }
    let layers = if dim == 3 { nz - 1 } else { 1 };
    for k in layers {
        for j in ny - 1 {
            for i in nx - 1 {
                let z = if dim == 3 { 0.5 * (lz[k] + lz[k + 1]) } else { 0.0 };
                let centre: f64[3] = [0.5 * (lx[i] + lx[i + 1]), 0.5 * (ly[j] + ly[j + 1]), z];
                guard inside(centre) else { continue; };
                // The VTK order: the face k counter-clockwise, then the face k + 1.
                var nodes: i64[..8] = [node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k)];
                if dim == 3 {
                    for n in [node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1), node(i, j + 1, k + 1)] {
                        nodes.push(n);
                    }
                }
                add_grid_cell(nodes[..]);
            }
        }
    }
}
