// The simplex mesh: P1 or P2 triangles in 2D, P1 or P2 tetrahedra in 3D.
//
// The mesh is unstructured: a list of node coordinates and a list of cells,
// each with `per_cell` local nodes in the sequence of shape.goose. The
// structured meshers below map a lattice of cells onto a box or onto a
// plate with a quarter hole, and cut each lattice cell into simplices.
// build_simplices takes the vertices and cells of another mesher, for
// example the tri3 and tet4 mesh of src/mesher/.
import std;
import os;
import math;
import dictionary;
namespace elements;
import .shape;
import .csg;

var dim = 0;
var order = 1;
var per_cell = 0;
var coords: f64[>..] = [];      // x, y and z of each node (z = 0 in 2D)
var conn: i64[>..] = [];        // the nodes of each cell
var bbox_lo: f64[3] = [0.0, 0.0, 0.0];
var bbox_hi: f64[3] = [0.0, 0.0, 0.0];
var cell_body: i64[>..] = [];   // the body of each cell of build_bodies (csg.goose), else empty
var node_body: i64[>..] = [];   // the body of each node, likewise
var body_shapes: i64[>..] = []; // the shape of each body

fn node_count() -> i64 { return coords.len / 3; }
fn cell_count() -> i64 { return if per_cell > 0 { conn.len / per_cell } else { 0 }; }

fn node_x(n: i64, t: i64) -> f64 {
    guard n >= 0 && n < node_count() && t >= 0 && t < 3 else { abort("node_x: n < node_count() and t < 3"); }
    return coords[3 * n + t];
}

// The tolerance of an equality on a coordinate: 1e-6 of the size of the
// bounding box (spec §8.2).
fn geometry_tolerance() -> f64 {
    var size = 0.0;
    for t in 3 { size = max(size, bbox_hi[t] - bbox_lo[t]); }
    return 1.0e-6 * size;
}

// --- the lattice -------------------------------------------------------------

// The lattice of nodes has order * n + 1 points along an axis of n cells,
// so a P2 midpoint node is a lattice point too. Node (i, j, k) has the
// index i + lat[0] (j + lat[1] k).
var lat: i64[3] = [1, 1, 1];
var on_lattice = false;     // node n is lattice point n: a mesh of build_lattice

fn lattice_node(i: i64, j: i64, k: i64) -> i64 {
    guard i >= 0 && i < lat[0] && j >= 0 && j < lat[1] && k >= 0 && k < lat[2] else { abort("lattice_node: (i, j, k) in the lattice"); }
    return i + lat[0] * (j + lat[1] * k);
}

// The kind of map from the lattice to the geometry.
let MAP_BOX = 0;
let MAP_HOLE = 1;

var map_kind = MAP_BOX;
var map_lo: f64[3] = [0.0, 0.0, 0.0];
var map_hi: f64[3] = [0.0, 0.0, 0.0];
var hole_r = 0.0;
var radial: f64[>..] = [];      // the parameter 0 to 1 of each radial vertex row

fn start_mesh(d: i64, ord: i64, lo: f64[3], hi: f64[3]) {
    guard (d == 2 || d == 3) && (ord == 1 || ord == 2) else { abort("start_mesh: d is 2 or 3 and ord 1 or 2"); }
    dim = d;
    order = ord;
    hex_cells = false;
    per_cell = node_count_of(d, ord);
    coords.clear();
    conn.clear();
    on_lattice = false;
    cell_body.clear();
    node_body.clear();
    body_shapes.clear();
    map_lo = lo;
    map_hi = hi;
    bbox_lo = lo;
    bbox_hi = hi;
    body_shape = -1;
}

// The position of lattice point (i, j, k), which may be a midpoint.
fn map_point(i: i64, j: i64, k: i64) -> f64[3] {
    guard i >= 0 && i < lat[0] && j >= 0 && j < lat[1] && k >= 0 && k < lat[2] && lat[0] > 1 && lat[1] > 1 && (map_kind == MAP_BOX || j / order < radial.len) else { abort("map_point: (i, j, k) in the lattice and a radial row for j"); }
    let s = i as f64 / (lat[0] - 1);
    let w = if lat[2] > 1 { k as f64 / (lat[2] - 1) } else { 0.0 };
    let z = map_lo[2] + w * (map_hi[2] - map_lo[2]);
    if map_kind == MAP_BOX {
        let t = j as f64 / (lat[1] - 1);
        return [map_lo[0] + s * (map_hi[0] - map_lo[0]), map_lo[1] + t * (map_hi[1] - map_lo[1]), z];
    }
    // The hole: s runs along the arc from the x axis to the y axis, and
    // along the outer edges x = hi and y = hi; j runs out from the arc.
    let phi = 0.5 * PI * s;
    let ax = map_lo[0] + hole_r * cos(phi);
    let ay = map_lo[1] + hole_r * sin(phi);
    let ox = if s <= 0.5 { map_hi[0] } else { map_hi[0] - (2.0 * s - 1.0) * (map_hi[0] - map_lo[0]) };
    let oy = if s <= 0.5 { map_lo[1] + 2.0 * s * (map_hi[1] - map_lo[1]) } else { map_hi[1] };
    let t = radial[j / order];
    return [ax + t * (ox - ax), ay + t * (oy - ay), z];
}

// Puts a node at each lattice point. A midpoint node is the midpoint of its
// edge, so each cell edge is straight, except on the arc of a hole, where
// the node is on the arc.
fn place_nodes() {
    guard lat[0] > 1 && lat[1] > 1 && lat[2] >= 1 && coords.len == 0 else { abort("place_nodes: a lattice of cells and no node yet"); }
    for k in lat[2] {
        for j in lat[1] {
            for i in lat[0] {
                let odd = i % order != 0 || j % order != 0 || k % order != 0;
                if !odd || (map_kind == MAP_HOLE && j == 0) {
                    let p = map_point(i, j, k);
                    coords.append(p);
                } else {
                    // The edges of the cut cells go up along each axis (see
                    // cut_cells), so the ends are -1 and +1 on each odd index.
                    let oi = i % 2;
                    let oj = j % 2;
                    let ok = k % 2;
                    let a = map_point(i - oi, j - oj, k - ok);
                    let b = map_point(i + oi, j + oj, k + ok);
                    for t in 3 { coords.push(0.5 * (a[t] + b[t])); }
                }
            }
        }
    }
}

// The signed size of a simplex from its vertex nodes.
fn signed_size(v: i64[4]) -> f64 {
    guard dim == 2 || dim == 3 else { abort("signed_size: dim is 2 or 3"); }
    var e: f64[9] = [0.0; 9];
    for r in dim {
        for t in 3 { e[3 * r + t] = node_x(v[r + 1], t) - node_x(v[0], t); }
    }
    if dim == 2 { return e[0] * e[4] - e[1] * e[3]; }
    return e[0] * (e[4] * e[8] - e[5] * e[7]) - e[1] * (e[3] * e[8] - e[5] * e[6])
        + e[2] * (e[3] * e[7] - e[4] * e[6]);
}

// Adds the simplex with the vertices at lattice offsets `off` (3 per
// vertex, in cell units) from lattice cell corner (ci, cj, ck).
fn add_simplex(ci: i64, cj: i64, ck: i64, off: i64[12]) {
    guard ci >= 0 && cj >= 0 && ck >= 0 && order * (ci + 1) < lat[0] && order * (cj + 1) < lat[1] && (dim == 2 || order * (ck + 1) < lat[2]) else { abort("add_simplex: cell (ci, cj, ck) in the lattice"); }
    var p: i64[12] = [0; 12];
    var v: i64[4] = [0; 4];
    for a in dim + 1 {
        p[3 * a] = order * (ci + off[3 * a]);
        p[3 * a + 1] = order * (cj + off[3 * a + 1]);
        p[3 * a + 2] = order * (ck + off[3 * a + 2]);
        v[a] = lattice_node(p[3 * a], p[3 * a + 1], p[3 * a + 2]);
    }
    if signed_size(v) < 0.0 {
        for t in 3 { swap(p[3 + t], p[6 + t]); }
        swap(v[1], v[2]);
    }
    for a in dim + 1 { conn.push(v[a]); }
    guard order == 2 else { return; }
    for e in edge_count(dim) {
        let a = edge_node(dim, e, 0);
        let b = edge_node(dim, e, 1);
        conn.push(lattice_node((p[3 * a] + p[3 * b]) / 2, (p[3 * a + 1] + p[3 * b + 1]) / 2,
                               (p[3 * a + 2] + p[3 * b + 2]) / 2));
    }
}

// Cuts each lattice cell into 2 triangles or 6 tetrahedra (Kuhn). Each
// tetrahedron goes from corner (0, 0, 0) to corner (1, 1, 1) by one step
// along each axis, so the diagonals of two cells agree on their shared face
// and the mesh conforms.
let KUHN: i64[3][6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

fn cut_cells(nx: i64, ny: i64, nz: i64) {
    guard nx > 0 && ny > 0 && nz > 0 && lat[0] == order * nx + 1 && lat[1] == order * ny + 1 else { abort("cut_cells: nx, ny and nz > 0 and the lattice of these cells"); }
    for ck in nz {
        for cj in ny {
            for ci in nx {
                if dim == 2 {
                    add_simplex(ci, cj, ck, [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0]);
                    add_simplex(ci, cj, ck, [0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0]);
                    continue;
                }
                for perm in KUHN {
                    var off: i64[12] = [0; 12];
                    for a in 1..4 {
                        for t in 3 { off[3 * a + t] = off[3 * (a - 1) + t]; }
                        off[3 * a + perm[a - 1]] = 1;
                    }
                    add_simplex(ci, cj, ck, off);
                }
            }
        }
    }
}

fn build_lattice(nx: i64, ny: i64, nz: i64) {
    guard nx > 0 && ny > 0 && nz > 0 && per_cell > 0 else { abort("build_lattice: nx, ny and nz > 0 after start_mesh"); }
    lat = [order * nx + 1, order * ny + 1, if dim == 3 { order * nz + 1 } else { 1 }];
    place_nodes();
    cut_cells(nx, ny, if dim == 3 { nz } else { 1 });
    on_lattice = true;
}

// --- the meshers -------------------------------------------------------------

// A box (3D) or a rectangle (2D, z ignored) with nx × ny × nz lattice cells.
fn build_box(d: i64, ord: i64, lo: f64[3], hi: f64[3], nx: i64, ny: i64, nz: i64) {
    guard nx > 0 && ny > 0 && nz > 0 && hi[0] > lo[0] && hi[1] > lo[1] && (d == 2 || hi[2] > lo[2]) else { abort("build_box: nx, ny and nz > 0 and hi > lo"); }
    start_mesh(d, ord, lo, hi);
    map_kind = MAP_BOX;
    build_lattice(nx, ny, nz);
}

// A box with nx × ny × nz hexahedra of order 1 (8 nodes, shape.goose): the
// lattice of the nodes is that of the simplex mesh, and each lattice cell is one cell.
fn build_hex_box(lo: f64[3], hi: f64[3], nx: i64, ny: i64, nz: i64) {
    guard nx > 0 && ny > 0 && nz > 0 && hi[0] > lo[0] && hi[1] > lo[1] && hi[2] > lo[2] else { abort("build_hex_box: nx, ny and nz > 0 and hi > lo"); }
    start_mesh(3, 1, lo, hi);
    hex_cells = true;
    per_cell = node_count_of(3, 1);
    map_kind = MAP_BOX;
    lat = [nx + 1, ny + 1, nz + 1];
    place_nodes();
    for ck in nz {
        for cj in ny {
            for ci in nx {
                for l in 8 { conn.push(lattice_node(ci + HEX_CORNERS[3 * l], cj + HEX_CORNERS[3 * l + 1], ck + HEX_CORNERS[3 * l + 2])); }
            }
        }
    }
    on_lattice = true;
}

// The radial parameter of each vertex row: cells that grow from `fine` at
// the hole by a factor `GRADE` up to `size`, then stay at `size`, scaled
// to fill the length `span`.
let GRADE = 1.25;

fn grade_rows(span: f64, fine: f64, size: f64) {
    guard span > 0.0 && fine > 0.0 && size > 0.0 else { abort("grade_rows: span, fine and size > 0"); }
    var cells: f64[>..] = [];
    var h = min(fine, size);
    var total = 0.0;
    while total < span * (1.0 - 1.0e-9) {
        cells.push(h);
        total += h;
        h = min(h * GRADE, size);
    }
    radial.clear();
    radial.push(0.0);
    var sum = 0.0;
    for c in cells {
        sum += c;
        radial.push(sum / total);
    }
}

// A rectangle (2D) or a slab (3D) from lo to hi, less a hole of radius r
// whose axis is z through the corner (lo.x, lo.y): the quarter model of a
// plate with a round hole. The cells have the size `fine` at the hole and
// grow to `size`.
fn build_hole(d: i64, ord: i64, lo: f64[3], hi: f64[3], r: f64, size: f64, fine: f64) {
    guard r > 0.0 && r < min(hi[0] - lo[0], hi[1] - lo[1]) else { abort("build_hole: 0 < r < the sides of the plate"); }
    start_mesh(d, ord, lo, hi);
    map_kind = MAP_HOLE;
    hole_r = r;
    let w = hi[0] - lo[0];
    let h = hi[1] - lo[1];
    grade_rows(sqrt(w * w + h * h) - r, fine, size);
    // The cells along the arc are of size `fine` at the hole and of at
    // most `size` on the outer edges. An even count puts a node on the
    // corner (hi.x, hi.y).
    var around = max(ceil(0.5 * PI * r / fine) as i64, ceil((w + h) / size) as i64);
    around += around % 2;
    let layers = if d == 3 { max(ceil((hi[2] - lo[2]) / size - 1.0e-9) as i64, 1) } else { 1 };
    build_lattice(around, radial.len - 1, layers);
}

// The cells `vertices` (d + 1 nodes each, in positive orientation) on the
// nodes `xyz` (x, y and z of each node). For order 2 it adds a node at the
// midpoint of each edge, once for the cells that share the edge.
fn build_simplices(d: i64, ord: i64, xyz: f64[:], vertices: i64[:]) {
    guard xyz.len >= 3 && xyz.len % 3 == 0 && vertices.len % (d + 1) == 0 && (vertices.len == 0 || (min(vertices) >= 0 && max(vertices) < xyz.len / 3)) else { abort("build_simplices: xyz holds the nodes and vertices d + 1 nodes of xyz for each cell"); }
    var lo: f64[3] = [0.0, 0.0, 0.0];
    var hi: f64[3] = [0.0, 0.0, 0.0];
    for t in 3 {
        lo[t] = xyz[t];
        hi[t] = xyz[t];
        for n in xyz.len / 3 {
            lo[t] = min(lo[t], xyz[3 * n + t]);
            hi[t] = max(hi[t], xyz[3 * n + t]);
        }
    }
    start_mesh(d, ord, lo, hi);
    coords.append(xyz);
    let nv = xyz.len / 3;
    var mid = dictionary<i64, i64> {};
    for c in vertices.len / (d + 1) {
        let v = vertices[(d + 1) * c..(d + 1) * (c + 1)];
        conn.append(v);
        guard ord == 2 else { continue; }
        for e in edge_count(d) {
            let a = v[edge_node(d, e, 0)];
            let b = v[edge_node(d, e, 1)];
            let key = min(a, b) * nv + max(a, b);
            let m = mid.get_or_insert(key, node_count());
            if m == node_count() {
                for t in 3 { coords.push(0.5 * (xyz[3 * a + t] + xyz[3 * b + t])); }
            }
            conn.push(m);
        }
    }
}

// --- the shared mesh of a study (spec §20.4) ---------------------------------

fn save_array<T>(dir: const u8[:], name: u8[:], a: const T[:]) -> bool { return write_file(str(dir, "/", name), to_bytes(a)); }

fn load_array<T>(dir: const u8[:], name: u8[:], into: T[>..]&) -> bool {
    var b: u8[>..] = [];
    guard read_file(str(dir, "/", name), b) else { return false; };
    let a, ok = from_bytes<T[>..]>(b);
    into.clear();
    into.append(a);
    return ok;
}

// Writes the mesh in the files of the directory `dir`.
fn save_mesh(dir: const u8[:]) -> bool {
    let sizes: i64[>..] = [dim, order, per_cell, body_shape];
    var bbox: f64[>..] = [];
    bbox.append(bbox_lo);
    bbox.append(bbox_hi);
    return save_array(dir, "sizes", sizes) && save_array(dir, "bbox", bbox) && save_array(dir, "coords", coords) &&
           save_array(dir, "conn", conn) && save_array(dir, "cell_body", cell_body) && save_array(dir, "node_body", node_body) &&
           save_array(dir, "body_shapes", body_shapes);
}

// Reads the mesh that save_mesh wrote in the directory `dir`.
fn load_mesh(dir: const u8[:]) -> bool {
    var sizes: i64[>..] = [];
    var bbox: f64[>..] = [];
    guard load_array(dir, "sizes", sizes) && sizes.len == 4 && load_array(dir, "bbox", bbox) && bbox.len == 6 &&
          load_array(dir, "coords", coords) && load_array(dir, "conn", conn) && load_array(dir, "cell_body", cell_body) &&
          load_array(dir, "node_body", node_body) && load_array(dir, "body_shapes", body_shapes) else { return false; };
    on_lattice = false;
    dim = sizes[0];
    order = sizes[1];
    per_cell = sizes[2];
    body_shape = sizes[3];
    hex_cells = dim == 3 && per_cell == 8;
    for t in 3 {
        bbox_lo[t] = bbox[t];
        bbox_hi[t] = bbox[3 + t];
    }
    return true;
}
