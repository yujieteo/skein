// Tests of the mesher. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
//
// Each case meshes a shape, prints its mesh_quality, and compares the area
// or the volume of the mesh with the exact value of the shape. The gate
// case is the plate with a hole of femahoy validation/cases/grid-plate-
// hole-ccx (examples/femahoy/grid_plate_hole_ccx.skn), in mm.
import std;
import math;
import dictionary;
import os;
import .mesher;
import .lint;
import .stl;
import syntax from "../syntax";   // the lint imports it
import check from "../check";     // the lint imports it
import runtime from "../runtime";   // the mesher imports it
import solvers from "../solvers";   // the runtime imports it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

// The mesh is conforming: each face (an edge in 2D) is in one cell or in
// two, and each face of one cell has its nodes on the surface of g.
fn expect_conforming(what: u8[:], g: i64) {
    let k = mesher::mesh_dim + 1;
    let n = mesher::point_count();
    var faces = dictionary<i64, i64> {};
    for c in mesher::cell_count() {
        for skip in k {
            var f: i64[..3] = [];
            for i in k {
                if i != skip { f.push(mesher::cell_node(c, i)); }
            }
            sort(f);
            let key = if k == 3 { f[0] * n + f[1] } else { (f[0] * n + f[1]) * n + f[2] };
            faces.update(key, 0) { it += 1; };
        }
    }
    let b = mesher::bounds(g);
    let tol = 1.0e-6 * sqrt((b[3] - b[0]) * (b[3] - b[0]) + (b[4] - b[1]) * (b[4] - b[1]) + (b[5] - b[2]) * (b[5] - b[2]));
    var bad = 0;
    var outer = 0;
    faces.each() { key, count =>
        if count > 2 { bad++; }
        if count == 1 {
            outer++;
            var rest = key;
            for i in k - 1 {
                let v = rest % n;
                rest /= n;
                let p = mesher::vertex(v);
                if abs(mesher::phi(g, p[0], p[1], p[2])) > tol { bad++; }
            }
        }
    };
    if bad > 0 { fail(str(what, ": ", bad, " faces are not conforming")); }
}

// Meshes g, prints the quality and the measure against the exact value,
// and checks both. Returns the relative error of the measure.
fn mesh_case(what: u8[:], g: i64, size: f64, exact: f64, tolerance: f64) -> f64 {
    let ok = if mesher::bounds(g)[5] > mesher::bounds(g)[2] {
        mesher::mesh_3d(g, size)
    } else {
        mesher::mesh_2d(g, size)
    };
    if !ok { fail(str(what, ": the mesher did not finish")); }
    let q = mesher::mesh_quality();
    let err = (q.measure - exact) / exact;
    print(what, ", size ", size, "\n", q, "\n  exact: ", exact, "\n  error: ", err, " (tolerance ", tolerance, ")");
    if !mesher::quality_passes(q) { fail(str(what, ": mesh_quality")); }
    if abs(err) > tolerance { fail(str(what, ": the measure is ", q.measure, ", not ", exact, " within ", tolerance)); }
    expect_conforming(what, g);
    return err;
}

// The area of the 2D plate with a hole: the plate less the polygon of the
// boundary edges on the circle. Each such edge is a chord that cuts off
// r^2 / 2 (theta - sin theta) of the hole. Thus this is exact for a mesh
// that fills the domain of its boundary edges, with no gap or overlap.
fn chord_area(w: f64, h: f64, cx: f64, cy: f64, r: f64) -> f64 {
    var area = w * h - PI * r * r;
    var edges = dictionary<i64, i64> {};
    let n = mesher::point_count();
    for c in mesher::cell_count() {
        for i in 3 {
            let a = mesher::cell_node(c, i);
            let b = mesher::cell_node(c, (i + 1) % 3);
            edges.update(min(a, b) * n + max(a, b), 0) { it += 1; };
        }
    }
    fn on_circle(v: i64) -> bool {
        let p = mesher::vertex(v);
        return abs(sqrt((p[0] - cx) * (p[0] - cx) + (p[1] - cy) * (p[1] - cy)) - r) <= 1.0e-9 * r;
    }
    edges.each() { key, count =>
        let a = key / n;
        let b = key % n;
        if count == 1 && on_circle(a) && on_circle(b) {
            let theta = 2.0 * asin(0.5 * mesher::norm3(mesher::sub3(mesher::vertex(a), mesher::vertex(b))) / r);
            area += 0.5 * r * r * (theta - sin(theta));
        }
    };
    return area;
}

// The shortest edge of the cells of a 2D mesh.
fn shortest_edge() -> f64 {
    var shortest = 1.0e300;
    for c in mesher::cell_count() {
        for i in 3 {
            let a = mesher::vertex(mesher::cell_node(c, i));
            let b = mesher::vertex(mesher::cell_node(c, (i + 1) % 3));
            shortest = min(shortest, mesher::norm3(mesher::sub3(a, b)));
        }
    }
    return shortest;
}

// The tolerances of the 3D cases. The boundary nodes are on the surface.
// A facet with edges of length h on a surface of radius R misses it by at
// most h^2 / (8 R), so a curved surface of area A loses at most
// A h^2 / (8 R). A chamfered sharp edge of length L loses at most a
// triangle with sides h, so L h^2 / 2. Both scale as h^2.
fn tolerance(curved_area: f64, radius: f64, edge_length: f64, size: f64, volume: f64) -> f64 {
    var loss = edge_length * size * size / 2.0;
    if radius > 0.0 { loss += curved_area * size * size / (8.0 * radius); }
    return loss / volume;
}

// The ids of the diagnostics of `skein lint --mesh` on `src`, one per line.
fn mesh_lint_ids(src: u8[:]) -> u8[>..] {
    check::lint_text("t.skn", src);
    mesher::mesh_lint();
    var out: u8[>..] = [];
    for d in check::diags {
        out.append(d.id);
        out.push('\n');
    }
    return out;
}

fn expect_mesh_lints(what: u8[:], src: u8[:], want: u8[:]) {
    let got = mesh_lint_ids(src);
    if got != want { fail(str("lint --mesh ", what, ": got [", got, "]")); }
}

// `mesher::lint(words)` exits 0 and gives exactly the ids `want`, one per line.
fn expect_lint_command(what: u8[:], words: const u8[:][:], want: u8[:]) {
    var copies: u8[..64][>..] = [];
    for w in words {
        var c: u8[..64] = [];
        c.append(w);
        copies.push(c);
    }
    var args: u8[:][>..] = [];
    for &w in copies { args.push(w[..]); }
    let status = mesher::lint(args[..]);
    var got: u8[>..] = [];
    for d in check::diags {
        got.append(d.id);
        got.push('\n');
    }
    if status != 0 || got != want { fail(str("skein lint ", what, ": the exit status is ", status, ", the ids [", got, "]")); }
}

// The area inside fill node g, by the Green theorem: the integral of x dy
// on each cubic span and on the closing line.
fn fill_area(g: i64) -> f64 {
    let n = mesher::shapes[g].b;
    var area = 0.0;
    for i in n - 1 {
        let c = mesher::span(g, i);
        for a in 4 {
            for b in 1..4 { area += c[a] * b * c[4 + b] / (a + b); }
        }
    }
    let p = mesher::fill_point(g, n - 1);
    let q = mesher::fill_point(g, 0);
    return area + 0.5 * (p[0] + q[0]) * (q[1] - p[1]);
}

// fill (spec §8.1): a half disc of radius 10, a spline through 13 points
// of its arc and the line of its diameter, and a plate of 4 by 4 less a
// quarter ellipse with the semi-axes sqrt(2) and 1/sqrt(2), a spline
// through the 61 points of shape_hole_biaxial.skn at power 2.
fn test_fill() {
    var arc: f64[>..] = [];
    for k in 13 {
        let a = PI * (k / 12.0 - 0.5);
        arc.push(10.0 * cos(a));
        arc.push(10.0 * sin(a));
    }
    let half = mesher::fill(arc[..]);
    // The spline passes through each point, is inside the disc between
    // the points, and the line closes it.
    for k in 13 {
        if abs(mesher::phi(half, arc[2 * k], arc[2 * k + 1], 0.0)) > 1.0e-9 { fail(str("fill: point ", k, " is not on the boundary")); }
    }
    if !(mesher::phi(half, 5.0, 0.0, 0.0) < -4.9 && mesher::phi(half, -1.0, 0.0, 0.0) > 0.99 && mesher::phi(half, 11.0, 0.0, 0.0) > 0.99) {
        fail("fill: the sign of phi of the half disc");
    }
    // The end spans differ most from the arc: 5e-3 of the radius. Each
    // boundary edge of length h is a chord, and cuts off h^3 / (12 r) of
    // the arc of length h: h^2 / (12 r) in all on each unit of length. The
    // tolerance is 2 times that on the arc of length pi r.
    let span_area = fill_area(half);
    if abs(span_area / (0.5 * PI * 100.0) - 1.0) > 2.0e-3 { fail(str("fill: the area of the spline is ", span_area)); }
    mesh_case("fill of a half disc", half, 0.5, span_area, 2.0 * PI * 10.0 * 0.25 / 120.0 / span_area);
    var quarter: f64[>..] = [];
    for k in 61 {
        let t = 0.5 * PI * k / 60.0;
        quarter.push(sqrt(2.0) * cos(t));
        quarter.push(sin(t) / sqrt(2.0));
    }
    let hole = mesher::union(mesher::fill(quarter[..]), mesher::polygon([0.0, 0.0, sqrt(2.0), 0.0, 0.0, 1.0 / sqrt(2.0)]));
    mesh_case("plate less a quarter ellipse", mesher::subtract(mesher::rect(0.0, 0.0, 4.0, 4.0), hole), 0.05, 16.0 - 0.25 * PI,
              1.0e-4);
}

// The lint rules that need the mesh (docs/checks.md §3) on a strip of
// 1 m by 50 mm, in plane stress, with a load on its tip.
fn test_lint() {
    fn strip(hole: u8[:], order: u8[:]) -> u8[>..] {
        return str("model t: 2d, SI\nscales 1 m, 1 MPa\nlet plate: Surface = rect((0, 0) m, (1, 0.05) m)", hole,
                   "\nlet root: Curve = boundary(plate, x == 0 m)\nlet tip: Curve = boundary(plate, x == 1 m)\n",
                   "material steel {\n    E: Pressure = 200 GPa\n    nu: Real = 0.3\n}\n",
                   "mesh m on plate {\n    size: 0.02 m\n    cell: tri\n}\n",
                   "physics s on plate: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plane_stress(1 mm)\n}\n",
                   "method s: cg { mesh: m, order: ", order, " }\n",
                   "step go: static {\n    bc root: u = 0\n    bc tip: sigma*n = (0, -1) MPa\n    output u\n}\n");
    }
    expect_mesh_lints("bending_lock", strip("", "1"), "bending_lock\n");
    expect_mesh_lints("order 2 does not lock", strip("", "2"), "");
    expect_mesh_lints("mesh_feature_small", strip(" - circle((0.5, 0.025) m, 5 mm)", "2"), "mesh_feature_small\n");
    // The reader fills a spline, and the last point equal to the first
    // closes it. It has no fill of another Curve.
    expect_mesh_lints("fill of a spline", strip(" - fill(spline((0.4, 0.015) m, (0.5, 0.01) m, (0.6, 0.025) m, (0.5, 0.04) m, (0.4, 0.015) m))", "2"), "");
    expect_mesh_lints("fill of a line", strip(" - fill(line((0.5, 0.01) m, (0.52, 0.04) m))", "2"), "verify_unsupported\n");
    expect_mesh_lints("fill of 2 points", strip(" - fill(spline((0.5, 0.01) m, (0.52, 0.04) m, (0.5, 0.01) m))", "2"), "verify_unsupported\n");
    // The copy of a pattern is a feature. The reader needs a whole copy
    // count of 1 or more, and has no pattern_polar of a 2D shape.
    expect_mesh_lints("mesh_feature_small of a pattern", strip(" - pattern(circle((0.2, 0.025) m, 5 mm), (0.2, 0) m, 3)", "2"), "mesh_feature_small\n");
    expect_mesh_lints("a copy count of 2.5", strip(" - pattern(circle((0.2, 0.025) m, 5 mm), (0.2, 0) m, 2.5)", "2"), "verify_unsupported\n");
    expect_mesh_lints("a 2D polar pattern", strip(" - pattern_polar(circle((0.5, 0.025) m, 5 mm), (0, 1), 2)", "2"), "verify_unsupported\n");
    // In a 3d model: a loft from a circle to a square, with a polar
    // pattern of pins.
    fn solid(body: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m\nlet s: Solid = ", body, "\nmesh m on s {\n    size: 0.5 m\n    cell: tet\n}\n");
    }
    expect_mesh_lints("a loft and a polar pattern",
                      solid("loft(circle((0, 0) m, 4 m), rect((-2, -2) m, (2, 2) m), 3 m) | pattern_polar(cylinder((4, 0, 0) m, (0, 0, 1), 1 m, 3 m), (0, 0, 1), 6)"), "");
    // cell: hex makes the grid on a union of boxes, with a plane on each
    // face of the boxes of the model, and needs a union of boxes. Only its
    // planes use the layer, so it is an unused_let.
    fn grid_model(region: u8[:], cell: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m\nlet s: Solid = ", region, "\nlet layer: Solid = box((0, 0, 0) m, (1, 1, 0.3) m)\n",
                   "mesh m on s {\n    size: 0.25 m\n    cell: ", cell, "\n}\n");
    }
    expect_mesh_lints("cell: hex on a box", grid_model("box((0, 0, 0) m, (1, 1, 1) m)", "hex"), "unused_let\n");
    if mesher::cell_count() != 80 { fail(str("lint --mesh cell: hex on a box: ", mesher::cell_count(), " cells, not 80")); }
    expect_mesh_lints("cell: hex on a cylinder", grid_model("cylinder((0, 0, 0) m, (0, 0, 1), 1 m, 1 m)", "hex"), "unused_let\nverify_unsupported\n");
    expect_mesh_lints("cell: quad in 3D", grid_model("box((0, 0, 0) m, (1, 1, 1) m)", "quad"), "unused_let\nverify_unsupported\n");
    // bending_lock is for triangles and tetrahedra only: an order-1 hex
    // grid on a thin bar does not lock.
    fn bar(cell: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m, 1 MPa\nlet bar: Solid = box((0, 0, 0) m, (1, 0.05, 0.05) m)\n",
                   "let root: Surface = boundary(bar, x == 0 m)\nlet tip: Surface = boundary(bar, x == 1 m)\n",
                   "material steel {\n    E: Pressure = 200 GPa\n    nu: Real = 0.3\n}\n",
                   "mesh m on bar {\n    size: 0.02 m\n    cell: ", cell, "\n}\n",
                   "physics s on bar: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n",
                   "method s: cg { mesh: m, order: 1 }\n",
                   "step go: static {\n    bc root: u = 0\n    bc tip: sigma*n = (0, -1, 0) MPa\n    output u\n}\n");
    }
    expect_mesh_lints("bending_lock on a thin tet bar", bar("tet"), "bending_lock\n");
    expect_mesh_lints("no bending_lock on a thin hex grid", bar("hex"), "");
    // `skein lint` applies the rules that need the mesh only with --mesh.
    let path = "build/mesh_lint_strip.skn";
    guard write_file(path, strip("", "1")) else {
        fail(str("cannot write ", path));
        return;
    }
    expect_lint_command("with no --mesh", ["skein", "lint", path], "");
    expect_lint_command("with --mesh", ["skein", "lint", "--mesh", path], "bending_lock\n");
}

// Meshes g with a grid, and checks the counts of cells and nodes, the
// measure and the quality.
fn grid_case(what: u8[:], g: i64, dim: i64, size: f64, planes: i64[:], cells: i64, nodes: i64, exact: f64) {
    if !mesher::mesh_grid(g, dim, size, planes) { fail(str(what, ": not a union of boxes")); }
    let q = mesher::mesh_quality();
    print(what, ", size ", size, "\n", q);
    if q.cells != cells || q.nodes != nodes { fail(str(what, ": ", q.cells, " cells and ", q.nodes, " nodes, not ", cells, " and ", nodes)); }
    if abs(q.measure - exact) > 1.0e-10 * exact { fail(str(what, ": the measure is ", q.measure, ", not ", exact)); }
    if !mesher::quality_passes(q) { fail(str(what, ": mesh_quality")); }
}

// The grid of quad and hex cells on a union of boxes (spec §9.2).
fn test_grid() {
    // The heated block of femahoy grid-heated-block-ccx, in mm: a plate
    // 100 x 100 x 10 and a block 30 x 30 x 20 on it. Along x and y, the
    // planes 35 and 65 give 4, 3 and 4 cells. Along z, 1 cell in the plate
    // and 2 in the block. The plate has 11 x 11 cells and 12 x 12 x 2
    // nodes, the block 3 x 3 x 2 cells and 4 x 4 x 2 more nodes.
    mesher::clear_shapes();
    let plate = mesher::box(0.0, 0.0, 0.0, 100.0, 100.0, 10.0);
    let top = mesher::box(35.0, 35.0, 10.0, 65.0, 65.0, 30.0);
    grid_case("grid on a plate and a block", mesher::union(plate, top), 3, 10.0, [], 139, 320, 118000.0);
    // A layer of 0.3 of a unit cube: the face of the layer is a plane, so
    // z has 2 cells below it and 3 above it.
    let cube = mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let layer = mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 0.3);
    grid_case("grid on a cube with a layer", cube, 3, 0.25, [layer], 80, 150, 1.0);
    let q = mesher::mesh_quality();
    if abs(q.max_aspect - 0.25 / 0.15) > 1.0e-9 { fail(str("grid on a cube with a layer: the aspect is ", q.max_aspect)); }
    // The 2D L of 2 rects: 2 cells along each axis at size 9, 3 quads in
    // the L, and 8 nodes.
    let l_rects = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::rect(0.0, 0.0, 2.0, 10.0));
    grid_case("2D grid on an L of 2 rects", l_rects, 2, 9.0, [], 3, 8, 36.0);
    if mesher::mesh_grid(mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0), 3, 0.5, []) {
        fail("grid on a cylinder: it is not a union of boxes");
    }
}

// `skein export` of the model `src` with the options `words`: the exit
// status is `want`. The result is the STL.
fn export_case(what: u8[:], src: u8[:], words: const u8[:][:], want: i64) -> u8[>..] {
    let path = "build/export_case.skn";
    var out: u8[>..] = [];
    guard write_file(path, src) else {
        fail(str("cannot write ", path));
        return out;
    }
    var copies: u8[..64][>..] = [];
    for w in ["skein", "export", path] { copies.push(w); }
    for w in words { copies.push(w); }
    var args: u8[:][>..] = [];
    for &w in copies { args.push(w[..]); }
    let status = mesher::export_stl(args[..], out);
    if status != want { fail(str("export ", what, ": the exit status is ", status, ", not ", want)); }
    return out;
}

// Reads the facets of `stl`. Each facet has its normal of unit length and
// in the direction of the right-hand rule of its vertices. The result is
// the facet count, the volume by the divergence theorem, and the sum of
// the area vectors, which is 0 on a closed surface.
fn read_stl(what: u8[:], stl: u8[:]) -> i64, f64, f64[3] {
    var count = 0;
    var volume = 0.0;
    var area: f64[3] = [0.0, 0.0, 0.0];
    var normal: f64[3] = [0.0, 0.0, 0.0];
    var p: f64[3][..3] = [];
    for line in split(stl, '\n') {
        let words = split(trim(line), ' ');
        guard words.len >= 4 && (words[0] == "vertex" || words[0] == "normal" || words[1] == "normal") else { continue; };
        let at = if words[0] == "vertex" { 1 } else { 2 };
        var v: f64[3] = [0.0, 0.0, 0.0];
        for j in 3 { v[j] = parse_flt(words[at + j]); }
        if words[0] != "vertex" {
            normal = v;
            continue;
        }
        p.push(v);
        guard p.len == 3 else { continue; };
        let a = mesher::cross3(mesher::sub3(p[1], p[0]), mesher::sub3(p[2], p[0]));
        if abs(mesher::norm3(normal) - 1.0) > 1.0e-9 || mesher::dot3(a, normal) <= 0.0 {
            fail(str("export ", what, ": facet ", count, " has the normal ", normal[0], " ", normal[1], " ", normal[2]));
        }
        volume += mesher::dot3(p[0], mesher::cross3(p[1], p[2])) / 6.0;
        for j in 3 { area[j] += 0.5 * a[j]; }
        count++;
        p.clear();
    }
    return count, volume, area;
}

// `skein export` (spec §22): the surface of the tetrahedra of a Solid,
// closed and with its normals out, and the triangles of a 2D Surface.
fn test_export() {
    fn solid(region: u8[:], mesh: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m\nlet L: Length = 2 m\nlet s: Solid = ", region,
                   "\nlet top: Surface = boundary(s, z == 1 m)\n", mesh);
    }
    let holed = solid("box((0, 0, 0) m, (2, 1, 1) m) - cylinder((1, 0.5, -1) m, (0, 0, 1), 0.3 m, 3 m)",
                      "mesh m on s {\n    size: 0.25 m\n}\n");
    let fine = export_case("a block with a hole", holed, ["--geometry", "s", "--format", "stl"], 0);
    let n_fine, volume, area = read_stl("a block with a hole", fine);
    let tets = mesher::mesh_quality().measure;
    print("export a block with a hole: ", n_fine, " facets, volume ", volume, ", tetrahedra ", tets);
    if !starts_with(fine, "solid s\n") || !ends_with(fine, "endsolid s\n") { fail("export a block with a hole: no solid and endsolid lines"); }
    if n_fine == 0 || abs(volume - tets) > 1.0e-9 * tets { fail(str("export a block with a hole: the volume ", volume, " is not ", tets)); }
    if mesher::norm3(area) > 1.0e-9 { fail(str("export a block with a hole: the surface is open, area sum ", mesher::norm3(area))); }
    // --size sets the triangles in place of the size of the mesh block.
    let coarse = export_case("--size 50 cm", holed, ["--geometry", "s", "--size", "50 cm"], 0);
    let n_coarse, coarse_volume, coarse_area = read_stl("--size 50 cm", coarse);
    if n_coarse == 0 || n_coarse >= n_fine { fail(str("export --size 50 cm: ", n_coarse, " facets, not fewer than ", n_fine)); }
    if mesher::norm3(coarse_area) > 1.0e-9 { fail("export --size 50 cm: the surface is open"); }
    let unmeshed = solid("box((0, 0, 0) m, (1, 1, 1) m)", "");
    let n_cube, cube_volume, cube_area = read_stl("a cube", export_case("a cube", unmeshed, ["--geometry", "s", "--size", "500mm"], 0));
    if n_cube == 0 || abs(cube_volume - 1.0) > 1.0e-12 { fail(str("export a cube: the volume is ", cube_volume, ", not 1")); }
    // A 2D Surface: its triangles in the plane z = 0, with the normal +z.
    let plate = "model t: 2d, SI\nscales 1 m\nlet p: Surface = rect((0, 0) m, (2, 1) m)\n";
    let n_plate, plate_volume, plate_area = read_stl("a 2D rect", export_case("a 2D rect", plate, ["--geometry", "p", "--size", "0.25 m"], 0));
    if n_plate == 0 || abs(plate_area[2] - 2.0) > 1.0e-12 { fail(str("export a 2D rect: the area is ", plate_area[2], ", not 2")); }
    // Faults of the command line: exit status 2 and no STL.
    export_case("an unknown let", holed, ["--geometry", "nope", "--size", "1 m"], 2);
    export_case("--format step", holed, ["--geometry", "s", "--format", "step"], 2);
    export_case("a size with no unit", holed, ["--geometry", "s", "--size", "3"], 2);
    export_case("no --geometry", holed, ["--size", "1 m"], 2);
    export_case("an unknown option", holed, ["--geometry", "s", "--sizes", "1 m"], 2);
    export_case("with no mesh block and no --size", unmeshed, ["--geometry", "s"], 2);
    // Faults of the model: exit status 1. A let that is not a shape, a
    // Surface of a Solid that the mesher does not read, and a check error.
    export_case("a Length", holed, ["--geometry", "L"], 1);
    export_case("a boundary", holed, ["--geometry", "top"], 1);
    export_case("a model with an error", str(holed, "let bad: Length = 2 s\n"), ["--geometry", "s"], 1);
}

// The interfaces of mesh_2d: a ground under a trapezoid lift, with the
// line y = 0 between them inside the body. No triangle crosses the line or
// the slopes of the lift, and the areas of the two parts are exact.
fn test_interfaces() {
    mesher::clear_shapes();
    let ground = mesher::rect(-6.0, -3.0, 6.0, 0.0);
    let lift = mesher::polygon([-4.0, 0.0, 4.0, 0.0, 3.0, 1.0, -3.0, 1.0]);
    let body = mesher::union(ground, lift);
    mesher::clear_interfaces();
    mesher::add_interface(lift);
    let done = mesher::mesh_2d(body, 0.7);
    mesher::clear_interfaces();
    if !done { fail("interfaces: the mesher did not finish"); }
    var area_lift = 0.0;
    var area = 0.0;
    var crossed = 0;
    for c in mesher::cell_count() {
        var inside = 0;
        var outside = 0;
        var centre: f64[3] = [0.0; 3];
        let a = mesher::vertex(mesher::cell_node(c, 0));
        let b = mesher::vertex(mesher::cell_node(c, 1));
        let d = mesher::vertex(mesher::cell_node(c, 2));
        let size = 0.5 * abs((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]));
        for i in 3 {
            let p = mesher::vertex(mesher::cell_node(c, i));
            let f = mesher::phi(lift, p[0], p[1], 0.0);
            if f < -1.0e-9 { inside++; }
            if f > 1.0e-9 { outside++; }
            for t in 3 { centre[t] += p[t] / 3.0; }
        }
        if inside > 0 && outside > 0 { crossed++; }
        if mesher::phi(lift, centre[0], centre[1], 0.0) < 0.0 { area_lift += size; }
        area += size;
    }
    print("interfaces: ", mesher::cell_count(), " triangles, ", crossed, " cross the lift, area of the lift ", area_lift);
    if crossed > 0 { fail(str("interfaces: ", crossed, " triangles cross the boundary of the lift")); }
    if abs(area_lift - 7.0) > 1.0e-9 || abs(area - 43.0) > 1.0e-9 { fail(str("interfaces: the areas ", area_lift, " and ", area, " are not 7 and 43")); }
}

fn main() {
    // The 2D plate: 100 x 60 mm with a hole of radius 10 mm at (30, 30).
    // The tolerance is the chord defect of the hole, (2 pi / 32)^2 / 6 of
    // its area for 32 boundary edges of 2 mm, over the plate area.
    mesher::clear_shapes();
    let plate2 = mesher::subtract(mesher::rect(0.0, 0.0, 100.0, 60.0), mesher::circle(30.0, 30.0, 10.0));
    let exact2 = 6000.0 - PI * 100.0;
    let chord_tol = (TAU / 32.0) * (TAU / 32.0) / 6.0 * PI * 100.0 / exact2;
    mesh_case("2D plate with a hole", plate2, 2.0, exact2, chord_tol);
    let polygon_area = chord_area(100.0, 60.0, 30.0, 30.0, 10.0);
    let q2 = mesher::mesh_quality();
    print("  area of the plate less the chord polygon of the hole: ", polygon_area,
          ", mesh error ", (q2.measure - polygon_area) / polygon_area);
    if abs(q2.measure - polygon_area) > 1.0e-10 * polygon_area {
        fail(str("2D plate: the mesh area ", q2.measure, " is not the polygon area ", polygon_area));
    }

    // A hole of radius 3 mm at size 5 mm: 4 chords, each with a sliver
    // between it and the arc that is deeper than the centroid of a small
    // triangle on the chord. The mesh must fill the chord polygon.
    let coarse = mesher::subtract(mesher::rect(0.0, 0.0, 20.0, 20.0), mesher::circle(10.0, 10.0, 3.0));
    if !mesher::mesh_2d(coarse, 5.0) { fail("2D coarse hole: the mesher did not finish"); }
    expect_conforming("2D coarse hole", coarse);
    let coarse_area = chord_area(20.0, 20.0, 10.0, 10.0, 3.0);
    let qc = mesher::mesh_quality();
    print("2D coarse hole, size 5\n", qc, "\n  area of the chord polygon: ", coarse_area);
    if abs(qc.measure - coarse_area) > 1.0e-10 * coarse_area {
        fail(str("2D coarse hole: the mesh area ", qc.measure, " is not the polygon area ", coarse_area));
    }

    // Edges of 2 primitives on one line: the L shape of a union, the
    // rectangle of an intersect, and an L of a rect and a polygon. At
    // size 3 the 2 edges on each line are cut at different points.
    let l_rects = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::rect(0.0, 0.0, 2.0, 10.0));
    mesh_case("2D L of 2 rects", l_rects, 3.0, 36.0, 1.0e-10);
    let overlap = mesher::intersect(mesher::rect(0.0, 0.0, 10.0, 10.0), mesher::rect(0.0, 0.0, 4.0, 6.0));
    mesh_case("2D intersect of 2 rects", overlap, 3.0, 24.0, 1.0e-10);
    let arm: f64[8] = [0.0, 0.0, 2.0, 0.0, 2.0, 10.0, 0.0, 10.0];
    let l_polygon = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::polygon(arm));
    mesh_case("2D L of a rect and a polygon", l_polygon, 3.0, 36.0, 1.0e-10);

    // Corners below MIN_ANGLE: wedges of 21.8, 16.7 and 20 degrees. At
    // the larger sizes, the 2 segments at the corner have different
    // lengths, and a midpoint split of each encroaches the other.
    let wide: f64[6] = [0.0, 0.0, 10.0, 0.0, 10.0, 4.0];
    let narrow: f64[6] = [0.0, 0.0, 10.0, 0.0, 10.0, 3.0];
    let twenty: f64[6] = [0.0, 0.0, 10.0, 0.0, 6.578, 2.394];
    for size in [1.0, 2.5] {
        mesh_case("2D wedge of 21.8 deg", mesher::polygon(wide), size, 20.0, 1.0e-10);
        mesh_case("2D wedge of 16.7 deg", mesher::polygon(narrow), size, 15.0, 1.0e-10);
    }
    mesh_case("2D wedge of 20 deg", mesher::polygon(twenty), 2.0, 0.5 * 10.0 * 2.394, 1.0e-10);

    // A rect with a round end, as for a slot: the circle touches the top
    // and the bottom edges of the rect at their ends.
    let slot = mesher::union(mesher::rect(0.0, 0.0, 20.0, 10.0), mesher::circle(20.0, 5.0, 5.0));
    for size in [0.7, 2.0, 4.5] {
        if !mesher::mesh_2d(slot, size) { fail(str("2D slot, size ", size, ": the mesher did not finish")); }
        expect_conforming("2D slot", slot);
        let q = mesher::mesh_quality();
        print("2D slot, size ", size, "\n", q);
        if !mesher::quality_passes(q) { fail(str("2D slot, size ", size, ": mesh_quality")); }
        if shortest_edge() < 0.1 * size { fail(str("2D slot, size ", size, ": the shortest edge is ", shortest_edge())); }
    }

    // A circle that cuts an edge of the plate at an angle that is not 90
    // deg: each of the 2 corners must be one point, so no edge is much
    // shorter than the size.
    for cy in [2.0, 5.0, 8.0] {
        let notch = mesher::subtract(mesher::rect(0.0, 0.0, 100.0, 60.0), mesher::circle(50.0, cy, 10.0));
        let what = str("2D notch at y = ", cy);
        if !mesher::mesh_2d(notch, 5.0) { fail(str(what, ": the mesher did not finish")); }
        expect_conforming(what, notch);
        if !mesher::quality_passes(mesher::mesh_quality()) { fail(str(what, ": mesh_quality")); }
        let shortest = shortest_edge();
        print(what, ", size 5: shortest edge ", shortest);
        if shortest < 0.1 * 5.0 { fail(str(what, ": the shortest edge is ", shortest)); }
    }

    // The femahoy plate with a hole and a pad, in mm. The sharp edges: the
    // 2 rims of the hole, and the edges of the pad, plate and pad faces on
    // the lattice planes, at most their perimeters.
    let bore = mesher::cylinder(30.0, 30.0, 0.0, 0.0, 0.0, 1.0, 10.0, 10.0);
    let plate = mesher::subtract(mesher::box(0.0, 0.0, 0.0, 100.0, 60.0, 10.0), bore);
    let body = mesher::union(plate, mesher::box(60.0, 20.0, 10.0, 90.0, 40.0, 20.0));
    let body_volume = 60000.0 - PI * 1000.0 + 6000.0;
    mesh_case("3D plate with a hole and a pad", body, 2.0, body_volume,
              tolerance(TAU * 100.0, 10.0, 2.0 * TAU * 10.0, 2.0, body_volume));

    // The cylinder at 2 sizes: the error must fall as h^2, so the order of
    // convergence must be near 2.
    let cyl = mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 10.0, 30.0);
    let cyl_volume = PI * 3000.0;
    let e2 = mesh_case("cylinder", cyl, 2.0, cyl_volume, tolerance(TAU * 300.0, 10.0, 2.0 * TAU * 10.0, 2.0, cyl_volume));
    let e1 = mesh_case("cylinder", cyl, 1.0, cyl_volume, tolerance(TAU * 300.0, 10.0, 2.0 * TAU * 10.0, 1.0, cyl_volume));
    let order = log2(e2 / e1);
    print("cylinder: order of convergence of the volume ", order);
    if order < 1.5 || order > 2.5 { fail(str("cylinder: the order of convergence is ", order, ", not 2")); }

    // A hexagon of radius 10 mm, extruded 20 mm: 6 sharp side edges of 20 mm
    // and 2 hexagon rims of 60 mm.
    let hexagon: f64[12] = [10.0, 0.0, 5.0, 8.660254037844386, -5.0, 8.660254037844386,
                            -10.0, 0.0, -5.0, -8.660254037844386, 5.0, -8.660254037844386];
    let prism = mesher::extrude(mesher::polygon(hexagon), 20.0);
    let prism_volume = 1.5 * sqrt(3.0) * 100.0 * 20.0;
    mesh_case("extruded hexagon", prism, 2.0, prism_volume, tolerance(0.0, 0.0, 240.0, 2.0, prism_volume));

    // A tube of radii 6 and 10 mm and length 30 mm, the revolve of a
    // rectangle: a full turn, and a quarter turn.
    let section = mesher::rect(6.0, 0.0, 10.0, 30.0);
    let tube = mesher::revolve(section, TAU);
    let tube_volume = PI * (100.0 - 36.0) * 30.0;
    mesh_case("revolved tube", tube, 1.0, tube_volume,
              tolerance(TAU * 16.0 * 30.0, 6.0, 4.0 * TAU * 8.0, 1.0, tube_volume));
    let quarter = mesher::revolve(section, 0.5 * PI);
    mesh_case("quarter tube", quarter, 1.0, 0.25 * tube_volume,
              tolerance(0.25 * TAU * 16.0 * 30.0, 6.0, 4.0 * 0.25 * TAU * 8.0 + 4.0 * 30.0 + 4.0 * 4.0, 1.0, 0.25 * tube_volume));

    // A loft from a circle of radius 10 mm to one of 5 mm, 20 mm up: a
    // frustum of a cone. A loft between 2 squares is a frustum of a
    // pyramid, as φ of a square centred on the origin is max(|x|, |y|) - a.
    let slant = sqrt(25.0 + 400.0);
    let cone = mesher::loft(mesher::circle(0.0, 0.0, 10.0), mesher::circle(0.0, 0.0, 5.0), 20.0);
    let cone_volume = PI * 20.0 / 3.0 * (100.0 + 50.0 + 25.0);
    mesh_case("lofted circles", cone, 1.0, cone_volume,
              tolerance(PI * 15.0 * slant, 5.0, TAU * 15.0, 1.0, cone_volume));
    let pyramid = mesher::loft(mesher::rect(-5.0, -5.0, 5.0, 5.0), mesher::rect(-2.0, -2.0, 2.0, 2.0), 20.0);
    let pyramid_volume = 20.0 / 3.0 * (100.0 + 16.0 + 40.0);
    mesh_case("lofted squares", pyramid, 1.0, pyramid_volume,
              tolerance(0.0, 0.0, 56.0 + 4.0 * sqrt(400.0 + 2.0 * 9.0), 1.0, pyramid_volume));

    // Patterns: 3 cubes in a row, 6 cylinders about the z axis, and 4
    // boxes about the y axis, each turned by 90 deg onto the next.
    let row = mesher::pattern(mesher::box(0.0, 0.0, 0.0, 4.0, 4.0, 4.0), 6.0, 0.0, 0.0, 3);
    mesh_case("pattern of 3 cubes", row, 1.0, 192.0, tolerance(0.0, 0.0, 144.0, 1.0, 192.0));
    let pins = mesher::pattern_polar(mesher::cylinder(8.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 5.0), 0.0, 0.0, 1.0, 6);
    let pins_volume = 6.0 * PI * 4.0 * 5.0;
    mesh_case("polar pattern of 6 cylinders", pins, 0.5, pins_volume,
              tolerance(6.0 * TAU * 2.0 * 5.0, 2.0, 6.0 * 2.0 * TAU * 2.0, 0.5, pins_volume));
    let cross = mesher::pattern_polar(mesher::box(3.0, -1.0, -1.0, 5.0, 1.0, 1.0), 0.0, 1.0, 0.0, 4);
    mesh_case("polar pattern of 4 boxes about y", cross, 0.5, 32.0, tolerance(0.0, 0.0, 4.0 * 24.0, 0.5, 32.0));

    // 2D patterns: 4 squares in a row, and a plate less a row of 4 holes.
    mesh_case("2D pattern of 4 squares", mesher::pattern(mesher::rect(0.0, 0.0, 2.0, 2.0), 3.0, 0.0, 0.0, 4), 1.0, 16.0, 1.0e-10);
    let holes = mesher::subtract(mesher::rect(0.0, 0.0, 40.0, 10.0), mesher::pattern(mesher::circle(5.0, 5.0, 2.0), 10.0, 0.0, 0.0, 4));
    let holes_area = 400.0 - 4.0 * PI * 4.0;
    mesh_case("2D plate with a row of 4 holes", holes, 0.5, holes_area,
              (TAU / 25.0) * (TAU / 25.0) / 6.0 * 4.0 * PI * 4.0 / holes_area);

    test_fill();
    test_interfaces();

    test_grid();
    test_lint();
    test_export();

    if failures > 0 { exit(1); }
    print("mesher: all tests pass");
}
