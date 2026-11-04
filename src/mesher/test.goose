import std;
import math;
import dictionary;
import os;
import .mesher;
import .lint;
import .stl;
import syntax from "../syntax";
import check from "../check";
import runtime from "../runtime";
import solvers from "../solvers";

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

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

fn tolerance(curved_area: f64, radius: f64, edge_length: f64, size: f64, volume: f64) -> f64 {
    var loss = edge_length * size * size / 2.0;
    if radius > 0.0 { loss += curved_area * size * size / (8.0 * radius); }
    return loss / volume;
}

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

fn test_fill() {
    var arc: f64[>..] = [];
    for k in 13 {
        let a = PI * (k / 12.0 - 0.5);
        arc.push(10.0 * cos(a));
        arc.push(10.0 * sin(a));
    }
    let half = mesher::fill(arc[..]);
    for k in 13 {
        if abs(mesher::phi(half, arc[2 * k], arc[2 * k + 1], 0.0)) > 1.0e-9 { fail(str("fill: point ", k, " is not on the boundary")); }
    }
    if !(mesher::phi(half, 5.0, 0.0, 0.0) < -4.9 && mesher::phi(half, -1.0, 0.0, 0.0) > 0.99 && mesher::phi(half, 11.0, 0.0, 0.0) > 0.99) {
        fail("fill: the sign of phi of the half disc");
    }
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

fn test_nurbs() {
    let r2 = sqrt(0.5);
    let ring: f64[27] = [10.0, 0.0, 1.0, 10.0, 10.0, r2, 0.0, 10.0, 1.0, -10.0, 10.0, r2, -10.0, 0.0, 1.0,
                         -10.0, -10.0, r2, 0.0, -10.0, 1.0, 10.0, -10.0, r2, 10.0, 0.0, 1.0];
    let disc = mesher::fill_nurbs(2, [0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0], ring[..]);
    if mesher::shapes[disc].b != 5 { fail(str("nurbs: the circle has ", mesher::shapes[disc].b - 1, " spans, not 4")); }
    var worst = 0.0;
    for i in mesher::shapes[disc].b - 1 {
        let c = mesher::span(disc, i);
        for k in 101 {
            let x, _, _ = mesher::cubic(c, 0, k / 100.0);
            let y, _, _ = mesher::cubic(c, 1, k / 100.0);
            worst = max(worst, abs(sqrt(x * x + y * y) - 10.0));
        }
    }
    if worst > 1.0e-12 { fail(str("nurbs: a point of the circle is ", worst, " from the radius")); }
    if abs(mesher::phi(disc, 0.0, 0.0, 0.0) + 10.0) > 1.0e-9 || abs(mesher::phi(disc, 12.0, 1.0, 0.0) - (sqrt(145.0) - 10.0)) > 1.0e-9 ||
       abs(mesher::phi(disc, 6.0, -7.0, 0.0) - (sqrt(85.0) - 10.0)) > 1.0e-9 {
        fail("nurbs: phi of the circle is not the distance to it");
    }
    if abs(nurbs_area(disc) / (PI * 100.0) - 1.0) > 1.0e-9 { fail(str("nurbs: the area of the circle is ", nurbs_area(disc))); }
    mesh_case("fill of a nurbs circle", disc, 0.5, PI * 100.0, 2.0 * PI * 10.0 * 0.25 / 120.0 / (PI * 100.0));
    let quarter = mesher::fill_nurbs(2, [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], [1.0, 0.0, 1.0, 1.0, 1.0, r2, 0.0, 1.0, 1.0]);
    if abs(nurbs_area(quarter) - (0.25 * PI - 0.5)) > 1.0e-9 { fail(str("nurbs: the area of the quarter segment is ", nurbs_area(quarter))); }
    if !(mesher::phi(quarter, 0.6, 0.6, 0.0) < 0.0 && mesher::phi(quarter, 0.4, 0.4, 0.0) > 0.0) { fail("nurbs: the sign of phi of the quarter segment"); }
    let corner = mesher::fill_nurbs(1, [0.0, 0.0, 1.0, 2.0, 2.0], [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0]);
    if mesher::shapes[corner].b != 3 { fail(str("nurbs: the polyline of degree 1 has ", mesher::shapes[corner].b - 1, " spans, not 2")); }
    if abs(nurbs_area(corner) - 0.5) > 1.0e-12 { fail(str("nurbs: the area of the triangle of degree 1 is ", nurbs_area(corner))); }
}

fn nurbs_area(g: i64) -> f64 {
    let n = mesher::shapes[g].b;
    var area = 0.0;
    for i in n - 1 {
        let c = mesher::span(g, i);
        for k in 201 {
            let x, _, _ = mesher::cubic(c, 0, k / 200.0);
            let _, dy, _ = mesher::cubic(c, 1, k / 200.0);
            let f = if k == 0 || k == 200 { 1.0 } else if k % 2 == 1 { 4.0 } else { 2.0 };
            area += f * x * dy / 600.0;
        }
    }
    let p = mesher::fill_point(g, n - 1);
    let q = mesher::fill_point(g, 0);
    return area + 0.5 * (p[0] + q[0]) * (q[1] - p[1]);
}

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
    expect_mesh_lints("fill of a spline", strip(" - fill(spline((0.4, 0.015) m, (0.5, 0.01) m, (0.6, 0.025) m, (0.5, 0.04) m, (0.4, 0.015) m))", "2"), "");
    expect_mesh_lints("fill of a line", strip(" - fill(line((0.5, 0.01) m, (0.52, 0.04) m))", "2"), "verify_unsupported\n");
    expect_mesh_lints("fill of a nurbs", strip(" - fill(nurbs(2, (0.51, 0.025) m, 1, (0.51, 0.035) m, 0.7071067811865476, (0.5, 0.035) m, 1, (0.49, 0.035) m, 0.7071067811865476, (0.49, 0.025) m, 1))", "2"), "mesh_feature_small\n");
    expect_mesh_lints("fill of a nurbs with knots that fall", strip(" - fill(nurbs(1, knots(0, 0, 2, 1, 1), (0.5, 0.01) m, 1, (0.52, 0.04) m, 1, (0.48, 0.04) m, 1))", "2"), "verify_unsupported\n");
    expect_mesh_lints("fill of a nurbs with an inner knot degree + 1 times", strip(" - fill(nurbs(2, knots(0, 0, 0, 1, 1, 1, 2, 2, 2), (0.51, 0.025) m, 1, (0.51, 0.035) m, 0.7071067811865476, (0.5, 0.035) m, 1, (0.5, 0.035) m, 1, (0.49, 0.035) m, 0.7071067811865476, (0.49, 0.025) m, 1))", "2"), "verify_unsupported\n");
    expect_mesh_lints("fill of 2 points", strip(" - fill(spline((0.5, 0.01) m, (0.52, 0.04) m, (0.5, 0.01) m))", "2"), "verify_unsupported\n");
    expect_mesh_lints("mesh_feature_small of a pattern", strip(" - pattern(circle((0.2, 0.025) m, 5 mm), (0.2, 0) m, 3)", "2"), "mesh_feature_small\n");
    expect_mesh_lints("a copy count of 2.5", strip(" - pattern(circle((0.2, 0.025) m, 5 mm), (0.2, 0) m, 2.5)", "2"), "verify_unsupported\n");
    expect_mesh_lints("a 2D polar pattern", strip(" - pattern_polar(circle((0.5, 0.025) m, 5 mm), (0, 1), 2)", "2"), "verify_unsupported\n");
    fn solid(body: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m\nlet s: Solid = ", body, "\nmesh m on s {\n    size: 0.5 m\n    cell: tet\n}\n");
    }
    expect_mesh_lints("a loft and a polar pattern",
                      solid("loft(circle((0, 0) m, 4 m), rect((-2, -2) m, (2, 2) m), 3 m) | pattern_polar(cylinder((4, 0, 0) m, (0, 0, 1), 1 m, 3 m), (0, 0, 1), 6)"), "");
    fn grid_model(region: u8[:], cell: u8[:]) -> u8[>..] {
        return str("model t: 3d, SI\nscales 1 m\nlet s: Solid = ", region, "\nlet layer: Solid = box((0, 0, 0) m, (1, 1, 0.3) m)\n",
                   "mesh m on s {\n    size: 0.25 m\n    cell: ", cell, "\n}\n");
    }
    expect_mesh_lints("cell: hex on a box", grid_model("box((0, 0, 0) m, (1, 1, 1) m)", "hex"), "unused_let\n");
    if mesher::cell_count() != 80 { fail(str("lint --mesh cell: hex on a box: ", mesher::cell_count(), " cells, not 80")); }
    expect_mesh_lints("cell: hex on a cylinder", grid_model("cylinder((0, 0, 0) m, (0, 0, 1), 1 m, 1 m)", "hex"), "unused_let\n");
    expect_mesh_lints("cell: quad in 3D", grid_model("box((0, 0, 0) m, (1, 1, 1) m)", "quad"), "unused_let\nverify_unsupported\n");
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
    let path = "build/mesh_lint_strip.skn";
    guard write_file(path, strip("", "1")) else {
        fail(str("cannot write ", path));
        return;
    }
    expect_lint_command("with no --mesh", ["skein", "lint", path], "");
    expect_lint_command("with --mesh", ["skein", "lint", "--mesh", path], "bending_lock\n");
}

fn grid_case(what: u8[:], g: i64, dim: i64, size: f64, planes: i64[:], cells: i64, nodes: i64, exact: f64) {
    if !mesher::mesh_grid(g, dim, size, planes) { fail(str(what, ": not a union of boxes")); }
    let q = mesher::mesh_quality();
    print(what, ", size ", size, "\n", q);
    if q.cells != cells || q.nodes != nodes { fail(str(what, ": ", q.cells, " cells and ", q.nodes, " nodes, not ", cells, " and ", nodes)); }
    if abs(q.measure - exact) > 1.0e-10 * exact { fail(str(what, ": the measure is ", q.measure, ", not ", exact)); }
    if !mesher::quality_passes(q) { fail(str(what, ": mesh_quality")); }
}

fn test_grid() {
    mesher::clear_shapes();
    let plate = mesher::box(0.0, 0.0, 0.0, 100.0, 100.0, 10.0);
    let top = mesher::box(35.0, 35.0, 10.0, 65.0, 65.0, 30.0);
    grid_case("grid on a plate and a block", mesher::union(plate, top), 3, 10.0, [], 139, 320, 118000.0);
    let cube = mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let layer = mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 0.3);
    grid_case("grid on a cube with a layer", cube, 3, 0.25, [layer], 80, 150, 1.0);
    let q = mesher::mesh_quality();
    if abs(q.max_aspect - 0.25 / 0.15) > 1.0e-9 { fail(str("grid on a cube with a layer: the aspect is ", q.max_aspect)); }
    let l_rects = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::rect(0.0, 0.0, 2.0, 10.0));
    grid_case("2D grid on an L of 2 rects", l_rects, 2, 9.0, [], 3, 8, 36.0);
    if mesher::mesh_grid(mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0), 3, 0.5, []) {
        fail("grid on a cylinder: it is not a union of boxes");
    }
}

fn expect_cells_conforming(what: u8[:], g: i64) {
    let n = mesher::point_count();
    let hex_faces: i64[4][6] = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
    var faces = dictionary<i64, i64> {};
    var nodes = dictionary<i64, i64> {};
    fn add(f: i64[..4]) {
        var v = f;
        sort(v);
        let key = if v.len == 2 { v[0] * n + v[1] } else { (v[0] * n + v[1]) * n + v[2] };
        faces.update(key, 0) { it += 1; };
        for x in f { nodes.insert(key * 4 + nodes_of(key), x); }
    }
    fn nodes_of(key: i64) -> i64 {
        var k = 0;
        while nodes.contains(key * 4 + k) { k++; }
        return k;
    }
    for c in mesher::cell_count() {
        let t = runtime::cell_type[c] as i64;
        if t == runtime::HEX8 {
            for f in 6 {
                var v: i64[..4] = [];
                for i in 4 { v.push(mesher::cell_node(c, hex_faces[f][i])); }
                add(v);
            }
        } else {
            let k = if t == runtime::QUAD4 { 4 } else { 3 };
            for i in k { add([mesher::cell_node(c, i), mesher::cell_node(c, (i + 1) % k)]); }
        }
    }
    let b = mesher::bounds(g);
    let tol = 1.0e-6 * sqrt((b[3] - b[0]) * (b[3] - b[0]) + (b[4] - b[1]) * (b[4] - b[1]) + (b[5] - b[2]) * (b[5] - b[2]));
    var bad = 0;
    faces.each() { key, count =>
        if count > 2 { bad++; }
        if count == 1 {
            for k in 4 {
                let v = nodes.get(key * 4 + k);
                guard v else { break; };
                let p = mesher::vertex(v);
                if abs(mesher::phi(g, p[0], p[1], p[2])) > tol { bad++; }
            }
        }
    };
    if bad > 0 { fail(str(what, ": ", bad, " faces are not conforming")); }
}

fn cells_case(what: u8[:], g: i64, dim: i64, size: f64, exact: f64, tolerance: f64) -> i64 {
    if !mesher::mesh_cells(g, dim, size, []) { fail(str(what, ": the mesher did not finish")); }
    let q = mesher::mesh_quality();
    var tris = 0;
    for c in mesher::cell_count() {
        if runtime::cell_type[c] as i64 == runtime::TRI3 { tris++; }
    }
    let err = (q.measure - exact) / exact;
    print(what, ", size ", size, "\n", q, "\n  triangles: ", tris, "\n  error: ", err, " (tolerance ", tolerance, ")");
    if !mesher::quality_passes(q) { fail(str(what, ": mesh_quality")); }
    if abs(err) > tolerance { fail(str(what, ": the measure is ", q.measure, ", not ", exact, " within ", tolerance)); }
    expect_cells_conforming(what, g);
    return tris;
}

fn test_cells() {
    mesher::clear_shapes();
    let ring = mesher::subtract(mesher::intersect(mesher::circle(0.0, 0.0, 2.0), mesher::rect(0.0, 0.0, 2.0, 2.0)),
                                mesher::circle(0.0, 0.0, 1.0));
    let ring_area = 16.0 * 3.0 * sin(PI / 32.0) / 2.0;
    cells_case("quarter ring", ring, 2, 0.2, ring_area, 1.0e-12);
    if mesher::cell_count() != 80 || mesher::point_count() != 102 { fail("quarter ring: not 80 cells and 102 nodes"); }
    let tube = mesher::intersect(mesher::subtract(mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0),
                                                  mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
                                 mesher::box(0.0, 0.0, 0.0, 2.0, 2.0, 1.0));
    cells_case("quarter tube", tube, 3, 0.2, ring_area, 1.0e-12);
    if mesher::cell_count() != 400 || mesher::point_count() != 612 { fail("quarter tube: not 400 cells and 612 nodes"); }
    let plate = mesher::subtract(mesher::rect(0.0, 0.0, 2.0, 1.0), mesher::circle(1.0, 0.5, 0.25));
    mesher::mesh_2d(plate, 0.05);
    let tri_area = mesher::mesh_quality().measure;
    let tri_count = mesher::cell_count();
    let tris = cells_case("Blossom-Quad on a plate with a hole", plate, 2, 0.05, tri_area, 1.0e-12);
    let quads = mesher::cell_count() - tris;
    if 2 * quads + tris != tri_count { fail(str("Blossom-Quad: ", quads, " quads and ", tris, " triangles from ", tri_count)); }
    if 10 * tris > mesher::cell_count() { fail(str("Blossom-Quad: ", tris, " of ", mesher::cell_count(), " cells are triangles")); }
    let rod = mesher::intersect(mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 3.0), mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 3.0));
    let rod_volume = 0.75 * PI;
    cells_case("swept hex grid on a quarter of a rod", rod, 3, 0.2, rod_volume, tolerance(0.5 * PI * 3.0, 1.0, 0.0, 0.2, rod_volume));
    for c in mesher::cell_count() {
        if runtime::cell_type[c] as i64 != runtime::HEX8 { fail("swept hex grid: a cell is not a hex"); }
    }
    let holed = mesher::subtract(mesher::box(0.0, 0.0, 0.0, 2.0, 1.0, 1.0), mesher::cylinder(0.0, 0.5, 0.5, 1.0, 0.0, 0.0, 0.25, 2.0));
    let block_volume = 2.0 - 2.0 * PI * 0.0625;
    cells_case("mapped hex grid on a block with a hole", holed, 3, 0.1, block_volume, tolerance(2.0 * TAU * 0.25, 0.25, 0.0, 0.1, block_volume));
    mesher::clear_mesh(2);
    let s0 = mesher::add_point(0.0, 0.0, 0.0);
    let s1 = mesher::add_point(1.0, 0.0, 0.0);
    let s2 = mesher::add_point(1.0, 1.0, 0.0);
    let s3 = mesher::add_point(0.0, 1.0, 0.0);
    mesher::add_grid_cell([s0, s1, s2, s3]);
    let r0 = mesher::add_point(2.0, 0.0, 0.0);
    let r1 = mesher::add_point(2.0 + 0.5 * sqrt(3.0), 0.5, 0.0);
    let r2 = mesher::add_point(3.0 + 0.5 * sqrt(3.0), 0.5, 0.0);
    let r3 = mesher::add_point(3.0, 0.0, 0.0);
    mesher::add_grid_cell([r0, r1, r2, r3]);
    mesher::add_grid_cell([s1, s1, r0, s2]);
    let bad = mesher::mesh_quality();
    if bad.inverted != 2 || bad.min_angle != 90.0 || bad.max_angle != 90.0 {
        fail(str("bad quads: ", bad.inverted, " inverted, angles ", bad.min_angle, " to ", bad.max_angle));
    }
}

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

fn integers(text: u8[:], marked: bool, list: i64[>..]&) {
    var i = 0;
    while i < text.len {
        guard syntax::is_digit(text[i]) && (!marked || (i > 0 && text[i - 1] == '#')) else { i++; continue; };
        var e = i;
        while e < text.len && syntax::is_digit(text[e]) { e++; }
        let v, good = parse_int(text[i..e]);
        list.push(v);
        i = e;
    }
}

fn check_step(what: u8[:], step: u8[:]) -> i64 {
    if !starts_with(step, "ISO-10303-21;\n") || !ends_with(step, "END-ISO-10303-21;\n") { fail(str("export ", what, ": no Part 21 frame")); }
    if find(step, "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF") < 0 { fail(str("export ", what, ": no AP242 schema")); }
    var defined: i64[>..] = [];
    var used: i64[>..] = [];
    var points = 0;
    var pnmax = -1;
    var corners: i64[>..] = [];
    for line in split(step, '\n') {
        guard line.len > 0 && line[0] == '#' else { continue; };
        let eq = position(line, '=');
        integers(line[..eq], false, defined);
        let body = line[eq + 1..];
        if starts_with(body, "COORDINATES_LIST") {
            var n: i64[>..] = [];
            integers(body[..find(body, ",(")], false, n);
            points = n[0];
            for real in split(body[find(body, ",(") + 1..], ',') {
                if position(real, '.') < 0 { fail(str("export ", what, ": the coordinate ", real, " has no point")); }
            }
        } else if starts_with(body, "TRIANGULATED_FACE") {
            let head = body[..find(body, ",(),$,(),(")];
            integers(body[find(body, ",$,(),(") + 7..], false, corners);
            integers(head, true, used);
            var n: i64[>..] = [];
            integers(head, false, n);
            pnmax = n[n.len - 1];
        } else {
            integers(body, true, used);
        }
    }
    for u in used { if !contains(defined[..], u) { fail(str("export ", what, ": #", u, " is used and not defined")); } }
    var top = 0;
    for c in corners { top = max(top, c); }
    if pnmax != points { fail(str("export ", what, ": pnmax ", pnmax, ", not the point count ", points)); }
    if points != corners.len || top != points || corners.len % 3 != 0 {
        fail(str("export ", what, ": ", corners.len, " corners, ", points, " points, the highest index ", top));
    }
    return corners.len / 3;
}

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
    let step = export_case("a block with a hole as STEP", holed, ["--geometry", "s", "--format", "step"], 0);
    let n_step = check_step("a block with a hole", step);
    if n_step != n_fine { fail(str("export STEP: ", n_step, " triangles, the STL has ", n_fine)); }
    var small: u8[>..] = [];
    mesher::step_real(small, 5.0e-5);
    if small[..] != "5.E-05" { fail(str("export STEP: 5e-05 is the REAL ", small, ", not 5.E-05")); }
    let coarse = export_case("--size 50 cm", holed, ["--geometry", "s", "--size", "50 cm"], 0);
    let n_coarse, coarse_volume, coarse_area = read_stl("--size 50 cm", coarse);
    if n_coarse == 0 || n_coarse >= n_fine { fail(str("export --size 50 cm: ", n_coarse, " facets, not fewer than ", n_fine)); }
    if mesher::norm3(coarse_area) > 1.0e-9 { fail("export --size 50 cm: the surface is open"); }
    let unmeshed = solid("box((0, 0, 0) m, (1, 1, 1) m)", "");
    let n_cube, cube_volume, cube_area = read_stl("a cube", export_case("a cube", unmeshed, ["--geometry", "s", "--size", "500mm"], 0));
    if n_cube == 0 || abs(cube_volume - 1.0) > 1.0e-12 { fail(str("export a cube: the volume is ", cube_volume, ", not 1")); }
    let plate = "model t: 2d, SI\nscales 1 m\nlet p: Surface = rect((0, 0) m, (2, 1) m)\n";
    let n_plate, plate_volume, plate_area = read_stl("a 2D rect", export_case("a 2D rect", plate, ["--geometry", "p", "--size", "0.25 m"], 0));
    if n_plate == 0 || abs(plate_area[2] - 2.0) > 1.0e-12 { fail(str("export a 2D rect: the area is ", plate_area[2], ", not 2")); }
    let plate_step = export_case("a 2D rect as STEP", plate, ["--geometry", "p", "--size", "0.25 m", "--format", "step"], 0);
    if check_step("a 2D rect", plate_step) != n_plate { fail("export STEP: a 2D rect has other triangles than its STL"); }
    export_case("an unknown let", holed, ["--geometry", "nope", "--size", "1 m"], 2);
    export_case("--format iges", holed, ["--geometry", "s", "--format", "iges"], 2);
    export_case("a size with no unit", holed, ["--geometry", "s", "--size", "3"], 2);
    export_case("no --geometry", holed, ["--size", "1 m"], 2);
    export_case("an unknown option", holed, ["--geometry", "s", "--sizes", "1 m"], 2);
    export_case("with no mesh block and no --size", unmeshed, ["--geometry", "s"], 2);
    export_case("a Length", holed, ["--geometry", "L"], 1);
    export_case("a boundary", holed, ["--geometry", "top"], 1);
    export_case("a model with an error", str(holed, "let bad: Length = 2 s\n"), ["--geometry", "s"], 1);
}

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

fn test_interfaces_3d() {
    mesher::clear_shapes();
    let slab = mesher::box(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0);
    let ball = mesher::sphere(0.0, 0.0, 0.0, 1.0);
    mesher::clear_interfaces();
    mesher::add_interface(ball);
    let done = mesher::mesh_3d(slab, 0.5);
    mesher::clear_interfaces();
    if !done { fail("interfaces 3d: the mesher did not finish"); }
    var volume_ball = 0.0;
    var volume = 0.0;
    var crossed = 0;
    var flat = 0;
    for c in mesher::cell_count() {
        var inside = 0;
        var outside = 0;
        var centre: f64[3] = [0.0; 3];
        let a = mesher::vertex(mesher::cell_node(c, 0));
        let u = mesher::sub3(mesher::vertex(mesher::cell_node(c, 1)), a);
        let v = mesher::sub3(mesher::vertex(mesher::cell_node(c, 2)), a);
        let w = mesher::sub3(mesher::vertex(mesher::cell_node(c, 3)), a);
        let size = mesher::dot3(u, mesher::cross3(v, w)) / 6.0;
        if size <= 0.0 { flat++; }
        for i in 4 {
            let p = mesher::vertex(mesher::cell_node(c, i));
            let f = mesher::phi(ball, p[0], p[1], p[2]);
            if f < -1.0e-6 { inside++; }
            if f > 1.0e-6 { outside++; }
            for t in 3 { centre[t] += p[t] / 4.0; }
        }
        if inside > 0 && outside > 0 { crossed++; }
        if mesher::phi(ball, centre[0], centre[1], centre[2]) < 0.0 { volume_ball += size; }
        volume += size;
    }
    print("interfaces 3d: ", mesher::cell_count(), " tetrahedra, ", crossed, " cross the ball, volume of the ball ", volume_ball);
    if crossed > 0 || flat > 0 { fail(str("interfaces 3d: ", crossed, " tetrahedra cross the ball, ", flat, " are not positive")); }
    if abs(volume - 64.0) > 1.0e-6 || abs(volume_ball - 4.18879) > 0.25 { fail(str("interfaces 3d: the volumes ", volume_ball, " and ", volume, " are not 4.19 and 64")); }
}

fn longest_near(d: f64, y0: f64, cx: f64, cy: f64, r: f64, circle: bool) -> f64 {
    fn near(p: f64[3]) -> bool {
        let e = if circle { abs(sqrt((p[0] - cx) * (p[0] - cx) + (p[1] - cy) * (p[1] - cy)) - r) } else { abs(p[1] - y0) };
        return e <= d + 1.0e-12;
    }
    var longest = 0.0;
    for c in mesher::cell_count() {
        for i in 3 {
            let a = mesher::vertex(mesher::cell_node(c, i));
            let b = mesher::vertex(mesher::cell_node(c, (i + 1) % 3));
            if near(a) && near(b) { longest = max(longest, mesher::norm3(mesher::sub3(a, b))); }
        }
    }
    return longest;
}

fn test_refine() {
    let square = mesher::rect(0.0, 0.0, 1.0, 1.0);
    mesher::clear_refines();
    mesher::refine_segment(0.0, 1.0, 1.0, 1.0, 1.0 / 128.0);
    let lid_ok = mesher::mesh_2d(square, 1.0 / 16.0);
    let lid = longest_near(0.0, 1.0, 0.0, 0.0, 0.0, false);
    let lid_cells = mesher::cell_count();
    print("2D refine on a segment: ", lid_cells, " cells, longest edge on the segment ", lid);
    if !lid_ok || !mesher::quality_passes(mesher::mesh_quality()) { fail("2D refine on a segment: the mesher or mesh_quality"); }
    expect_conforming("2D refine on a segment", square);
    if lid > 1.01 / 128.0 { fail(str("2D refine on a segment: an edge on the segment is ", lid)); }
    if lid_cells > 20000 { fail(str("2D refine on a segment: ", lid_cells, " cells")); }

    let channel = mesher::subtract(mesher::rect(-2.0, -2.0, 6.0, 2.0), mesher::circle(0.0, 0.0, 0.05));
    mesher::clear_refines();
    mesher::refine_circle(0.0, 0.0, 0.05, 0.0025);
    let cyl_ok = mesher::mesh_2d(channel, 0.1);
    let at_cyl = longest_near(0.0, 0.0, 0.0, 0.0, 0.05, true);
    let far = longest_near(0.5, 2.0, 0.0, 0.0, 0.0, false);
    let cyl_cells = mesher::cell_count();
    print("2D refine on a circle: ", cyl_cells, " cells, longest edge on the circle ", at_cyl, ", near y = 2 ", far);
    if !cyl_ok || !mesher::quality_passes(mesher::mesh_quality()) { fail("2D refine on a circle: the mesher or mesh_quality"); }
    expect_conforming("2D refine on a circle", channel);
    if at_cyl > 1.01 * 0.0025 { fail(str("2D refine on a circle: an edge on the circle is ", at_cyl)); }
    if far < 0.05 { fail(str("2D refine on a circle: the far edges are ", far, ", not near the size 0.1")); }
    if cyl_cells > 60000 { fail(str("2D refine on a circle: ", cyl_cells, " cells")); }
    mesher::clear_refines();
}

fn edges_in_size(size: f64) -> f64 {
    let n = mesher::point_count();
    var seen = dictionary<i64, bool> {};
    var total = 0;
    var inside = 0;
    for c in mesher::cell_count() {
        for i in 4 {
            for j in i + 1..4 {
                let a = mesher::cell_node(c, i);
                let b = mesher::cell_node(c, j);
                guard !seen.get(min(a, b) * n + max(a, b)) else { continue; }
                seen.insert(min(a, b) * n + max(a, b), true);
                let pa = mesher::vertex(a);
                let pb = mesher::vertex(b);
                let w = mesher::wanted_size(0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2]), 0.0, size);
                let r = mesher::norm3(mesher::sub3(pa, pb)) / w;
                total++;
                if r >= 0.5 && r <= 1.5 { inside++; }
            }
        }
    }
    return inside * 1.0 / total;
}

fn uncovered(p: f64[3], q: f64[3], tol: f64) -> f64 {
    let u = mesher::sub3(q, p);
    let l = mesher::norm3(u);
    fn along(x: f64[3]) -> f64, bool {
        let t = mesher::dot3(mesher::sub3(x, p), u) / (l * l);
        let off = mesher::norm3(mesher::sub3(x, [p[0] + t * u[0], p[1] + t * u[1], p[2] + t * u[2]]));
        return t, off <= tol;
    }
    var spans: f64[>..] = [];
    for c in mesher::cell_count() {
        for i in 4 {
            for j in i + 1..4 {
                let ta, on_a = along(mesher::vertex(mesher::cell_node(c, i)));
                let tb, on_b = along(mesher::vertex(mesher::cell_node(c, j)));
                if on_a && on_b { spans.append([min(ta, tb), max(ta, tb)]); }
            }
        }
    }
    var missed = 0;
    for k in 1000 {
        let t = (k + 0.5) / 1000.0;
        var hit = false;
        for s in spans.len / 2 {
            if spans[2 * s] <= t && t <= spans[2 * s + 1] { hit = true; }
        }
        if !hit { missed++; }
    }
    return missed / 1000.0;
}

fn test_refine_clip() {
    let b = mesher::cylinder(0.0, 0.0, -20.0, 0.0, 0.0, 1.0, 20.0, 20.0);
    mesher::clear_refine();
    mesher::add_refine_clip(b, 0.02, true, 2, [0.0, 0.0], 1, true, [0.0, 0.0, 1.0]);
    let rim = mesher::wanted_size(1.0, 0.0, 0.0, 0.0, 4.0);
    let off_ring = mesher::wanted_size(5.0, 0.0, 0.0, 0.0, 4.0);
    let below = mesher::wanted_size(1.0, 0.0, -5.0, 0.0, 4.0);
    mesher::clear_refine();
    print("refine clip: wanted size ", rim, " at the rim, ", off_ring, " 4 mm off the ring, ", below, " 5 mm below");
    if abs(rim - 0.02) > 1.0e-12 || off_ring < 1.0 || below < 1.0 { fail("refine clip: small at the rim, and larger away from the plane and the ring"); }
}

fn test_stuffing() {
    let dent = mesher::sphere(0.5, 0.5, 1.0, 0.3);
    let dimple = mesher::subtract(mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0), dent);
    let dimple_volume = 1.0 - 2.0 / 3.0 * PI * 0.027;
    mesh_case("dimple", dimple, 0.1, dimple_volume, tolerance(2.0 * PI * 0.09, 0.3, 4.0 + 2.0 * PI * 0.3, 0.1, dimple_volume));

    var cells = 0;
    for grade in [1.1, 1.3, 2.0] {
        mesher::clear_refine();
        mesher::set_grade(grade);
        mesher::add_refine(dent, 0.02, false);
        let what = str("dimple_ref, grade ", grade);
        if grade == 1.3 {
            mesh_case(what, dimple, 0.1, dimple_volume, tolerance(2.0 * PI * 0.09, 0.3, 4.0 + 2.0 * PI * 0.3, 0.02, dimple_volume));
        } else if !mesher::mesh_3d(dimple, 0.1) {
            fail(str(what, ": the mesher did not finish"));
        }
        let q = mesher::mesh_quality();
        print(what, ": ", q.cells, " cells");
        if cells > 0 && q.cells >= cells { fail(str(what, ": ", q.cells, " cells, not fewer than ", cells)); }
        cells = q.cells;
        if q.min_angle < 10.0 || q.max_angle > 165.0 { fail(str(what, ": dihedral angles outside 10 to 165 deg")); }
        if grade == 1.3 {
            let share = edges_in_size(0.1);
            print("  edges in [0.5, 1.5] of the wanted size: ", share);
            if share < 0.9 { fail(str(what, ": ", share, " of the edges in [0.5, 1.5] of the wanted size, not 0.9")); }
            improved_case(what, dimple, 0.1);
        }
    }
    mesher::clear_refine();

    let lblock = mesher::subtract(mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0), mesher::box(0.43, 0.43, -1.0, 2.0, 2.0, 2.0));
    mesh_case("lblock", lblock, 0.1, 0.6751, 1.0e-9);
    if mesher::mesh_quality().min_angle < 10.0 { fail("lblock: a dihedral angle below 10 deg"); }
    let edges: f64[6][7] = [[0.43, 1.0, 0.0, 0.43, 1.0, 1.0], [1.0, 0.43, 0.0, 1.0, 0.43, 1.0], [0.43, 0.43, 0.0, 0.43, 0.43, 1.0],
                            [0.43, 0.43, 1.0, 0.43, 1.0, 1.0], [0.43, 0.43, 1.0, 1.0, 0.43, 1.0], [0.43, 0.43, 0.0, 0.43, 1.0, 0.0],
                            [0.43, 0.43, 0.0, 1.0, 0.43, 0.0]];
    for e in edges {
        let missed = uncovered([e[0], e[1], e[2]], [e[3], e[4], e[5]], 1.0e-9 * 0.1);
        if missed > 0.0 { fail(str("lblock: ", missed, " of the sharp edge ", e, " is off the edges of the mesh")); }
    }

    let ball = mesher::sphere(0.0, 0.0, 0.0, 1.0);
    if !mesher::mesh_3d(ball, 0.05) { fail("ball, size 0.05: the mesher did not finish"); }
    improved_case("ball, size 0.05", ball, 0.05);
}

fn improved_case(what: u8[:], g: i64, size: f64) {
    let q = mesher::mesh_quality();
    var off = 0;
    for v in mesher::point_count() {
        let p = mesher::vertex(v);
        if mesher::on_surface[v] && abs(mesher::phi(g, p[0], p[1], p[2])) > 1.0e-9 * size { off++; }
    }
    if off > 0 { fail(str(what, ": ", off, " nodes of the surface are off the surface")); }
    mesher::improve_cells = false;
    let done = mesher::mesh_3d(g, size);
    mesher::improve_cells = true;
    let stuffed = mesher::mesh_quality();
    print(what, ": max_aspect ", stuffed.max_aspect, " stuffed, ", q.max_aspect, " improved");
    if !done || q.max_aspect > stuffed.max_aspect { fail(str(what, ": the improvement made the worst cell worse")); }
    if q.max_aspect > 1.0 / 0.35 { fail(str(what, ": a cell with gamma ", 1.0 / q.max_aspect, ", below 0.35")); }
}

fn circle_in_edges(what: u8[:], cx: f64, cy: f64, z: f64, r: f64, size: f64) {
    let n = mesher::point_count();
    var edges = dictionary<i64, i64> {};
    for c in mesher::cell_count() {
        for i in 4 {
            for j in i + 1..4 {
                let a = mesher::cell_node(c, i);
                let b = mesher::cell_node(c, j);
                edges.update(min(a, b) * n + max(a, b), 0) { it += 1; };
            }
        }
    }
    var angle: f64[>..] = [];
    var ids: i64[>..] = [];
    for v in n {
        let p = mesher::vertex(v);
        if abs(sqrt((p[0] - cx) * (p[0] - cx) + (p[1] - cy) * (p[1] - cy)) - r) < 1.0e-9 * size && abs(p[2] - z) < 1.0e-9 * size {
            angle.push(atan2(p[1] - cy, p[0] - cx));
            ids.push(v);
        }
    }
    let m = ids.len;
    print(what, ": ", m, " nodes on the circle");
    if m < 8 { fail(str(what, ": ", m, " nodes on the sharp circle")); return; }
    for i in m {
        for j in i + 1..m {
            if angle[j] < angle[i] {
                let ta = angle[i];
                angle[i] = angle[j];
                angle[j] = ta;
                let ti = ids[i];
                ids[i] = ids[j];
                ids[j] = ti;
            }
        }
    }
    var gap = 0.0;
    for i in m {
        let a = ids[i];
        let b = ids[(i + 1) % m];
        if !edges.get(min(a, b) * n + max(a, b)) { fail(str(what, ": no edge from the circle node ", i, " to the next")); return; }
        gap = max(gap, mesher::norm3(mesher::sub3(mesher::vertex(a), mesher::vertex(b))));
    }
    print("  longest edge on the circle ", gap, ", chord error ", gap * gap / (8.0 * r) / size, " h");
    if gap > 1.6 * size { fail(str(what, ": an edge of ", gap, " on the circle, over 1.6 * size")); }
}

fn edge_case(what: u8[:], g: i64, size: f64, cx: f64, cy: f64, z: f64, r: f64, must_find: bool, gate_time: bool) {
    mesher::keep_edges = false;
    let t0 = clock();
    let plain_ok = mesher::mesh_3d(g, size);
    let t1 = clock();
    let plain = mesher::mesh_quality();
    mesher::keep_edges = true;
    let t2 = clock();
    let ok = mesher::mesh_3d(g, size);
    let t3 = clock();
    let q = mesher::mesh_quality();
    if !plain_ok || !ok { fail(str(what, ": the mesher did not finish")); return; }
    print(what, ": ", q.cells, " cells (", plain.cells, " stuffed)");
    print("  time ", t3 - t2, " (", t1 - t0, " stuffed)");
    print("  dihedral ", q.min_angle, " to ", q.max_angle, ", curves ", mesher::curves_found, ", missed ", mesher::curves_missed);
    if must_find && (mesher::curves_found < 1 || mesher::curves_missed > 0) {
        fail(str(what, ": ", mesher::curves_found, " curves, ", mesher::curves_missed, " not on a chain"));
    }
    if q.min_angle < 5.0 || q.max_angle > 175.0 { fail(str(what, ": dihedral angles outside 5 to 175 deg")); }
    if q.cells > 1.5 * plain.cells { fail(str(what, ": ", q.cells, " cells, over 1.5 times ", plain.cells)); }
    if gate_time && (t3 - t2) > 3.0 * (t1 - t0) { fail(str(what, ": time ", t3 - t2, ", over 3 times ", t1 - t0)); }
    if mesher::curves_missed == 0 { circle_in_edges(what, cx, cy, z, r, size); }
}

fn test_edges() {
    let box = mesher::box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let dimple = mesher::subtract(box, mesher::sphere(0.5, 0.5, 1.0, 0.3));
    edge_case("dimple rim", dimple, 0.1, 0.5, 0.5, 1.0, 0.3, true, true);
    let z0 = 0.5 + (0.25 + 0.45 * 0.45 - 0.09) / (2.0 * 0.5);
    let r0 = sqrt(0.45 * 0.45 - (z0 - 0.5) * (z0 - 0.5));
    let cup = mesher::subtract(mesher::sphere(0.5, 0.5, 0.5, 0.45), mesher::sphere(0.5, 0.5, 1.0, 0.3));
    edge_case("cup rim", cup, 0.06, 0.5, 0.5, z0, r0, true, false);
    for size in [0.08, 0.13] {
        edge_case(str("dimple rim, size ", size), dimple, size, 0.5, 0.5, 1.0, 0.3, false, false);
    }
    let lblock = mesher::subtract(box, mesher::box(0.43, 0.43, -1.0, 2.0, 2.0, 2.0));
    mesher::keep_edges = true;
    if !mesher::mesh_3d(lblock, 0.1) { fail("lblock: the mesher did not finish"); }
    print("lblock: curves ", mesher::curves_found, ", missed ", mesher::curves_missed);
    if mesher::curves_missed > 0 { fail("lblock: a curve not on a chain"); }
}

fn expect_layers(what: u8[:], g: i64, wall: i64, first: f64, outer_area: f64) {
    let n = mesher::point_count();
    var edges = dictionary<i64, i64> {};
    var total = 0.0;
    for c in mesher::cell_count() {
        let p = mesher::vertex(mesher::cell_node(c, 0));
        let q = mesher::vertex(mesher::cell_node(c, 1));
        let r = mesher::vertex(mesher::cell_node(c, 2));
        let a2 = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
        if a2 <= 0.0 { fail(str(what, ": triangle ", c, " has the area ", 0.5 * a2)); }
        total += 0.5 * a2;
        for i in 3 {
            let a = mesher::cell_node(c, i);
            let b = mesher::cell_node(c, (i + 1) % 3);
            edges.update(min(a, b) * n + max(a, b), 0) { it += 1; };
        }
    }
    var boundary = 0.0;
    var bad = 0;
    for c in mesher::cell_count() {
        for i in 3 {
            let a = mesher::cell_node(c, i);
            let b = mesher::cell_node(c, (i + 1) % 3);
            let have = edges.get(min(a, b) * n + max(a, b));
            var count = 0;
            if have { count = have; }
            let p = mesher::vertex(a);
            let q = mesher::vertex(b);
            if count > 2 { bad++; }
            if count == 1 {
                boundary += 0.5 * (p[0] * q[1] - q[0] * p[1]);
                if abs(mesher::phi(g, p[0], p[1], 0.0)) > 1.0e-6 || abs(mesher::phi(g, q[0], q[1], 0.0)) > 1.0e-6 { bad++; }
            }
        }
    }
    if bad > 0 { fail(str(what, ": ", bad, " edges are hanging or off the boundary")); }
    if abs(total - boundary) > 1.0e-9 * total { fail(str(what, ": the area ", total, " is not the enclosed area ", boundary)); }
    if abs(total - outer_area) > 5.0e-4 * outer_area { fail(str(what, ": the area ", total, " is not ", outer_area)); }
    var wrong = 0;
    var walls = 0;
    for v in n {
        let p = mesher::vertex(v);
        if abs(mesher::phi(wall, p[0], p[1], 0.0)) > 1.0e-9 { continue; }
        walls++;
        var near = 1.0e300;
        for u in n {
            if u == v { continue; }
            let q = mesher::vertex(u);
            near = min(near, sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1])));
        }
        if abs(near - first) > 1.0e-9 * first { wrong++; }
    }
    print(what, ": ", n, " nodes, ", mesher::cell_count(), " cells, ", walls, " wall nodes, ", wrong, " with a first layer not ", first);
    if wrong > 0 { fail(str(what, ": the first layer is not ", first, " at ", wrong, " wall nodes")); }
}

fn test_layers() {
    var xy: f64[>..] = [];
    for s in 41 {
        let k = if s <= 20 { s } else { 40 - s };
        let x = 0.5 * (1.0 + cos(PI * (1.0 - k / 20.0)));
        let y = 0.6 * (0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
        xy.push(x);
        xy.push(if s == 0 || s == 40 { 0.0 } else if s <= 20 { y } else { -y });
    }
    mesher::clear_shapes();
    mesher::clear_refines();
    let foil = mesher::fill(xy[..]);
    let far = mesher::circle(0.5, 0.0, 50.0);
    let g = mesher::subtract(far, foil);
    mesher::refine_segment(0.0, 0.0, 1.0, 0.0, 0.05);
    mesher::refine_circle(0.5, 0.0, 50.0, 2.0);
    let t0 = clock();
    if !mesher::mesh_2d_layers(g, foil, 5.0, 1.0e-4, 1.2, 30) { fail("layers on a foil: the mesher did not finish"); }
    print("layers on a foil: ", clock() - t0, " s, ", mesher::bl_lk, " layers");
    let q = mesher::mesh_quality();
    print(q);
    expect_layers("layers on a foil", g, foil, 1.0e-4, PI * 2500.0 - 0.0817);

    for y in [0.0, -1.7e-18, 1.0e-18] {
        if mesher::phi(foil, -1.0e-5, y, 0.0) <= 0.0 { fail(str("sharp tip: (-1e-5, ", y, ") is inside the foil")); }
        if mesher::phi(foil, 1.00001, y, 0.0) <= 0.0 { fail(str("sharp tip: (1.00001, ", y, ") is inside the foil")); }
        if mesher::phi(foil, 0.5, y, 0.0) >= 0.0 { fail(str("sharp tip: (0.5, ", y, ") is outside the foil")); }
    }
    mesher::clear_refines();
    mesher::refine_segment(0.0, 0.0, 1.0, 0.0, 0.05);
    let box_g = mesher::subtract(mesher::rect(-2.0, -2.0, 3.0, 2.0), foil);
    if !mesher::mesh_2d_layers(box_g, foil, 1.0, 1.0e-5, 1.25, 40) { fail("layers on a sharp tip: the mesher did not finish"); }
    print("layers on a sharp tip: ", mesher::bl_lk, " layers");
    expect_layers("layers on a sharp tip", box_g, foil, 1.0e-5, 20.0 - 0.0817);

    mesher::clear_shapes();
    mesher::clear_refines();
    let cyl = mesher::circle(0.0, 0.0, 0.5);
    let tube = mesher::subtract(mesher::circle(0.0, 0.0, 20.0), cyl);
    mesher::refine_circle(0.0, 0.0, 0.5, 0.05);
    mesher::refine_circle(0.0, 0.0, 20.0, 1.0);
    if !mesher::mesh_2d_layers(tube, cyl, 2.0, 1.0e-3, 1.15, 40) { fail("layers on a cylinder: the mesher did not finish"); }
    print("layers on a cylinder: ", mesher::bl_lk, " layers");
    expect_layers("layers on a cylinder", tube, cyl, 1.0e-3, PI * (400.0 - 0.25));
    mesher::clear_refines();
}

fn main() {
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

    let coarse = mesher::subtract(mesher::rect(0.0, 0.0, 20.0, 20.0), mesher::circle(10.0, 10.0, 3.0));
    if !mesher::mesh_2d(coarse, 5.0) { fail("2D coarse hole: the mesher did not finish"); }
    expect_conforming("2D coarse hole", coarse);
    let coarse_area = chord_area(20.0, 20.0, 10.0, 10.0, 3.0);
    let qc = mesher::mesh_quality();
    print("2D coarse hole, size 5\n", qc, "\n  area of the chord polygon: ", coarse_area);
    if abs(qc.measure - coarse_area) > 1.0e-10 * coarse_area {
        fail(str("2D coarse hole: the mesh area ", qc.measure, " is not the polygon area ", coarse_area));
    }

    let l_rects = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::rect(0.0, 0.0, 2.0, 10.0));
    mesh_case("2D L of 2 rects", l_rects, 3.0, 36.0, 1.0e-10);
    let overlap = mesher::intersect(mesher::rect(0.0, 0.0, 10.0, 10.0), mesher::rect(0.0, 0.0, 4.0, 6.0));
    mesh_case("2D intersect of 2 rects", overlap, 3.0, 24.0, 1.0e-10);
    let arm: f64[8] = [0.0, 0.0, 2.0, 0.0, 2.0, 10.0, 0.0, 10.0];
    let l_polygon = mesher::union(mesher::rect(0.0, 0.0, 10.0, 2.0), mesher::polygon(arm));
    mesh_case("2D L of a rect and a polygon", l_polygon, 3.0, 36.0, 1.0e-10);

    let wide: f64[6] = [0.0, 0.0, 10.0, 0.0, 10.0, 4.0];
    let narrow: f64[6] = [0.0, 0.0, 10.0, 0.0, 10.0, 3.0];
    let twenty: f64[6] = [0.0, 0.0, 10.0, 0.0, 6.578, 2.394];
    for size in [1.0, 2.5] {
        mesh_case("2D wedge of 21.8 deg", mesher::polygon(wide), size, 20.0, 1.0e-10);
        mesh_case("2D wedge of 16.7 deg", mesher::polygon(narrow), size, 15.0, 1.0e-10);
    }
    mesh_case("2D wedge of 20 deg", mesher::polygon(twenty), 2.0, 0.5 * 10.0 * 2.394, 1.0e-10);

    let slot = mesher::union(mesher::rect(0.0, 0.0, 20.0, 10.0), mesher::circle(20.0, 5.0, 5.0));
    for size in [0.7, 2.0, 4.5] {
        if !mesher::mesh_2d(slot, size) { fail(str("2D slot, size ", size, ": the mesher did not finish")); }
        expect_conforming("2D slot", slot);
        let q = mesher::mesh_quality();
        print("2D slot, size ", size, "\n", q);
        if !mesher::quality_passes(q) { fail(str("2D slot, size ", size, ": mesh_quality")); }
        if shortest_edge() < 0.1 * size { fail(str("2D slot, size ", size, ": the shortest edge is ", shortest_edge())); }
    }

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

    let bore = mesher::cylinder(30.0, 30.0, 0.0, 0.0, 0.0, 1.0, 10.0, 10.0);
    let plate = mesher::subtract(mesher::box(0.0, 0.0, 0.0, 100.0, 60.0, 10.0), bore);
    let body = mesher::union(plate, mesher::box(60.0, 20.0, 10.0, 90.0, 40.0, 20.0));
    let body_volume = 60000.0 - PI * 1000.0 + 6000.0;
    mesh_case("3D plate with a hole and a pad", body, 2.0, body_volume,
              tolerance(TAU * 100.0, 10.0, 2.0 * TAU * 10.0, 2.0, body_volume));

    let cyl = mesher::cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 10.0, 30.0);
    let cyl_volume = PI * 3000.0;
    let e2 = mesh_case("cylinder", cyl, 2.0, cyl_volume, tolerance(TAU * 300.0, 10.0, 2.0 * TAU * 10.0, 2.0, cyl_volume));
    let e1 = mesh_case("cylinder", cyl, 1.0, cyl_volume, tolerance(TAU * 300.0, 10.0, 2.0 * TAU * 10.0, 1.0, cyl_volume));
    let order = log2(e2 / e1);
    print("cylinder: order of convergence of the volume ", order);
    if order < 1.5 || order > 2.5 { fail(str("cylinder: the order of convergence is ", order, ", not 2")); }

    let hexagon: f64[12] = [10.0, 0.0, 5.0, 8.660254037844386, -5.0, 8.660254037844386,
                            -10.0, 0.0, -5.0, -8.660254037844386, 5.0, -8.660254037844386];
    let prism = mesher::extrude(mesher::polygon(hexagon), 20.0);
    let prism_volume = 1.5 * sqrt(3.0) * 100.0 * 20.0;
    mesh_case("extruded hexagon", prism, 2.0, prism_volume, tolerance(0.0, 0.0, 240.0, 2.0, prism_volume));

    let section = mesher::rect(6.0, 0.0, 10.0, 30.0);
    let tube = mesher::revolve(section, TAU);
    let tube_volume = PI * (100.0 - 36.0) * 30.0;
    mesh_case("revolved tube", tube, 1.0, tube_volume,
              tolerance(TAU * 16.0 * 30.0, 6.0, 4.0 * TAU * 8.0, 1.0, tube_volume));
    let quarter = mesher::revolve(section, 0.5 * PI);
    mesh_case("quarter tube", quarter, 1.0, 0.25 * tube_volume,
              tolerance(0.25 * TAU * 16.0 * 30.0, 6.0, 4.0 * 0.25 * TAU * 8.0 + 4.0 * 30.0 + 4.0 * 4.0, 1.0, 0.25 * tube_volume));

    let slant = sqrt(25.0 + 400.0);
    let cone = mesher::loft(mesher::circle(0.0, 0.0, 10.0), mesher::circle(0.0, 0.0, 5.0), 20.0);
    let cone_volume = PI * 20.0 / 3.0 * (100.0 + 50.0 + 25.0);
    mesh_case("lofted circles", cone, 1.0, cone_volume,
              tolerance(PI * 15.0 * slant, 5.0, TAU * 15.0, 1.0, cone_volume));
    let pyramid = mesher::loft(mesher::rect(-5.0, -5.0, 5.0, 5.0), mesher::rect(-2.0, -2.0, 2.0, 2.0), 20.0);
    let pyramid_volume = 20.0 / 3.0 * (100.0 + 16.0 + 40.0);
    mesh_case("lofted squares", pyramid, 1.0, pyramid_volume,
              tolerance(0.0, 0.0, 56.0 + 4.0 * sqrt(400.0 + 2.0 * 9.0), 1.0, pyramid_volume));

    let row = mesher::pattern(mesher::box(0.0, 0.0, 0.0, 4.0, 4.0, 4.0), 6.0, 0.0, 0.0, 3);
    mesh_case("pattern of 3 cubes", row, 1.0, 192.0, tolerance(0.0, 0.0, 144.0, 1.0, 192.0));
    let pins = mesher::pattern_polar(mesher::cylinder(8.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 5.0), 0.0, 0.0, 1.0, 6);
    let pins_volume = 6.0 * PI * 4.0 * 5.0;
    mesh_case("polar pattern of 6 cylinders", pins, 0.5, pins_volume,
              tolerance(6.0 * TAU * 2.0 * 5.0, 2.0, 6.0 * 2.0 * TAU * 2.0, 0.5, pins_volume));
    let cross = mesher::pattern_polar(mesher::box(3.0, -1.0, -1.0, 5.0, 1.0, 1.0), 0.0, 1.0, 0.0, 4);
    mesh_case("polar pattern of 4 boxes about y", cross, 0.5, 32.0, tolerance(0.0, 0.0, 4.0 * 24.0, 0.5, 32.0));

    mesh_case("2D pattern of 4 squares", mesher::pattern(mesher::rect(0.0, 0.0, 2.0, 2.0), 3.0, 0.0, 0.0, 4), 1.0, 16.0, 1.0e-10);
    let holes = mesher::subtract(mesher::rect(0.0, 0.0, 40.0, 10.0), mesher::pattern(mesher::circle(5.0, 5.0, 2.0), 10.0, 0.0, 0.0, 4));
    let holes_area = 400.0 - 4.0 * PI * 4.0;
    mesh_case("2D plate with a row of 4 holes", holes, 0.5, holes_area,
              (TAU / 25.0) * (TAU / 25.0) / 6.0 * 4.0 * PI * 4.0 / holes_area);

    test_fill();
    test_layers();
    test_nurbs();
    test_interfaces();
    test_interfaces_3d();
    test_refine();
    test_refine_clip();
    test_stuffing();
    test_edges();

    test_grid();
    test_cells();
    test_lint();
    test_export();

    if failures > 0 { exit(1); }
    print("mesher: all tests pass");
}
