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
import mesher from "../mesher";

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
                if i != skip { f.push(mesher::cells[k * c + i]); }
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
            let a = mesher::cells[3 * c + i];
            let b = mesher::cells[3 * c + (i + 1) % 3];
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

    if failures > 0 { exit(1); }
    print("mesher: all tests pass");
}
