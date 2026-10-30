// mesh_quality (docs/checks.md): the worst cells of the mesh. In 2D, the
// minimum and the maximum angle of a triangle. In 3D, the minimum and the
// maximum dihedral angle of a tetrahedron. The aspect ratio is the
// longest edge over the inradius, scaled so that an equilateral cell has 1.
import std;
import math;
namespace mesher;
import .mesh;

struct Quality {
    cells: i64,
    nodes: i64,
    min_angle: f64,     // degrees
    max_angle: f64,     // degrees
    max_aspect: f64,
    measure: f64,       // the area in 2D, the volume in 3D
    inverted: i64,      // cells of zero or negative measure
}

// The limits of docs/checks.md: 3D 5 to 175 degrees, 2D at least 10.
fn quality_passes(q: Quality) -> bool {
    if mesh_dim == 2 { return q.inverted == 0 && q.min_angle >= 10.0; }
    return q.inverted == 0 && q.min_angle >= 5.0 && q.max_angle <= 175.0;
}

fn vertex(n: i64) -> f64[3] { return [points[3 * n], points[3 * n + 1], points[3 * n + 2]]; }

fn sub3(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
fn dot3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
fn cross3(a: f64[3], b: f64[3]) -> f64[3] {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
fn norm3(a: f64[3]) -> f64 { return sqrt(dot3(a, a)); }

fn mesh_quality() -> Quality {
    var q = Quality { cells: cell_count(), nodes: point_count(), min_angle: 180.0, max_angle: 0.0,
                      max_aspect: 0.0, measure: 0.0, inverted: 0 };
    let k = mesh_dim + 1;
    for c in cell_count() {
        var v: f64[3][4] = default<f64[3][4]>();
        for i in k { v[i] = vertex(cells[k * c + i]); }
        var longest = 0.0;
        for i in k {
            for j in i + 1..k { longest = max(longest, norm3(sub3(v[i], v[j]))); }
        }
        if mesh_dim == 2 {
            let area = 0.5 * cross3(sub3(v[1], v[0]), sub3(v[2], v[0]))[2];
            q.measure += area;
            if area <= 0.0 { q.inverted++; continue; }
            var perimeter = 0.0;
            for i in 3 {
                let a = sub3(v[(i + 1) % 3], v[i]);
                let b = sub3(v[(i + 2) % 3], v[i]);
                let angle = degrees(acos(clamp(dot3(a, b) / (norm3(a) * norm3(b)), -1.0, 1.0)));
                q.min_angle = min(q.min_angle, angle);
                q.max_angle = max(q.max_angle, angle);
                perimeter += norm3(a);
            }
            let inradius = 2.0 * area / perimeter;
            q.max_aspect = max(q.max_aspect, longest / (2.0 * sqrt(3.0) * inradius));
        } else {
            let volume = dot3(sub3(v[1], v[0]), cross3(sub3(v[2], v[0]), sub3(v[3], v[0]))) / 6.0;
            q.measure += volume;
            if volume <= 0.0 { q.inverted++; continue; }
            // n[i] is the outward normal of the face opposite vertex i, with
            // the length of twice the face area.
            var n: f64[3][4] = default<f64[3][4]>();
            var faces = 0.0;
            for i in 4 {
                let a = v[(i + 1) % 4];
                let b = v[(i + 2) % 4];
                let d = v[(i + 3) % 4];
                var f = cross3(sub3(b, a), sub3(d, a));
                if dot3(f, sub3(v[i], a)) > 0.0 { f = [-f[0], -f[1], -f[2]]; }
                n[i] = f;
                faces += 0.5 * norm3(f);
            }
            for i in 4 {
                for j in i + 1..4 {
                    let cosine = clamp(-dot3(n[i], n[j]) / (norm3(n[i]) * norm3(n[j])), -1.0, 1.0);
                    let angle = degrees(acos(cosine));
                    q.min_angle = min(q.min_angle, angle);
                    q.max_angle = max(q.max_angle, angle);
                }
            }
            let inradius = 3.0 * volume / faces;
            q.max_aspect = max(q.max_aspect, longest / (2.0 * sqrt(6.0) * inradius));
        }
    }
    return q;
}

// The mesh_quality record, in the TOON form of the run report.
fn format(out: u8[>..]&, q: Quality) {
    let angle = if mesh_dim == 2 { "angle" } else { "dihedral" };
    format(out, "mesh_quality: ", if quality_passes(q) { "pass" } else { "FAIL" },
        "\n  cells: ", q.cells, "\n  nodes: ", q.nodes, "\n  min_", angle, ": ");
    format_flt(out, q.min_angle, 2);
    format(out, " deg\n  max_", angle, ": ");
    format_flt(out, q.max_angle, 2);
    out.append(" deg\n  max_aspect: ");
    format_flt(out, q.max_aspect, 3);
    format(out, "\n  inverted: ", q.inverted, "\n  ", if mesh_dim == 2 { "area" } else { "volume" }, ": ", q.measure);
}
