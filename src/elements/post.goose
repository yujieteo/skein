// Values of the solution at a point, and VTU output.
//
// A point value is the mean over the cells that hold the point. At a node
// this is the nodal average of the cell values, which is the usual
// recovery of a stress or a flux at a node.
import std;
import math;
import os;
namespace elements;
import .shape;
import .mesh;
import .assemble;

let FIELD_VALUE = 0;    // component comp of u, or T
let FIELD_STRESS = 1;   // sigma at comp = 3 i + j
let FIELD_FLUX = 2;     // the heat flux -k grad(T) at comp

// The barycentric coordinates of p in cell c, from its vertices.
fn barycentric(c: i64, p: f64[3]) -> f64[4] {
    var e: f64[9] = [0.0; 9];
    var b: f64[3] = [0.0; 3];
    let v0 = cell_node(c, 0);
    for r in dim {
        let v = cell_node(c, r + 1);
        for t in dim { e[3 * t + r] = node_x(v, t) - node_x(v0, t); }
        b[r] = p[r] - node_x(v0, r);
    }
    var xi: f64[3] = [0.0; 3];
    if dim == 2 {
        let det = e[0] * e[4] - e[1] * e[3];
        xi[0] = (b[0] * e[4] - e[1] * b[1]) / det;
        xi[1] = (e[0] * b[1] - b[0] * e[3]) / det;
    } else {
        let det = e[0] * (e[4] * e[8] - e[5] * e[7]) - e[1] * (e[3] * e[8] - e[5] * e[6])
            + e[2] * (e[3] * e[7] - e[4] * e[6]);
        // Cramer's rule: column r of e replaced by b.
        for r in 3 {
            var m = e;
            for t in 3 { m[3 * t + r] = b[t]; }
            xi[r] = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
        }
    }
    var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
    for r in dim {
        l[r + 1] = xi[r];
        l[0] -= xi[r];
    }
    return l;
}

// The field in cell c at reference point xi.
fn cell_field(c: i64, xi: f64[3], what: i64, comp: i64) -> f64 {
    cell_geometry(c, xi);
    if what == FIELD_VALUE {
        var v = 0.0;
        if comp < ncomp {
            for l in per_cell { v += shp[l] * sol[ncomp * cell_node(c, l) + comp]; }
        }
        return v;
    }
    // du[3 a + t] = du_a / dx_t
    var du: f64[9] = [0.0; 9];
    for l in per_cell {
        let n = cell_node(c, l);
        for a in ncomp {
            for t in dim { du[3 * a + t] += sol[ncomp * n + a] * grad[3 * l + t]; }
        }
    }
    if what == FIELD_FLUX { return if comp < dim { -conductivity * du[comp] } else { 0.0 }; }
    let i = comp / 3;
    let j = comp % 3;
    var tr = 0.0;
    for a in dim { tr += du[4 * a]; }
    if i >= dim || j >= dim {
        if i == j && plane == PLANE_STRAIN { return lame_3d * tr; }
        return 0.0;
    }
    var s = shear * (du[3 * i + j] + du[3 * j + i]);
    if i == j { s += lame * tr; }
    return s;
}

// The field at point (x, y, z), and whether a cell holds the point. A point
// a little outside the mesh, as on the arc of a hole between two nodes,
// takes the nearest cell.
fn field_at(what: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard cell_count() > 0 && sol.len > 0 else { return 0.0, false; }
    let p: f64[3] = [x, y, z];
    let tol = 1.0e-9;
    var sum = 0.0;
    var count = 0;
    var best = -1;
    var best_l = -1.0e300;
    for c in cell_count() {
        let l = barycentric(c, p);
        var lowest = l[0];
        for a in 1..dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
        }
        if lowest >= -tol {
            sum += cell_field(c, [l[1], l[2], l[3]], what, comp);
            count++;
        }
    }
    if count > 0 { return sum / count, true; }
    let l = barycentric(best, p);
    return cell_field(best, [l[1], l[2], l[3]], what, comp), best_l >= -0.1;
}

// --- VTU -----------------------------------------------------------------------

// VTK cell types: the triangle, the quadratic triangle, the tetrahedron and
// the quadratic tetrahedron. Their node sequence is the one of shape.goose.
fn vtk_type() -> i64 {
    if dim == 2 { return if order == 1 { 5 } else { 22 }; }
    return if order == 1 { 10 } else { 24 };
}

fn write_vtu_file(path: const u8[:], field: const u8[:]) -> bool {
    let nodes = node_count();
    let cells = cell_count();
    var out: u8[>..] = [];
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", cells, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, node_x(n, 0), " ", node_x(n, 1), " ", node_x(n, 2), "\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cells {
        for l in per_cell { format(out, cell_node(c, l), if l + 1 < per_cell { " " } else { "\n" }); }
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cells { format(out, per_cell * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cells { format(out, vtk_type(), "\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    let comps = if ncomp == 1 { 1 } else { 3 };
    format(out, "<DataArray type=\"Float64\" Name=\"", field, "\" NumberOfComponents=\"", comps,
           "\" format=\"ascii\">\n");
    for n in nodes {
        for a in comps {
            let v = if a < ncomp { sol[ncomp * n + a] } else { 0.0 };
            format(out, v, if a + 1 < comps { " " } else { "\n" });
        }
    }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
