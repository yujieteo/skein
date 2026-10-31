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
import .stages;

let FIELD_VALUE = 0;    // component comp of u, or T
let FIELD_STRESS = 1;   // sigma at comp = 3 i + j
let FIELD_FLUX = 2;     // the heat flux -k grad(T) at comp

// The barycentric coordinates of p in cell c, from its vertices.
fn barycentric(c: i64, p: f64[3]) -> f64[4] {
    guard c >= 0 && c < cell_count() && (dim == 2 || dim == 3) else { abort("barycentric: c < cell_count()"); }
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
    guard what >= FIELD_VALUE && what <= FIELD_FLUX && comp >= 0 && comp < 9 && ncomp <= 3 && sol.len == dof_total() else { abort("cell_field: a field, 0 <= comp < 9 and a solution of each dof"); }
    cell_geometry(c, xi);
    use_law(c);
    if what == FIELD_VALUE {
        var v = 0.0;
        if comp < ncomp {
            for l in per_cell { v += shp[l] * sol[ncomp * cell_node(c, l) + comp]; }
        }
        return v;
    }
    // du[3 a + t] = du_a / dx_t
    // A born cell is strain-free at its birth (stages.goose).
    let b = born.len > 0 && born[c] != 0;
    var du: f64[9] = [0.0; 9];
    for l in per_cell {
        let n = cell_node(c, l);
        for a in ncomp {
            let v = sol[ncomp * n + a] - (if b { birth_u[(per_cell * c + l) * ncomp + a] } else { 0.0 });
            for t in dim { du[3 * a + t] += v * grad[3 * l + t]; }
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
    guard what >= FIELD_VALUE && what <= FIELD_FLUX && comp >= 0 && comp < 9 else { abort("field_at: a field and 0 <= comp < 9"); }
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

// The largest von Mises stress at the nodes of the region, from the
// stress of each cell at the node. 0 if the region holds no node.
fn region_max_von_mises() -> f64 {
    guard physics == ELASTIC && sol.len == dof_total() else { abort("region_max_von_mises: an elastic model and a solution"); }
    var most = 0.0;
    for c in cell_count() {
        for l in per_cell {
            let n = cell_node(c, l);
            if !on_region(n) { continue; }
            let b = barycentric(c, [node_x(n, 0), node_x(n, 1), node_x(n, 2)]);
            var s: f64[9] = [0.0; 9];
            for k in 9 { s[k] = cell_field(c, [b[1], b[2], b[3]], FIELD_STRESS, k); }
            let d = (s[0] - s[4]) * (s[0] - s[4]) + (s[4] - s[8]) * (s[4] - s[8]) + (s[8] - s[0]) * (s[8] - s[0]);
            let t = s[1] * s[1] + s[5] * s[5] + s[2] * s[2];
            most = max(most, sqrt(0.5 * d + 3.0 * t));
        }
    }
    return most;
}

// --- VTU -----------------------------------------------------------------------

// VTK cell types: the triangle, the quadratic triangle, the tetrahedron and
// the quadratic tetrahedron. Their node sequence is the one of shape.goose.
fn vtk_type() -> i64 {
    if dim == 2 { return if order == 1 { 5 } else { 22 }; }
    return if order == 1 { 10 } else { 24 };
}

fn write_vtu_file(path: const u8[:], field: const u8[:]) -> bool {
    guard ncomp >= 1 && sol.len == dof_total() else { abort("write_vtu_file: a solution of each dof"); }
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
