// Values of the solution at a point, and VTU output.
import std;
import math;
import os;
namespace flow;
import .mesh;
import .cell;
import .solve;
import .fsi;

let VALUE_VX = 0;
let VALUE_VY = 1;
let VALUE_P = 2;
let VALUE_K = 3;    // k and w of a turbulent flow
let VALUE_W = 4;
let VALUE_UX = 5;   // the displacement of the solid (fsi.goose)
let VALUE_UY = 6;

// The barycentric coordinates of (x, y) in cell c.
fn barycentric(c: i64, x: f64, y: f64) -> f64[3] {
    let g = geometry(c);
    var l: f64[3] = [0.0; 3];
    for k in 3 {
        // L_k is 0 at the next vertex i and has the gradient gl[k].
        let i = cell_node(c, (k + 1) % 3);
        l[k] = g.gl[2 * k] * (x - coords[2 * i]) + g.gl[2 * k + 1] * (y - coords[2 * i + 1]);
    }
    return l;
}

// The value `what` (a VALUE_ constant, scaled) at (x, y), and whether a
// cell holds the point.
fn value_at(what: i64, x: f64, y: f64) -> f64, bool {
    if coupled && (what == VALUE_UX || what == VALUE_UY) { return solid_value_at(what - VALUE_UX, x, y); }
    guard what >= VALUE_VX && what <= (if turbulent { VALUE_W } else { VALUE_P }) && cell_count() > 0 && sol.len == dof_count() else { abort("value_at: a value, and a solve"); }
    let tol = 1.0e-9;
    for c in cell_count() {
        let l = barycentric(c, x, y);
        if l[0] < -tol || l[1] < -tol || l[2] < -tol { continue; }
        if what == VALUE_P {
            var p = 0.0;
            for k in 3 { p += l[k] * sol[pdof(cell_node(c, k))]; }
            return p, true;
        }
        if what >= VALUE_K {
            var t = 0.0;
            for k in 3 { t += l[k] * sol[turbulence_dof(cell_node(c, k), what - VALUE_K)]; }
            return t, true;
        }
        let s = shape(geometry(c), l);
        var v = 0.0;
        for a in per_cell { v += s.n[a] * sol[vdof(cell_node(c, a), what)]; }
        return v, true;
    }
    return 0.0, false;
}

// Writes the mesh, v and p, and k and w of a turbulent flow, to a VTU
// file. The values are multiplied by the references: x_ref for the
// coordinates, v_ref and p_ref for the fields, v_ref² for k and
// v_ref/x_ref for w, so that the file is dimensional (spec §21.4). An edge node takes the
// mean pressure of the ends of its edge.
fn write_vtu(path: const u8[:], x_ref: f64, v_ref: f64, p_ref: f64) -> bool {
    guard sol.len == dof_count() else { abort("write_vtu: a solve"); }
    let nodes = node_count();
    var pressure: f64[>..] = [];
    for n in nodes { pressure.push(if n < vertex_count { sol[pdof(n)] } else { 0.0 }); }
    for c in cell_count() {
        for e in per_cell - 3 {
            let a = cell_node(c, EDGES[2 * e]);
            let b = cell_node(c, EDGES[2 * e + 1]);
            pressure[cell_node(c, 3 + e)] = 0.5 * (sol[pdof(a)] + sol[pdof(b)]);
        }
    }
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", cell_count(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, x_ref * coords[2 * n], " ", x_ref * coords[2 * n + 1], " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cell_count() {
        for l in per_cell { format(out, cell_node(c, l), if l + 1 < per_cell { " " } else { "\n" }); }
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cell_count() { format(out, per_cell * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    // VTK_TRIANGLE (5) or VTK_QUADRATIC_TRIANGLE (22), with the edge nodes
    // in the sequence of mesh.goose.
    for c in cell_count() { format(out, if order == 1 { 5 } else { 22 }, "\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    out.append("<DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, v_ref * sol[vdof(n, 0)], " ", v_ref * sol[vdof(n, 1)], " 0\n"); }
    out.append("</DataArray>\n<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n");
    for n in nodes { format(out, p_ref * pressure[n], "\n"); }
    if turbulent {
        out.append("</DataArray>\n<DataArray type=\"Float64\" Name=\"k\" format=\"ascii\">\n");
        for n in nodes { format(out, v_ref * v_ref * sol[turbulence_dof(n, 0)], "\n"); }
        out.append("</DataArray>\n<DataArray type=\"Float64\" Name=\"w\" format=\"ascii\">\n");
        for n in nodes { format(out, v_ref / x_ref * sol[turbulence_dof(n, 1)], "\n"); }
    }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
