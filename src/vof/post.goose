import std;
import math;
import os;
namespace vof;
import .mesh;

fn locate(x: f64, y: f64) -> i64 {
    let tol = tolerance();
    guard x >= x0 - tol && x <= x0 + nx * hx + tol && y >= y0 - tol && y <= y0 + ny * hy + tol else { return -1; };
    let i = clamp(floor((x - x0) / hx) as! i64, 0, nx - 1);
    let j = clamp(floor((y - y0) / hy) as! i64, 0, ny - 1);
    return i + nx * j;
}

let FIELD_V = 1;
let FIELD_P = 2;
let FIELD_ALPHA = 4;

fn write_vtu(path: const u8[:], mask: i64) -> bool {
    guard alpha.len == cell_count() else { abort("write_vtu: a mesh and a state"); }
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", (nx + 1) * (ny + 1), "\" NumberOfCells=\"", cell_count(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for j in ny + 1 {
        for i in nx + 1 { format(out, x0 + i * hx, " ", y0 + j * hy, " 0\n"); }
    }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cell_count() {
        let p = c % nx + (nx + 1) * (c / nx);
        format(out, p, " ", p + 1, " ", p + nx + 2, " ", p + nx + 1, "\n");
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cell_count() { format(out, 4 * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cell_count() { out.append("9\n"); }
    out.append("</DataArray>\n</Cells>\n<CellData>\n");
    if (mask & FIELD_V) != 0 {
        out.append("<DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for c in cell_count() { format(out, vel[2 * c], " ", vel[2 * c + 1], " 0\n"); }
        out.append("</DataArray>\n");
    }
    if (mask & FIELD_P) != 0 { scalar_array(out, "p", pres); }
    if (mask & FIELD_ALPHA) != 0 { scalar_array(out, "alpha", alpha); }
    out.append("</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}

fn scalar_array(out: u8[>..]&, name: u8[:], q: f64[:]) {
    format(out, "<DataArray type=\"Float64\" Name=\"", name, "\" format=\"ascii\">\n");
    for v in q { format(out, v, "\n"); }
    out.append("</DataArray>\n");
}
