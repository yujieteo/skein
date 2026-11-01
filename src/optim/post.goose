// The results of the package optim: the VTU file of u at the nodes and
// of the density of each cell.
import std;
import os;
namespace optim;
import .q4;

// Writes u (as `name`) and the density `rho` to the ASCII VTU file `path`.
// False if the file cannot be written.
fn write_vtu(path: const u8[:], name: const u8[:]) -> bool {
    guard u.len == dof_count() else { abort("write_vtu: a solve"); }
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", node_count(), "\" NumberOfCells=\"", cell_count(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    var xy: f64[>..] = [];
    for n in 2 * node_count() { xy.push(0.0); }
    for i in nx + 1 {
        for j in ny + 1 {
            xy[2 * node(i, j)] = x0 + i * dx;
            xy[2 * node(i, j) + 1] = y0 + j * dy;
        }
    }
    for n in node_count() { format(out, xy[2 * n], " ", xy[2 * n + 1], " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cell_count() {
        let n = cell_nodes(c);
        format(out, n[0], " ", n[1], " ", n[2], " ", n[3], "\n");
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cell_count() { format(out, 4 * (c + 1), "\n"); }
    // VTK_QUAD (9).
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cell_count() { out.append("9\n"); }
    format(out, "</DataArray>\n</Cells>\n<PointData>\n<DataArray type=\"Float64\" Name=\"", name, "\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in node_count() { format(out, u[2 * n], " ", u[2 * n + 1], " 0\n"); }
    out.append("</DataArray>\n</PointData>\n<CellData>\n<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n");
    for c in cell_count() { format(out, density[c], "\n"); }
    out.append("</DataArray>\n</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
