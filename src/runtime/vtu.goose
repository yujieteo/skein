// VTU output (VTK XML unstructured grid, ASCII), which ParaView reads.
// Each 27-node cell is written as its 8 linear sub-hexahedra, so the file
// holds all nodes with no cell type of higher order.
import std;
import os;
namespace runtime;
import .fem;
import .mesh;
import .solve;

let VTK_HEXAHEDRON = 12;

fn write_vtu_file(path: const u8[:]) -> bool {
    let nodes = node_count();
    let sub = (lattice[0] - 1) * (lattice[1] - 1) * (lattice[2] - 1);
    var out: u8[>..] = [];
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", sub, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes {
        format(out, coords[3 * n], " ", coords[3 * n + 1], " ", coords[3 * n + 2], "\n");
    }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for k in lattice[2] - 1 {
        for j in lattice[1] - 1 {
            for i in lattice[0] - 1 {
                for dk in 2 {
                    format(out, node_at(i, j, k + dk), " ", node_at(i + 1, j, k + dk), " ",
                        node_at(i + 1, j + 1, k + dk), " ", node_at(i, j + 1, k + dk), " ");
                }
                out.push('\n');
            }
        }
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in sub { format(out, 8 * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in sub { format(out, VTK_HEXAHEDRON, "\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData Vectors=\"u\">\n");
    out.append("<DataArray type=\"Float64\" Name=\"u\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes {
        format(out, disp[3 * n], " ", disp[3 * n + 1], " ", disp[3 * n + 2], "\n");
    }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
