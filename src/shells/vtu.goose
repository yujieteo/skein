// VTU output (VTK XML unstructured grid, ASCII), which ParaView reads: the
// shell cells as quads and the beams as lines, with the displacement and
// the rotation of each node as point arrays.
import std;
import os;
namespace shells;
import .model;

fn write_vtu(path: const u8[:], w: const u8[:], theta: const u8[:]) -> bool {
    let nodes = node_count();
    var out: u8[>..] = [];
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", shells.len + beams.len, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, coords[3 * n], " ", coords[3 * n + 1], " ", coords[3 * n + 2], "\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in shells { format(out, c.n[0], " ", c.n[1], " ", c.n[2], " ", c.n[3], "\n"); }
    for b in beams { format(out, b.n[0], " ", b.n[1], "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    var offset = 0;
    for c in shells.len {
        offset += 4;
        format(out, offset, "\n");
    }
    for b in beams.len {
        offset += 2;
        format(out, offset, "\n");
    }
    // VTK types 9 (quad) and 3 (line).
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in shells.len { out.append("9\n"); }
    for b in beams.len { out.append("3\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    for f in 2 {
        format(out, "<DataArray type=\"Float64\" Name=\"", if f == 0 { w } else { theta },
            "\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for n in nodes {
            format(out, disp[6 * n + 3 * f], " ", disp[6 * n + 3 * f + 1], " ", disp[6 * n + 3 * f + 2], "\n");
        }
        out.append("</DataArray>\n");
    }
    out.append("</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
