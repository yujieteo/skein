// VTU output (VTK XML unstructured grid, ASCII), which ParaView reads: the
// shell cells as quads (4 quads for each 9-node cell) and the beams as lines,
// with the displacement and the rotation of each node as point arrays, and
// the point and cell arrays that the caller gives as XML.
import std;
import os;
namespace shells;
import .model;

fn write_vtu(path: const u8[:], w: const u8[:], theta: const u8[:], points: const u8[:], cells: const u8[:]) -> bool {
    guard disp.len == dof_count() else { abort("write_vtu: a solution for each dof"); }
    let nodes = node_count();
    var out: u8[>..] = [];
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    var quads = 0;
    for c in shells { quads += if c.count == 9 { 4 } else { 1 }; }
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", quads + beams.len, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, coords[3 * n], " ", coords[3 * n + 1], " ", coords[3 * n + 2], "\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in shells { write_quads(out, c); }
    for b in beams { format(out, b.n[0], " ", b.n[1], "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    var offset = 0;
    for c in quads {
        offset += 4;
        format(out, offset, "\n");
    }
    for b in beams.len {
        offset += 2;
        format(out, offset, "\n");
    }
    // VTK types 9 (quad) and 3 (line).
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in quads { out.append("9\n"); }
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
    out.append(points);
    out.append("</PointData>\n");
    if cells.len > 0 {
        out.append("<CellData>\n");
        out.append(cells);
        out.append("</CellData>\n");
    }
    out.append("</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}

// Writes the connectivity of shell cell c: one quad, or 4 quads for a
// 9-node cell, so that each viewer of linear cells draws it. LATTICE9
// gives the node of place (i, j) of the 3 × 3 nodes of the cell.
let LATTICE9 = [0, 4, 1, 7, 8, 5, 3, 6, 2];

fn write_quads(out: u8[>..]&, c: Shell) {
    if c.count == 4 {
        format(out, c.n[0], " ", c.n[1], " ", c.n[2], " ", c.n[3], "\n");
        return;
    }
    for q in 4 {
        let i = q % 2;
        let j = q / 2;
        format(out, c.n[LATTICE9[i + 3 * j]], " ", c.n[LATTICE9[i + 1 + 3 * j]], " ",
               c.n[LATTICE9[i + 1 + 3 * (j + 1)]], " ", c.n[LATTICE9[i + 3 * (j + 1)]], "\n");
    }
}
