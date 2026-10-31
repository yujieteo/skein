// VTU output (VTK XML unstructured grid, ASCII), which ParaView reads.
// Each cell is written with the VTK type of its cell type, and its nodes
// in the order of its registry entry. A 27-node cell is written as its 8
// linear sub-hexahedra, so the file holds all nodes with no cell type of
// higher order. Each field of point_fields (fem.goose) is a point array.
import std;
import os;
namespace runtime;
import .fem;
import .solve;

// The VTK cell type of each cell type (fem.goose): triangle, quadratic
// triangle, tetra, quadratic tetra, hexahedron, hexahedron (the sub-cells
// of hex27), quad, line, hexahedron (hex8i) and quad.
let VTK_TYPES: i64[10] = [5, 22, 10, 24, 12, 12, 9, 3, 12, 9];

// The local nodes of sub-hexahedron (i, j, k) of a hex27 cell, in the VTK
// order of a hexahedron.
fn hex27_sub_cell(out: u8[>..]&, c: i64, i: i64, j: i64, k: i64) {
    for dk in 2 {
        let z = 9 * (k + dk);
        format(out, cell_node(c, i + 3 * j + z), " ", cell_node(c, i + 1 + 3 * j + z), " ",
            cell_node(c, i + 1 + 3 * (j + 1) + z), " ", cell_node(c, i + 3 * (j + 1) + z), " ");
    }
    out.push('\n');
}

fn write_vtu_file(path: const u8[:]) -> bool {
    let nodes = node_count();
    var cells = 0;
    for t in cell_type { cells += if t as i64 == HEX27 { 8 } else { 1 }; }
    var out: u8[>..] = [];
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", cells, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes {
        format(out, coords[3 * n], " ", coords[3 * n + 1], " ", coords[3 * n + 2], "\n");
    }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cell_count() {
        if cell_type[c] as i64 == HEX27 {
            for k in 2 {
                for j in 2 {
                    for i in 2 { hex27_sub_cell(out, c, i, j, k); }
                }
            }
            continue;
        }
        for e in cell_start[c]..cell_start[c + 1] { format(out, cell_nodes[e], " "); }
        out.push('\n');
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    var offset = 0;
    for c in cell_count() {
        let t = cell_type[c] as i64;
        let parts = if t == HEX27 { 8 } else { 1 };
        let size = if t == HEX27 { 8 } else { CELL_NODES[t] };
        for p in parts {
            offset += size;
            format(out, offset, "\n");
        }
    }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cell_count() {
        let t = cell_type[c] as i64;
        for p in (if t == HEX27 { 8 } else { 1 }) { format(out, VTK_TYPES[t], "\n"); }
    }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
    for f in point_fields {
        format(out, "<DataArray type=\"Float64\" Name=\"", f.name, "\" NumberOfComponents=\"", f.count,
            "\" format=\"ascii\">\n");
        for n in nodes {
            for k in f.count {
                if k > 0 { out.push(' '); }
                format(out, disp[dofs_per_node * n + f.first + k]);
            }
            out.push('\n');
        }
        out.append("</DataArray>\n");
    }
    out.append("</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
