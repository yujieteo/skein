// The results of a conservation-law state: the integral and the mean of a
// field on a region of cells, and a VTU file.
import std;
import math;
import os;
namespace cons;
import .mesh;

// The integral of F(u, x, y) on the cells of region r, with the value at
// the centre of each cell. F reads the N values of the cell in u.
fn integral<F>(r: i64) -> f64 {
    guard r >= 0 && r < region_count && sol.len == dof_count() else { abort("integral: a region and a state"); }
    var total = 0.0;
    for c in cell_total {
        if cell_region[r * cell_total + c] == 0 { continue; }
        total += F(sol[c * comps..(c + 1) * comps], centre_x[c], centre_y[c]) * area[c];
    }
    return total;
}

// The mean of F on the cells of region r, weighted by the area.
fn mean<F>(r: i64) -> f64 {
    let a = region_area(r);
    return integral(r, F) / a;
}

// The area of region r.
fn region_area(r: i64) -> f64 {
    guard r >= 0 && r < region_count else { abort("region_area: a region"); }
    var total = 0.0;
    for c in cell_total { if cell_region[r * cell_total + c] != 0 { total += area[c]; } }
    return total;
}

// The words of a text of names that ';' separates.
fn field_names(names: const u8[:]) -> u8[..64][>..] {
    var words: u8[..64][>..] = [];
    var w: u8[..64] = [];
    for c in names {
        if c == ';' {
            words.push(w);
            w = [];
        } else {
            w.push(c);
        }
    }
    words.push(w);
    return words;
}

// Writes the cells to a VTU file with the field F(u, x, y, i, k) of
// component k of each of the fields that names holds (separated by ';').
// counts gives the components of each field. A field of 2 components is a
// vector.
fn write_vtu<F>(path: const u8[:], names: const u8[:], counts: i64[:]) -> bool {
    let words = field_names(names);
    guard sol.len == dof_count() && words.len == counts.len else { abort("write_vtu: a state, and a count for each name"); }
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", vertex_total, "\" NumberOfCells=\"", cell_total, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for v in vertex_total { format(out, vx[v], " ", vy[v], " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cell_total { format(out, tri[3 * c], " ", tri[3 * c + 1], " ", tri[3 * c + 2], "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cell_total { format(out, 3 * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cell_total { out.append("5\n"); }
    out.append("</DataArray>\n</Cells>\n<CellData>\n");
    for f, i in words {
        let m = if counts[i] == 2 { 3 } else { counts[i] };
        format(out, "<DataArray type=\"Float64\" Name=\"", f, "\" NumberOfComponents=\"", m, "\" format=\"ascii\">\n");
        for c in cell_total {
            for k in m {
                let v = if k < counts[i] { F(sol[c * comps..(c + 1) * comps], centre_x[c], centre_y[c], i, k) } else { 0.0 };
                format(out, v, if k + 1 < m { " " } else { "\n" });
            }
        }
        out.append("</DataArray>\n");
    }
    out.append("</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
