// Values of the solution at a point, the L2 error against a field, and
// VTU output.
import std;
import math;
import os;
namespace grid;
import .mesh;
import .dg;

// The cell that holds (x, y), and the reference point there, or -1.
fn locate(x: f64, y: f64) -> i64, f64[2] {
    let tol = tolerance();
    let i = floor((x - x0) / hx) as! i64;
    let j = floor((y - y0) / hy) as! i64;
    let ci = clamp(i, 0, nx - 1);
    let cj = clamp(j, 0, ny - 1);
    let c = ci + nx * cj;
    let p = centre(c);
    let xi = 2.0 * (x - p[0]) / hx;
    let eta = 2.0 * (y - p[1]) / hy;
    if abs(xi) > 1.0 + 2.0 * tol / hx || abs(eta) > 1.0 + 2.0 * tol / hy { return -1, [0.0, 0.0]; }
    return c, [clamp(xi, -1.0, 1.0), clamp(eta, -1.0, 1.0)];
}

// The value of u in cell c at the reference point (ξ, η): the value of
// the cell (fv), or the sum of the basis (dg).
fn value_in(c: i64, xi: f64, eta: f64) -> f64 {
    guard method == FV || method == DG else { abort("value_in: a method"); }
    if method == FV { return sol[c]; }
    let n = per_cell();
    let f = basis(xi, eta);
    var v = 0.0;
    for a in n { v += sol[c * n + a] * f.v[a]; }
    return v;
}

// The value of u at (x, y), and whether a cell holds the point.
fn value_at(x: f64, y: f64) -> f64, bool {
    guard sol.len == dof_count() && method >= 0 else { abort("value_at: a method and a state"); }
    let c, p = locate(x, y);
    if c < 0 { return 0.0, false; }
    return value_in(c, p[0], p[1]), true;
}

// The L2 norm of u - F(x, y, t) at the time of the state, and the L2 norm
// of F. fv takes the value at the centre of each cell, which is second
// order, and dg the Gauss rule of p + 2 points on each axis.
fn l2_error<F>() -> f64, f64 {
    guard sol.len == dof_count() && method >= 0 else { abort("l2_error: a method and a state"); }
    var e = 0.0;
    var f2 = 0.0;
    for c in cell_count() {
        if method == FV {
            let p = centre(c);
            let v = F(p[0], p[1], time);
            e += (sol[c] - v) * (sol[c] - v) * hx * hy;
            f2 += v * v * hx * hy;
            continue;
        }
        for qj in points() {
            for qi in points() {
                let p = cell_point(c, gauss_x(qi), gauss_x(qj));
                let w = gauss_w(qi) * gauss_w(qj) * 0.25 * hx * hy;
                let v = F(p[0], p[1], time);
                let f = vol_basis[qi + points() * qj];
                var d = -v;
                for a in per_cell() { d += sol[c * per_cell() + a] * f.v[a]; }
                e += w * d * d;
                f2 += w * v * v;
            }
        }
    }
    return sqrt(e), sqrt(f2);
}

// Writes the cells and u to a VTU file: 4 points for each cell, so that a
// dg field keeps its jumps, and u of the cell at each of them.
fn write_vtu(path: const u8[:], name: const u8[:]) -> bool {
    guard sol.len == dof_count() && method >= 0 else { abort("write_vtu: a method and a state"); }
    let cells = cell_count();
    let corners: f64[8] = [-1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0];
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", 4 * cells, "\" NumberOfCells=\"", cells, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for c in cells {
        for k in 4 {
            let p = cell_point(c, corners[2 * k], corners[2 * k + 1]);
            format(out, p[0], " ", p[1], " 0\n");
        }
    }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cells { format(out, 4 * c, " ", 4 * c + 1, " ", 4 * c + 2, " ", 4 * c + 3, "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cells { format(out, 4 * (c + 1), "\n"); }
    // VTK_QUAD (9).
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cells { out.append("9\n"); }
    format(out, "</DataArray>\n</Cells>\n<PointData>\n<DataArray type=\"Float64\" Name=\"", name, "\" format=\"ascii\">\n");
    for c in cells {
        for k in 4 { format(out, value_in(c, corners[2 * k], corners[2 * k + 1]), "\n"); }
    }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
