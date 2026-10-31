// Values of the solution at a point, the stress, the L2 error against a
// field, the reactions, and VTU output. rkpm gives u at x as Σ Ψ_j(x) d_j,
// and its gradient as the smoothed gradient of the part of a nodal domain
// that holds x: the strain of the stiffness. Reason: on the plate of
// samples/05_meshless_plate.skn with a smaller support, the gradient of
// Σ Ψ_j d_j at the hole was 10% below the Kirsch stress, and the smoothed
// gradient 5%. gfd
// gives the values at the points, with the linear interpolation of each
// triangle between them.
import std;
import math;
import os;
namespace meshless;
import .cloud;
import .rk;
import .model;
import .rkpm;

// The triangle that holds (x, y), and its barycentric coordinates there,
// or -1.
fn locate(x: f64, y: f64) -> i64, f64[3] {
    let tol = tolerance();
    for t in triangle_count() {
        let a = tri[3 * t];
        let b = tri[3 * t + 1];
        let c = tri[3 * t + 2];
        let s = twice_area(a, b, c);
        let la = ((px(b) - x) * (py(c) - y) - (py(b) - y) * (px(c) - x)) / s;
        let lb = ((px(c) - x) * (py(a) - y) - (py(c) - y) * (px(a) - x)) / s;
        let lc = 1.0 - la - lb;
        let e = tol / sqrt(abs(s));
        if la >= -e && lb >= -e && lc >= -e { return t, [la, lb, lc]; }
    }
    return -1, [0.0, 0.0, 0.0];
}

// Component c of u at (x, y), and its gradient, in the triangle t with
// the barycentric coordinates l. False if rkpm has no shape functions
// there.
fn field_in(c: i64, x: f64, y: f64, t: i64, l: f64[3]) -> f64, f64[2], bool {
    if method == GFD {
        var v = 0.0;
        for k in 3 { v += l[k] * sol[ncomp * tri[3 * t + k] + c]; }
        return v, [0.0, 0.0], true;
    }
    if !shape_at(x, y) { return 0.0, [0.0, 0.0], false; }
    var v = 0.0;
    for j, i in near { v += psi[i] * sol[ncomp * j + c]; }
    var k = 0;
    for m in 3 { if l[m] > l[k] { k = m; } }
    let r = 3 * t + k;
    var g: f64[2] = [0.0, 0.0];
    for i in rec_start[r]..rec_start[r + 1] {
        let d = sol[ncomp * item[i] + c];
        g[0] += item_dx[i] * d;
        g[1] += item_dy[i] * d;
    }
    return v, g, true;
}

// Component c of u at (x, y), and whether the cloud holds the point.
fn value_at(c: i64, x: f64, y: f64) -> f64, bool {
    guard method >= 0 && c >= 0 && c < ncomp && sol.len == dof_count() else { abort("value_at: a method and a component of the unknown"); }
    let t, l = locate(x, y);
    if t < 0 { return 0.0, false; }
    let v, g, found = field_in(c, x, y, t, l);
    return v, found;
}

// Component (i, j) of the stress at (x, y), with 0, 1 and 2 for x, y and
// z, and whether the cloud holds the point.
fn stress_at(i: i64, j: i64, x: f64, y: f64) -> f64, bool {
    guard method == RKPM && physics == ELASTIC && i >= 0 && i < 3 && j >= 0 && j < 3 else { abort("stress_at: rkpm, the elastic physics, and 0 <= i, j < 3"); }
    let t, l = locate(x, y);
    if t < 0 { return 0.0, false; }
    let ux, gx, fx = field_in(0, x, y, t, l);
    let uy, gy, fy = field_in(1, x, y, t, l);
    if !fx { return 0.0, false; }
    let exx = gx[0];
    let eyy = gy[1];
    if i == 2 && j == 2 { return if plane_stress { 0.0 } else { lame_3d * (exx + eyy) }, true; }
    if i == 2 || j == 2 { return 0.0, true; }
    if i != j { return shear * (gx[1] + gy[0]), true; }
    return lame * (exx + eyy) + 2.0 * shear * (if i == 0 { exx } else { eyy }), true;
}

// The L2 norm of u - F(x, y, t), and the L2 norm of F, by the Gauss rule
// of 3 points on each triangle.
fn l2_error<F>() -> f64, f64 {
    guard method >= 0 && physics == SCALAR && sol.len == dof_count() else { abort("l2_error: a method and the scalar physics"); }
    var e = 0.0;
    var f2 = 0.0;
    for t in triangle_count() {
        let w = twice_area(tri[3 * t], tri[3 * t + 1], tri[3 * t + 2]) / 6.0;
        for k in 3 {
            var l: f64[3] = [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];
            l[k] = 2.0 / 3.0;
            var x = 0.0;
            var y = 0.0;
            for m in 3 {
                x += l[m] * px(tri[3 * t + m]);
                y += l[m] * py(tri[3 * t + m]);
            }
            let u, g, found = field_in(0, x, y, t, l);
            let v = F(x, y, time);
            e += w * (u - v) * (u - v);
            f2 += w * v * v;
        }
    }
    return sqrt(e), sqrt(f2);
}

// The reaction of each component on the selected edges (rkpm).
fn reaction() -> f64[2] {
    guard method == RKPM && selected.len == edge_count() else { abort("reaction: rkpm"); }
    let r, size = rkpm_reaction(selected);
    return r;
}

// The reaction of each component on all edges, the sum of the sizes of
// its parts, and the load of each component (rkpm).
fn balance() -> f64[2], f64, f64[2] {
    guard method == RKPM else { abort("balance: rkpm"); }
    var all: u8[>..] = [];
    for e in edge_count() { all.push(1); }
    var r, size = rkpm_reaction(all);
    var f: f64[2] = [0.0, 0.0];
    for e in edge_count() {
        let n, l = edge_normal(e);
        for c in ncomp {
            if edge_fix[ncomp * e + c] != FREE { continue; }
            let v = edge_load[ncomp * e + c] * l * thickness;
            f[c] += v;
            size += abs(v);
        }
    }
    return r, size, f;
}

// Writes the points, the triangles and u at each point to a VTU file.
fn write_vtu(path: const u8[:], name: const u8[:]) -> bool {
    guard method >= 0 && sol.len == dof_count() else { abort("write_vtu: a method and a state"); }
    let n = point_count();
    var out: u8[>..] = [];
    out.append("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    format(out, "<Piece NumberOfPoints=\"", n, "\" NumberOfCells=\"", triangle_count(), "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for p in n { format(out, px(p), " ", py(p), " 0\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for t in triangle_count() { format(out, tri[3 * t], " ", tri[3 * t + 1], " ", tri[3 * t + 2], "\n"); }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for t in triangle_count() { format(out, 3 * (t + 1), "\n"); }
    // VTK_TRIANGLE (5).
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for t in triangle_count() { out.append("5\n"); }
    format(out, "</DataArray>\n</Cells>\n<PointData>\n<DataArray type=\"Float64\" Name=\"", name, "\" NumberOfComponents=\"",
           if ncomp == 1 { 1 } else { 3 }, "\" format=\"ascii\">\n");
    for p in n {
        let t = point_tris[point_start[p]];
        var l: f64[3] = [0.0, 0.0, 0.0];
        for k in 3 { if tri[3 * t + k] == p { l[k] = 1.0; } }
        for c in ncomp {
            let v, g, found = field_in(c, px(p), py(p), t, l);
            format(out, v, if c + 1 < ncomp { " " } else if ncomp == 1 { "\n" } else { " 0\n" });
        }
    }
    out.append("</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    return write_file(path, out);
}
