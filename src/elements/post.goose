import std;
import math;
import os;
namespace elements;
import .shape;
import .mesh;
import .assemble;
import .stages;
import .field;
import .edge;

let FIELD_VALUE = 0;
let FIELD_STRESS = 1;
let FIELD_FLUX = 2;

fn cube_coordinates(c: i64, p: f64[3]) -> f64[4] {
    var xi: f64[3] = [0.5; 3];
    var n: f64[20] = [0.0; 20];
    var dn: f64[60] = [0.0; 60];
    for it in 8 {
        shape(3, 1, xi, n, dn);
        var r: f64[3] = [0.0; 3];
        var j: f64[9] = [0.0; 9];
        for l in 8 {
            let v = cell_node(c, l);
            for t in 3 {
                r[t] += n[l] * node_x(v, t);
                for a in 3 { j[3 * t + a] += dn[3 * l + a] * node_x(v, t); }
            }
        }
        for t in 3 { r[t] -= p[t]; }
        let det = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6]) + j[2] * (j[3] * j[7] - j[4] * j[6]);
        guard det > 0.0 else { break; };
        for a in 3 {
            var m = j;
            for t in 3 { m[3 * t + a] = r[t]; }
            xi[a] -= (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
        }
    }
    return [min(1.0 - xi[0], min(1.0 - xi[1], 1.0 - xi[2])), xi[0], xi[1], xi[2]];
}

fn newton_refine(c: i64, p: f64[3], guess: f64[3]) -> f64[3] {
    var xi = guess;
    let affine = xi;
    var n: f64[20] = [0.0; 20];
    var dn: f64[60] = [0.0; 60];
    for it in 6 {
        shape(dim, 2, xi, n, dn);
        var r: f64[3] = [0.0; 3];
        var j: f64[9] = [0.0; 9];
        for l in per_cell {
            let v = cell_node(c, l);
            for t in dim {
                r[t] += n[l] * node_x(v, t);
                for a in dim { j[3 * t + a] += dn[3 * l + a] * node_x(v, t); }
            }
        }
        for t in dim { r[t] -= p[t]; }
        var step: f64[3] = [0.0; 3];
        if dim == 2 {
            let det = j[0] * j[4] - j[1] * j[3];
            guard abs(det) > 0.0 else { break; };
            step[0] = (r[0] * j[4] - j[1] * r[1]) / det;
            step[1] = (j[0] * r[1] - r[0] * j[3]) / det;
        } else {
            let det = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6]) + j[2] * (j[3] * j[7] - j[4] * j[6]);
            guard abs(det) > 0.0 else { break; };
            for a in 3 {
                var m = j;
                for t in 3 { m[3 * t + a] = r[t]; }
                step[a] = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
            }
        }
        for a in dim { xi[a] -= step[a]; }
    }
    for a in dim { if !(xi[a] > -0.5 && xi[a] < 1.5) { xi = affine; } }
    shape(dim, 2, xi, n, dn);
    var back: f64[9] = [0.0; 9];
    for l in per_cell {
        for t in dim { for a in dim { back[3 * t + a] += dn[3 * l + a] * node_x(cell_node(c, l), t); } }
    }
    let folded = if dim == 2 { back[0] * back[4] - back[1] * back[3] <= 0.0 } else {
        back[0] * (back[4] * back[8] - back[5] * back[7]) - back[1] * (back[3] * back[8] - back[5] * back[6])
            + back[2] * (back[3] * back[7] - back[4] * back[6]) <= 0.0 };
    var miss = 0.0;
    var span = 0.0;
    for t in dim {
        var x = 0.0;
        for l in per_cell { x += n[l] * node_x(cell_node(c, l), t); }
        miss = max(miss, abs(x - p[t]));
        span = max(span, abs(node_x(cell_node(c, 1), t) - node_x(cell_node(c, 0), t)));
    }
    if folded || miss > 1.0e-6 * span { xi = affine; }
    return xi;
}

fn barycentric(c: i64, p: f64[3]) -> f64[4] {
    guard c >= 0 && c < cell_count() && (dim == 2 || dim == 3) else { abort("barycentric: c < cell_count()"); }
    if hex_cells { return cube_coordinates(c, p); }
    var e: f64[9] = [0.0; 9];
    var b: f64[3] = [0.0; 3];
    let v0 = cell_node(c, 0);
    for r in dim {
        let v = cell_node(c, r + 1);
        for t in dim { e[3 * t + r] = node_x(v, t) - node_x(v0, t); }
        b[r] = p[r] - node_x(v0, r);
    }
    var xi: f64[3] = [0.0; 3];
    if dim == 2 {
        let det = e[0] * e[4] - e[1] * e[3];
        xi[0] = (b[0] * e[4] - e[1] * b[1]) / det;
        xi[1] = (e[0] * b[1] - b[0] * e[3]) / det;
    } else {
        let det = e[0] * (e[4] * e[8] - e[5] * e[7]) - e[1] * (e[3] * e[8] - e[5] * e[6])
            + e[2] * (e[3] * e[7] - e[4] * e[6]);
        for r in 3 {
            var m = e;
            for t in 3 { m[3 * t + r] = b[t]; }
            xi[r] = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
        }
    }
    if order == 2 { xi = newton_refine(c, p, xi); }
    var l: f64[4] = [1.0, 0.0, 0.0, 0.0];
    for r in dim {
        l[r + 1] = xi[r];
        l[0] -= xi[r];
    }
    return l;
}

fn cell_field(c: i64, xi: f64[3], what: i64, comp: i64) -> f64 {
    guard what >= FIELD_VALUE && what <= FIELD_FLUX && comp >= 0 && comp < 9 && ncomp <= 3 && sol.len == dof_total() else { abort("cell_field: a field, 0 <= comp < 9 and a solution of each dof"); }
    if what == FIELD_VALUE { shape(dim, order, xi, shp, dref); } else { cell_geometry(c, xi); }
    use_law(c);
    if what == FIELD_VALUE {
        var v = 0.0;
        if comp < ncomp {
            for l in per_cell { v += shp[l] * shown(ncomp * cell_node(c, l) + comp); }
        }
        return v;
    }
    let b = born.len > 0 && born[c] != 0;
    var du: f64[9] = [0.0; 9];
    for l in per_cell {
        let n = cell_node(c, l);
        for a in ncomp {
            let v = sol[ncomp * n + a] - (if b { birth_u[(per_cell * c + l) * ncomp + a] } else { 0.0 });
            for t in dim { du[3 * a + t] += v * grad[3 * l + t]; }
        }
    }
    if what == FIELD_FLUX { return if comp < dim { -conductivity_at(c) * du[comp] } else { 0.0 }; }
    let i = comp / 3;
    let j = comp % 3;
    var tr = 0.0;
    for a in dim { tr += du[4 * a]; }
    if i >= dim || j >= dim {
        if aniso_on && plane == PLANE_STRAIN { return aniso_stress(du, i, j); }
        if i == j && plane == PLANE_STRAIN { return lame_3d * tr; }
        return 0.0;
    }
    var s = shear * (du[3 * i + j] + du[3 * j + i]);
    if i == j { s += lame * tr; }
    if aniso_on { return aniso_stress(du, i, j); }
    if i == j && ft_kind.len > 0 { s -= strain_stress(c); }
    return s;
}

fn field_at(what: i64, comp: i64, x: f64, y: f64, z: f64) -> f64, bool {
    guard what >= FIELD_VALUE && what <= FIELD_FLUX && comp >= 0 && comp < 9 else { abort("field_at: a field and 0 <= comp < 9"); }
    guard cell_count() > 0 && sol.len > 0 else { return 0.0, false; }
    let p: f64[3] = [x, y, z];
    let tol = 1.0e-9;
    var sum = 0.0;
    var count = 0;
    var best = -1;
    var best_l = -1.0e300;
    for c in cell_count() {
        if !is_live(c) { continue; }
        let l = barycentric(c, p);
        var lowest = l[0];
        for a in 1..dim + 1 { lowest = min(lowest, l[a]); }
        if lowest > best_l {
            best_l = lowest;
            best = c;
        }
        if lowest >= -tol {
            sum += cell_field(c, [l[1], l[2], l[3]], what, comp);
            count++;
        }
    }
    if count > 0 { return sum / count, true; }
    guard best >= 0 else { return 0.0, false; };
    let l = barycentric(best, p);
    return cell_field(best, [l[1], l[2], l[3]], what, comp), best_l >= -0.1;
}

fn region_max_von_mises() -> f64 {
    guard physics == ELASTIC && sol.len == dof_total() else { abort("region_max_von_mises: an elastic model and a solution"); }
    var most = 0.0;
    for c in cell_count() {
        if !is_live(c) { continue; }
        for l in per_cell {
            let n = cell_node(c, l);
            if !on_region(n) { continue; }
            let b = barycentric(c, [node_x(n, 0), node_x(n, 1), node_x(n, 2)]);
            var s: f64[9] = [0.0; 9];
            for k in 9 { s[k] = cell_field(c, [b[1], b[2], b[3]], FIELD_STRESS, k); }
            let d = (s[0] - s[4]) * (s[0] - s[4]) + (s[4] - s[8]) * (s[4] - s[8]) + (s[8] - s[0]) * (s[8] - s[0]);
            let t = s[1] * s[1] + s[5] * s[5] + s[2] * s[2];
            most = max(most, sqrt(0.5 * d + 3.0 * t));
        }
    }
    return most;
}

fn region_integral(code: i64[:], consts: f64[:]) -> f64 {
    guard physics >= 0 && sol.len == dof_total() && ncomp <= 3 else { abort("region_integral: a physics and its solution"); }
    let fd = dim - 1;
    let fnodes = node_count_of(fd, order);
    var fn_shape: f64[20] = [0.0; 20];
    var total = 0.0;
    var any = false;
    for c in cell_count() {
        for f in facet_count(dim) {
            if !region_facet(c, f) { continue; }
            any = true;
            for qp in quad_count(fd) {
                let dw = facet_point(c, f, qp, fn_shape);
                var x: f64[3] = [0.0; 3];
                var u: f64[3] = [0.0; 3];
                for l in fnodes {
                    let n = cell_node(c, facet_node(dim, f, l));
                    for t in 3 { x[t] += fn_shape[l] * node_x(n, t); }
                    for k in ncomp { u[k] += fn_shape[l] * shown(ncomp * n + k); }
                }
                total += eval_facet_code(code, consts, x, u) * dw;
            }
        }
    }
    return if any { total } else { 0.0 / 0.0 };
}

fn region_extreme(comp: i64, sign: f64) -> f64 {
    guard sol.len == dof_total() && comp >= 0 && comp < ncomp else { abort("region_extreme: a solution and a component of it"); }
    var most = 0.0 / 0.0;
    for n in node_count() {
        guard on_region(n) else { continue; };
        let v = sign * shown(ncomp * n + comp);
        if !(most >= v) { most = v; }
    }
    return most;
}

fn region_max_stress(i: i64, j: i64, sign: f64) -> f64 {
    guard physics == ELASTIC && sol.len == dof_total() && i >= 0 && i < 3 && j >= 0 && j < 3 else { abort("region_max_stress: an elastic model, a solution and 0 <= i, j < 3"); }
    var sum: f64[>..] = [];
    var count: f64[>..] = [];
    for n in node_count() {
        sum.push(0.0);
        count.push(0.0);
    }
    for c in cell_count() {
        if !is_live(c) { continue; }
        for l in per_cell {
            let n = cell_node(c, l);
            if !on_region(n) { continue; }
            let b = barycentric(c, [node_x(n, 0), node_x(n, 1), node_x(n, 2)]);
            sum[n] += sign * cell_field(c, [b[1], b[2], b[3]], FIELD_STRESS, 3 * i + j);
            count[n] += 1.0;
        }
    }
    var most = 0.0 / 0.0;
    for n in node_count() {
        if count[n] > 0.0 && !(most >= sum[n] / count[n]) { most = sum[n] / count[n]; }
    }
    return most;
}

fn vtk_type() -> i64 {
    if dim == 2 { return if order == 1 { 5 } else { 22 }; }
    if hex_cells { return if order == 2 { 25 } else { 12 }; }
    return if order == 1 { 10 } else { 24 };
}

fn vtu_mesh(out: u8[>..]&) {
    let nodes = node_count();
    let cells = cell_count();
    format(out, """
        <?xml version="1.0"?>
        <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
        """, "\n");
    format(out, "<Piece NumberOfPoints=\"", nodes, "\" NumberOfCells=\"", cells, "\">\n");
    out.append("<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes { format(out, node_x(n, 0), " ", node_x(n, 1), " ", node_x(n, 2), "\n"); }
    out.append("</DataArray>\n</Points>\n<Cells>\n");
    out.append("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for c in cells {
        for l in per_cell { format(out, cell_node(c, l), if l + 1 < per_cell { " " } else { "\n" }); }
    }
    out.append("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for c in cells { format(out, per_cell * (c + 1), "\n"); }
    out.append("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for c in cells { format(out, vtk_type(), "\n"); }
    out.append("</DataArray>\n</Cells>\n<PointData>\n");
}

let VTU_END = "</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";

fn write_vtu_file(path: const u8[:], field: const u8[:]) -> bool {
    if physics == EDGE { return write_edge_vtu(path, field); }
    guard ncomp >= 1 && sol.len == dof_total() else { abort("write_vtu_file: a solution of each dof"); }
    var out: u8[>..] = [];
    vtu_mesh(out);
    let comps = if ncomp == 1 { 1 } else { 3 };
    format(out, "<DataArray type=\"Float64\" Name=\"", field, "\" NumberOfComponents=\"", comps,
           "\" format=\"ascii\">\n");
    for n in node_count() {
        for a in comps {
            let v = if a < ncomp { shown(ncomp * n + a) } else { 0.0 };
            format(out, v, if a + 1 < comps { " " } else { "\n" });
        }
    }
    out.append("</DataArray>\n");
    out.append(VTU_END);
    return write_file(path, out);
}
