import std;
namespace elements;
import .shape;
import .mesh;
import .assemble;

fn elastic_cube(c: i64) {
    guard c >= 0 && c < cell_count() && physics == ELASTIC && dim == 3 && hex_cells && ncomp == 3 else { abort("elastic_cube: a cell, the elastic physics and hexahedra"); }
    let m = 3 * per_cell;
    for &v in ke { v = 0.0; }
    for &v in fe { v = 0.0; }
    for qp in quad_count(3) {
        let xi, w = quad_point(3, qp);
        let dv = w * cell_geometry(c, xi) * thickness;
        for i in per_cell {
            for p in 3 { fe[3 * i + p] += volume_force[p] * shp[i] * dv; }
            for j in per_cell {
                var dot = 0.0;
                for t in 3 { dot += grad[3 * i + t] * grad[3 * j + t]; }
                for p in 3 {
                    for q in 3 {
                        var v = lame * grad[3 * i + p] * grad[3 * j + q] + shear * grad[3 * i + q] * grad[3 * j + p];
                        if p == q { v += shear * dot; }
                        ke[m * (3 * i + p) + 3 * j + q] += v * dv;
                    }
                }
            }
        }
    }
}

fn heat_cube(c: i64) {
    guard c >= 0 && c < cell_count() && physics == HEAT && dim == 3 && hex_cells && ncomp == 1 else { abort("heat_cube: a cell, the heat physics and hexahedra"); }
    for &v in ke { v = 0.0; }
    for &v in fe { v = 0.0; }
    for qp in quad_count(3) {
        let xi, w = quad_point(3, qp);
        let dv = w * cell_geometry(c, xi) * thickness;
        let kq = conductivity_at(c);
        for i in per_cell {
            fe[i] += source * shp[i] * dv;
            for j in per_cell {
                var dot = 0.0;
                for t in 3 { dot += grad[3 * i + t] * grad[3 * j + t]; }
                ke[per_cell * i + j] += kq * dot * dv;
            }
        }
    }
}
