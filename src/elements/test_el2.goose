import std;
import math;
namespace elements;
import .api;
import .mesh;
import .assemble;
import .shape;

fn el2_check(what: u8[:], ok: bool) -> i64 {
    if ok { return 0; }
    print("FAIL ", what);
    return 1;
}

fn el2_midpoints_straight() -> i64 {
    mesh_box(3, 2, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);
    var worst = 0.0;
    for c in cell_count() {
        for e in 6 {
            let a = cell_node(c, edge_node(3, e, 0));
            let b = cell_node(c, edge_node(3, e, 1));
            let m = cell_node(c, 4 + e);
            for t in 3 { worst = max(worst, abs(node_x(m, t) - 0.5 * (node_x(a, t) + node_x(b, t)))); }
        }
    }
    return el2_check("mesh_box order 2: each midpoint node is the middle of its edge", worst < 1.0e-12);
}

fn el2_cells_for() -> i64 {
    return el2_check("cells_for rounds up to a cell size of at most size",
                 cells_for(1.0, 0.25) == 4 && cells_for(1.0, 0.3) == 4 && cells_for(0.1, 1.0) == 1);
}

fn el2_geometry_tolerance() -> i64 {
    mesh_box(3, 1, 0.0, 0.0, 0.0, 40.0, 10.0, 10.0, 2, 1, 1);
    return el2_check("geometry_tolerance is 1e-6 of the largest side of the box",
                 abs(geometry_tolerance() - 40.0e-6) < 1.0e-12);
}

fn el2_hole_corner_node() -> i64 {
    mesh_hole(2, 1, 0.0, 0.0, 0.0, 40.0, 30.0, 0.0, 1.0, 4.0, 0.5);
    var found = false;
    for n in node_count() {
        if abs(node_x(n, 0) - 40.0) < 1.0e-9 && abs(node_x(n, 1) - 30.0) < 1.0e-9 { found = true; }
    }
    return el2_check("mesh_hole has a node on the corner (hi.x, hi.y)", found);
}

fn el2_failures() -> i64 {
    return el2_midpoints_straight() + el2_cells_for() + el2_geometry_tolerance() + el2_hole_corner_node();
}
