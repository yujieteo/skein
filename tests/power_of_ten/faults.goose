// Calls of the runtime with a broken precondition, for the test of rule 5
// of the Power of Ten in scripts/power_of_ten.sh. The argument names the
// call. Each call must stop the program with the abort of its precondition.
import std;
import os;
import runtime from "../../src/runtime";
import elements from "../../src/elements";
import shells from "../../src/shells";
import flow from "../../src/flow";
import solvers from "../../src/solvers";
import mesher from "../../src/mesher";

fn main() {
    let a = args();
    guard a.len == 2 else { exit(2); }
    let call = a[1];
    if call == "runtime_material" {
        runtime::material(210000.0, 0.5);
    } else if call == "runtime_cells_for" {
        runtime::cells_for(10.0, 0.0);
    } else if call == "runtime_mesh_box" {
        runtime::mesh_box(0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 1, 1, 1);
    } else if call == "runtime_fix" {
        runtime::fix("x0");
    } else if call == "runtime_solve_static" {
        runtime::material(210000.0, 0.3);
        runtime::solve_static();
    } else if call == "elements_elastic" {
        elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
        elements::elastic(210000.0, 0.5);
    } else if call == "elements_fix" {
        elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
        elements::heat(50.0, 0.0, 1.0);
        elements::plane(0, 0.0);
        elements::fix(1, 0.0);
    } else if call == "elements_traction" {
        elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
        elements::plane(0, 0.0);
        elements::traction(1.0, 0.0, 0.0);
    } else if call == "elements_solve" {
        elements::mesh_box(3, 1, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
        elements::solve();
    } else if call == "shells_laminate_ply" {
        shells::laminate_ply(0.0, 0.0, 1.0, 1.0, 1.0, 0.3, 1.0, 1.0);
    } else if call == "shells_fix" {
        shells::fix(0, 0);
    } else if call == "shells_force" {
        shells::force(0, 0, 1.0);
    } else if call == "shells_mesh_cylinder" {
        shells::mesh_cylinder(1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 0, 1.0, 0.0, 0.0);
    } else if call == "flow_mesh" {
        flow::mesh_from_arrays([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], [0, 1, 3], 2);
    } else if call == "flow_fluid" {
        flow::fluid(1.0, 0.0, true, flow::SUPG);
    } else if call == "flow_fix" {
        flow::fix_velocity(0, 1.0);
    } else if call == "flow_solve" {
        flow::solve_steady(1.0e-8);
    } else if call == "flow_value" {
        flow::value_at(flow::VALUE_P, 0.0, 0.0);
    } else if call == "solvers_multiply" {
        var y: f64[2] = [0.0, 0.0];
        solvers::multiply([0, 1], [0], [1.0], [1.0], y);
    } else if call == "solvers_ilu0_apply" {
        var z: f64[1] = [0.0];
        solvers::ilu0_apply([0, 1], [0], [1.0, 2.0], [0], [1.0], z);
    } else {
        exit(2);
    }
}
