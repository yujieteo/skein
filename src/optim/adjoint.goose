import std;
import math;
namespace optim;
import .q4;
import .band;
import solvers from "../solvers";

var objective_comp = -1;
var objective_at: f64[2] = [0.0, 0.0];
var objective_scale = 0.0;
var adjoint: f64[>..] = [];

fn objective(comp: i64, x: f64, y: f64, a: f64) -> bool {
    guard nx > 0 && (comp == X || comp == Y) else { abort("objective: a mesh and a component X or Y"); }
    objective_comp = comp;
    objective_at = [x, y];
    objective_scale = a;
    let _, _, found = point_weights(x, y);
    return found;
}

fn objective_value() -> f64 {
    guard objective_comp >= 0 else { abort("objective_value: an objective"); }
    let v, _ = value_at(objective_comp, objective_at[0], objective_at[1]);
    return objective_scale * v;
}

fn solve_adjoint() {
    guard objective_comp >= 0 && u.len == dof_count() else { abort("solve_adjoint: an objective and a solve"); }
    let n, w, _ = point_weights(objective_at[0], objective_at[1]);
    var g: f64[>..] = [];
    solvers::zeros(g, dof_count());
    for k in 4 { g[2 * n[k] + objective_comp] += objective_scale * w[k]; }
    for d in g.len { if fixed[d] != 0 { g[d] = 0.0; } }
    if adjoint.len != g.len { solvers::zeros(adjoint, g.len); }
    substitute(g, adjoint);
}

fn cell_sensitivity(c: i64) -> f64 {
    let x = density[c];
    let dm = young * (1.0 - E_MIN) * penalty * pow(x, penalty - 1.0);
    return -dm * cell_product(k0, c, adjoint, u);
}

fn input_gradient() -> f64[>..] {
    guard adjoint.len == dof_count() else { abort("input_gradient: solve_adjoint"); }
    let dk = cell_matrix(d_by_nu(1.0, poisson));
    var ku = 0.0;
    var dnu = 0.0;
    for c in cell_count() {
        ku += cell_modulus(c) * cell_product(k0, c, adjoint, u);
        dnu -= cell_modulus(c) * cell_product(dk, c, adjoint, u);
    }
    var f: f64[>..] = [];
    solvers::zeros(f, dof_count());
    for l in loads { if !l.point { add_load(l, 1.0, f); } }
    var out: f64[>..] = [-ku / young, dnu, (solvers::dot(adjoint, f) - ku) / thickness];
    for l in loads {
        for c in 2 {
            var unit = l;
            unit.f = [0.0, 0.0];
            unit.f[c] = 1.0;
            for &v in f { v = 0.0; }
            add_load(unit, 1.0, f);
            out.push(solvers::dot(adjoint, f));
        }
    }
    return out;
}
