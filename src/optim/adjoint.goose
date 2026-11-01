// The objective of the package optim and its adjoint. The objective is
// J = a u_c(p): component c of u at a point p, times a. A study adds the
// terms that hold no u (gen.goose). The adjoint λ solves K λ = dJ/du = a w
// on the free dofs, with the weights w of p in its cell, and λ = 0 on the
// fixed dofs. With the residual r = f - K u of a value v that the step
// reads (E, nu, t, a load), dJ/dv = λᵀ ∂r/∂v at the u of the solve: one
// solve of the adjoint for each design gives the derivative by each value
// and by the density of each cell. Reason: a finite difference needs one
// solve for each value, and a grid of 1200 cells has 1200 densities.
import std;
import math;
namespace optim;
import .q4;
import .band;
import solvers from "../solvers";

var objective_comp = -1;
var objective_at: f64[2] = [0.0, 0.0];
var objective_scale = 0.0;
var adjoint: f64[>..] = [];   // λ

// The objective a u_comp(x, y). False if the point is outside the grid.
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

// λ = K⁻¹ dJ/du, with the factor of the last solve.
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

// dJ/dx of the density x of cell c: -λᵀ (dK/dx) u, with
// dK/dx = p x^(p - 1) (1 - E_MIN) E k0 on the cell.
fn cell_sensitivity(c: i64) -> f64 {
    let x = density[c];
    let dm = young * (1.0 - E_MIN) * penalty * pow(x, penalty - 1.0);
    return -dm * cell_product(k0, c, adjoint, u);
}

// dJ/dv = λᵀ ∂r/∂v of each value v that the step reads: E, nu, t, then
// the 2 components of each load. K is linear in E and t, the tractions in
// t, and the loads in their components. ∂K/∂nu comes from ∂D/∂nu.
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
