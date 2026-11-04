import std;
namespace runtime;
import .solve;
import .fem;
import .mesh;
import solvers from "../solvers";

let GMRES_RESTART = 300;

var direct = false;

fn solvers_solve(solver: u8[:], preconditioner: u8[:], b: f64[:], x: f64[>..]&, tol: f64,
                 max_it: i64) -> bool, bool, i64 {
    if solver == "cg" && preconditioner == "amg" && direct { return direct_solve(b, x, tol, max_it); }
    if solver == "cg" && preconditioner == "amg" { return amg_solve(b, x, tol, max_it); }
    guard solver == "gmres" && preconditioner == "ilu0" else { return false, false, 0; }
    var a: f64[>..] = [];
    for v in val { a.push(v); }
    var rhs: f64[>..] = [];
    for v in b { rhs.push(v); }
    solvers::constrain(row_start, col, a, fixed, rhs);
    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    guard solvers::ilu0(row_start, col, a, lu, diag) else {
        write_stderr("runtime: ILU(0) has a zero pivot\n");
        return true, false, 0;
    }
    let its, converged = solvers::gmres(row_start, col, a, rhs, x, tol, GMRES_RESTART, max_it) { r, z =>
        solvers::ilu0_apply(row_start, col, lu, diag, r, z);
    };
    if !converged {
        write_stderr(str("runtime: gmres with ilu0 did not converge in ", its,
                         " iterations with a restart of ", GMRES_RESTART, "\n"));
    }
    return true, converged, its;
}

fn amg_precondition(r: f64[:], z: f64[:]) {
    solvers::amg_apply(r, z);
    for d in z.len { if fixed[d] != 0 { z[d] = 0.0; } }
}

fn amg_solve(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, bool, i64 {
    var a: f64[>..] = [];
    a.append(val);
    for d in fixed.len {
        for e in row_start[d]..row_start[d + 1] {
            if col[e] != d && (fixed[d] != 0 || fixed[col[e]] != 0) { a[e] = 0.0; }
        }
    }
    let key = solvers::warm_key();
    if !solvers::amg_reuse(key, row_start, col, a, fixed) {
        guard rigid_setup(a) else {
            write_stderr("runtime: the coarsest level of AMG is singular\n");
            return true, false, 0;
        }
        solvers::amg_keep(key, row_start, col, fixed);
    }
    if solvers::warm_start(key, x, fixed) {
        var ax = zeros(x.len);
        free_operator(x, ax);
        let xax = dot(x, ax);
        let alpha = if xax > 0.0 { dot(x, b) / xax } else { 0.0 };
        for d in x.len { if fixed[d] == 0 { x[d] *= alpha; } }
    }
    let converged, its = cg(b, x, tol, max_it, free_operator, amg_precondition);
    if converged { solvers::warm_keep(key, x); }
    return true, converged, its;
}

fn direct_solve(b: f64[:], x: f64[>..]&, tol: f64, max_it: i64) -> bool, bool, i64 {
    let key = solvers::warm_key();
    guard !solvers::nested_reuse(key, row_start, col, val, fixed, b, x) else { return true, true, 1; }
    guard solvers::nested_factor(row_start, col, val, fixed) else {
        write_stderr(str("runtime: the direct solve has a pivot with no stiffness at dof ", solvers::nested_singular(), ", AMG with CG takes the solve\n"));
        solvers::nested_free();
        check_iterative_memory(col.len);
        return amg_solve(b, x, tol, max_it);
    }
    solvers::nested_keep(key, row_start, col, val, fixed);
    solvers::nested_solve(b, x, fixed);
    return true, true, 1;
}

fn rigid_setup(a: f64[:]) -> bool {
    let modes, nb = solvers::near_null_space(mode_kind(), coords, 3, dofs_per_node, fixed);
    return solvers::amg_setup_modes(row_start, col, a, dofs_per_node, modes, nb);
}

fn mode_kind() -> i64 { return if dofs_per_node == 3 { solvers::ELASTIC_MODES } else { solvers::SCALAR_MODES }; }

fn check_memory(entries: i64) {
    direct = preconditioner_name == "amg" && solver_name == "cg" && node_place.len == node_count() &&
             solvers::nested_plan(lattice, node_place, dofs_per_node, row_start, col) && solvers::nested_pays(2);
    if direct {
        solvers::memory_check("lattice ldlt", solvers::nested_entries(), dof_total(), 2, 0, "");
        return;
    }
    check_iterative_memory(entries);
}

fn check_iterative_memory(entries: i64) {
    let modes = if preconditioner_name == "amg" { solvers::mode_count(mode_kind(), coords, 3, dofs_per_node) } else { 0 };
    let vectors = if solver_name == "gmres" { GMRES_RESTART + 8 } else { 7 + 2 * modes };
    let bending = if modes > 6 { " bending" } else { "" };
    let fix = if preconditioner_name == "jacobi" { "" } else { ", or write solve { linear: cg, precond: jacobi } in the step: it needs less memory, but it can need many more iterations" };
    let mesh = 8 * (coords.len + cell_start.len + cell_nodes.len + cell_young.len + cell_poisson.len) + cell_type.len;
    solvers::memory_check(str("box ", solver_name, " ", preconditioner_name, bending), entries, dof_total(), vectors, mesh, fix);
}
