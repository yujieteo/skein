// Newton's method for F(u) = 0, with a backtracking line search.
// The caller gives the residual R(u, f), which sets f = F(u), and the
// Jacobian J(u, val), which sets the values of dF/du on a fixed CSR
// pattern. Each step solves J du = -F by GMRES with ILU(0).
import std;
import math;
namespace solvers;
import .sparse;
import .ilu;
import .krylov;

let NEWTON_RESTART = 50;            // m of GMRES(m) in each step
let NEWTON_LINEAR_MAX = 5000;       // the GMRES iterations in each step
let NEWTON_SUFFICIENT = 1.0e-4;     // the Armijo constant of the line search
let NEWTON_HALVINGS = 20;           // the most halvings of a step

// Starts at `u`, and stops when the norm of F(u) is at most `tol`.
// `history` gets the norm of F at the start and after each step.
// `linear_rtol` is the tolerance of GMRES relative to the norm of F: a
// small value keeps the quadratic convergence. Returns the number of
// steps and true if it converged.
fn newton<R, J>(row_start: i64[:], col: i32[:], val: f64[:], u: f64[:], tol: f64,
                max_steps: i64, linear_rtol: f64,
                history: f64[>..]&) -> i64, bool {
    guard row_start.len == u.len + 1 && col.len == val.len && max_steps > 0 else { abort("newton: row_start has u.len + 1 entries, col and val one size, max_steps > 0"); }
    let n = u.len;
    var f: f64[>..] = [];
    var rhs: f64[>..] = [];
    var du: f64[>..] = [];
    var trial: f64[>..] = [];
    var lu: f64[>..] = [];
    var diag: i64[>..] = [];
    zeros(f, n);
    zeros(rhs, n);
    zeros(du, n);
    zeros(trial, n);
    history.clear();
    R(u, f);
    var fnorm = norm(f);
    history.push(fnorm);
    var step = 0;
    while fnorm > tol {
        if step >= max_steps { return step, false; }
        J(u, val);
        if !ilu0(row_start, col, val, lu, diag) { return step, false; }
        for d in n {
            rhs[d] = -f[d];
            du[d] = 0.0;
        }
        let its, solved = gmres(row_start, col, val, rhs, du, linear_rtol, NEWTON_RESTART,
                                NEWTON_LINEAR_MAX) { r, z => ilu0_apply(row_start, col, lu, diag, r, z); };
        if !solved { return step, false; }
        var lambda = 1.0;
        var halvings = 0;
        while true {
            for d in n { trial[d] = u[d] + lambda * du[d]; }
            R(trial, f);
            let tn = norm(f);
            if tn <= (1.0 - NEWTON_SUFFICIENT * lambda) * fnorm { fnorm = tn; break; }
            halvings++;
            if halvings > NEWTON_HALVINGS { return step, false; }
            lambda *= 0.5;
        }
        for d in n { u[d] = trial[d]; }
        step++;
        history.push(fnorm);
    }
    return step, true;
}
