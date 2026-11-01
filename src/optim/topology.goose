// The topology loop of femahoy gap G1, as the 88 line code with the
// sensitivity filter: E. Andreassen, A. Clausen, M. Schevenels, B. S.
// Lazarov and O. Sigmund, "Efficient topology optimization in MATLAB
// using 88 lines of code", Struct. Multidisc. Optim. 43 (2011) 1-16.
//
// Each cell has a density x, and the modulus E (E_MIN + (1 - E_MIN) x^p)
// (SIMP). Each iteration solves the design, gives the sensitivities of
// the objective by the adjoint (adjoint.goose), filters them, and moves
// each density by the optimality criteria update. A bisection finds the
// multiplier that keeps the volume at its limit. The loop stops when no
// density changes by more than STOP_CHANGE. The densities start at the
// volume limit. The update needs dJ/dx <= 0, as for a compliance.
import std;
import math;
namespace optim;
import .q4;
import .band;
import .adjoint;

let MOVE = 0.2;              // the largest change of a density in an iteration
let STOP_CHANGE = 0.01;

// The history of the loop: the objective, the volume fraction and the
// largest change of a density in each iteration.
var history_objective: f64[>..] = [];
var history_volume: f64[>..] = [];
var history_change: f64[>..] = [];

// The weights H(e, f) = max(0, r - |centre e - centre f|) of the filter,
// of the cells f in the square of ±reach cells about e.
fn filter_weight(e: i64, f: i64, r: f64) -> f64 {
    let ddx = (e / ny - f / ny) * dx;
    let ddy = (e % ny - f % ny) * dy;
    return max(0.0, r - sqrt(ddx * ddx + ddy * ddy));
}

// The sensitivity filter: dc(e) = Σ H(e, f) x(f) dc(f) / (max(1e-3, x(e)) Σ H(e, f)).
fn filtered(dc: f64[:], r: f64) -> f64[>..] {
    let reach_x = max(0, ceil(r / dx) as! i64 - 1);
    let reach_y = max(0, ceil(r / dy) as! i64 - 1);
    var out: f64[>..] = [];
    for e in cell_count() {
        let i = e / ny;
        let j = e % ny;
        var sum = 0.0;
        var weights = 0.0;
        for k in max(0, i - reach_x)..min(nx, i + reach_x + 1) {
            for m in max(0, j - reach_y)..min(ny, j + reach_y + 1) {
                let h = filter_weight(e, k * ny + m, r);
                sum += h * density[k * ny + m] * dc[k * ny + m];
                weights += h;
            }
        }
        out.push(sum / (max(1.0e-3, density[e]) * weights));
    }
    return out;
}

// The optimality criteria update of the densities, with the volume limit
// `volume` (a fraction). It gives the new densities.
fn update(dc: f64[:], volume: f64) -> f64[>..] {
    var lo = 0.0;
    var hi = 1.0e9;
    var x: f64[>..] = [];
    for e in cell_count() { x.push(0.0); }
    while (hi - lo) / (lo + hi) > 1.0e-3 {
        let mid = 0.5 * (hi + lo);
        var sum = 0.0;
        for e in cell_count() {
            let old = density[e];
            let b = old * sqrt(max(0.0, -dc[e]) / mid);
            x[e] = max(0.0, max(old - MOVE, min(1.0, min(old + MOVE, b))));
            sum += x[e];
        }
        if sum > volume * cell_count() { lo = mid; } else { hi = mid; }
    }
    return x;
}

// Runs the loop with the volume limit `volume`, the power `p` and the
// filter radius `r`, for at most `most` iterations. At the end, the
// densities and u are those of the last solve. False if a solve finds
// that the conditions do not hold the body.
fn topology(volume: f64, p: f64, r: f64, most: i64) -> bool {
    guard volume > 0.0 && volume <= 1.0 && p >= 1.0 && r > 0.0 && most >= 1 && objective_comp >= 0 else { abort("topology: 0 < volume <= 1, p >= 1, r > 0, most >= 1 and an objective"); }
    penalty = p;
    for &x in density { x = volume; }
    history_objective.clear();
    history_volume.clear();
    history_change.clear();
    var dc: f64[>..] = [];
    for e in cell_count() { dc.push(0.0); }
    for it in most {
        guard solve() else { return false; }
        solve_adjoint();
        for e in cell_count() { dc[e] = cell_sensitivity(e); }
        let next = update(filtered(dc, r), volume);
        var change = 0.0;
        var sum = 0.0;
        for e in cell_count() {
            change = max(change, abs(next[e] - density[e]));
            sum += density[e];
        }
        history_objective.push(objective_value());
        history_volume.push(sum / cell_count());
        history_change.push(change);
        guard change > STOP_CHANGE && it + 1 < most else { return true; }
        for e in cell_count() { density[e] = next[e]; }
    }
    return true;
}
