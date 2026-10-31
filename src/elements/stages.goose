// The materials of the cells, and birth and death (spec §11, §14.3).
//
// A cell is in a shape g of src/mesher/ if the centroid of its vertices is
// in g. The mesh of mesh_shape keeps the boundary of each shape of
// mesher::add_interface as edges in 2D, so each cell is then in or out of
// the shape as a whole.
//
// Law 0 is the material of set_elastic and body_force. add_law gives the
// cells in a shape another one.
//
// A dead cell adds no stiffness and no load. A node that no live cell holds
// is held at its value of the last step: fixed is 2 there, and the
// reactions and the check `equilibrium` skip it. A born cell is
// strain-free at its birth (spec §14.3): it stores the unknown of its nodes
// then, u_b, and its stress is C (strain(u) - strain(u_b)). The assembly
// adds K_e u_b to the load of its nodes for this. The new nodes that no
// live cell held take the value of the nearest node that one held.
import std;
import math;
namespace elements;
import mesher from "../mesher";
import .mesh;
import .assemble;

struct Law {
    lame: f64,
    lame_3d: f64,
    shear: f64,
    force: f64[3],
}

var laws: Law[>..] = [];
var cell_law: i64[>..] = [];    // the law of each cell; empty for law 0 in each cell
var dead: u8[>..] = [];         // 1 for each dead cell; empty if each cell is live
var born: u8[>..] = [];         // 1 for each cell with a birth; empty if none
var birth_u: f64[>..] = [];     // u_b of each node of each cell, per_cell ncomp values
var held: f64[>..] = [];        // the solution of the last step

// Clears the laws and the live set, for a new physics.
fn clear_stages() {
    laws.clear();
    cell_law.clear();
    dead.clear();
    born.clear();
    birth_u.clear();
    held.clear();
}

fn law_now() -> Law { return Law { lame: lame, lame_3d: lame_3d, shear: shear, force: volume_force }; }

// Sets the moduli and the body force of the assembly to the law of cell c.
fn use_law(c: i64) {
    guard cell_law.len > 0 else { return; }
    let w = laws[cell_law[c]];
    lame = w.lame;
    lame_3d = w.lame_3d;
    shear = w.shear;
    volume_force = w.force;
}

// Sets them to law 0 again.
fn use_law_0() {
    guard laws.len > 0 else { return; }
    lame = laws[0].lame;
    lame_3d = laws[0].lame_3d;
    shear = laws[0].shear;
    volume_force = laws[0].force;
}

fn cell_in(c: i64, g: i64) -> bool {
    var p: f64[3] = [0.0, 0.0, 0.0];
    for l in dim + 1 {
        for t in 3 { p[t] += node_x(cell_node(c, l), t) / (dim + 1); }
    }
    return mesher::phi(g, p[0], p[1], p[2]) < 0.0;
}

fn is_live(c: i64) -> bool { return dead.len == 0 || dead[c] == 0; }

// Gives the cells in shape g the elastic law of E, nu and the body force
// (bx, by, bz), with the plane assumption and the thickness of law 0. It
// gives the number of these cells.
fn add_law(g: i64, E: f64, nu: f64, b: f64[3]) -> i64 {
    guard physics == ELASTIC && E > 0.0 && nu > -1.0 && nu < 0.5 else { abort("add_law: an elastic physics, E > 0 and -1 < nu < 0.5"); }
    if laws.len == 0 {
        laws.push(law_now());
        for c in cell_count() { cell_law.push(0); }
    }
    let l3 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    laws.push(Law { lame: if plane == PLANE_STRESS { E * nu / (1.0 - nu * nu) } else { l3 }, lame_3d: l3,
                    shear: E / (2.0 * (1.0 + nu)), force: b });
    var count = 0;
    for c in cell_count() {
        if cell_in(c, g) {
            cell_law[c] = laws.len - 1;
            count++;
        }
    }
    return count;
}

// Kills the live cells in shape g. It gives the number of cells in g.
fn kill_cells(g: i64) -> i64 {
    if dead.len == 0 { for c in cell_count() { dead.push(0); } }
    var count = 0;
    for c in cell_count() {
        if cell_in(c, g) {
            dead[c] = 1;
            count++;
        }
    }
    return count;
}

// Whether a live cell holds each node.
fn live_nodes() -> u8[>..] {
    var live: u8[>..] = [];
    for n in node_count() { live.push(0); }
    for c in cell_count() {
        if !is_live(c) { continue; }
        for l in per_cell { live[cell_node(c, l)] = 1; }
    }
    return live;
}

// Gives birth to the dead cells in shape g, strain-free on the solution of
// the last step. It gives the number of cells in g.
fn birth_cells(g: i64) -> i64 {
    guard dofs_sized() else { abort("birth_cells: the arrays of the dofs"); }
    if dead.len == 0 { for c in cell_count() { dead.push(0); } }
    if born.len == 0 {
        for c in cell_count() { born.push(0); }
        for i in cell_count() * per_cell * ncomp { birth_u.push(0.0); }
    }
    if held.len != sol.len {
        held.clear();
        held.append(sol);
    }
    let before = live_nodes();
    // The nodes of the new cells that a live cell held: the start of the
    // nearest-node search.
    var new_cells: i64[>..] = [];
    var shared: i64[>..] = [];
    var seen: u8[>..] = [];
    for n in node_count() { seen.push(0); }
    var count = 0;
    for c in cell_count() {
        guard cell_in(c, g) else { continue; };
        count++;
        guard dead[c] != 0 else { continue; };
        new_cells.push(c);
        for l in per_cell {
            let n = cell_node(c, l);
            if before[n] != 0 && seen[n] == 0 { shared.push(n); }
            seen[n] = 1;
        }
    }
    for c in new_cells {
        dead[c] = 0;
        born[c] = 1;
        for l in per_cell {
            let n = cell_node(c, l);
            if before[n] != 0 || seen[n] != 1 || shared.len == 0 { continue; }
            var best = shared[0];
            var best_d = 1.0e300;
            for m in shared {
                var d = 0.0;
                for t in 3 { d += (node_x(n, t) - node_x(m, t)) * (node_x(n, t) - node_x(m, t)); }
                if d < best_d {
                    best_d = d;
                    best = m;
                }
            }
            for p in ncomp { held[ncomp * n + p] = held[ncomp * best + p]; }
            seen[n] = 2;
        }
    }
    for c in new_cells {
        for l in per_cell {
            for p in ncomp { birth_u[(per_cell * c + l) * ncomp + p] = held[ncomp * cell_node(c, l) + p]; }
        }
    }
    return count;
}

// Holds each dof of a node that no live cell holds at its value of the
// last step, with a 1 on the diagonal of its row.
fn hold_dead_nodes() {
    guard dead.len > 0 else { return; }
    let live = live_nodes();
    for n in node_count() {
        if live[n] != 0 { continue; }
        for p in ncomp {
            let d = ncomp * n + p;
            fixed[d] = 2;
            prescribed[d] = if held.len == sol.len { held[d] } else { 0.0 };
            val[entry(d, d)] = 1.0;
        }
    }
}
