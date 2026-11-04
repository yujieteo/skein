import std;
import math;
namespace elements;
import mesher from "../mesher";
import .mesh;
import .assemble;

struct Law {
    cond: f64,
    heat: f64,
    lame: f64,
    lame_3d: f64,
    shear: f64,
    force: f64[3],
    ani: bool,
    cm: f64[36],
    s0: f64[6],
}

var laws: Law[>..] = [];
var cell_law: i64[>..] = [];
var dead: u8[>..] = [];
var born: u8[>..] = [];
var birth_u: f64[>..] = [];
var held: f64[>..] = [];
var origin: f64[>..] = [];

fn clear_stages() {
    laws.clear();
    cell_law.clear();
    dead.clear();
    born.clear();
    birth_u.clear();
    held.clear();
    origin.clear();
}

fn shown(d: i64) -> f64 {
    return if origin.len == sol.len { sol[d] - origin[d] } else { sol[d] };
}

fn law_now() -> Law { return Law { cond: conductivity, heat: source, lame: lame, lame_3d: lame_3d, shear: shear, force: volume_force,
                                          ani: aniso_on, cm: aniso_c, s0: aniso_s0 }; }

fn use_law(c: i64) {
    guard cell_law.len > 0 else { return; }
    let w = laws[cell_law[c]];
    conductivity = w.cond;
    source = w.heat;
    lame = w.lame;
    lame_3d = w.lame_3d;
    shear = w.shear;
    volume_force = w.force;
    use_aniso(w);
}

fn use_aniso(w: Law) {
    if w.ani { set_aniso(w.cm, w.s0); } else { aniso_on = false; }
}

fn use_law_0() {
    guard laws.len > 0 else { return; }
    conductivity = laws[0].cond;
    source = laws[0].heat;
    lame = laws[0].lame;
    lame_3d = laws[0].lame_3d;
    shear = laws[0].shear;
    volume_force = laws[0].force;
    use_aniso(laws[0]);
}

fn cell_in(c: i64, g: i64) -> bool {
    var p: f64[3] = [0.0, 0.0, 0.0];
    let corners = if hex_cells { 8 } else { dim + 1 };
    for l in corners {
        for t in 3 { p[t] += node_x(cell_node(c, l), t) / corners; }
    }
    return mesher::phi(g, p[0], p[1], p[2]) < 0.0;
}

fn is_live(c: i64) -> bool { return dead.len == 0 || dead[c] == 0; }

fn add_law(g: i64, E: f64, nu: f64, b: f64[3]) -> i64 {
    guard physics == ELASTIC && E > 0.0 && nu > -1.0 && nu < 0.5 else { abort("add_law: an elastic physics, E > 0 and -1 < nu < 0.5"); }
    if laws.len == 0 {
        laws.push(law_now());
        for c in cell_count() { cell_law.push(0); }
    }
    let l3 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    laws.push(Law { cond: conductivity, heat: source, lame: if plane == PLANE_STRESS { E * nu / (1.0 - nu * nu) } else { l3 }, lame_3d: l3,
                    shear: E / (2.0 * (1.0 + nu)), force: b, ani: false, cm: [0.0; 36], s0: [0.0; 6] });
    var count = 0;
    for c in cell_count() {
        if cell_in(c, g) {
            cell_law[c] = laws.len - 1;
            count++;
        }
    }
    return count;
}

fn add_aniso_law(g: i64, c: f64[36], s0: f64[6], b: f64[3]) -> i64 {
    guard physics == ELASTIC else { abort("add_aniso_law: an elastic physics"); }
    if laws.len == 0 {
        laws.push(law_now());
        for k in cell_count() { cell_law.push(0); }
    }
    laws.push(Law { cond: conductivity, heat: source, lame: lame, lame_3d: lame_3d, shear: shear, force: b, ani: true, cm: c, s0: s0 });
    var count = 0;
    for k in cell_count() {
        if cell_in(k, g) {
            cell_law[k] = laws.len - 1;
            count++;
        }
    }
    return count;
}

fn add_heat_law(g: i64, k: f64, q: f64) -> i64 {
    guard physics == HEAT && k > 0.0 && curve_k.len == 0 else { abort("add_heat_law: a heat physics with a constant conductivity, and k > 0"); }
    if laws.len == 0 {
        laws.push(law_now());
        for c in cell_count() { cell_law.push(0); }
    }
    laws.push(Law { cond: k, heat: q, lame: lame, lame_3d: lame_3d, shear: shear, force: volume_force, ani: false, cm: [0.0; 36], s0: [0.0; 6] });
    var count = 0;
    for c in cell_count() {
        if cell_in(c, g) {
            cell_law[c] = laws.len - 1;
            count++;
        }
    }
    return count;
}

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

fn live_nodes() -> u8[>..] {
    var live: u8[>..] = [];
    for n in node_count() { live.push(0); }
    for c in cell_count() {
        if !is_live(c) { continue; }
        for l in per_cell { live[cell_node(c, l)] = 1; }
    }
    return live;
}

fn adopt_values(new_cells: i64[>..]&, shared: i64[>..]&, before: u8[>..]&, seen: u8[>..]&) {
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
}

fn record_birth(new_cells: i64[>..]&, before: u8[>..]&) {
    if origin.len != sol.len {
        origin.clear();
        for d in sol.len { origin.push(0.0); }
    }
    for c in new_cells {
        for l in per_cell {
            let n = cell_node(c, l);
            if before[n] == 0 { for p in ncomp { origin[ncomp * n + p] = held[ncomp * n + p]; } }
        }
    }
    for c in new_cells {
        for l in per_cell {
            for p in ncomp { birth_u[(per_cell * c + l) * ncomp + p] = held[ncomp * cell_node(c, l) + p]; }
        }
    }
}

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
    adopt_values(new_cells, shared, before, seen);
    record_birth(new_cells, before);
    return count;
}

fn hold_dead_nodes() {
    guard dead.len > 0 else { return; }
    let live = live_nodes();
    for n in node_count() {
        if live[n] != 0 { continue; }
        for p in ncomp {
            let d = ncomp * n + p;
            fixed[d] = FIXED_HELD;
            prescribed[d] = if held.len == sol.len { held[d] } else { 0.0 };
            val[entry(d, d)] = 1.0;
        }
    }
}
