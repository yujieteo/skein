// `skein lint` (docs/command_line.md §4). With `--mesh`, it makes the mesh
// of each `mesh` block, as `skein verify` does, and applies the lint rules
// that need the mesh (docs/checks.md §3): `bending_lock` and
// `mesh_feature_small`. The mesher makes a uniform mesh, so the local mesh
// size is the `size` of the block. The other rules are in src/check/lint.goose.
import std;
namespace mesher;
import syntax from "../syntax";
import check from "../check";
import .geometry;
import .verify;

// `skein lint`. `args[1]` is "lint". The result is the exit status
// (command_line.md §10).
fn lint(args: u8[:][:]) -> i64 {
    let status, json, mesh = check::start(args);
    if status >= 0 { return status; }
    if mesh { mesh_lint(); }
    return check::finish(json);
}

// Applies the rules that need the mesh to the model of `check::lint_text`,
// when it has no error. An error of the mesher goes into `check::diags`.
fn mesh_lint() {
    guard check::error_count() == 0 else { return; };
    let prog = check::model_prog;
    let dim = read_geometry(prog);
    var records: MeshRecord[>..] = [];
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        if syntax::nodes[item].kind == syntax::N_MESH {
            if !mesh_block(item, dim, records) { break; }
            mesh_rules(item, records[records.len - 1], dim);
        }
        item = syntax::nodes[item].next;
    }
    check::take_syntax_errors();
    check::sort_diags();
}

// The smallest feature of shape `g`: the shortest side of the bounding box
// of each box, cylinder, rect, circle, polygon, extrude and revolve in it.
recursive fn smallest_feature(g: i64) -> f64 {
    let s = shapes[g];
    if s.kind == UNION || s.kind == INTERSECT || s.kind == SUBTRACT {
        return min(smallest_feature(s.a), smallest_feature(s.b));
    }
    var f = if s.kind == EXTRUDE || s.kind == REVOLVE { smallest_feature(s.a) } else { 1.0e300 };
    let b = bounds(g);
    for c in 3 {
        if b[3 + c] > b[c] { f = min(f, b[3 + c] - b[c]); }
    }
    return f;
}

// The rules on mesh block `n`, with its record `r`.
//
// `mesh_feature_small`: the smallest feature of the region is under the
// mesh size.
//
// `bending_lock`: an order-1 displacement on the mesh, whose region is thin:
// the shortest side of its bounding box is under 4 cells, and the longest
// side is at least 4 times the shortest, so a load bends it.
fn mesh_rules(n: i64, r: MeshRecord, dim: i64) {
    let unit = length_name;
    let f = smallest_feature(r.g);
    if f < r.size {
        check::warn(syntax::nodes[n].tok, "mesh_feature_small",
                    str("a feature of the region of mesh '", node_text(n), "' is ", check::approx(f), " ", unit,
                        ", under the mesh size ", check::approx(r.size), " ", unit),
                    str("size: ", check::approx(f), " ", unit, " or less"));
    }
    let b = bounds(r.g);
    var thickness = 1.0e300;
    var length = 0.0;
    for c in dim {
        thickness = min(thickness, b[3 + c] - b[c]);
        length = max(length, b[3 + c] - b[c]);
    }
    guard thickness < 4.0 * r.size && length >= 4.0 * thickness else { return; };
    for p, i in check::physics {
        let u = check::displacement_of(i);
        let m = check::method_of(i);
        guard check::mesh_of(i) == n && u >= 0 && check::order_of(m, u) == 1.0 else { continue; };
        let order = check::setting(syntax::nodes[m].first, "order");
        check::warn(check::first_tok(order), "bending_lock",
                    str("the displacement '", check::tok_text(check::syms[u].tok), "' has order-1 ", r.cell,
                        " cells, and the region of mesh '", node_text(n), "' is ", check::approx(thickness), " ", unit,
                        " thick, under 4 cells of ", check::approx(r.size), " ", unit),
                    str("order: 2, or size: ", check::approx(thickness / 4.0), " ", unit, " or less"));
    }
}
