// The stages of the work of a case that a param reaches (spec §20.4). The
// study runner reuses the stages that no varied param reaches. A later
// `skein check --study` can call it to report the stages of each param.
import std;
namespace check;
import syntax from "../syntax";

let STAGE_POST = 0;    // the expect and report expressions
let STAGE_SOLVE = 1;   // the assembly and the solve: a material, a physics, a bc or a load
let STAGE_MESH = 2;    // the geometry and the mesh, and each stage after them

// Whether the tree `n` uses one of the names `names`.
recursive fn uses_name(n: i64, names: u8[..64][:]) -> bool {
    let k = syntax::nodes[n];
    if k.kind == syntax::N_NAME && position(names, syntax::text(k.tok)) >= 0 { return true; }
    var c = k.first;
    while c >= 0 {
        if uses_name(c, names) { return true; }
        c = syntax::nodes[c].next;
    }
    return false;
}

fn declares_value(kind: i64) -> bool {
    return kind == syntax::N_LET || kind == syntax::N_PARAM || kind == syntax::N_FN || kind == syntax::N_TABLE;
}

// The stage of an item other than a step that uses a param.
fn item_stage(kind: i64) -> i64 {
    if kind == syntax::N_EXPECT { return STAGE_POST; }
    if kind == syntax::N_MATERIAL || kind == syntax::N_PHYSICS || kind == syntax::N_TEMPLATE || kind == syntax::N_BC ||
       kind == syntax::N_MATASSIGN || kind == syntax::N_CONTACT || kind == syntax::N_INTERFACE {
        return STAGE_SOLVE;
    }
    return STAGE_MESH;
}

// The first stage that the param `name` of program `prog` reaches, or -1
// when no item uses it. The param reaches each `let`, `fn` or `table` whose
// value uses it or a declaration that it reaches.
fn param_stage(prog: i64, name: u8[:]) -> i64 {
    var reached: u8[..64][>..] = [];
    var w: u8[..64] = [];
    w.append(name[..min(name.len, 64)]);
    reached.push(w);
    var grown = true;
    while grown {
        grown = false;
        var c = syntax::nodes[prog].first;
        while c >= 0 {
            if declares_value(syntax::nodes[c].kind) {
                let d = syntax::text(syntax::nodes[c].tok);
                if position(reached, d) < 0 && uses_name(c, reached) {
                    var v: u8[..64] = [];
                    v.append(d[..min(d.len, 64)]);
                    reached.push(v);
                    grown = true;
                }
            }
            c = syntax::nodes[c].next;
        }
    }
    var stage = -1;
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        let kind = syntax::nodes[c].kind;
        if kind == syntax::N_STEP {
            var s = syntax::nodes[c].first;
            while s >= 0 {
                if uses_name(s, reached) { stage = max(stage, if syntax::nodes[s].kind == syntax::N_REPORT { STAGE_POST } else { STAGE_SOLVE }); }
                s = syntax::nodes[s].next;
            }
        } else if !declares_value(kind) && kind != syntax::N_STUDY && uses_name(c, reached) {
            stage = max(stage, item_stage(kind));
        }
        c = syntax::nodes[c].next;
    }
    return stage;
}
