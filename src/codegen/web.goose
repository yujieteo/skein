// The browser solve (implementation.md §8): the part of `skein build
// --target wasm` that also runs in the page.
//
// The program of a model holds its values as number literals, so a new
// value of a param needs a new program. The page has no Goose compiler.
// Thus `skein build` generates the program at the default values and at
// the ends of the range of each param, and the literals that change
// become slots: `web_f[k]` for an f64 and `web_i[k]` for an i64. The
// WebAssembly module holds the parser, the checker and the generator too.
// For the values of the page, it generates the program again, takes the
// literals of the slots from it, and runs the program with them. A text
// that differs in another place than a slot stops the run: the module
// cannot run that program.
import std;
namespace codegen;
import syntax from "../syntax";
import check from "../check";
import fmt.fmt;
// elements_gen and grid_gen come from the importer: src/codegen/run.goose,
// or the program of the page, which imports elements.gen and grid.gen
// before the packages elements and grid (page.goose).
import .codegen;
import .study;
import .split;
import .slots;

// The results of a run in the page go here: the program of the page and
// the programs that `skein build` compares must have the same strings.
fn web_results_dir(model: u8[:]) -> u8[>..] { return str(".skein/results/", model); }

var web_grid = false;   // the model of the last web_program is an fv or dg model

// Whether the last web_program used the simplex generator of
// src/elements/gen.goose, and not the hex generator of codegen.goose.
var web_simplex = false;

// The program of the model `text` of `path` with the values `sets`
// (name=value, as `skein run --set` takes them), in `program`. Only the
// simplex models of src/elements/gen.goose, the fv and dg models of
// src/grid/gen.goose and the hex models of codegen.goose. False, with the
// errors on stderr, if the model does not check or the generator does not
// take it.
fn web_program(path: u8[:], text: u8[:], sets: const u8[:][:], program: u8[>..]&) -> bool {
    var set_text: u8[>..] = [];
    guard set_source(path, text, sets, set_text) else { return false; };
    check::check_text(path, set_text);
    if check::error_count() > 0 {
        check::finish(false);
        return false;
    }
    let prog = syntax::parse_model(path, set_text);
    web_grid = prog >= 0 && grid_gen::claims(prog);
    web_simplex = prog >= 0 && !web_grid && elements_gen::claims(prog);
    guard prog >= 0 && (if web_grid { grid_gen::read(prog) } else if web_simplex { elements_gen::read(prog) } else { read(prog) }) else {
        fmt::write_errors();
        return false;
    };
    if web_grid && grid_gen::method != grid_gen::FV && grid_gen::method != grid_gen::DG {
        write_stderr("skein build: --target wasm takes the fv and dg models of the grid generator, not rkpm and gfd\n");
        return false;
    }
    program.clear();
    if web_grid {
        grid_gen::generate("./grid", web_results_dir(model_name), program);
    } else if web_simplex {
        elements_gen::generate("./elements", web_results_dir(model_name), no_share(), program);
    } else {
        generate("./runtime", web_results_dir(model_name), no_share(), program);
    }
    split_long_functions(program);
    return true;
}

// The dofs of the model of the last web_program.
fn web_dof_estimate() -> i64 { return if web_grid { grid_gen::dof_estimate() } else if web_simplex { elements_gen::dof_estimate() } else { dof_estimate() }; }
