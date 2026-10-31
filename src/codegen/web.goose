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
// elements_gen comes from the importer: src/codegen/run.goose, or the
// program of the page, which imports elements.gen before the package
// elements (page.goose).
import .codegen;
import .study;
import .split;
import .slots;

// The results of a run in the page go here: the program of the page and
// the programs that `skein build` compares must have the same strings.
fn web_results_dir(model: u8[:]) -> u8[>..] { return str(".skein/results/", model); }

// The program of the model `text` of `path` with the values `sets`
// (name=value, as `skein run --set` takes them), in `program`. Only the
// simplex models of src/elements/gen.goose. False, with the errors on
// stderr, if the model does not check or the generator does not take it.
fn web_program(path: u8[:], text: u8[:], sets: const u8[:][:], program: u8[>..]&) -> bool {
    var set_text: u8[>..] = [];
    guard set_source(path, text, sets, set_text) else { return false; };
    check::check_text(path, set_text);
    if check::error_count() > 0 {
        check::finish(false);
        return false;
    }
    let prog = syntax::parse_model(path, set_text);
    if prog >= 0 && !elements_gen::claims(prog) {
        write_stderr("skein build: --target wasm takes only the models of the simplex generator (tri or tet cells)\n");
        return false;
    }
    guard prog >= 0 && elements_gen::read(prog) else {
        fmt::write_errors();
        return false;
    };
    program.clear();
    elements_gen::generate("./elements", web_results_dir(model_name), no_share(), program);
    split_long_functions(program);
    return true;
}
