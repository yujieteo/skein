import std;
namespace codegen;
import syntax from "../syntax";
import check from "../check";
import fmt.fmt;
import .codegen;
import .study;
import .split;
import .slots;

fn web_results_dir(model: u8[:]) -> u8[>..] { return str(".skein/results/", model); }

var web_grid = false;

var web_simplex = false;

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

fn web_dof_estimate() -> i64 { return if web_grid { grid_gen::dof_estimate() } else if web_simplex { elements_gen::dof_estimate() } else { dof_estimate() }; }
