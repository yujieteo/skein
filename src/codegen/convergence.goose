// `skein verify convergence model.skn [--levels n]` (docs/checks.md §4.3,
// command_line.md §5): solve the model on n meshes, each with half the
// size of the one before, and give the observed order of the L2 error of
// each `expect u == g`.
//
//     .skein/build/<model>/convergence/<level>/main.goose    the program of a level
//     .skein/results/<model>/convergence/<level>/            its results
//
// The fv and dg models of src/grid/gen.goose run. The order of an expect
// is log2 of the ratio of its errors on the two finest meshes. It fails
// when it is lower than the expected order of the method by more than 0.3.
// The id is `convergence`.
import std;
import os;
import math;
namespace codegen;
import syntax from "../syntax";
import fmt.fmt;
import grid.gen;
import .run;

let CONVERGENCE_USAGE = "usage: skein verify convergence <file.skn> [--levels n]\n";
let ORDER_MARGIN = 0.3;   // checks.md §4.3

// The L2 errors of an expect on each level.
struct Errors {
    name: u8[..64],
    cells: i64[..8],    // --levels is 8 or less
    error: f64[..8],
}

// Adds the rows of the record l2 of the program output `out` to `errors`.
fn read_errors(out: u8[:], errors: Errors[>..]&) {
    var inside = false;
    for line in split(out, '\n') {
        if starts_with(line, "l2[") {
            inside = true;
            continue;
        }
        if !inside || !starts_with(line, "  ") { inside = false; continue; }
        let fields = split(trim(line), ',');
        guard fields.len == 4 else { continue; };
        let cells, good_cells = parse_int(fields[1]);
        let error, good_error = parse_flt(fields[2]);
        guard good_cells && good_error else { continue; };
        var e = -1;
        for x, i in errors { if x.name == fields[0] { e = i; } }
        if e < 0 {
            errors.push(Errors { name: [], cells: [], error: [] });
            e = errors.len - 1;
            set64(errors[e].name, fields[0]);
        }
        errors[e].cells.push(cells);
        errors[e].error.push(error);
    }
}

// Builds and runs the program of the model on the mesh of `level`, and
// adds its L2 errors. False, with an error on stderr, if it cannot.
fn run_level(path: u8[:], text: u8[:], level: i64, errors: Errors[>..]&) -> bool {
    grid_gen::level = level;
    let prog = syntax::parse_model(path, text);
    let lane = if prog >= 0 { read_model(prog) } else { -1 };
    guard lane >= 0 else {
        fmt::write_errors();
        return false;
    };
    guard lane == LANE_GRID else {
        write_stderr("skein verify convergence: the model has no method fv or dg. Only these run yet.\n");
        return false;
    };
    let build_dir = str(".skein/build/", model_name, "/convergence/", level);
    let results_dir = str(".skein/results/", model_name, "/convergence/", level);
    guard make_dirs(build_dir) && make_dirs(results_dir) else {
        write_stderr("skein verify convergence: cannot make the directory .skein\n");
        return false;
    };
    var program: u8[>..] = [];
    guard write_program(lane, results_dir, no_share(), program) else { return false; };
    let source = str(build_dir, "/main.goose");
    let exe = str(build_dir, "/", model_name);
    let out = str(build_dir, "/out.toon");
    guard write_file(source, program) else {
        write_stderr(str("skein verify convergence: cannot write '", source, "'\n"));
        return false;
    };
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    guard shell(str(shell_quoted(goose), " build -q ", shell_quoted(source), " -o ", shell_quoted(exe), " 1>&2")) == 0 else {
        write_stderr(str("skein verify convergence: internal fault: goose cannot build '", source, "'. Please report it.\n"));
        return false;
    };
    // A failed expect on a coarse mesh is not a fault here: exit status 1.
    let status = shell(str(shell_quoted(str("./", exe)), " >", shell_quoted(out)));
    var text_out: u8[>..] = [];
    guard (status == 0 || status == 1) && read_file(out, text_out) else {
        write_stderr(str("skein verify convergence: the program of level ", level, " stopped with status ", status, "\n"));
        return false;
    };
    for line in split(text_out, '\n') {
        guard !(starts_with(line, "  solve,") && find(line, ",fail,") >= 0) else {
            write_stderr(str("skein verify convergence: the solve of level ", level, " did not converge: ", trim(line), "\n"));
            return false;
        };
    }
    read_errors(text_out, errors);
    return true;
}

// `skein verify convergence`. `args` are the words after `convergence`.
fn convergence(args: const u8[:][:]) -> i64 {
    var levels = 3;
    var at = -1;
    var i = 0;
    while i < args.len {
        if args[i] == "--levels" && i + 1 < args.len {
            let k, good = parse_int(args[i + 1]);
            guard good && k >= 2 && k <= 8 else {
                write_stderr(str("skein verify convergence: --levels takes a whole number from 2 to 8\n", CONVERGENCE_USAGE));
                return 2;
            };
            levels = k;
            i += 2;
        } else if at < 0 && args[i].len > 0 && args[i][0] != '-' {
            at = i;
            i++;
        } else {
            write_stderr(str("skein verify convergence: give one model file\n", CONVERGENCE_USAGE));
            return 2;
        }
    }
    guard at >= 0 else {
        write_stderr(str("skein verify convergence: give one model file\n", CONVERGENCE_USAGE));
        return 2;
    };
    let path = args[at];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein verify convergence: cannot read '", path, "'\n"));
        return 1;
    };
    var errors: Errors[>..] = [];
    for level in levels { guard run_level(path, text, level, errors) else { return 1; }; }
    let want, why = grid_gen::expected_order();
    return print_convergence(errors, levels, want, why);
}

// Prints the table of the errors and the check of each expect, and gives
// the exit status.
fn print_convergence(errors: Errors[:], levels: i64, want: f64, why: u8[:]) -> i64 {
    var table: u8[>..] = [];
    var checks: u8[>..] = [];
    var rows = 0;
    var failed = errors.len == 0;
    for e in errors {
        var order = 0.0;
        for k in e.error.len {
            if k > 0 { order = log(e.error[k - 1] / e.error[k]) / log(2.0); }
            format(table, "  ", e.name, ",", k, ",", e.cells[k], ",", e.error[k], ",", if k > 0 { str(order) } else { str("-") }, "\n");
            rows++;
        }
        let pass = e.error.len == levels && order >= want - ORDER_MARGIN;
        if !pass { failed = true; }
        format(checks, "  post,convergence,", if pass { "pass" } else { "fail" }, ",\"", e.name, ": order ", order,
               "\",\">= ", want - ORDER_MARGIN, " (", why, ")\"\n");
    }
    if errors.len == 0 {
        checks.append("  post,convergence,fail,\"no expect u == g\",\"an expect that compares two fields\"\n");
    }
    var s: u8[>..] = [];
    format(s, "convergence[", rows, "]{expect,level,cells,error,order}:\n", table);
    format(s, "verify[", max(1, errors.len), "]{stage,id,result,value,limit}:\n", checks);
    write_stdout(s);
    return if failed { 1 } else { 0 };
}
