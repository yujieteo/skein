import std;
import os;
import syntax from "./syntax";
import fmt.fmt;
import cli.card;
import codegen.run;
import codegen.build;
import codegen.convergence;
import check from "./check";
import md from "./md";
import report from "./report";
import mesher.verify;
import mesher.lint;
import mesher.stl;
import mesher from "./mesher";
import runtime from "./runtime";
import solvers from "./solvers";

let USAGE = """
    usage: skein <command> [options] [files]

    commands:
      build   build the page of a model that meshes and solves in the browser (--target wasm)
      check   parse and type-check a model
      export  write a let shape as an STL surface on standard output
      fmt     rewrite files in the standard layout
      help    print this text (help std: the standard library, help --llm: the language card)
      lint    check a model, then give the lint warnings
      run     build and run a model (linear static or steady heat: a box, a quarter of a tube, a plate with a hole, a cylindrical shell or beams; --study: each case of a study)
      md      write a model as a Markdown document
      report  run a model, then write one HTML file with its source, the solver tables and a viewer of the results
      verify  check a model, mesh it and report mesh_quality
    """;

fn main() {
    let args = args();
    guard args.len > 1 else {
        write_stderr(str(USAGE, "\n"));
        exit(2);
    }
    let command = args[1];
    if command == "help" && args.len > 2 && args[2] == "--llm" {
        print(cli::CARD);
        exit(0);
    }
    if command == "help" && args.len > 2 && args[2] == "std" {
        print(check::std_text());
        exit(0);
    }
    if command == "help" || command == "--help" || command == "-h" {
        print(USAGE);
        exit(0);
    }
    if command == "fmt" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(fmt::run(rest));
    }
    if command == "run" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(codegen::run(rest));
    }
    if command == "build" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(codegen::build(rest));
    }
    if command == "check" { exit(check::command(args)); }
    if command == "lint" { exit(mesher::lint(args)); }
    if command == "verify" && args.len > 2 && args[2] == "convergence" {
        var rest: const u8[:][>..] = [];
        for a in args[3..] { rest.push(a); }
        exit(codegen::convergence(rest));
    }
    if command == "verify" { exit(mesher::verify(args)); }
    if command == "export" { exit(mesher::export_command(args)); }
    if command == "md" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(md::run(rest));
    }
    if command == "report" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(report::run(rest) { prog => codegen::page_takes(prog) });
    }
    write_stderr(str("skein: unknown command '", command, "'\n", USAGE, "\n"));
    exit(2);
}
