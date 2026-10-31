// The `skein` command line (docs/command_line.md). Each command goes to the
// module that owns it. The commands that have no module yet print that they
// are not implemented.
import std;
import os;
import syntax from "./syntax";
import fmt.fmt;
import cli.card;
import codegen.run;
import codegen.build;
import check from "./check";
import md from "./md";
import report from "./report";   // skein report
import mesher.verify;   // mesher lane: skein verify
import mesher.lint;     // skein lint, with the rules that need the mesh
import mesher.stl;      // skein export
import mesher from "./mesher";   // skein-flow: skein run meshes a flow model
import runtime from "./runtime";   // mesher lane: the mesh type of the mesher
import solvers from "./solvers";   // the runtime imports it

let USAGE = """
    usage: skein <command> [options] [files]

    commands:
      build   build the page of a model that meshes and solves in the browser (--target wasm)
      check   parse and type-check a model
      export  write a let shape as an STL surface on standard output
      fmt     rewrite files in the standard layout
      help    print this text (help std: the standard library, help --llm: the language card)
      lint    check a model, then give the lint warnings
      run     build and run a model (linear static or steady heat: a box, a quarter of a tube, a plate with a hole or a cylindrical shell; --study: each case of a study)
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
    if command == "verify" { exit(mesher::verify(args)); }   // mesher lane
    if command == "export" { exit(mesher::export_command(args)); }
    if command == "md" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(md::run(rest));
    }
    if command == "report" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(report::run(rest));
    }
    write_stderr(str("skein: unknown command '", command, "'\n", USAGE, "\n"));
    exit(2);
}
