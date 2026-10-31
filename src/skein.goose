// The `skein` command line (docs/command_line.md). Each command goes to the
// module that owns it. The commands that have no module yet print that they
// are not implemented.
import std;
import os;
import syntax from "./syntax";
import fmt.fmt;
import cli.card;
import codegen.run;
import check from "./check";
import md from "./md";
import mesher.verify;   // mesher lane: skein verify
import mesher.lint;     // skein lint, with the rules that need the mesh
import runtime from "./runtime";   // mesher lane: the mesh type of the mesher
import solvers from "./solvers";   // the runtime imports it

let USAGE = """
    usage: skein <command> [options] [files]

    commands:
      check   parse and type-check a model
      fmt     rewrite files in the standard layout
      help    print this text
      lint    check a model, then give the lint warnings
      run     build and run a model (linear static or steady heat: a box, a quarter of a tube, a plate with a hole or a cylindrical shell)
      md      write a model as a Markdown document
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
    if command == "check" { exit(check::command(args)); }
    if command == "lint" { exit(mesher::lint(args)); }
    if command == "verify" { exit(mesher::verify(args)); }   // mesher lane
    if command == "md" {
        var rest: const u8[:][>..] = [];
        for a in args[2..] { rest.push(a); }
        exit(md::run(rest));
    }
    write_stderr(str("skein: unknown command '", command, "'\n", USAGE, "\n"));
    exit(2);
}
