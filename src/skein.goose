// The `skein` command line (docs/command_line.md). Each command goes to the
// module that owns it. The commands that have no module yet print that they
// are not implemented.
import std;
import os;
import syntax from "./syntax";
import fmt.fmt;
import cli.card;
import codegen.run;
import check.check;

let USAGE = """
    usage: skein <command> [options] [files]

    commands:
      check   parse and type-check a model
      fmt     rewrite files in the standard layout
      help    print this text
      run     build and run a model (linear static or steady heat, on a box or a plate with a hole)
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
    write_stderr(str("skein: unknown command '", command, "'\n", USAGE, "\n"));
    exit(2);
}
