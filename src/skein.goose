// The `skein` command line (docs/command_line.md). Each command goes to the
// module that owns it. The commands that have no module yet print that they
// are not implemented.
import std;
import os;
import syntax.tree;

let USAGE = """
    usage: skein <command> [options] [files]

    commands:
      check   parse and type-check a model
      fmt     rewrite files in the standard layout
      help    print this text
    """;

fn main() {
    let args = args();
    guard args.len > 1 else {
        write_stderr(str(USAGE, "\n"));
        exit(2);
    }
    let command = args[1];
    if command == "help" || command == "--help" || command == "-h" {
        print(USAGE);
        exit(0);
    }
    if command == "check" || command == "fmt" {
        write_stderr(str("skein: ", command, " is not implemented yet\n"));
        exit(3);
    }
    write_stderr(str("skein: unknown command '", command, "'\n", USAGE, "\n"));
    exit(2);
}
