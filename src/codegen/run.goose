// `skein run` (docs/command_line.md §5) for the first slice: parse the
// model, generate a Goose program, build it with goose to a native
// executable and run it. The program prints the reports in TOON and writes
// the VTU files.
//
//     .skein/build/<model>/main.goose    the generated program
//     .skein/build/<model>/<model>       the executable
//     .skein/results/<model>/<step>.vtu  the results of each step
import std;
import os;
namespace codegen;
import syntax from "../syntax";
import fmt.fmt;
import .codegen;
import elements.gen;   // skein-elements: the simplex models
import shells.gen;     // skein-shells: the shell models

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

let USAGE = "usage: skein run <file.skn>\n";

// Runs `command` in the shell, and gives its exit status.
fn shell(command: u8[:]) -> i64 {
    var c: u8[>..] = [];
    c.append(command);
    c.push(0);
    let status = system(c[0]);
    if status < 0 || (status & 0x7f) != 0 { return 3; }
    return (status >> 8) & 0xff;
}

// `s` in single quotes for the shell.
fn shell_quoted(s: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    out.push('\'');
    for c in s {
        if c == '\'' { out.append("'\\''"); } else { out.push(c); }
    }
    out.push('\'');
    return out;
}

fn make_dirs(path: u8[:]) -> bool {
    for c, i in path {
        if c == '/' && i > 0 { make_dir(path[..i]); }
    }
    return make_dir(path);
}

// `skein run`. `args` are the words after `run`. It gives the exit status
// (command_line.md §10).
fn run(args: const u8[:][:]) -> i64 {
    guard args.len == 1 && args[0].len > 0 && args[0][0] != '-' else {
        write_stderr(str("skein run: give one model file\n", USAGE));
        return 2;
    }
    let path = args[0];
    var text: u8[>..] = [];
    if !read_file(path, text) {
        write_stderr(str("skein run: cannot read '", path, "'\n"));
        return 1;
    }
    let prog = fmt::parse_text(path, text);
    // skein-elements: a simplex model goes to src/elements/gen.goose.
    // skein-shells: a shell model goes to src/shells/gen.goose.
    let shell_model = prog >= 0 && shells_gen::claims(prog);
    let simplex = prog >= 0 && !shell_model && elements_gen::claims(prog);
    if prog < 0 || !(if shell_model { shells_gen::read(prog) } else if simplex { elements_gen::read(prog) } else { read(prog) }) {
        fmt::write_errors(path);
        return 1;
    }
    let build_dir = str(".skein/build/", model_name);
    let results_dir = str(".skein/results/", model_name);
    if !make_dirs(build_dir) || !make_dirs(results_dir) {
        write_stderr("skein run: cannot make the directory .skein\n");
        return 1;
    }
    // The skein executable is build/skein in the repository, so the runtime
    // is in src/runtime beside it.
    let runtime_path = str(resource_dir(), if shell_model { "../src/shells" } else if simplex { "../src/elements" } else { "../src/runtime" });
    var program: u8[>..] = [];
    if shell_model {
        shells_gen::generate(runtime_path, results_dir, program);
    } else if simplex {
        elements_gen::generate(runtime_path, results_dir, program);
    } else {
        generate(runtime_path, results_dir, program);
    }
    let source = str(build_dir, "/main.goose");
    let exe = str(build_dir, "/", model_name);
    if !write_file(source, program) {
        write_stderr(str("skein run: cannot write '", source, "'\n"));
        return 1;
    }
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    // goose writes its progress on stdout, and stdout is for the TOON.
    let build = str(shell_quoted(goose), " build -q ", shell_quoted(source), " -o ", shell_quoted(exe), " 1>&2");
    if shell(build) != 0 {
        write_stderr(str("skein run: internal fault: goose cannot build '", source, "'. Please report it.\n"));
        return 3;
    }
    return shell(shell_quoted(str("./", exe)));
}
