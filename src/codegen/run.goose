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
import check from "../check";
import .study;

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

let USAGE = "usage: skein run <file.skn> [--study <name> [--jobs k] [--format toon|json]]\n";

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

// The generators of `read_model`.
let LANE_BOX = 0;
let LANE_SIMPLEX = 1;
let LANE_SHELL = 2;

// Reads program `prog` with the generator that claims it. It gives that
// generator, or -1 on an error, and `syntax::errors` holds the errors.
// skein-elements: a simplex model goes to src/elements/gen.goose.
// skein-shells: a shell model goes to src/shells/gen.goose.
fn read_model(prog: i64) -> i64 {
    if shells_gen::claims(prog) { return if shells_gen::read(prog) { LANE_SHELL } else { -1 }; }
    if elements_gen::claims(prog) { return if elements_gen::read(prog) { LANE_SIMPLEX } else { -1 }; }
    return if read(prog) { LANE_BOX } else { -1 };
}

// Writes the program of the model that `read_model` read with generator
// `lane`. The program puts its results in `results_dir`.
fn write_program(lane: i64, results_dir: u8[:], program: u8[>..]&) {
    // The skein executable is build/skein in the repository, so the runtime
    // is in src/runtime beside it.
    let runtime_path = str(resource_dir(), if lane == LANE_SHELL { "../src/shells" } else if lane == LANE_SIMPLEX { "../src/elements" } else { "../src/runtime" });
    if lane == LANE_SHELL {
        shells_gen::generate(runtime_path, results_dir, program);
    } else if lane == LANE_SIMPLEX {
        elements_gen::generate(runtime_path, results_dir, program);
    } else {
        generate(runtime_path, results_dir, program);
    }
}

// `skein run`. `args` are the words after `run`. It gives the exit status
// (command_line.md §10).
fn run(args: const u8[:][:]) -> i64 {
    // The indices in `args` of the model file and of the study name.
    var at = -1;
    var study_at = -1;
    var jobs = 0;
    var json = false;
    var i = 0;
    while i < args.len {
        let a = args[i];
        let more = i + 1 < args.len;
        if a == "--study" && more {
            study_at = i + 1;
        } else if a == "--jobs" && more {
            let k, good = parse_int(args[i + 1]);
            guard good && k >= 1 else {
                write_stderr(str("skein run: --jobs takes a whole number of 1 or more\n", USAGE));
                return 2;
            };
            jobs = k;
        } else if a == "--format" && more && (args[i + 1] == "toon" || args[i + 1] == "json") {
            json = args[i + 1] == "json";
        } else if a.len > 0 && a[0] != '-' && at < 0 {
            at = i;
            i++;
            continue;
        } else {
            write_stderr(str("skein run: give one model file\n", USAGE));
            return 2;
        }
        i += 2;
    }
    guard at >= 0 && (study_at >= 0 || (jobs == 0 && !json)) else {
        write_stderr(str("skein run: give one model file, and --jobs and --format only with --study\n", USAGE));
        return 2;
    };
    let path = args[at];
    var text: u8[>..] = [];
    if !read_file(path, text) {
        write_stderr(str("skein run: cannot read '", path, "'\n"));
        return 1;
    }
    if study_at >= 0 { return run_study(path, text, args[study_at], jobs, json); }
    let prog = syntax::parse_model(path, text);
    let lane = if prog >= 0 { read_model(prog) } else { -1 };
    if lane < 0 {
        fmt::write_errors();
        return 1;
    }
    let build_dir = str(".skein/build/", model_name);
    let results_dir = str(".skein/results/", model_name);
    if !make_dirs(build_dir) || !make_dirs(results_dir) {
        write_stderr("skein run: cannot make the directory .skein\n");
        return 1;
    }
    var program: u8[>..] = [];
    write_program(lane, results_dir, program);
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

// `skein run --study` (spec §20, command_line.md §5, §6.5): the study
// runner, with cases first (spec §20.5). The model is checked once. The
// generator writes one program for each case (study.goose), and contained
// Goose workers build and run the cases, `--jobs` at a time. The collector
// writes one table of the cases, in case order:
//
//     .skein/build/<model>/<study>/<case>/main.goose    the program of a case
//     .skein/studies/<model>/<study>/<case>/             the results of a case
//     .skein/studies/<model>/<study>.toon                the table
//
// Phase 2 meshes and assembles each case again (spec §20.4 gives the reuse).

// A case record. The scheduler gives it to the worker as its argument, so
// that the abort of a worker names its case.
struct CaseJob {
    case: i64,
    build: u8[..4096],   // the shell command that builds the program
    run: u8[..4096],     // the shell command that runs it
    out: u8[..1024],     // the file of its standard output
    err: u8[..1024],     // the file of its standard error
}

// A result record.
struct CaseResult {
    case: i64,
    built: bool,
    status: i64,
    seconds: f64,
    err: u8[..1024],   // the start of its standard error, or the abort message
    out: u8[>..],      // its standard output
}

thread_fn case_worker(job: CaseJob) {
    var r = CaseResult { case: job.case, built: shell(job.build) == 0, status: 0, seconds: 0.0, out: [], err: [] };
    if r.built {
        let start = clock();
        r.status = shell(job.run);
        r.seconds = clock() - start;
        read_file(job.out, r.out);
    }
    var err: u8[>..] = [];
    read_file(job.err, err);
    r.err.append(err[..min(err.len, 1024)]);
    qput(r);
}

// Runs `jobs` with `k` contained workers at a time, and gives the end of
// each case, in case order.
fn run_cases(jobs: CaseJob[:], k: i64) -> CaseEnd[>..] {
    outs.clear();
    var ends: CaseEnd[>..] = [];
    var ids: i64[>..] = [];
    var next = 0;
    while next < jobs.len && next < k {
        ids.push(thread_spawn_contained(case_worker, jobs[next]));
        next++;
    }
    var ended = 0;
    for i in 0..jobs.len {
        var e = CaseEnd { ended: true, built: false, status: 0, seconds: 0.0, err: [], out_lo: 0, out_hi: 0 };
        var why: u8[>..] = [];
        e.ended = thread_wait(ids[i], why);
        if e.ended { ended++; } else { set1024(e.err, str("the worker stopped: ", why)); }
        ends.push(e);
        if next < jobs.len {
            ids.push(thread_spawn_contained(case_worker, jobs[next]));
            next++;
        }
    }
    for i in 0..ended {
        let r = qget<CaseResult>();
        var e = CaseEnd { ended: true, built: r.built, status: r.status, seconds: r.seconds, err: r.err, out_lo: outs.len, out_hi: 0 };
        outs.append(r.out);
        e.out_hi = outs.len;
        ends[r.case] = e;
    }
    return ends;
}

// `skein run <path> --study <study>`, with the model `text`. `jobs` is the
// `--jobs` value, or 0 for the default. The result is the exit status.
fn run_study(path: u8[:], text: u8[:], study: u8[:], jobs: i64, json: bool) -> i64 {
    check::check_text(path, text);
    if check::error_count() > 0 { return check::finish(json); }
    let prog = syntax::parse_model(path, text);
    guard prog >= 0 else {
        fmt::write_errors();
        return 1;
    };
    let s = find_study(prog, study);
    guard s >= 0 else {
        write_stderr(str("skein run: the model has no study '", study, "'\n"));
        return 2;
    };
    if !study_cases(prog, s) {
        fmt::write_errors();
        return 1;
    }
    // The program of each case.
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    var case_jobs: CaseJob[>..] = [];
    var model: u8[>..] = [];
    for c in 0..case_count {
        let cp = syntax::parse_model(path, case_source(text, c));
        let lane = if cp >= 0 { read_model(cp) } else { -1 };
        guard lane >= 0 else {
            fmt::write_errors();
            return 1;
        };
        model.clear();
        model.append(model_name);
        let build_dir = str(".skein/build/", model, "/", study, "/", c + 1);
        let case_dir = str(".skein/studies/", model, "/", study, "/", c + 1);
        if !make_dirs(build_dir) || !make_dirs(case_dir) {
            write_stderr("skein run: cannot make the directory .skein\n");
            return 1;
        }
        var program: u8[>..] = [];
        write_program(lane, case_dir, program);
        let source = str(build_dir, "/main.goose");
        let exe = str(build_dir, "/", model);
        if !write_file(source, program) {
            write_stderr(str("skein run: cannot write '", source, "'\n"));
            return 1;
        }
        let out = str(case_dir, "/run.toon");
        let err = str(case_dir, "/run.err");
        var job = CaseJob { case: c, build: [], run: [], out: [], err: [] };
        format(job.build, shell_quoted(goose), " build -q ", shell_quoted(source), " -o ", shell_quoted(exe), " >", shell_quoted(err), " 2>&1");
        format(job.run, shell_quoted(str("./", exe)), " >", shell_quoted(out), " 2>", shell_quoted(err));
        set1024(job.out, out);
        set1024(job.err, err);
        case_jobs.push(job);
    }
    // The default `k` is min(H, number of cases) (spec §20.5).
    let k = min(if jobs > 0 { jobs } else { max(1, hardware_threads()) }, case_count);
    let ends = run_cases(case_jobs, k);
    // The tokens of the cases replaced those of the model: read it again
    // for the columns and the place of the study.
    let p = syntax::parse_model(path, text);
    var printed: u8[>..] = [];
    var table: u8[>..] = [];
    let status = study_table(study, tok(find_study(p, study)), k, p, ends, json, printed, table);
    let file = str(".skein/studies/", model, "/", study, ".toon");
    if !write_file(file, table) {
        write_stderr(str("skein run: cannot write '", file, "'\n"));
        return 1;
    }
    write_stdout(printed);
    return status;
}
