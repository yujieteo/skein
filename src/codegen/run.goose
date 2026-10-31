// `skein run` (docs/command_line.md §5) for the first slice: parse the
// model, generate a Goose program, build it with goose to a native
// executable and run it. The program prints the reports in TOON and writes
// the VTU files.
//
//     .skein/build/<model>/main.goose    the generated program
//     .skein/build/<model>/<model>       the executable
//     .skein/results/<model>/<step>.vtu  the results of each step
//     .skein/build/<model>/set/<words>/  the program of a run with --set
//     .skein/results/<model>/set/<words>/ the results of a run with --set
import std;
import os;
namespace codegen;
import syntax from "../syntax";
import fmt.fmt;
import .codegen;
import elements.gen;   // skein-elements: the simplex models
import shells.gen;     // skein-shells: the shell models
import flow.gen;       // skein-flow: the flow models
import check from "../check";
import .study;
import .split;
import .status;

import solvers from "../solvers";   // the memory budget of a study

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

let USAGE = "usage: skein run <file.skn> [--set name=value]... [--study <name> [--jobs k] [--format toon|json]]\n";

// Runs `command` in the shell, and gives its exit status.
fn shell(command: u8[:]) -> i64 {
    var c: u8[>..] = [];
    c.append(command);
    c.push(0);
    let status = system(c[0]);
    if status < 0 || (status & 0x7f) != 0 { return 3; }
    return (status >> 8) & 0xff;
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
let LANE_FLOW = 3;

// Reads program `prog` with the generator that claims it. It gives that
// generator, or -1 on an error, and `syntax::errors` holds the errors.
// skein-elements: a simplex model goes to src/elements/gen.goose.
// skein-shells: a shell model goes to src/shells/gen.goose.
// skein-flow: a flow model goes to src/flow/gen.goose.
fn read_model(prog: i64) -> i64 {
    if shells_gen::claims(prog) { return if shells_gen::read(prog) { LANE_SHELL } else { -1 }; }
    if flow_gen::claims(prog) { return if flow_gen::read(prog) { LANE_FLOW } else { -1 }; }
    if elements_gen::claims(prog) { return if elements_gen::read(prog) { LANE_SIMPLEX } else { -1 }; }
    return if read(prog) { LANE_BOX } else { -1 };
}

// The dofs of the mesh of the model that `read_model` read with generator
// `lane`.
fn dof_estimate(lane: i64) -> i64 {
    if lane == LANE_SHELL { return shells_gen::dof_estimate(); }
    if lane == LANE_FLOW { return flow_gen::dof_estimate(); }
    if lane == LANE_SIMPLEX { return elements_gen::dof_estimate(); }
    return codegen::dof_estimate();
}

// The steps of the model that `read_model` read with generator `lane`.
// Each step solves.
fn step_count(lane: i64) -> i64 {
    if lane == LANE_SHELL { return shells_gen::steps.len; }
    if lane == LANE_FLOW { return flow_gen::steps.len; }
    if lane == LANE_SIMPLEX { return elements_gen::steps.len; }
    return steps.len;
}

// Writes the program of the model that `read_model` read with generator
// `lane`. The program puts its results in `results_dir`. In a study, it
// writes or reads the stages of `share`. The shell and flow generators
// share none. For a flow model, it also meshes the domain and writes the
// mesh in `results_dir` (src/flow/gen.goose). False, with an error on
// stderr, if it cannot.
fn write_program(lane: i64, results_dir: u8[:], share: Share, program: u8[>..]&) -> bool {
    // The skein executable is build/skein in the repository, so the runtime
    // is in src/runtime beside it.
    let runtime_path = str(resource_dir(), if lane == LANE_SHELL { "../src/shells" } else if lane == LANE_SIMPLEX { "../src/elements" }
                                           else if lane == LANE_FLOW { "../src/flow" } else { "../src/runtime" });
    if lane == LANE_FLOW {
        guard flow_gen::write_mesh(results_dir) else { return false; }
        flow_gen::generate(runtime_path, results_dir, program);
    } else if lane == LANE_SHELL {
        shells_gen::generate(runtime_path, results_dir, program);
    } else if lane == LANE_SIMPLEX {
        elements_gen::generate(runtime_path, results_dir, share, program);
    } else {
        generate(runtime_path, results_dir, share, program);
    }
    split_long_functions(program);
    let at = unchecked_status(program);
    if at > 0 {
        write_stderr(str("skein run: internal fault: line ", at, " of the program drops the status of a solver call. Please report it.\n"));
        exit(3);
    }
    return true;
}

// `skein run`. `args` are the words after `run`. It gives the exit status
// (command_line.md §10).
fn run(args: const u8[:][:]) -> i64 {
    // The indices in `args` of the model file and of the study name.
    var at = -1;
    var study_at = -1;
    var jobs = 0;
    var json = false;
    var sets: const u8[:][>..] = [];
    var i = 0;
    while i < args.len {
        let a = args[i];
        let more = i + 1 < args.len;
        if a == "--study" && more {
            study_at = i + 1;
        } else if a == "--set" && more {
            sets.push(args[i + 1]);
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
    if sets.len > 0 {
        var set_text: u8[>..] = [];
        guard set_source(path, text, sets, set_text) else { return 2; };
        text.clear();
        text.append(set_text);
        check::check_text(path, text);
        if check::error_count() > 0 { return check::finish(false); }
    }
    if study_at >= 0 { return run_study(path, text, args[study_at], jobs, json); }
    let prog = syntax::parse_model(path, text);
    let lane = if prog >= 0 { read_model(prog) } else { -1 };
    if lane < 0 {
        fmt::write_errors();
        return 1;
    }
    let build_dir = set_path(".skein/build/", model_name, sets);
    let results_dir = results_path(model_name, sets);
    if !make_dirs(build_dir) || !make_dirs(results_dir) {
        write_stderr("skein run: cannot make the directory .skein\n");
        return 1;
    }
    var program: u8[>..] = [];
    if !write_program(lane, results_dir, no_share(), program) { return 1; }
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
// runner (spec §20.5). The model is checked once. The generator writes
// one program for each case (study.goose), and contained Goose workers
// build and run the cases, at most `--jobs` at a time, each when the
// memory ledger admits it (study.goose).
// The collector writes one table of the cases, in case order:
//
//     .skein/build/<model>/<study>/<case>/main.goose    the program of a case
//     .skein/studies/<model>/<study>/<case>/             the results of a case
//     .skein/studies/<model>/<study>/shared/             the shared stages
//     .skein/studies/<model>/<study>.toon                the table
//
// Case 1 runs first and writes the stages that no varied param reaches
// (spec §20.4), and each other case reads them. When the varied params
// reach the physics or the loads, the mesh is shared. When they reach only
// expect and report expressions, the mesh and the solution of each step
// are shared, and each other case runs only the post-processing. The shell
// generator shares no stage, because its mesh holds the section of each
// cell.

// A case record. The scheduler gives it to the worker as its argument, so
// that the abort of a worker names its case.
struct CaseJob {
    case: i64,
    threads: i64,        // the threads of the case (spec §20.5)
    budget: i64,         // the memory budget of the study in bytes, or 0 (SKEIN_MEMORY_BUDGET)
    build: u8[..4096],   // the shell command that builds the program
    exe: u8[..1024],     // the program, quoted for the shell
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

// The shell command that runs the program of `job` with the environment
// `env` and the threads and the budget of the job, its standard output
// to `out`.
fn case_command(job: CaseJob, env: u8[:], out: u8[:]) -> u8[>..] {
    var s = str(env, "SKEIN_THREADS=", job.threads, " ");
    if job.budget > 0 { format(s, "SKEIN_MEMORY_BUDGET=", job.budget, " "); }
    format(s, job.exe, " >", shell_quoted(out), " 2>", shell_quoted(job.err));
    return s;
}

thread_fn case_worker(job: CaseJob) {
    var r = CaseResult { case: job.case, built: shell(job.build) == 0, status: 0, seconds: 0.0, out: [], err: [] };
    if r.built {
        let start = clock();
        r.status = shell(case_command(job, "", job.out));
        r.seconds = clock() - start;
        read_file(job.out, r.out);
    }
    var err: u8[>..] = [];
    read_file(job.err, err);
    r.err.append(err[..min(err.len, 1024)]);
    qput(r);
}

// The end of a task: its case, its measured peak in bytes (or 0), and the
// abort message of a contained worker that aborted.
struct TaskEnd {
    case: i64,
    ended: bool,
    peak: f64,
    why: u8[..1024],
}

// Runs case_worker on `job` as a contained worker, and gives its end with
// the peak_mb of its records.
thread_fn task_worker(job: CaseJob) {
    var why: u8[>..] = [];
    var e = TaskEnd { case: job.case, ended: thread_wait(thread_spawn_contained(case_worker, job), why), peak: 0.0, why: [] };
    var out: u8[>..] = [];
    if read_file(job.out, out) { e.peak = 1.0e6 * max(0.0, toon_number(out, "peak_mb:")); }
    set1024(e.why, why);
    qput(e);
}

// Runs `jobs` in case order, at most `limit` at a time, each when the
// ledger admits its estimate from `dofs[c]` (spec §20.5). A case that
// starts when no case runs, and that no second case can join, has all
// hardware threads and runs alone. Each other case has 1 thread. It writes the
// ledger on stderr, learns from the peak of each case, and gives the end
// of each case, in case order. It appends the outputs to `outs`.
fn run_cases(jobs: CaseJob[:], dofs: i64[:], limit: i64, l: Ledger&) -> CaseEnd[>..] {
    var ends: CaseEnd[>..] = [];
    var need: f64[>..] = [];
    for i in jobs.len {
        ends.push(CaseEnd { ended: true, built: false, status: 0, seconds: 0.0, err: [], out_lo: 0, out_hi: 0 });
        need.push(0.0);
    }
    var ids: i64[>..] = [];
    var next = 0;
    var running = 0;
    var ended = 0;
    var wide = false;   // a case on more than 1 thread runs
    loop {
        while next < jobs.len && running < limit && !wide {
            let c = jobs[next].case;
            let n = case_need(l, dofs[c]);
            guard running == 0 || admits(l, n, running) else { break; };
            var job = jobs[next];
            job.threads = 1;
            var two = l;
            two.used += n;
            if running == 0 && (next + 1 == jobs.len || limit < 2 || !admits(two, case_need(l, dofs[jobs[next + 1].case]), 1)) {
                job.threads = max(1, hardware_threads());
                wide = true;
            }
            need[next] = n;
            l.used += n;
            running++;
            l.most = max(l.most, running);
            write_stderr(str("skein run: case ", c + 1, " starts (threads: ", job.threads, "): about ", mb(n), " MB, ", mb(l.used), " of ", mb(l.budget),
                             " MB of the budget in use\n"));
            ids.push(thread_spawn(task_worker, job));
            next++;
        }
        guard running > 0 else { break; };
        let e = qget<TaskEnd>();
        let i = e.case - jobs[0].case;
        running--;
        l.used -= need[i];
        learn(l, e.case, dofs[e.case], e.peak);
        wide = false;
        write_stderr(str("skein run: case ", e.case + 1, " ended: peak ", mb(e.peak), " MB\n"));
        if e.ended { ended++; } else {
            ends[i].ended = false;
            set1024(ends[i].err, str("the worker stopped: ", e.why));
        }
    }
    for id in ids { thread_wait(id); }
    for i in 0..ended {
        let r = qget<CaseResult>();
        var e = CaseEnd { ended: true, built: r.built, status: r.status, seconds: r.seconds, err: r.err, out_lo: outs.len, out_hi: 0 };
        outs.append(r.out);
        e.out_hi = outs.len;
        ends[r.case - jobs[0].case] = e;
    }
    return ends;
}

fn mb(bytes: f64) -> i64 { return round(bytes / 1.0e6) as! i64; }

// Builds the program of `job` and runs it with SKEIN_MEMORY_PROBE, so that
// it prints the estimate of its first solve on 1 thread and stops. It
// gives the estimate in bytes, or 0. The job then needs no build.
fn probe_case(job: CaseJob&) -> f64 {
    guard shell(job.build) == 0 else { return 0.0; };
    set4096(job.build, "true");
    let probe = str(job.out, ".probe");
    guard shell(case_command(job, "SKEIN_MEMORY_PROBE=1 ", probe)) == 0 else { return 0.0; };
    var out: u8[>..] = [];
    read_file(probe, out);
    return max(0.0, toon_number(out, "estimate_bytes:"));
}

fn set4096(into: u8[..4096]&, s: u8[:]) {
    into = [];
    into.append(s[..min(s.len, 4096)]);
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
    var stage = -1;
    for p in varied { stage = max(stage, check::param_stage(prog, p.name)); }
    // The program of each case.
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    var case_jobs: CaseJob[>..] = [];
    var model: u8[>..] = [];
    var share = no_share();
    var dofs: i64[>..] = [];
    var most_steps = 0;
    for c in 0..case_count {
        let cp = syntax::parse_model(path, case_source(text, c));
        let lane = if cp >= 0 { read_model(cp) } else { -1 };
        guard lane >= 0 else {
            fmt::write_errors();
            return 1;
        };
        model.clear();
        model.append(model_name);
        if c == 0 && lane != LANE_SHELL && lane != LANE_FLOW {
            share.stages = check::STAGE_MESH - max(stage, check::STAGE_POST);
            format(share.dir, ".skein/studies/", model, "/", study, "/shared");
            shell(str("rm -rf ", shell_quoted(share.dir)));
            if share.stages > 0 && !make_dirs(share.dir) {
                write_stderr("skein run: cannot make the directory .skein\n");
                return 1;
            }
        }
        share.reads = c > 0;
        dofs.push(dof_estimate(lane));
        most_steps = max(most_steps, step_count(lane));
        let build_dir = str(".skein/build/", model, "/", study, "/", c + 1);
        let case_dir = str(".skein/studies/", model, "/", study, "/", c + 1);
        if !make_dirs(build_dir) || !make_dirs(case_dir) {
            write_stderr("skein run: cannot make the directory .skein\n");
            return 1;
        }
        var program: u8[>..] = [];
        if !write_program(lane, case_dir, share, program) { return 1; }
        let source = str(build_dir, "/main.goose");
        let exe = str(build_dir, "/", model);
        if !write_file(source, program) {
            write_stderr(str("skein run: cannot write '", source, "'\n"));
            return 1;
        }
        let out = str(case_dir, "/run.toon");
        let err = str(case_dir, "/run.err");
        var job = CaseJob { case: c, threads: 1, budget: 0, build: [], exe: [], out: [], err: [] };
        format(job.build, shell_quoted(goose), " build -q ", shell_quoted(source), " -o ", shell_quoted(exe), " >", shell_quoted(err), " 2>&1");
        job.exe.append(shell_quoted(str("./", exe)));
        set1024(job.out, out);
        set1024(job.err, err);
        case_jobs.push(job);
    }
    // At most min(H, number of cases) cases at a time, or --jobs, and each
    // when the ledger admits it (spec §20.5). A case has 1 thread when a
    // second case can run beside it. Reason: on sample 01, 16 threads took
    // 25% more memory than 1 thread, and the same time.
    let h = max(1, hardware_threads());
    let limit = min(if jobs > 0 { jobs } else { h }, case_count);
    let free = solvers::free_bytes();
    var l = new_ledger(solvers::BUDGET_SHARE * free);
    l.partial = most_steps > 1;
    for &job in case_jobs { job.budget = l.budget as! i64; }
    outs.clear();
    var first = if share.stages > 0 { 1 } else { 0 };
    if l.budget > 0.0 && limit > 1 && first == 0 {
        // The cases with the fewest and the most dofs, which read no stage
        // of case 1, and the larger bytes for each dof of the two. Reason:
        // the bytes for each dof change with the dofs, in the two
        // directions: the RSS before the solve and the GMRES restart of a
        // wave solve do not grow with the dofs, and the skyline profile of a
        // shell grows faster. On a PML duct study, the per-dof bytes of the
        // largest case gave the smaller cases 30% to 39% less than their peaks.
        var small = 0;
        var big = 0;
        for d, c in dofs {
            if d < dofs[small] { small = c; }
            if d > dofs[big] { big = c; }
        }
        for c in [small, big] {
            if c == big && big == small { break; }
            let estimate = probe_case(case_jobs[c]);
            if estimate > 0.0 && estimate / dofs[c] > l.per_dof {
                l.per_dof = estimate / dofs[c];
                l.source = c;
            }
        }
    }
    if l.budget > 0.0 {
        write_stderr(str("skein run: memory budget ", mb(l.budget), " MB (", solvers::BUDGET_SHARE, " of ", mb(free),
                         " MB free: MemAvailable and the cgroup limit)"));
        if l.source >= 0 { write_stderr(str(", case ", l.source + 1, " estimates ", mb(l.per_dof * dofs[l.source]), " MB for ", dofs[l.source], " dofs")); }
        write_stderr("\n");
    }
    var ends = run_cases(case_jobs[..first], dofs, 1, l);
    ends.append(run_cases(case_jobs[first..], dofs, limit, l));
    // The tokens of the cases replaced those of the model: read it again
    // for the columns and the place of the study.
    let p = syntax::parse_model(path, text);
    let hash = syntax::source_hash();
    var printed: u8[>..] = [];
    var table: u8[>..] = [];
    let status = study_table(study, tok(find_study(p, study)), l.most, p, ends, json, printed, table);
    let file = str(".skein/studies/", model, "/", study);
    if !write_file(str(file, ".toon"), table) || !write_file(str(file, ".hash"), hash) {
        write_stderr(str("skein run: cannot write '", file, ".toon'\n"));
        return 1;
    }
    write_stdout(printed);
    return status;
}
