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
import grid.gen;       // skein-p4-grid: the fv and dg models
import optim.gen;      // skein-p5-optim: the quad models of linear_elastic, and topology
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
let LANE_GRID = 4;
let LANE_OPTIM = 5;

// Reads program `prog` with the generator that claims it. It gives that
// generator, or -1 on an error, and `syntax::errors` holds the errors.
// skein-elements: a simplex model goes to src/elements/gen.goose.
// skein-shells: a shell model goes to src/shells/gen.goose.
// skein-flow: a flow model goes to src/flow/gen.goose.
// skein-p4-grid: an fv or dg model goes to src/grid/gen.goose.
// skein-p5-optim: a quad model of linear_elastic goes to src/optim/gen.goose.
fn read_model(prog: i64) -> i64 {
    if shells_gen::claims(prog) { return if shells_gen::read(prog) { LANE_SHELL } else { -1 }; }
    // Edge and face elements go to the simplex generator, before flow claims a q and a p.
    if elements_gen::family_claims(prog) { return if elements_gen::read(prog) { LANE_SIMPLEX } else { -1 }; }
    if flow_gen::claims(prog) { return if flow_gen::read(prog) { LANE_FLOW } else { -1 }; }
    if grid_gen::claims(prog) { return if grid_gen::read(prog) { LANE_GRID } else { -1 }; }
    if optim_gen::claims(prog) { return if optim_gen::read(prog) { LANE_OPTIM } else { -1 }; }
    if elements_gen::claims(prog) { return if elements_gen::read(prog) { LANE_SIMPLEX } else { -1 }; }
    return if read(prog) { LANE_BOX } else { -1 };
}

// The dofs of the mesh of the model that `read_model` read with generator
// `lane`.
fn dof_estimate(lane: i64) -> i64 {
    if lane == LANE_SHELL { return shells_gen::dof_estimate(); }
    if lane == LANE_FLOW { return flow_gen::dof_estimate(); }
    if lane == LANE_GRID { return grid_gen::dof_estimate(); }
    if lane == LANE_OPTIM { return optim_gen::dof_estimate(); }
    if lane == LANE_SIMPLEX { return elements_gen::dof_estimate(); }
    return codegen::dof_estimate();
}

// The steps of the model that `read_model` read with generator `lane`.
// Each step solves.
fn step_count(lane: i64) -> i64 {
    if lane == LANE_SHELL { return shells_gen::steps.len; }
    if lane == LANE_FLOW { return flow_gen::steps.len; }
    if lane == LANE_GRID { return grid_gen::steps.len; }
    if lane == LANE_OPTIM { return optim_gen::steps.len; }
    if lane == LANE_SIMPLEX { return elements_gen::steps.len; }
    return steps.len;
}

// Writes the program of the model that `read_model` read with generator
// `lane`. The program puts its results in `results_dir`. In a study, it
// writes or reads the stages of `share`. The shell, flow and grid
// generators share none. For a flow model, it also meshes the domain and writes the
// mesh in `results_dir` (src/flow/gen.goose). False, with an error on
// stderr, if it cannot.
fn write_program(lane: i64, results_dir: u8[:], share: Share, program: u8[>..]&) -> bool {
    // The skein executable is build/skein in the repository, so the runtime
    // is in src/runtime beside it.
    let runtime_path = str(resource_dir(), if lane == LANE_SHELL { "../src/shells" } else if lane == LANE_SIMPLEX { "../src/elements" }
                                           else if lane == LANE_FLOW { "../src/flow" } else if lane == LANE_GRID { "../src/grid" }
                                           else if lane == LANE_OPTIM { "../src/optim" } else { "../src/runtime" });
    if lane == LANE_FLOW {
        guard flow_gen::write_mesh(results_dir) else { return false; }
        flow_gen::generate(runtime_path, results_dir, program);
    } else if lane == LANE_GRID {
        grid_gen::generate(runtime_path, results_dir, program);
    } else if lane == LANE_OPTIM {
        optim_gen::generate(runtime_path, results_dir, program);
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
    if name(syntax::nodes[s].first) == "minimize" { return run_minimize(path, text, prog, s, study, json); }
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
        optim_gen::topology_study = if cp >= 0 && name(syntax::nodes[find_study(cp, study)].first) == "topology" { find_study(cp, study) } else { -1 };
        let lane = if cp >= 0 { read_model(cp) } else { -1 };
        guard optim_gen::topology_study < 0 || lane == LANE_OPTIM || lane < 0 else {
            write_stderr("skein run: a topology study runs on a quad mesh of a rect with linear_elastic (src/optim/gen.goose)\n");
            return 1;
        };
        guard lane >= 0 else {
            fmt::write_errors();
            return 1;
        };
        model.clear();
        model.append(model_name);
        if c == 0 && lane != LANE_SHELL && lane != LANE_FLOW && lane != LANE_GRID && lane != LANE_OPTIM {
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
    let status = study_table(study, tok(find_study(p, study)), l.most, p, ends, -1, "", json, printed, table);
    let file = str(".skein/studies/", model, "/", study);
    if !write_file(str(file, ".toon"), table) || !write_file(str(file, ".hash"), hash) {
        write_stderr(str("skein run: cannot write '", file, ".toon'\n"));
        return 1;
    }
    write_stdout(printed);
    return status;
}

// --- minimize with gradients (spec §20.2) -----------------------------------------

// The texts of the param values z, each in 0..1 of its range.
fn range_texts(z: f64[:]) -> u8[..256][>..] {
    var out: u8[..256][>..] = [];
    for x, k in z {
        var t: u8[..256] = [];
        set256(t, str(range_lo[k] + x * (range_hi[k] - range_lo[k]), if range_unit[k].len > 0 { " " } else { "" }, range_unit[k]));
        out.push(t);
    }
    return out;
}

// Reads the model `text` with the param values z and `tail` on the optim
// lane: the design values and the values that a param must not change.
// False on an error, with the message written.
fn read_design(path: u8[:], text: u8[:], z: f64[:], tail: u8[:]) -> f64[>..], u8[>..], bool {
    let cp = syntax::parse_model(path, values_source(text, range_texts(z), tail));
    let lane = if cp >= 0 { read_model(cp) } else { -1 };
    guard lane == LANE_OPTIM else {
        if lane < 0 { fmt::write_errors(); } else {
            write_stderr("skein run: minimize with method: gradient runs on a quad mesh of a rect with linear_elastic (src/optim/gen.goose)\n");
        }
        return [], [], false;
    };
    return optim_gen::design_values(), optim_gen::frozen(), true;
}

// dv/dp of each design value v by each param p at z, by central
// differences of the values that the generator reads, in the order
// chain[i * params + k]. A difference of values that the model text gives
// is exact for a value linear in p. False if a param changes the mesh, a
// region, a fixed value or the point of the objective.
fn design_chain(path: u8[:], text: u8[:], z: f64[:]) -> f64[>..], bool {
    let h = 1.0e-4;
    let np = z.len;
    var rows: f64[>..] = [];   // rows[k * nv + i] = dv_i/dp_k
    var frozen: u8[>..] = [];
    for k in np {
        var zp: f64[>..] = [];
        zp.append(z);
        zp[k] += h;
        let vp, fp, okp = read_design(path, text, zp, "");
        zp[k] -= 2.0 * h;
        let vm, fm, okm = read_design(path, text, zp, "");
        guard okp && okm else { return [], false; };
        if k == 0 { frozen.append(fp); }
        guard fp == frozen && fm == frozen else {
            write_stderr(str("skein run: the param '", varied[k].name, "' changes the mesh, a region, a fixed value or the point of the objective: method gradient varies E, nu, t, a load or the objective (src/optim/gen.goose)\n"));
            return [], false;
        };
        for i in vp.len { rows.push((vp[i] - vm[i]) / (2.0 * h * (range_hi[k] - range_lo[k]))); }
    }
    let nv = rows.len / np;
    var chain: f64[>..] = [];
    for i in nv { for k in np { chain.push(rows[k * nv + i]); } }
    return chain, true;
}

// The objective `objective` and the gradient dJ/dz in the records of a
// case. ran is false if a record is missing: the program of a case with a
// fault prints no gradient, and a failed expect is not a fault.
fn case_gradient(e: CaseEnd, objective: u8[:]) -> f64, f64[>..], bool {
    var j = 0.0;
    var found = false;
    var g: f64[>..] = [];
    for k in varied.len { g.push(0.0); }
    var got = 0;
    var section: u8[..16] = [];
    for line in split(outs[e.out_lo..e.out_hi], "\n") {
        if line.len > 0 && line[0] != ' ' {
            section = [];
            if starts_with(line, "report[") { section.append("report"); }
            if starts_with(line, "gradient[") { section.append("gradient"); }
            continue;
        }
        let f = toon_fields(line);
        if section == "report" && f.len == 3 && starts_with(f[1], str(objective, " =")) {
            let parts = split(f[2], " ");
            let v, good = parse_flt(parts[0]);
            if good { j = v; found = true; }
        }
        if section == "gradient" && f.len == 2 {
            for w, k in varied {
                let v, good = parse_flt(f[1]);
                if w.name == f[0] && good { g[k] = v * (range_hi[k] - range_lo[k]); got++; }
            }
        }
    }
    return j, g, found && got == varied.len && e.ended && e.built;
}

// Builds and runs case c of the minimize loop at z, with `tail` at the end
// of the model text. It gives the end of the case, J and dJ/dz. ran is
// false if the case gave no J or no gradient, and ok is false on an error
// of the model, with the message written.
fn gradient_case(path: u8[:], text: u8[:], study: u8[:], c: i64, z: f64[:], tail: u8[:], l: Ledger&) -> CaseEnd, f64, f64[>..], bool, bool {
    var none = CaseEnd { ended: false, built: false, status: 0, seconds: 0.0, err: [], out_lo: 0, out_hi: 0 };
    let chain, good = design_chain(path, text, z);
    guard good else { return none, 0.0, [], false, false; };
    let _, _, read = read_design(path, text, z, tail);
    guard read else { return none, 0.0, [], false, false; };
    optim_gen::chain.clear();
    optim_gen::chain.append(chain);
    let build_dir = str(".skein/build/", model_name, "/", study, "/", c + 1);
    let case_dir = str(".skein/studies/", model_name, "/", study, "/", c + 1);
    var program: u8[>..] = [];
    guard make_dirs(build_dir) && make_dirs(case_dir) && write_program(LANE_OPTIM, case_dir, no_share(), program) else {
        write_stderr("skein run: cannot make the directory .skein\n");
        return none, 0.0, [], false, false;
    };
    let source = str(build_dir, "/main.goose");
    guard write_file(source, program) else {
        write_stderr(str("skein run: cannot write '", source, "'\n"));
        return none, 0.0, [], false, false;
    };
    var job = CaseJob { case: c, threads: 1, budget: l.budget as! i64, build: [], exe: [], out: [], err: [] };
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    format(job.build, shell_quoted(goose), " build -q ", shell_quoted(source), " -o ", shell_quoted(str(build_dir, "/", model_name)), " >",
           shell_quoted(str(case_dir, "/run.err")), " 2>&1");
    job.exe.append(shell_quoted(str("./", build_dir, "/", model_name)));
    set1024(job.out, str(case_dir, "/run.toon"));
    set1024(job.err, str(case_dir, "/run.err"));
    var dofs: i64[>..] = [];
    for i in c + 1 { dofs.push(optim_gen::dof_estimate()); }
    let ends = run_cases([job], dofs, 1, l);
    let j, g, ran = case_gradient(ends[0], optim_gen::gradient_objective);
    return ends[0], j, g, ran, true;
}

fn largest(v: f64[:]) -> f64 {
    var m = 0.0;
    for x in v { m = max(m, abs(x)); }
    return m;
}

fn dot_of(a: f64[:], b: f64[:]) -> f64 {
    var s = 0.0;
    for x, i in a { s += x * b[i]; }
    return s;
}

// The name of the objective of minimize `pr`: the first argument, or the
// argument `objective`.
fn objective_name(pr: i64) -> const u8[:] {
    var a = syntax::nodes[pr].first;
    while a >= 0 {
        let aux = syntax::nodes[a].aux;
        if aux < 0 || syntax::text(aux) == "objective" { return name(syntax::nodes[a].first); }
        a = syntax::nodes[a].next;
    }
    return "";
}

// `skein run --study` of a minimize study with method: gradient: a
// projected gradient on z, the param values in 0..1 of their ranges. The
// first step moves the param of the largest derivative by `step` of its
// range, and each next step is that of Barzilai and Borwein, s.s / s.y. A
// step that does not lower J by 1e-4 of its first-order estimate (Armijo)
// is halved, after the cut of the bounds. The loop stops when a step moves no param by xtol of its
// range, or after max_cases cases. Each case solves the step and its
// adjoint once. The expects of the study then run in the best case.
fn run_minimize(path: u8[:], text: u8[:], prog: i64, s: i64, study: u8[:], json: bool) -> i64 {
    set64(optim_gen::gradient_objective, objective_name(syntax::nodes[s].first));
    optim_gen::gradient_params.clear();
    for w in varied { optim_gen::gradient_params.push(w.name); }
    let tail = str(case_expects[expect_at[0]..expect_at[1]]);
    var l = new_ledger(solvers::BUDGET_SHARE * solvers::free_bytes());
    outs.clear();
    case_values.clear();
    var ends: CaseEnd[>..] = [];
    var z: f64[>..] = [];
    z.append(range_start);
    let e0, j0, g0, ran0, ok0 = gradient_case(path, text, study, 0, z, "", l);
    guard ok0 else { return 1; };
    case_values.append(range_texts(z));
    ends.push(e0);
    var stop: u8[>..] = str(if ran0 { "max_cases" } else { "case 1 failed" });
    var best = 0;
    var j = j0;
    var g = g0;
    var alpha = first_step / max(largest(g), 1.0e-300);
    while ran0 && ends.len < max_cases {
        var trial: f64[>..] = [];
        var d: f64[>..] = [];
        for x, k in z {
            trial.push(min(1.0, max(0.0, x - alpha * g[k])));
            d.push(trial[k] - x);
        }
        if largest(d) < xtol { stop = str("xtol"); break; }
        let e, jt, gt, ran, ok = gradient_case(path, text, study, ends.len, trial, "", l);
        guard ok else { return 1; };
        case_values.append(range_texts(trial));
        ends.push(e);
        if !ran { stop = str("case ", ends.len, " failed"); break; }
        if jt <= j + 1.0e-4 * dot_of(g, d) {
            var y: f64[>..] = [];
            for x, k in gt { y.push(x - g[k]); }
            let sy = dot_of(d, y);
            alpha = if sy > 0.0 { dot_of(d, d) / sy } else { 2.0 * alpha };
            z = copy(trial);
            j = jt;
            g = copy(gt);
            best = ends.len - 1;
            if largest(d) < xtol { stop = str("xtol"); break; }
        } else {
            // Half the step that the bounds left. Reason: a long step
            // that ends on a bound gives the same trial for each halving.
            alpha = 0.5 * min(alpha, largest(d) / max(largest(g), 1.0e-300));
        }
    }
    if tail.len > 0 {
        let e, _, _, _, ok = gradient_case(path, text, study, best, z, tail, l);
        guard ok else { return 1; };
        ends[best] = e;
    }
    optim_gen::gradient_objective = [];
    case_count = ends.len;
    let p = syntax::parse_model(path, text);
    let hash = syntax::source_hash();
    var printed: u8[>..] = [];
    var table: u8[>..] = [];
    let status = study_table(study, tok(find_study(p, study)), 1, p, ends, best, stop, json, printed, table);
    let file = str(".skein/studies/", model_name, "/", study);
    if !write_file(str(file, ".toon"), table) || !write_file(str(file, ".hash"), hash) {
        write_stderr(str("skein run: cannot write '", file, ".toon'\n"));
        return 1;
    }
    write_stdout(printed);
    return status;
}
