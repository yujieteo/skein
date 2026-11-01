// The memory of a case (spec §20.5). Each solution path checks, after the
// sparsity pattern of its matrix and before the values, that the case
// fits in the free memory. The estimate counts the entries of the pattern
// that the program really built, so it needs no table of cell types or
// physics: the element, the components and the coupling of the fields
// are all in the count.
//
//     peak = RSS now + entries × (bytes of the path + bytes of the threads)
//                    + dofs × (bytes of the path + vectors × 8)
//                    + threads × (THREAD_BYTES + 2 × mesh)
//
// The RSS now holds the mesh and the pattern. `vectors` is the number of
// vectors of one value for each dof that the solver of the path holds
// (the basis of GMRES, the blocks of LOBPCG), which the caller knows. The
// constants of each path are the measured peak less the RSS at the check
// and less the vectors, for each entry. `mesh` is the bytes of the globals
// that each worker of the assembly copies at its spawn (the mesh), which
// the caller knows, or 0 if the path has no such workers.
//
// The free memory (the budget) is 0.9 × min(MemAvailable, memory.max −
// memory.current of each cgroup of the process), read once. In a study,
// the runner reads it once for all cases and gives it in
// SKEIN_MEMORY_BUDGET. Linux only: on another OS the budget is 0, and the
// check passes. With SKEIN_MEMORY_PROBE, the program prints the estimate
// of its first check and stops: the study runner probes one case so.
import std;
import os;
namespace solvers;

// The memory model of one solution path: the bytes of the peak for each
// entry of the pattern and for each dof, and the bytes more for each entry
// when the case has more than one thread (the rows of solvers::rows_start,
// and the cell matrices of the assembly workers).
struct MemoryPath {
    name: u8[..32],
    entry_bytes: f64,
    dof_bytes: f64,
    thread_bytes: f64,
}

// Each solution path that a generated program can take, and its model.
// A new path needs a line here: memory_check aborts on a name that is not
// here, and scripts/generated_c.sh checks that each solve of each model
// that skein run can run made a check.
//
// Measured on 2026-10-31 with SKEIN_THREADS=1 and 16: (peak - RSS at the
// check - 8 dofs vectors) / entries. The AMG paths add 10%. Reason: the
// peak of the AMG setup for each entry changed by 8% between three mesh
// sizes of sample 01. The thread bytes of a path with solvers::rows_start
// are its copy of the matrix (12 bytes) and the assembly. The models:
//
//     box cg jacobi       closed cylinder, 65,835 dofs: 0.2, and 16.3 more on 16 threads
//     box cg amg          closed cylinder 2 mm, 143,787 dofs: 57.9, and 0.8 more. Measured
//                         again on 2026-11-01: the 6 rigid-body modes of the AMG setup
//                         (#58) took the peak from 38.1.
//     box gmres ilu0      closed cylinder 5 mm, 10,251 dofs: 17.0 (the copy of K and its ILU(0)), and 0
//     simplex cg amg      sample 01, 319,320 dofs (elastic, 79 entries for each dof): 47.0,
//                         and 15.4 more. Measured again on 2026-11-01 with steady heat (a
//                         scalar field, 27 entries for each dof), which took 78 at 137,781
//                         dofs and 80 at 1,257,525 dofs: the memory of each dof is a larger
//                         part of the peak when each dof has fewer entries. The fit of
//                         sample 01 and the heat at 1,257,525 dofs is 28.5 for each entry
//                         and 1468 for each dof.
//     simplex modal       frequency cantilever 1.5 mm, 54,675 dofs: 66.3, and 1.0 more
//     simplex buckling    buckling column 3 mm, 32,805 dofs: 36.6, and 0 more. The static
//                         step before it leaves pages that the buckling step uses again.
//     simplex central_difference  bar_wave_tri at 0.5 mm and 0.25 mm, 32,882 and
//                         129,762 dofs, on 2026-11-01: 8.1 and 7.5, and 12 more on 16
//                         threads (the copy of the matrix of rows_start). On 16 threads the
//                         estimate is 245 MB and the peak 170 MB.
//     shell skyline ldlt  thin cylinder 38.1 mm, 23,760 dofs: 8.2 (one value for each
//                         entry of the profile), and no threads
//     simplex wave        PML duct at 6 mm and 4 mm, 35,070 and 76,602 dofs: 17.6 and 9.1
//                         (lu and the ILU(0)), and -2.5 and 14.4 more on 16 threads. The
//                         highest of each, measured on 2026-11-01. The Krylov vectors of
//                         GMRES are most of the peak, and the caller counts them.
//     flow                sample 02 at a/128, 102,675 dofs: 24.2 (val, lu and pre). Not
//                         measured on 16 threads: the 16 of simplex cg amg, which also
//                         has rows_start.
//     grid                the examples/grid models at the sizes 0.004, 0.008 and 0.002 m, on
//                         2026-11-01, with no output: fv gmres -7.7 (62,500 dofs), dg order 2
//                         gmres -0.9 (140,625 dofs), fv ssprk3 0.3 (250,000 dofs). fv and dg
//                         size the matrix and the vectors of the state before the check, so
//                         the RSS at the check holds them. No threads.
//     meshless skyline lu sample 05, 14,342 dofs, 15,389,813 entries of the profile, on
//                         2026-11-01: 14.5 (U and L, 8 bytes each). The check is before the
//                         factors are sized; the check of the solve counts 0 entries. No
//                         threads.
//                         The Nédélec elements take this path too: the coax line of
//                         examples/femahoy/grid_coax_tem.skn, 104,856 dofs, on 2026-11-01:
//                         estimate 613 MB and peak 588 MB. The AMG of the gradient
//                         correction (src/elements/edge.goose) is on the nodes, and the
//                         estimate holds it.
//     simplex mixed cg    examples/elements/darcy_annulus_tri.skn, 80,804 dofs, on
//                         2026-11-01: estimate 55 MB and peak 28 MB with 16 bytes for each
//                         entry (val and the Schur product). No threads.
//     radiation dense     counted, not measured: the n × n arrays F, M, X and Q of the
//                         facets of an enclosure (src/elements/radiation.goose). No threads.
//     q4 band cholesky    counted, not measured: the band of the Cholesky factor, 8 bytes
//                         for each entry (src/optim/band.goose). No threads.
let MEMORY_PATHS = [
    MemoryPath { name: "box cg jacobi", entry_bytes: 1.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "box cg amg", entry_bytes: 64.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "box gmres ilu0", entry_bytes: 18.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex cg amg", entry_bytes: 32.0, dof_bytes: 1600.0, thread_bytes: 16.0 },
    MemoryPath { name: "simplex modal lobpcg amg", entry_bytes: 73.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex buckling lobpcg amg", entry_bytes: 40.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex central_difference", entry_bytes: 9.0, dof_bytes: 0.0, thread_bytes: 12.0 },
    MemoryPath { name: "shell skyline ldlt", entry_bytes: 9.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex wave gmres ilu0", entry_bytes: 18.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "flow newton gmres ilu0", entry_bytes: 25.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "grid gmres ilu0", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "grid ssprk3", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "meshless skyline lu", entry_bytes: 16.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex mixed cg", entry_bytes: 16.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "radiation dense", entry_bytes: 32.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "q4 band cholesky", entry_bytes: 8.0, dof_bytes: 0.0, thread_bytes: 0.0 },
];

// The bytes of each worker thread of a case: its stacks. Measured: 16
// workers of the assembly of the closed cylinder (box lane, 143,787 dofs)
// took 23 MB more than 1 thread. Each worker also holds its spawn argument
// and its globals, two copies of `mesh`. Measured on 2026-11-01: linear
// steady heat (conduction slab, 915,756 tets, 34 MB of mesh) took 57 MB
// more for each thread.
let THREAD_BYTES = 1.5e6;

// The fraction of the free memory that the cases may take. Reason: the
// desktop and the other programs grow too.
let BUDGET_SHARE = 0.9;

var budget = -1.0;            // the budget in bytes, 0 if it is not known, or -1 before the first read
var budget_from_env = false;  // SKEIN_MEMORY_BUDGET set the budget
var largest_estimate = 0.0;   // the largest estimate of a check, in bytes

// The value in kB of the line `key` of a /proc file such as
// "MemAvailable:   1234 kB", or -1.
fn proc_kb(text: u8[:], key: u8[:]) -> f64 {
    for line in split(text, '\n') {
        if !starts_with(line, key) { continue; }
        var field = trim(line[key.len..]);
        if ends_with(field, " kB") { field = field[..field.len - 3]; }
        let v, good = parse_int(field);
        return if good { v * 1.0 } else { -1.0 };
    }
    return -1.0;
}

// The bytes of a cgroup file that holds one number, or -1 for `max` or
// no file.
fn cgroup_bytes(path: u8[:]) -> f64 {
    var s: u8[>..] = [];
    guard read_file(path, s) else { return -1.0; };
    let v, good = parse_int(trim(s));
    return if good { v * 1.0 } else { -1.0 };
}

// The free bytes of this machine for a new case: MemAvailable, and at most
// memory.max − memory.current of the cgroup of the process and of each
// cgroup above it. 0 if /proc/meminfo is not there.
fn free_bytes() -> f64 {
    var s: u8[>..] = [];
    guard read_file("/proc/meminfo", s) else { return 0.0; };
    var free = 1024.0 * proc_kb(s, "MemAvailable:");
    guard free > 0.0 else { return 0.0; };
    s.clear();
    guard read_file("/proc/self/cgroup", s) && starts_with(s, "0::/") else { return free; };
    var dir = trim(s[3..]);
    loop {
        let limit = cgroup_bytes(str("/sys/fs/cgroup", dir, "/memory.max"));
        let now = cgroup_bytes(str("/sys/fs/cgroup", dir, "/memory.current"));
        if limit > 0.0 && now >= 0.0 { free = min(free, max(0.0, limit - now)); }
        let up = rfind(dir, "/");
        if up <= 0 { break; }
        dir = dir[..up];
    }
    return free;
}

// The budget of a case in bytes: SKEIN_MEMORY_BUDGET in a study, or
// BUDGET_SHARE of the free bytes now plus the RSS that this case holds
// already. 0 if it is not known.
fn memory_budget() -> f64 {
    if budget < 0.0 {
        var s: u8[>..] = [];
        let free = free_bytes();
        budget = if free > 0.0 { BUDGET_SHARE * free + resident_bytes() } else { 0.0 };
        if env("SKEIN_MEMORY_BUDGET", s) {
            let v, good = parse_int(s);
            if good { budget = v * 1.0; budget_from_env = true; }
        }
    }
    return budget;
}

// The value in bytes of the line `key` of /proc/self/status, or 0 if
// /proc is not there.
fn status_bytes(key: u8[:]) -> f64 {
    var s: u8[>..] = [];
    guard read_file("/proc/self/status", s) else { return 0.0; };
    return 1024.0 * max(0.0, proc_kb(s, key));
}

// The resident bytes of this process now.
fn resident_bytes() -> f64 { return status_bytes("VmRSS:"); }

// The peak resident bytes of this process: VmHWM, the ru_maxrss of
// getrusage on Linux. Reason: a binding of libc.so.6 stops the WebAssembly
// page of skein build from starting, since emcc has no dlopen of it.
fn peak_bytes() -> f64 { return status_bytes("VmHWM:"); }

// The memory model of path `name`. It aborts if the path has none.
fn memory_path(name: u8[:]) -> MemoryPath {
    for p in MEMORY_PATHS {
        if p.name == name { return p; }
    }
    abort(str("memory_path: the solution path '", name, "' has no memory model in MEMORY_PATHS"));
}

// The bytes above the RSS of a solve on path `p` with a pattern of
// `entries` entries, `vectors` vectors of `dofs` values, `mesh` bytes that
// each worker copies and `threads` threads.
fn path_bytes(p: MemoryPath, entries: i64, dofs: i64, vectors: i64, mesh: i64, threads: i64) -> f64 {
    guard entries >= 0 && dofs >= 0 && vectors >= 0 && mesh >= 0 && threads >= 1 else { abort("path_bytes: entries, dofs, vectors and mesh >= 0, and threads >= 1"); }
    let bytes = p.entry_bytes * entries + (p.dof_bytes + 8.0 * vectors) * dofs;
    guard threads > 1 else { return bytes; };
    return bytes + p.thread_bytes * entries + (THREAD_BYTES + 2.0 * mesh) * threads;
}

// The peak bytes of a solve on path `name`, from the RSS now (path_bytes).
fn memory_estimate(name: u8[:], entries: i64, dofs: i64, vectors: i64, mesh: i64) -> f64 {
    return resident_bytes() + path_bytes(memory_path(name), entries, dofs, vectors, mesh, case_threads());
}

var memory_checks = 0;   // the checks of the run

// Checks, after the pattern of a solve on path `name` and before its
// values, that the case fits in the budget (memory_estimate). If it does
// not on its threads but does on 1 thread, the case goes on with 1 thread.
// If it does not on 1 thread, it stops the program with the 1-thread
// estimate, the budget and `fix`, the lower-memory choice of the path or "". It also asks the kernel to kill
// this case first if the memory runs out (oom_score_adj 1000, which a
// process may set on itself).
fn memory_check(name: u8[:], entries: i64, dofs: i64, vectors: i64, mesh: i64, fix: u8[:]) {
    var need = memory_estimate(name, entries, dofs, vectors, mesh);
    memory_checks++;
    var probe: u8[>..] = [];
    if env("SKEIN_MEMORY_PROBE", probe) {
        write_stdout(str("estimate_bytes: ", need as! i64, "\n"));
        exit(0);
    }
    write_file("/proc/self/oom_score_adj", "1000");
    let have = memory_budget();
    let one = resident_bytes() + path_bytes(memory_path(name), entries, dofs, vectors, mesh, 1);
    if have > 0.0 && need > have && one <= have {
        write_stderr(str("skein run: the solve runs on 1 thread, since ", case_threads(), " threads need about ", mb(need), " MB and ", mb(have), " MB is free\n"));
        threads = 1;
        need = one;
    }
    largest_estimate = max(largest_estimate, need);
    if have <= 0.0 || need <= have { return; }
    let source = if budget_from_env { str("from SKEIN_MEMORY_BUDGET") } else { str(BUDGET_SHARE, " of MemAvailable and of the cgroup limit") };
    write_stderr(str("skein run: the solve needs about ", mb(one), " MB of memory, and ", mb(have),
                     " MB is free for it (", source, "). ",
                     "Close other programs or use a coarser mesh", fix, ".\n"));
    exit(1);
}

fn mb(bytes: f64) -> i64 { return round(bytes / 1.0e6) as! i64; }

// The memory records of the run, in TOON: the checks, the largest
// estimate and the measured peak, in MB.
fn memory_records() -> u8[>..] {
    return str("memory_checks: ", memory_checks, "\nestimate_mb: ", mb(largest_estimate), "\npeak_mb: ", mb(peak_bytes()), "\n");
}
