import std;
import os;
namespace solvers;

struct MemoryPath {
    name: u8[..32],
    entry_bytes: f64,
    dof_bytes: f64,
    thread_bytes: f64,
}

let MEMORY_PATHS = [
    MemoryPath { name: "box cg jacobi", entry_bytes: 1.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "box cg amg", entry_bytes: 64.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "box cg amg bending", entry_bytes: 95.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "box gmres ilu0", entry_bytes: 18.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex cg amg", entry_bytes: 32.0, dof_bytes: 1600.0, thread_bytes: 16.0 },
    MemoryPath { name: "simplex elastic cg amg", entry_bytes: 86.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "simplex modal lobpcg amg", entry_bytes: 73.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex buckling lobpcg amg", entry_bytes: 40.0, dof_bytes: 0.0, thread_bytes: 1.0 },
    MemoryPath { name: "simplex central_difference", entry_bytes: 9.0, dof_bytes: 0.0, thread_bytes: 12.0 },
    MemoryPath { name: "simplex newton cg amg", entry_bytes: 32.0, dof_bytes: 1600.0, thread_bytes: 16.0 },
    MemoryPath { name: "lattice ldlt", entry_bytes: 9.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "shell skyline ldlt", entry_bytes: 9.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex newton skyline ldlt", entry_bytes: 16.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex harmonic skyline ldlt", entry_bytes: 32.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex wave gmres ilu0", entry_bytes: 18.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "simplex wave lu", entry_bytes: 24.0, dof_bytes: 150.0, thread_bytes: 0.0 },
    MemoryPath { name: "flow newton gmres ilu0", entry_bytes: 25.0, dof_bytes: 0.0, thread_bytes: 16.0 },
    MemoryPath { name: "grid gmres ilu0", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "grid ssprk3", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "cons ssprk3", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "vof cg amg", entry_bytes: 64.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "meshless skyline lu", entry_bytes: 16.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "simplex mixed cg", entry_bytes: 16.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "radiation dense", entry_bytes: 32.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "q4 band cholesky", entry_bytes: 8.0, dof_bytes: 0.0, thread_bytes: 0.0 },
    MemoryPath { name: "ode backward differences", entry_bytes: 0.0, dof_bytes: 0.0, thread_bytes: 0.0 },
];

let THREAD_BYTES = 1.5e6;

let BUDGET_SHARE = 0.9;

var budget = -1.0;
var budget_from_env = false;
var largest_estimate = 0.0;

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

fn cgroup_bytes(path: u8[:]) -> f64 {
    var s: u8[>..] = [];
    guard read_file(path, s) else { return -1.0; };
    let v, good = parse_int(trim(s));
    return if good { v * 1.0 } else { -1.0 };
}

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

fn status_bytes(key: u8[:]) -> f64 {
    var s: u8[>..] = [];
    guard read_file("/proc/self/status", s) else { return 0.0; };
    return 1024.0 * max(0.0, proc_kb(s, key));
}

fn resident_bytes() -> f64 { return status_bytes("VmRSS:"); }

fn peak_bytes() -> f64 { return status_bytes("VmHWM:"); }

fn memory_path(name: u8[:]) -> MemoryPath {
    for p in MEMORY_PATHS {
        if p.name == name { return p; }
    }
    abort(str("memory_path: the solution path '", name, "' has no memory model in MEMORY_PATHS"));
}

fn path_bytes(p: MemoryPath, entries: i64, dofs: i64, vectors: i64, mesh: i64, threads: i64) -> f64 {
    guard entries >= 0 && dofs >= 0 && vectors >= 0 && mesh >= 0 && threads >= 1 else { abort("path_bytes: entries, dofs, vectors and mesh >= 0, and threads >= 1"); }
    let bytes = p.entry_bytes * entries + (p.dof_bytes + 8.0 * vectors) * dofs;
    guard threads > 1 else { return bytes; };
    return bytes + p.thread_bytes * entries + (THREAD_BYTES + 2.0 * mesh) * threads;
}

fn memory_estimate(name: u8[:], entries: i64, dofs: i64, vectors: i64, mesh: i64) -> f64 {
    return resident_bytes() + path_bytes(memory_path(name), entries, dofs, vectors, mesh, case_threads());
}

var memory_checks = 0;

fn memory_fits(name: u8[:], entries: i64, dofs: i64, vectors: i64, mesh: i64) -> bool {
    let have = memory_budget();
    return have <= 0.0 || resident_bytes() + path_bytes(memory_path(name), entries, dofs, vectors, mesh, 1) <= have;
}

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

fn memory_records() -> u8[>..] {
    return str("memory_checks: ", memory_checks, "\nestimate_mb: ", mb(largest_estimate), "\npeak_mb: ", mb(peak_bytes()), "\n");
}
