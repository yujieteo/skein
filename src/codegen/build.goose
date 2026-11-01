// `skein build <file.skn> --target wasm` (docs/command_line.md §7): the
// page of a model that meshes and solves in the browser (implementation.md
// §8). It writes, in .skein/build/<model>/web/:
//
//     main.goose      the program of the page (web.goose)
//     main.c          its C, from goose
//     profile.c       main.c with the profile of the data stacks
//     native.toon     the run of profile.c at the default values
//     corner.toon     its run at the corner of the ranges with the most dofs
//     wasm.c          main.c with the regions of WebAssembly
//     <model>.html    the page, with the module in it
//
// A Goose data stack is a region of reserved address space with a guard
// gap behind it. WebAssembly memory has no reservation and no guard
// pages, so each stack gets a fixed region of memory64. Its size is the
// larger peak of the stack in two native runs, at the default values and
// at the corner of the param ranges with the most dofs, times 1.5 (1 MB
// at least). A canary gap behind each region stops the run with an error
// if a stack overruns it. The memory of the module grows to the regions
// and 1.5 times the larger heap peak of the two runs. A model that needs
// more than the browser gives (MAX_MEMORY_MB) gets the page without the
// module, with the native results and its size.
import std;
import os;
import math;
namespace codegen;
import syntax from "../syntax";
import elements.gen;
import .codegen;
import .study;
import .run;
import .web;
import .page;

let BUILD_USAGE = "usage: skein build <file.skn> --target wasm\n";

// The controls of program `prog`: each `param` with a range of two
// literals in the unit of its default value. False, with an error, for
// another range.
fn controls(prog: i64, out: Control[>..]&) -> bool {
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        let range = syntax::child_of_kind(c, syntax::N_RANGE);
        if kind(c) == syntax::N_PARAM && range >= 0 {
            let v, v_unit, v_ok = literal(syntax::child(c, 1));
            let lo, lo_unit, lo_ok = literal(syntax::nodes[range].first);
            let hi, hi_unit, hi_ok = literal(syntax::nodes[syntax::nodes[range].first].next);
            guard v_ok && lo_ok && hi_ok && v_unit == lo_unit && v_unit == hi_unit else {
                unsupported(tok(c), "a param of the page needs a literal default and range in one unit",
                            "write the range in the unit of the default, for example L: Length = 200 mm in 100 mm..400 mm");
                return false;
            };
            var k = Control { name: [], unit: [], lo: lo, hi: hi, value: v };
            set64(k.name, name(c));
            set64(k.unit, v_unit);
            out.push(k);
        }
        c = syntax::nodes[c].next;
    }
    return true;
}

fn set_word(k: Control, v: f64) -> u8[>..] { return str(k.name, "=", v, if k.unit.len > 0 { " " } else { "" }, k.unit); }

// The literals of `base` that change at the ends of the ranges of
// `params`, in `slots`, and the corner of the ranges with the most dofs,
// as words of --set, in `corner`: the end of each range with more dofs
// than the other end. False, with an error on stderr, if a value changes
// the program in another place.
fn find_slots(path: u8[:], text: u8[:], base: u8[:], params: Control[:], slots: i64[>..]&, corner: u8[..256][>..]&) -> bool {
    let lb = literals(base);
    var changed: bool[>..] = [];
    for l in lb { changed.push(false); }
    var program: u8[>..] = [];
    for k in params {
        var most = -1;
        var end: u8[..256] = [];
        for v in [k.lo, k.hi] {
            let word = set_word(k, v);
            var words: const u8[:][>..] = [];
            words.push(word);
            guard web_program(path, text, words, program) else { return false; };
            let dofs = web_dof_estimate();
            if dofs >= most {
                most = dofs;
                set256(end, word);
            }
            let lp = literals(program);
            guard same_text(base, lb, program, lp) else {
                write_stderr(str("skein build: ", word, " changes the structure of the program: --target wasm cannot make '",
                                 k.name, "' a control\n"));
                return false;
            };
            for l, i in lp {
                guard program[l.lo..l.hi] != base[lb[i].lo..lb[i].hi] else { continue; };
                guard l.float == lb[i].float else {
                    write_stderr(str("skein build: ", word, " changes the type of a value of the program\n"));
                    return false;
                };
                changed[i] = true;
            }
        }
        corner.push(end);
    }
    slots.clear();
    for c, i in changed {
        if c { slots.push(i); }
    }
    return true;
}

// Appends `s` as a Goose string literal.
fn goose_string(out: u8[>..]&, s: u8[:]) {
    out.push('"');
    for c in s {
        if c == '"' || c == '\\' { out.push('\\'); out.push(c); }
        else if c == '\n' { out.append("\\n"); }
        else if c == '\t' { out.append("\\t"); }
        else if c == '\r' { out.append("\\r"); }
        else { out.push(c); }
    }
    out.push('"');
}

// The solves that the program of the page calls through src/gpu/: the
// GPU path when the page has WebGPU, else the solve of the package.
let GPU_CALLS = [
    ["elements::solve()", "gpu::solve()"],
    ["elements::central_difference(", "gpu::central_difference("],
    ["grid::solve_transient(", "gpu::grid_transient("],
    ["grid::solve_static(", "gpu::grid_static("],
    ["elements::modal(", "gpu::modal("],
    ["elements::buckling(", "gpu::buckling("],
    ["meshless::solve_static(", "gpu::meshless_static("],
];

// `line` with each call of GPU_CALLS replaced.
fn gpu_line(line: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    var at = 0;
    while at < line.len {
        var hit = false;
        for c in GPU_CALLS {
            if !hit && line.len - at >= c[0].len && line[at..at + c[0].len] == c[0] {
                out.append(c[1]);
                at += c[0].len;
                hit = true;
            }
        }
        if !hit {
            out.push(line[at]);
            at++;
        }
    }
    return out;
}

// The program of the page: the program `base` with its slots, after the
// model, the text of the default program and the slots (web.goose). Each
// step tells src/gpu/ its name, for the record of the devices.
fn page_program(path: u8[:], text: u8[:], base: u8[:], slots: i64[:]) -> u8[>..] {
    var out: u8[>..] = [];
    format(out, "// Generated by skein build --target wasm from the model ", model_name, ". Do not edit it: edit the model.\n");
    out.append(PAGE_IMPORTS);
    out.append("let web_path = ");
    goose_string(out, path);
    out.append(";\nlet web_model = ");
    goose_string(out, text);
    out.append(";\nlet web_base = ");
    goose_string(out, base);
    out.append(";\nlet web_slots = \"");
    for s, i in slots { format(out, if i > 0 { "," } else { "" }, s); }
    out.append("\";\nvar web_f: f64[>..] = [];\nvar web_i: i64[>..] = [];\n");
    let body = slot_program(base, slots);
    var at = 0;
    while at < body.len {
        var end = at;
        while end < body.len && body[end] != '\n' { end++; }
        let line = body[at..end];
        if line == "fn main() {" { out.append("fn web_run() {"); }
        else if line == "    print_records();" { out.append("    print_records();\n    web_quality();"); }
        else if line.len > 12 && line[..8] == "fn step_" && line[line.len - 4..] == "() {" {
            format(out, line, "\n    gpu::step_name(\"", line[8..line.len - 4], "\");");
        } else if !(line.len > 7 && line[..7] == "import ") { out.append(gpu_line(line)); }
        out.push('\n');
        at = end + 1;
    }
    out.append(PAGE_MAIN);
    return out;
}

// The numbers of the line `<key>a,b,...` of `err`: the data stack
// peaks (MB) of the native run for `skein-profile: `, and its peak
// resident MB for `skein-rss: `.
fn profile_numbers(err: u8[:], key: u8[:], out: i64[>..]&) -> bool {
    var at = -1;
    for i in 0..err.len - key.len + 1 {
        if at < 0 && err[i..i + key.len] == key { at = i + key.len; }
    }
    guard at >= 0 else { return false; };
    var start = at;
    var i = at;
    while i <= err.len {
        if i == err.len || err[i] == ',' || err[i] == '\n' {
            let v, good = parse_int(err[start..i]);
            guard good else { return false; };
            out.push(v);
            if i == err.len || err[i] == '\n' { return true; }
            start = i + 1;
        }
        i++;
    }
    return true;
}

// `text` with the first `old` replaced by `new`. False if `text` has no
// `old`: the Goose runtime changed, and skein build needs an update.
fn replace_once(text: u8[>..]&, old: u8[:], new: u8[:]) -> bool {
    for i in 0..text.len - old.len + 1 {
        if text[i..i + old.len] == old {
            var out: u8[>..] = [];
            out.append(text[..i]);
            out.append(new);
            out.append(text[i + old.len..]);
            text.clear();
            text.append(out);
            return true;
        }
    }
    write_stderr(str("skein build: internal fault: the C of goose has no '", old, "'. Please report it.\n"));
    return false;
}

// main.c with a call of `hook` at each end of the program, and the C `tail`.
fn hooked_c(c: u8[:], hook: u8[:], tail: u8[:], out: u8[>..]&) -> bool {
    out.clear();
    out.append(c);
    let end = "static GS_NORETURN void gs_end(int code) {\n    if (!gs_end_first()) gs_park();\n";
    guard replace_once(out, end, str("static void ", hook, "(void);\n", end, "    ", hook, "();\n")) &&
          replace_once(out, "    main_g(0);\n", str("    main_g(0);\n    ", hook, "();\n")) else { return false; };
    out.append(tail);
    return true;
}

// The C compiler of the profile run.
fn c_compiler() -> u8[>..] {
    var cc: u8[>..] = [];
    if !env("CC", cc) || cc.len == 0 { cc.append("cc"); }
    return cc;
}

// `skein build`. `args` are the words after `build`. It gives the exit
// status (command_line.md §10).
fn build(args: const u8[:][:]) -> i64 {
    guard args.len == 3 && args[1] == "--target" && args[2] == "wasm" && args[0].len > 0 && args[0][0] != '-' else {
        write_stderr(str("skein build: give one model file and --target wasm\n", BUILD_USAGE));
        return 2;
    };
    let path = args[0];
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        write_stderr(str("skein build: cannot read '", path, "'\n"));
        return 1;
    };
    var base: u8[>..] = [];
    var none: const u8[:][>..] = [];
    guard web_program(path, text, none, base) else { return 1; };
    var params: Control[>..] = [];
    guard controls(syntax::parse_model(path, text), params) else {
        fmt::write_errors();
        return 1;
    };
    var model: u8[>..] = [];
    model.append(model_name);
    var slots: i64[>..] = [];
    var corner: u8[..256][>..] = [];
    guard find_slots(path, text, base, params, slots, corner) else { return 1; };
    let dir = str(".skein/build/", model, "/web");
    guard make_dirs(str(dir, "/.skein/results/", model)) else {
        write_stderr("skein build: cannot make the directory .skein\n");
        return 1;
    };
    guard write_file(str(dir, "/main.goose"), page_program(path, text, base, slots)) else {
        write_stderr(str("skein build: cannot write '", dir, "/main.goose'\n"));
        return 1;
    };
    return build_page(dir, model, params, corner);
}

// Links each part of the compiler in `dir`, beside main.goose: the program
// imports them as the compiler does (src/skein.goose).
fn link_sources(dir: u8[:]) -> bool {
    let src = str(resource_dir(), "../src/");
    for part in ["syntax", "check", "codegen", "elements", "shells", "flow", "grid", "meshless", "gpu", "mesher", "runtime", "solvers", "fmt"] {
        guard shell(str("ln -sfn ", shell_quoted(str(src, part)), " ", shell_quoted(str(dir, "/", part)))) == 0 else { return false; };
    }
    return true;
}

// Runs the profile of `dir` with the words of --set `sets`. Its standard
// output goes to `out` in `dir`. It gives the peak of each data stack in
// `peaks` (MB), the seconds of the run, and the peak MB of the heap: the
// peak resident MB less the data stacks. False, with an error, if the run
// fails or gives no profile.
fn profile_run(dir: u8[:], sets: u8[..256][:], out: u8[:], peaks: i64[>..]&) -> bool, f64, i64 {
    var command = str("cd ", shell_quoted(dir), " && SKEIN_THREADS=1 ./profile");
    for w in sets { format(command, " ", shell_quoted(w)); }
    format(command, " > ", out, " 2> profile.err");
    write_stderr(str("skein build: the native run", if sets.len > 0 { " at" } else { " at the default values" }));
    for w in sets { write_stderr(str(" ", w)); }
    write_stderr("\n");
    let start = clock();
    let status = shell(command);
    let seconds = clock() - start;
    var err: u8[>..] = [];
    read_file(str(dir, "/profile.err"), err);
    guard status == 0 else {
        write_stderr(str("skein build: the native run failed with exit status ", status, ":\n", err));
        return false, seconds, 0;
    };
    var rss: i64[>..] = [];
    guard profile_numbers(err, "skein-profile: ", peaks) && profile_numbers(err, "skein-rss: ", rss) && rss.len == 1 else {
        write_stderr(str("skein build: the native run gave no profile:\n", err));
        return false, seconds, 0;
    };
    var heap = rss[0];
    for p in peaks { heap -= p; }
    return true, seconds, max(heap, 0);
}

// Builds main.goose of `dir` to C, runs the profile at the default values
// and at `corner`, and writes the page. goose gets the absolute path of
// main.goose, so that each link resolves to one path.
fn build_page(dir: u8[:], model: u8[:], params: Control[:], corner: u8[..256][:]) -> i64 {
    var goose: u8[>..] = [];
    if !env("GOOSE", goose) || goose.len == 0 { goose.append("goose"); }
    let q = shell_quoted(dir);
    guard link_sources(dir) && shell(str("cd ", q, " && ", shell_quoted(goose), " --standalone -q -o main.c \"$PWD/main.goose\" 1>&2")) == 0 else {
        write_stderr(str("skein build: internal fault: goose cannot build '", dir, "/main.goose'. Please report it.\n"));
        return 3;
    };
    var c: u8[>..] = [];
    var hooked: u8[>..] = [];
    guard read_file(str(dir, "/main.c"), c) && hooked_c(c, "skein_profile", str(PROFILE_C, GPU_STUB_C), hooked) &&
          write_file(str(dir, "/profile.c"), hooked) else { return 3; };
    guard shell(str(shell_quoted(c_compiler()), " -O2 -w -o ", q, "/profile ", q, "/profile.c -lm -lpthread 1>&2")) == 0 else {
        write_stderr("skein build: the C compiler cannot build the profile (set CC)\n");
        return 1;
    };
    var peaks: i64[>..] = [];
    var more: i64[>..] = [];
    let none: u8[..256][>..] = [];
    let ran, seconds, heap = profile_run(dir, none, "native.toon", peaks);
    var heap_corner = 0;
    if ran && corner.len > 0 {
        let corner_ran, _, h = profile_run(dir, corner, "corner.toon", more);
        guard corner_ran else { return 1; };
        heap_corner = h;
    }
    guard ran && (more.len == 0 || more.len == peaks.len) else { return 1; };
    var table: u8[>..] = [];
    // The memory of the module: the regions, 1.5 times the larger heap
    // peak, and 64 MB for the loader and the files, in 64 MB steps.
    var total = ceil(1.5 * max(heap, heap_corner)) as! i64 + 64;
    var largest = 1;
    for p, i in peaks {
        let mb = max(1, ceil(1.5 * max(p, if more.len > 0 { more[i] } else { 0 })) as! i64);
        format(table, if i > 0 { "," } else { "" }, mb);
        total += mb;
        largest = max(largest, mb);
    }
    total = (total + 63) / 64 * 64;
    var glue: u8[>..] = [];
    if total <= MAX_MEMORY_MB {
        guard wasm_module(dir, model, c, table, largest, total, glue) else { return 1; };
    }
    var native: u8[>..] = [];
    read_file(str(dir, "/native.toon"), native);
    var wgsl: u8[>..] = [];
    var gpu_js: u8[>..] = [];
    guard read_file(str(resource_dir(), "../src/gpu/kernels.wgsl"), wgsl) && read_file(str(resource_dir(), "../src/gpu/gpu.js"), gpu_js) else {
        write_stderr("skein build: internal fault: cannot read src/gpu/kernels.wgsl and src/gpu/gpu.js. Please report it.\n");
        return 3;
    };
    let page = str(dir, "/", model, ".html");
    guard write_file(page, page_html(model, params, native, glue, total, seconds, wgsl, gpu_js)) else {
        write_stderr(str("skein build: cannot write '", page, "'\n"));
        return 1;
    };
    write_stdout(str("page: ", page, "\nmemory_mb: ", total, "\nsolve_in_page: ", if glue.len > 0 { "yes" } else { "no" }, "\n"));
    return 0;
}

// Compiles main.c with the region `table` (MB, comma separated) to the
// WebAssembly module and its loader, one file of JavaScript, in `glue`.
// The memory of the module starts at 64 MB and grows to `memory_mb` at
// most. A growth that the browser cannot give stops the run with an error.
// The C stack is 8 MB, as natively. Reason: the default of 64 KB overflowed
// in the run of mms_fv_transient and broke the data with no error. JSPI
// suspends the module while the GPU read-back waits. Reason: Asyncify
// made emcc take 5 min in place of 80 s.
fn wasm_module(dir: u8[:], model: u8[:], c: u8[:], table: u8[:], largest: i64, memory_mb: i64, glue: u8[>..]&) -> bool {
    var wasm: u8[>..] = [];
    var patched: u8[>..] = [];
    patched.append(c);
    guard replace_once(patched, "GS_API uint8_t *gs_reserve_region(gs_stack *stack) {", "static uint8_t *gs_reserve_region_native(gs_stack *stack) {") &&
          hooked_c(patched, "skein_check_regions", str("\nstatic const int64_t skein_region_mb[] = {", table, "};\n", WASM_C, GPU_C), wasm) &&
          write_file(str(dir, "/wasm.c"), wasm) else { return false; };
    let q = shell_quoted(dir);
    let js = str(dir, "/", model, ".js");
    let emcc = str("emcc -O3 -w ", q, "/wasm.c -o ", shell_quoted(js), " -DGS_STACK_RESERVE='(", largest, "ull<<20)' -DGS_STACK_GAP=65536",
                   " -sMEMORY64=1 -sMODULARIZE=1 -sEXPORT_NAME=createSkein -sENVIRONMENT=worker,node -sSINGLE_FILE=1",
                   " -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sMAXIMUM_MEMORY=", memory_mb, "MB -sEXIT_RUNTIME=1 -sFORCE_FILESYSTEM=1",
                   " -sEXPORTED_RUNTIME_METHODS=ENV,FS -sSTACK_SIZE=8MB -sJSPI=1 -sJSPI_IMPORTS=skein_gpu_js_read 1>&2");
    write_stderr("skein build: emcc\n");
    guard shell(emcc) == 0 && read_file(js, glue) else {
        write_stderr("skein build: emcc cannot build the module: put emcc on PATH\n");
        return false;
    };
    return true;
}
