// The texts of `skein build --target wasm` (build.goose): the parts of the
// program of the page, the C of the profile and of the regions of
// WebAssembly, and the solve of the page. The page is the page of
// `skein report` (src/report).
import std;
namespace codegen;

// The largest memory64 of V8, less the heap of the module.
let MAX_MEMORY_MB = 15 * 1024;

// A param with a range: a control of the page.
struct Control {
    name: u8[..64],
    unit: u8[..64],
    lo: f64,
    hi: f64,
    value: f64,
}

// The imports of the program of the page. The modules come before the
// packages, as in src/skein.goose: a package and a module of the same
// directory need that sequence.
let PAGE_IMPORTS = """
    import std;
    import os;
    import codegen.codegen;
    import codegen.web;
    import elements.gen;
    import grid.gen;
    import mesher.verify;
    import fmt.fmt;
    import syntax from "./syntax";
    import check from "./check";
    import mesher from "./mesher";
    import runtime from "./runtime";
    import solvers from "./solvers";
    import elements from "./elements";
    import grid from "./grid";
    import meshless from "./meshless";
    import gpu from "./gpu";

    """;

// The end of the program of the page. The arguments are the words of
// `--set`. web_quality prints the device record of src/gpu/ and, for a
// simplex mesh, the mesh_quality record after the records.
let PAGE_MAIN = """

    fn web_quality() {
        write_stdout(gpu::records());
        if mesher::cell_count() == 0 { return; }
        var s: u8[>..] = [];
        mesher::format(s, mesher::mesh_quality());
        s.push('\n');
        write_stdout(s);
    }

    fn main() {
        var sets: const u8[:][>..] = [];
        let words = args();
        for a in words[1..] { sets.push(a); }
        var text: u8[>..] = [];
        guard codegen::web_program(web_path, web_model, sets, text) &&
              codegen::slot_values(web_base, codegen::parse_slots(web_slots), text, web_f, web_i) else { exit(2); };
        web_run();
    }

    """;

// The profile of the native run, on stderr: the resident MB of each data
// stack region at the end, and the peak resident MB of the process.
let PROFILE_C = """

    #include <sys/mman.h>
    #include <sys/resource.h>
    #include <unistd.h>
    static void skein_profile(void) {
        long page = sysconf(_SC_PAGESIZE);
        fprintf(stderr, "skein-profile: ");
        for (long i = 0; i < gs_nregions; i++) {
            gs_region *r = &gs_regions[i];
            size_t n = r->usable / (size_t)page, top = 0;
            unsigned char *v = malloc(n);
            if (v && mincore(r->base, r->usable, v) == 0) {
                for (size_t k = 0; k < n; k++) if (v[k] & 1) top = k + 1;
            }
            free(v);
            fprintf(stderr, "%s%zu", i ? "," : "", (top * (size_t)page + (1u << 20) - 1) >> 20);
        }
        struct rusage use;
        getrusage(RUSAGE_SELF, &use);
        fprintf(stderr, "\nskein-rss: %ld\n", (use.ru_maxrss + 1023) / 1024);
        fflush(stderr);
    }

    """;

// The externs of src/gpu/gpu.goose in the native profile run: no GPU, so
// each solve takes its CPU path.
let GPU_STUB_C = """

    int32_t skein_gpu_ready(void) { return 0; }
    void skein_gpu_reset(void) {}
    int32_t skein_gpu_buffer(int64_t bytes) { (void)bytes; return -1; }
    void skein_gpu_write(int32_t buf, int64_t at, sl_f32 data) { (void)buf; (void)at; (void)data; }
    void skein_gpu_write_ints(int32_t buf, int64_t at, sl_i32 data) { (void)buf; (void)at; (void)data; }
    int32_t skein_gpu_bind(sl_u8 kernel, sl_i32 pairs) { (void)kernel; (void)pairs; return -1; }
    void skein_gpu_dispatch(int32_t group, int64_t threads) { (void)group; (void)threads; }
    void skein_gpu_submit(void) {}
    int32_t skein_gpu_read(int32_t buf, int64_t at, sl_f32 out) { (void)buf; (void)at; (void)out; return 0; }

    """;

// The externs of src/gpu/gpu.goose in the module: calls of the methods of
// Module.skeinGpu (src/gpu/gpu.js), which the worker of the page sets when
// it has a WebGPU device. The read-back waits for the GPU: JSPI suspends
// the module during the wait (skein_gpu_js_read in JSPI_IMPORTS). The bytes are copied, since Goose packs its arrays
// without alignment.
let GPU_C = """

    #include <emscripten.h>
    EM_JS(int32_t, skein_gpu_js_ready, (void), { return Module.skeinGpu ? 1 : 0; });
    EM_JS(void, skein_gpu_js_reset, (void), { Module.skeinGpu.reset(); });
    EM_JS(int32_t, skein_gpu_js_buffer, (int64_t bytes), { return Module.skeinGpu.buffer(Number(bytes)); });
    EM_JS(void, skein_gpu_js_write, (int32_t buf, int64_t at, const void *p, int64_t bytes), {
        Module.skeinGpu.write(buf, Number(at), HEAPU8.slice(Number(p), Number(p) + Number(bytes)));
    });
    EM_JS(int32_t, skein_gpu_js_bind, (const void *k, int64_t kn, const void *p, int64_t n), {
        const name = new TextDecoder().decode(HEAPU8.slice(Number(k), Number(k) + Number(kn)));
        const pairs = new Int32Array(HEAPU8.slice(Number(p), Number(p) + 4 * Number(n)).buffer);
        return Module.skeinGpu.bind(name, Array.from(pairs));
    });
    EM_JS(void, skein_gpu_js_dispatch, (int32_t group, int64_t threads), { Module.skeinGpu.dispatch(group, Number(threads)); });
    EM_JS(void, skein_gpu_js_submit, (void), { Module.skeinGpu.submit(); });
    EM_ASYNC_JS(int32_t, skein_gpu_js_read, (int32_t buf, int64_t at, void *p, int64_t bytes), {
        const out = await Module.skeinGpu.read(buf, Number(at), Number(bytes));
        if (!out) return 0;
        HEAPU8.set(out, Number(p));
        return 1;
    });

    int32_t skein_gpu_ready(void) { return skein_gpu_js_ready(); }
    void skein_gpu_reset(void) { skein_gpu_js_reset(); }
    int32_t skein_gpu_buffer(int64_t bytes) { return skein_gpu_js_buffer(bytes); }
    void skein_gpu_write(int32_t buf, int64_t at, sl_f32 data) { skein_gpu_js_write(buf, 4 * at, data.data, 4 * data.len); }
    void skein_gpu_write_ints(int32_t buf, int64_t at, sl_i32 data) { skein_gpu_js_write(buf, 4 * at, data.data, 4 * data.len); }
    int32_t skein_gpu_bind(sl_u8 kernel, sl_i32 pairs) { return skein_gpu_js_bind(kernel.data, kernel.len, pairs.data, pairs.len); }
    void skein_gpu_dispatch(int32_t group, int64_t threads) { skein_gpu_js_dispatch(group, threads); }
    void skein_gpu_submit(void) { skein_gpu_js_submit(); }
    int32_t skein_gpu_read(int32_t buf, int64_t at, sl_f32 out) { return skein_gpu_js_read(buf, 4 * at, out.data, 4 * out.len); }

    """;

// The data stack regions of WebAssembly: region i has skein_region_mb[i]
// MB of memory64, from malloc (whose memory nothing touched before), and
// a canary gap behind it. Each end of the program checks the gaps.
let WASM_C = """
    #define SKEIN_CANARY 0xA5
    #define SKEIN_GAP (64u << 10)

    static void skein_check_regions(void) {
        for (long i = 0; i < gs_nregions; i++) {
            const uint8_t *g = gs_regions[i].base + gs_regions[i].usable;
            for (size_t k = 0; k < SKEIN_GAP; k++) {
                if (g[k] != SKEIN_CANARY) {
                    fprintf(stderr, "skein web: data stack %ld overran its %lld MB: the results are not valid\n",
                            i, (long long)(gs_regions[i].usable >> 20));
                    fflush(NULL);
                    _Exit(4);
                }
            }
        }
    }

    GS_API uint8_t *gs_reserve_region(gs_stack *stack) {
        long i = gs_nregions;
        long count = (long)(sizeof skein_region_mb / sizeof skein_region_mb[0]);
        if (i >= gs_regions_cap || i >= count) gs_panic("skein web: more data stack regions than the profile");
        size_t usable = (size_t)skein_region_mb[i] << 20;
        uint8_t *p = (uint8_t *)malloc(usable + SKEIN_GAP);
        if (!p) {
            fprintf(stderr, "skein web: the browser gives no more memory for data stack %ld (%lld MB)\n",
                    i, (long long)skein_region_mb[i]);
            fflush(NULL);
            _Exit(4);
        }
        memset(p + usable, SKEIN_CANARY, SKEIN_GAP);
        gs_regions[i].base = p;
        gs_regions[i].usable = usable;
        gs_regions[i].size = usable + SKEIN_GAP;
        gs_regions[i].stack = stack;
        gs_regions[i].committed = p;
        gs_regions[i].reuse_lo = gs_regions[i].reuse_hi = NULL;
        gs_nregions = i + 1;
        return p;
    }

    """;

// Appends `s` as a JSON string, with `<` escaped so that it cannot end
// the script element.
fn json_string(out: u8[>..]&, s: u8[:]) {
    out.push('"');
    for c in s {
        if c == '"' || c == '\\' { out.push('\\'); out.push(c); }
        else if c == '\n' { out.append("\\n"); }
        else if c == '<' { out.append("\\u003c"); }
        else if c < 32 { out.append("\\u00"); out.push("0123456789abcdef"[c >> 4]); out.push("0123456789abcdef"[c & 15]); }
        else { out.push(c); }
    }
    out.push('"');
}

fn base64(out: u8[>..]&, s: u8[:]) {
    let digits = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    var i = 0;
    while i < s.len {
        let n = min(3, s.len - i);
        var v: i64 = 0;
        for k in n {
            let b: i64 = s[i + k];
            v |= b << (16 - 8 * k);
        }
        out.push(digits[(v >> 18) & 63]);
        out.push(digits[(v >> 12) & 63]);
        out.push(if n > 1 { digits[(v >> 6) & 63] } else { '=' });
        out.push(if n > 2 { digits[v & 63] } else { '=' });
        i += 3;
    }
}

// The solve of a page, for the page of `skein build` and for the report
// of `skein report` (src/report): the module `glue` of the model `model`
// in base64 (empty if the model is too large), the controls PARAMS of
// `params`, the memory MEMORY_MB of the run, `memory_mb`, the kernels
// `wgsl` and the WebGPU calls `gpu_js` of src/gpu/, and skeinSolve. It is
// two script elements.
fn solve_html(model: u8[:], params: Control[:], glue: u8[:], memory_mb: i64, wgsl: u8[:], gpu_js: u8[:]) -> u8[>..] {
    var out: u8[>..] = "<script id=\"skein-module\" type=\"application/octet-stream\">";
    base64(out, glue);
    out.append("</script>\n<script>\n\"use strict\";\nconst PARAMS = [");
    for k, i in params {
        out.append(if i > 0 { ", " } else { "" });
        out.append("{\"name\": ");
        json_string(out, k.name);
        out.append(", \"unit\": ");
        json_string(out, k.unit);
        format(out, ", \"lo\": ", k.lo, ", \"hi\": ", k.hi, ", \"value\": ", k.value, "}");
    }
    out.append("];\nconst RESULTS = ");
    json_string(out, str(".skein/results/", model));
    format(out, ";\nconst MEMORY_MB = ", memory_mb, ";\nconst MAX_MEMORY_MB = ", MAX_MEMORY_MB, ";\nconst WGSL = ");
    json_string(out, wgsl);
    out.append(";\nconst GPU_JS = ");
    json_string(out, gpu_js);
    out.append(";\n");
    out.append(SOLVE_JS);
    out.append("</script>\n");
    return out;
}

// The solve of the page (solve_html). skeinSolve(values, options) takes
// the value of each param (in its unit) and gives the run. Each run has a
// worker and a module of its own, so that a new value stops the run
// before it, and the memory of a run goes back when it ends. With WebGPU
// the worker gives the module a GPU device (src/gpu/gpu.js): the explicit
// steps always use it, and the implicit solves use it when
// options.implicitGpu is on.
let SOLVE_JS = """
    // The loader of the module, in UTF-8 bytes.
    const MODULE = Uint8Array.from(atob(document.getElementById("skein-module").textContent.trim()), (c) => c.charCodeAt(0));

    // The message of a model that is too large for the memory of the browser.
    function tooLarge() {
      return "This model is too large for the memory of the browser: it needs " + MEMORY_MB + " MB, and the page has " + MAX_MEMORY_MB + " MB.";
    }

    // The message of a run that failed: the memory limit for a run that
    // the browser cannot give the memory, else the last line of its log.
    // A run with no log and no result file says so.
    function failure(r) {
      if (/memory|OOM|RangeError|allocat/i.test(r.err)) {
        return "The browser could not give this run the memory that it can need, up to " + MEMORY_MB + " MB. Close other tabs or programs " +
          "and solve again, or run the model natively with skein run" + r.sets.map((w) => ' --set "' + w + '"').join("") + ".";
      }
      const last = r.err.trim().split("\n").pop();
      return "The solve in the page failed with exit " + r.code + (last ? ": " + last : Object.keys(r.files).length ? "." : ": no result file");
    }

    // The files of the records `results[n]:` of a run.
    function resultPaths(toon) {
      const lines = toon.split("\n"), paths = [];
      let i = lines.findIndex((l) => l.startsWith("results["));
      if (i < 0) return paths;
      for (i++; i < lines.length && lines[i].startsWith("  "); i++) paths.push(lines[i].trim());
      return paths;
    }

    // The worker of one run: it runs the module with the words of --set,
    // with a GPU device if e.data.gpu.
    function workerMain() {
      self.onmessage = async (e) => {
        const out = [], err = [];
        let mod = null, sent = false;
        const send = (code) => {
          if (sent) return;
          sent = true;
          const files = {};
          for (const p of resultPaths(out.join("\n"))) {
            try { files[p] = mod.FS.readFile(p, { encoding: "utf8" }); } catch (x) { err.push("cannot read " + p); }
          }
          postMessage({ code, out: out.join("\n"), err: err.join("\n"), files });
        };
        let gpu = null;
        if (e.data.gpu && self.navigator.gpu) {
          try {
            const adapter = await navigator.gpu.requestAdapter({ powerPreference: "high-performance" });
            if (adapter) gpu = skeinGpu(await adapter.requestDevice(), WGSL);
          } catch (x) { err.push("no GPU device: " + x); }
        }
        createSkein({
          arguments: e.data.args,
          skeinGpu: gpu,
          print: (s) => out.push(s),
          printErr: (s) => err.push(s),
          preRun: [(m) => {
            mod = m;
            m.ENV.SKEIN_THREADS = "1";
            if (e.data.implicitGpu) m.ENV.SKEIN_GPU_IMPLICIT = "1";
            m.FS.mkdirTree(e.data.results);
          }],
          onExit: (code) => send(code),
          onAbort: (why) => { err.push("abort: " + why); send(-1); },
        }).catch((x) => { err.push(String(x)); send(-1); });
      };
    }

    // WebGPU in the page: the worker asks for its own adapter, and a run
    // without one takes the CPU path of each step.
    const HAS_GPU = !!(self.navigator && navigator.gpu);

    let running = null;

    // The solve of the page: values holds a number for each param, in the
    // unit of the param (the default for a param that it does not hold).
    // It gives { code, out, err, files, seconds, sets }: code 0 is a run
    // whose checks pass, out its records in TOON, files its VTU files.
    // options.gpu (default: WebGPU is present) gives the run a GPU, and
    // options.implicitGpu puts the implicit solves on it too.
    function skeinSolve(values, options = {}) {
      if (running) { running.worker.terminate(); running.resolve({ code: -2, out: "", err: "stopped by a new run", files: {}, sets: running.sets }); }
      const sets = PARAMS.map((p) => p.name + "=" + (values[p.name] ?? p.value) + (p.unit ? " " + p.unit : ""));
      const source = [MODULE, "\nconst WGSL = " + JSON.stringify(WGSL) + ";\n" + GPU_JS + "\n" + resultPaths.toString() + "\n(" + workerMain.toString() + ")();"];
      const worker = new Worker(URL.createObjectURL(new Blob(source, { type: "text/javascript;charset=utf-8" })));
      const start = performance.now();
      return new Promise((resolve) => {
        running = { worker, resolve, sets };
        const end = (r) => {
          worker.terminate();
          if (running && running.worker === worker) running = null;
          resolve({ ...r, seconds: (performance.now() - start) / 1000, sets });
        };
        worker.onmessage = (e) => end(e.data);
        worker.onerror = (e) => end({ code: -1, out: "", err: String(e.message), files: {} });
        worker.postMessage({ args: sets, results: RESULTS, gpu: options.gpu ?? HAS_GPU, implicitGpu: !!options.implicitGpu });
      });
    }
    """;
