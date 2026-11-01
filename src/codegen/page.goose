// The texts of `skein build --target wasm` (build.goose): the parts of the
// program of the page, the C of the profile and of the regions of
// WebAssembly, and the page.
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
// `--set`. web_quality prints the device record of src/gpu/ and the
// mesh_quality record of the mesh after the records.
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

// The page of `model`, with its controls `params`, the native run
// `native` of `seconds`, the module `glue` (empty if the model is too
// large), the memory of its regions, and the kernels and the WebGPU calls
// of src/gpu/.
fn page_html(model: u8[:], params: Control[:], native: u8[:], glue: u8[:], memory_mb: i64, seconds: f64, wgsl: u8[:], gpu_js: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    var data: u8[>..] = [];
    data.append("const MODEL = ");
    json_string(data, model);
    data.append(";\nconst PARAMS = [");
    for k, i in params {
        data.append(if i > 0 { ", " } else { "" });
        data.append("{\"name\": ");
        json_string(data, k.name);
        data.append(", \"unit\": ");
        json_string(data, k.unit);
        format(data, ", \"lo\": ", k.lo, ", \"hi\": ", k.hi, ", \"value\": ", k.value, "}");
    }
    data.append("];\nconst NATIVE = ");
    json_string(data, native);
    format(data, ";\nconst MEMORY_MB = ", memory_mb, ";\nconst MAX_MEMORY_MB = ", MAX_MEMORY_MB, ";\nconst NATIVE_SECONDS = ", seconds, ";\nconst WGSL = ");
    json_string(data, wgsl);
    data.append(";\nconst GPU_JS = ");
    json_string(data, gpu_js);
    data.append(";\n");
    var at = 0;
    while at < PAGE_HTML.len {
        if PAGE_HTML.len - at >= 6 && PAGE_HTML[at..at + 6] == "@DATA@" {
            out.append(data);
            at += 6;
        } else if PAGE_HTML.len - at >= 8 && PAGE_HTML[at..at + 8] == "@MODULE@" {
            base64(out, glue);
            at += 8;
        } else if PAGE_HTML.len - at >= 7 && PAGE_HTML[at..at + 7] == "@MODEL@" {
            out.append(model);
            at += 7;
        } else {
            out.push(PAGE_HTML[at]);
            at++;
        }
    }
    return out;
}

// The page. skeinSolve(values, options) is the solve of the page: it takes
// the value of each param (in its unit) and gives the run. Each run has a
// worker and a module of its own, so that a new value stops the run
// before it, and the memory of a run goes back when it ends. With WebGPU
// the worker gives the module a GPU device (src/gpu/gpu.js): the explicit
// steps always use it, and the implicit solves use it when the control
// "Implicit solves on the GPU" is on.
let PAGE_HTML = """
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>@MODEL@ in the browser</title>
    <style>
    :root { color-scheme: light dark; --bg: #ffffff; --fg: #1b1b1b; --muted: #606060; --line: #c8c8c8; --bad: #b00020; }
    @media (prefers-color-scheme: dark) { :root { --bg: #161616; --fg: #e8e8e8; --muted: #a0a0a0; --line: #444444; --bad: #ff6b6b; } }
    body { background: var(--bg); color: var(--fg); font: 15px/1.45 system-ui, sans-serif; margin: 0 auto; max-width: 960px; padding: 16px; }
    pre { overflow-x: auto; border: 1px solid var(--line); padding: 8px; font-size: 12px; }
    .param { display: grid; grid-template-columns: 5em 1fr 10em; gap: 8px; align-items: center; margin: 6px 0; }
    .param input[type=number] { width: 7em; }
    canvas { display: block; width: 100%; aspect-ratio: 2 / 1; border: 1px solid var(--line); }
    #state { color: var(--muted); }
    #error { color: var(--bad); font-weight: 600; }
    </style>
    </head>
    <body>
    <h1>@MODEL@</h1>
    <p id="note"></p>
    <div id="params"></div>
    <p><button id="solve">Solve in the page</button> <label><input type="checkbox" id="implicit-gpu"> Implicit solves on the GPU</label> <span id="state"></span></p>
    <p id="gpu-note"></p>
    <p id="error"></p>
    <canvas id="view" width="1200" height="600"></canvas>
    <p id="summary"></p>
    <table id="devices"></table>
    <h2>Run in the page</h2>
    <pre id="out">No run yet.</pre>
    <h2>Native run (skein run, default values)</h2>
    <pre id="native"></pre>
    <script id="skein-module" type="application/octet-stream">@MODULE@</script>
    <script>
    "use strict";
    @DATA@
    const RESULTS = ".skein/results/" + MODEL;
    // The loader of the module, in UTF-8 bytes.
    const MODULE = Uint8Array.from(atob(document.getElementById("skein-module").textContent.trim()), (c) => c.charCodeAt(0));

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

    // The arrays of a VTU file of skein (ascii).
    function parseVtu(text) {
      const after = (tag, from) => {
        const a = text.indexOf(tag, from);
        const b = text.indexOf(">", a) + 1;
        return { body: text.slice(b, text.indexOf("</DataArray>", b)), head: text.slice(a, b), end: b };
      };
      const nums = (s) => Float64Array.from(s.trim().split(/\s+/), Number);
      const pts = after("<DataArray", text.indexOf("<Points>"));
      const conn = after('Name="connectivity"', 0), offs = after('Name="offsets"', 0), types = after('Name="types"', 0);
      const field = after("<DataArray", text.indexOf("<PointData>"));
      const comps = Number((/NumberOfComponents="(\d+)"/.exec(field.head) || [0, 1])[1]);
      return { xyz: nums(pts.body), conn: nums(conn.body), offs: nums(offs.body), types: nums(types.body), field: nums(field.body), comps };
    }

    // The boundary triangles of the corner nodes of the cells.
    function surface(v) {
      const tris = [], count = new Map();
      for (let c = 0; c < v.types.length; c++) {
        const o = c > 0 ? v.offs[c - 1] : 0, t = v.types[c];
        const n = Array.from(v.conn.slice(o, o + 4));
        if (t === 5 || t === 22) { tris.push([n[0], n[1], n[2]]); continue; }
        if (t === 9 || t === 23 || t === 28) { tris.push([n[0], n[1], n[2]], [n[0], n[2], n[3]]); continue; }
        for (const f of [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]) {
          const tri = [n[f[0]], n[f[1]], n[f[2]]], key = tri.slice().sort((a, b) => a - b).join(",");
          count.set(key, count.has(key) ? null : tri);
        }
      }
      for (const tri of count.values()) if (tri) tris.push(tri);
      return tris;
    }

    // Draws the deformed mesh, coloured by the size of the field.
    function draw(v) {
      const cv = document.getElementById("view"), g = cv.getContext("2d");
      const nn = v.xyz.length / 3, size = new Float64Array(nn);
      let most = 0;
      for (let i = 0; i < nn; i++) {
        let s = 0;
        for (let a = 0; a < v.comps; a++) s += v.field[v.comps * i + a] ** 2;
        size[i] = Math.sqrt(s);
        most = Math.max(most, size[i]);
      }
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      for (let i = 0; i < nn; i++) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], v.xyz[3 * i + a]); hi[a] = Math.max(hi[a], v.xyz[3 * i + a]); }
      const diag = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
      const scale = v.comps === 3 && most > 0 ? 0.05 * diag / most : 0;
      const yaw = -0.6, pitch = 0.45, cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
      const px = new Float64Array(nn), py = new Float64Array(nn), pz = new Float64Array(nn);
      for (let i = 0; i < nn; i++) {
        const p = [0, 1, 2].map((a) => v.xyz[3 * i + a] - 0.5 * (lo[a] + hi[a]) + (v.comps === 3 ? scale * v.field[3 * i + a] : 0));
        const x = cy * p[0] - sy * p[1], y = sy * p[0] + cy * p[1];
        px[i] = x; py[i] = cp * p[2] - sp * y; pz[i] = sp * p[2] + cp * y;
      }
      const k = 0.85 * Math.min(cv.width, 2 * cv.height) / diag;
      const tris = surface(v).map((t) => [t, pz[t[0]] + pz[t[1]] + pz[t[2]]]).sort((a, b) => b[1] - a[1]);
      g.clearRect(0, 0, cv.width, cv.height);
      g.lineWidth = 0.4;
      g.strokeStyle = "rgba(0, 0, 0, 0.45)";
      for (const [t] of tris) {
        const f = most > 0 ? (size[t[0]] + size[t[1]] + size[t[2]]) / (3 * most) : 0;
        g.fillStyle = "hsl(" + Math.round(240 * (1 - f)) + ", 70%, 55%)";
        g.beginPath();
        g.moveTo(cv.width / 2 + k * px[t[0]], cv.height / 2 - k * py[t[0]]);
        g.lineTo(cv.width / 2 + k * px[t[1]], cv.height / 2 - k * py[t[1]]);
        g.lineTo(cv.width / 2 + k * px[t[2]], cv.height / 2 - k * py[t[2]]);
        g.closePath();
        g.fill();
        g.stroke();
      }
      return { nodes: nn, cells: v.types.length, most, scale };
    }

    // The rows of the record `name[n]{...}:` of a run.
    function recordRows(toon, name) {
      const lines = toon.split("\n"), rows = [];
      let i = lines.findIndex((l) => l.startsWith(name + "["));
      if (i < 0) return rows;
      for (i++; i < lines.length && lines[i].startsWith("  "); i++) rows.push(lines[i].trim().split(","));
      return rows;
    }

    // The table of the device of each step of a run.
    function showDevices(toon) {
      const t = document.getElementById("devices");
      t.textContent = "";
      const rows = recordRows(toon, "device");
      if (!rows.length) return;
      for (const r of [["step", "solver", "device"], ...rows]) {
        const tr = t.insertRow();
        for (const c of r) tr.insertCell().textContent = c;
      }
    }

    function solveRow(toon) {
      const lines = toon.split("\n"), i = lines.findIndex((l) => l.startsWith("solve["));
      return i >= 0 && i + 1 < lines.length ? lines[i + 1].trim().split(",") : null;
    }

    const state = document.getElementById("state"), error = document.getElementById("error");

    // The message of a run that failed.
    function failure(r) {
      if (/memory|OOM|RangeError|allocat/i.test(r.err)) {
        return "The browser could not give this run the memory that it can need, up to " + MEMORY_MB + " MB. Close other tabs or programs " +
          "and solve again, or run the model natively with" + r.sets.map((w) => ' --set "' + w + '"').join("") + ".";
      }
      return "The run failed with exit " + r.code + ". The log below gives the cause.";
    }
    const values = {};
    for (const p of PARAMS) {
      values[p.name] = p.value;
      const row = document.createElement("div");
      row.className = "param";
      row.innerHTML = "<label></label><input type=range step=any><span><input type=number step=any> </span>";
      row.querySelector("label").textContent = p.name;
      row.querySelector("span").append(p.unit);
      const range = row.querySelector("input[type=range]"), num = row.querySelector("input[type=number]");
      for (const e of [range, num]) { e.min = p.lo; e.max = p.hi; e.value = p.value; }
      range.oninput = () => { num.value = range.value; };
      const set = (x) => { values[p.name] = Math.min(p.hi, Math.max(p.lo, Number(x))); range.value = num.value = values[p.name]; run(); };
      range.onchange = () => set(range.value);
      num.onchange = () => set(num.value);
      document.getElementById("params").append(row);
    }

    async function run() {
      if (!MODULE.length) return;
      state.textContent = "meshing and solving in the page...";
      error.textContent = "";
      const r = await skeinSolve(values, { implicitGpu: document.getElementById("implicit-gpu").checked });
      if (r.code === -2) return;
      showDevices(r.out);
      document.getElementById("out").textContent = r.sets.join("  ") + "\n\n" + r.out + (r.err ? "\n" + r.err : "");
      state.textContent = "exit " + r.code + " in " + r.seconds.toFixed(1) + " s";
      if (r.code !== 0) error.textContent = failure(r);
      const vtu = Object.values(r.files)[0];
      if (vtu) {
        const d = draw(parseVtu(vtu));
        document.getElementById("summary").textContent = d.cells + " cells, " + d.nodes + " nodes, largest field size " + d.most +
          (d.scale ? ", displacement drawn " + d.scale.toPrecision(3) + " times" : "");
      }
    }

    // WebGPU in the page: the worker asks for its own adapter, and a run
    // without one takes the CPU path of each step.
    const HAS_GPU = !!(self.navigator && navigator.gpu);
    document.getElementById("implicit-gpu").disabled = !HAS_GPU;
    document.getElementById("gpu-note").textContent = HAS_GPU ?
      "WebGPU is present: the explicit steps run on the GPU, in f32. The implicit solves run on the CPU, unless you turn on the GPU for them." :
      "This browser has no WebGPU: all steps run on the CPU.";
    document.getElementById("native").textContent = NATIVE;
    document.getElementById("solve").onclick = run;
    const row = solveRow(NATIVE);
    if (!MODULE.length) {
      document.getElementById("solve").disabled = true;
      document.getElementById("note").textContent = "This model is too large for the memory of the browser: it needs " +
        MEMORY_MB + " MB, and the page has " + MAX_MEMORY_MB + " MB." + (row ? " The native run has " + row[1] + " cells and " + row[2] +
        " dofs, and took " + NATIVE_SECONDS.toFixed(1) + " s." : "") + " The page shows the native results.";
    } else {
      document.getElementById("note").textContent = "Each control meshes and solves the model again in the page (WebAssembly, one thread, up to " +
        MEMORY_MB + " MB of memory). The native run at the default values took " + NATIVE_SECONDS.toFixed(1) + " s.";
    }
    </script>
    </body>
    </html>
    """;
