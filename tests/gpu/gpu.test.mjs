// Checks of the GPU path of the page of `skein build --target wasm`
// (src/gpu/) against the CPU path of the same page. Each test builds the
// page of a model, takes the module, the kernels and gpu.js out of it,
// and runs the module in Node twice: with no GPU, and with a WebGPU device
// of the webgpu package (Dawn). With no GPU adapter the tests skip, but
// SKEIN_GPU=1 makes that a failure. Run with `./build.sh gpu`.
import { test } from "node:test";
import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { createRequire } from "node:module";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const root = join(dirname(fileURLToPath(import.meta.url)), "..", "..");

// Dawn ends the device when the GC frees the GPU object, so `webgpu`
// keeps it. Each GPU run asks for a new adapter: Dawn crashed on the
// device of an adapter from the start of the file.
let webgpu = null, adapter = null, why = "";
try {
  const { create, globals } = await import("webgpu");
  Object.assign(globalThis, globals);
  webgpu = create([]);
  adapter = await webgpu.requestAdapter({ powerPreference: "high-performance" });
  if (!adapter) why = "no WebGPU adapter";
} catch (e) {
  why = `no WebGPU: ${e.message}`;
}
const skip = !adapter && !process.env.SKEIN_GPU && why;

// The module and the data of the page of `model`.
function page(model) {
  const dir = mkdtempSync(join(tmpdir(), "skein-gpu-"));
  execFileSync(join(root, "build", "skein"), ["build", join(root, model), "--target", "wasm"], { cwd: dir, stdio: "pipe" });
  const name = /^model ([A-Za-z0-9_]+):/m.exec(readFileSync(join(root, model), "utf8"))[1];
  const html = readFileSync(join(dir, ".skein", "build", name, "web", name + ".html"), "utf8");
  const data = html.slice(html.indexOf('"use strict";') + 13, html.indexOf("const RESULTS"));
  const p = new Function(data + "return { MODEL, PARAMS, WGSL, GPU_JS };")();
  const b64 = /<script id="skein-module"[^>]*>([^<]*)<\/script>/.exec(html)[1].trim();
  writeFileSync(join(dir, "skein.cjs"), Buffer.from(b64, "base64"));
  p.createSkein = createRequire(import.meta.url)(join(dir, "skein.cjs"));
  p.skeinGpu = new Function(p.GPU_JS + "\nreturn skeinGpu;")();
  rmSync(dir, { recursive: true });
  return p;
}

// One run of the page at the default values: on the GPU if `gpu`, with
// the implicit solves on it too if `implicit`. A run that takes more than
// `seconds` fails. Reason: a WebAssembly fault after a GPU read-back left
// a run of sample 01 waiting with no end.
async function run(p, gpu, implicit, seconds) {
  const out = [], err = [];
  let mod = null, dev = null;
  if (gpu) dev = p.skeinGpu(await (await webgpu.requestAdapter({ powerPreference: "high-performance" })).requestDevice(), p.WGSL);
  const keep = setInterval(() => {}, 1000);   // Dawn does not keep Node alive while the GPU works
  let timer = null;
  const code = await new Promise((resolve) => {
    timer = setTimeout(() => { err.push(`no exit in ${seconds} s`); resolve(-3); }, 1000 * seconds);
    p.createSkein({
      arguments: p.PARAMS.map((k) => k.name + "=" + k.value + (k.unit ? " " + k.unit : "")),
      skeinGpu: dev,
      print: (s) => out.push(s),
      printErr: (s) => err.push(s),
      preRun: [(m) => {
        mod = m;
        m.ENV.SKEIN_THREADS = "1";
        if (implicit) m.ENV.SKEIN_GPU_IMPLICIT = "1";
        m.FS.mkdirTree(".skein/results/" + p.MODEL);
      }],
      onExit: resolve,
      onAbort: (w) => { err.push("abort: " + w); resolve(-1); },
    }).catch((x) => { err.push(String(x)); resolve(-1); });
  });
  clearInterval(keep);
  clearTimeout(timer);
  const toon = out.join("\n");
  assert.equal(code, 0, toon + "\n" + err.join("\n"));
  if (dev) assert.deepEqual(dev.errors, []);
  const rows = (name) => {
    const l = toon.split("\n"), i = l.findIndex((s) => s.startsWith(name + "["));
    const r = [];
    for (let k = i + 1; i >= 0 && k < l.length && l[k].startsWith("  "); k++) r.push(l[k].trim().split(","));
    return r;
  };
  const vtu = mod.FS.readFile(rows("results")[0][0], { encoding: "utf8" });
  const a = vtu.indexOf("<DataArray", vtu.indexOf("<PointData>"));
  const body = vtu.slice(vtu.indexOf(">", a) + 1, vtu.indexOf("</DataArray>", a));
  return { devices: rows("device").map((r) => r[2]), field: Float64Array.from(body.trim().split(/\s+/), Number) };
}

// The largest difference of the two fields, relative to the largest value.
function difference(a, b) {
  assert.equal(a.length, b.length);
  let d = 0, most = 0;
  for (let i = 0; i < a.length; i++) { d = Math.max(d, Math.abs(a[i] - b[i])); most = Math.max(most, Math.abs(a[i])); }
  return d / most;
}

// The GPU works in f32, so the explicit steps agree with the CPU (f64) to
// the rounding of f32 over the steps. The implicit solve corrects its
// f32 passes in f64 to the tolerance of the CPU solve.
// `seconds` is the limit of each run, about 10 times the run on the
// desktop.
const cases = [
  { model: "examples/elements/bar_wave_tri.skn", implicit: false, tolerance: 1e-3, seconds: 60 },
  { model: "examples/grid/mms_fv_transient.skn", implicit: false, tolerance: 1e-3, seconds: 60 },
  { model: "examples/grid/mms_dg_transient.skn", implicit: false, tolerance: 1e-3, seconds: 60 },
  { model: "examples/elements/kirsch_plate_tri.skn", implicit: true, tolerance: 1e-6, seconds: 60 },
];

for (const c of cases) {
  test(`${c.model}: the GPU run matches the CPU run within ${c.tolerance}`, { skip }, async () => {
    assert.ok(adapter, why);
    const p = page(c.model);
    const cpu = await run(p, false, false, c.seconds);
    const gpu = await run(p, true, c.implicit, c.seconds);
    assert.ok(cpu.devices.length > 0 && cpu.devices.every((d) => d === "cpu"), cpu.devices.join());
    assert.ok(gpu.devices.length > 0 && gpu.devices.every((d) => d === "gpu"), gpu.devices.join());
    const d = difference(cpu.field, gpu.field);
    console.log(`${c.model}: ${gpu.devices.length} solve(s) on the GPU, difference ${d.toExponential(3)}`);
    assert.ok(d <= c.tolerance, `difference ${d}`);
  });
}
