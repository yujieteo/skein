// Checks of the explicit solves in the page of `skein build --target wasm`,
// in a real headless Chromium with WebGPU (not Node with Dawn). For each
// model the script opens the page, solves it twice with skeinSolve of the
// page (once on the CPU, once with the GPU), and checks:
//   - the GPU run reports only device "gpu", and the CPU run only "cpu"
//     (a GPU run that falls back to the CPU is a failure);
//   - the page checks the expect values of the model (verify rows) on both
//     runs, and each of them passes;
// The reference (below, with its source) only feeds the error column of
// the table. The page verdict checks the expect values.
// Chromium is driven over the DevTools protocol (global WebSocket of Node).
// The WebGPU adapter must match SKEIN_ADAPTER (default "nvidia").
//
//   node tests/gpu/browser.mjs [model-substring ...]
//
// Environment:
//   SKEIN_CHROMIUM  the browser (default chromium)
//   SKEIN_ADAPTER   regular expression for the adapter (default nvidia)
//   SKEIN_PAGES     directory that keeps the built pages (default a new
//                   temporary directory); pages in it are not built again
//   SKEIN_OUT       directory for table.md and table.json (default: print)
// On NixOS, Chromium and Dawn find the Nvidia driver only with
// LD_LIBRARY_PATH=<vulkan-loader>/lib:/run/opengl-driver/lib and
// VK_ICD_FILENAMES=/run/opengl-driver/share/vulkan/icd.d/nvidia_icd.json.
import { execFile, spawn } from "node:child_process";
import { existsSync, mkdirSync, mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { promisify } from "node:util";
import { fileURLToPath } from "node:url";

const root = join(dirname(fileURLToPath(import.meta.url)), "..", "..");
const run = promisify(execFile);

// ref: the reference of each expect id, in the unit of the page value.
// pick: "expect" gives the value of the verify row, "report" the value of
// the first report row (for an expect of a field).
// seconds: the limit of one run (the CPU run of the wasm module is slow).
const cases = [
  // Graff 1975, 1D wave theory: tip = S*L/E = 100 MPa*100 mm/200 GPa.
  { model: "examples/elements/bar_wave_tri.skn", seconds: 120, refs: { tip_u: 0.05 } },
  // -V0*t_end = -227 m/s*5 us.
  { model: "examples/elements/taylor_bar_tet.skn", seconds: 900, refs: { top_free: -1.135 } },
  // Half of the static value: 0.5*S*L/E.
  { model: "examples/elements/bar_contact_tri.skn", seconds: 120, refs: { mid_b: 0.025 } },
  // The page gives the L2 error against the manufactured solution
  // exp(-t/T) sin sin u0, so the reference is 0 and the error is absolute.
  { model: "examples/grid/mms_fv_transient.skn", seconds: 120, refs: { mms: 0 } },
  { model: "examples/grid/mms_dg_transient.skn", seconds: 120, refs: { mms: 0 } },
  // femahoy validation/cases/dynamics-bar-step: 2*F*L/(E*A), Graff 1975.
  { model: "examples/femahoy/dynamics_bar_step.skn", seconds: 300, refs: { u_end: 0.019047619047619049 } },
  // femahoy validation/cases/taylor-bar-impact, DYNA3D/ccx reference.
  { model: "examples/femahoy/taylor_bar_impact.skn", seconds: 1800, refs: { length: -10.74, radius: 3.82 } },
];

const only = process.argv.slice(2);
const chosen = cases.filter((c) => only.length === 0 || only.some((o) => c.model.includes(o)));

// --- the pages ---------------------------------------------------------

const pages = process.env.SKEIN_PAGES || mkdtempSync(join(tmpdir(), "skein-browser-"));
mkdirSync(pages, { recursive: true });
const nameOf = (c) => /^model ([A-Za-z0-9_]+):/m.exec(readFileSync(join(root, c.model), "utf8"))[1];

async function build(c) {
  const name = nameOf(c), dir = join(pages, name);
  const html = join(dir, ".skein", "build", name, "web", name + ".html");
  if (!existsSync(html)) {
    mkdirSync(dir, { recursive: true });
    await run(join(root, "build", "skein"), ["build", join(root, c.model), "--target", "wasm"], { cwd: dir, maxBuffer: 1 << 28 });
  }
  return html;
}

// --- the DevTools protocol ---------------------------------------------

const userData = mkdtempSync(join(tmpdir(), "skein-chromium-"));
const chromium = spawn(process.env.SKEIN_CHROMIUM || "chromium", [
  "--headless=new", "--remote-debugging-port=0", `--user-data-dir=${userData}`,
  "--enable-unsafe-webgpu", "--enable-features=Vulkan", "--ignore-gpu-blocklist",
  "--use-angle=vulkan", "--disable-vulkan-surface", "--no-sandbox", "about:blank",
], { stdio: ["ignore", "ignore", "ignore"] });
const stop = () => { try { chromium.kill(); } catch {} };
process.on("exit", stop);

let ws = null, nextId = 1;
const waiting = new Map();
async function connect() {
  const file = join(userData, "DevToolsActivePort");
  for (let i = 0; i < 100 && !existsSync(file); i++) await new Promise((r) => setTimeout(r, 100));
  const [port, path] = readFileSync(file, "utf8").trim().split("\n");
  ws = new WebSocket(`ws://127.0.0.1:${port}${path}`);
  await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
  ws.onmessage = (m) => {
    const d = JSON.parse(m.data), w = waiting.get(d.id);
    if (w) { waiting.delete(d.id); d.error ? w.reject(new Error(JSON.stringify(d.error))) : w.resolve(d.result); }
  };
}
function cdp(method, params = {}, sessionId) {
  const id = nextId++;
  return new Promise((resolve, reject) => {
    waiting.set(id, { resolve, reject });
    ws.send(JSON.stringify({ id, method, params, sessionId }));
  });
}

// Evaluates `expression` in the page of `sessionId`; waits at most `seconds`.
async function evaluate(sessionId, expression, seconds) {
  let timer;
  const limit = new Promise((_, reject) => { timer = setTimeout(() => reject(new Error(`no result in ${seconds} s`)), 1000 * seconds); });
  try {
    const r = await Promise.race([cdp("Runtime.evaluate", { expression, awaitPromise: true, returnByValue: true }, sessionId), limit]);
    if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description || r.exceptionDetails.text);
    return r.result.value;
  } finally { clearTimeout(timer); }
}

// --- the records of a run ----------------------------------------------

function csv(line) {
  const f = [];
  for (const m of line.matchAll(/"([^"]*)"|([^,]+)/g)) f.push(m[1] ?? m[2]);
  return f;
}
function rows(toon, name) {
  const l = toon.split("\n"), r = [];
  for (let i = 0; i < l.length; i++) {
    if (!l[i].startsWith(name + "[")) continue;
    for (let k = i + 1; k < l.length && l[k].startsWith("  "); k++) r.push(csv(l[k].trim()));
  }
  return r;
}
// The post expect rows: id -> { result, value (number), unit, limit }.
function expects(toon) {
  const out = {};
  for (const r of rows(toon, "verify")) {
    if (r[0] !== "post") continue;
    const m = /^([A-Za-z0-9_]+): (?:L2 error )?(\S+)(?: (.*))?$/.exec(r[3]);
    if (m) out[m[1]] = { result: r[2], value: Number(m[2]), unit: m[3] ?? "", limit: r[4] };
  }
  return out;
}
const devices = (toon) => rows(toon, "device").map((r) => r[2]);

// --- the run -----------------------------------------------------------

const results = [];
let failures = 0;
function fail(row, why) { row.failures.push(why); failures++; }

async function main() {
  console.error(`building ${chosen.length} page(s) in ${pages}`);
  const queue = [...chosen], htmls = new Map();
  await Promise.all([0, 1, 2].map(async () => {
    for (let c; (c = queue.shift()); ) htmls.set(c.model, await build(c));
  }));
  await connect();
  const version = (await cdp("Browser.getVersion")).product;
  let adapter = null;

  for (const c of chosen) {
    const row = { model: c.model, failures: [], runs: {} };
    results.push(row);
    const { targetId } = await cdp("Target.createTarget", { url: "file://" + htmls.get(c.model) });
    const { sessionId } = await cdp("Target.attachToTarget", { targetId, flatten: true });
    try {
      for (let i = 0; i < 300 && !(await evaluate(sessionId, "typeof skeinSolve === 'function'", 10)); i++) await new Promise((r) => setTimeout(r, 100));
      if (!adapter) {
        // The GPU process may need a moment to start.
        for (let i = 0; i < 20 && !adapter; i++) {
          adapter = await evaluate(sessionId, `(async () => {
            const a = navigator.gpu && await navigator.gpu.requestAdapter({ powerPreference: "high-performance" });
            if (!a) return null;
            const i = a.info;
            return { vendor: i.vendor, architecture: i.architecture, device: i.device, description: i.description, fallback: a.isFallbackAdapter === true };
          })()`, 60);
          if (!adapter) await new Promise((r) => setTimeout(r, 500));
        }
        // The name of the GPU, from the GPU process of the browser (the
        // adapter info of WebGPU has no description in this Chromium).
        if (adapter) adapter.name = (await cdp("SystemInfo.getInfo")).gpu.devices.map((d) => d.deviceString).join("; ");
        const text = adapter ? JSON.stringify(adapter) : "none";
        console.error(`WebGPU adapter: ${text}`);
        if (!adapter || adapter.fallback || !new RegExp(process.env.SKEIN_ADAPTER || "nvidia", "i").test(text)) {
          throw new Error(`the WebGPU adapter is not ${process.env.SKEIN_ADAPTER || "nvidia"}: ${text}`);
        }
      }
      for (const mode of ["cpu", "gpu"]) {
        const r = await evaluate(sessionId, `(async () => {
          const r = await skeinSolve({}, { gpu: ${mode === "gpu"} });
          return { code: r.code, out: r.out, err: r.err, seconds: r.seconds };
        })()`, c.seconds);
        const dev = devices(r.out);
        row.runs[mode] = { code: r.code, seconds: r.seconds, devices: [...new Set(dev)], expects: expects(r.out), err: r.err };
        if (r.code !== 0) fail(row, `${mode}: exit ${r.code}: ${r.err.trim().split("\n").pop()}`);
        if (dev.length === 0 || dev.some((d) => d !== mode)) fail(row, `${mode}: devices ${dev.join(",") || "none"}`);
        if (/no GPU device/.test(r.err)) fail(row, `${mode}: ${r.err}`);
      }
      for (const [id, ref] of Object.entries(c.refs)) {
        for (const mode of ["cpu", "gpu"]) {
          const e = row.runs[mode].expects[id];
          if (!e) { fail(row, `${mode}: no expect ${id} in the page`); continue; }
          if (e.result !== "pass") fail(row, `${mode}: expect ${id} ${e.result}`);
        }
      }
    } catch (x) {
      fail(row, String(x.message || x));
    } finally {
      await cdp("Target.closeTarget", { targetId }).catch(() => {});
    }
    console.error(`${row.failures.length ? "FAIL" : "ok  "} ${c.model} ${row.failures.join("; ")}`);
    if (!adapter) break;
  }

  const table = render(version, adapter, chosen, results);
  console.log(table);
  if (process.env.SKEIN_OUT) {
    mkdirSync(process.env.SKEIN_OUT, { recursive: true });
    writeFileSync(join(process.env.SKEIN_OUT, "table.json"), JSON.stringify({ version, adapter, results }, null, 1));
    writeFileSync(join(process.env.SKEIN_OUT, "table.md"), table);
  }
}

const g = (v) => (Number.isFinite(v) ? Number(v.toPrecision(6)).toString() : "-");
function render(version, adapter, cs, rs) {
  const l = [];
  l.push(`Browser: ${version}`);
  l.push(`WebGPU adapter: ${adapter ? `${adapter.vendor} ${adapter.architecture} ${adapter.device} ${adapter.description} (${adapter.name})`.replace(/\s+/g, " ").trim() : "none"}`, "");
  l.push("| model | expect | unit | reference | GPU value | CPU value | GPU error | CPU error | page GPU | page CPU | GPU s | CPU s | result |");
  l.push("|---|---|---|---|---|---|---|---|---|---|---|---|---|");
  for (const [i, r] of rs.entries()) {
    const ids = Object.keys(cs[i].refs);
    for (const id of ids) {
      const gp = r.runs.gpu?.expects[id], cp = r.runs.cpu?.expects[id], ref = cs[i].refs[id];
      const err = (e) => (!e ? "-" : ref === 0 ? `${g(Math.abs(e.value))} ${e.unit} (abs)` : `${g((100 * Math.abs(e.value - ref)) / Math.abs(ref))} %`);
      const secs = (m) => (r.runs[m] ? r.runs[m].seconds.toFixed(1) : "-");
      l.push(`| ${r.model.replace(/^examples\//, "")} | ${id} | ${gp?.unit || cp?.unit || ""} | ${g(ref)} | ${g(gp?.value)} | ${g(cp?.value)} | ${err(gp)} | ${err(cp)} | ${gp?.result ?? "-"} | ${cp?.result ?? "-"} | ${secs("gpu")} | ${secs("cpu")} | ${r.failures.length ? "FAIL: " + r.failures.join("; ").replace(/\|/g, "/") : "pass"} |`);
    }
  }
  return l.join("\n") + "\n";
}

try {
  await main();
} catch (x) {
  console.error(String(x.stack || x));
  failures++;
} finally {
  stop();
}
process.exit(failures ? 1 : 0);
