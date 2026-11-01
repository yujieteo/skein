// The viewer of `skein report` (src/report/report.goose). Plain JavaScript
// and WebGL 2, with no library, so that the page opens offline.
//
// The page holds the data in the script element #skein-data: the source of
// the model, the TOON of check, verify and run, and for each case the run
// TOON and its result files. Each result file is a small VTU, gzip in
// base64, in the script element #vtu-<n>. The viewer also opens a VTU file
// of another solver (the button "Open VTU"). That view has no source links.
//
// The page can also hold the solve of `skein build --target wasm`
// (solve_html of src/codegen/page.goose): PARAMS, MODULE and skeinSolve.
// Then each param of PARAMS has a slider over its range. A slider at the
// values of a stored case shows that case at once. A slider that stops at
// other values meshes and solves in the page, and the solve becomes a case.
'use strict';

// --- TOON --------------------------------------------------------------------------

// The values of a TOON row. A value in quotes can hold commas and \" escapes.
function toonRow(s) {
  const out = [];
  let i = 0;
  while (i <= s.length) {
    if (s[i] === '"') {
      let v = '';
      i++;
      while (i < s.length && s[i] !== '"') {
        if (s[i] === '\\' && i + 1 < s.length) i++;
        v += s[i++];
      }
      out.push(v);
      i = s.indexOf(',', i);
      if (i < 0) break;
      i++;
    } else {
      const j = s.indexOf(',', i);
      out.push(s.slice(i, j < 0 ? s.length : j).trim());
      if (j < 0) break;
      i = j + 1;
    }
  }
  return out;
}

// The tables (name[n]{cols}:), lists (name[n]:) and values (key: value) of
// a TOON text.
function parseToon(text) {
  const out = { tables: {}, lists: {}, values: {} };
  const lines = (text || '').split('\n');
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    let m = /^(\w+)\[(\d+)\]\{(.*)\}:\s*$/.exec(line);
    if (m) {
      const rows = [];
      for (let k = 0; k < +m[2] && i + 1 < lines.length; k++) rows.push(toonRow(lines[++i].trim()));
      out.tables[m[1]] = { cols: toonRow(m[3]), rows };
      continue;
    }
    m = /^(\w+)\[(\d+)\]:\s*$/.exec(line);
    if (m) {
      const rows = [];
      for (let k = 0; k < +m[2] && i + 1 < lines.length; k++) rows.push(lines[++i].trim());
      out.lists[m[1]] = rows;
      continue;
    }
    m = /^(\w+):\s*(.*)$/.exec(line);
    if (m) out.values[m[1]] = m[2];
  }
  return out;
}

// --- VTU ---------------------------------------------------------------------------

const TYPED = {
  Float32: Float32Array, Float64: Float64Array, Int8: Int8Array, UInt8: Uint8Array,
  Int16: Int16Array, UInt16: Uint16Array, Int32: Int32Array, UInt32: Uint32Array,
  Int64: BigInt64Array, UInt64: BigUint64Array,
};

async function inflate(bytes, format) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream(format));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

function fromBase64(s) {
  const bin = atob(s.replace(/\s+/g, ''));
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

function attr(tag, key) {
  const m = new RegExp('\\s' + key + '="([^"]*)"').exec(tag);
  return m ? m[1] : null;
}

// The values of a binary block at `at` of `bytes`: a byte count (or the
// block table of vtkZLibDataCompressor) in header ints, then the data.
async function binaryBlock(bytes, at, head, zlib) {
  const hs = head === 'UInt64' ? 8 : 4;
  const dv = new DataView(bytes.buffer, bytes.byteOffset);
  const int = (k) => (hs === 8 ? Number(dv.getBigUint64(at + k * 8, true)) : dv.getUint32(at + k * 4, true));
  if (!zlib) {
    const n = int(0);
    return bytes.slice(at + hs, at + hs + n);
  }
  const blocks = int(0);
  let from = at + hs * (3 + blocks);
  const parts = [];
  for (let b = 0; b < blocks; b++) {
    const size = int(3 + b);
    parts.push(await inflate(bytes.slice(from, from + size), 'deflate'));
    from += size;
  }
  const out = new Uint8Array(parts.reduce((s, p) => s + p.length, 0));
  let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}

// The arrays of a VTU file (ASCII, binary or raw appended data, with or
// without zlib): points, cells, and the point and cell data.
// `part` allows a file with no points or no cells: the fields part of a
// result of a report.
async function parseVtu(bytes, part = false) {
  const text = new TextDecoder('latin1').decode(bytes);
  const vtk = /<VTKFile[^>]*>/.exec(text);
  if (!vtk || attr(vtk[0], 'type') !== 'UnstructuredGrid') throw new Error('not a VTU file of an UnstructuredGrid');
  if ((attr(vtk[0], 'byte_order') || 'LittleEndian') !== 'LittleEndian') throw new Error('a big endian VTU file');
  const head = attr(vtk[0], 'header_type') || 'UInt32';
  const zlib = (attr(vtk[0], 'compressor') || '').includes('ZLib');
  const app = text.indexOf('<AppendedData');
  const appStart = app < 0 ? -1 : text.indexOf('_', app) + 1;
  const appBase64 = app >= 0 && attr(text.slice(app, text.indexOf('>', app) + 1), 'encoding') === 'base64';
  const xml = app < 0 ? text : text.slice(0, app);
  const out = { points: null, conn: null, offsets: null, types: null, pointData: [], cellData: [] };
  const re = /<(Points|Cells|PointData|CellData)[\s>]|<DataArray[^>]*?(\/?)>/g;
  let section = '';
  let m;
  while ((m = re.exec(xml))) {
    if (m[1]) { section = m[1]; continue; }
    const tag = m[0];
    const Type = TYPED[attr(tag, 'type')];
    if (!Type) continue;
    const format = attr(tag, 'format') || 'ascii';
    let values;
    if (format === 'ascii') {
      const end = xml.indexOf('</DataArray>', re.lastIndex);
      values = xml.slice(re.lastIndex, end).trim().split(/\s+/).filter((s) => s.length).map(Number);
    } else {
      let raw;
      if (format === 'binary') {
        const end = xml.indexOf('</DataArray>', re.lastIndex);
        raw = await binaryBlock(fromBase64(xml.slice(re.lastIndex, end)), 0, head, zlib);
      } else if (appBase64) {
        const from = appStart + +attr(tag, 'offset');
        const end = text.slice(from).search(/[^A-Za-z0-9+/=]/);
        raw = await binaryBlock(fromBase64(text.slice(from, end < 0 ? text.length : from + end)), 0, head, zlib);
      } else {
        raw = await binaryBlock(bytes, appStart + +attr(tag, 'offset'), head, zlib);
      }
      const typed = new Type(raw.buffer, raw.byteOffset, raw.byteLength / Type.BYTES_PER_ELEMENT);
      values = Type === BigInt64Array || Type === BigUint64Array ? Array.from(typed, Number) : typed;
    }
    const name = attr(tag, 'Name') || '';
    const comps = +(attr(tag, 'NumberOfComponents') || 1);
    const scale = +(attr(tag, 'Scale') || 1);
    if (section === 'Points') out.points = Float64Array.from(values);
    else if (section === 'Cells' && name === 'connectivity') out.conn = Int32Array.from(values);
    else if (section === 'Cells' && name === 'offsets') out.offsets = Int32Array.from(values);
    else if (section === 'Cells' && name === 'types') out.types = Uint8Array.from(values);
    else if (section === 'PointData') out.pointData.push({ name, comps, data: Float64Array.from(values, (v) => v * scale) });
    else if (section === 'CellData') out.cellData.push({ name, comps, data: Float64Array.from(values, (v) => v * scale) });
  }
  if (!part && (!out.points || !out.conn || !out.offsets || !out.types)) throw new Error('the VTU file has no points or no cells');
  return out;
}

// --- the mesh ----------------------------------------------------------------------

// The linear type and the corner count of each VTK cell type.
const LINEAR = { 21: [3, 2], 22: [5, 3], 23: [9, 4], 28: [9, 4], 24: [10, 4], 25: [12, 8], 29: [12, 8], 26: [13, 6] };
const FACES = {
  10: [[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]],
  12: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
  13: [[0, 1, 2], [3, 5, 4], [0, 3, 4, 1], [1, 4, 5, 2], [2, 5, 3, 0]],
  14: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
};
const TETS = {
  10: [[0, 1, 2, 3]],
  12: [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]],
  13: [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]],
  14: [[0, 1, 2, 4], [0, 2, 3, 4]],
};

// The mesh of a parsed VTU: the corner nodes of each cell, the boundary
// faces as triangles, and the tetrahedra for the cut plane.
function buildMesh(vtu) {
  const cells = [];
  let from = 0;
  for (let c = 0; c < vtu.offsets.length; c++) {
    const to = vtu.offsets[c];
    const lin = LINEAR[vtu.types[c]] || [vtu.types[c], to - from];
    cells.push({ type: lin[0], nodes: Array.from(vtu.conn.slice(from, from + Math.min(lin[1], to - from))) });
    from = to;
  }
  const faces = new Map();
  const tris = [];      // node triples
  const triCell = [];   // the cell of each triangle
  const addFace = (nodes, c) => {
    for (let k = 1; k + 1 < nodes.length; k++) { tris.push(nodes[0], nodes[k], nodes[k + 1]); triCell.push(c); }
  };
  cells.forEach((cell, c) => {
    const fs = FACES[cell.type];
    if (!fs) {
      if (cell.type === 5 || cell.type === 9) addFace(cell.nodes, c);
      return;
    }
    for (const f of fs) {
      const nodes = f.map((k) => cell.nodes[k]);
      const key = nodes.slice().sort((a, b) => a - b).join(',');
      if (faces.has(key)) faces.delete(key); else faces.set(key, [nodes, c]);
    }
  });
  const boundary = [];
  cells.forEach((cell, c) => { if (cell.type === 5 || cell.type === 9) boundary.push(cell.nodes); });
  for (const [nodes, c] of faces.values()) { addFace(nodes, c); boundary.push(nodes); }
  const p = vtu.points;
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < p.length; i += 3) {
    for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[i + k]); hi[k] = Math.max(hi[k], p[i + k]); }
  }
  const size = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  // the line cells (VTK 3): beams
  const beams = [];
  cells.forEach((cell, c) => { if (cell.type === 3) beams.push(c); });
  return { vtu, cells, faces: boundary, beams, tris: Int32Array.from(tris), triCell: Int32Array.from(triCell), lo, hi, size,
    fields: fieldsOf(vtu), disp: displacement(vtu), regions: new Map() };
}

// The fields of a VTU that the viewer shows: each scalar, and the size and
// the components of each vector and tensor. von_mises comes from a tensor
// of 6 or 9 components.
function fieldsOf(vtu) {
  const out = [];
  const add = (arr, cell) => {
    const n = arr.data.length / arr.comps;
    const comp = (k) => () => { const v = new Float64Array(n); for (let i = 0; i < n; i++) v[i] = arr.data[i * arr.comps + k]; return v; };
    if (arr.comps === 1) { out.push({ name: arr.name, base: arr.name, cell, values: () => arr.data }); return; }
    if (arr.comps === 3) {
      out.push({ name: '|' + arr.name + '|', base: arr.name, cell, values: () => {
        const v = new Float64Array(n);
        for (let i = 0; i < n; i++) v[i] = Math.hypot(arr.data[3 * i], arr.data[3 * i + 1], arr.data[3 * i + 2]);
        return v;
      } });
      ['x', 'y', 'z'].forEach((c, k) => out.push({ name: arr.name + '.' + c, base: arr.name, cell, values: comp(k) }));
      return;
    }
    if (arr.comps === 6 || arr.comps === 9) {
      const idx = arr.comps === 6 ? [0, 1, 2, 4, 5, 3] : [0, 4, 8, 5, 2, 1];   // xx yy zz yz xz xy (VTK order of 6 is xx yy zz xy yz xz)
      out.push({ name: 'von_mises(' + arr.name + ')', base: arr.name, cell, values: () => {
        const v = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          const s = idx.map((k) => arr.data[i * arr.comps + k]);
          v[i] = Math.sqrt(0.5 * ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) + 3 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2));
        }
        return v;
      } });
      ['xx', 'yy', 'zz', 'yz', 'xz', 'xy'].forEach((c, k) => out.push({ name: arr.name + '.' + c, base: arr.name, cell, values: comp(idx[k]) }));
      return;
    }
    for (let k = 0; k < arr.comps; k++) out.push({ name: arr.name + '[' + k + ']', base: arr.name, cell, values: comp(k) });
  };
  vtu.pointData.forEach((a) => add(a, false));
  vtu.cellData.forEach((a) => add(a, true));
  return out;
}

// The displacement of each node: the point vector u, or else the first
// point vector, or null.
function displacement(vtu) {
  const v = vtu.pointData.find((a) => a.name === 'u' && a.comps === 3) || vtu.pointData.find((a) => a.comps === 3);
  return v ? v.data : null;
}

// The triangles of the cut of the tetrahedra of `mesh` by the plane
// n·x = d. Each vertex is [edge a, edge b, t]: the point a + t (b - a).
function cutPlane(mesh, n, d) {
  const p = mesh.vtu.points;
  const out = [];
  const outCell = [];
  const s = (i) => n[0] * p[3 * i] + n[1] * p[3 * i + 1] + n[2] * p[3 * i + 2] - d;
  const edge = (a, b, sa, sb) => [a, b, sa / (sa - sb)];
  mesh.cells.forEach((cell, c) => {
    const tets = TETS[cell.type];
    if (!tets) return;
    for (const t of tets) {
      const q = t.map((k) => cell.nodes[k]);
      const v = q.map(s);
      const below = [0, 1, 2, 3].filter((k) => v[k] < 0);
      const above = [0, 1, 2, 3].filter((k) => v[k] >= 0);
      if (below.length === 0 || above.length === 0) continue;
      const e = (i, j) => edge(q[i], q[j], v[i], v[j]);
      if (below.length === 1 || above.length === 1) {
        const one = below.length === 1 ? below[0] : above[0];
        const rest = below.length === 1 ? above : below;
        out.push(e(one, rest[0]), e(one, rest[1]), e(one, rest[2]));
        outCell.push(c);
      } else {
        const [a, b] = below;
        const [c1, d1] = above;
        out.push(e(a, c1), e(a, d1), e(b, d1), e(a, c1), e(b, d1), e(b, c1));
        outCell.push(c, c);
      }
    }
  });
  return { verts: out, cells: outCell };
}

// --- the source --------------------------------------------------------------------

// The tokens of a Skein expression: numbers, names, and operators.
function tokens(s) {
  const out = [];
  const re = /\s*(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|[A-Za-z_]\w*|==|!=|<=|>=|&&|\|\||[-+*/^(),.|&<>:\[\]=])/y;
  let m;
  while ((m = re.exec(s))) out.push(m[1]);
  return out;
}

// The factor of each unit to the unit system of the model: length to the
// length of the header, an angle to radians. Other units are 1: the
// regions need only lengths and angles.
function unitFactor(word, metre) {
  const len = { m: 1, cm: 0.01, mm: 0.001, um: 1e-6, nm: 1e-9, km: 1000 };
  if (word in len) return len[word] * metre;
  if (word === 'deg') return Math.PI / 180;
  return 1;
}

const UNIT_WORDS = new Set(['m', 'cm', 'mm', 'um', 'nm', 'km', 'deg', 'rad', 'K', 'degC', 'N', 'kN', 'MN', 'Pa', 'kPa', 'MPa',
  'GPa', 'W', 'J', 's', 'ms', 'kg', 'g', 'Hz', 'V', 'A', 'mol']);

// A parser of the expressions of the geometry, as a tree of nodes. The
// operators have the precedence of spec §3: `|| && | &`, then comparisons,
// then `+ -`, `* /`, unary `-` and `^`.
function parseExpr(toks, metre) {
  let i = 0;
  const peek = () => toks[i];
  const take = (t) => { if (toks[i] === t) { i++; return true; } return false; };
  const unit = () => {
    // a unit after a value: mm, K^-1, W/(m*K): only the length factor matters
    let f = 1;
    while (i < toks.length && UNIT_WORDS.has(toks[i])) {
      f *= unitFactor(toks[i++], metre);
      if (take('^')) { const neg = take('-'); const e = +toks[i++]; f = Math.pow(f, neg ? -e : e); }
      if (peek() === '/' && UNIT_WORDS.has(toks[i + 1])) { i++; f /= unitFactor(toks[i++], metre); }
    }
    return f;
  };
  function primary() {
    const t = toks[i++];
    if (t === undefined) throw new Error('end of the expression');
    if (/^[\d.]/.test(t)) return { k: 'num', v: +t * unit() };
    if (t === '(') {
      const items = [expr()];
      while (take(',')) items.push(expr());
      if (!take(')')) throw new Error('expected )');
      const f = unit();
      return items.length === 1 ? (f === 1 ? items[0] : { k: 'mul', a: items[0], b: { k: 'num', v: f } }) : { k: 'vec', items, f };
    }
    if (t === '-') return { k: 'neg', a: power() };
    if (/^[A-Za-z_]/.test(t)) {
      if (take('(')) {
        const args = [];
        if (!take(')')) {
          do {
            if (/^[A-Za-z_]/.test(peek()) && toks[i + 1] === ':') i += 2;   // a named argument
            args.push(expr());
          } while (take(','));
          if (!take(')')) throw new Error('expected )');
        }
        return { k: 'call', name: t, args };
      }
      let node = { k: 'name', name: t };
      while (peek() === '.' && /^[A-Za-z_]/.test(toks[i + 1])) { i++; node = { k: 'member', a: node, name: toks[i++] }; }
      return node;
    }
    throw new Error('unexpected ' + t);
  }
  function power() {
    const a = primary();
    if (take('^')) return { k: 'pow', a, b: power() };
    return a;
  }
  function product() {
    let a = power();
    while (peek() === '*' || peek() === '/') { const op = toks[i++]; a = { k: op === '*' ? 'mul' : 'div', a, b: power() }; }
    return a;
  }
  function sum() {
    let a = product();
    while (peek() === '+' || peek() === '-') { const op = toks[i++]; a = { k: op === '+' ? 'add' : 'sub', a, b: product() }; }
    return a;
  }
  function compare() {
    let a = sum();
    while (['==', '!=', '<', '<=', '>', '>='].includes(peek())) { const op = toks[i++]; a = { k: 'cmp', op, a, b: sum() }; }
    return a;
  }
  function expr() {
    let a = compare();
    while (['||', '&&', '|', '&'].includes(peek())) { const op = toks[i++]; a = { k: 'bool', op, a, b: compare() }; }
    return a;
  }
  const out = expr();
  return out;
}

// A shape: inside(p, tol) is true for a point in the shape grown by tol.
function boxShape(a, b) {
  const lo = [0, 1, 2].map((k) => Math.min(a[k], b[k]));
  const hi = [0, 1, 2].map((k) => Math.max(a[k], b[k]));
  return { inside: (p, t) => [0, 1, 2].every((k) => p[k] >= lo[k] - t && p[k] <= hi[k] + t) };
}

function cylinderShape(c, axis, r, h) {
  const l = Math.hypot(axis[0], axis[1], axis[2]) || 1;
  const n = axis.map((x) => x / l);
  return {
    inside: (p, t) => {
      const d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
      const s = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
      const q = Math.hypot(d[0] - s * n[0], d[1] - s * n[1], d[2] - s * n[2]);
      return s >= -t && s <= h + t && q <= r + t;
    },
  };
}

// p turned by -angle about the axis through the origin (Rodrigues).
function unturn(p, axis, angle) {
  const l = Math.hypot(axis[0], axis[1], axis[2]) || 1;
  const k = axis.map((x) => x / l);
  const c = Math.cos(-angle);
  const s = Math.sin(-angle);
  const dot = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
  const cr = [k[1] * p[2] - k[2] * p[1], k[2] * p[0] - k[0] * p[2], k[0] * p[1] - k[1] * p[0]];
  return [0, 1, 2].map((j) => p[j] * c + cr[j] * s + k[j] * dot * (1 - c));
}

// The value of node `n`: a number, a vector (array), a shape, a region of
// the boundary or a point. `env.point` is the point of the x, y and z of a
// condition.
function evaluate(n, env) {
  const ev = (a) => evaluate(a, env);
  const num = (v) => (Array.isArray(v) ? v[0] : v);
  switch (n.k) {
    case 'num': return n.v;
    case 'vec': return n.items.map((x) => num(ev(x)) * n.f);
    case 'neg': { const a = ev(n.a); return Array.isArray(a) ? a.map((x) => -x) : -a; }
    case 'add': case 'sub': {
      const a = ev(n.a); const b = ev(n.b);
      if (n.k === 'sub' && a && a.inside && b && b.inside) {
        return { inside: (p, t) => a.inside(p, t) && !b.inside(p, -t) };
      }
      const s = n.k === 'add' ? 1 : -1;
      if (Array.isArray(a)) return a.map((x, k) => x + s * b[k]);
      return a + s * b;
    }
    case 'mul': { const a = ev(n.a); const b = ev(n.b); if (Array.isArray(a)) return a.map((x) => x * b); if (Array.isArray(b)) return b.map((x) => x * a); return a * b; }
    case 'div': { const a = ev(n.a); const b = ev(n.b); return Array.isArray(a) ? a.map((x) => x / b) : a / b; }
    case 'pow': return Math.pow(ev(n.a), ev(n.b));
    case 'cmp': {
      const a = num(ev(n.a)); const b = num(ev(n.b)); const t = env.tol;
      return { '==': Math.abs(a - b) <= t, '!=': Math.abs(a - b) > t, '<': a < b + t, '<=': a <= b + t, '>': a > b - t, '>=': a >= b - t }[n.op];
    }
    case 'bool': {
      const a = ev(n.a); const b = ev(n.b);
      if (a && a.inside) {
        if (n.op === '|') return { inside: (p, t) => a.inside(p, t) || b.inside(p, t) };
        return { inside: (p, t) => a.inside(p, t) && b.inside(p, t) };
      }
      return n.op === '||' || n.op === '|' ? a || b : a && b;
    }
    case 'name': {
      if (env.point && ['x', 'y', 'z'].includes(n.name)) return env.point['xyz'.indexOf(n.name)];
      if (n.name === 'pi') return Math.PI;
      const v = env.lets.get(n.name);
      if (v === undefined) throw new Error('unknown name ' + n.name);
      return v(env);
    }
    case 'member': {
      const a = ev(n.a);
      const k = 'xyz'.indexOf(n.name);
      if (Array.isArray(a) && k >= 0) return a[k];
      throw new Error('unknown member ' + n.name);
    }
    case 'call': return call(n, env);
    default: throw new Error('unknown node');
  }
}

function call(n, env) {
  const a = n.args.map((x) => (x.k === 'cmp' || x.k === 'bool' ? x : evaluate(x, env)));
  switch (n.name) {
    case 'box': return boxShape(a[0], a[1]);
    case 'rect': return boxShape([a[0][0], a[0][1], -Infinity], [a[1][0], a[1][1], Infinity]);
    case 'cylinder': return cylinderShape(a[0], a[1], a[2], a[3]);
    case 'circle': return cylinderShape([a[0][0], a[0][1], -1e30], [0, 0, 1], a[1], 2e30);
    case 'sphere': return { inside: (p, t) => Math.hypot(p[0] - a[0][0], p[1] - a[0][1], p[2] - a[0][2]) <= a[1] + t };
    case 'move': return { inside: (p, t) => a[0].inside([p[0] - a[1][0], p[1] - a[1][1], p[2] - (a[1][2] || 0)], t) };
    case 'rotate': return { inside: (p, t) => a[0].inside(unturn(p, a[1], a[2]), t) };
    case 'point': return { point: a[0] };
    case 'abs': return Math.abs(a[0]);
    case 'min': return Math.min(...a);
    case 'max': return Math.max(...a);
    case 'sqrt': return Math.sqrt(a[0]);
    case 'boundary': return { boundary: a[0], on: n.args[1] || null, env };
    default: throw new Error('the region function ' + n.name + ' is not in the viewer');
  }
}

// The items of the model source: each line that starts at column 0, with
// its continuation lines, and each statement inside a block.
function sourceItems(source) {
  const lines = source.split('\n');
  const items = [];
  for (let i = 0; i < lines.length; i++) {
    const raw = lines[i];
    const t = raw.trim();
    if (!t || t.startsWith('//') || t === '}') continue;
    const word = /^(\w+)/.exec(t);
    if (!word) continue;
    const indent = raw.length - raw.trimStart().length;
    let text = t;
    let end = i;
    // continuation lines are more indented, or follow a line that ends in an operator or a comma
    while (end + 1 < lines.length) {
      const next = lines[end + 1];
      const ni = next.length - next.trimStart().length;
      const open = /[|&,(+\-*/=]$/.test(lines[end].trim()) || (lines[end].match(/\(/g) || []).length > (lines[end].match(/\)/g) || []).length;
      if (!next.trim() || (!open && !(ni > indent && !/[{}]\s*$/.test(lines[end].trim()) && word[1] !== 'table')) || next.trim().startsWith('//')) break;
      end++;
      text += ' ' + next.trim();
    }
    items.push({ line: i + 1, last: end + 1, word: word[1], text, indent });
    i = end;
  }
  return items;
}

// The doc lines (///) above line `line`.
function docAbove(lines, line) {
  const out = [];
  for (let i = line - 2; i >= 0 && lines[i].trim().startsWith('///'); i--) out.unshift(lines[i].trim().slice(3).trim());
  return out.join(' ');
}

// The model of the source: the lets, params, regions, steps, outputs,
// expects and reports, with their lines.
function readModel(source) {
  const lines = source.split('\n');
  const items = sourceItems(source);
  const header = items.find((it) => it.word === 'model');
  const metre = header && /mm_/.test(header.text) ? 1000 : 1;
  const model = { lines, items, metre, lets: new Map(), letNodes: new Map(), params: [], regions: new Map(), steps: [], outputs: [],
    expects: [], reports: [], materials: [], meshes: [], physics: [], methods: [], bcs: [], solves: [], studies: [] };
  let step = null;
  let study = null;
  for (const it of items) {
    const t = it.text;
    let m;
    if (it.indent === 0) { step = null; study = null; }
    if ((m = /^(let|param)\s+(\w+)\s*:\s*([\w<>/*^() ]+?)\s*=\s*(.*?)(?:\s+in\s+(.*))?$/.exec(t))) {
      const [, kind, name, type, value, range] = m;
      let node = null;
      try { node = parseExpr(tokens(value), metre); } catch (e) { node = null; }
      model.letNodes.set(name, node);
      if (kind === 'param') model.params.push({ name, type, value, range, line: it.line, doc: docAbove(lines, it.line) });
      if (['Solid', 'Surface', 'Point', 'Curve'].includes(type.trim())) model.regions.set(name, { name, type: type.trim(), line: it.line, node });
      continue;
    }
    if ((m = /^step\s+(\w+)/.exec(t))) { step = m[1]; model.steps.push({ name: m[1], line: it.line }); continue; }
    if ((m = /^study\s+(\w+)/.exec(t))) { study = m[1]; model.studies.push({ name: m[1], line: it.line }); continue; }
    if ((m = /^bc\s+(\w+)\s*:\s*(.*)$/.exec(t))) { model.bcs.push({ region: m[1], text: m[2], step, line: it.line }); continue; }
    if ((m = /^output\s+(.*)$/.exec(t))) {
      for (const f of splitTop(m[1])) model.outputs.push({ expr: f, step, line: it.line });
      continue;
    }
    if ((m = /^expect\s+(\w+)\s*:\s*(.*)$/.exec(t))) { model.expects.push({ name: m[1], text: m[2], step, study, line: it.line }); continue; }
    if ((m = /^report\s+(.*)$/.exec(t))) { model.reports.push({ expr: m[1], step, line: it.line }); continue; }
    if ((m = /^solve\s*\{(.*)\}/.exec(t))) { model.solves.push({ text: m[1].trim(), step, line: it.line }); continue; }
    if ((m = /^mesh\s+(\w+)/.exec(t))) { model.meshes.push({ name: m[1], line: it.line }); continue; }
    if ((m = /^physics\s+(\w+)/.exec(t))) { model.physics.push({ name: m[1], line: it.line }); continue; }
    if ((m = /^method\s+(\w+)/.exec(t))) { model.methods.push({ name: m[1], line: it.line }); continue; }
    if ((m = /^material\s+(.*)$/.exec(t)) && it.indent > 0) {
      for (const part of splitTop(m[1])) {
        const pm = /^(\w+)(?:\s+on\s+(.*))?$/.exec(part.trim());
        if (pm) model.materials.push({ name: pm[1], region: pm[2] || null, line: it.line });
      }
    }
  }
  return model;
}

// The parts of `s` that commas at the top level separate.
function splitTop(s) {
  const out = [];
  let depth = 0;
  let from = 0;
  for (let i = 0; i < s.length; i++) {
    if (s[i] === '(' || s[i] === '[') depth++;
    if (s[i] === ')' || s[i] === ']') depth--;
    if (s[i] === ',' && depth === 0) { out.push(s.slice(from, i).trim()); from = i + 1; }
  }
  out.push(s.slice(from).trim());
  return out.filter((x) => x.length);
}

// The environment of the lets for the param values `values` (name → text).
function letEnv(model, values, tol) {
  const env = { lets: new Map(), tol, point: null };
  const cache = new Map();
  for (const [name, node] of model.letNodes) {
    let n = node;
    if (values && values[name] !== undefined) {
      try { n = parseExpr(tokens(values[name]), model.metre); } catch (e) { n = node; }
    }
    env.lets.set(name, () => {
      if (!n) throw new Error('cannot read ' + name);
      if (!cache.has(name)) cache.set(name, evaluate(n, { lets: env.lets, tol: env.tol, point: null }));
      return cache.get(name);
    });
  }
  return env;
}

// The node of `mesh` nearest to the point q.
function nearestNode(mesh, q) {
  const p = mesh.vtu.points;
  let best = -1;
  let dist = Infinity;
  for (let i = 0; i < p.length / 3; i++) {
    const d = Math.hypot(p[3 * i] - q[0], p[3 * i + 1] - q[1], p[3 * i + 2] - (q[2] || 0));
    if (d < dist) { dist = d; best = i; }
  }
  return best;
}

// The cells, the boundary triangles or the node of region `r` on `mesh`,
// for the param values `values`. It gives { tris, cells, nodes } as sets,
// or throws when the viewer cannot read the region.
function regionOf(model, mesh, r, values) {
  const tol = 1e-5 * mesh.size;
  const env = letEnv(model, values, tol);
  const v = env.lets.get(r.name)(env);
  const p = mesh.vtu.points;
  const centre = (nodes) => {
    const c = [0, 0, 0];
    for (const i of nodes) for (let k = 0; k < 3; k++) c[k] += p[3 * i + k] / nodes.length;
    return c;
  };
  const out = { tris: new Set(), cells: new Set(), nodes: new Set() };
  if (v && v.point) {
    out.nodes.add(nearestNode(mesh, v.point));
    return out;
  }
  if (v && v.inside) {
    mesh.cells.forEach((cell, c) => { if (v.inside(centre(cell.nodes), tol)) out.cells.add(c); });
    for (let t = 0; t < mesh.triCell.length; t++) if (out.cells.has(mesh.triCell[t])) out.tris.add(t);
    return out;
  }
  if (v && v.boundary) {
    const body = v.boundary;
    const on = v.on;
    let shape = null;
    if (on && on.k !== 'cmp' && on.k !== 'bool') shape = evaluate(on, env);
    else if (on && on.k === 'bool') {
      try { const s = evaluate(on, env); if (s && s.inside) shape = s; } catch (e) { shape = null; }
    }
    for (let t = 0; t < mesh.triCell.length; t++) {
      const c = centre([mesh.tris[3 * t], mesh.tris[3 * t + 1], mesh.tris[3 * t + 2]]);
      if (body.inside && !body.inside(c, tol)) continue;
      let good = true;
      if (shape) good = shape.inside(c, tol) && !shape.inside(c, -tol);
      else if (on) good = !!evaluate(on, { lets: env.lets, tol, point: c });
      if (good) out.tris.add(t);
    }
    return out;
  }
  throw new Error('the viewer cannot show this region');
}

// --- WebGL ---------------------------------------------------------------------------

const VERTEX = `#version 300 es
in vec3 aPos; in vec3 aDisp; in float aVal; in float aMark;
uniform mat4 uMvp; uniform float uScale; uniform float uPointSize;
out vec3 vPos; out vec3 vWorld; out float vVal; out float vMark; flat out int vId;
void main() {
  vec3 p = aPos + uScale * aDisp;
  vPos = aPos; vWorld = p; vVal = aVal; vMark = aMark; vId = gl_VertexID / 3;
  gl_Position = uMvp * vec4(p, 1.0);
  gl_PointSize = uPointSize;
}`;

const FRAGMENT = `#version 300 es
precision highp float;
in vec3 vPos; in vec3 vWorld; in float vVal; in float vMark; flat in int vId;
uniform float uMin; uniform float uMax; uniform int uBands; uniform vec4 uClip; uniform int uClipOn;
uniform int uMode; uniform vec3 uColor; uniform vec3 uEye; uniform int uIdBase;
out vec4 outColor;
vec3 turbo(float t) {
  const vec4 r4 = vec4(0.13572138, 4.61539260, -42.66032258, 132.13108234);
  const vec4 g4 = vec4(0.09140261, 2.19418839, 4.84296658, -14.18503333);
  const vec4 b4 = vec4(0.10667330, 12.64194608, -60.58204836, 110.36276771);
  const vec2 r2 = vec2(-152.94239396, 59.28637943);
  const vec2 g2 = vec2(4.27729857, 2.82956604);
  const vec2 b2 = vec2(-89.90310912, 27.34824973);
  vec4 v4 = vec4(1.0, t, t * t, t * t * t);
  vec2 v2 = v4.zw * v4.z;
  return vec3(dot(v4, r4) + dot(v2, r2), dot(v4, g4) + dot(v2, g2), dot(v4, b4) + dot(v2, b2));
}
void main() {
  if (uClipOn == 1 && dot(uClip.xyz, vPos) > uClip.w) discard;
  if (uMode == 2) {
    int id = vId + uIdBase;
    outColor = vec4(float(id & 255) / 255.0, float((id >> 8) & 255) / 255.0, float((id >> 16) & 255) / 255.0, 1.0);
    return;
  }
  vec3 base = uColor;
  if (uMode == 0) {
    float t = clamp((vVal - uMin) / max(uMax - uMin, 1e-30), 0.0, 1.0);
    if (uBands > 0) t = (floor(t * float(uBands) * 0.9999) + 0.5) / float(uBands);
    base = turbo(t);
  }
  if (vMark > 0.5) base = mix(base, vec3(1.0, 0.1, 0.75), 0.65);
  if (uMode == 3) { outColor = vec4(base, 1.0); return; }
  vec3 n = normalize(cross(dFdx(vWorld), dFdy(vWorld)));
  float light = 0.4 + 0.6 * abs(dot(n, normalize(uEye)));
  outColor = vec4(base * light, 1.0);
}`;

function compile(gl, type, text) {
  const s = gl.createShader(type);
  gl.shaderSource(s, text);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}

function mat4Mul(a, b) {
  const o = new Float32Array(16);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    let s = 0;
    for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k];
    o[i * 4 + j] = s;
  }
  return o;
}

function perspective(fov, aspect, near, far) {
  const f = 1 / Math.tan(fov / 2);
  return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1, 0, 0, 2 * far * near / (near - far), 0]);
}

function lookAt(eye, at, up) {
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const norm = (a) => { const l = Math.hypot(...a) || 1; return a.map((x) => x / l); };
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const z = norm(sub(eye, at));
  const x = norm(cross(up, z));
  const y = cross(z, x);
  return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot(x, eye), -dot(y, eye), -dot(z, eye), 1]);
}

// A WebGL 2 view of one mesh: the surface, the cut, the mesh lines and the
// marked points.
class View {
  constructor(canvas) {
    this.canvas = canvas;
    const gl = canvas.getContext('webgl2', { antialias: true, preserveDrawingBuffer: true });
    if (!gl) throw new Error('this browser has no WebGL 2');
    this.gl = gl;
    const prog = gl.createProgram();
    gl.attachShader(prog, compile(gl, gl.VERTEX_SHADER, VERTEX));
    gl.attachShader(prog, compile(gl, gl.FRAGMENT_SHADER, FRAGMENT));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    this.prog = prog;
    this.loc = {};
    for (const a of ['aPos', 'aDisp', 'aVal', 'aMark']) this.loc[a] = gl.getAttribLocation(prog, a);
    for (const u of ['uMvp', 'uScale', 'uPointSize', 'uMin', 'uMax', 'uBands', 'uClip', 'uClipOn', 'uMode', 'uColor', 'uEye', 'uIdBase']) {
      this.loc[u] = gl.getUniformLocation(prog, u);
    }
    this.parts = {};
    this.yaw = -0.9; this.pitch = 0.5; this.dist = 1; this.target = [0, 0, 0];
    this.scale = 0; this.range = [0, 1]; this.bands = 0; this.clip = null; this.lines = false;
    this.pick = gl.createFramebuffer();
  }

  // Sets the vertex arrays of part `name`: positions, displacements, values
  // and marks, one entry for each vertex.
  setPart(name, pos, disp, val, mark, mode) {
    const gl = this.gl;
    const old = this.parts[name];
    if (old) for (const b of old.buffers) gl.deleteBuffer(b);
    const vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    const buffers = [];
    const put = (loc, data, size) => {
      const b = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, b);
      gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(loc);
      gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0);
      buffers.push(b);
    };
    put(this.loc.aPos, pos, 3);
    put(this.loc.aDisp, disp, 3);
    put(this.loc.aVal, val, 1);
    put(this.loc.aMark, mark, 1);
    gl.bindVertexArray(null);
    this.parts[name] = { vao, buffers, count: pos.length / 3, mode };
  }

  fit(lo, hi) {
    this.target = [0, 1, 2].map((k) => (lo[k] + hi[k]) / 2);
    this.dist = 2.2 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
    this.size = this.dist;
  }

  eye() {
    const c = Math.cos(this.pitch);
    return [this.target[0] + this.dist * c * Math.cos(this.yaw), this.target[1] + this.dist * c * Math.sin(this.yaw),
      this.target[2] + this.dist * Math.sin(this.pitch)];
  }

  matrix() {
    const w = this.canvas.width; const h = this.canvas.height;
    const proj = perspective(0.6, w / Math.max(h, 1), this.dist * 0.01, this.dist * 10);
    return mat4Mul(proj, lookAt(this.eye(), this.target, [0, 0, 1]));
  }

  draw(mode) {
    const gl = this.gl;
    const dpr = window.devicePixelRatio || 1;
    const w = Math.round(this.canvas.clientWidth * dpr); const h = Math.round(this.canvas.clientHeight * dpr);
    if (this.canvas.width !== w || this.canvas.height !== h) { this.canvas.width = w; this.canvas.height = h; }
    gl.viewport(0, 0, w, h);
    const bg = getComputedStyle(this.canvas).getPropertyValue('--stage').trim() || '#f4f4f2';
    const c = parseInt(bg.slice(1), 16);
    if (mode === 'pick') gl.clearColor(0, 0, 0, 1); else gl.clearColor((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255, 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    gl.useProgram(this.prog);
    const L = this.loc;
    gl.uniformMatrix4fv(L.uMvp, false, this.matrix());
    gl.uniform1f(L.uScale, this.scale);
    gl.uniform1f(L.uMin, this.range[0]);
    gl.uniform1f(L.uMax, this.range[1]);
    gl.uniform1i(L.uBands, this.bands);
    gl.uniform3fv(L.uEye, this.eye().map((x, k) => x - this.target[k]));
    gl.uniform1i(L.uClipOn, this.clip ? 1 : 0);
    if (this.clip) {
      // The cut removes the half on the side of the eye, so the cut face shows.
      const [nx, ny, nz, d] = this.clip; const e = this.eye();
      gl.uniform4fv(L.uClip, nx * e[0] + ny * e[1] + nz * e[2] > d ? this.clip : [-nx, -ny, -nz, -d]);
    }
    gl.uniform1f(L.uPointSize, 12 * dpr);
    let base = 1;
    for (const name of ['surface', 'cut', 'beams', 'lines', 'points']) {
      const p = this.parts[name];
      if (!p || !p.count || (name === 'lines' && (!this.lines || mode === 'pick'))) continue;
      if ((name === 'points' || name === 'beams') && mode === 'pick') continue;
      gl.uniform1i(L.uClipOn, this.clip && name !== 'cut' && name !== 'points' ? 1 : 0);
      gl.uniform1i(L.uMode, mode === 'pick' ? 2 : p.mode);
      gl.uniform1i(L.uIdBase, base);
      gl.uniform3fv(L.uColor, name === 'lines' ? [0.15, 0.15, 0.15] : name === 'points' ? [1, 0.1, 0.75] : [0.8, 0.8, 0.8]);
      if (name === 'surface' || name === 'cut') { gl.enable(gl.POLYGON_OFFSET_FILL); gl.polygonOffset(1, 1); } else gl.disable(gl.POLYGON_OFFSET_FILL);
      gl.bindVertexArray(p.vao);
      gl.drawArrays(name === 'lines' || name === 'beams' ? gl.LINES : name === 'points' ? gl.POINTS : gl.TRIANGLES, 0, p.count);
      base += p.count / 3;
    }
    gl.bindVertexArray(null);
  }

  // The part and the triangle at the pixel (x, y) of the canvas, or null.
  pickAt(x, y) {
    const gl = this.gl;
    const dpr = window.devicePixelRatio || 1;
    const w = this.canvas.width; const h = this.canvas.height;
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, w, h, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
    const depth = gl.createRenderbuffer();
    gl.bindRenderbuffer(gl.RENDERBUFFER, depth);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT24, w, h);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pick);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depth);
    this.draw('pick');
    const px = new Uint8Array(4);
    gl.readPixels(Math.round(x * dpr), Math.round(h - y * dpr), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.deleteTexture(tex);
    gl.deleteRenderbuffer(depth);
    this.draw();
    const id = px[0] + (px[1] << 8) + (px[2] << 16) - 1;
    if (id < 0) return null;
    const surface = this.parts.surface ? this.parts.surface.count / 3 : 0;
    return id < surface ? { part: 'surface', tri: id } : { part: 'cut', tri: id - surface };
  }

  // The ray of the pixel (x, y): its origin and direction.
  ray(x, y) {
    const inv = invert(this.matrix());
    const nx = 2 * x / this.canvas.clientWidth - 1; const ny = 1 - 2 * y / this.canvas.clientHeight;
    const at = (z) => {
      const v = [0, 1, 2, 3].map((r) => inv[r] * nx + inv[4 + r] * ny + inv[8 + r] * z + inv[12 + r]);
      return [v[0] / v[3], v[1] / v[3], v[2] / v[3]];
    };
    const a = at(-1); const b = at(1);
    return [a, [b[0] - a[0], b[1] - a[1], b[2] - a[2]]];
  }
}

function invert(m) {
  const inv = new Float32Array(16);
  const a = Array.from(m);
  inv[0] = a[5] * a[10] * a[15] - a[5] * a[11] * a[14] - a[9] * a[6] * a[15] + a[9] * a[7] * a[14] + a[13] * a[6] * a[11] - a[13] * a[7] * a[10];
  inv[4] = -a[4] * a[10] * a[15] + a[4] * a[11] * a[14] + a[8] * a[6] * a[15] - a[8] * a[7] * a[14] - a[12] * a[6] * a[11] + a[12] * a[7] * a[10];
  inv[8] = a[4] * a[9] * a[15] - a[4] * a[11] * a[13] - a[8] * a[5] * a[15] + a[8] * a[7] * a[13] + a[12] * a[5] * a[11] - a[12] * a[7] * a[9];
  inv[12] = -a[4] * a[9] * a[14] + a[4] * a[10] * a[13] + a[8] * a[5] * a[14] - a[8] * a[6] * a[13] - a[12] * a[5] * a[10] + a[12] * a[6] * a[9];
  inv[1] = -a[1] * a[10] * a[15] + a[1] * a[11] * a[14] + a[9] * a[2] * a[15] - a[9] * a[3] * a[14] - a[13] * a[2] * a[11] + a[13] * a[3] * a[10];
  inv[5] = a[0] * a[10] * a[15] - a[0] * a[11] * a[14] - a[8] * a[2] * a[15] + a[8] * a[3] * a[14] + a[12] * a[2] * a[11] - a[12] * a[3] * a[10];
  inv[9] = -a[0] * a[9] * a[15] + a[0] * a[11] * a[13] + a[8] * a[1] * a[15] - a[8] * a[3] * a[13] - a[12] * a[1] * a[11] + a[12] * a[3] * a[9];
  inv[13] = a[0] * a[9] * a[14] - a[0] * a[10] * a[13] - a[8] * a[1] * a[14] + a[8] * a[2] * a[13] + a[12] * a[1] * a[10] - a[12] * a[2] * a[9];
  inv[2] = a[1] * a[6] * a[15] - a[1] * a[7] * a[14] - a[5] * a[2] * a[15] + a[5] * a[3] * a[14] + a[13] * a[2] * a[7] - a[13] * a[3] * a[6];
  inv[6] = -a[0] * a[6] * a[15] + a[0] * a[7] * a[14] + a[4] * a[2] * a[15] - a[4] * a[3] * a[14] - a[12] * a[2] * a[7] + a[12] * a[3] * a[6];
  inv[10] = a[0] * a[5] * a[15] - a[0] * a[7] * a[13] - a[4] * a[1] * a[15] + a[4] * a[3] * a[13] + a[12] * a[1] * a[7] - a[12] * a[3] * a[5];
  inv[14] = -a[0] * a[5] * a[14] + a[0] * a[6] * a[13] + a[4] * a[1] * a[14] - a[4] * a[2] * a[13] - a[12] * a[1] * a[6] + a[12] * a[2] * a[5];
  inv[3] = -a[1] * a[6] * a[11] + a[1] * a[7] * a[10] + a[5] * a[2] * a[11] - a[5] * a[3] * a[10] - a[9] * a[2] * a[7] + a[9] * a[3] * a[6];
  inv[7] = a[0] * a[6] * a[11] - a[0] * a[7] * a[10] - a[4] * a[2] * a[11] + a[4] * a[3] * a[10] + a[8] * a[2] * a[7] - a[8] * a[3] * a[6];
  inv[11] = -a[0] * a[5] * a[11] + a[0] * a[7] * a[9] + a[4] * a[1] * a[11] - a[4] * a[3] * a[9] - a[8] * a[1] * a[7] + a[8] * a[3] * a[5];
  inv[15] = a[0] * a[5] * a[10] - a[0] * a[6] * a[9] - a[4] * a[1] * a[10] + a[4] * a[2] * a[9] + a[8] * a[1] * a[6] - a[8] * a[2] * a[5];
  const d = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12] || 1;
  return inv.map((x) => x / d);
}

// --- the page ----------------------------------------------------------------------

function el(tag, props, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props || {})) {
    if (k === 'class') e.className = v;
    else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else if (v !== undefined && v !== null && v !== false) e.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) if (kid !== null && kid !== undefined) e.append(kid instanceof Node ? kid : String(kid));
  return e;
}

function fmt(x) {
  if (!isFinite(x)) return String(x);
  const a = Math.abs(x);
  return a !== 0 && (a < 1e-3 || a >= 1e5) ? x.toExponential(4) : (+x.toPrecision(5)).toString();
}

// The number of a value text of the model, such as "150 um", or NaN.
function valueOf(model, text) {
  try { const v = evaluate(parseExpr(tokens(text), model.metre), { lets: new Map(), tol: 0, point: null }); return Array.isArray(v) ? NaN : v; } catch (e) { return NaN; }
}

const D = JSON.parse(document.getElementById('skein-data').textContent);
const M = readModel(D.source || '');
const S = {
  caseIndex: 0, values: {}, step: 0, field: null, scale: 1, bands: 0, cut: 'off', cutAt: 0.5, mark: null, probe: null,
  mesh: null, foreign: null, lines: false,
};
const meshes = new Map();
let view = null;

const checkToon = parseToon(D.check);
const verifyToon = parseToon(D.verify);
const tableToon = parseToon(D.table);
const caseToons = (D.cases || []).map((c) => parseToon(c.run));
const studyTable = tableToon.tables.cases || null;
const studyParams = studyTable ? studyTable.cols.slice(1, studyTable.cols.indexOf('status')) : [];
// The solve of the page, or false.
const SOLVE = typeof skeinSolve === 'function' && MODULE.length > 0;
const solveParams = SOLVE ? PARAMS.filter((p) => M.params.some((x) => x.name === p.name)) : [];

// The param values of case i, name → text: the values of a solve in the
// page, else the study values and the default of each other param of PARAMS.
function caseValues(i) {
  if (D.cases[i] && D.cases[i].values) return { ...D.cases[i].values };
  const out = {};
  for (const p of solveParams) out[p.name] = M.params.find((x) => x.name === p.name).value;
  if (!studyTable || !studyTable.rows[i]) return out;
  studyParams.forEach((p, k) => { out[p] = studyTable.rows[i][k + 1]; });
  return out;
}

function sameValue(a, b) {
  const x = valueOf(M, a), y = valueOf(M, b);
  return Math.abs(x - y) <= 1e-9 * Math.max(Math.abs(x), Math.abs(y), 1e-300);
}

// The case of the param values `values` (name → value text): a stored case
// or a solve in the page, or null.
async function caseFor(values) {
  const i = D.cases.findIndex((c, k) => {
    const v = caseValues(k);
    return Object.keys(values).every((p) => !(p in v) || sameValue(v[p], values[p]));
  });
  return i < 0 ? null : D.cases[i];
}

// The value of the text `t` in the unit of the param p of PARAMS.
function inUnit(p, t) { return valueOf(M, t) / (p.unit ? valueOf(M, '1 ' + p.unit) : 1); }

let solves = 0;
// Whether a solve in the page puts the implicit solves on the GPU too.
let implicitGpu = false;

// Meshes and solves the param values `values` in the page. It gives the
// new case, or null if the solve fails or a newer solve stops it.
async function solveCase(values) {
  const note = document.getElementById('note');
  note.textContent = 'meshing and solving ' + solveParams.map((p) => p.name + ' = ' + values[p.name]).join(', ') + ' in the page';
  document.body.dataset.solve = 'running';
  const nums = {};
  for (const p of solveParams) nums[p.name] = inUnit(p, values[p.name]);
  const r = await skeinSolve(nums, { implicitGpu });
  if (r.code === -2) return null;
  document.body.dataset.solve = String(r.code);
  const results = Object.entries(r.files).map(([path, text]) => ({ path, text }));
  if (!results.length) {
    note.textContent = failure(r);
    return null;
  }
  note.textContent = '';
  const c = { id: 'page ' + ++solves, values, run: r.out, err: r.err, results, seconds: r.seconds, code: r.code, failure: r.code === 0 ? null : failure(r) };
  D.cases.push(c);
  caseToons.push(parseToon(r.out));
  return c;
}

// The parsed VTU of blob k.
async function blobVtu(k) {
  return parseVtu(await inflate(fromBase64(document.getElementById('vtu-' + k).textContent), 'gzip'), true);
}

// The mesh of result r: its mesh blob with its fields blob, decoded and
// built once. The cache keeps the last 12 meshes. Reason: a study of 100
// cases of 50000 nodes needs about 1 GB for all its meshes. A result of a
// solve in the page has the VTU text, and keeps its mesh.
const geometries = new Map();
function meshOf(r) {
  if ('text' in r) {
    if (!r.built) r.built = parseVtu(new TextEncoder().encode(r.text)).then(buildMesh);
    return r.built;
  }
  if (!meshes.has(r.blob)) {
    if (!geometries.has(r.mesh)) geometries.set(r.mesh, blobVtu(r.mesh));
    meshes.set(r.blob, (async () => {
      const [geo, data] = await Promise.all([geometries.get(r.mesh), blobVtu(r.blob)]);
      return buildMesh({ ...geo, pointData: data.pointData, cellData: data.cellData });
    })());
    if (meshes.size > 12) meshes.delete(meshes.keys().next().value);
  }
  return meshes.get(r.blob);
}

function stepName(path) { return path.replace(/^.*\//, '').replace(/\.vtu$/, ''); }

function currentCase() { return D.cases[S.caseIndex]; }

// --- the source pane -----------------------------------------------------------------

function goLine(line) {
  const e = document.getElementById('L' + line);
  if (!e) return;
  e.scrollIntoView({ block: 'center', behavior: 'smooth' });
  e.classList.remove('flash');
  void e.offsetWidth;
  e.classList.add('flash');
}

function lineLink(line, text) {
  return el('a', { href: '#L' + line, class: 'line-link', 'data-line': line, onclick: (ev) => { ev.preventDefault(); goLine(line); } }, text || 'line ' + line);
}

// The field of the results that output expression `expr` names, or null.
function outputField(mesh, expr) {
  if (!mesh) return null;
  const f = mesh.fields;
  return f.find((x) => x.name === expr) || f.find((x) => x.name === '|' + expr + '|') ||
    f.find((x) => x.name === expr.replace(/\s/g, '')) || null;
}

// The field of an expect or report expression, such as `u.z at p`.
function exprField(mesh, text) {
  const m = /([A-Za-z_]\w*(?:\.\w+)?(?:\([^)]*\))?)\s+at\b/.exec(text);
  return m ? outputField(mesh, m[1]) : null;
}

function renderSource() {
  const pre = document.getElementById('source');
  pre.textContent = '';
  const tokensOf = new Map();
  for (const t of D.tokens || []) tokensOf.set(t[0], (tokensOf.get(t[0]) || []).concat([t]));
  const fieldLines = new Map();
  for (const o of M.outputs) fieldLines.set(o.line, (fieldLines.get(o.line) || []).concat([o.expr]));
  M.lines.forEach((text, i) => {
    const line = i + 1;
    const row = el('div', { id: 'L' + line, class: 'src-line' }, el('span', { class: 'ln' }, String(line)));
    const code = el('span', { class: 'code' });
    const ci = text.indexOf('//');
    const body = ci >= 0 ? text.slice(0, ci) : text;
    // The words and links of a piece of text.
    const words = (piece, into) => {
      const re = /[A-Za-z_]\w*/g;
      let at = 0;
      let m;
      while ((m = re.exec(piece))) {
        into.append(piece.slice(at, m.index));
        const w = m[0];
        if (M.regions.has(w)) {
          into.append(el('a', { href: '#', class: 'region', 'data-region': w, onclick: (ev) => { ev.preventDefault(); markRegion(w, line); } }, w));
        } else if (M.params.some((p) => p.name === w) && studyParams.includes(w)) {
          into.append(el('a', { href: '#', class: 'param', 'data-param': w, onclick: (ev) => { ev.preventDefault(); focusSlider(w); } }, w));
        } else {
          into.append(w);
        }
        at = m.index + w.length;
      }
      into.append(piece.slice(at));
    };
    // The token classes of the line, from the lexer of the compiler (token_classes in report.goose).
    let at = 0;
    for (const [, col, len, cls] of tokensOf.get(line) || []) {
      if (col < at || col + len > body.length) continue;
      words(body.slice(at, col), code);
      const tk = el('span', { class: 'tk-' + cls });
      words(body.slice(col, col + len), tk);
      code.append(tk);
      at = col + len;
    }
    words(body.slice(at), code);
    if (ci >= 0) code.append(el('span', { class: text.slice(ci).startsWith('///') ? 'doc' : 'comment' }, text.slice(ci)));
    row.append(code);
    const actions = sourceActions(line);
    if (actions.length) row.append(el('span', { class: 'actions' }, actions));
    pre.append(row);
  });
}

// The buttons after line `line`: show the field of an output, the point of
// an expect, the row of a report.
function sourceActions(line) {
  const out = [];
  for (const o of M.outputs.filter((x) => x.line === line)) {
    out.push(el('button', { class: 'button chip', 'data-output': o.expr, onclick: () => showOutput(o.expr) }, 'show ' + o.expr));
  }
  for (const e of M.expects.filter((x) => x.line === line)) {
    out.push(el('button', { class: 'button chip', 'data-expect': e.name, onclick: () => showCheck('expect', e.name, e.text) }, 'check ' + e.name));
  }
  for (const r of M.reports.filter((x) => x.line === line)) {
    out.push(el('button', { class: 'button chip', 'data-report': r.expr, onclick: () => showCheck('report', r.expr, r.expr) }, 'report'));
  }
  for (const st of M.steps.filter((x) => x.line === line)) {
    out.push(el('button', { class: 'button chip', 'data-step': st.name, onclick: () => selectStep(st.name) }, 'show step'));
  }
  for (const m of M.materials.filter((x) => x.line === line && !x.region)) {
    out.push(el('button', { class: 'button chip', onclick: () => markRegion(null, line) }, 'material ' + m.name));
  }
  return out;
}

// Shows the field of output `expr` and goes to its control.
function showOutput(expr) {
  const f = outputField(S.mesh, expr);
  const note = document.getElementById('note');
  if (!f) { note.textContent = expr + ' is not in the result files of this step'; return; }
  S.field = f.name;
  note.textContent = '';
  refresh();
}

// Marks the point of an expect or a report, shows its field, and marks
// its row in the tables.
function showCheck(kind, name, text) {
  const f = exprField(S.mesh, text);
  if (f) S.field = f.name;
  const at = /\bat\s+(\w+|\([^)]*\)\s*\w*)/.exec(text);
  S.mark = null;
  if (at && M.regions.has(at[1])) S.mark = { name: at[1] };
  else if (at) S.mark = { point: at[1] };
  refresh();
  document.querySelectorAll('tr.flash').forEach((r) => r.classList.remove('flash'));
  const rows = document.querySelectorAll(`tr[data-${kind}="${CSS.escape(name)}"]`);
  rows.forEach((r) => r.classList.add('flash'));
  if (rows[0]) rows[0].scrollIntoView({ block: 'nearest', behavior: 'smooth' });
}

function markRegion(name, line) {
  S.mark = name ? { name } : { all: true };
  document.querySelectorAll('.src-line.active').forEach((e) => e.classList.remove('active'));
  const e = document.getElementById('L' + line);
  if (e) e.classList.add('active');
  refresh();
}

// --- the controls and the sliders ------------------------------------------------------

function select(id, options, value, onchange) {
  const s = el('select', { id, onchange: (ev) => onchange(ev.target.value) });
  for (const o of options) s.append(el('option', { value: o, selected: o === value }, o));
  return s;
}

function renderControls() {
  const c = document.getElementById('controls');
  c.textContent = '';
  document.getElementById('note').textContent = '';
  const mesh = S.mesh;
  if (!mesh) return;
  const fields = mesh.fields.map((f) => f.name);
  const results = S.foreign ? [] : currentCase().results;
  c.append(...[
    el('label', {}, 'Field ', select('field', fields, S.field, (v) => { S.field = v; refresh(); })),
    results.length > 1 ? el('label', {}, 'Step ', select('step', results.map((r) => stepName(r.path)), stepName(results[S.step].path),
      (v) => selectStep(v))) : null,
    el('label', {}, 'Deform ', el('input', { id: 'scale', type: 'number', min: 0, step: 'any', value: fmt(S.scale), class: 'num',
      onchange: (ev) => { S.scale = Math.max(0, +ev.target.value || 0); draw(); } }),
    mesh.disp ? el('button', { class: 'button', onclick: () => { S.scale = autoScale(mesh); renderControls(); draw(); } }, 'auto') : null),
    el('label', {}, 'Cut ', select('cut', ['off', 'x', 'y', 'z'], S.cut, (v) => { S.cut = v; refresh(); }),
      el('input', { id: 'cut-at', type: 'range', min: 0, max: 1000, value: Math.round(S.cutAt * 1000),
        oninput: (ev) => { S.cutAt = ev.target.value / 1000; refresh(); } })),
    el('label', {}, el('input', { type: 'checkbox', id: 'bands', checked: S.bands > 0, onchange: (ev) => { S.bands = ev.target.checked ? 12 : 0; draw(); } }), ' bands'),
    el('label', {}, el('input', { type: 'checkbox', id: 'edges', checked: S.lines, onchange: (ev) => { S.lines = ev.target.checked; draw(); } }), ' mesh'),
    el('button', { class: 'button', onclick: () => { view.fit(mesh.lo, mesh.hi); draw(); } }, 'fit'),
    el('label', { class: 'button file' }, 'Open VTU', el('input', { type: 'file', accept: '.vtu', id: 'open-vtu', onchange: openFile }))].filter(Boolean));
}

// The deformation factor that moves the largest displacement by 1/10 of the
// size of the mesh.
function autoScale(mesh) {
  const u = mesh.disp;
  if (!u) return 0;
  let m = 0;
  for (let i = 0; i < u.length; i++) m = Math.max(m, Math.abs(u[i]));
  return m > 0 ? +(0.1 * mesh.size / m).toPrecision(2) : 0;
}

function selectStep(name) {
  const k = currentCase().results.findIndex((r) => stepName(r.path) === name);
  if (k >= 0) { S.step = k; load(); }
}

// A slider for each param of the study and of PARAMS. S.values holds the
// values of the sliders.
function renderSliders() {
  const box = document.getElementById('sliders');
  box.textContent = '';
  const names = studyParams.concat(solveParams.map((p) => p.name).filter((n) => !studyParams.includes(n)));
  if (!names.length) return;
  const study = M.studies.find((x) => x.name === D.study);
  box.append(studyTable ? el('h2', {}, 'Study ', D.study, ' ', study ? lineLink(study.line) : null) : el('h2', {}, 'Params'));
  if (SOLVE) {
    box.append(el('p', { class: 'note' }, 'A slider at the values of a stored case shows that case. A slider that stops at other values meshes and solves the model in the page.'));
    box.append(el('p', { class: 'note' }, HAS_GPU ?
      'WebGPU is present: the explicit steps run on the GPU, in f32. The implicit solves run on the CPU, unless you turn on the GPU for them.' :
      'This browser has no WebGPU: all steps run on the CPU.'));
    box.append(el('label', {}, el('input', { type: 'checkbox', id: 'implicit-gpu', checked: implicitGpu, disabled: !HAS_GPU,
      onchange: (ev) => { implicitGpu = ev.target.checked; } }), ' Implicit solves on the GPU'));
  }
  else if (typeof tooLarge === 'function') box.append(el('p', { class: 'note' }, tooLarge() + ' The sliders move only between the stored cases.'));
  for (const p of names) {
    const seen = [];
    for (let i = 0; i < D.cases.length; i++) {
      const t = caseValues(i)[p];
      if (t !== undefined && !seen.some((x) => sameValue(x, t))) seen.push(t);
    }
    seen.sort((a, b) => valueOf(M, a) - valueOf(M, b));
    const param = M.params.find((x) => x.name === p);
    const control = solveParams.find((x) => x.name === p);
    const out = el('output', { id: 'value-' + p }, S.values[p]);
    const set = (t, solve) => { S.values[p] = t; out.textContent = t; showValues(solve); };
    let input;
    if (control) {
      // A value near a stored value takes the stored value.
      const at = (x) => seen.find((t) => Math.abs(inUnit(control, t) - x) <= 0.01 * (control.hi - control.lo)) ||
        fmt(x) + (control.unit ? ' ' + control.unit : '');
      input = el('input', { type: 'range', id: 'slider-' + p, 'data-param': p, min: control.lo, max: control.hi, step: 'any',
        value: inUnit(control, S.values[p]), oninput: (ev) => set(at(+ev.target.value), false), onchange: (ev) => set(at(+ev.target.value), true) });
    } else {
      input = el('input', { type: 'range', id: 'slider-' + p, 'data-param': p, min: 0, max: seen.length - 1, step: 1,
        value: Math.max(0, seen.findIndex((x) => sameValue(x, S.values[p]))),
        oninput: (ev) => set(seen[+ev.target.value], false), onchange: (ev) => set(seen[+ev.target.value], true) });
    }
    box.append(el('div', { class: 'slider' },
      el('label', { for: 'slider-' + p }, el('b', {}, p), ' = ', out, ' ', param ? lineLink(param.line) : null),
      input,
      param && param.doc ? el('div', { class: 'doc' }, param.doc) : null));
  }
}

// The params of the values `values` that the page solve cannot change, and
// that are not at their default.
function fixedParams(values) {
  return Object.keys(values).filter((p) => !solveParams.some((q) => q.name === p) &&
    !M.params.some((x) => x.name === p && sameValue(x.value, values[p])));
}

// Shows the case of the values of the sliders. Without a case, `solve`
// meshes and solves them in the page.
async function showValues(solve) {
  const values = { ...S.values };
  let c = await caseFor(values);
  const fixed = fixedParams(values);
  const inPage = SOLVE && fixed.length === 0;
  if (!c && solve && inPage) {
    c = await solveCase(values);
    if (JSON.stringify(values) !== JSON.stringify(S.values)) return;    // the sliders moved on
  }
  if (c) {
    if (D.cases.indexOf(c) !== S.caseIndex) selectCase(D.cases.indexOf(c));
  } else if (SOLVE && !inPage) {
    document.getElementById('note').textContent = 'no stored case has these values, and the page solves only with the default of ' + fixed.join(', ');
  } else if (!SOLVE) {
    document.getElementById('note').textContent = 'no stored case has these values' + (typeof tooLarge === 'function' ? '. ' + tooLarge() : '');
  } else if (!solve) {
    document.getElementById('note').textContent = 'no stored case has these values: release the slider to solve them in the page';
  }
}

function focusSlider(p) {
  const s = document.getElementById('slider-' + p);
  if (s) { s.scrollIntoView({ block: 'center', behavior: 'smooth' }); s.focus(); }
}

function selectCase(i) {
  if (i < 0 || i >= D.cases.length) return;
  // The sliders stay if they have the values of the case. Reason: a new
  // slider element stops the drag of the old one.
  const v = caseValues(i);
  const same = Object.keys(v).every((p) => p in S.values && sameValue(S.values[p], v[p]));
  S.caseIndex = i;
  S.values = v;
  S.step = Math.max(0, Math.min(S.step, D.cases[i].results.length - 1));
  S.foreign = null;
  if (!same) renderSliders();
  renderTables();
  load();
}

// --- the tables --------------------------------------------------------------------

// A table of TOON rows. `link(col, value, row)` gives the cell content.
function table(title, t, link, rowAttrs) {
  if (!t || !t.rows.length) return null;
  return el('div', { class: 'table' }, el('h3', {}, title),
    el('table', {}, el('thead', {}, el('tr', {}, t.cols.map((c) => el('th', {}, link ? link(c, null, null) : c)))),
      el('tbody', {}, t.rows.map((r, i) => el('tr', rowAttrs ? rowAttrs(r, i) : {}, r.map((v, k) => el('td', {}, link ? link(t.cols[k], v, r) : v)))))));
}

function itemLine(kind, name) {
  const list = { step: M.steps, mesh: M.meshes, physics: M.physics, param: M.params, expect: M.expects, method: M.methods };
  const it = (list[kind] || []).find((x) => x.name === name);
  return it ? it.line : null;
}

function reportLine(step, expr) {
  const r = M.reports.find((x) => x.expr === expr && (!x.step || x.step === step)) || M.reports.find((x) => x.expr === expr);
  return r ? r.line : null;
}

// The source line of a verify row: an expect, or the step of a check.
function verifyLine(row, cols) {
  const id = row[cols.indexOf('id')];
  const value = row[cols.indexOf('value')] || '';
  if (id === 'expect') { const n = value.split(':')[0]; return itemLine('expect', n); }
  if (id === 'mesh_quality') return M.meshes.length ? M.meshes[0].line : null;
  return M.steps.length ? M.steps[0].line : null;
}

function withLine(text, line) { return line ? el('span', {}, text, ' ', lineLink(line, '↗')) : text; }

function renderTables() {
  const box = document.getElementById('tables');
  box.textContent = '';
  const add = (...kids) => box.append(...kids.filter((k) => k !== null));
  const toon = S.foreign ? null : caseToons[S.caseIndex];
  if (S.foreign) { add(el('p', { class: 'note' }, 'This view shows ', S.foreign, ', a VTU file of another solver: it has no source links and no solver tables.')); return; }
  const c = D.cases[S.caseIndex];
  add(el('h2', {}, c.values ? `Solve in the page (${c.seconds.toFixed(1)} s, exit ${c.code})` : studyTable ? `Case ${c.id}` : 'Run'));
  if (c.failure) add(el('p', { class: 'note' }, c.failure));
  const v = caseValues(S.caseIndex);
  if (Object.keys(v).length) add(el('p', {}, Object.entries(v).map(([p, t]) => `${p} = ${t}`).join(', ')));
  // the solver
  const method = checkToon.tables.method;
  add(table('Method', method, (c, v, r) => (v === null ? c : c === 'physics' ? withLine(v, itemLine('physics', v)) :
    c === 'mesh' ? withLine(v, itemLine('mesh', v)) : v)));
  if (!method) {
    const lines = M.methods.map((m) => M.items.find((it) => it.line === m.line));
    if (lines.length) add(el('div', { class: 'table' }, el('h3', {}, 'Method'), lines.map((it) => el('p', {}, el('code', {}, it.text), ' ', lineLink(it.line, '↗')))));
  }
  const solveText = M.solves.length ? M.solves.map((s) => el('p', {}, el('code', {}, 'solve { ' + s.text + ' }'), ' ', lineLink(s.line, '↗')))
    : [el('p', {}, 'The model has no solve item, so the solver uses the defaults of spec §14.4: CG for a symmetric positive definite system, else GMRES, with the AMG preconditioner.')];
  add(el('div', { class: 'table' }, el('h3', {}, 'Linear solver'), solveText));
  add(table('Mesh (of the default param values)', verifyToon.tables.mesh, (c, v) => (v === null ? c : c === 'name' ? withLine(v, itemLine('mesh', v)) : v)));
  add(table('Solve', toon.tables.solve, (c, v) => (v === null ? c : c === 'step' ? withLine(v, itemLine('step', v)) : v)));
  add(table('Devices', toon.tables.device, (c, v) => (v === null ? c : c === 'step' ? withLine(v, itemLine('step', v)) : v)));
  const verify = { cols: ['stage', 'id', 'result', 'value', 'limit'], rows: [] };
  if (verifyToon.tables.verify) verify.rows.push(...verifyToon.tables.verify.rows);
  if (toon.tables.verify) verify.rows.push(...toon.tables.verify.rows);
  add(table('Verification', verify, (c, v, r) => {
    if (v === null) return c;
    if (c === 'result') return el('span', { class: v === 'pass' ? 'pass' : 'fail' }, v);
    if (c === 'id') return withLine(v, verifyLine(r, verify.cols));
    return v;
  }, (r) => {
    const id = r[1];
    return id === 'expect' ? { 'data-expect': (r[3] || '').split(':')[0] } : { 'data-check': id };
  }));
  add(table('Reports', toon.tables.report, (c, v, r) => (v === null ? c : c === 'expr' ? withLine(v, reportLine(r[0], v)) : v),
    (r) => ({ 'data-report': r[1] })));
  if (c.values) add(el('div', { class: 'table' }, el('h3', {}, 'Records of the run'), el('pre', {}, c.run + (c.err ? '\n' + c.err : ''))));
  if (studyTable) {
    add(table('Study table', studyTable, (c, v) => {
      if (v !== null) return c === 'status' || c === 'verify' ? el('span', { class: v === 'ok' || v === 'pass' ? 'pass' : 'fail' }, v) : v;
      const line = itemLine('param', c) || itemLine('expect', c) || (M.reports.find((x) => x.expr.startsWith(c + ' ') || x.expr === c) || {}).line;
      return withLine(c, line);
    }, (r, i) => ({ class: i === S.caseIndex ? 'current' : '', 'data-case': i, onclick: () => selectCase(i) })));
    add(expectTable());
  }
  const notes = M.lines.map((t, i) => [t.trim(), i + 1]).filter(([t]) => t.startsWith('///'));
  if (notes.length) {
    add(el('div', { class: 'table' }, el('h3', {}, 'Notes and assumptions of the model'),
      el('ul', {}, notes.map(([t, line]) => el('li', {}, t.slice(3).trim(), ' ', lineLink(line, '↗'))))));
  }
}

// The table of each expect of the study: the value of each case and its
// limit (spec §15.3).
function expectTable() {
  const names = M.expects.filter((e) => e.study === D.study).map((e) => e.name);
  if (!names.length) return null;
  const t = { cols: ['case'].concat(studyParams, ['expect', 'value', 'expected', 'result']), rows: [] };
  D.cases.forEach((c, i) => {
    const v = caseToons[i].tables.verify;
    if (!v) return;
    for (const r of v.rows) {
      const name = (r[3] || '').split(':')[0];
      if (r[1] !== 'expect' || !names.includes(name)) continue;
      const vals = caseValues(i);
      t.rows.push([String(c.id)].concat(studyParams.map((p) => vals[p]), [name, (r[3] || '').slice(name.length + 1).trim(), r[4], r[2]]));
    }
  });
  return table('Measured and expected values of the study', t, (c, v) => {
    if (v === null) return withLine(c, c === 'expect' ? itemLine('expect', names[0]) : itemLine('param', c));
    return c === 'result' ? el('span', { class: v === 'pass' ? 'pass' : 'fail' }, v) : v;
  }, (r) => ({ 'data-case': +r[0] - 1, class: +r[0] - 1 === S.caseIndex ? 'current' : '' }));
}

// --- the view ----------------------------------------------------------------------

// The marked triangles, cells and nodes of S.mark on `mesh`.
function marked(mesh) {
  const none = { tris: new Set(), cells: new Set(), nodes: new Set() };
  if (!S.mark || S.foreign) return none;
  if (S.mark.all) return { tris: new Set(Array.from({ length: mesh.triCell.length }, (_, i) => i)), cells: new Set(), nodes: new Set() };
  const note = document.getElementById('note');
  try {
    if (S.mark.point) {
      const env = letEnv(M, caseValues(S.caseIndex), 1e-5 * mesh.size);
      const p = evaluate(parseExpr(tokens(S.mark.point), M.metre), env);
      return { tris: new Set(), cells: new Set(), nodes: new Set([nearestNode(mesh, p.point || p)]) };
    }
    const key = S.caseIndex + ':' + S.mark.name;
    if (!mesh.regions.has(key)) mesh.regions.set(key, regionOf(M, mesh, M.regions.get(S.mark.name), caseValues(S.caseIndex)));
    if (note) note.textContent = '';
    return mesh.regions.get(key);
  } catch (e) {
    if (note) note.textContent = 'cannot show ' + (S.mark.name || S.mark.point) + ': ' + e.message;
    return none;
  }
}

function currentField(mesh) {
  return mesh.fields.find((f) => f.name === S.field) || mesh.fields[0] || null;
}

// Builds the vertex arrays of the surface, the cut, the mesh lines and the
// marked points, and draws them.
function refresh() {
  const mesh = S.mesh;
  if (!mesh || !view) return;
  const field = currentField(mesh);
  if (field) S.field = field.name;
  const vals = field ? field.values() : null;
  const p = mesh.vtu.points;
  const u = mesh.disp;
  const mark = marked(mesh);
  const valAt = (node, cell) => (vals ? (field.cell ? vals[cell] : vals[node]) : 0);
  let lo = Infinity; let hi = -Infinity;
  if (vals) for (let i = 0; i < vals.length; i++) { if (vals[i] < lo) lo = vals[i]; if (vals[i] > hi) hi = vals[i]; }
  if (!(hi > lo)) { hi = lo + (Math.abs(lo) || 1) * 1e-9; }
  view.range = [lo, hi];
  const nt = mesh.triCell.length;
  const pos = new Float32Array(9 * nt); const disp = new Float32Array(9 * nt);
  const val = new Float32Array(3 * nt); const mk = new Float32Array(3 * nt);
  for (let t = 0; t < nt; t++) {
    const c = mesh.triCell[t];
    const m = mark.tris.has(t) || mark.cells.has(c) ? 1 : 0;
    for (let j = 0; j < 3; j++) {
      const n = mesh.tris[3 * t + j];
      for (let k = 0; k < 3; k++) { pos[9 * t + 3 * j + k] = p[3 * n + k]; disp[9 * t + 3 * j + k] = u ? u[3 * n + k] : 0; }
      val[3 * t + j] = valAt(n, c);
      mk[3 * t + j] = m;
    }
  }
  view.setPart('surface', pos, disp, val, mk, vals ? 0 : 1);
  // the beams, in the colours of the field
  const bp = []; const bd = []; const bv = []; const bm = [];
  for (const c of mesh.beams) {
    const m = mark.cells.has(c) ? 1 : 0;
    for (const n of mesh.cells[c].nodes) {
      for (let k = 0; k < 3; k++) { bp.push(p[3 * n + k]); bd.push(u ? u[3 * n + k] : 0); }
      bv.push(valAt(n, c));
      bm.push(m);
    }
  }
  view.setPart('beams', Float32Array.from(bp), Float32Array.from(bd), Float32Array.from(bv), Float32Array.from(bm), vals ? 0 : 1);
  // the cut
  view.clip = null;
  S.cutMesh = null;
  if (S.cut !== 'off') {
    const a = 'xyz'.indexOf(S.cut);
    const n = [0, 0, 0]; n[a] = 1;
    const d = mesh.lo[a] + S.cutAt * (mesh.hi[a] - mesh.lo[a]);
    const cut = cutPlane(mesh, n, d);
    S.cutMesh = cut;
    const nv = cut.verts.length;
    const cp = new Float32Array(3 * nv); const cd = new Float32Array(3 * nv); const cv = new Float32Array(nv); const cm = new Float32Array(nv);
    cut.verts.forEach(([ea, eb, s], i) => {
      const cell = cut.cells[Math.floor(i / 3)];
      for (let k = 0; k < 3; k++) {
        cp[3 * i + k] = p[3 * ea + k] + s * (p[3 * eb + k] - p[3 * ea + k]);
        cd[3 * i + k] = u ? u[3 * ea + k] + s * (u[3 * eb + k] - u[3 * ea + k]) : 0;
      }
      cv[i] = field && field.cell ? vals[cell] : valAt(ea, cell) + s * (valAt(eb, cell) - valAt(ea, cell));
      cm[i] = mark.cells.has(cell) ? 1 : 0;
    });
    view.setPart('cut', cp, cd, cv, cm, vals ? 0 : 1);
    view.clip = [n[0], n[1], n[2], d];
  } else {
    view.setPart('cut', new Float32Array(0), new Float32Array(0), new Float32Array(0), new Float32Array(0), 0);
  }
  // the mesh lines: the edges of the boundary faces
  const seen = new Set();
  const lp = []; const ld = [];
  for (const f of mesh.faces) {
    for (let j = 0; j < f.length; j++) {
      const a = f[j]; const b = f[(j + 1) % f.length];
      const key = a < b ? a * 1e7 + b : b * 1e7 + a;
      if (seen.has(key)) continue;
      seen.add(key);
      for (const n of [a, b]) for (let k = 0; k < 3; k++) { lp.push(p[3 * n + k]); ld.push(u ? u[3 * n + k] : 0); }
    }
  }
  view.setPart('lines', Float32Array.from(lp), Float32Array.from(ld), new Float32Array(lp.length / 3), new Float32Array(lp.length / 3), 1);
  // the marked nodes and the probe
  const pts = []; const pd = [];
  for (const n of mark.nodes) for (let k = 0; k < 3; k++) { pts.push(p[3 * n + k]); pd.push(u ? u[3 * n + k] : 0); }
  if (S.probe) for (let k = 0; k < 3; k++) { pts.push(S.probe.at[k]); pd.push(S.probe.disp[k]); }
  view.setPart('points', Float32Array.from(pts), Float32Array.from(pd), new Float32Array(pts.length / 3), new Float32Array(pts.length / 3), 3);
  renderLegend(field, lo, hi);
  const fs = document.getElementById('field');
  if (fs && field) fs.value = field.name;
  draw();
}

function draw() {
  if (!view) return;
  view.scale = S.scale;
  view.bands = S.bands;
  view.lines = S.lines;
  view.draw();
}

function renderLegend(field, lo, hi) {
  const box = document.getElementById('legend');
  box.textContent = '';
  if (!field) return;
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((t) => el('span', {}, fmt(lo + t * (hi - lo))));
  box.append(el('div', { class: 'legend-name' }, field.name), el('div', { class: 'bar' }), el('div', { class: 'ticks' }, ticks));
}

// Loads the mesh of the current case and step, and draws it.
// A newer load stops an older load that still decodes its mesh.
let loads = 0;

async function load() {
  const c = currentCase();
  if (!c || !c.results.length) { document.getElementById('note').textContent = 'this case has no result file'; return; }
  const r = c.results[Math.min(S.step, c.results.length - 1)];
  const first = !S.mesh;
  const me = ++loads;
  document.body.dataset.ready = 'no';
  const mesh = await meshOf(r);
  if (me !== loads) return;
  S.mesh = mesh;
  S.probe = null;
  document.getElementById('probe').textContent = '';
  if (first) { view.fit(S.mesh.lo, S.mesh.hi); S.scale = autoScale(S.mesh); }
  renderControls();
  refresh();
  document.body.dataset.case = String(c.id);
  document.body.dataset.ready = 'yes';
}

// The probe at the pixel (x, y): the point, the value of the field, and
// the regions of the model that hold it.
function probe(x, y) {
  const hit = view.pickAt(x, y);
  const box = document.getElementById('probe');
  if (!hit || (hit.part === 'cut' && !S.cutMesh)) { S.probe = null; box.textContent = ''; refresh(); return; }
  const mesh = S.mesh;
  const t = hit.tri;
  const p = mesh.vtu.points; const u = mesh.disp;
  // Each corner of the triangle as [node, weight] pairs: a node of the
  // surface, or a point on an edge of a cell for the cut.
  const corners = [0, 1, 2].map((j) => {
    if (hit.part === 'surface') return [[mesh.tris[3 * t + j], 1]];
    const [ea, eb, f] = S.cutMesh.verts[3 * t + j];
    return [[ea, 1 - f], [eb, f]];
  });
  const sum = (corner, arr, n, k) => corner.reduce((acc, [node, f]) => acc + f * arr[n * node + k], 0);
  const at = (corner) => [0, 1, 2].map((k) => sum(corner, p, 3, k) + S.scale * (u ? sum(corner, u, 3, k) : 0));
  const [a, b, c] = corners.map(at);
  const [o, d] = view.ray(x, y);
  const e1 = [0, 1, 2].map((k) => b[k] - a[k]); const e2 = [0, 1, 2].map((k) => c[k] - a[k]);
  const cross = (q, r) => [q[1] * r[2] - q[2] * r[1], q[2] * r[0] - q[0] * r[2], q[0] * r[1] - q[1] * r[0]];
  const dot = (q, r) => q[0] * r[0] + q[1] * r[1] + q[2] * r[2];
  const h = cross(d, e2); const det = dot(e1, h) || 1e-30;
  const s = [0, 1, 2].map((k) => o[k] - a[k]);
  const w1 = Math.min(1, Math.max(0, dot(s, h) / det));
  const w2 = Math.min(1 - w1, Math.max(0, dot(d, cross(s, e1)) / det));
  const w = [1 - w1 - w2, w1, w2];
  const field = currentField(mesh);
  const vals = field ? field.values() : null;
  const cell = hit.part === 'surface' ? mesh.triCell[t] : S.cutMesh.cells[t];
  const value = vals ? (field.cell ? vals[cell] : w.reduce((acc, wi, j) => acc + wi * sum(corners[j], vals, 1, 0), 0)) : null;
  const point = [0, 1, 2].map((k) => w.reduce((acc, wi, j) => acc + wi * sum(corners[j], p, 3, k), 0));
  const disp = [0, 1, 2].map((k) => w.reduce((acc, wi, j) => acc + wi * (u ? sum(corners[j], u, 3, k) : 0), 0));
  S.probe = { at: point, disp, value };
  box.textContent = '';
  box.append(el('div', {}, `(${point.map(fmt).join(', ')})`));
  if (field) box.append(el('div', {}, el('b', {}, field.name), ' = ', el('span', { id: 'probe-value' }, fmt(value))));
  const regions = S.foreign ? null : probeRegions(mesh, hit.part === 'surface' ? t : -1, cell);
  if (regions) box.append(regions);
  refresh();
}

// The links to the regions that hold triangle t (cell c).
function probeRegions(mesh, t, c) {
  const links = [];
  for (const r of M.regions.values()) {
    if (r.type === 'Point') continue;
    const key = S.caseIndex + ':' + r.name;
    try {
      if (!mesh.regions.has(key)) mesh.regions.set(key, regionOf(M, mesh, r, caseValues(S.caseIndex)));
    } catch (e) { continue; }
    const set = mesh.regions.get(key);
    if (set.tris.has(t) || set.cells.has(c)) links.push(lineLink(r.line, r.name));
  }
  return links.length ? el('div', { class: 'regions' }, 'in ', links.flatMap((l, i) => (i ? [', ', l] : [l]))) : null;
}

function bindMouse(canvas) {
  let drag = null;
  canvas.addEventListener('pointerdown', (ev) => {
    drag = { x: ev.clientX, y: ev.clientY, moved: false, pan: ev.button === 2 || ev.shiftKey };
    canvas.setPointerCapture(ev.pointerId);
  });
  canvas.addEventListener('pointermove', (ev) => {
    if (!drag) return;
    const dx = ev.clientX - drag.x; const dy = ev.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
    drag.x = ev.clientX; drag.y = ev.clientY;
    if (drag.pan) {
      const f = view.dist / canvas.clientHeight * 0.6;
      const cy = Math.cos(view.yaw); const sy = Math.sin(view.yaw);
      const sp = Math.sin(view.pitch); const cp = Math.cos(view.pitch);
      const right = [-sy, cy, 0]; const up = [-cy * sp, -sy * sp, cp];
      for (let k = 0; k < 3; k++) view.target[k] += f * (-dx * right[k] + dy * up[k]);
    } else {
      view.yaw -= dx * 0.008;
      view.pitch = Math.max(-1.55, Math.min(1.55, view.pitch + dy * 0.008));
    }
    draw();
  });
  canvas.addEventListener('pointerup', (ev) => {
    if (drag && !drag.moved) { const r = canvas.getBoundingClientRect(); probe(ev.clientX - r.left, ev.clientY - r.top); }
    drag = null;
  });
  canvas.addEventListener('contextmenu', (ev) => ev.preventDefault());
  canvas.addEventListener('wheel', (ev) => { ev.preventDefault(); view.dist *= Math.exp(ev.deltaY * 0.001); draw(); }, { passive: false });
}

// Opens a VTU file of another solver. Its view has no source links.
async function openFile(ev) {
  const file = ev.target.files[0];
  if (!file) return;
  ++loads;
  try {
    let bytes = new Uint8Array(await file.arrayBuffer());
    if (bytes[0] === 0x1f && bytes[1] === 0x8b) bytes = await inflate(bytes, 'gzip');
    S.mesh = buildMesh(await parseVtu(bytes));
    S.foreign = file.name;
    S.mark = null; S.probe = null; S.field = null;
    view.fit(S.mesh.lo, S.mesh.hi);
    S.scale = autoScale(S.mesh);
    renderControls();
    renderTables();
    refresh();
  } catch (e) {
    document.getElementById('note').textContent = 'cannot open ' + file.name + ': ' + e.message;
  }
}

function renderHeader() {
  document.getElementById('title').textContent = D.model || 'Skein report';
  const intro = [];
  for (const line of M.lines) {
    const t = line.trim();
    if (t.startsWith('//')) intro.push(t.replace(/^\/+\s?/, ''));
    else if (t) break;
  }
  document.getElementById('intro').textContent = intro.join(' ');
  document.getElementById('file').textContent = D.file ? `${D.file}: the results of skein run` + (D.study ? `, study ${D.study} (${D.cases.length} cases)` : '') : '';
}

async function main() {
  renderHeader();
  renderSource();
  if (studyTable) {
    // the first case with the default value of each varied param, or case 1
    const defaults = {};
    for (const p of M.params) defaults[p.name] = p.value;
    const c = await caseFor(defaults);
    S.caseIndex = c ? D.cases.indexOf(c) : 0;
  }
  S.values = caseValues(S.caseIndex);
  renderSliders();
  renderTables();
  try {
    view = new View(document.getElementById('canvas'));
    bindMouse(document.getElementById('canvas'));
    window.addEventListener('resize', draw);
    await load();
  } catch (e) {
    document.getElementById('stage').append(el('p', { class: 'fail' }, 'The viewer cannot start: ' + e.message));
    document.body.dataset.ready = 'error';
  }
}

// The functions that a test or a browser solve uses.
window.skein = { caseFor, solveCase, selectCase, parseVtu, parseToon, probe, get view() { return view; }, state: S, data: D, model: M };
main();
