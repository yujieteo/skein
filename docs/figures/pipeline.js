// The path of a Skein model from source to a native executable or a browser page.
// Drawn with the schematic-emil style of the schematic-figure skill (emil.js is that style's style.js).
// Render: put the schematic-emil style.js (as emil.js) and its fonts/ folder next to this file, then run
//   node <schematic-figure>/foundation/tools/render.mjs --root <dir> --fig pipeline.js
// Hero: none, a software pipeline has no physical object.
import * as E from './emil.js';

export const figure = { width: 80, height: 85.5, name: 'pipeline', format: 'journal' };

export default function draw(F, L) {
  const { r, line, text } = E.useEmil(F, L, 'web');
  const m = L.mm, root = F.root, G = E.gray, B = E.hue.blue;
  const skein = { fill: B[3], stroke: B[7], r: r.card, width: line.thin };   // Skein code, in src/
  const tool = { fill: G[100], stroke: G[700], r: r.card, width: line.thin }; // an external tool
  const out = { fill: G[100], stroke: G[400], r: r.card, width: line.thin };  // an output file

  const name = (x, y, s, o = {}) => L.text(root, m(x), m(y), s, { size: text.label, weight: E.weight.name, color: G[1200], ...o });
  const note = (x, y, s, o = {}) => L.text(root, m(x), m(y), s, { size: text.small, color: G[1100], ...o });
  const card = (x, y, w, h, s) => L.rect(root, m(x), m(y), m(w), m(h), s);
  const arrow = (a, b, o = {}) => L.arrow(root, a.map(m), b.map(m), { stroke: G[900], width: line.normal, ...o });

  // main column x 2..45, centre 23.5; library column x 54..78
  const cx = 23.5;

  // 1. the model
  card(2, 2, 43, 9, out);
  name(cx, 6.5, 'model.skn', { anchor: 'middle', valign: 'middle' });
  arrow([cx, 11], [cx, 15.5]);

  // 2. the front end
  card(2, 15.5, 43, 13, skein);
  name(5, 19.5, 'Front end', { valign: 'middle' });
  note(5, 23.6, 'src/syntax: parse', { valign: 'middle' });
  note(5, 26.4, 'src/check: types, units', { valign: 'middle' });
  arrow([cx, 28.5], [cx, 33]);

  // 3. the code generator
  card(2, 33, 43, 10, skein);
  name(5, 37, 'Code generator', { valign: 'middle' });
  note(5, 40.3, 'src/codegen', { valign: 'middle' });
  arrow([cx, 43], [cx, 51.5]);
  note(cx + 2, 47.2, 'main.goose', { valign: 'middle', color: B[11] });

  // 4. the solver library, which main.goose imports
  card(54, 15.5, 24, 46, skein);
  name(66, 19.5, 'Library', { anchor: 'middle', valign: 'middle' });
  const parts = ['runtime', 'solvers', 'elements', 'mesher', 'shells', 'transient', 'flow', 'cons', 'meshless', 'optim'];
  note(57, 23.6, 'src/', { valign: 'middle' });
  parts.forEach((p, i) => note(58.5, 27 + i * 3.35, p, { valign: 'middle' }));

  // 5. Goose
  card(2, 51.5, 43, 10, tool);
  name(5, 55.5, 'Goose', { valign: 'middle' });
  note(5, 58.8, 'writes C', { valign: 'middle' });
  arrow([54, 57.5], [45, 57.5]);
  note(49.5, 54.6, 'imports', { anchor: 'middle', valign: 'middle' });

  // 6. the two targets
  L.polyline(root, [[cx, 61.5], [cx, 63.5], [19, 63.5], [59, 63.5], [cx, 63.5]].map(([x, y]) => [m(x), m(y)]), { stroke: G[900], width: line.normal, fill: 'none' });
  arrow([19, 63.5], [19, 67.5]);
  arrow([59, 63.5], [59, 67.5]);
  card(2, 67.5, 36, 15.5, out);
  name(4.5, 71.5, 'Native executable', { valign: 'middle' });
  note(4.5, 75.2, 'skein run', { valign: 'middle', color: B[11] });
  note(4.5, 78.4, 'TOON tables, VTU files', { valign: 'middle' });
  card(40, 67.5, 38, 15.5, out);
  name(42.5, 71.5, 'HTML page', { valign: 'middle' });
  note(42.5, 75.2, 'skein build --target wasm', { valign: 'middle', color: B[11] });
  note(42.5, 78.4, 'C compiled by emcc', { valign: 'middle' });

  E.emilTheme(F);
}
