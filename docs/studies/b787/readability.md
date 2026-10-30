# Readability: Skein and the Abaqus deck for the B787 study

Two writers modelled the same [problem sheet](problem.md) separately: once
in Skein ([`skein/`](skein/)) and once as an Abaqus deck with a mesh script
([`abaqus/`](abaqus/)). This page compares the two.

## Result

- **Size.** The Skein model is 287 code lines and 5,719 tokens. The Abaqus
  input is 164 deck lines and 349 script lines, 9,329 tokens in all. The
  script generates 55,836 more lines (1,067,377 tokens).
- **Edits.** For the 4 edits, Skein changed 8 to 25 lines (added plus
  removed). The Abaqus input changed 6 to 72 hand-written lines, and 0 to
  98,196 generated lines. The window row is the largest edit in both
  formats: 25 Skein lines, and 72 lines of `mesh.py`.
- **Check.** `skein check` stops each Skein file with 15 errors. Each error
  comes from one of the gaps G3, G31, G32 and G33 that line 1 of the file
  names. No edit added a new error id. The edits found 2 new gaps, G34 and
  G35. The [coverage record](../../femahoy_coverage.md#gaps-of-the-787-study)
  holds them. Before, N1 of `skein/notes.md` became G31. N2 and N3 stay on that page.
- **LLM runs.** All 3 Skein runs gave a result, and 1 of 3 Abaqus-path runs
  gave a result. But no run modelled the barrel and also gave a correct
  result. Each Skein run first wrote a barrel of 3D solids that passed
  `skein check`. Then `skein run` (one box only) made each run reduce the
  model to a flat strip of skin. The strip results agree with their hand
  values, and `skein run` checks each `expect` of the strip. 2
  Abaqus-path runs sent a Python mesh script to ccx in each round. The
  third run solved a coarse barrel deck, but the pressure acts
  inward, and the deck fails C3, C4 and C5. Nothing in that path told the
  model.
- **Solvers.** `skein run` meshes only one box, so no Skein model of the
  barrel runs. CalculiX (`fh run`) cannot read the hand-written Abaqus deck
  (`XSYMM`, `*BEAM GENERAL SECTION` with `*CENTROID`, `*SUBMODEL`). Thus
  neither hand-written model gives a result on this machine.

## Sources

The models are at `origin/main` `2de58a4`. Each number comes from one of
these commits and commands.

| Item | Source |
|---|---|
| Skein model | `skein/fuselage.skn`, `skein/door_corner.skn` at `2de58a4` (PR 12, gap marks of PR 15, C2 fix of PR 24) |
| Abaqus deck and mesh script | `abaqus/*.inp`, `abaqus/mesh.py` at `2de58a4` (PR 13, C2 fix of PR 24) |
| `skein` | `./build.sh` at `2de58a4` |
| Generated mesh | `python3 abaqus/mesh.py`. The SHA-256 values agree with [`abaqus/notes.md`](abaqus/notes.md). |
| Size | [`compare/size.py`](compare/size.py), tiktoken 0.14.0, `o200k_base` |
| Edits | [`compare/edits/`](compare/edits/), `git diff --no-index --numstat` against the files at `2de58a4` |
| LLM runs | [`compare/llm_run.py`](compare/llm_run.py) on `scripts/llm/workflow.py` at `fe4fe23` (PR 14), `skein` built at `7ffbeb6`, sheet at `5d032bf` |

## 1. `skein check` on the Skein model

`build/skein check skein/fuselage.skn` and `build/skein check
skein/door_corner.skn` each print `summary: errors=15 lints=0`. Each error
is under a `// gap` comment (PR 15).

| Error | Gap | `fuselage.skn` | `door_corner.skn` |
|---|---|---|---|
| `reissner_mindlin` is not a template (`unknown_name`) | G3 | 1 | 1 |
| `N`, `eps_m` are not declared, `N*n` is not a flux, `global.e_hoop` is not a field (`unknown_name`, `bc_form`) | G33 | 7 | 8 |
| `timoshenko` is not a template, `Nx` is not declared (`unknown_name`) | G31 | 4 | 1 |
| a `Surface` or `Curve` made with `&` is a `Solid` (`geometry_dim`) | G32 | 3 | 5 |
| Total | | 15 | 15 |

The gaps are in the
[femahoy coverage](../../femahoy_coverage.md#gaps-of-the-787-study). G32 is
a defect of the checker as much as of the spec: `geom_op_ty` in
`src/check/expr.goose` gives `a & b` the higher dimension.

## 2. Size

Lines are non-blank lines. Code lines also drop the full-line comments
(`//`, `**`, `#`). Command: `compare/size.py FILE...`.

| Part | Files | Lines | Code lines | Tokens |
|---|---|---|---|---|
| Skein model | `fuselage.skn`, `door_corner.skn` | 391 | 287 | 5,719 |
| Abaqus deck, hand-written | `global.inp`, `door_corner.inp`, `properties.inp` | 240 | 164 | 3,639 |
| Mesh script | `mesh.py` | 352 | 349 | 5,690 |
| Generated mesh and loads | `*_mesh.inp`, `*_loads.inp` | 55,847 | 55,836 | 1,067,377 |

The Skein files repeat about 165 lines, because the parser rejects a `use`
library file (see [`skein/notes.md`](skein/notes.md)). The Abaqus deck has
72 ply lines out of 164 code lines.

The Skein model also holds what the deck keeps in comments. The deck gives
each check value as a number in a comment, for example `C1: SF2 / t at P in
48.8 to 61.3 MPa`. The Skein model computes the same value in an `expect`.

## 3. LLM runs

Each run is a new `claude-opus-5-5` session with no tools. The input is
one line plus the whole sheet:

> Model the Boeing 787 forward fuselage pressurization of this problem
> sheet, and report the outputs and the checks it asks for.

The Skein path also gets `skein help --llm` first. The driver feeds the
output of a failed `skein check`, a failed `skein run` or a failed `fh run`
back to the model, for 5 rounds at most. Command:
`compare/llm_run.py skein|inp --out DIR`. The outputs of each round are in
[`compare/llm/`](compare/llm/), with `results-*.toon`.

| Run | Rounds | Check pass at round | Result | Tokens in | Tokens out | Wall time |
|---|---|---|---|---|---|---|
| Skein 1 | 5 | 2 | strip, reaction −28.29 kN, 5 `expect` pass | 283,521 | 71,726 | 673 s |
| Skein 2 | 4 | 2 | strip, reaction −28.29 kN, 1 `expect` pass | 219,128 | 71,391 | 678 s |
| Skein 3 | 5 | 2 | strip, reaction −28.29 kN, 3 `expect` pass | 295,332 | 79,582 | 1,058 s |
| Abaqus 1 | 5 | never | none | 492,306 | 164,668 | 1,385 s |
| Abaqus 2 | 5 | 5 | coarse barrel, fails C3, C4, C5 | 544,186 | 230,109 | 2,055 s |
| Abaqus 3 | 5 | never | none | 524,611 | 172,331 | 1,483 s |

For the Abaqus path, check and run are one `fh run`. No run passed on the
first round.

**Skein runs.** In round 1, each model was a quarter barrel of 3D solids
(314 to 399 lines). The models said that Skein has no shell, laminate or
beam template (G3, G31), and used solids in place of them. The round 1
errors were `unit_fractional` (`(y^2 + z^2)^0.5`, where `sqrt` is
available), `duplicate_name` (`t`, `shell`), and `geometry_dim`. 4 of the
`geometry_dim` errors in run 2 are G32. Each model passed `skein check` in
round 2. Then `skein run` stopped with `the solid 'skin' is not one box`.
In the next rounds, each model became a flat strip of skin at P (84 to 101
lines), with a comment that names what it leaves out. The strip reaction
−28.29 kN agrees with p_U R / 2 × 0.2013 m. The hoop displacement
0.0925 mm agrees with the hand value in each model.

**Abaqus runs.** The driver asks for one deck with the mesh in it. In 14 of
15 rounds, the model sent a Python script that writes the deck. `fh run`
gave the script to ccx, and ccx stopped on signal 11. `fh` then said "too
little memory", so the model looked for the wrong cause. ccx itself works
on this machine: a deck with one cube solves in 0.01 s. In round 5 of run
2, the model wrote a coarse deck (2,678 lines, one element for each
stringer pitch, square door corners, isotropic shells, no submodel). It
solves in 1.03 s. The reaction sums in `compare/llm/inp/inp-2/job.dat` are:

| Check | Expected | Run 2 | Pass |
|---|---|---|---|
| C3, sum of x | −1273.0 kN | −1327.6 kN (4.3 %) | no |
| C4, sum of z | 0 | −138.7 kN | no |
| C5, sum of y | −3427.1 kN | +3030.8 kN | no |

The sign of C5 shows that the pressure acts inward. The deck has no
assertion, so nothing in this path finds the error.

**Gap ids hit by the runs:** G3 and G31 (all 3 Skein runs used solids in
place of shells and beams), and G32 (Skein run 2). The `skein run` limit
(one box) stopped each barrel. That is not a gap of the spec.

**Limits of these runs.**

- The runs used the sheet at `5d032bf`, before the C2 fix. No run got as
  far as C2.
- The one line above is not the line in sheet section 8. Each run also
  saw section 8, because it is part of the sheet.
- An earlier set of runs used `workflow.py` before PR 14, which does not
  feed `skein run` errors back. That set is not in this table.

## 4. Edits

Each edit starts from a copy of both models at `2de58a4` and is the
smallest correct change. Each edit also updates the comments and the check
values that the change makes wrong. The diffs are in
[`compare/edits/`](compare/edits/), one patch for each edit and format. To
apply one, run
`git apply -p3 --directory=docs/studies/b787/skein docs/studies/b787/compare/edits/e1_pitch.skein.patch`
(use `abaqus` for an Abaqus patch).

| Edit | Change |
|---|---|
| E1 frame pitch | 0.6096 m to 0.4064 m (16 in), 15 frames in the half length. P stays mid-bay. |
| E2 skin layup | skin [45/0/−45/90]3s to [45/0/−45/90]2s (16 plies, 3.056 mm). The pad-up stays 48 plies. |
| E3 window row | 8 windows 0.272 m × 0.467 m, corner radius 0.100 m, centred at θ = 75° and x = 2 to 9 frame pitches. Each window cuts the skin and the stringers. The pane load goes to the window edge by the rule of the door, with x − x_w in place of x. |
| E4 pressure sweep | 3 cases: 65.0 kPa (limit), 86.5 kPa (1.33 × limit, sheet 6.2) and 97.5 kPa (ultimate) |

Diff lines (+/−), from `git diff --no-index --numstat` of each file. The
generated lines are the 4 generated files after `python3 mesh.py`.

| Edit | Skein | Abaqus deck | `mesh.py` | Other hand-written | Generated |
|---|---|---|---|---|---|
| E1 frame pitch | +5 −5 | +1 −1 | +2 −2 | | +33,695 −28,986 |
| E2 skin layup | +4 −4 | +8 −16 | | | 0 |
| E3 window row | +21 −4 | | +52 −20 | | +47,680 −50,516 |
| E4 pressure sweep | +18 −4 | +40 −0 | +10 −5 | `.gitignore` +1 −1 | +5 −5, and 4 new files of +1,286 |

Gap ids and findings for each edit:

- **E1.** Skein: no new error. The pitch is one `let` in each file, and the
  frame count is a literal (`10`) in each file. The `20-bay` comment also
  changes. Abaqus: `mesh.py` holds the pitch and the count. The C1 range in
  a `global.inp` comment (48.8 MPa) is a number that the writer computed by
  hand. It changes to 44.3 MPa. In Skein, the `expect` computes it.
- **E2.** Skein: no new error. `pad_lam` had `h = 2*t_skin`, so a thinner
  skin also made the pad-up thinner. The edit writes `48*ply_t`. A change of
  the stacking sequence of the skin only is not possible: the layup is an
  argument of the physics, so both materials take one layup (G3). Abaqus:
  13 of the 24 ply lines change, and 3 check values in comments change
  (C2 becomes 42.0 MPa).
- **E3.** Skein: 2 more errors, both G33 (`N` is not declared, `N*n` is not
  a flux). New gap **G35**: a `pattern` gives no value for each copy, so the
  load cannot name the centre of each window. The edit computes x − x_w
  with `atan2(sin(a), cos(a))`, which is correct but hard to read. The
  submodel does not change, because no window is in its region. Abaqus: the
  deck does not change, but `mesh.py` must handle many cutouts. The edit
  generalizes the distance function and the corner snap, adds the window
  edges to the grid, and adds the window load. The first try failed the
  assertion `check point P is not a node`, because the new grid lines moved
  the nodes. The script printed the load resultant Fy = 3427.08 kN, as
  before, so the window load is correct. The smallest shell corner angle
  falls from 15.8° to 11.4°, and no shell is inverted.
- **E4.** Skein: `param` and a 4-line `study` in each file. 2 gaps: **G10**
  (each `expect` runs once, with the default value) and new gap **G34**
  (`results()` takes no study case, so the submodel cannot take the global
  case with the same pressure). Abaqus: 2 more steps in each deck, and
  `mesh.py` writes 1 load file for each step. The `.gitignore` pattern
  `*_loads.inp` does not match the new files, so it changes too. The deck
  gives no expected value for steps 2 and 3. A comment says that they scale
  with p.

The generated files change by 62,681 lines for E1 and 98,196 lines for E3.
No reviewer can read that diff. Thus a review of the Abaqus model must read
`mesh.py`, and must trust that the script and the deck agree.

## 5. Defects found

These are defects in the inputs to this study. This page does not change
the models or the sheet.

1. **C2 of the sheet was wrong.** The earlier sheet gave C2 = 25.2 MPa,
   which leaves out the Poisson term of the skin. PR 24 corrects it to
   28.7 MPa (sheet 7.3), in the sheet, the deck comment and the Skein
   `expect`. The Skein `expect` writes ν as the literal `0.326`, not as the
   material value.
2. **Two versions of the sheet.** PR 12 (the Skein model) held an older
   `problem.md`, and PR 13 held the newer one. Main now has the PR 13
   version. The Skein model computes s_c from the door size, so the newer
   value agrees with it.
3. **The `fh run` message.** When ccx gets a file that is not a deck, it
   stops on signal 11, and `fh` says "too little memory". The LLM runs
   followed this wrong cause for many rounds.
4. **The Abaqus-path driver.** `workflow.py` takes the last code block only.
   The model wrote a script in place of a deck, and the driver gave it to
   ccx without a test of the file kind.

## 6. Reproduce

```sh
./build.sh
build/skein check docs/studies/b787/skein/fuselage.skn
python3 docs/studies/b787/compare/size.py docs/studies/b787/skein/*.skn
docs/studies/b787/compare/llm_run.py skein --out /tmp/b787-llm
docs/studies/b787/compare/llm_run.py inp --out /tmp/b787-llm
```

`size.py` needs `tiktoken`. `llm_run.py` needs the `claude` CLI, and `fh`
for the `inp` path.
