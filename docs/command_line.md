# The skein command line

<details>
<summary>Table of contents</summary>

- [1. Commands](#1-commands)
- [2. Checking a model](#2-checking-a-model)
- [3. Formatting](#3-formatting)
- [4. Lint](#4-lint)
- [5. Running and verifying](#5-running-and-verifying)
- [6. What skein prints](#6-what-skein-prints)
- [7. Building](#7-building)
- [8. The language server](#8-the-language-server)
- [9. Files](#9-files)
- [10. Exit status](#10-exit-status)
- [11. The Markdown document](#11-the-markdown-document)

</details>

This page describes the current commands and the planned interface.
The command table states which commands exist. Notes in each section
separate current support from the design. The
[roadmap](roadmap.md) gives the planned phases.

`lint` checks the model first. `run` supports a limited set of native
models. It reports the verification checks that its solver implements.

---

## 1. Commands

| Command | What it does | Status |
|---|---|---|
| `skein check model.skn` | parses and type-checks. It creates no files and performs no solve. | Current |
| `skein fmt model.skn…` | rewrites the files in the standard layout (§3) | Current |
| `skein lint model.skn` | `check`, then the lint rules | Current |
| `skein verify model.skn` | `check`, the mesh, and the pre-solve verification. No solve. | Current: mesh quality only |
| `skein verify convergence model.skn` | the convergence study (checks.md §4.3) | Current: `fv`, `dg`, and scalar `rkpm` and `gfd` models |
| `skein build model.skn --target wasm` | writes the page of the model, which meshes and solves it in the browser | Current: `--target wasm`, simplex and hex models |
| `skein run model.skn` | builds and runs the model, with the checks of its solver | Current: limited models |
| `skein export model.skn` | writes the geometry of a part to STL or STEP (§7.1) | Current |
| `skein view model.skn` | opens the results in the WebGL viewer | Planned |
| `skein report model.skn` | runs the model, then writes one HTML file with the source, the solver tables and a WebGL viewer of the results (§12) | Current |
| `skein md model.skn` | `check`, then writes the model as a Markdown document (§11) | Current |
| `skein lsp` | serves an editor over the Language Server Protocol | Planned |
| `skein help [command]` | prints help | Current: general help only |
| `skein help --llm` | prints the language card for an LLM: the grammar, the keywords, the standard library and one example | Current |
| `skein help std` | prints the Skein source of the standard library (spec §4.2) | Current |

Output format:

| Option | Meaning |
|---|---|
| `--format toon\|json` | output format for `check`, `lint` and `run --study` (§6). Default `toon`. |

---

## 2. Checking a model

```
skein check samples/01_cantilever.skn
```

`check` applies the rules of the [check catalogue](checks.md#2-errors).
It continues after an error where it can, so that one run reports
as many errors as possible. It also prints the parts of the boundary that
get the default natural condition (spec §12.3), because a missing `bc` is
the most common modelling fault.

The current `check` command accepts one model file and `--format toon|json`.
It checks `param` declarations at their default values.
The options below are planned:

| Option | Meaning |
|---|---|
| `--groups` | print the dimensionless groups of the model (spec §21.3) |
| `--study <name>` | check each case of the study, and print the stages that each parameter reaches (spec §20.4). No solve. |
| `--set name=value` | give a `param` a value for this command (spec §20.1). The value has a unit and must be in the range of the `param`. |

`check` needs no mesh. Its target is under 50 ms for each of the samples.
No measured run confirms this target yet.
Reason: an LLM runs `check` after each edit.

---

## 3. Formatting

`skein fmt` gives each model one layout. Two models with the same meaning
and the same comments have the same text after `fmt`. The layout keeps each
change local, so that the diff of a small change is small (spec §3.2).

| Option | Meaning |
|---|---|
| `--check` | lists the files that `fmt` would change, and exits 1 if any |
| no file, or `-` | formats standard input to standard output |

### 3.1 Format rules

1. Indentation is 4 spaces. Skein has no tabs.
2. A line is at most 100 characters. A longer argument list or study list
   gets one item on each line, indented by 4 spaces, with a comma after
   each item, also after the last. A longer expression breaks after a
   binary operator (spec §2.2). Thus a new item adds one line and changes
   no other.
3. One space goes on each side of `=`, `==`, `!=`, `<`, `<=`, `>`, `>=`,
   `+`, `-`, `|`, `&`, `||`, `&&`, and after `:` and `,`.
4. No space goes around `*`, `/` and `^`. Thus `2*mu*eps` and `kg/m^3`
   look as in a textbook.
5. One space goes between a number and its unit: `210 GPa`, `(0, 0, 1) mm`.
6. The header is one line. The `use` lines come after it.
7. One blank line goes between two top-level items. A run of `bc` lines
   has no blank lines in it. In a run of `let` and `param` lines, `fmt`
   keeps one blank line that the author wrote between two lines. `fmt`
   collapses two or more blank lines to one.
8. In a block, the statements keep their order. `fmt` never reorders items,
   because order carries meaning (declaration before use, step order).
9. `fmt` never aligns columns, in a table or anywhere else. A table row is
   its values with `, ` between them. Reason: alignment changes the nearby
   lines when one value gets longer.
10. A number keeps its digits. `fmt` does not change `1.0` to `1`.
11. Comments keep their text and their line. A trailing comment gets 2
    spaces before `//`.
12. `fmt` changes only the lines that break a rule. It does not reflow a
    paragraph of a comment or a line that already obeys the rules.

---

## 4. Lint

```
skein lint samples/03_birth_death.skn
```

`lint` does `check` and then the rules of [checks.md §3](checks.md#3-lint-rules).
Rules that need the mesh, for example `bending_lock` and
`mesh_feature_small`, run only with `--mesh`, because the mesh takes time.

---

## 5. Running and verifying

Current `run` accepts one model file. With `--study`, it also accepts
`--jobs` and `--format`. The other options below are planned.
Current `verify` accepts only one model file and reports `mesh_quality`.
It does not support `--mesh`. `verify convergence` accepts one model file
and `--levels`, and runs the `fv` and `dg` models and the scalar `rkpm`
and `gfd` models (checks.md §4.3).

The examples and option table below describe the planned interface:

```
skein run samples/02_cavity.skn
skein verify convergence examples/grid/mms_dg_static.skn --levels 3
skein run samples/06_cylinder_sweep.skn --study re_sweep --jobs 3
```

| Option | Commands | Meaning |
|---|---|---|
| `--target wasm` | `build` | the page that solves in the browser (§7) |
| `--threads n` | `run` | worker threads for assembly and matrix products. Default: all hardware threads. With `--study`, a case has all threads when no second case can run beside it, and 1 thread otherwise (spec §20.5). |
| `--levels n` | `verify convergence` | number of meshes. Default 3. |
| `--mesh` | `lint`, `verify` | make the mesh and run the checks that need it |
| `--study <name>` | `run` | run each case of the study (spec §20) |
| `--set name=value` | `run`, `verify`, `build` | give a `param` a value for this command |
| `--jobs k` | `run --study` | the most cases that run at the same time. Default: `min(H, number of cases)` for `H` hardware threads. A case starts only when its memory estimate fits in the budget beside the running cases, and the table gives the most cases that ran at a time (spec §20.5). |
| `-o dir` | `run`, `build` | the output directory. Default `.skein/`. |

Current `run` writes VTU results under `.skein/results/<model>/`.
The solver determines each file name.
A run with `--set` writes them under `.skein/results/<model>/set/<words>/`,
so that it does not write over the results of the default values. A
submodel (spec §15.1) reads the state of the other model from this
directory, and runs the other model first if the state is not there.
`run --study` writes the results of each case under
`.skein/studies/<model>/<study>/<case>/`, and the table of all cases to
`.skein/studies/<model>/<study>.toon` (§6.5).
Case 1 writes the stages that no varied parameter reaches (spec §20.4) to
`.skein/studies/<model>/<study>/shared/`, and the other cases read them:
the mesh, and also the solution of each step when the varied parameters
reach only `expect` and `report` expressions. The shell and flow
generators share no stages.
For a flow model, `run` makes the mesh itself and writes it to
`flow.xy` and `flow.tri` in the results directory of the model or the
case. The generated program reads the mesh from these files.
Each case gets its threads in the environment variable `SKEIN_THREADS`.
The program of a model reads it for its assembly and its matrix products. Without
it, the program uses all hardware threads.
Each case also gets the memory budget of the study in bytes in
`SKEIN_MEMORY_BUDGET`, in place of the budget that the program reads
itself (spec §20.5).

---

## 6. What skein prints

The default output is [TOON](https://github.com/toon-format/toon). Reason:
TOON uses fewer tokens than JSON for tables of uniform records, and an LLM
reads it directly. `--format json` gives the same records as JSON.

### 6.1 Diagnostics

```toon
diagnostics[2]{class,id,file,line,col,msg,fix}:
  error,unit_mismatch,beam.skn,31,30,"Pressure + Density in eq of physics solid","rho*g has Force/Volume; did you mean rho*(0, 0, -g)?"
  lint,unused_material,beam.skn,18,10,"material aluminium is not assigned","assign it in a physics block or remove it"
summary: errors=1 lints=1
```

- `line` and `col` start at 1. `col` counts Unicode code points.
- `fix` is a suggested change in words or in Skein text. It is never empty.
- The list keeps the order of the file.

### 6.2 The check report

`check` prints the default natural conditions after the diagnostics:

```toon
free_boundary[1]{physics,field,region,condition}:
  solid,u,"boundary(beam) - root - tip","sigma*n = 0"
```

Then it prints the method that it chose for each physics with no `method`
item (spec §13). A physics with a `method` item has no row:

```toon
method[1]{physics,kind,mesh,order,stabilize}:
  solid,cg,m,2,none
```

### 6.3 Verification

```toon
verify[6]{stage,id,result,value,limit}:
  pre,rigid_body,pass,0,0
  pre,mesh_quality,pass,"12.4 deg","5 deg"
  post,residual,pass,3.1e-9,1e-7
  post,equilibrium,pass,"2.0e-6 N","2.0e-3 N"
  post,energy_balance,pass,4.2e-10,1e-6
  post,expect,pass,"kirsch: 301.2 MPa","300 MPa within 5%"
```

### 6.4 Reports

Each `report` line in a step gives one record:

```toon
report[1]{step,expr,value}:
  load,"reaction(root)","(-0.0, 0.0, 2000.0) N"
```

### 6.5 Studies

`run --study` gives one record for each case (spec §20.6):

```toon
study: re_sweep
jobs: 3
cases[3]{case,Re,status,verify,cd,time_s}:
  1,10,ok,pass,2.81,4.2
  2,20,ok,pass,2.05,4.6
  3,40,ok,pass,1.54,5.1
```

A case that fails gets the status `failed` and the diagnostic
[`case_failed`](checks.md). The other cases continue.

A `minimize` study also gives its best case and the reason that its loop
stopped (`xtol`, `max_cases` or `case n failed`). Only the verification
records of the best case set the exit status. With `method: gradient`,
the records of each case also hold the derivative of the objective by
each param, in the unit of the objective for each unit of the param:

```toon
study: best
jobs: 1
best: 10
stop: xtol
cases[10]{case,stiff,status,verify,cost,optimum,time_s}:
  1,50.0 MPa,ok,pass,52.000000000000085 N,"",0.01
  ...
  10,10.001682370269508 MPa,ok,pass,20.000000282989255 N,pass,0.01
```

```toon
gradient[1]{param,value}:
  stiff,0.0003363891618904935
```

### 6.6 Memory

`run` ends its records with the memory of the run (spec §20.5): the
number of memory checks (one before the assembly of each solve), the
largest estimate, and the measured peak, in MB. When a check finds that
the solve fits only on 1 thread, the run goes on with 1 thread and
writes one line on stderr.

```toon
memory_checks: 1
estimate_mb: 1688
peak_mb: 1611
```

---

## 7. Building

Current `skein run` writes `.skein/build/<model>/main.goose` and calls
Goose to build a native executable. A run with `--set` builds in
`.skein/build/<model>/set/<words>/`, so that runs with other values do
not share the program. Set `GOOSE` to select the compiler.
`./build.sh` builds the Skein compiler itself.

```
skein build samples/01_cantilever.skn --target wasm
```

`build --target wasm` writes one HTML page,
`.skein/build/<model>/web/<model>.html`, that meshes and solves the model
in the browser (implementation.md §8). The page has a control for each
`param` with a range of literals, for example `param L: Length = 200 mm
in 100 mm..400 mm`. When a control changes, the page meshes the model
again, solves it with the new values, and shows the records of the run
in TOON, the `mesh_quality` record of a `tri` or `tet` mesh and the deformed mesh. The run in the
page gives the same records and VTU file as `skein run --set` with the
same values. The page also shows the records of a native run at the
default values.

It needs Goose, a C compiler for the native runs (`cc`, or set `CC`) and
`emcc` on `PATH`. It runs the model natively twice: at the default values
and at the corner of the param ranges with the most dofs. A native run
that fails stops `build` with an error. The peak memory of these runs
sets the memory limit of the page: the author sets nothing. The memory of
a run starts at 64 MB and grows to that limit. If the browser cannot give
the memory, the run stops and the page shows an error with the limit. A
model that needs more than 15 GB gets a page with the native results only,
and the page gives the size of the model.

Current limits:

- only the models of the simplex generator (`tri` or `tet` cells) and
  of the hex generator (a box, a sector or a union of boxes, `hex`
  cells). Sample 01 (the cantilever) and the Hisada interposer model
  (`examples/package/interposer_warpage.skn`, `hex27` cells) are the
  validated models of the page;
- a value must not change the structure of the program, for example the
  number of steps or of shapes. `build` stops with an error for such a
  `param`;
- the page runs on one thread, in a browser with WebAssembly memory64
  (tested in Chromium 156), from a `file://` or an `http://` address;
- no `use` files.

The page has a function `skeinSolve(values)` for other pages: `values`
holds a number for each param, in the unit of its default value. It gives
a promise of `{ code, out, err, files, seconds, sets }`: the exit status
(section 10), the records in TOON, standard error, the VTU text of each
result file, the time and the words of `--set`.

The planned `build --target native` writes the executable of `run`.

### 7.1 Export

```
skein export samples/01_cantilever.skn --geometry beam --format stl
```

`export` writes one `let` shape as an STL file, or as a tessellated STEP AP242 file with `--format step` (spec §22). It needs
no solve. It checks the model, meshes the `let` shape that `--geometry`
names, and writes the file to standard output:

```
skein export samples/01_cantilever.skn --geometry beam > beam.stl
```

The STL is ASCII. A `Solid` gives the closed surface of its tetrahedra,
with each normal out of the solid. A 2D `Surface` gives its triangles in
the plane z = 0. The coordinates are in the length unit of the header (mm
for the mm systems, m for SI). An error of the model goes to standard
error.

| Option | Meaning |
|---|---|
| `--geometry <name>` | the geometry to write. Required. |
| `--format stl\|step` | `stl` is the default. `step` writes STEP AP242 with tessellated geometry: one `TESSELLATED_SHELL` with the same triangles as the STL, in the length unit of the header. |
| `--size length` | the size of the triangles, for example `--size 2mm`. Default: the `size` of the first `mesh` block on the geometry, else of the first `mesh` block. |

`export` writes geometry only. Skein does not import a mesh or a CAD file
(decision 8).

---

## 8. The language server

The language server is planned. The current executable rejects `lsp`.

`skein lsp --stdio` serves an editor. It gives:

- the diagnostics of `check` and `lint` after each change;
- the type and the dimension of the name under the cursor, for example
  `sigma: Tensor2<Pressure>   [M L^-1 T^-2]`;
- completion of names, units, settings and assumptions;
- `fmt` on save;
- go to definition, across a `results` file.

---

## 9. Files

| Path | Contents |
|---|---|
| `*.skn` | Skein source |
| `.skein/build/` | generated Goose, C and native executables. `.skein/build/<model>/web/` holds the page of `build --target wasm`. |
| `.skein/results/<model>/` | VTU result files. A single run prints its TOON log to standard output. |
| `.skein/results/<model>/<step>.state_sol` | the state of a step, for submodels |
| `.skein/studies/<model>/<study>/` | one directory of results for each case |
| `.skein/studies/<model>/<study>.toon` | the table of the study |
| `.skein/report/<model>/` | the TOON of `check`, `verify` and `run`, and the small VTU files of a report (§12) |

---

## 10. Exit status

| Status | Meaning |
|---|---|
| 0 | success |
| 1 | an error, a failed verification check, or a failed `expect` |
| 2 | wrong use of the command line |
| 3 | an internal fault of `skein`. Please report it. |
| 4 | in the page: the model needs more memory than the page gives it |

---

## 11. The Markdown document

```
skein md samples/01_cantilever.skn > cantilever.md
```

`md` checks the model, then writes one Markdown document to standard
output. A person reads it as the write-up of the model. The math is TeX
inside `$...$` and `$$...$$`, which GitHub and the common Markdown readers
render. The tables are Markdown tables. The document contains:

- the comment at the top of the file, as the introduction;
- the space, the unit system and the `scales` line;
- tables of the `param` and `let` values and of each material, with their
  types and units;
- the geometry and the regions in words;
- the mesh and method settings;
- each physics: its unknowns, its material, its named assumptions, and its
  `let` and `eq` statements as equations in strong form. A physics that
  uses a template also shows the statements of the template;
- the boundary conditions as equations;
- the steps in order;
- the checks: each `expect` with its expected value and its limit.

The items are in the order of the file, with the names of the file, and
each item gives its line, so that a reader can find each part in the
source. Reason: the document is for a review of the model, and the
reviewer changes the source.

The document contains no results and no date. The same model always gives
the same document.

The math has no backslash before a punctuation character, so it has no
`\,` and no `\_`. Reason: GitHub applies the Markdown escapes inside
`$...$`, and the other readers do not. A name such as `h_air` is h with
the subscript "air".

A model with errors gives no document: `md` writes the diagnostics of
`check` to standard error and exits 1.

---

## 12. The HTML report

```
skein report samples/01_cantilever.skn --study sizes -o cantilever.html
```

`report` runs `check`, `verify` and `run` on the model, then writes one
HTML file. The file opens offline in a browser with WebGL2: it holds the
viewer, the source, the tables and the results. Without `-o`, the file is
`<model>.html`. Without `--study`, the report uses the first study of the
model, or a single run if the model has no study. With `--no-run`,
`report` uses the last run of the study in `.skein/studies/` and solves
nothing. Reason: the interposer sweep of 72 cases runs for more than an
hour. `skein run --study` writes a hash of the source of the model and of
its `use` files in `.skein/studies/<model>/<study>.hash`. If the hash is
missing or is not the hash of the current source, `report --no-run` stops
with exit status 2 and tells you to run `report` without `--no-run`.

The page has two panes:

- the source of the model. A region name marks its region in the viewer.
  A param name moves to its study slider. The chips at the end of an
  `output`, `expect`, `report`, `step` or `material` line show its field,
  point, step or region.
- the viewer of the results: the contours of a field, the deformed shape
  with a scale, a cut plane, bands, the mesh lines and a probe. A click
  on the mesh gives the value at that point and the regions there. Each
  step of a case is in the step list.

Under the viewer are the solver tables: the method of each physics, the
linear solver and its preconditioner, the mesh and its quality, the
solve of each step (cells, dofs, iterations, time, convergence), the
verification checks with their values, limits and results, the reports,
and the study table. Each row links to its line of the source.

With a study, each param of the grid has a slider. A slider moves
between the stored cases. The page does not solve: it shows only the
cases of the native run.

The input "Open VTU" opens a VTU file of another solver in the same
viewer, without source links.

To keep the file small, each VTU result keeps only the corner nodes of
each cell. The mesh is Float32, and cases with the same mesh share it.
The fields are Int16, scaled to the largest value of each field: the
viewer shows about 4 significant digits, and the tables give the exact
values of the run. gzip compresses each part. Reason: the VTU files of
sample 01 and its 16 cases are 230 MB as text. `report` prints the path,
the size in bytes and the number of cases:

```
report: cantilever.html
bytes: 9123456
cases: 16
```

The exit status is that of `run` (§10).

### 12.1 The look of the page

The page has the look of the site teoyujie.org (its rule 6, `REQUIREMENTS.md`
in the repository yujieteo/site), light only: a near-white background,
16 px Inter for the text, IBM Plex Mono for the source, the tables and the
legend, quiet underlined links. There are no themes and no dark mode.
Reason: the captain asked on 2026-10-31 that the HTML reports obey the same
principles as his website.

`src/report/kowalski.css` is the file `public/shared/kowalski.css` of the site
as it is. `report` puts it before `viewer.css` in the page. It has the tokens,
the `.button` look of Sonner, the press (`scale(0.97)`, 150 ms), the ease-out
curve, the focus ring, the disabled look, the tap area of 44 px, the touch
baseline and the reduced-motion rule. Each rule has its source link in that
file (the skills repository of Emil Kowalski, Sonner, animations.dev,
aiforui.dev). Hover styles apply only under `(hover: hover) and (pointer:
fine)`. Reason: the captain asked on 2026-10-31 that the controls be exactly
those of Emil Kowalski and of the site. Every button of the report (`auto`,
`fit`, `Open VTU`, the chips of the source) is a `.button`. The file input is
hidden and its label is a `.button`. The selects, the number input, the range
slider and the checkboxes take the rules of the site tool
`public/play/goosestrike/style.css` (kowalski.css has no rule for them): the
`--bg` fill, the `--line` border, the 44 px height, 16 px text, `accent-color:
var(--ink)` and the outline of 2 px in `--accent`. The buttons of the controls
have the 44 px height of the site tool buttons.

The code colours are those of the code blocks of emilkowal.ski, which use the
Shiki theme `github-light` (read from the computed styles of
https://emilkowal.ski/ui/building-a-toast-component): keywords `#d73a49`,
numbers `#005cc5`, strings `#032f62`, types `#6f42c1`, units `#e36209`,
comments `#6a737d`. Pass is `#22863a` and fail is `#d73a49` of the same theme.
`report` takes the token classes from the lexer of the compiler
(`token_classes` in `report.goose`), not from a regex. Reason: the captain
asked that the colour coding of Skein code match Kowalski's recommendations.
He publishes no other recommendation: the colours above are those of his site.

The terminal has no colour. Reason: Emil Kowalski publishes no guidance on
terminal colours (checked on his site, his skills repository, aiforui.dev,
Sonner), and the CLI prints only TOON, which stays plain for the tools that
read it.

The report is one template for every model: no model name is in `page.html`,
`viewer.css`, `viewer.js` or `report.goose`.

The two panes keep their layout, and the column of 692 px does not apply.
Reason: the source and the viewer are side by side, and one column of 692 px
has no room for both. The flash that marks the target of a jump lasts 1.6 s,
not under 300 ms. Reason: the eye must find the line after the scroll.
The colours of the contour bar are data, so they stay.

The two fonts are in `src/report/fonts/` with their SIL OFL 1.1 licence
files. `report` writes them into the file as woff2 in base64, so the file
needs no network. Reason: the report opens offline. Together they cost
87 628 bytes of woff2 and 116 840 bytes in the file.
