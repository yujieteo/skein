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
| `skein verify convergence model.skn` | the convergence study (checks.md §4.3) | Planned |
| `skein build model.skn` | writes the Goose source and builds a native executable or a WebAssembly module | Planned |
| `skein run model.skn` | builds and runs the model, with the checks of its solver | Current: limited models |
| `skein export model.skn` | writes the geometry of a part to STL or STEP (§7.1) | Current: STL only |
| `skein view model.skn` | opens the results in the WebGL viewer | Planned |
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
It does not support `convergence` or `--mesh`.

The examples and option table below describe the planned interface:

```
skein run samples/02_cavity.skn
skein verify convergence samples/05_meshless_plate.skn --levels 3
skein run samples/06_cylinder_sweep.skn --study re_sweep --jobs 3
```

| Option | Commands | Meaning |
|---|---|---|
| `--target native\|wasm` | `run`, `build` | native executable (default) or WebAssembly |
| `--threads n` | `run` | worker threads for assembly and solvers. Default: all hardware threads. With `--study`, the default is `max(1, H / k)` for `H` hardware threads and `k` jobs (spec §20.5). |
| `--levels n` | `verify convergence` | number of meshes. Default 3. |
| `--mesh` | `lint`, `verify` | make the mesh and run the checks that need it |
| `--study <name>` | `run` | run each case of the study (spec §20) |
| `--set name=value` | `run`, `verify`, `build` | give a `param` a value for this command |
| `--jobs k` | `run --study` | cases that run at the same time. Default: `min(H, number of cases)` for `H` hardware threads. The runner lowers `k` when `k` cases do not fit in the memory, and the table gives the `k` that it used (spec §20.5). |
| `-o dir` | `run`, `build` | the output directory. Default `.skein/`. |

Current `run` writes VTU results under `.skein/results/<model>/`.
The solver determines each file name.
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
The program of a model reads it for its assembly and its solves. Without
it, the program uses all hardware threads.

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

---

## 7. Building

The `skein build` command and its options below are planned.
Current `skein run` writes `.skein/build/<model>/main.goose` and calls
Goose to build a native executable. Set `GOOSE` to select the compiler.
`./build.sh` builds the Skein compiler itself.

```
skein build samples/01_cantilever.skn --target wasm
```

`build` writes `.skein/build/<model>.goose`, then calls the Goose compiler.
For `native`, Goose writes C and builds it with the C compiler that it
finds. For `wasm`, Goose writes C and `skein` compiles it with `emcc`, in
the same way as earlier Goose projects that run in the browser.

| Option | Meaning |
|---|---|
| `--emit goose` | write the Goose source and stop |
| `--emit c` | write the C and stop |
| `--goose path` | the Goose compiler to use |

The WebAssembly target compiles the C from Goose with emcc (spec §18.1).

### 7.1 Export

```
skein export samples/01_cantilever.skn --geometry beam --format stl
```

`export` writes one `let` shape as an STL file (spec §22). It needs
no solve. It checks the model, meshes the `let` shape that `--geometry`
names, and writes the STL to standard output:

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
| `--format stl` | current format, and the default. STEP remains planned work. |
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
| `.skein/build/` | generated Goose, C and native executables |
| `.skein/results/<model>/` | VTU result files. A single run prints its TOON log to standard output. |
| `.skein/state/<model>/<step>.sks` | planned state files for submodels |
| `.skein/studies/<model>/<study>/` | one directory of results for each case |
| `.skein/studies/<model>/<study>.toon` | the table of the study |

---

## 10. Exit status

| Status | Meaning |
|---|---|
| 0 | success |
| 1 | an error, a failed verification check, or a failed `expect` |
| 2 | wrong use of the command line |
| 3 | an internal fault of `skein`. Please report it. |

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
