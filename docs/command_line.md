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

</details>

`skein` is one tool for all the work on a model. Each command does the
checks of the commands before it: `lint` does `check` first, and `run` does
`check`, the pre-solve verification, the solve and the post-solve
verification.

---

## 1. Commands

| Command | What it does |
|---|---|
| `skein check model.skn` | parses and type-checks. It writes nothing and runs nothing. |
| `skein fmt model.skn…` | rewrites the files in the standard layout (§3) |
| `skein lint model.skn` | `check`, then the lint rules |
| `skein verify model.skn` | `check`, the mesh, and the pre-solve verification. No solve. |
| `skein verify convergence model.skn` | the convergence study (checks.md §4.3) |
| `skein build model.skn` | writes the Goose source and builds a native executable or a WebAssembly module |
| `skein run model.skn` | builds and runs the model, with all verification |
| `skein export model.skn` | writes the geometry of a part to STL or STEP (§7.1) |
| `skein view model.skn` | opens the results in the WebGL viewer |
| `skein lsp` | serves an editor over the Language Server Protocol |
| `skein help [command]` | prints help |
| `skein help --llm` | prints the language card for an LLM: the grammar, the keywords, the standard library and one example, in about 2000 tokens |
| `skein help std` | prints the Skein source of the standard library (spec §4.2) |

Global options:

| Option | Meaning |
|---|---|
| `--format toon\|json` | the output format (§6). Default `toon`. |

---

## 2. Checking a model

```
skein check samples/01_cantilever.skn
```

`check` applies the rules of the [specification](skein_spec.md#19-static-semantics-checklist)
in order. It continues after an error where it can, so that one run reports
as many errors as possible. It also prints the parts of the boundary that
get the default natural condition (spec §12.3), because a missing `bc` is
the most common modelling fault.

| Option | Meaning |
|---|---|
| `--groups` | print the dimensionless groups of the model (spec §21.3) |
| `--study <name>` | check each case of the study, and print the stages that each parameter reaches (spec §20.4). No solve. |
| `--set name=value` | give a `param` a value for this command (spec §20.1). The value has a unit and must be in the range of the `param`. |

`check` needs no mesh, so it is fast. Its target is under 50 ms for each of
the samples. Reason: an LLM runs `check` after each edit.

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
   has no blank lines in it. In a run of `let` lines, `fmt` keeps one blank
   line that the author wrote between two lines. `fmt` collapses two or
   more blank lines to one.
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

```
skein run samples/02_cavity.skn
skein run samples/03_birth_death.skn --from lift_3
skein verify convergence samples/05_meshless_plate.skn --levels 3
skein run samples/06_cylinder_sweep.skn --study re_sweep --jobs 3
skein verify samples/01_cantilever.skn --dfm
```

| Option | Commands | Meaning |
|---|---|---|
| `--from <step>` | `run` | restart from the state at the end of the step before `<step>` (spec §14.7) |
| `--target native\|wasm` | `run`, `build` | native executable (default) or WebAssembly |
| `--threads n` | `run` | worker threads for assembly and solvers. Default: all hardware threads. With `--study`, the default is `max(1, H / k)` for `H` hardware threads and `k` jobs (spec §20.5). |
| `--levels n` | `verify convergence` | number of meshes. Default 3. |
| `--mesh` | `lint`, `verify` | make the mesh and run the checks that need it |
| `--dfm` | `verify` | run the DFM checks of the `manufacture` items ([checks.md §4.2a](checks.md#42a-design-for-manufacture)) |
| `--study <name>` | `run` | run each case of the study (spec §20) |
| `--set name=value` | `run`, `verify`, `build` | give a `param` a value for this command |
| `--jobs k` | `run --study` | cases that run at the same time. Default: `min(H, number of cases)` for `H` hardware threads (spec §20.5). |
| `-o dir` | `run`, `build` | the output directory. Default `.skein/`. |

`run` writes the results to `.skein/results/<model>/<step>.vtu` (one file
per output time) and the state to `.skein/state/<model>/<step>.sks`.
`run --study` writes the results of each case under
`.skein/studies/<model>/<study>/<case>/`, and the table of all cases to
`.skein/studies/<model>/<study>.toon` (§6.5).

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
cases[3]{case,Re,status,verify,cd,time_s}:
  1,10,ok,pass,2.81,4.2
  2,20,ok,pass,2.05,4.6
  3,40,ok,pass,1.54,5.1
```

A case that fails gets the status `failed` and the diagnostic
[`case_failed`](checks.md). The other cases continue.

---

## 7. Building

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

`export` writes one `geometry` value as a part file (spec §22.5). It needs
no solve.

| Option | Meaning |
|---|---|
| `--geometry <name>` | the geometry to write. Required. |
| `--format stl\|step` | STL is a triangle mesh (phase 2). STEP is tessellated in phase 4 and exact in phase 5. |
| `--size length` | the size of the triangles. Default: the mesh size of the model. |

`export` writes geometry only. Skein does not import a mesh or a CAD file
(decision 8).

---

## 8. The language server

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
| `.skein/build/` | generated Goose and C |
| `.skein/results/<model>/` | VTU result files and the TOON run log |
| `.skein/state/<model>/<step>.sks` | the state at the end of each step, for restarts and submodels |
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
