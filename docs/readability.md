# Readability: Skein and femahoy

This page compares one model in Skein and in the two femahoy inputs, with
numbers. The model is the femahoy case `validation/cases/grid-cantilever`:
a steel cantilever 100 × 10 × 10 mm, fixed at x = 0, with a shear force of
100 N on its free end. The page also gives the results of the end-to-end
LLM workflow (§4).

**Status:** complete. §4.1 gives 3 runs for each path of the LLM workflow,
with `skein check` and `skein run` for Skein.

---

## 1. The model in each format

The three inputs give the same model: the same beam, mesh size, material,
support and load. Each input reports the tip deflection and the root
reaction.

- **Skein:** one file.
- **femahoy CalculiX deck:** `beam.geo` (gmsh) makes the mesh, and `job.inp`
  holds the model. The deck takes the mesh with `*INCLUDE`. This is the
  layout of the femahoy cases that use `ccx`, for example
  `frequency-cantilever`.
- **femahoy grid input:** `grid.toml`, the input of this femahoy case. The
  built-in grid solver of femahoy solves it.

<table>
<tr><th>Skein: cantilever.skn</th><th>femahoy CalculiX: beam.geo and job.inp</th></tr>
<tr><td valign="top"><pre>
// A steel cantilever 100 x 10 x 10 mm with a shear force of 100 N on its free end.

model cantilever: 3d, mm_N_MPa_s

let L: Length = 100 mm
let b: Length = 10 mm
let H: Length = 10 mm
let P: Force = 100 N

scales H, P/(b*H)

let beam: Solid = box((0, 0, 0) mm, (L, b, H))
let root: Surface = boundary(beam, x == 0 mm)
let tip: Surface = boundary(beam, x == L)

material steel {
    E: Pressure = 210 GPa
    nu: Real = 0.3
}

mesh m on beam {
    size: 2.5 mm
    cell: hex
}

physics solid on beam: linear_elastic solves u: Vector&lt;Length&gt; {
    material steel
}

bc root: u = 0

step load: static {
    bc tip: sigma*n = (0, 0, -P/(b*H))
    output u, sigma
    report tip_uz = u.z at (L, b/2, H/2)
    report reaction(root)
}
</pre></td><td valign="top"><pre>
// Steel cantilever 100 x 10 x 10 mm: a mapped grid of 40 x 4 x 4 cells.
L = 100; B = 10; H = 10;
Point(1) = {0, 0, 0}; Point(2) = {L, 0, 0}; Point(3) = {L, B, 0}; Point(4) = {0, B, 0};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4}; Plane Surface(1) = {1};
Transfinite Curve{1, 3} = 41; Transfinite Curve{2, 4} = 5;
Transfinite Surface{1}; Recombine Surface{1};
out[] = Extrude {0, 0, H} { Surface{1}; Layers{4}; Recombine; };
Physical Volume("BEAM") = {out[1]};
Physical Surface("ROOT") = {out[5]};
Physical Surface("TIP") = {out[3]};
Mesh.ElementOrder = 2; Mesh.SecondOrderIncomplete = 1;
Mesh.SaveGroupsOfNodes = 1;
</pre><pre>
** Steel cantilever 100 x 10 x 10 mm, 100 N shear on the free end. Units: SI-mm.
*INCLUDE, INPUT=mesh.inp
*NODE, NSET=LOADPT
900001, 100.0, 5.0, 5.0
*RIGID BODY, NSET=TIP, REF NODE=900001
*MATERIAL, NAME=STEEL
*ELASTIC
210000.0, 0.3
*SOLID SECTION, ELSET=BEAM, MATERIAL=STEEL
*BOUNDARY
ROOT, 1, 3
*STEP, NAME=LOAD
*STATIC
*CLOAD
900001, 3, -100.0
*NODE PRINT, NSET=LOADPT
U
*NODE PRINT, NSET=ROOT, TOTALS=ONLY
RF
*NODE FILE
U
*EL FILE
S
*END STEP
</pre></td></tr>
</table>

<details>
<summary>femahoy grid input: grid.toml</summary>

```toml
# A steel beam 100 x 10 x 10 mm, clamped at x = 0, with a shear force of 100 N on its free end.
units = "SI-mm"

[grid]
max_cell = 2.5

[material.steel]
youngs_modulus = 210e3
poisson_ratio = 0.3

[part.beam]
material = "steel"
x = [0.0, 100.0]
y = [0.0, 10.0]
z = [0.0, 10.0]

[face.root]
side = "x-"

[face.tip]
side = "x+"

[[load]]
ux = 0.0
uy = 0.0
uz = 0.0
face = "root"

[[load]]
traction = [0.0, 0.0, -1.0]
face = "tip"

[step]
type = "static"

[probe.tip]
point = [100.0, 5.0, 5.0]
```

</details>

### 1.1 Size

| Input | Files | Lines | Lines without comments and blank lines | Bytes |
|---|---|---|---|---|
| Skein | 1 | 37 | 27 | 730 |
| femahoy CalculiX | 2 (+1 made) | 37 (13 + 24) | 35 (12 + 23) | 1121 (674 + 447) |
| femahoy grid | 1 | 37 | 27 | 491 |

The CalculiX input also needs the mesh file. These commands make it:

```sh
gmsh beam.geo -3 -format inp -o mesh.inp
# ccx refuses the 2D faces that gmsh writes for the physical surfaces: remove them.
awk '/^\*ELEMENT, type=CPS8/{s=1;next} /^\*ELSET,ELSET=(ROOT|TIP)/{s=1;next} /^\*/{s=0} !s' mesh.inp > m && mv m mesh.inp
```

`mesh.inp` has 5399 lines: 3665 nodes and 640 C3D20 elements. A person
does not write it, but a deck that an LLM writes in one file must hold it
(§4). The femahoy cases do the same step with a `mesh.py` of 51 lines.

### 1.2 What each input states, and what it leaves out

| Fact | Skein | femahoy CalculiX | femahoy grid |
|---|---|---|---|
| units of each value | in the value (`210 GPa`) | a comment (`Units: SI-mm`) | one line (`units = "SI-mm"`) |
| the load as a force | `P: Force = 100 N` | `-100.0` on node 900001, DOF 3 | a traction of 1 MPa: a person divides 100 N by the face area |
| the tip face | `boundary(beam, x == L)` | a physical surface of `beam.geo`, by surface number `out[3]` | `side = "x+"` |
| how the load reaches the tip face | uniform traction | a rigid tip section (`*RIGID BODY`) | uniform traction |
| the element | default: order 2 on hex cells | `Mesh.ElementOrder = 2` in `beam.geo`, C3D20 | the grid solver chooses |
| the dimensions of the model | type-checked by `skein check` | not checked | not checked |

### 1.3 The results

The Skein file ran with `skein run` on main. The two femahoy inputs ran
with `fh run` (femahoy and CalculiX on this computer). All runs are from
2026-10-30.

| Input | Mesh and element | Tip deflection u_z | Root reaction F_z |
|---|---|---|---|
| Skein | 40 × 4 × 4 hex27 (order 2), 19683 unknowns, 542 CG iterations | −0.190502 mm | 100.000 N |
| femahoy grid | grid solver, `max_cell` 2.5 mm | −0.190113 mm | 100.000 N |
| femahoy CalculiX | 40 × 4 × 4 C3D20, rigid tip | −0.190408 mm | 100.000 N |
| Euler-Bernoulli, PL³/(3EI) | | −0.190476 mm | 100 N |
| Timoshenko, the reference of the femahoy case | | −0.191933 mm | 100 N |

The Skein tip deflection is 0.75 % below the Timoshenko reference, 0.20 %
above the femahoy grid solver and 0.01 % above Euler-Bernoulli. The
equilibrium check of `skein run` passed: the reaction equals the applied
load to 2.4e-8 N.

---

## 2. Three edit tasks

Each task starts from the inputs of §1. Each answer passed its checker:
`skein check` for Skein, and `fh run` for the femahoy inputs. The diff
counts come from `git diff --no-index --numstat` over the files that a
person writes. The made mesh files are not in the counts.

| Task | Request |
|---|---|
| E1 | change the tip load from 100 N to 150 N |
| E2 | add a second static step in which the tip load is 200 N |
| E3 | sweep the beam height H over 10, 15 and 20 mm |

### 2.1 Diff size

Lines added and lines removed. "Files" is the number of files that the
change adds or changes.

| Task | Skein | femahoy CalculiX | femahoy grid |
|---|---|---|---|
| E1 | +1 −1, 1 file | +1 −1, 1 file | +1 −1, 1 file |
| E2 | +7 −0, 1 file | +13 −0, 1 file | +37 −0, 1 new file |
| E3 | +5 −1, 1 file | +74 −0, 4 new files, and 2 more mesh runs | +74 −0, 2 new files |

For E3, each new femahoy file is a copy of the base file. The copies differ
from the base in these lines:

| Input | Lines that differ in each copy | What the person computes by hand |
|---|---|---|
| femahoy CalculiX | 2 in `beam.geo`, 3 in `job.inp` | the point of the load node (H/2) |
| femahoy grid | 4 in `grid.toml` | the traction 100 N / (10 mm × H), and the probe point (H/2) |
| Skein | none: no copy | none: the traction `-P/(b*H)` and the report point `H/2` follow `H` |

Femahoy has no parameter sweep. A grid input has one `[step]`, so E2 in the
grid input is a second case. CalculiX and gmsh have no parameters.

### 2.2 The diffs

E1, Skein and CalculiX:

```diff
-let P: Force = 100 N
+let P: Force = 150 N
```

```diff
-900001, 3, -100.0
+900001, 3, -150.0
```

E2, Skein:

```diff
+
+step load2: static {
+    bc tip: sigma*n = (0, 0, -2*P/(b*H))
+    output u, sigma
+    report tip_uz = u.z at (L, b/2, H/2)
+    report reaction(root)
+}
```

E2, CalculiX:

```diff
+*STEP, NAME=LOAD2
+*STATIC
+*CLOAD
+900001, 3, -200.0
+*NODE PRINT, NSET=LOADPT
+U
+*NODE PRINT, NSET=ROOT, TOTALS=ONLY
+RF
+*NODE FILE
+U
+*EL FILE
+S
+*END STEP
```

E3, Skein:

```diff
-let H: Length = 10 mm
+param H: Length = 10 mm in 5 mm..50 mm
@@
+
+study height: grid {
+    H: [10, 15, 20] mm
+}
```

E3, CalculiX: `beam_h15.geo` and `job_h15.inp` against the base files. The
files for 20 mm are the same with 20 and 10.0.

```diff
-// Steel cantilever 100 x 10 x 10 mm: a mapped grid of 40 x 4 x 4 cells.
-L = 100; B = 10; H = 10;
+// Steel cantilever 100 x 10 x 15 mm: a mapped grid of 40 x 4 x 4 cells.
+L = 100; B = 10; H = 15;
```

```diff
-** Steel cantilever 100 x 10 x 10 mm, 100 N shear on the free end. Units: SI-mm.
-*INCLUDE, INPUT=mesh.inp
+** Steel cantilever 100 x 10 x 15 mm, 100 N shear on the free end. Units: SI-mm.
+*INCLUDE, INPUT=mesh_h15.inp
@@
-900001, 100.0, 5.0, 5.0
+900001, 100.0, 5.0, 7.5
```

E3, femahoy grid: `h15/grid.toml` against the base file.

```diff
-# A steel beam 100 x 10 x 10 mm, clamped at x = 0, with a shear force of 100 N on its free end.
+# A steel beam 100 x 10 x 15 mm, clamped at x = 0, with a shear force of 100 N on its free end.
@@
-z = [0.0, 10.0]
+z = [0.0, 15.0]
@@
-traction = [0.0, 0.0, -1.0]
+traction = [0.0, 0.0, -0.666667]
@@
-point = [100.0, 5.0, 5.0]
+point = [100.0, 5.0, 7.5]
```

### 2.3 The results of the edits

The femahoy answers ran with `fh run`, and the Skein answers ran with
`skein run` on main. The results show that each edit gives the model that
the task asks for.

| Task | Skein u_z, F_z | femahoy CalculiX u_z, F_z | femahoy grid u_z, F_z |
|---|---|---|---|
| E1 | −0.285753 mm, 150.000 N | −0.285612 mm, 150.000 N | −0.285169 mm, 150.000 N |
| E2, step 2 | −0.381004 mm, 200.000 N | −0.380815 mm, 200.000 N | −0.380225 mm, 200.000 N |
| E3, H = 15 mm | not in this slice | −0.056905 mm, 100.000 N | −0.056823 mm, 100.000 N |
| E3, H = 20 mm | not in this slice | −0.024300 mm, 100.000 N | −0.024266 mm, 100.000 N |

`skein run` refuses the `study` of E3 with the error `unsupported`, because
this slice runs one linear static model. The E3 file passes `skein check`.

---

## 3. Repeat the counts

1. Write the base files of §1 and the edited files of §2 into a folder
   outside the repository.
2. Run `git diff --no-index --numstat BASE EDITED` for each task.
3. Run `fh run` on each femahoy answer, and `skein check` on each Skein
   answer.

---

## 4. The LLM workflow

`scripts/llm/workflow.py` gives a fresh model session one request and
records each round:

> a steel cantilever 1 m x 50 mm x 100 mm, fixed at one end, 1 kN tip load, report the tip deflection

| Path | The model gets | Check | Run | Final number |
|---|---|---|---|---|
| `skein` | the card of `skein help --llm`, the request, "Write the complete Skein model." | `skein check` | `skein run` | the first `report` value in a length unit |
| `inp` | the request, and "Write the complete CalculiX input deck", with the mesh in the deck | `fh run` (ccx reads the deck) | the same `fh run` | the largest displacement component in the `.dat` file |

If the check or the run fails, the driver gives the model its output and
"Fix the model." (or "Fix the deck."). For Skein, the run is `skein run`
after a passed `skein check`. For a deck, `fh run` is the check and the run.
The driver stops after the run passes or after round 5. Each run is one new
`claude -p` session with no tools, no settings and no memory. The repair
rounds resume that session. The input tokens count the whole context of each
round, which includes the card.

```sh
./build.sh
scripts/llm/workflow.py skein --out /tmp/llm-runs
scripts/llm/workflow.py inp --out /tmp/llm-runs
```

The output folder gets the answer, the checker output and the `skein run`
output of each round, and
`results-skein.toon` or `results-inp.toon`.
`--request plate` picks the other request of [llm_test.md §3](llm_test.md#3-the-requests).
With `--request`, the results file is `results-<path>-<request>.toon`,
for example `results-skein-plate.toon`.

### 4.1 Results

Euler-Bernoulli for the request: PL³/(3EI) = 0.381 mm, with I = 50 × 100³ / 12 mm⁴.
Model: claude-opus-5-5, 2026-10-30. The round limit is 5.
For the femahoy case of §1, `skein run` gives −0.190502 mm and 100.000 N
on 40 × 4 × 4 hex27 cells (§1.3).

"First-try pass" means that the check passed in round 1.

| Path | Run | First-try pass | Repair rounds | Input tokens | Output tokens | Wall time | Tip deflection | Reaction | Mesh and element |
|---|---|---|---|---|---|---|---|---|---|
| skein | 1 | yes | 1 | 7823 | 1736 | 53.1 s | −0.382161 mm | 1000.000 N | 10 mm hex, order 2 |
| skein | 2 | yes | 0 | 3339 | 1127 | 55.1 s | −0.382161 mm | 1000.000 N | 10 mm hex, order 2 |
| skein | 3 | yes | 0 | 3339 | 998 | 44.3 s | −0.382161 mm | 1000.000 N | 10 mm hex, order 2 |
| inp | 1 | yes | 0 | 656 | 2706 | 26.1 s | −0.379948 mm | 1000 N | 10 B32R beam elements |
| inp | 2 | yes | 0 | 656 | 3364 | 32.6 s | −0.379948 mm | 1000 N | 10 B32R beam elements |
| inp | 3 | yes | 0 | 656 | 6538 | 52.7 s | −0.379456 mm | 1000 N | 40 C3D8I, 126 nodes |

Each Skein solve had 139293 unknowns and took 1223 CG iterations.

An earlier version of the driver gave only the check output back to the
model, not the run output. With it, `skein run` gave a number in 1 of 3
Skein runs (same model, same day):

| Path | Run | First-try pass | Repair rounds | Input tokens | Output tokens | Wall time | Tip deflection |
|---|---|---|---|---|---|---|---|
| skein | 1 | yes | 0 | 3337 | 1055 | 11.5 s | none: `skein run` refused the model |
| skein | 2 | yes | 0 | 3337 | 1125 | 13.1 s | none: `skein run` refused the model |
| skein | 3 | yes | 0 | 3337 | 1070 | 45.9 s | −0.382161 mm |

- The Skein input tokens include the card of `skein help --llm` (about
  2700 tokens).
- The Skein result is 0.3 % above Euler-Bernoulli, and the equilibrium
  check of `skein run` passed. Shear adds about 0.8 % to Euler-Bernoulli.
- `skein run` refuses a report that this slice does not compute, with the
  error `unsupported`. The card lists `min`, `mean` and `max` as report
  functions, and the models used them: `min(u.z, tip)` and `mean(u.z, tip)`
  with the earlier driver, and `max(von_mises(sigma), beam)` in run 1 above.
  In run 1, one repair round removed it.
- Inp runs 1 and 2 used beam elements, not solid elements. Inp run 3 wrote
  its 126 nodes and 40 elements by hand in 234 lines.

The Skein runs above used main with `skein run`, merged with the checker
PR, because `skein check` is not on main yet. The earlier runs and the inp
runs used the same code before the merge to main.
