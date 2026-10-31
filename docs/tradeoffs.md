# Skein tradeoffs

<details>
<summary>Table of contents</summary>

- [0. The goals, in priority order](#0-the-goals-in-priority-order)
- [1. How to read this document](#1-how-to-read-this-document)
- [2. The design against the workloads](#2-the-design-against-the-workloads)
- [3. Static analysis: what Skein proves and what it cannot](#3-static-analysis-what-skein-proves-and-what-it-cannot)
- [4. Parametric parts, text-to-CAD and export](#4-parametric-parts-text-to-cad-and-export)
- [5. Design for manufacture](#5-design-for-manufacture)
- [6. Reliability: the Power of Ten](#6-reliability-the-power-of-ten)
- [7. Parameters and concurrency](#7-parameters-and-concurrency)
- [8. Strict dimensions and dimensionless work](#8-strict-dimensions-and-dimensionless-work)
- [9. LLM-native and AXI-friendly](#9-llm-native-and-axi-friendly)
- [10. The web](#10-the-web)
- [11. Summary of the gaps](#11-summary-of-the-gaps)

</details>

Skein must handle a wide range of engineering simulation work with one
language. Each design choice helps some workloads and costs others. This
document tests the choices against the workloads and records each cost. The
captain answered the questions of [grilling.md](grilling.md) on 2026-10-30,
and this document follows the answers.

---

## 0. The goals, in priority order

Skein has too many goals to build them all at once. This table ranks them.
Simulation is the core. A **protected** goal is one that no other goal may
weaken. A **layered** goal builds on the core, comes later, and can be cut
if it costs too much. The captain set this ranking in
[grilling.md round 0](grilling.md#round-0-the-ranking-of-the-goals).

| Rank | Goal | Class | Syntax cost | Compiler size | Build time | Runtime cost | Phase |
|---|---|---|---|---|---|---|---|
| 1 | Simulation core: strong-form physics, methods, mesher, solvers | protected | the base language | the largest part: lowering, methods, mesher, solvers | the most | none extra, it is the work | 2 to 5 |
| 2 | Dimensions done right: strict units, static analysis and verification by default, with dimensional analysis and parameter sweeps where they count | protected | a type on each declaration, a unit on each literal, a `scales` line, `param` and `study` | the checker (about a third of the front end), the π basis, the scaled SIR, the study runner | small | the post-solve checks cost a few percent of a solve. Scaled solves are equal or faster. | 1 (checks, `param`, `check --groups`), 2 (verification, cases-first studies), 3 (scaled SIR) |
| 3 | LLM- and AXI-native syntax and output, with small human-readable diffs | protected | it shapes all syntax: types only at declarations, TOON, ids with fixes, one item on each line, no alignment | small: the fix strategies, the TOON writer, the language card | none | none | 1 |
| 4 | Concurrency: one work-stealing task pool | layered, high value | none | the task pool in the runtime | small | memory grows with the number of cases in flight | 3 |
| 5 | Web output: the tool chain and the viewer in the browser | layered | none | the browser shell and the viewer | a second target: Goose to C, then `emcc` | the browser limits memory | 4 |
| 6 | Power of Ten reliability | layered, mostly free | none: most rules hold by construction | tests, emitter checks, the kernel splitter | small | none | 2 onwards, as tests |
| 7 | Parametric parts and text-to-CAD | layered | `shell`, `loft`, `pattern` | the STL and tessellated STEP writers. No B-rep kernel. | small | none for simulation | 2 (STL), 4 (tessellated STEP) |
| 8 | Design for manufacture | cut on 2026-10-31 (§5) | none | none | none | none | — |

Rank 2 merges the units, the static analysis and the dimensional analysis
into one goal. Reason, in the words of the captain: dimensional analysis
and units are subtle, and they need "strong typing for real problems
requiring dimensions, but yet … dimensional analysis / parametric sweep
when it counts".

Rules that follow from the ranking:

- A layered goal may add syntax, but it must not change the meaning of a
  model that does not use it.
- When two goals conflict, the higher rank wins. Example: a geometry
  operation that the mesher cannot mesh with good cells is not added for
  text-to-CAD (rank 7 against rank 1).
- The LLM-native goal (rank 3) is protected because it is the reason for
  the language, but it does not overrule correctness: units stay strict
  (rank 2) even when they cost tokens (G0.3).
- If the schedule is late, cut the tessellated STEP export (rank 7) first.
  STL export stays. Rank 6 costs almost nothing, so
  it stays.

---

## 1. How to read this document

Each section gives a choice, what it gives, what it costs, and what the
spec does about the cost. "Gap" means that the spec today cannot express
the workload well. Each gap has a proposed fix and a question in
[grilling.md](grilling.md).

The core choices that each section tests:

| Choice | Spec |
|---|---|
| C1. Physics in strong form, weak form derived | §10, §17 |
| C2. Every type explicit, every quantity with a dimension | §5, §6 |
| C3. Declarative: no loops, no mutable state, functions of one expression only | §1, §7.4 |
| C4. One model per file plus `use` libraries, one namespace | §1, §4 |
| C5. Geometry from distance functions, built-in mesher, no mesh import | §8, §9 |
| C6. Method-neutral SIR, methods chosen per physics | §13, §16 |
| C7. Code generation to Goose, all solvers in Goose | §18 |
| C8. Verification on by default | checks.md |

---

## 2. The design against the workloads

| Workload | What works | Cost or gap | What the spec does |
|---|---|---|---|
| Linear structures (static, modal, buckling) | C1 with `hooke`, `cg`, all procedures, the `linear_elastic` template | none of note | — |
| Periodic structures | `symmetry` with `bloch` (spec §12.6) | no rotational group: a rotor is a full model | accepted. The femahoy cells need only `bloch`. Add a group when a real model needs it. |
| Nonlinear structures (plasticity, hyperelasticity, large deformation) | `plastic`, `diff(W, C)`, `finite_strain` (§10.6) | a user law with many state variables needs a return map that the compiler must derive. A general return map for any yield function is hard. | phase 5. Start with von Mises and Drucker-Prager, and give a clear error for a law the compiler cannot map. |
| Shells and beams | lower-dimension regions (§8.3) | the strong form of a shell (Reissner-Mindlin, Kirchhoff-Love) is long. Users do not write it from memory. | standard-library templates `reissner_mindlin` and `timoshenko`, in Skein source (spec §4.2, G4.1), phase 5 |
| Contact | `contact` (§12.4) | not a PDE: an inequality. It does not fit the strong form, so it is a separate item. | accepted cost |
| Heat transfer | C1 with `fourier`, `transient` | radiation between surfaces needs view factors, a non-local term | a `radiation` item with view factors, phase 5 (G4.2) |
| Incompressible CFD | `supg_pspg`, `vms`, `fv` with Rhie-Chow | none of note for laminar flow | — |
| Turbulent CFD | `turbulent(model)` assumption | RANS models add equations with wall functions. LES needs `vms`. Each model is many lines. | the standard-library template `k_omega_sst` (spec §4.2, G4.1), phase 5 |
| Compressible CFD with shocks | `dg`, `fv` with Roe flux | shocks need limiters or artificial viscosity, and the spec has no setting for them | `limiter` for `fv` (spec §13, G4.3). `dg` has no shock capturing until a model needs it. |
| Multiphysics (thermal-structural) | two physics that share fields (§10.4), `monolithic` or `staggered` | none of note | — |
| FSI | `interface` with ALE (§10.7) | mesh motion and large deformation of the fluid mesh are hard to make robust | phase 5, accepted cost |
| Electromagnetics | `curl` exists | `curl curl` problems need H(curl) (Nédélec) elements. The method table has no element family setting. Time-harmonic problems need complex numbers. | `family: nedelec` in `cg`, phase 5 (G4.4). `Complex<D>` and the `harmonic` procedure (spec §5.1, §14.2, G4.5). |
| Acoustics (Helmholtz) | the equation is easy to write | frequency domain needs complex numbers and absorbing boundaries | `Complex<D>` and `harmonic` (G4.5) |
| Porous media (Darcy, consolidation) | mixed forms with two unknowns | a mixed Darcy form needs H(div) (Raviart-Thomas) elements | `family: raviart_thomas`, phase 5 (G4.4) |
| ODE systems (control, lumped models, chemistry) | global unknowns (§10.5), `ode` method | stiff chemistry needs adaptive implicit time steps with error control. The spec has fixed `dt` only. | a fixed `dt`. Adaptive steps come back when a real model needs them (G4.6). |
| Phase field, reaction-diffusion | C1, `transient` | fourth-order terms (Cahn-Hilliard) need a mixed form or C1 elements | accepted. Write the mixed form with two unknowns. |
| Free surface (level set, VOF) | a level set is a PDE | reinitialization is a process, not an equation | out of scope for v1. No `process` statement until a real model shows which processes are necessary (G4.7). |
| Meshless (large deformation, fracture) | `rkpm`, `gfd` | fragmentation needs particle methods with contact. That is closer to DEM, which is not a PDE. | out of scope for v1 |
| Parameter studies | `param`, `study` (§20) | each case is a whole run | §20.4 reuses the mesh. Each case runs, with no reuse of similar cases (G2.3). |
| Optimization | `study … minimize` (§20.2) | gradient-free methods need many cases. Shape and topology optimization need adjoint gradients. | adjoint gradients from the SIR, phase 5 (G3.3) |
| Uncertainty quantification | `study … sample` with `lhs` or `sobol` | no statistics of the results in the tool | the table is TOON. Other tools compute the statistics (G3.4). |
| General PDEs with no physics | a `nondim` model (§6.5) | no dimension checks | accepted cost, the check report states it |

---

## 3. Static analysis: what Skein proves and what it cannot

The declarative form (C3) and explicit types (C2) let `check` prove these
properties before any run:

| Property | How | Spec |
|---|---|---|
| no memory fault in the model | no pointers, no allocation, no loops in the language | §1 |
| dimensional consistency of each equation, condition and output | dimension algebra | §6 |
| tensor ranks match | type rules | §7.2 |
| each region exists and is inside its field | CSG tree and sampling | §8 |
| each unknown has one physics and one equation | the unknowns are declared in `solves`, and the equations are counted | §10 |
| assumptions do not conflict | assumption table | §11 |
| conditions do not conflict and are on the boundary | condition rules | §12 |
| the method can discretize the physics, with a stable pair | method rules | §13 |
| each study value is in the range of its parameter | ranges | §20 |
| the scales span the dimensions | rank of the dimension matrix | §21 |
| which parts of the model each parameter changes | dependency analysis | §20.4 |

What static analysis cannot prove, and what Skein does instead:

| Not provable before the run | Instead |
|---|---|
| the mesh resolves the solution | `verify convergence` |
| the nonlinear solve converges | the `solver_diverged` run error, `ramp` |
| the time step is stable for a nonlinear problem | `dt_unstable` lint is an estimate. The post-solve energy check finds growth. |
| the model is the right model of the real part | `expect` against tests or benchmarks. This is validation, not verification. |

Cost of C3: a function is one expression. A law that needs a loop or a
local iteration must come from the standard library or the built-ins
(G1.2).

---

## 4. Parametric parts, text-to-CAD and export

Recent text-to-CAD tools have an LLM write a short parametric script, in
the style of CadQuery or build123d, from a sentence. Skein geometry has the
same shape (spec §8, §20): named parameters drive a CSG tree.

What Skein adds to that trend:

- **One set of parameters** drives the part, the simulation and the study.
  A sweep of the wall thickness changes the part, the mesh and the stress
  in one run.
- **Types and units** on each dimension of the part, so an LLM that writes
  `2` for a wall gets an error, not a part 2 m thick.
- **Ranges** on parameters, so a study cannot make a part with a hole wider
  than its wall.

What it costs:

| Cost | Why | Proposal |
|---|---|---|
| distance-function CSG is weaker than a B-rep kernel for some features | sketches with constraints, variable fillets and text are hard on distance functions | `shell(a, t)`, `loft` and `pattern` are distance functions (spec §8.1). No constraint sketches (G1.1, G5.1). |
| no exact STEP export | a CSG tree of distance functions has no exact faces and edges, and a B-rep kernel is a project of its own size | phase 4 exports STEP AP242 with tessellated geometry, which most CAM and inspection tools read. Skein has no exact B-rep (G5.2). |
| no import of CAD | decision 8 | accepted. An LLM can rewrite a simple part as Skein. A complex existing part is out of scope. |

---

## 5. Design for manufacture

The spec cut DFM on 2026-10-31: the `manufacture` and `tolerance` items,
`tol` on a `param`, and the DFM checks. Reason: no sample and no femahoy
case used them, and rank 8 is the first goal to cut.

---

## 6. Reliability: the Power of Ten

The "Power of Ten" rules (Holzmann, 2006, NASA/JPL) are ten rules for
safety-critical C code. The table maps each rule to Skein source and to the
Goose code that the compiler writes and the runtime.

| # | Rule | Skein source | Generated Goose and runtime |
|---|---|---|---|
| 1 | simple control flow, no goto, no recursion | **by construction**: no control flow at all | the emitter writes no `recursive fn`. The runtime has none either: the octree and AMG use loops over levels with a fixed maximum. A test fails if a `recursive` keyword is in the runtime or the output. |
| 2 | a fixed upper bound on each loop | **by construction**: no loops | each loop is over an array whose length is fixed before the solve, or an iteration with a `max_iter` that always has a value (§14.4) |
| 3 | no dynamic memory after initialization | **by construction**: no allocation | Goose has no heap. The runtime sizes each array after the mesh is made and before the solve. The solve does not grow an array. A test runs each sample with a check that no array grows after the start of the solve. |
| 4 | functions of at most 60 lines | **by construction**: a function is one expression (spec §7.4) | runtime functions obey it. The emitter splits each long kernel into functions of at most 60 lines, by blocks of the cell matrix (spec §18.3). A test checks the length of each generated function (G6.4). |
| 5 | at least 2 assertions per function | verification checks are on by default (checks.md) | each runtime function asserts its preconditions. Each kernel asserts its array sizes once at entry. |
| 6 | data in the smallest scope | **the one deviation**: one namespace per model (spec §4), and the standard library and the `use` files enter it with no prefix. Reason: a model is short, a prefix costs tokens on each use, and a clash of two names is a `check` error, so no name is ambiguous. | Goose scopes free data at scope exit |
| 7 | check each return value | no calls that return status | each solver call returns a status that the plan checks. An unchecked status is a compile error in the emitter. |
| 8 | limited preprocessor | **by construction**: none | Goose has no preprocessor |
| 9 | limited pointers | **by construction**: none | Goose has no raw pointers. References are checked by the Goose compiler. |
| 10 | all warnings on, static analysis clean | `skein check` and `skein lint` with no warnings on each sample and benchmark | the Goose compiler with `--debug` on the test suite, and `checkfree fn` on the hot kernels, which proves that no bounds check is left |

Rules 1, 2, 3, 4, 8 and 9 hold for Skein source by construction: the
language cannot express a violation. Rule 6 is the one deviation, for the
one namespace. Rules 5, 7 and 10 are enforced by tests and by the
emitter.

---

## 7. Parameters and concurrency

| Choice | Gives | Costs |
|---|---|---|
| only `param` can change (§20.1) | a reviewer finds each input. The dependency analysis is exact. | a user must mark a value `param` before a study can vary it |
| Go-style workers over typed queues (§20.3) | no shared memory, so no data races. A failed case does not stop the study. The same design works on native threads and Web Workers. | each worker holds its own copy of the shared mesh. Memory grows with the number of workers. §20.5 limits the workers by memory. |
| cases first in phase 2, one work-stealing task pool in phase 3 (§20.5) | cases scale almost perfectly. The pool keeps each thread busy also when the cases are few and large, with no `--jobs` from the user. | each queue needs atomics, or a lock that costs time for each task (G3.2b). |

---

## 8. Strict dimensions and dimensionless work

Strict dimensions (C2) catch the most common modelling fault. But engineers,
mostly in fluids, think in dimensionless groups. The spec keeps both:

| Need | How | Spec |
|---|---|---|
| vary Re, not ρ, U, L, μ | `param Re: Real` and a dimensional `let mu = rho*U*D/Re` | §21.1 |
| the solver in scaled variables | a `scales` line, required in each model with dimensions (G2.2) | §21.2 |
| know the groups of a model | `check --groups`, Buckingham π | §21.3 |
| a model with no physical units at all | a `nondim` model, for example `model m: 2d, nondim` | §6.5 |

Costs:

- A `nondim` model has no dimension checks. It is a deliberate choice in
  the header, and the report says so.
- The π basis is not unique. The tool prefers the `scales` as repeated
  variables and names known forms, but a user can see a group in a
  different form than the one in the textbook.
- A compressible flow needs more scales (a temperature, a speed of sound).
  The `scales_incomplete` error lists the missing dimension.
- Each model has one more line, the `scales` line. An LLM writes it easily,
  and it makes the size of each unknown explicit.
- A study runs each case, also when two cases have equal groups. A reuse of
  similar cases needs every dimensional constant of the model in the test,
  and a missed one gives a wrong result with no error (G2.3).

---

## 9. LLM-native and AXI-friendly

| Choice | Gives | Costs |
|---|---|---|
| types only at declarations | equations stay short and look like a textbook | each `let` has a type, which costs tokens |
| TOON tables in source and output | fewer tokens than JSON for uniform rows | a second syntax inside the language |
| stable ids with a fix in each diagnostic | an LLM repairs a model from the output alone | each new rule needs a fix strategy |
| `allow` only in the file | an LLM cannot turn off a check from the command line to make a run pass | a user must edit the file |
| strong form, with standard-library templates for the common physics | an LLM knows textbook equations well, and a template removes them from the common models | a Skein answer holds the equations, so the model gets a card of about 2700 tokens in each round, and a CalculiX deck needs none ([readability.md §4.1](readability.md#41-results)) |
| one item on each line, no column alignment, a comma after each item of a long list (spec §3.2) | a change to one load, one step or one study parameter is a diff of a few lines that a person can review | `fmt` output is less compact than an aligned table, and rows of a TOON table do not line up |
| `skein help --llm`, a language card of about 2700 tokens | the context of each LLM session is small | one more text to keep in step with the spec |

The input of each Skein round stays larger than for a deck, because of the
card. The captain chose total cost as the measure (G6.1): the tokens of all
rounds until the run passes. The [LLM test](llm_test.md) measures it in phase 1, with the
diff size of three edit tasks.

---

## 10. The web

| Choice | Gives | Costs |
|---|---|---|
| the whole tool chain in Goose. Goose writes C, and `emcc` compiles the C to WebAssembly, as in earlier Goose projects (spec §18.1). | `check`, `fmt`, the mesher, the solver and the viewer run in a browser with no server | the browser limits memory, so large 3D models run only on the native target |
| studies on Web Workers | a study uses all cores of the machine in the browser | each worker loads the module and its memory |

---

## 11. Summary of the gaps

The captain decided each gap on 2026-10-30 ([grilling.md](grilling.md#answers)).

| Gap | Workloads | Decision | Phase | Question |
|---|---|---|---|---|
| standard-library physics | shells, beams, turbulence | physics templates in Skein source (spec §4.2) | 2 to 5 | G4.1 |
| radiation | heat transfer | a `radiation` item with view factors | 5 | G4.2 |
| shock capturing | compressible CFD | `limiter` for `fv` | 4 | G4.3 |
| element families | electromagnetics, porous media | `family: nedelec \| raviart_thomas` | 5 | G4.4 |
| number systems | electromagnetics, acoustics | `Real` and `Complex`, and the `harmonic` procedure (spec §5.1). `Quaternion` and `Dual` cut on 2026-10-31: no model used them. | 4 | G4.5, G4.5b |
| symmetry groups | periodic structures | `bloch` (spec §12.6). `cyclic` and `dihedral` cut on 2026-10-31: no model used them. | 4 | G4.5b |
| adaptive time steps | stiff ODEs | cut on 2026-10-31: no model used `adapt` | — | G4.6 |
| processes in a step | level-set reinitialization | out of scope for v1 | — | G4.7 |
| adjoint gradients | shape and topology optimization | adjoints from the SIR | 5 | G3.3 |
| exact STEP | export for CAM | no exact B-rep. STL and tessellated STEP only. | — | G5.2 |
