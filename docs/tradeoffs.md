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
open questions go to the captain in [grilling.md](grilling.md).

---

## 0. The goals, in priority order

Skein has too many goals to build them all at once. This table ranks them.
Simulation is the core. A goal above the line is **protected**: no later
goal may weaken it. A goal below the line is **layered**: it builds on the
core, it comes later, and it can be cut if it costs too much. The captain
confirms or changes this ranking in [grilling.md round 0](grilling.md#round-0-the-ranking-of-the-goals).

| Rank | Goal | Class | Syntax cost | Compiler size | Build time | Runtime cost | Phase |
|---|---|---|---|---|---|---|---|
| 1 | Simulation core: strong-form physics, methods, mesher, solvers | protected | the base language | the largest part: lowering, methods, mesher, solvers | the most | none extra, it is the work | 2 to 5 |
| 2 | Units, static analysis and verification by default | protected | a type on each declaration, a unit on each literal | the checker, about a third of the front end | small | the post-solve checks cost a few percent of a solve | 1 (checks), 2 (verification) |
| 3 | LLM- and AXI-native syntax and output | protected | it shapes all syntax: types only at declarations, TOON, ids with fixes | small: the fix strategies and the TOON writer | none | none | 1 |
| 4 | Parameters, studies and concurrency | layered, high value | `param`, `study`, `report name =` | the study runner and the dependency analysis | small | memory grows with the number of workers | 1 (`param` in `check`), 2 (studies) |
| 5 | Dimensional analysis: groups and scales | layered | `scales`, `units nondim` | the π basis and the scaled SIR | small | scaled solves are equal or faster | 1 (`check --groups`), 3 (scaled SIR) |
| 6 | Web output: the tool chain and the viewer in the browser | layered, depends on Goose | none | the browser shell and the viewer | a second target | the browser limits memory | 4 |
| 7 | Power of Ten reliability | layered, mostly free | none: most rules hold by construction | tests and emitter checks | small | none | 2 onwards, as tests |
| 8 | Parametric parts and text-to-CAD | layered | `tol`, more geometry operations | export writers. A B-rep kernel for exact STEP. | medium for exact STEP | none for simulation | 2 (STL), 4 (tessellated STEP), 5 (exact STEP) |
| 9 | Design for manufacture | layered, can be cut | `manufacture`, `tolerance` | the DFM checks on distance functions | small | geometry checks only, fast | 5 |

Rules that follow from the ranking:

- A layered goal may add syntax, but it must not change the meaning of a
  model that does not use it. A model with no `study`, `scales`,
  `manufacture` or `tolerance` item is the same model as before.
- When two goals conflict, the higher rank wins. Example: a geometry
  operation that the mesher cannot mesh with good cells is not added for
  text-to-CAD (rank 8 against rank 1).
- The LLM-native goal (rank 3) is protected because it is the reason for
  the language, but it does not overrule correctness: units stay strict
  (rank 2) even when they cost tokens.
- Ranks 8 and 9 are the first to cut if the schedule is late. Rank 7 costs
  almost nothing, so it stays.

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
| C3. Declarative: no loops, no functions, no mutable state | §1, §7.4 |
| C4. One model per file, one namespace | §1, §4 |
| C5. Geometry from distance functions, built-in mesher, no mesh import | §8, §9 |
| C6. Method-neutral SIR, methods chosen per physics | §13, §16 |
| C7. Code generation to Goose, all solvers in Goose | §18 |
| C8. Verification on by default | checks.md |

---

## 2. The design against the workloads

| Workload | What works | Cost or gap | Fix in the spec, or proposal |
|---|---|---|---|
| Linear structures (static, modal, buckling) | C1 with `hooke`, `cg`, all procedures | none of note | — |
| Nonlinear structures (plasticity, hyperelasticity, large deformation) | `plastic`, `diff(W, C)`, `finite_strain` (§10.6) | a user law with many state variables needs a return map that the compiler must derive. A general return map for any yield function is hard. | phase 5. Start with von Mises and Drucker-Prager, and give a clear error for a law the compiler cannot map. |
| Shells and beams | lower-dimension regions (§8.3) | the strong form of a shell (Reissner-Mindlin, Kirchhoff-Love) is long. Users do not write it from memory. | **Gap.** Proposal: standard-library physics for shells and beams, written in Skein. Question G4.1. |
| Contact | `contact` (§12.4) | not a PDE: an inequality. It does not fit the strong form, so it is a separate item. | accepted cost |
| Heat transfer | C1 with `fourier`, `transient` | radiation between surfaces needs view factors, a non-local term | **Gap.** Proposal: `radiation` item, phase 5. Question G4.2. |
| Incompressible CFD | `supg_pspg`, `vms`, `fv` with Rhie-Chow | none of note for laminar flow | — |
| Turbulent CFD | `turbulent(model)` assumption | RANS models add equations with wall functions. LES needs `vms`. Each model is many lines. | standard-library physics for k-ω SST and Spalart-Allmaras. Question G4.1. |
| Compressible CFD with shocks | `dg`, `fv` with Roe flux | shocks need limiters or artificial viscosity, and the spec has no setting for them | **Gap.** Proposal: `limiter` (`minmod`, `barth_jespersen`) for `fv`, and `shock_capturing` (entropy viscosity) for `dg`. Question G4.3. |
| Multiphysics (thermal-structural) | two physics that share fields (§10.4), `monolithic` or `staggered` | none of note | — |
| FSI | `interface` with ALE (§10.7) | mesh motion and large deformation of the fluid mesh are hard to make robust | phase 5, accepted cost |
| Electromagnetics | `curl` exists | `curl curl` problems need H(curl) (Nédélec) elements. The method table has no element family setting. Time-harmonic problems need complex numbers. | **Gap.** Proposal: `family: nedelec` in `cg`, and a `Complex<D>` type. Questions G4.4 and G4.5. |
| Acoustics (Helmholtz) | the equation is easy to write | frequency domain needs complex numbers and absorbing boundaries | **Gap.** Same `Complex<D>` proposal. Question G4.5. |
| Porous media (Darcy, consolidation) | mixed forms with two unknowns | a mixed Darcy form needs H(div) (Raviart-Thomas) elements | **Gap.** Same element family proposal. Question G4.4. |
| ODE systems (control, lumped models, chemistry) | global unknowns (§10.5), `ode` method | stiff chemistry needs adaptive implicit time steps with error control. The spec has fixed `dt` only. | **Gap.** Proposal: `transient(adapt: { rtol, atol })` with BDF of variable order. Question G4.6. |
| Phase field, reaction-diffusion | C1, `transient` | fourth-order terms (Cahn-Hilliard) need a mixed form or C1 elements | accepted. Write the mixed form with two unknowns. |
| Free surface (level set, VOF) | a level set is a PDE | reinitialization is a process, not an equation | **Gap.** Question G4.7. |
| Meshless (large deformation, fracture) | `rkpm`, `gfd` | fragmentation needs particle methods with contact. That is closer to DEM, which is not a PDE. | out of scope for v1 |
| Parameter studies | `param`, `study` (§20) | each case is a whole run | §20.4 reuses the mesh |
| Optimization | `study … minimize` (§20.2) | gradient-free methods need many cases. Shape and topology optimization need adjoint gradients. | **Gap** for gradients. Proposal: adjoints from the SIR, phase 5. Question G3.3. |
| Uncertainty quantification | `study … sample` with `lhs` or `sobol` | no statistics of the results in the tool | the table is TOON. Compute the statistics outside, or add `mean`, `std` and quantiles to the study output. Question G3.4. |
| General PDEs with no physics | `units nondim` (§6.5) | no dimension checks | accepted cost, the check report states it |

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
| each unknown has one physics and one equation | counts | §10 |
| assumptions do not conflict | assumption table | §11 |
| conditions do not conflict and are on the boundary | condition rules | §12 |
| the method can discretize the physics, with a stable pair | method rules | §13 |
| each study value is in the range of its parameter | ranges | §20 |
| the scales span the dimensions | rank of the dimension matrix | §21 |
| which parts of the model each parameter changes | dependency analysis | §20.4 |

What static analysis cannot prove, and what Skein does instead:

| Not provable before the run | Instead |
|---|---|
| the mesh resolves the solution | `verify convergence`, `adapt` |
| the nonlinear solve converges | `newton_converged` check, `ramp` |
| the time step is stable for a nonlinear problem | `dt_unstable` lint is an estimate. The post-solve energy check finds growth. |
| the model is the right model of the real part | `expect` against tests or benchmarks. This is validation, not verification. |

Cost of C3: no user functions. A user who needs a helper for a long law
must repeat it, or use the standard library. Question G1.2.

---

## 4. Parametric parts, text-to-CAD and export

Recent text-to-CAD tools have an LLM write a short parametric script, in
the style of CadQuery or build123d, from a sentence. Skein geometry has the
same shape (spec §22.1): named parameters drive a CSG tree.

What Skein adds to that trend:

- **One set of parameters** drives the part, the simulation and the study.
  A sweep of the wall thickness changes the part, the mesh, the stress and
  the manufacturing checks in one run.
- **Types and units** on each dimension of the part, so an LLM that writes
  `2` for a wall gets an error, not a part 2 m thick.
- **Ranges** on parameters, so a study cannot make a part with a hole wider
  than its wall.

What it costs:

| Cost | Why | Proposal |
|---|---|---|
| distance-function CSG is weaker than a B-rep kernel for some features | sketches with constraints, lofts, shells of constant thickness, variable fillets and text are hard on distance functions | add `shell(s, t)` (an offset) and `loft` as distance functions. Keep constraint sketches out of v1. Question G5.1. |
| exact STEP export needs a B-rep | a CSG tree of distance functions has no exact faces and edges | phase 4 exports STEP AP242 with tessellated geometry, which most CAM and inspection tools read. An exact B-rep for the analytic primitives (planes, cylinders, spheres, cones, tori) is phase 5. Question G5.2. |
| no import of CAD | decision 8 | accepted. An LLM can rewrite a simple part as Skein. A complex existing part is out of scope. |

---

## 5. Design for manufacture

The checks of spec §22.3 run on the parametric part with the process that
the file declares. They are in the `verify` family, so they report the same
way as the other checks, and `allow` with a reason turns one off.

| Choice | Gives | Costs |
|---|---|---|
| checks on the distance function | wall thickness and draft come directly from φ and ∇φ, with no B-rep | the result has the accuracy of the sampling. Thin features near the size of the sample can hide. The checks sample at a fraction of `min_wall`. |
| the process in the file | the checks know the rules, and a reviewer sees them | a part made by two processes (cast, then machined) needs two `manufacture` items on one region. The spec allows it. |
| tolerance stack-ups by derivatives | exact sensitivities from the expression, worst case and RSS | linear only. A stack-up through a nonlinear geometry relation needs a sampling method. Proposal: `method: monte_carlo` with `study … sample`. Question G5.3. |
| DFM on each study case | a sweep shows which designs cannot be made | the checks cost time for each case. They are geometry-only, so they are fast compared with a solve. |

What the DFM checks do not cover: cost, cycle time, machine choice and
fixturing. These are out of scope.

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
| 4 | functions of at most 60 lines | **by construction**: no functions | runtime functions obey it. Generated kernels are straight-line code and can be longer. This is a stated exception: they have no branches and come from one template each. |
| 5 | at least 2 assertions per function | verification checks are on by default (checks.md) | each runtime function asserts its preconditions. Each kernel asserts its array sizes once at entry. |
| 6 | data in the smallest scope | **partly**: one namespace per model (spec §4) is a deliberate exception. A model is short, and one namespace keeps `sigma` usable in a step. | Goose scopes free data at scope exit |
| 7 | check each return value | no calls that return status | each solver call returns a status that the plan checks. An unchecked status is a compile error in the emitter. |
| 8 | limited preprocessor | **by construction**: none | Goose has no preprocessor |
| 9 | limited pointers | **by construction**: none | Goose has no raw pointers. References are checked by the Goose compiler. |
| 10 | all warnings on, static analysis clean | `skein check` and `skein lint` with no warnings on each sample and benchmark | the Goose compiler with `--debug` on the test suite, and `checkfree fn` on the hot kernels, which proves that no bounds check is left |

Rules 1, 2, 3, 8 and 9 hold for Skein source by construction: the language
cannot express a violation. Rule 4 has a stated exception for generated
kernels, and rule 6 for the one namespace. Rules 5, 7 and 10 are
enforced by tests and by the emitter.

---

## 7. Parameters and concurrency

| Choice | Gives | Costs |
|---|---|---|
| only `param` can change (§20.1) | a reviewer finds each input. The dependency analysis is exact. | a user must mark a value `param` before a study can vary it |
| Go-style workers over typed queues (§20.3) | no shared memory, so no data races. A failed case does not stop the study. The same design works on native threads and Web Workers. | each worker holds its own copy of the shared mesh. Memory grows with the number of workers. §20.5 limits the workers by memory. |
| parallel cases before parallel assembly (§20.5) | cases scale almost perfectly | one large case on a machine with many threads gets threads only when it is alone |

---

## 8. Strict dimensions and dimensionless work

Strict dimensions (C2) catch the most common modelling fault. But engineers,
mostly in fluids, think in dimensionless groups. The spec keeps both:

| Need | How | Spec |
|---|---|---|
| vary Re, not ρ, U, L, μ | `param Re: Real` and a dimensional `let mu = rho*U*D/Re` | §21.1 |
| the solver in scaled variables | `scales D, U, rho` | §21.2 |
| know the groups of a model | `check --groups`, Buckingham π | §21.3 |
| a model with no physical units at all | `units nondim` | §6.5 |

Costs:

- `units nondim` removes all dimension checks. It is a deliberate choice in
  the header, and the report says so.
- The π basis is not unique. The tool prefers the `scales` as repeated
  variables and names known forms, but a user can see a group in a
  different form than the one in the textbook.
- A compressible flow needs more scales (a temperature, a speed of sound).
  The `scales_incomplete` error lists the missing dimension.

---

## 9. LLM-native and AXI-friendly

| Choice | Gives | Costs |
|---|---|---|
| types only at declarations | equations stay short and look like a textbook | each `let` has a type, which costs tokens |
| TOON tables in source and output | fewer tokens than JSON for uniform rows | a second syntax inside the language |
| stable ids with a fix in each diagnostic | an LLM repairs a model from the output alone | each new rule needs a fix strategy |
| `allow` only in the file | an LLM cannot turn off a check from the command line to make a run pass | a user must edit the file |
| strong form | an LLM knows textbook equations well | the reference answers are 1.4 to 2.9 times larger than APDL decks in bytes ([llm_test.md](llm_test.md) §5) |

The size cost is the largest open question. The [LLM test](llm_test.md)
measures it in phase 1. Question G6.1.

---

## 10. The web

| Choice | Gives | Costs |
|---|---|---|
| the whole tool chain in Goose, compiled to WebAssembly | `check`, `fmt`, the mesher, the solver and the viewer run in a browser with no server | Goose on WebAssembly is a proposal today (spec §18.1). The browser limits memory, so large 3D models run only on the native target. |
| studies on Web Workers | a study uses all cores of the machine in the browser | each worker loads the module and its memory |

---

## 11. Summary of the gaps

| Gap | Workloads | Proposal | Question |
|---|---|---|---|
| standard-library physics | shells, beams, turbulence | physics templates in Skein source | G4.1 |
| radiation | heat transfer | a `radiation` item | G4.2 |
| shock capturing | compressible CFD | `limiter`, `shock_capturing` | G4.3 |
| element families | electromagnetics, porous media | `family: nedelec \| raviart_thomas` | G4.4 |
| complex numbers | electromagnetics, acoustics | `Complex<D>` and a `harmonic` procedure | G4.5 |
| adaptive time steps | stiff ODEs | `transient(adapt:)` with variable-order BDF | G4.6 |
| processes in a step | level-set reinitialization | a `process` statement in a step | G4.7 |
| adjoint gradients | shape and topology optimization | adjoints from the SIR | G3.3 |
| exact STEP | export for CAM | B-rep of the analytic primitives | G5.2 |
