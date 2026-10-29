# Design decisions

The design of Skein comes from a review with the captain on 2026-10-29.
This index links each of the 17 decisions to the parts of the
documentation that implement it.

| # | Decision | Where |
|---|---|---|
| 1 | The name is Skein. The repository stays private until the captain approves the spec. The licence is Apache 2.0. | [README](../README.md), [LICENSE](../LICENSE) |
| 2 | General PDEs and ODEs from the start. FEM and CFD are the main uses. | spec [§10](skein_spec.md#10-fields-materials-and-physics), [§10.5](skein_spec.md#105-odes-and-global-unknowns) |
| 3 | One language for geometry, mesh or point cloud, physics, assumptions, methods, steps, solvers and output. | spec [§1](skein_spec.md#1-core-concepts) |
| 4 | Physics in strong form with constitutive laws and named assumptions. The compiler derives the weak form. | spec [§10.3](skein_spec.md#103-physics), [§11](skein_spec.md#11-modelling-assumptions), [§17.1](skein_spec.md#171-continuous-galerkin-cg) |
| 5 | A method-neutral core with no mesh. CG, DG, FV and meshless methods. Explicit and implicit time schemes. Build order CG, FV, DG, meshless. | spec [§13](skein_spec.md#13-methods), [§14.2](skein_spec.md#142-procedures), [§16](skein_spec.md#16-the-intermediate-form-sir), [§17](skein_spec.md#17-lowering-by-method), [roadmap](roadmap.md) |
| 6 | Rust-like syntax. Every type explicit, no inference. Declarative and memory-safe. Units in the type system. Equations as in Modelica. Easy for LLMs: types only at declarations, TOON tables. | spec [§1](skein_spec.md#1-core-concepts), [§3](skein_spec.md#3-grammar), [§5](skein_spec.md#5-types), [§5.6](skein_spec.md#56-tables) |
| 7 | Dimension types, literals with units, declared unit systems. Compile-time checks of units, rank, regions, assumptions and boundary conditions. | spec [§6](skein_spec.md#6-dimensions-and-units), [§19](skein_spec.md#19-static-semantics-checklist), [checks.md §2](checks.md#2-errors) |
| 8 | Geometry from distance functions (CSG), lower-dimension entities, a feature-aware mesher in the browser, NURBS later. No mesh import. No Gmsh code. | spec [§8](skein_spec.md#8-geometry-and-regions), [§9](skein_spec.md#9-meshes-and-point-clouds), [implementation §7](implementation.md#7-the-mesher) |
| 9 | Advanced analysis in the spec from the start: birth and death, submodels, CFD stabilization, contact, adaptivity, restarts, coupling, modal, buckling, plasticity, hyperelasticity, large deformation. | spec [§10.4](skein_spec.md#104-coupling), [§10.6](skein_spec.md#106-inelastic-materials), [§10.7](skein_spec.md#107-fsi-interfaces), [§12.4](skein_spec.md#124-contact), [§14](skein_spec.md#14-steps), [§15.2](skein_spec.md#152-submodels), [§17.2](skein_spec.md#172-stabilization) |
| 10 | Ordered Abaqus-style steps with a procedure, loads, birth and death, and output. Restart from a named step. | spec [§14](skein_spec.md#14-steps), [§14.7](skein_spec.md#147-restarts) |
| 11 | The compiler writes Goose. Goose builds C and WebAssembly. All solvers in Goose. | spec [§18](skein_spec.md#18-code-generation-to-goose), [implementation §6](implementation.md#6-the-runtime) |
| 12 | The compiler is written in Goose, so the toolchain runs in the browser. | spec [§18.5](skein_spec.md#185-the-compiler-in-goose), [implementation §8](implementation.md#8-the-browser) |
| 13 | One tool `skein` with `check`, `fmt`, `lint`, `verify`, `build`, `run`. TOON output, stable ids, line and column, a fix. A language server. | [command_line.md](command_line.md), [checks.md §1](checks.md#1-how-diagnostics-work) |
| 14 | Verification by default, before and after the solve. Convergence on request. `allow <check>: "<reason>"` only in the file. | [checks.md §4](checks.md#4-verification-checks) |
| 15 | The lint rules: unused parts, loads on regions with no cells, poor elements, unstable time steps, missing stabilization. | [checks.md §3](checks.md#3-lint-rules) |
| 16 | Results in a WebGL viewer and as VTU for ParaView. | [implementation §9](implementation.md#9-results-and-the-viewer), [command_line.md §1](command_line.md#1-commands) |
| 17 | A verification suite of manufactured solutions and published benchmarks, not tied to one company. | [implementation §10](implementation.md#10-testing-and-the-verification-suite) |

## Later steers

The captain added these goals on 2026-10-29, after the 17 decisions. The
[tradeoffs](tradeoffs.md#0-the-goals-in-priority-order) rank them below the
simulation core. The [grilling record](grilling.md) holds the open choices.

| Steer | Where |
|---|---|
| Tradeoffs against diverse problems | [tradeoffs §2](tradeoffs.md#2-the-design-against-the-workloads), [§11](tradeoffs.md#11-summary-of-the-gaps) |
| Parameter sweeps that run at the same time, in the style of Go | spec [§20](skein_spec.md#20-parameters-and-studies), [tradeoffs §7](tradeoffs.md#7-parameters-and-concurrency) |
| Dimensionless work: groups, scales, Buckingham π | spec [§21](skein_spec.md#21-dimensional-analysis-and-scaling), [tradeoffs §8](tradeoffs.md#8-strict-dimensions-and-dimensionless-work) |
| A grilling step with a design tree and recommended answers | [grilling.md](grilling.md) |
| LLM-native, AXI-friendly, runs on the web | [tradeoffs §9](tradeoffs.md#9-llm-native-and-axi-friendly), [§10](tradeoffs.md#10-the-web) |
| Parametric text-to-CAD, STEP and STL export, no mesh import | spec [§22.1](skein_spec.md#221-parametric-parts), [§22.5](skein_spec.md#225-export-for-manufacture), [tradeoffs §4](tradeoffs.md#4-parametric-parts-text-to-cad-and-export) |
| DFM checks with the process in the file | spec [§22.2](skein_spec.md#222-manufacture) to [§22.4](skein_spec.md#224-tolerances-and-stack-ups), [checks.md §4.2a](checks.md#42a-design-for-manufacture) |
| The Power of Ten rules for Skein and the generated Goose | [tradeoffs §6](tradeoffs.md#6-reliability-the-power-of-ten) |
| Too many goals: rank them, with simulation as the core | [tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order), [grilling round 0](grilling.md#round-0-the-ranking-of-the-goals) |

## Choices made in the spec

These choices are not in the 17 decisions. The spec gives a reason for each.

| Choice | Reason | Where |
|---|---|---|
| A newline ends a statement. No semicolons. | a semicolon costs one token on each line and tells the parser nothing | spec [§2.2](skein_spec.md#22-lines-and-statements) |
| Tables omit the TOON row count `[N]` | the count changes with each new row, and the checker finds a short row without it | spec [§3](skein_spec.md#3-grammar) |
| Declaration before use | one pass for names, and each error at the line of the use | spec [§3.1](skein_spec.md#31-order-of-items) |
| One namespace per model | a step can use `sigma` from a physics block with no prefix | spec [§4](skein_spec.md#4-names-and-scope) |
| `0` is zero of any type | `u = 0` is the most common condition | spec [§5.7](skein_spec.md#57-assignability) |
| No type aliases and no user functions | an alias puts one fact in two places, and a function body needs loops | spec [§6.6](skein_spec.md#66-what-skein-leaves-out), [§7.4](skein_spec.md#74-functions) |
| Diagnostic ids are names, not numbers | an LLM and a reviewer read `rigid_body` without a lookup | [checks.md](checks.md) |
| `allow` needs a reason in the file, with no command-line switch | a switch lets an LLM turn off a check to make a run pass | [checks.md §1](checks.md#1-how-diagnostics-work) |
| Dead cells leave the system, with no small stiffness factor | no ill conditioning from the dead cells | spec [§14.3](skein_spec.md#143-birth-and-death) |
