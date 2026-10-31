# Skein roadmap

Each phase ends with a gate: a result that a person can measure. The next
phase starts only after the gate passes.

| Phase | Contents | Gate |
|---|---|---|
| 0 | this specification and the six samples | the captain approves the spec |
| 1 | front end: lexer, parser, type and unit checker, `fmt` with the layout rules for small diffs, `check`, diagnostics in TOON, `param` and `--set`, `fn`, the required `scales` line, `check --groups`, `check --study`, `skein help --llm` | the six samples pass `check` and `fmt --check`. Each femahoy model passes `check`, or it has a `// blocked-by:` mark that names the gaps that block it ([femahoy_coverage.md](femahoy_coverage.md#spec-gaps)). Phase 1 closes the gaps G6, G7, G9, G12, G16, G29 and G30. Each error in [checks.md §2](checks.md#2-errors) has a test. The [LLM test](llm_test.md) passes on pass rate and total cost against a CalculiX deck that femahoy runs (`scripts/llm/workflow.py`), and records the diff size of its edit tasks. |
| 2 | first solver build: `cg` for linear static and transient, the mesher (2D and 3D), CG, GMRES, ILU(0), Newton, VTU output, pre- and post-solve verification, `lint`, the `linear_elastic` and `heat` templates, the `method` defaults, `use` files, `shell`, `loft` and `pattern`, the study runner with cases first (§20), STL export, the kernel splitter (at most 60 lines for each generated function) and the other Power of Ten tests | sample 01 runs, and a `--set` study of it runs on all threads. The patch test and the plate with a hole pass. The femahoy gaps G5, G10, G11, G14, G15, G17, G18, G19, G20, G21, G24, G25 and G28 are closed, and their `// blocked-by:` marks are removed. |
| 3 | stabilized Stokes and Navier-Stokes (SUPG/PSPG), birth and death, submodels, AMG, LOBPCG for modal and buckling, the scaled SIR (§21), the `stokes` and `navier_stokes` templates, the task pool (§20.5) | samples 02, 03, 04 and 06 run. The cavity matches Ghia et al. within the `expect` limit. The femahoy gap G36 is closed, and its `// blocked-by:` marks are removed. |
| 4 | `fv`, then `dg`, then the meshless methods `rkpm` and `gfd`. `verify convergence`. Limiters. `Complex`, the `harmonic` procedure, the `bloch` symmetry group. The WebGL viewer. The browser target: Goose to C, then `emcc`. Tessellated STEP export. | sample 05 runs. Each method passes its manufactured solution test. The femahoy gaps G2, G8 and G26 are closed, and their `// blocked-by:` marks are removed. |
| 5 | plasticity, hyperelasticity, finite strain, contact, FSI, turbulence, NURBS geometry, the shell, beam and turbulence templates, radiation with view factors, Nédélec and Raviart-Thomas elements, `minimize` with adjoint gradients. If the schedule is late, cut the tessellated STEP export of phase 4 first ([tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order)). | the benchmarks of [implementation.md §10](implementation.md#10-testing-and-the-verification-suite) pass. The femahoy gaps G1, G3, G4, G13, G22, G23 and G27 are closed, and their `// blocked-by:` marks are removed. |

The build order of the methods is `cg`, `fv`, `dg`, then meshless
(decision 5 in [decisions.md](decisions.md)).

## The end-to-end slice

One thin slice of phase 2 lands before phase 1 is complete. Reason: a Skein
model must give a measured result early. The slice is the package `runtime`
in `src/runtime/`. It has a box mesher with 27-node hexahedra (Q2), CG with
a Jacobi preconditioner, reactions and VTU output. The package `codegen` in
`src/codegen/` gives `skein run` for one linear static box: it generates a
Goose program that calls `runtime`. Phase 2 widens the slice to the
full runtime of [spec §18.4](skein_spec.md#184-the-runtime).

## Phase 1 in detail

Phase 1 adds no solver, except the end-to-end slice. Reason: the language is the product, and the
[LLM test](llm_test.md) can show a bad syntax before the solver costs time.

Phase 1 delivers:

1. the lexer and the parser of spec §2 and §3, with error recovery;
2. name resolution, types, dimensions and units (spec §4 to §7);
3. the geometry types and the region rules that need no mesh (spec §8);
4. the rules for fields, materials, physics, assumptions, conditions,
   methods and steps (spec §10 to §15), except the ones that need a mesh;
5. `skein fmt` and `skein check`, with TOON and JSON output;
6. a test for each error id.

## Phase 3 status

Phase 3 is in three tasks: flow, solid, and the mesher with the task pool.

| Item | Status |
|---|---|
| stabilized Stokes and Navier-Stokes (SUPG/PSPG) | done: `src/flow/`, steady, 2D, P1/P1 and P2/P1 triangles |
| the scaled SIR (§21) | done for the flow models: `src/flow/gen.goose` |
| the `stokes` and `navier_stokes` templates | done |
| the task pool (§20.5) | done for the cases of a study, with the memory ledger (`src/codegen/run.goose`), and for the threads in a case (`src/solvers/pool.goose`) |
| birth and death | done: `birth` and `death` regions of a step, `src/elements/gen.goose` |
| submodels | done: `results(file, step:, case:)`, the cut values and `submodel_cut`, `src/elements/source.goose`. The femahoy gap G34 is closed. |
| AMG | done: the default preconditioner of the elements solver, with the rigid-body modes for the eigen solves, `src/solvers/amg.goose` |
| LOBPCG for modal and buckling | done: `modal` and `buckling` steps and `eigenvalue(k)`, `src/elements/eigen.goose` |

Gate: `scripts/generated_c.sh` runs samples 02 and 06. Sample 02 passes
each `expect` against Ghia et al. within the limit of 0.02 m/s. The gap
G36 is closed. Samples 03 and 04 run with `skein run`, and each
`expect` passes. The gap G34 is closed.

## Phase 4 status

Each phase 4 lane changes only the rows of its own items.

| Item | Status |
|---|---|
| the browser target: Goose to C, then `emcc` | done for the simplex models: `skein build --target wasm` writes a page that meshes and solves in the browser, on the CPU with the implicit solvers (`src/codegen/build.goose`, `src/codegen/web.goose`). Sample 01 in the page gives the records and the VTU file of `skein run --set`, bit for bit. The GPU and the explicit solver are phase 4.5. |
| the WebGL viewer | done as `skein report` ([command_line.md §12](command_line.md#12-the-html-report)): one offline HTML file with the source, the solver tables and the results of each study case, `src/report/`. The page of `skein build --target wasm` solves in the browser. |
| the interposer warpage of Hisada et al. 2012 | the model is in [interposer_warpage.skn](../examples/package/interposer_warpage.skn), with the unknown RDL thickness, bump layer and chip offset as params. Its study `sweep` of 72 cases keeps the order of Figure 7, but at the best set only 3 of 9 cases are within 25 % (ratios 0.54 to 0.98), so its `expect` fails. |

## Risks

| Risk | Effect | What to do |
|---|---|---|
| Goose is young | the compiler or the runtime can find Goose faults | report them upstream, with a small test |
| a 3D mesher with sharp features and a quality guarantee is a research problem | phase 2 takes longer | start with 2D and with boxes, and measure each mesh (`mesh_quality`) |
| the LLM test fails | the syntax must change | phase 1 is cheap to change because it has no solver |
