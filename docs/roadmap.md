# Skein roadmap

Each phase ends with a gate: a result that a person can measure. The next
phase starts only after the gate passes.

| Phase | Contents | Gate |
|---|---|---|
| 0 | this specification and the six samples | the captain approves the spec |
| 1 | front end: lexer, parser, type and unit checker, `fmt` with the layout rules for small diffs, `check`, diagnostics in TOON, `param` and `--set`, `fn`, the required `scales` line, `check --groups`, `check --study`, `skein help --llm` | the six samples pass `check` and `fmt --check`. Each femahoy model passes `check`, or it has a `// blocked-by:` mark that names the gaps that block it ([femahoy_coverage.md](femahoy_coverage.md#spec-gaps)). Phase 1 closes the gaps G6, G7, G9, G12, G16, G29 and G30. Each error in [checks.md §2](checks.md#2-errors) has a test. The [LLM test](llm_test.md) passes on total cost, and records the diff size of its edit tasks. |
| 2 | first solver build: `cg` for linear static and transient, the mesher (2D and 3D), CG, GMRES, ILU(0), Newton, VTU output, pre- and post-solve verification, `lint`, the `linear_elastic` and `heat` templates, the `method` defaults, `use` files, `shell`, `loft` and `pattern`, the study runner with cases first (§20), STL export, the kernel splitter (at most 60 lines for each generated function) and the other Power of Ten tests | sample 01 runs, and a `--set` study of it runs on all threads. The patch test and the plate with a hole pass. The femahoy gaps G5, G10, G11, G14, G15, G17, G18, G19, G20, G21, G24, G25 and G28 are closed, and their `// blocked-by:` marks are removed. |
| 3 | stabilized Stokes and Navier-Stokes (SUPG/PSPG), birth and death, submodels, AMG, LOBPCG for modal and buckling, restarts, the scaled SIR (§21), the `stokes` and `navier_stokes` templates, adaptive time steps with variable-order BDF, the work-stealing task pool (§20.5) | samples 02, 03, 04 and 06 run. The cavity matches Ghia et al. within the `expect` limit. |
| 4 | `fv`, then `dg`, then the meshless methods `rkpm` and `gfd`. `verify convergence`. Limiters and entropy viscosity. `Complex`, `Quaternion` and `Dual`, the `harmonic` procedure, the symmetry groups. The WebGL viewer. The browser target: Goose to C, then `emcc`. Tessellated STEP export. | sample 05 runs. Each method passes its manufactured solution test. The femahoy gaps G2, G8 and G26 are closed, and their `// blocked-by:` marks are removed. |
| 5 | plasticity, hyperelasticity, finite strain, contact, mesh adaptivity, FSI, turbulence, NURBS geometry, the shell, beam and turbulence templates, radiation with view factors, Nédélec and Raviart-Thomas elements, `minimize` with adjoint gradients, DFM checks and tolerance stack-ups with Monte Carlo (§22). If the schedule is late, cut DFM first, then the tessellated STEP export of phase 4 ([tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order)). | the benchmarks of [implementation.md §10](implementation.md#10-testing-and-the-verification-suite) pass. The femahoy gaps G1, G3, G4, G13, G22, G23 and G27 are closed, and their `// blocked-by:` marks are removed. |

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

## Risks

| Risk | Effect | What to do |
|---|---|---|
| Goose may have no atomic operations | the task pool of phase 3 needs atomic queues | check at the start of phase 3. A lock on each queue gives the same results, but slower. |
| Goose is young | the compiler or the runtime can find Goose faults | report them upstream, with a small test |
| a 3D mesher with sharp features and a quality guarantee is a research problem | phase 2 takes longer | start with 2D and with boxes, and measure each mesh (`mesh_quality`) |
| the LLM test fails | the syntax must change | phase 1 is cheap to change because it has no solver |
