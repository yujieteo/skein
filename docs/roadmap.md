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
| 5 | plasticity, hyperelasticity, finite strain, contact, FSI, turbulence, NURBS geometry, the shell, beam and turbulence templates, radiation with view factors, Nédélec and Raviart-Thomas elements, `minimize` with adjoint gradients. The tessellated STEP export of phase 4 is done, so it is no longer a cut candidate. The next cut if the schedule is late is not decided, and the captain decides it ([tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order)). | the benchmarks of [implementation.md §10](implementation.md#10-testing-and-the-verification-suite) pass. The femahoy gaps G1, G3, G4, G13, G22, G23 and G27 are closed, and their `// blocked-by:` marks are removed. |

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
| `Complex`, the `harmonic` procedure, the `bloch` symmetry group, `pml` | done for a `Complex` scalar unknown: `eq div(grad(u)) + c*u = s`, solved by GMRES with ILU(0) at each frequency (`src/elements/wave.goose`). [harmonic_duct_tri.skn](../examples/elements/harmonic_duct_tri.skn) and [bloch_cell_tri.skn](../examples/elements/bloch_cell_tri.skn) match their closed forms within 0.2 %, and the PML of [pml_duct_tri.skn](../examples/elements/pml_duct_tri.skn) reflects less than 1e-4 of the wave. A `Complex` vector unknown, `rms` and a step other than `harmonic` with a `Complex` unknown do not run yet. |
| the femahoy gaps G2, G8 and G26 | closed: `harmonic(damping:)`, `rms(f, psd)`, the built-in `omega` and the method setting `pml` pass `check` ([skein_spec.md §13](skein_spec.md#13-methods), [§14.2](skein_spec.md#142-procedures), [§14.5](skein_spec.md#145-output-and-report)). Their `// blocked-by:` marks are removed. |
| the interposer warpage of Hisada et al. 2012 | the model is in [interposer_warpage.skn](../examples/package/interposer_warpage.skn), with the unknown RDL thickness, bump layer and chip offset as params. Its study `sweep` of 72 cases keeps the order of Figure 7, but at the best set only 3 of 9 cases are within 25 % (ratios 0.54 to 0.98), so its `expect` fails. |
| `fv` | done for one scalar unknown, `ct ddt(u) + div(a*u) - d*div(grad(u)) = f`, on a `rect` of `quad` cells, `src/grid/` ([skein_spec.md §17.3](skein_spec.md#173-finite-volume-fv)). `static` steps and `transient` steps with `ssprk3`. The manufactured solutions of [mms_fv_static.skn](../examples/grid/mms_fv_static.skn) and [mms_fv_transient.skn](../examples/grid/mms_fv_transient.skn) converge with the L2 orders 2.40 and 1.98. `pv_coupling` and flow do not run yet. |
| `dg` | done for the same equation and grid, `order` 1 to 3, `src/grid/` ([§17.4](skein_spec.md#174-discontinuous-galerkin-dg)). The manufactured solution of [mms_dg_static.skn](../examples/grid/mms_dg_static.skn) converges with the L2 order 3.16 for `order: 2`. The tests of `src/grid/test.goose` measure about 2, 3 and 4 for the orders 1, 2 and 3. |
| limiters | done for `fv`: `minmod`, `van_leer` and `venkatakrishnan`. On a square wave, `minmod` and `van_leer` make no overshoot, and `venkatakrishnan` less than half the overshoot of `none` (`src/grid/test.goose`). |
| `verify convergence` | done for the `fv` and `dg` models, from the L2 error of each `expect` of two fields ([checks.md §4.3](checks.md#43-on-request)). The order of a `report` does not run yet. |
| `rkpm` | done in `src/meshless/` ([skein_spec.md §17.5](skein_spec.md#175-meshless)) for one scalar unknown, `div(a*u) - d*div(grad(u)) = f`, and for 2D linear elasticity: the linear basis, nodal integration on the quadrilaterals of each nodal domain, and Nitsche conditions. The manufactured solution of [mms_rkpm.skn](../examples/meshless/mms_rkpm.skn) converges with the L2 order 2.02 on a cloud that is not a grid, and the patch test of `src/meshless/test.goose` is exact to 1e-15. `basis: quadratic` and `integrate: gauss` do not run yet. |
| `gfd` | done in `src/meshless/` for the same scalar equation: a weighted fit of order 2 on 12 points. The manufactured solution of [mms_gfd.skn](../examples/meshless/mms_gfd.skn) converges with the L2 order 2.03. Elasticity does not run with `gfd` yet. |
| sample 05 | done: [05_meshless_plate.skn](../samples/05_meshless_plate.skn) runs with `rkpm` on the 7,171 nodes of a P1 mesh of the plate, in 12 s. Its `expect` passes: σxx at the top of the hole is 298 MPa, 0.6 % from the Kirsch value 3S. Poisson disk points do not run yet. |
| tessellated STEP AP242 export | done: `skein export --format step` writes a `TESSELLATED_SHELL` with the triangles of the STL export, `src/mesher/stl.goose`. No STEP reader runs on the desktop, so the test in `src/mesher/test.goose` is a schema check: the Part 21 frame, the `#n` references, the index range and the triangle count against the STL. |

## Phase 5 status

Each phase 5 lane changes only the rows of its own items.

| Item | Status |
|---|---|
| the shell template | done for a composite `reissner_mindlin` shell on a part of the side of a cylinder, `src/shells/`: MITC4 quads for `order: 1` and MITC9 quads for `order: 2`, with a nested dissection solve for order 2. The shell can have cuts (`skin - door_box(0 m)`), a material on each part (`material skin_lam on plain, pad_lam on pad`) and `refine` curves. [thin_cylinder_shell.skn](../examples/shells/thin_cylinder_shell.skn) and [scordelis_lo_roof.skn](../examples/femahoy/scordelis_lo_roof.skn) run. The 787 global model [fuselage.skn](studies/b787/skein/fuselage.skn) (order 2, 714,834 dofs, 1.9 GB) passes C1 to C5, and its submodel [door_corner.skn](studies/b787/skein/door_corner.skn) passes C6. |
| the beam template | done for `timoshenko` beams of 2 nodes on `line`s, rings and their cuts, alone or as offset beams tied to the shell, with a material on each region. The benchmarks of [implementation.md §10](implementation.md#10-testing-and-the-verification-suite) pass: the tip deflection of [cantilever_beam.skn](../examples/shells/cantilever_beam.skn) is P L³/(3 E I) + P L/(k G A) within 0.1 % for 4 lengths, the first two frequencies of [beam_modes.skn](../examples/shells/beam_modes.skn) are within 1 % of the Euler-Bernoulli values, and the load factor of [beam_buckling.skn](../examples/shells/beam_buckling.skn) is within 1 % of the Euler load. In [stiffened_barrel.skn](../examples/shells/stiffened_barrel.skn), the hoop and axial stresses of a skin with stringers and frames are within the sheet values of the 787 study. Modes and buckling of a shell do not run. |
| the femahoy gaps G3, G31 and G33 | closed in the spec ([§4.3](skein_spec.md#43-shells-beams-and-laminates), [§12.1](skein_spec.md#121-form), [§15.1](skein_spec.md#151-results)), and `skein run` runs them in the two 787 models: plies, two laminates, offset beams with the tie, `N`, `eps_m`, `Nx`, the edge load `N*n`, and `global.w`, `global.theta` and `global.e_hoop` of the global model in the submodel. A beam force `F_b*n = P` on a `Point` is in §12.1. |
| the femahoy gaps G32 and G35 | closed: `a & b` of a `Surface` or a `Curve` and a `Solid` has the lower dimension ([§8.4](skein_spec.md#84-region-algebra)), and `copy(r)` is the index of the copy of a pattern ([§8.1](skein_spec.md#81-distance-functions)). The 787 models pass `check` with no `// blocked-by:` mark, and `skein run` cuts the skin, the stringers and the frames with the door and the pad-up solids. In [window_barrel.skn](../examples/shells/window_barrel.skn), `copy(windows)` gives the window at each point of the window edges, and the y reaction is the pressure resultant of the closed barrel within 0.01 %. |
| radiation with view factors | done for 2D steady heat on `cg`: the natural bc `dot(q, n) = radiation(eps)` makes the faces of the bc one gray diffuse enclosure, with view factors by crossed strings and a check of their row sums (`src/elements/radiation.goose`). [radiation_annulus_tri.skn](../examples/elements/radiation_annulus_tri.skn) matches the closed form of two concentric tubes within 0.001 %. 3D and transient steps do not run yet. |
| turbulence and the turbulence template | done: the `k_omega_sst` template (Menter 1994) on P1 triangles with SUPG/PSPG, pseudo-transient continuation, and automatic wall values of k and w (`src/flow/rans.goose`). A `bloch` group with `k: 0` makes a rect periodic. The centre velocity of [channel_sst_tri.skn](../examples/flow/channel_sst_tri.skn) at Re_tau = 395 is 20.86 u_tau, 3.8 % above the DNS of Moser et al. (1999). Inlet values of k and w do not run yet. |
| FSI | done for a steady flow and a `linear_elastic` solid in plane strain: an `interface` with `mesh_motion: harmonic`, a partitioned solve with Aitken relaxation, and one mesh of the two shapes with shared interface nodes (`src/flow/fsi.goose`). [turek_hron_fsi1_tri.skn](../examples/flow/turek_hron_fsi1_tri.skn) gives drag 14.298 N/m (+0.02 %), lift 0.758 N/m (-0.7 %) and uy at point A 0.816 mm (-0.6 %). ux at A is 4.5 % high, because the small strain does not shorten the flag as it bends. Its `expect` waits for a `finite_strain` St. Venant-Kirchhoff solid, which is not merged. Transient FSI does not run yet. |
| the femahoy gap G1: topology optimization | closed: the study procedure `topology` ([skein_spec.md §20.2](skein_spec.md#202-study)) runs SIMP, the sensitivity filter and the optimality criteria update with a volume limit, on a `rect` of `quad` cells, `src/optim/`. The derivatives come from one adjoint solve for each iteration, and `src/optim/test.goose` checks them against central differences. The study `shape` of [topology_mbb_beam.skn](../examples/femahoy/topology_mbb_beam.skn) gives c = 216.814 in 106 iterations (Andreassen et al. 2011: 216.81). With `threshold: 0.5`, its study `solid` gives c = 193.274 for the solid of the cells at 0.5 or more (femahoy: 193.27). A solid with a smooth boundary and its STEP file do not run yet. |
| NURBS geometry | done in 2D: `nurbs(d, knots(k…), p1, w1, …)` of degree 1 to 3 ([skein_spec.md §8.1](skein_spec.md#81-distance-functions)), and `fill` of it in the mesher, as rational Bézier spans from knot insertion (`src/mesher/geometry.goose`). In `src/mesher/test.goose`, each point of the rational circle is on the circle within 1e-12 of the radius, and its area is within 1e-9 of πr². [kirsch_nurbs_tri.skn](../examples/elements/kirsch_nurbs_tri.skn) has a NURBS hole and gives 298.99 MPa at the top of the hole (Kirsch: 300 MPa). NURBS surfaces and a NURBS curve in 3D do not run yet. |
| `minimize` with adjoint gradients | done for the models of the optim lane: `minimize(r, method: gradient)` ([skein_spec.md §20.2](skein_spec.md#202-study)), a projected gradient with the steps of Barzilai and Borwein, in `src/codegen/run.goose`. Each case solves the step and its adjoint once, and prints dJ/dp for each param, through E, nu, t, the loads and the objective. [stiffness_bar.skn](../examples/optim/stiffness_bar.skn) finds E = 10.0017 MPa in 10 cases (closed form: 10 MPa). On a model of 4 params, each derivative agrees with a central difference to 1e-7 relative. `nelder_mead`, `bobyqa`, and params of the geometry do not run yet. |

No femahoy gap belongs to the radiation, turbulence and FSI items, and
implementation.md §10 has no benchmark for them.

## Risks

| Risk | Effect | What to do |
|---|---|---|
| Goose is young | the compiler or the runtime can find Goose faults | report them upstream, with a small test |
| a 3D mesher with sharp features and a quality guarantee is a research problem | phase 2 takes longer | start with 2D and with boxes, and measure each mesh (`mesh_quality`) |
| the LLM test fails | the syntax must change | phase 1 is cheap to change because it has no solver |
