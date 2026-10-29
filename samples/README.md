# Skein samples

These six models show the main parts of the language. Read the
[tutorial](../docs/tutorial.md) first. The samples are the reference
answers of the [LLM test](../docs/llm_test.md).

No compiler exists yet, so the samples do not run. Phase 1 of the
[roadmap](../docs/roadmap.md) makes each of them pass `skein check` and
`skein fmt --check`.

```
skein check 01_cantilever.skn
skein run 01_cantilever.skn
```

| Sample | What it shows |
|---|---|
| [01_cantilever](01_cantilever.skn) | A 3D steel cantilever with a hole: CSG geometry, `boundary` selections, a local mesh refinement, linear elasticity in strong form, quadratic tetrahedra, a static step, a reaction report. |
| [02_cavity](02_cavity.skn) | The lid-driven cavity at Re = 100: steady incompressible Navier-Stokes, two unknowns, SUPG/PSPG stabilization, a pressure reference at a point, a TOON table of the Ghia et al. (1982) data, and an `expect` for each row. |
| [03_birth_death](03_birth_death.skn) | An embankment built in three lifts and a berm removed: plane strain, gravity, three materials on one physics, top-level conditions, five steps with `birth` and `death`. |
| [04_submodel](04_submodel.skn) | A submodel of the hole in sample 01: `results` from another model, cut-boundary displacements, a fine local mesh, a peak stress report. |
| [05_meshless_plate](05_meshless_plate.skn) | A plate with a hole without a mesh: a point cloud, RKPM with stabilized nodal integration, Nitsche boundary conditions, plane stress with thickness, symmetry, and an `expect` against the Kirsch solution. |
| [06_cylinder_sweep](06_cylinder_sweep.skn) | Flow past a cylinder swept over the Reynolds number: a `param` with a range, `scales`, a viscosity computed from `Re`, a named `report` of the drag coefficient, and a `grid` study that runs the cases at the same time (spec §20, §21). |

The Ghia et al. values in sample 02 are from Table I of the paper, for
Re = 100. Check them against the paper before you use them as a benchmark.
