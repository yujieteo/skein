# Skein checks: errors, lint rules and verification

<details>
<summary>Table of contents</summary>

- [1. How diagnostics work](#1-how-diagnostics-work)
- [2. Errors](#2-errors)
- [3. Lint rules](#3-lint-rules)
- [4. Verification checks](#4-verification-checks)
- [5. Run errors](#5-run-errors)

</details>

This is the catalogue of every diagnostic that the Skein tools give. Each
diagnostic has one **id**, a snake_case name. An id never changes and is
never used again for a different rule. The tools print the id, the file,
the line, the column, a message and a suggested fix. The
[command line](command_line.md#6-what-skein-prints) gives the output format.

---

## 1. How diagnostics work

| Class | Tool that gives it | Effect | Can `allow` turn it off? |
|---|---|---|---|
| `error` | `skein check` (and all tools after it) | the program is not valid | no |
| `lint` | `skein lint` | the program is valid but probably wrong or slow | yes |
| `verify` | `skein run` and `skein verify` | a check on the model before the solve, or on the results after it | yes |
| `run` | `skein run` | the run stopped | no |

`skein run` does `check`, the pre-solve verification, the solve and the
post-solve verification, in that sequence. It does not need a separate
command for verification. A failed `verify` check stops the run before the
solve, or gives a non-zero exit status after it.

To turn off a lint rule or a verification check, write `allow` in the
model, with a reason:

```skein
allow rigid_body: "the springs in the step remove the rigid-body motion; the check sees only essential conditions"
```

An `allow` with an empty reason, or with an id that does not exist, is an
error [`allow_reason`]. The `run` report lists each `allow` and its reason.

Reason for the rule: a check that a user can turn off with a flag on the
command line is a check that an LLM turns off to make a run pass. A reason
in the file stays in the review.

---

## 2. Errors

`skein check` gives these. The column "Spec" points into the
[specification](skein_spec.md).

### 2.1 Syntax and names

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `syntax` | §3 | the text does not match the grammar | `let E Pressure = 1 Pa` | the expected token, for example "insert `:`" |
| `unknown_name` | §4 | a name is not declared | `hooke(E, nu, strain(v))` with no `v` | the nearest declared name |
| `use_before_decl` | §3.1 | a name is used above its declaration | `bc tip: …` above `let tip` | move the declaration up |
| `duplicate_name` | §4 | two declarations have one name, also a name from the standard library, a `use` file or Appendix A | two `let sigma`, or `let t` (time) | rename one |
| `use_item` | §4.1 | a `use` file does not exist, or holds a header, a `use` line or an item other than `material`, `param`, `let` or `fn` | a `physics` in a library file | move the item to the model |
| `unknown_unit` | §6.3 | a unit name does not exist | `5 mms` | the nearest unit name |
| `table_row_width` | §5.6 | a table row has the wrong number of values | 3 values for 2 columns | the column list |
| `header_order` | §3 | the header line is missing or not first | | `model name: 3d, SI` |

### 2.2 Types and units

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `type_mismatch` | §5.7 | the type of the value is not the declared type | `let count: Int = 2.5` | the type of the value |
| `unit_mismatch` | §6.1 | a sum, comparison or equation joins two dimensions | `E + rho` | the two dimensions, and the unit that would make them equal |
| `rank_mismatch` | §7.2 | a sum or an equation joins two ranks | `div(sigma) = 0 Pa` | the rank of each side |
| `ambiguous_product` | §7.2 | `*` between two vectors | `v*n` | `dot(v, n)` |
| `unit_fractional` | §6.1 | `sqrt` or `^` gives an exponent that is not a multiple of 1/2 | `sqrt(sqrt(E))` | the exponents |
| `unit_transcendental` | §6.1 | `exp`, `ln` or a trigonometric function other than `atan2` takes a quantity with a dimension | `exp(T)` | divide by a reference value |
| `unit_affine` | §6.3 | `degC` in a compound unit | `1.2 W/(m*degC)` | `W/(m*K)` |
| `sum_no_column` | §5.6 | the argument of `sum` uses no table column | `sum(2 m)` | `sum(t.col)` |
| `sum_mixed_tables` | §5.6 | the argument of `sum` uses the columns of more than one table | `sum(a.x*b.y)` | `sum(a.x)*sum(b.y)` |
| `not_constant` | §7.6 | a varying value where a constant must be | `size: 1 mm*(1 + x/L)` | `refine` for local size |

### 2.3 Geometry and regions

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `geometry_dim` | §5.3 | the geometry type does not match the dimension | `let s: Solid = rect(…)` | the type that matches |
| `geometry_space` | §8.1 | a 2D function in a `3d` model, or the reverse | `circle` in 3D | `cylinder` or `sphere` |
| `empty_region` | §8.2 | a selection has no points | `boundary(beam, x == 300 mm)` with `L = 200 mm` | the bounding box of the geometry |
| `region_outside` | §8.4 | a region is not inside the region of its field | `bc` on a face of another body | the region of the field |
| `crack_tip_form` | §9.1 | `crack_tip` is not in the mesh of a `2d` model, its region is not a `Point`, its size is not a `Length`, or the displacement on the mesh does not have order 2 | `crack_tip front: 0.5 mm` with a `Curve` | a `Point` at the tip |

### 2.4 Fields, materials and physics

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `template_eq` | §4.2 | a physics that uses a template adds an `eq` | `physics s on b: linear_elastic solves u: … { eq … }` | write the physics in full, from the source of the template |
| `eq_count` | §10.3 | the number of `eq` lines is not the number of unknowns | `solves v, p` with one `eq` | the missing equation's unknown |
| `global_eq` | §10.5 | the `eq` of a global unknown uses a field outside `mean` and `integral` | `eq mean(sigma.zz, b) + sigma.zz = S0` | the field to put in `mean(f, r)` or `integral(f, r)`. A field PDE with a mean term is not valid |
| `material_cover` | §10.3 | a point of the region has no material, or two | | the uncovered region |
| `material_prop_missing` | §4 | a physics uses a property that one of its materials lacks | `rho` used, `steel` has no `rho` | the material to change |
| `material_prop_type` | §4 | two materials give one property with different types | | the two types |

### 2.5 Assumptions

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `assumption_unknown` | §11.1 | an assumption name does not exist | `assume plain_stress` | `plane_stress` |
| `assumption_conflict` | §11.1 | two assumptions conflict | `small_strain, finite_strain` | remove one |
| `assumption_invalid` | §11.1 | the condition of an assumption fails | `newtonian` with a stress nonlinear in the rate of strain | the term that breaks it |
| `assumption_space` | §11.1 | an assumption needs another space | `plane_stress` in a `3d` model | a `2d` model |

### 2.6 Boundary conditions

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `bc_form` | §12.1 | the left side is not an unknown or a flux | `bc tip: sigma = 0` | `sigma*n = 0` |
| `bc_conflict` | §12.3 | two conditions fight on one region and component | `u = 0` and `sigma*n = …` on one face | the two lines |
| `rigid_form` | §12.1 | `rigid(p)` is not the whole right side of an essential condition on a displacement, `p` is not a `Point`, or a `moment` is not on the point of an earlier `rigid` | `bc f: u.x = rigid(c)` | `bc f: u = rigid(c)` |
| `bc_off_boundary` | §12.3 | a `bc` region is not on the boundary of its field's region | an interior surface | `boundary(…)` |

### 2.7 Methods and steps

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `method_missing` | §13 | a physics has no `method`, and the default cannot choose one mesh | two meshes on the region | a `method` item that names the mesh |
| `method_setting` | §13 | a setting does not belong to the method kind, or a `family` or `quadrature` word is not a value of the setting | `points:` in `cg` | the valid settings or values |
| `family_cell` | §13 | `serendipity` without order 2, `incompatible` without order 1 or with `quadrature: reduced`, or either family on cells that are not `quad` or `hex` | `family: incompatible` with `order: 2` | `order: 1`, or `cell: hex` in the mesh |
| `method_physics` | §13 | the method cannot discretize the physics | `gfd` with a `plastic` block | a method that can |
| `inf_sup` | §13 | an unstable velocity-pressure pair with no stabilization | `order: { v: 1, p: 1 }`, `stabilize: none` | `stabilize: supg_pspg` or `order: { v: 2, p: 1 }` |
| `order_too_low` | §13 | the cell order is too low for the derivatives | `order: 1` with `lap(lap(w))` | the minimum order |
| `procedure_args` | §14.2 | a procedure argument is missing or not valid | `transient(end: 1 s)` with no `dt` | the argument list |
| `birth_not_meshed` | §14.3 | a `birth` region is not a part of the mesh region | | the mesh region |
| `allow_reason` | §19 | `allow` has an empty reason or an unknown id | | the list of ids |

### 2.8 Results and expectations

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `results_missing` | §15.1 | the file or the step of a `results` does not exist | | the steps in that file |
| `submodel_outside` | §15.2 | the cut region is not inside the global field region | | the bounding box of the global region |
| `expect_form` | §15.3 | the two sides of an `expect` are not comparable | a field and a point value | `at` |
| `eigenvalue_step` | §14.5 | `eigenvalue(k)` is not in a `report` or an `expect`, has no `modal` or `buckling` step, or `k` is outside `1..count` of the step | `report eigenvalue(1)` in a `static` step, `report eigenvalue(13)` in `modal(count: 12)` | a `modal` or `buckling` step, or a `k` in `1..count` |
| `when_step` | §14.5 | `when(f, value)` is not in a `report` of a `transient` step, or in an `expect` whose last step is `transient` | `report when(front, 1 m)` in a `static` step | a `transient` step |

### 2.9 Parameters, scales, studies and manufacture

| Id | Spec | Rule | Example | Suggested fix |
|---|---|---|---|---|
| `param_range` | §20.1 | a default, `--set` or study value is outside the range of the parameter | `param Re: Real = 2000 in 1..1000` | the range |
| `study_param` | §20.2 | a study varies a name that is not a `param` | `D: [0.1, 0.2] m` with `let D` | change `let` to `param` |
| `study_objective` | §20.2 | `minimize` names no scalar `report` of the last step | `minimize(drag)` with no `report drag = …` | the named reports |
| `scales_dependent` | §21.2 | two scales have dependent dimensions | `scales D, a` (both `Length`) | remove one |
| `scales_incomplete` | §21.2 | the scales do not make a dimension that the model uses | `scales D, U` in a model with a density | the missing base dimension |
| `scales_missing` | §21.2 | a model with dimensions has no `scales` line | | a `scales` line from the sizes of the model, for example `scales L, E` |
| `symmetry_map` | §12.6 | the transformation of a `symmetry` item does not map `low` onto `high` | `n: 24` on a sector of 30 degrees | the angle of the sector |
| `nondim_unit` | §6.5 | a literal has a unit in a `nondim` model | `model m: 2d, nondim` and `1 m` | remove the unit |
| `process_setting` | §22.2 | a setting does not belong to the process | `draft:` in `cnc_mill` | the valid settings |

---

## 3. Lint rules

`skein lint` gives these. Each one is a warning. The program still runs.

| Id | Rule | Why it matters | Suggested fix |
|---|---|---|---|
| `unused_let` | a `let` has no use | dead text costs review time and tokens | remove it |
| `unused_material` | a material is not assigned | the user probably forgot to assign it | assign it or remove it |
| `unused_field` | an unknown is solved but never output, reported or used | the solve costs time for nothing | output it or remove it |
| `load_no_cells` | a `bc` names a region with no active cells in that step | the load does nothing, for example on a dead region | move the `bc` to a later step |
| `bending_lock` | order-1 triangles or tetrahedra in a solid with a thickness under 4 cells in a bending load | linear simplices are too stiff in bending | `order: 2` or a finer mesh |
| `volumetric_lock` | an order-1 displacement method with `nu` > 0.49 and no pressure unknown | the solution locks and is too stiff | `order: 2` or a mixed form |
| `dt_unstable` | an explicit scheme with `dt` above the estimated stable limit (CFL for flow, the highest frequency for solids) | the run diverges | the estimated limit |
| `stabilization_missing` | `stabilize: none` where the cell Peclet or Reynolds number is above 1 | the solution oscillates | `stabilize: supg` or a finer mesh |
| `mesh_feature_small` | a feature of the geometry is smaller than the local mesh size | the mesh does not resolve it | a `refine` with the feature size |

---

## 4. Verification checks

`skein run` and `skein verify` give these. Each check prints its id, `pass`
or `fail`, the measured value and the limit.

### 4.1 Before the solve

The type and unit rules, the completeness of the boundary conditions and the
assumption rules are errors of `skein check` (§2). These checks need the
mesh or the assembled operator:

| Id | Test | Limit |
|---|---|---|
| `rigid_body` | each rigid-body mode of each solid region (3 in 2D, 6 in 3D) is removed by the essential conditions or by a coupling | no free mode |
| `pressure_level` | an incompressible flow with essential velocity on all of the boundary fixes the pressure at a point or by its mean | one reference |
| `mesh_quality` | the worst cells of each mesh: the minimum and the maximum dihedral angle in 3D, the minimum angle in 2D | 3D: 5° to 175°. 2D: at least 10°. |
| `points_quality` | the fill distance divided by the separation distance of each point cloud | at most 4 |
| `material_admissible` | `E` > 0, −1 < `nu` < 0.5, `rho` > 0, `mu` > 0, `k` > 0, where the material gives them | the physical range |

### 4.2a Design for manufacture

`skein verify --dfm` runs these checks on each `manufacture` region (spec
§22.3). They need the geometry only, so they run with no solve, and on each
case of a study.

| Id | Test | Limit |
|---|---|---|
| `dfm_wall` | local wall thickness | `min_wall` to `max_wall` |
| `dfm_draft` | draft angle of each face along the pull direction | at least `draft` |
| `dfm_undercut` | each surface is visible from a pull or approach direction | no hidden surface |
| `dfm_inner_radius` | radius of each concave edge | at least `tool_radius` |
| `dfm_hole` | hole diameter and depth ratio | `min_hole`, `max_depth_ratio` |
| `dfm_overhang` | angle of each downward face to the build direction | at most `max_overhang` |
| `dfm_bend` | bend radius and flange length | `min_bend_radius`, `min_flange` |
| `dfm_stack` | the range of a `tolerance` stack-up by its `method` | inside the limits of the item |

### 4.2 After the solve

| Id | Test | Limit |
|---|---|---|
| `residual` | the relative residual of the final linear or nonlinear iteration | at most 10 × `rtol` |
| `newton_converged` | Newton met `rtol` in each increment, with no cut of the increment left | converged |
| `equilibrium` | the sum of the reactions and the applied loads, for each static step of a solid | 10⁻⁶ × the largest applied load |
| `energy_balance` | linear static: the external work equals 2 × the strain energy. Transient: kinetic energy + strain energy − external work + dissipation stays constant. | 10⁻⁶ relative (static), 10⁻³ relative (transient) |
| `conservation` | the net mass flux through the boundary of an incompressible flow, divided by the total flux | `fv`: 10⁻⁸. `cg`, `dg`: 10⁻³. |
| `submodel_cut` | the resultant force on the cut of a submodel against the global stress resultant on the same surface | 2 % |
| `expect` | each `expect` in the model (§15.3 of the spec) | the `within` value |

### 4.3 On request

`skein verify convergence model.skn --levels 3` solves the model on
`levels` meshes. Each mesh halves the size of the one before. The tool
reports the observed order of convergence:

- of each `expect` that compares two field expressions, from the L2 error;
- of each `report` value, from the Richardson estimate of three levels.

It compares the observed order with the expected order of the method (for
example 3 in L2 for `cg` with `order: 2` and a smooth solution) and fails
when the observed order is lower by more than 0.3. The id is `convergence`.

---

## 5. Run errors

These stop a run. `allow` cannot turn them off.

| Id | Cause | Suggested fix |
|---|---|---|
| `probe_outside` | an `at` point is outside the region of the field | the bounding box of the region |
| `restart_stale` | `--from` with a source that changed before that step | run from the first step |
| `mesh_failed` | the mesher could not mesh a region | the smallest feature and the size that it needs |
| `solver_diverged` | the linear or the nonlinear solver did not converge | a smaller increment, `ramp`, or another preconditioner |
| `out_of_memory` | the model does not fit the memory of the target | a coarser mesh, or the native target |
| `case_failed` | one case of a study stopped. The study continues and the table gives the case this status with the message. | the message of the case |
