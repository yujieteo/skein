# Femahoy coverage

[Femahoy](https://github.com/yujieteo/femahoy) is the requirements corpus
of Skein. Each of its validation cases, and the fixture
`fixtures/results/cantilever-2step`, is a Skein model in
[examples/femahoy](../examples/femahoy/). Line 1 of each model names the
femahoy input that it comes from. These models and the
[samples](../samples/README.md) are the acceptance tests of phase 1 of the
[roadmap](roadmap.md). The phase 1 gate: each model passes `skein check`,
or line 2 of the model is a mark that names the open gaps that block it, in
the form `// blocked-by: G10 G18`. Each later phase removes the marks of
the gaps that it closes.

This page maps each femahoy feature to the Skein construct that states it,
and to the section of the [specification](skein_spec.md). A femahoy feature
that the merged specification cannot state is a **spec gap**. The language
does not change for a gap without a decision of the captain.

Femahoy paths on this page are relative to the root of the femahoy
repository.

## Spec gaps

Each gap has the femahoy evidence and the smallest spec change that closes
it. The models mark each gap with a comment `// gap Gn:` at the place
where the feature is missing.

The captain decided the phase of each gap on 2026-10-30. The Phase column
gives the roadmap phase that closes the gap:

- Phase 1 closes G6, G7, G9, G12, G16, G29 and G30. A separate change to
  the specification states them. This page does not change the
  specification.
- Each other gap closes in the phase that builds its solver or its physics.
  A gap that blocks a case closes when its phase adds the model of the case.
- The shell, beam and laminate lane closes G3 before phase 5. Spec
  [§4.3](skein_spec.md#43-shells-beams-and-laminates) states the shell template, and `src/shells/` solves it. The new model
  [scordelis_lo_roof.skn](../examples/femahoy/scordelis_lo_roof.skn) and
  the thin plate of
  [laminate_pagano.skn](../examples/femahoy/laminate_pagano.skn) use it.

The `// blocked-by:` mark of a model names the gaps of the model that stay
open after phase 1.

### Gaps that block a case

Skein cannot state the core physics of these cases, so they have no model.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G1 | Topology optimization: a density for each cell, SIMP, a sensitivity filter, an optimality criteria update, a volume limit. Then a solid from the density at 0.5. | `validation/cases/topology-mbb-beam/derivation.md:24`, `validation/cases/topology-mbb-beam-cad/derivation.md:23` | `minimize` (§20.2) varies only a few `param` values, and a `param` is a constant (§20.1). Geometry comes only from CSG of constants (§8). | A study procedure `topology(objective, region, volume, penalty, filter)`, and a geometry function `level(f, c)`: the region where a `Results` field `f` is `c` or more. | 5 |
| G2 | Random vibration: a modal damping ratio, a PSD input and the RMS response over a frequency band. | `validation/cases/random-vibration-cantilever/job.inp:14`, `:25`, `check.toml:9` | `harmonic` (§14.2) has no damping ratio for each mode, and no report function integrates over the frequencies of a step (§14.5). Also needs G8 and G9. | A `harmonic` argument `damping: Real`, and a report function `rms(f, psd)` over the frequencies of the step. | 4 |
| G3 | Shells (`S8R`, `S6`), also a composite shell section with plies. | `validation/cases/scordelis-lo-roof/job.inp:26`, `validation/cases/laminate-pagano/job.inp:16` | §4.2 names `reissner_mindlin(t)`, but gives no source, no types of `w` and `theta`, and no load argument. Its parameter `t` is also the built-in time. The thin plate of `laminate-pagano` is missing from its model. | Give the template source in §4.2, with a parameter `h: Length`, a load argument, typed unknowns and an optional `layup` table of plies (thickness, angle, material). | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) |
| G4 | A perfect electric conductor on a curved wall: only the tangential part of the field is zero. | `validation/cases/grid-coax-tem/coax/grid.toml:1` | The left side of an essential `bc` is an unknown or a component of it (§12.1). A component is along x, y or z, not along the wall. | Accept `cross(n, u)` as the left side of an essential condition. | 5 |
| G5 | A least-squares curve fit to 40 points. | `validation/cases/fit-nist-pontius/pontius/fit.toml:4` | A table feeds only `interp` and the rows of an `expect` (§5.6, §15.3). No construct sums an expression over the rows. `minimize` has no gradients, so it cannot get to a relative error of 1e-10. | A function `sum(expr)` over the rows of a table. A `nondim` model can then state the normal equations. Or: state that a curve fit is not a Skein model. | 2 |

### Gaps that make a model fail `check`

The model states the physics, but the merged spec rejects the construct at
the `// gap` comment. The model fails `check` until the spec change of
phase 1 closes the gap.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G9 | The frequencies of a harmonic step: one frequency, many frequencies, or a band of N points. A model with one frequency fails `check`. A case with many frequencies or a band also loses the checks at the other frequencies. | `validation/cases/line-wire-skin/wire/grid.toml:27`, `validation/cases/line-bars-getdp/bars/grid.toml:35`, `validation/cases/random-vibration-cantilever/job.inp:24` | `harmonic(freq:)` needs a list (§14.2), but a procedure argument is an `expr` (§3), and `expr` has no list. | Let a procedure argument be a `list` (§3), as in a study. | 1 |
| G30 | A power with an exponent that is not an integer constant: the Ogden terms with alpha = 1.3, 5.0 and -2.0, the points of a superellipse hole. | `validation/cases/hyperelastic-tube/job.inp:12`, `validation/cases/shape-hole-biaxial/L40/mesh.geo:15` | §6.1 needs an integer constant `k` in `a^k`. The models write `jac^(-al/3)` and `cos(k*dth)^(2/power)`. | In §6.1, let `k` be any constant `Real` when `a` is `Real`. | 1 |

### Gaps that remove a check

The model states the physics, but it cannot state some reference values of
`check.toml`. The model reports a related value, or it leaves a `// gap`
comment.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G6 | The frequency of mode k (`FREQ`) and the buckling factor of mode k (`BUCKLE`). | `validation/cases/frequency-cantilever/check.toml:9`, `validation/cases/buckling-column/check.toml:7`, `validation/cases/grid-free-block-modes-ccx/check.toml:13`, `validation/cases/grid-cavity-te101/check.toml:10` | `modal` and `buckling` give eigenpairs (§14.2), but no expression names an eigenvalue (§14.5, Appendix A). | A built-in `eigenvalue(k: Int)`, for `report` and `expect`: a `Frequency` after `modal`, a `Real` load factor after `buckling`. | 1 |
| G7 | The real part, the imaginary part, the magnitude and the phase of a complex result (S-parameters, impedance, R and L of a line). | `validation/cases/grid-waveguide-te10/check.toml:10`, `validation/cases/line-bars-getdp/check.toml:11`, `validation/cases/acoustics-pulsating-sphere/check.toml:17`, `validation/cases/team7-eddy-current/50hz/case.pro:89` | §5.1 defines `Complex<D>` and `im`, but Appendix A has no function that takes a part of a complex value. | Add `re`, `imag`, `arg` and `conj` to Appendix A, and accept `Complex` in `abs`. | 1 |
| G8 | The angular frequency of the step in a load or a wall impedance: a base acceleration, a surface impedance of copper. | `validation/cases/random-vibration-cantilever/job.inp:27`, `validation/cases/grid-cavity-te101/cavity/grid.toml:20` | No built-in name gives the frequency of a `harmonic` or `modal` step (§14.2, Appendix A). | A built-in constant `omega`, the angular frequency of the current `harmonic` frequency or `modal` eigenvalue. | 4 |
| G10 | A check after an earlier step, a check for each case of a parameter sweep, and a check of the optimum of a shape loop. | `validation/cases/hyperelastic-tube/check.toml:18`, `validation/cases/cavity-devahldavis/check.toml:40`, `validation/cases/naca0012-ladson/check.toml:19`, `validation/cases/shape-hole-biaxial/check.toml:23` | §15.3 runs each `expect` once, after the last step, with the default `param` values. | Let `expect` stand in a `step` (checked at the end of the step) and in a `study` block (checked for each case, and for the best case of `minimize`). | 2 |
| G11 | The time at which a value gets to a level: the front of a dam break. | `validation/cases/dam-break-vof/check.toml:19`, `validation/cases/dam-break-martin-moyce/check.toml:17` | Reports give values at the end of a step only (§14.5, §15.3). | A report function `when(expr, value)`: the first time at which `expr` gets to `value`. | 2 |
| G12 | An angle from a ratio: the shock angle. | `validation/cases/wedge-oblique-shock/check.toml:19` | Appendix A has `tan` but no inverse. The model checks `tan(beta)`. | Add `atan` and `atan2` to Appendix A. | 1 |
| G13 | The directivity and the efficiency of an antenna from a near-to-far-field transform. | `validation/cases/dipole-king-middleton/omega12.5/openems.py:65` | Report functions act on regions of the mesh only (§14.5). | A report function `far_field(E, H, surface, direction)`. | 5 |
| G14 | A stress intensity factor in MPa·mm^0.5. | `validation/cases/fracture-centre-crack/check.toml:7` | A dimension has integer exponents, and `sqrt` needs even exponents (§6.1). The model checks the shape factor `K_I/(S*sqrt(pi*a))`, with the same tolerance. | Allow exponents in halves (§6.1). | 2 |
| G15 | The effective Young's modulus of a periodic cell under a mean stress. | `validation/cases/grid-fibre-composite-cell/check.toml:31` | The mean strain must be a global unknown that a mean stress fixes, but `mean` is valid only in `output`, `report` and `expect` (§14.5). | Allow `mean(f, r)` and `integral(f, r)` in the `eq` of a global unknown (§10.5). | 2 |
| G16 | A concentrated force on a node (`*CLOAD`). | `validation/cases/topology-mbb-beam/60x20/job.inp:2517`, `fixtures/results/cantilever-2step/cantilever-2step.inp` | §8.3 says that a `Point` can be the region of a point load, but §12.1 gives no left side for a force at a point: `n` has no value there. `cantilever_2step.skn` uses a traction with the same resultant. | On a `Point` region, read a natural condition `F*n = P` as the force `P`, with the dimension of `F` times area. | 1 |
| G17 | A rigid constraint of a boundary (`*EQUATION`): the face stays straight and can turn. | `validation/cases/contact-cattaneo/mesh.inp:13915` | §12 has only essential and natural conditions. The model puts tractions with the same force and moment on the face. | A built-in `rigid(p)` for the right side of an essential `bc`. A `bc` on the `Point` `p` takes the force and the moment. | 2 |

### Gaps that change the discretization or the physics

The model states the problem, but the result can differ from the femahoy
result, or the mesh or the solve is less efficient.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G18 | Element families: reduced integration (`C3D20R`, `C3D8R`, `CPE8R`, `S8R`), serendipity cells (`C3D20`, `CPE8`) and incompatible modes (`C3D8I`, the hex8 of the grid solver). | `validation/cases/frequency-cantilever/mesh.inp:2120`, `validation/cases/dynamics-bar-step/mesh.inp:2028`, `validation/cases/grid-cantilever/derivation.md:3` | The setting `quadrature` has no values, and `family` is `lagrange`, `nedelec` or `raviart_thomas` (§13). The models use `lagrange` of order 2. | Give `quadrature` the values `full` and `reduced`, and add `serendipity` and `incompatible` to `family`. | 2 |
| G19 | A structured hex grid through layers and blocks, with a grid plane on each material face. | `validation/cases/warpage-2p5d-frame/mesh.inp:5462`, `validation/cases/grid-heated-block-ccx/block/grid.toml:16`, `validation/cases/grid-layered-wall/wall/grid.toml:18` | `hex` and `quad` need one `box` or `rect` (or, for `hex`, a quarter of a tube) with one `size` (§9.1). §9.2 does not put a grid plane on an inner material face. The models use tetrahedra, or a size that puts the planes on the faces. | Allow `hex` and `quad` on a union of axis-aligned boxes or rects, with a grid plane on each face of each box. | 2 |
| G20 | Quarter-point cells at a crack tip. | `validation/cases/fracture-centre-crack/derivation.md:38` | No mesh setting moves the mid-side nodes (§9.1). The model uses a fine mesh. | A mesh entry `crack_tip r: size`. | 2 |
| G21 | The stable increment of an explicit step, chosen by the solver. | `validation/cases/dynamics-bar-step/job.inp:14` | `transient` needs `dt` (§14.2). The model gives a `dt` below the stable limit. | Accept `dt: stable` with an explicit scheme. | 2 |
| G22 | Finite-strain J2 plasticity with a multiplicative split (plastic strain up to 3.24). | `validation/cases/taylor-bar-impact/job.inp:12` | §10.6 states `plastic` with an additive plastic strain. Under `finite_strain` this is a different law at large strain. | One rule in §10.6: under `finite_strain`, the `strain` state is `Fp`, the elastic part is `defgrad(u)*inv(Fp)`, and the return map is the exponential map. | 5 |
| G23 | The stick stiffness of penalty friction. It changes the reference traction by 0.34 % to 3.4 %. | `validation/cases/contact-cattaneo/job.inp:15` | `contact` has only `friction`, `enforce` and `stiffness` (§12.4). | A `contact` setting `stick_stiffness`. | 5 |
| G24 | Smooth closed shapes from points or parameters: the airfoil splines, a superellipse hole. | `validation/cases/naca0012-ladson/a4.04/system/blockMeshDict:67`, `validation/cases/shape-hole-biaxial/L40/mesh.geo:9` | `spline` gives a `Curve`, and no function makes a `Surface` from a closed curve (§8.1). The models use polygons. | A geometry function `fill(c: Curve)` for a closed curve. | 2 |
| G25 | The stop rule and the first simplex of Nelder-Mead. The case tolerance comes from the stop rule. | `validation/cases/shape-hole-biaxial/derivation.md:50` | `minimize` has only `objective`, `method` and `max_cases` (§20.2). | `minimize` arguments `xtol` and `step`. | 2 |
| G26 | An open boundary: a PML, or a map of a shell to infinity. | `validation/cases/team7-eddy-current/50hz/case.pro:50`, `validation/cases/dipole-king-middleton/omega12.5/openems.py:39`, `validation/cases/optics-mie-cylinder/res40/meep_case.py:38` | §8 and §13 have no absorbing layer and no infinite map. The models put a first-order absorbing `bc`, or a zero field, on a far boundary. | A method setting `pml: { region, thickness }`. | 4 |
| G27 | A tree gauge for edge elements in a region with no conductor. | `validation/cases/team7-eddy-current/50hz/case.pro:36` | `family: nedelec` has no gauge setting (§13). | A method setting `gauge: tree`. | 5 |
| G28 | A material property that is a function of the unknown of the physics that uses it: a conductivity of T. | `validation/cases/grid-radiating-fin-ccx/fin/grid.toml:9` | §10.2 allows a property of a temperature field, but §3.1 needs each name before its use, and the `material` comes before the `physics` that declares `T`. The model puts the property in a physics `let`. | In §3.1, let a material property name an unknown of a physics that assigns the material. | 2 |
| G29 | A temperature below 0 °C, for example a sink at -5 °C. | `validation/cases/grid-layered-wall/wall/grid.toml:49` | `-` is a unary operator (§2.4) and `5 degC` is 278.15 K (§6.3), so `-5 degC` is -278.15 K. The model gives the values in K. | A `-` directly before a literal with `degC` is part of the literal. | 1 |

### Gaps of the 787 study

The [787 fuselage study](studies/b787/skein/notes.md) is not a femahoy
case, but its two models use the same marks. Its evidence paths are
relative to the root of this repository. The composite skin of the models
uses the shell template of G3 (`docs/studies/b787/skein/fuselage.skn:150`, `docs/studies/b787/skein/door_corner.skn:158`).

| Gap | Study feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G31 | Offset beams on a shell: the stringers and the frames, with a torsion constant, a section axis, an offset and the axial force. | `docs/studies/b787/skein/fuselage.skn:166`, `:184`, `:204`, `docs/studies/b787/skein/door_corner.skn:175` | §4.2 gives `timoshenko(A, I_y, I_z)` with no source, no torsion constant, no section axis, no offset and no axial force. No rule ties a beam to the shell points on its line. | `timoshenko` arguments `J`, `y_axis` and `offset`, and a template `let Nx` (axial force). In §12.1, let the right side of an essential `bc` be an unknown of another physics on the same points (a tie). | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) and [§12.1](skein_spec.md#121-form) |
| G32 | A shell or a beam line cut by a solid: `skin & door_box(pad_w)`, `stringers & left`. | `docs/studies/b787/skein/fuselage.skn:61`, `:70`, `:86`, `docs/studies/b787/skein/door_corner.skn:65`, `:75`, `:78`, `:91` | §5.3 and §8.4 do not give the dimension of `a & b` when `a` is a `Surface` or a `Curve` and `b` is a `Solid`. The checker gives the higher dimension, a `Solid`. | In §8.4: `a & b` has the lower dimension of `a` and `b`, and `a - b` has the dimension of `a`. | not set |
| G33 | Shell resultants: the force resultant `N`, the membrane strain `eps_m`, an edge load `N*n`, and a strain of the global model in a submodel. | `docs/studies/b787/skein/fuselage.skn:160`, `:201`, `docs/studies/b787/skein/door_corner.skn:168`, `:203`, `:208`, `:214`, `:217` | The shell template of G3 declares no resultant. A natural `bc` needs a flux of an `eq` (§12.1). §15.1 gives a `Results` value the unknowns of the other model, not its physics `let` fields. | The G3 template gives `let N` and `let eps_m`, and its `eq` has `div(N)`, so `N*n` is a natural condition. In §15.1, a `Results` value also gives the physics `let` fields. | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) and [§15.1](skein_spec.md#151-results) |
| G34 | A submodel for each case of a parameter sweep: the door corner takes the global case with the same `p_U`. | `docs/studies/b787/compare/edits/e4_sweep.skein.patch:15` | §15.1 gives `results(file, step:)` no study case, so a `Results` value refers to one state only. | A `case:` argument of `results`, which names a `param` value of the other model. | not set |
| G35 | A value for each copy of a `pattern`: the centre of each window, for the window load. | `docs/studies/b787/compare/edits/e3_windows.skein.patch:53` | §8.1 gives `pattern(a, v, n)` the type of `a`, with no index and no value for each copy. The model computes x - x_w with `atan2`. | `pattern` gives an index `k` of the copy, and the position of each copy, for use in a `bc`. | not set |

### Points that the spec does not make clear

These are not gaps, because the models can state the features. But the
models chose one meaning, and the spec must confirm it.

- `expect` uses `mean`, `max`, `min`, `reaction` and `integral` in many
  models. §15.3 does not list the functions that an `expect` may use.
- `orthotropic(R)` changes `hooke` (§11.1), but `hooke` takes one `E` and
  one `nu`. The orthotropic models write the stiffness out in full.
- §9.2 does not say that the mesh follows the faces between the materials of
  one physics, that two touching solids of a union share nodes, or that the
  two meshes of an `interface` conform.
- `buckling` (§14.2) does not say which load the factor multiplies, or if
  the geometric stiffness needs `finite_strain`.
- `interp` (Appendix A) does not say what it gives outside the table.
- `principal` (Appendix A) has no return type.
- `static(ramp:)` (§14.2) does not say if it ramps the full load or the
  change of load from the step before.
- `within p%` of a `degC` value is a percent of the value in K. The models
  give an absolute tolerance in K.
- The literal `2.3e-6/K` does not parse as a unit, because a unit
  expression starts with a unit name (§6.3). The models write `2.3e-6 K^-1`.
- A `report` name is not in the namespace (§4). Some models use the same
  name for a `report` and an `expect`, which gives two study columns with
  one name (§20.6).
- A `bc` on a surface inside a region (the feed of the dipole) needs the
  rule of §8.3 or §12.3, and neither rule makes the case clear.
- §4.2 puts the `let` names of a template in the namespace. It does not say
  if a model that does not use the template may declare `sigma`.

## Feature coverage

One row for each femahoy feature. "Cases" names the femahoy cases that use
the feature. The model of a case has the same name, with `_` in place of
`-`. A gap in the Spec column links to the gap tables above.

### Units and geometry

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `** units: SI-mm` and `SI-m` | each case | header unit system `mm_N_MPa_s` or `SI` | [§6.5](skein_spec.md#65-unit-systems) |
| A model with no units | shape-hole-biaxial, fit-nist-pontius | header `nondim` | [§6.5](skein_spec.md#65-unit-systems) |
| A mesh in the input (`*NODE`, `*ELEMENT`, `blockMeshDict`, `.msh`) | each mesh case | the geometry again as CSG. Skein does not import meshes (decision 8). | [§8.1](skein_spec.md#81-distance-functions), [§9](skein_spec.md#9-meshes-and-point-clouds) |
| A STEP part (`bar.step`) | bar-pressure, topology-mbb-beam-cad | the same part as CSG | [§8.1](skein_spec.md#81-distance-functions) |
| Box, plate, bar, slab, layer | most cases | `box`, `rect` | [§8.1](skein_spec.md#81-distance-functions) |
| Hole, tube, rod, ball, spherical shell | grid-lame-cylinder, plasticity-sphere, contact-hertz, acoustics-pulsating-sphere | `cylinder`, `circle`, `sphere`, `-`, `&` | [§8.1](skein_spec.md#81-distance-functions), [§8.4](skein_spec.md#84-region-algebra) |
| A quarter or an eighth of a symmetric part | submodel-kirsch-hole, plasticity-sphere, warpage-2p5d-frame | `&` with a `box`, and roller `bc` lines | [§8.4](skein_spec.md#84-region-algebra), [§12.1](skein_spec.md#121-form) |
| A sloped fillet, a coil with round corners | warpage-2p5d-frame, team7-eddy-current | `&` with a turned `box`, a union of `box` and `cylinder` with `pattern` | [§8.1](skein_spec.md#81-distance-functions) |
| Polygon footprints of the grid solver | grid-kirsch-hole, grid-plate-hole-ccx | `polygon`, `box`, or the true circle | [§8.1](skein_spec.md#81-distance-functions) |
| Grid parts with `priority` over other parts | grid-fibre-composite-cell, grid-waveguide-dielectric | `-` and `material a on r1, b on r2` | [§8.4](skein_spec.md#84-region-algebra), [§10.3](skein_spec.md#103-physics) |
| Node sets, element sets and surfaces (`*NSET`, `*ELSET`, `*SURFACE`), face sets of the grid solver | each case | `boundary(g, s)` with a condition of `x`, `y`, `z`, or with a shape | [§8.2](skein_spec.md#82-regions) |
| A point (a node set of 1 node, a probe) | grid-lame-cylinder, warpage-2p5d-frame, grid-heated-block-stress-ccx | `point(p)` of type `Point` | [§8.3](skein_spec.md#83-lower-dimension-entities) |
| A curved shell surface | scordelis-lo-roof | `boundary(cylinder(...), ...)` of type `Surface` | [§8.3](skein_spec.md#83-lower-dimension-entities), [§4.3](skein_spec.md#43-shells-beams-and-laminates) |
| Spline airfoils, a superellipse hole | naca0012-ladson, naca23012-flap, shape-hole-biaxial | `polygon` through the points | gaps [G24](#gaps-that-change-the-discretization-or-the-physics) and [G30](#gaps-that-make-a-model-fail-check) |
| A solid from a density field | topology-mbb-beam-cad | none | gap [G1](#gaps-that-block-a-case) |

### Elements and meshes

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `C3D10`, and the 10-node tetrahedra of the grid solver | bar-pressure, buckling-column, contact-hertz, plasticity-sphere, heat-fin, cantilever-2step, grid-kirsch-hole, grid-lame-cylinder, grid-round-cantilever, grid-plate-hole-stress-ccx | `mesh m on r { cell: tet }`, default `cg` with `order: 2` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| Linear tetrahedra of the grid solver (heat, conduction) | grid-pipe-wall, grid-plate-hole-ccx, grid-coax-capacitor, grid-holm-spot, grid-fibre-composite-cell | `cell: tet`, `method p: cg { order: 1 }` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| `C3D8`, and the hex8 heat cells of the grid solver | birth-death-column, grid-cooling-slab, grid-heated-bar, grid-layered-wall, grid-radiating-fin-ccx, grid-radiation-slab, grid-resistor-bar, grid-joule-bar, grid-conductivity-slab | `cell: hex` on a `box`, `order: 1` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| `C3D8I`, and the hex8 with incompatible modes of the grid solver | grid-cantilever, grid-bimetal-strip, grid-frequency-cantilever, grid-free-block-modes-ccx, turek-hron-fsi1 | `cell: hex` with `order: 2` in its place | gap [G18](#gaps-that-change-the-discretization-or-the-physics) |
| `C3D20`, `C3D20R`, `C3D8R`, `CPE8R`, `S8R` (serendipity, reduced integration) | frequency-cantilever, laminate-pagano, dynamics-bar-step, taylor-bar-impact, submodel-kirsch-hole, hyperelastic-tube | Lagrange cells of order 1 or 2 | gap [G18](#gaps-that-change-the-discretization-or-the-physics) |
| A hex grid through layers or a union of blocks | warpage-2p5d-frame, warpage-thick-copper, grid-heated-block-ccx, grid-layered-wall, grid-laminate-cell | `cell: tet` with `refine`, or one `box` | gap [G19](#gaps-that-change-the-discretization-or-the-physics) |
| Plane strain elements (`CPE4`, `CPE6`, `CPE8`, `CPE8R`) | submodel-kirsch-hole, fracture-centre-crack, contact-cattaneo, hyperelastic-tube, topology-mbb-beam | header `2d`, `cell: tri`, `assume plane_strain(th)` | [§3](skein_spec.md#3-grammar), [§9.1](skein_spec.md#91-mesh), [§11.1](skein_spec.md#111-assumption-table) |
| Plane stress elements (`CPS6`) | shape-hole-biaxial | `assume plane_stress(th)` | [§11.1](skein_spec.md#111-assumption-table), [§11.2](skein_spec.md#112-how-plane_stress-reduces-a-3d-law) |
| Quarter-point cells at a crack tip | fracture-centre-crack | `refine tip: size` | gap [G20](#gaps-that-change-the-discretization-or-the-physics) |
| Shell elements (`S8R`, `S6`), a composite shell section | scordelis-lo-roof, laminate-pagano | the `reissner_mindlin` template on a `Surface`, with a `layup` table for a composite section | [§4.3](skein_spec.md#43-shells-beams-and-laminates), [§8.3](skein_spec.md#83-lower-dimension-entities) |
| First-order edge elements | grid-waveguide-te10, grid-waveguide-dielectric, grid-waveguide-step, grid-cavity-te101, team7-eddy-current, dipole-king-middleton | `method p: cg { order: 1, family: nedelec }` | [§13](skein_spec.md#13-methods) |
| Triangles of the line solver | line-coax, line-stripline, line-wire-skin, line-bars-getdp | header `2d`, `cell: tri`, `order: 1` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| OpenFOAM cells and schemes (`Gauss linear`, `linearUpwind`, `vanLeer`, Rhie-Chow) | cavity-devahldavis, cht-heated-channel, dam-break-vof, naca0012-ladson, naca23012-flap, wedge-oblique-shock, turek-hron-fsi1 | `method p: fv { gradient, flux, limiter, pv_coupling }` | [§13](skein_spec.md#13-methods), [§17.3](skein_spec.md#173-finite-volume-fv) |
| Boundary layer cells at a wall | naca0012-ladson, naca23012-flap | `layers wall: { first, ratio, count }` | [§9.1](skein_spec.md#91-mesh), [§9.2](skein_spec.md#92-mesh-generation) |
| Mesh size, local size and growth (`max_cell`, `max_growth`, local fields of gmsh) | most cases | `size`, `refine r: size`, `grade` | [§9.1](skein_spec.md#91-mesh) |
| A count of cells through a layer (`min_cells`) | grid-bimetal-strip, grid-layered-wall | `size` or `refine` that gives the count | [§9.1](skein_spec.md#91-mesh) |
| SPH particles | dam-break-martin-moyce | `points` and `method p: rkpm` in a Lagrangian frame | [§9.3](skein_spec.md#93-points), [§13](skein_spec.md#13-methods), [§17.5](skein_spec.md#175-meshless) |
| FDTD grids of openEMS and Meep | dipole-king-middleton, optics-mie-cylinder | a frequency-domain `cg` model for each frequency | [§13](skein_spec.md#13-methods), [§14.2](skein_spec.md#142-procedures) |
| Rigid bodies and joints (Chrono) | slider-crank | global unknowns and `method p: ode` | [§10.5](skein_spec.md#105-odes-and-global-unknowns), [§13](skein_spec.md#13-methods) |

### Materials

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*ELASTIC` (isotropic), `*SOLID SECTION` | most structural cases | `material m { E: Pressure, nu: Real }`, `physics ... : linear_elastic` or `hooke(E, nu, strain(u))` | [§10.2](skein_spec.md#102-material), [§4.2](skein_spec.md#42-the-standard-library), Appendix [A](skein_spec.md#appendix-a-built-in-names) |
| `*DENSITY` | frequency-cantilever, dynamics-bar-step, birth-death-column, taylor-bar-impact | `rho: Density` and `rho*ddt2(u)` | [§10.2](skein_spec.md#102-material), [§7.3](skein_spec.md#73-differential-operators) |
| `*ELASTIC, TYPE=ENGINEERING CONSTANTS`, `*ORIENTATION`, orthotropic grid materials | laminate-pagano, grid-constrained-block | the orthotropic law in full with `outer`, one material for each ply | [§7.2](skein_spec.md#72-arithmetic-on-tensors), [§10.2](skein_spec.md#102-material) |
| `*EXPANSION` with a reference temperature | warpage-2p5d-frame, warpage-thick-copper, grid-bimetal-strip, grid-constrained-block, grid-heated-block-stress-ccx, grid-plate-hole-stress-ccx | `alpha: Expansion` and `strain(u) - alpha*(T - T0)*I` | [§10.4](skein_spec.md#104-coupling) |
| `*PLASTIC` (J2, linear hardening) | plasticity-sphere, taylor-bar-impact | `state`, `plastic { yield: von_mises(sigma) - ..., strain, hardening }` | [§10.6](skein_spec.md#106-inelastic-materials), gap [G22](#gaps-that-change-the-discretization-or-the-physics) |
| `*HYPERELASTIC` (Ogden) | hyperelastic-tube | a strain energy `W` and `S = 2*diff(W, C)` | [§10.6](skein_spec.md#106-inelastic-materials), gap [G30](#gaps-that-make-a-model-fail-check) |
| St. Venant-Kirchhoff under `NLGEOM` | turek-hron-fsi1, buckling-column | `assume finite_strain` and `hooke` of the Green-Lagrange strain | [§11.1](skein_spec.md#111-assumption-table) |
| `*CONDUCTIVITY`, `*SPECIFIC HEAT`, grid heat materials | heat-fin, grid-cooling-slab and the other heat cases | `k: Conductivity`, `cp: SpecificHeat`, `q = fourier(k, T)` | [§6.2](skein_spec.md#62-named-dimension-types), [§10.2](skein_spec.md#102-material) |
| A conductivity table of the temperature | grid-radiating-fin-ccx, grid-conductivity-slab | `table` and `interp(..., T)` in a physics `let` | [§5.6](skein_spec.md#56-tables), Appendix [A](skein_spec.md#appendix-a-built-in-names), gap [G28](#gaps-that-change-the-discretization-or-the-physics) |
| Electric conductivity, permittivity, permeability, loss | grid-resistor-bar, grid-coax-capacitor, grid-waveguide-dielectric, line-*, team7-eddy-current | properties with dimension expressions, for example `A/(V*m)` | [§6.2](skein_spec.md#62-named-dimension-types), [§6.3](skein_spec.md#63-units), [§10.2](skein_spec.md#102-material) |
| Viscosity, density, Prandtl number of a fluid | cavity-devahldavis, cht-heated-channel, naca0012-ladson, turek-hron-fsi1 | `mu: DynViscosity`, `rho: Density`, `k`, `cp` | [§10.2](skein_spec.md#102-material) |
| Boussinesq buoyancy | cavity-devahldavis | the body force `rho*(1 - beta*(T - T0))*g` in the momentum equation | [§10.3](skein_spec.md#103-physics), [§11.1](skein_spec.md#111-assumption-table) |
| Perfect gas | wedge-oblique-shock | `gamma: Real` and a gas constant `let` | [§10.2](skein_spec.md#102-material) |
| Two phases with a volume fraction (VOF) | dam-break-vof | mixture laws of `alpha` in physics `let` items | [§10.3](skein_spec.md#103-physics) |
| Tait equation of state | dam-break-martin-moyce | `p = B*((1/J)^7 - 1)` | [§10.3](skein_spec.md#103-physics) |
| A material on each element set or grid part | most multi-material cases | `material a on r1, b on r2` | [§10.3](skein_spec.md#103-physics) |
| A homogenized material of a periodic cell | grid-laminate-cell, grid-fibre-composite-cell | one physics for each load case, `symmetry s: bloch { k: 0 }`, `mean` in `expect` | [§12.6](skein_spec.md#126-symmetry-groups), [§14.5](skein_spec.md#145-output-and-report), gap [G15](#gaps-that-remove-a-check) |
| Rigid body mass and inertia | slider-crank | `let` constants in the kinetic energy | [§10.5](skein_spec.md#105-odes-and-global-unknowns) |

### Steps

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*STATIC` | most structural cases | `step s: static` | [§14.2](skein_spec.md#142-procedures) |
| `*STATIC` with increments and `NLGEOM` | hyperelastic-tube, plasticity-sphere, contact-cattaneo | `static(ramp: n)` and `assume finite_strain` | [§14.2](skein_spec.md#142-procedures), [§11.1](skein_spec.md#111-assumption-table) |
| Two or more steps, with the boundary conditions kept | cantilever-2step, contact-hertz, contact-cattaneo, hyperelastic-tube | steps in the order of the file, `release r` | [§14.1](skein_spec.md#141-form), [§12.2](skein_spec.md#122-where-bc-may-stand) |
| `*FREQUENCY`, a modes step of the grid solver | frequency-cantilever, grid-frequency-cantilever, grid-free-block-modes-ccx | `modal(count: n)` | [§14.2](skein_spec.md#142-procedures), gap [G6](#gaps-that-remove-a-check) |
| Modes near a frequency | grid-cavity-te101 | `modal(count: 1, shift: f)` | [§14.2](skein_spec.md#142-procedures) |
| `*BUCKLE` after a preload | buckling-column | `static` then `buckling(count: n)` | [§14.2](skein_spec.md#142-procedures) |
| `*DYNAMIC, EXPLICIT` | dynamics-bar-step, taylor-bar-impact | `transient(end, dt, scheme: central_difference)` or `ssprk3` | [§14.2](skein_spec.md#142-procedures), gap [G21](#gaps-that-change-the-discretization-or-the-physics) |
| `*HEAT TRANSFER, STEADY STATE`, a steady grid heat step | heat-fin and the grid heat cases | `static` | [§14.2](skein_spec.md#142-procedures) |
| A transient heat step (backward Euler) | grid-cooling-slab | `transient(end, dt, scheme: backward_euler)` | [§14.2](skein_spec.md#142-procedures) |
| `*MODEL CHANGE, REMOVE` and `ADD=STRAIN FREE` | birth-death-column | `death r` and `birth r` | [§14.3](skein_spec.md#143-birth-and-death) |
| `*STEADY STATE DYNAMICS` with modal damping | random-vibration-cantilever | none | gap [G2](#gaps-that-block-a-case) |
| One frequency of an electromagnetic or acoustic problem | grid-waveguide-*, acoustics-pulsating-sphere, dipole-king-middleton, optics-mie-cylinder, team7-eddy-current, line-wire-skin, line-bars-getdp | `harmonic(freq: [f])` with `Complex` unknowns | [§5.1](skein_spec.md#51-scalar-types), [§14.2](skein_spec.md#142-procedures), gap [G9](#gaps-that-make-a-model-fail-check) |
| Many frequencies in one step | line-wire-skin, line-bars-getdp, random-vibration-cantilever | none | gap [G9](#gaps-that-make-a-model-fail-check) |
| DC conduction and electrostatic steps | grid-resistor-bar, grid-holm-spot, grid-coax-capacitor, line-coax, line-stripline | a potential physics `eq div(J) = 0` in a `static` step | [§10.3](skein_spec.md#103-physics), [§14.2](skein_spec.md#142-procedures) |
| Joule heat as the source of a heat step | grid-joule-bar | 2 coupled physics | [§10.4](skein_spec.md#104-coupling) |
| A thermal stress step after a steady heat step | grid-heated-block-stress-ccx | 2 physics, `couple: staggered` | [§10.4](skein_spec.md#104-coupling), [§13](skein_spec.md#13-methods) |
| A steady CFD solve (SIMPLE) | cavity-devahldavis, cht-heated-channel, naca0012-ladson, naca23012-flap, turek-hron-fsi1 | `static` | [§14.2](skein_spec.md#142-procedures) |
| A transient CFD solve | dam-break-vof, wedge-oblique-shock | `transient(end, dt, scheme)` | [§14.2](skein_spec.md#142-procedures) |
| An implicit time scheme of a mechanism (HHT) | slider-crank | `transient(..., scheme: generalized_alpha(rho_inf))` | [§14.2](skein_spec.md#142-procedures) |
| Solver settings (Newton tolerance) | grid-radiation-slab, grid-radiating-fin-ccx | `solve { nonlinear: newton, rtol }` | [§14.4](skein_spec.md#144-solve) |
| A submodel: the global model runs first | submodel-kirsch-hole | `let global: Results = results(file, step:)` | [§15.1](skein_spec.md#151-results), [§15.2](skein_spec.md#152-submodels) |
| Fluid-structure coupling (preCICE) | turek-hron-fsi1 | `interface fsi: wet_fluid, wet_solid { mesh_motion: harmonic }` | [§10.7](skein_spec.md#107-fsi-interfaces) |
| Conjugate heat transfer of 2 regions | cht-heated-channel | 2 physics and `bc` lines on the shared boundary | [§10.4](skein_spec.md#104-coupling), [§12.1](skein_spec.md#121-form) |
| Case folders for each parameter value | cavity-devahldavis, naca0012-ladson, naca23012-flap, dipole-king-middleton, optics-mie-cylinder | `param` and `study s: grid` | [§20.1](skein_spec.md#201-param), [§20.2](skein_spec.md#202-study), gap [G10](#gaps-that-remove-a-check) |
| Shape optimization (Nelder-Mead) | shape-hole-biaxial | `study s: minimize(objective, method: nelder_mead, max_cases)` | [§20.2](skein_spec.md#202-study), gap [G25](#gaps-that-change-the-discretization-or-the-physics) |
| Topology optimization (SIMP) | topology-mbb-beam, topology-mbb-beam-cad | none | gap [G1](#gaps-that-block-a-case) |
| A least-squares curve fit | fit-nist-pontius | none | gap [G5](#gaps-that-block-a-case) |

### Loads

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*DSLOAD` and `*DLOAD` pressure | bar-pressure, buckling-column, plasticity-sphere, dynamics-bar-step, fracture-centre-crack | `bc r: sigma*n = -p*n` | [§12.1](skein_spec.md#121-form) |
| A traction vector on a face | grid-cantilever, grid-kirsch-hole, grid-round-cantilever | `bc r: sigma*n = (...)` | [§12.1](skein_spec.md#121-form) |
| A traction that varies along a face | laminate-pagano, contact-cattaneo | `sigma*n` with a function of `x`, `y` | [§7.6](skein_spec.md#76-coordinates-and-time), [§12.1](skein_spec.md#121-form) |
| `*CLOAD` at a node | cantilever-2step, topology-mbb-beam, contact-cattaneo | a traction with the same resultant | gap [G16](#gaps-that-remove-a-check) |
| `*DLOAD, GRAV` | birth-death-column, scordelis-lo-roof | the template argument `b: rho*(0, 0, -g)` | [§4.2](skein_spec.md#42-the-standard-library), [§4.3](skein_spec.md#43-shells-beams-and-laminates) |
| Gravity on a fluid | cavity-devahldavis, dam-break-vof, dam-break-martin-moyce | a body force term in `eq` | [§10.3](skein_spec.md#103-physics) |
| `*TEMPERATURE` and a uniform temperature of the grid solver | warpage-2p5d-frame, warpage-thick-copper, grid-bimetal-strip, grid-constrained-block | a constant `T` in the thermal strain | [§10.4](skein_spec.md#104-coupling) |
| A temperature field from a heat solve | grid-heated-block-stress-ccx | the unknown `T` of the heat physics in the thermal strain | [§10.4](skein_spec.md#104-coupling) |
| A volume heat source | grid-heated-bar, grid-heated-block-ccx, grid-plate-hole-ccx | `eq div(q) = q_v` | [§10.3](skein_spec.md#103-physics) |
| `*FILM` convection | heat-fin, grid-heated-block-ccx, grid-layered-wall, grid-pipe-wall | `bc r: dot(q, n) = h*(T - T_inf)` | [§12.1](skein_spec.md#121-form) |
| Radiation to the ambient | grid-radiating-fin-ccx, grid-radiation-slab | `bc r: dot(q, n) = eps*sb*(T^4 - T_amb^4)` | [§12.1](skein_spec.md#121-form) |
| A heat flux on a face | cht-heated-channel | `bc r: dot(q, n) = -q_h` | [§12.1](skein_spec.md#121-form) |
| `*AMPLITUDE` of a sudden load | dynamics-bar-step | a constant `bc` in a `transient` step | [§12.2](skein_spec.md#122-where-bc-may-stand), [§14.2](skein_spec.md#142-procedures) |
| `*INITIAL CONDITIONS` of temperature or velocity | heat-fin, grid-cooling-slab, taylor-bar-impact | the initial value in `solves`, a velocity unknown | [§10.1](skein_spec.md#101-unknown-fields), [§14.1](skein_spec.md#141-form) |
| A water column at the start (`setFields`) | dam-break-vof | an initial value of `alpha` from `interp` | [§10.1](skein_spec.md#101-unknown-fields), Appendix [A](skein_spec.md#appendix-a-built-in-names) |
| Voltage, current and charge on terminals | grid-resistor-bar, grid-joule-bar, grid-holm-spot, grid-coax-capacitor, line-coax | `bc r: phi = V0` | [§12.1](skein_spec.md#121-form) |
| A coil current density | team7-eddy-current | a varying `let js` | [§7.6](skein_spec.md#76-coordinates-and-time) |
| Wave ports (TE10, TEM) | grid-waveguide-te10, grid-waveguide-dielectric, grid-waveguide-step | a Robin `bc` on `Hs*n` | [§12.1](skein_spec.md#121-form) |
| A source velocity on a surface (acoustics) | acoustics-pulsating-sphere | `bc r: dot(grad(p), n) = im*omega*rho*u_s` | [§12.1](skein_spec.md#121-form) |
| A plane wave (Meep source) | optics-mie-cylinder | a source term of the scattered field | [§10.3](skein_spec.md#103-physics) |
| A delta-gap feed | dipole-king-middleton | an essential `bc` on the gap surface | [§8.3](skein_spec.md#83-lower-dimension-entities), [§12.1](skein_spec.md#121-form) |
| `*BASE MOTION` with a PSD amplitude | random-vibration-cantilever | none | gaps [G2](#gaps-that-block-a-case), [G8](#gaps-that-remove-a-check) |
| A motor at a constant speed | slider-crank | `let th = w*t` | [§7.6](skein_spec.md#76-coordinates-and-time) |

### Boundary conditions

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*BOUNDARY` fixed | most structural cases | `bc r: u = 0` | [§12.1](skein_spec.md#121-form) |
| `*BOUNDARY` on 1 DOF, symmetry planes, rollers | most structural cases | `bc r: u.x = 0` | [§12.1](skein_spec.md#121-form) |
| A support at a point | grid-lame-cylinder, warpage-2p5d-frame, grid-heated-block-stress-ccx | `bc p: u.z = 0` on a `Point` | [§8.3](skein_spec.md#83-lower-dimension-entities), [§12.1](skein_spec.md#121-form) |
| A displacement value | contact-hertz, shape-hole-biaxial | `bc r: u.z = -0.075 mm` | [§12.1](skein_spec.md#121-form) |
| A displacement that varies on a face (Kirsch field, `*TRANSFORM` radial) | submodel-kirsch-hole, hyperelastic-tube | a varying `let` or a vector expression of `x`, `y` | [§7.6](skein_spec.md#76-coordinates-and-time), [§12.1](skein_spec.md#121-form) |
| `*BOUNDARY, SUBMODEL` | submodel-kirsch-hole | `bc cut: u = global.u` | [§12.5](skein_spec.md#125-values-from-other-models), [§15.2](skein_spec.md#152-submodels) |
| A condition that a later step removes | contact-cattaneo, contact-hertz | `release r` | [§12.2](skein_spec.md#122-where-bc-may-stand) |
| A condition at the model level and in a step | frequency-cantilever, birth-death-column | a top-level `bc`, or a `bc` in a step | [§12.2](skein_spec.md#122-where-bc-may-stand) |
| A fixed temperature | heat-fin and the grid heat cases | `bc r: T = 100 degC` | [§6.3](skein_spec.md#63-units), [§12.1](skein_spec.md#121-form) |
| An adiabatic face, a free face, a magnetic wall | each case | no `bc`: the natural condition with zero right side | [§12.3](skein_spec.md#123-completeness) |
| `*EQUATION`: a face stays straight and can turn | contact-cattaneo | tractions with the same force and moment | gap [G17](#gaps-that-remove-a-check) |
| `*CONTACT PAIR`, `*SURFACE INTERACTION`, `*FRICTION` | contact-hertz, contact-cattaneo | `contact c: master, slave { friction, enforce: penalty, stiffness }` | [§12.4](skein_spec.md#124-contact), gap [G23](#gaps-that-change-the-discretization-or-the-physics) |
| A periodic cell | grid-laminate-cell, grid-fibre-composite-cell | `symmetry s: bloch { low, high, period, k: 0 }` | [§12.6](skein_spec.md#126-symmetry-groups) |
| A free body with no support | grid-free-block-modes-ccx | no `bc`, `allow rigid_body: "reason"` | [checks.md](checks.md) |
| No slip, inlet velocity, outlet pressure, slip wall, symmetry plane | cavity-devahldavis, cht-heated-channel, naca0012-ladson, wedge-oblique-shock, turek-hron-fsi1 | `bc r: v = ...`, `sigma*n = -p*n`, `dot(v, n) = 0` | [§12.1](skein_spec.md#121-form) |
| A pressure reference (`pRefCell`) | cavity-devahldavis | `bc p0: p = 0` on a `Point` | [§8.3](skein_spec.md#83-lower-dimension-entities), [§12.3](skein_spec.md#123-completeness) |
| The coupled temperature of 2 regions | cht-heated-channel | `bc r: Tf = Ts` and equal fluxes | [§12.1](skein_spec.md#121-form) |
| A perfect electric conductor on a flat wall | grid-waveguide-*, grid-cavity-te101 | the tangential components: `bc r: E.y = 0` | [§12.1](skein_spec.md#121-form) |
| A perfect electric conductor on a curved wall | grid-coax-tem | none | gap [G4](#gaps-that-block-a-case) |
| A surface impedance wall | grid-cavity-te101 | none | gap [G8](#gaps-that-remove-a-check) |
| An open boundary (PML, a map to infinity) | team7-eddy-current, dipole-king-middleton, optics-mie-cylinder | a first-order absorbing `bc` on a far boundary | gap [G26](#gaps-that-change-the-discretization-or-the-physics) |
| An exact outgoing wave condition | acoustics-pulsating-sphere | a Robin `bc` with `im` | [§5.1](skein_spec.md#51-scalar-types), [§12.1](skein_spec.md#121-form) |
| Joint constraints of a mechanism | slider-crank | `eq` lines of the loop closure on global unknowns | [§10.5](skein_spec.md#105-odes-and-global-unknowns) |

### Outputs and checks

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*NODE FILE`, `*EL FILE`, VTU output, output at each increment | each case | `output u, sigma every k` | [§14.5](skein_spec.md#145-output-and-report) |
| `*NODE PRINT` U or NT of 1 node, a probe of the grid solver | most cases | `expect e: u.z at (x, y, z) == ref within tol` | [§7.5](skein_spec.md#75-at), [§15.3](skein_spec.md#153-expect) |
| `*EL PRINT` S at an integration point | submodel-kirsch-hole, laminate-pagano | `sigma.xx at (...)` | [§7.5](skein_spec.md#75-at), [§15.3](skein_spec.md#153-expect) |
| `op = "mean"` and `op = "max"` over a set | bar-pressure, heat-fin, taylor-bar-impact, contact-hertz | `mean(f, r)`, `max(f, r)` | [§14.5](skein_spec.md#145-output-and-report), [§15.3](skein_spec.md#153-expect) |
| `RF` with `op = "sum"` | bar-pressure, hyperelastic-tube, grid-bimetal-strip | `reaction(r)` | [§14.5](skein_spec.md#145-output-and-report) |
| `CSTR` contact pressure | contact-hertz, contact-cattaneo | `max(-sigma.zz, face)` | [§14.5](skein_spec.md#145-output-and-report) |
| `FREQ` and `BUCKLE` | frequency-cantilever, buckling-column, grid-frequency-cantilever, grid-free-block-modes-ccx, grid-cavity-te101 | none | gap [G6](#gaps-that-remove-a-check) |
| `A_RMS` | random-vibration-cantilever | none | gap [G2](#gaps-that-block-a-case) |
| Heat flow through a face (`energy.txt`) | grid-pipe-wall, grid-radiation-slab | `integral(dot(q, n), r)` | [§14.5](skein_spec.md#145-output-and-report) |
| Resistance, conductance, capacitance of terminals | grid-resistor-bar, grid-holm-spot, grid-coax-capacitor | a ratio of `integral(dot(J, n), r)` and the voltage | [§14.5](skein_spec.md#145-output-and-report) |
| S-parameters, impedance, R and L of a line | grid-waveguide-*, line-*, dipole-king-middleton | projection integrals and a complex `expect` | [§14.5](skein_spec.md#145-output-and-report), [§15.3](skein_spec.md#153-expect), gap [G7](#gaps-that-remove-a-check) |
| Lift, drag and moment coefficients (`forceCoeffs`) | naca0012-ladson, naca23012-flap, turek-hron-fsi1 | `report cl = integral(dot(sigma*n, e_l), wall)/...` | [§14.5](skein_spec.md#145-output-and-report) |
| Nusselt number, peak velocities, mean temperatures | cavity-devahldavis, cht-heated-channel | `report` with `mean`, `max`, `integral` | [§14.5](skein_spec.md#145-output-and-report) |
| A curvature from a fit of the top face | grid-bimetal-strip | `integral` moments of `u.z` | [§14.5](skein_spec.md#145-output-and-report) |
| A stress intensity factor | fracture-centre-crack | a dimensionless shape factor | gap [G14](#gaps-that-remove-a-check) |
| The time at which a front gets to a position | dam-break-vof, dam-break-martin-moyce | `report` of the front position | gap [G11](#gaps-that-remove-a-check) |
| A shock angle | wedge-oblique-shock | `tan` of the angle | gap [G12](#gaps-that-remove-a-check) |
| Directivity and efficiency of an antenna | dipole-king-middleton | none | gap [G13](#gaps-that-remove-a-check) |
| A check after an earlier step | hyperelastic-tube | `report` in the step | gap [G10](#gaps-that-remove-a-check) |
| The compliance of a design | topology-mbb-beam | none | gap [G1](#gaps-that-block-a-case) |
| `rel_tol` and `abs_tol` of `check.toml` | each case | `within p%` or `within` an absolute value | [§15.3](skein_spec.md#153-expect) |

## Case inventory

The corpus is each validation case of femahoy and its result fixture. The
femahoy tests run each validation case (`fh validate`), and the viewer
tests read the fixture.

| Femahoy case | Skein model |
|---|---|
| `fixtures/results/cantilever-2step` | [cantilever_2step.skn](../examples/femahoy/cantilever_2step.skn) |
| `validation/cases/acoustics-pulsating-sphere` | [acoustics_pulsating_sphere.skn](../examples/femahoy/acoustics_pulsating_sphere.skn) |
| `validation/cases/bar-pressure` | [bar_pressure.skn](../examples/femahoy/bar_pressure.skn) |
| `validation/cases/birth-death-column` | [birth_death_column.skn](../examples/femahoy/birth_death_column.skn) |
| `validation/cases/buckling-column` | [buckling_column.skn](../examples/femahoy/buckling_column.skn) |
| `validation/cases/cavity-devahldavis` | [cavity_devahldavis.skn](../examples/femahoy/cavity_devahldavis.skn) |
| `validation/cases/cht-heated-channel` | [cht_heated_channel.skn](../examples/femahoy/cht_heated_channel.skn) |
| `validation/cases/contact-cattaneo` | [contact_cattaneo.skn](../examples/femahoy/contact_cattaneo.skn) |
| `validation/cases/contact-hertz` | [contact_hertz.skn](../examples/femahoy/contact_hertz.skn) |
| `validation/cases/dam-break-martin-moyce` | [dam_break_martin_moyce.skn](../examples/femahoy/dam_break_martin_moyce.skn) |
| `validation/cases/dam-break-vof` | [dam_break_vof.skn](../examples/femahoy/dam_break_vof.skn) |
| `validation/cases/dipole-king-middleton` | [dipole_king_middleton.skn](../examples/femahoy/dipole_king_middleton.skn) |
| `validation/cases/dynamics-bar-step` | [dynamics_bar_step.skn](../examples/femahoy/dynamics_bar_step.skn) |
| `validation/cases/fit-nist-pontius` | none: gap [G5](#gaps-that-block-a-case) |
| `validation/cases/fracture-centre-crack` | [fracture_centre_crack.skn](../examples/femahoy/fracture_centre_crack.skn) |
| `validation/cases/frequency-cantilever` | [frequency_cantilever.skn](../examples/femahoy/frequency_cantilever.skn) |
| `validation/cases/grid-bimetal-strip` | [grid_bimetal_strip.skn](../examples/femahoy/grid_bimetal_strip.skn) |
| `validation/cases/grid-cantilever` | [grid_cantilever.skn](../examples/femahoy/grid_cantilever.skn) |
| `validation/cases/grid-cavity-te101` | [grid_cavity_te101.skn](../examples/femahoy/grid_cavity_te101.skn) |
| `validation/cases/grid-coax-capacitor` | [grid_coax_capacitor.skn](../examples/femahoy/grid_coax_capacitor.skn) |
| `validation/cases/grid-coax-tem` | none: gap [G4](#gaps-that-block-a-case) |
| `validation/cases/grid-conductivity-slab` | [grid_conductivity_slab.skn](../examples/femahoy/grid_conductivity_slab.skn) |
| `validation/cases/grid-constrained-block` | [grid_constrained_block.skn](../examples/femahoy/grid_constrained_block.skn) |
| `validation/cases/grid-cooling-slab` | [grid_cooling_slab.skn](../examples/femahoy/grid_cooling_slab.skn) |
| `validation/cases/grid-fibre-composite-cell` | [grid_fibre_composite_cell.skn](../examples/femahoy/grid_fibre_composite_cell.skn) |
| `validation/cases/grid-free-block-modes-ccx` | [grid_free_block_modes_ccx.skn](../examples/femahoy/grid_free_block_modes_ccx.skn) |
| `validation/cases/grid-frequency-cantilever` | [grid_frequency_cantilever.skn](../examples/femahoy/grid_frequency_cantilever.skn) |
| `validation/cases/grid-heated-bar` | [grid_heated_bar.skn](../examples/femahoy/grid_heated_bar.skn) |
| `validation/cases/grid-heated-block-ccx` | [grid_heated_block_ccx.skn](../examples/femahoy/grid_heated_block_ccx.skn) |
| `validation/cases/grid-heated-block-stress-ccx` | [grid_heated_block_stress_ccx.skn](../examples/femahoy/grid_heated_block_stress_ccx.skn) |
| `validation/cases/grid-holm-spot` | [grid_holm_spot.skn](../examples/femahoy/grid_holm_spot.skn) |
| `validation/cases/grid-joule-bar` | [grid_joule_bar.skn](../examples/femahoy/grid_joule_bar.skn) |
| `validation/cases/grid-kirsch-hole` | [grid_kirsch_hole.skn](../examples/femahoy/grid_kirsch_hole.skn) |
| `validation/cases/grid-lame-cylinder` | [grid_lame_cylinder.skn](../examples/femahoy/grid_lame_cylinder.skn) |
| `validation/cases/grid-laminate-cell` | [grid_laminate_cell.skn](../examples/femahoy/grid_laminate_cell.skn) |
| `validation/cases/grid-layered-wall` | [grid_layered_wall.skn](../examples/femahoy/grid_layered_wall.skn) |
| `validation/cases/grid-pipe-wall` | [grid_pipe_wall.skn](../examples/femahoy/grid_pipe_wall.skn) |
| `validation/cases/grid-plate-hole-ccx` | [grid_plate_hole_ccx.skn](../examples/femahoy/grid_plate_hole_ccx.skn) |
| `validation/cases/grid-plate-hole-stress-ccx` | [grid_plate_hole_stress_ccx.skn](../examples/femahoy/grid_plate_hole_stress_ccx.skn) |
| `validation/cases/grid-radiating-fin-ccx` | [grid_radiating_fin_ccx.skn](../examples/femahoy/grid_radiating_fin_ccx.skn) |
| `validation/cases/grid-radiation-slab` | [grid_radiation_slab.skn](../examples/femahoy/grid_radiation_slab.skn) |
| `validation/cases/grid-resistor-bar` | [grid_resistor_bar.skn](../examples/femahoy/grid_resistor_bar.skn) |
| `validation/cases/grid-round-cantilever` | [grid_round_cantilever.skn](../examples/femahoy/grid_round_cantilever.skn) |
| `validation/cases/grid-waveguide-dielectric` | [grid_waveguide_dielectric.skn](../examples/femahoy/grid_waveguide_dielectric.skn) |
| `validation/cases/grid-waveguide-step` | [grid_waveguide_step.skn](../examples/femahoy/grid_waveguide_step.skn) |
| `validation/cases/grid-waveguide-te10` | [grid_waveguide_te10.skn](../examples/femahoy/grid_waveguide_te10.skn) |
| `validation/cases/heat-fin` | [heat_fin.skn](../examples/femahoy/heat_fin.skn) |
| `validation/cases/hyperelastic-tube` | [hyperelastic_tube.skn](../examples/femahoy/hyperelastic_tube.skn) |
| `validation/cases/laminate-pagano` | [laminate_pagano.skn](../examples/femahoy/laminate_pagano.skn) |
| `validation/cases/line-bars-getdp` | [line_bars_getdp.skn](../examples/femahoy/line_bars_getdp.skn) |
| `validation/cases/line-coax` | [line_coax.skn](../examples/femahoy/line_coax.skn) |
| `validation/cases/line-stripline` | [line_stripline.skn](../examples/femahoy/line_stripline.skn) |
| `validation/cases/line-wire-skin` | [line_wire_skin.skn](../examples/femahoy/line_wire_skin.skn) |
| `validation/cases/naca0012-ladson` | [naca0012_ladson.skn](../examples/femahoy/naca0012_ladson.skn) |
| `validation/cases/naca23012-flap` | [naca23012_flap.skn](../examples/femahoy/naca23012_flap.skn) |
| `validation/cases/optics-mie-cylinder` | [optics_mie_cylinder.skn](../examples/femahoy/optics_mie_cylinder.skn) |
| `validation/cases/plasticity-sphere` | [plasticity_sphere.skn](../examples/femahoy/plasticity_sphere.skn) |
| `validation/cases/random-vibration-cantilever` | none: gap [G2](#gaps-that-block-a-case) |
| `validation/cases/scordelis-lo-roof` | [scordelis_lo_roof.skn](../examples/femahoy/scordelis_lo_roof.skn) |
| `validation/cases/shape-hole-biaxial` | [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn) |
| `validation/cases/slider-crank` | [slider_crank.skn](../examples/femahoy/slider_crank.skn) |
| `validation/cases/submodel-kirsch-hole` | [submodel_kirsch_hole.skn](../examples/femahoy/submodel_kirsch_hole.skn), [submodel_kirsch_hole_global.skn](../examples/femahoy/submodel_kirsch_hole_global.skn) |
| `validation/cases/taylor-bar-impact` | [taylor_bar_impact.skn](../examples/femahoy/taylor_bar_impact.skn) |
| `validation/cases/team7-eddy-current` | [team7_eddy_current.skn](../examples/femahoy/team7_eddy_current.skn) |
| `validation/cases/topology-mbb-beam` | none: gap [G1](#gaps-that-block-a-case) |
| `validation/cases/topology-mbb-beam-cad` | none: gap [G1](#gaps-that-block-a-case) |
| `validation/cases/turek-hron-fsi1` | [turek_hron_fsi1.skn](../examples/femahoy/turek_hron_fsi1.skn) |
| `validation/cases/warpage-2p5d-frame` | [warpage_2p5d_frame.skn](../examples/femahoy/warpage_2p5d_frame.skn) |
| `validation/cases/warpage-thick-copper` | [warpage_thick_copper.skn](../examples/femahoy/warpage_thick_copper.skn) |
| `validation/cases/wedge-oblique-shock` | [wedge_oblique_shock.skn](../examples/femahoy/wedge_oblique_shock.skn) |
