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
  specification. The models use them:
  - G6: [§14.5](skein_spec.md#145-output-and-report) states `eigenvalue(k)`.
    The expectations of [frequency_cantilever.skn](../examples/femahoy/frequency_cantilever.skn),
    [grid_frequency_cantilever.skn](../examples/femahoy/grid_frequency_cantilever.skn),
    [grid_free_block_modes_ccx.skn](../examples/femahoy/grid_free_block_modes_ccx.skn),
    [buckling_column.skn](../examples/femahoy/buckling_column.skn) and
    [grid_cavity_te101.skn](../examples/femahoy/grid_cavity_te101.skn) use it.
  - G7: [Appendix A](skein_spec.md#appendix-a-built-in-names) states `re`,
    `imag`, `abs`, `arg` and `conj` of a `Complex` value. The expectations
    of the waveguide, line, acoustics, eddy current and Mie models check
    each part of the femahoy reference.
  - G9: [§3](skein_spec.md#3-grammar) lets a procedure argument be a
    list. The harmonic models write `freq: [f0]`.
  - G12: [Appendix A](skein_spec.md#appendix-a-built-in-names) states
    `atan` and `atan2`. [wedge_oblique_shock.skn](../examples/femahoy/wedge_oblique_shock.skn) checks the
    shock angle.
  - G16: [§12.1](skein_spec.md#121-form) reads `F*n = P` on a `Point` as a
    force. [cantilever_2step.skn](../examples/femahoy/cantilever_2step.skn) puts a force on each of the 9 tip
    nodes.
  - G29: [§6.3](skein_spec.md#63-units) reads `-5 degC` as one literal.
    [grid_layered_wall.skn](../examples/femahoy/grid_layered_wall.skn) uses it.
  - G30: [§6.1](skein_spec.md#61-dimensions) allows a constant `Real`
    exponent. [hyperelastic_tube.skn](../examples/femahoy/hyperelastic_tube.skn) and
    [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn) use it.
- Each other gap closes in the phase that builds its solver or its physics.
  A gap that blocks a case closes when its phase adds the model of the case.
- The shell, beam and laminate lane closes G3 before phase 5. Spec
  [§4.3](skein_spec.md#43-shells-beams-and-laminates) states the shell template, and `src/shells/` solves it. The new model
  [scordelis_lo_roof.skn](../examples/femahoy/scordelis_lo_roof.skn) and
  the thin plate of
  [laminate_pagano.skn](../examples/femahoy/laminate_pagano.skn) use it.
  Phase 5 runs G3, G31 and G33 in the two 787 models with `skein run`:
  two laminates, offset beams of two sections, order 2, and the global
  values in the submodel ([notes](studies/b787/skein/notes.md#run-status)).
- Phase 2 closes G5. Spec [§5.6](skein_spec.md#56-tables) states `sum`,
  and the new model
  [fit_nist_pontius.skn](../examples/femahoy/fit_nist_pontius.skn) uses it.
- Phase 2 closes G10. Spec [§15.3](skein_spec.md#153-expect) lets an
  `expect` stand in a step and in a study. The step `inflate_25` of
  [hyperelastic_tube.skn](../examples/femahoy/hyperelastic_tube.skn) checks
  the force at 2.5. The studies of
  [cavity_devahldavis.skn](../examples/femahoy/cavity_devahldavis.skn),
  [naca0012_ladson.skn](../examples/femahoy/naca0012_ladson.skn),
  [naca23012_flap.skn](../examples/femahoy/naca23012_flap.skn) and
  [dipole_king_middleton.skn](../examples/femahoy/dipole_king_middleton.skn)
  check each case, and the `minimize` study of
  [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn)
  checks its best case.
- Phase 2 closes G11. Spec [§14.5](skein_spec.md#145-output-and-report)
  states `when`, and the expectations of
  [dam_break_vof.skn](../examples/femahoy/dam_break_vof.skn) and
  [dam_break_martin_moyce.skn](../examples/femahoy/dam_break_martin_moyce.skn)
  use it.
- Phase 2 closes G14. Spec [§6.1](skein_spec.md#61-dimensions) allows
  exponents in halves, and the expectation of
  [fracture_centre_crack.skn](../examples/femahoy/fracture_centre_crack.skn)
  compares `K_I` in `MPa*mm^0.5`.
- Phase 2 closes G15. Spec [§10.5](skein_spec.md#105-odes-and-global-unknowns)
  allows `mean` and `integral` in the `eq` of a global unknown, and
  [grid_fibre_composite_cell.skn](../examples/femahoy/grid_fibre_composite_cell.skn)
  uses it.
- Phase 2 closes G17. Spec [§12.1](skein_spec.md#121-form) states
  `rigid(p)` and `moment`, and the top of
  [contact_cattaneo.skn](../examples/femahoy/contact_cattaneo.skn) uses
  them.
- Phase 2 closes G18. Spec [§13](skein_spec.md#13-methods) gives `family`
  the values `serendipity` and `incompatible`, and `quadrature` the values
  `full` and `reduced`. The grid models use `family: incompatible`, and
  `skein run` runs
  [grid_cantilever.skn](../examples/femahoy/grid_cantilever.skn) on 8-node
  hexahedra with incompatible modes. The elements lane (`src/elements/`)
  runs the same family in a `modal` step:
  [grid_frequency_cantilever.skn](../examples/femahoy/grid_frequency_cantilever.skn)
  gives 419.52 Hz (+0.42 %), and
  [grid_free_block_modes_ccx.skn](../examples/femahoy/grid_free_block_modes_ccx.skn)
  gives its 6 elastic modes within 0.7 % of CalculiX, after 6 rigid modes at 0 Hz.
  [frequency_cantilever.skn](../examples/femahoy/frequency_cantilever.skn)
  needs `family: serendipity` with `quadrature: reduced`, the 20-node
  hexahedron. It does not run yet.
- Phase 2 runs the hex8 heat cells in the elements lane (`cell: hex` on a
  `box`, `order: 1`): [grid_heated_bar.skn](../examples/femahoy/grid_heated_bar.skn),
  [grid_conductivity_slab.skn](../examples/femahoy/grid_conductivity_slab.skn)
  (k of T), [grid_radiation_slab.skn](../examples/femahoy/grid_radiation_slab.skn)
  and [grid_radiating_fin_ccx.skn](../examples/femahoy/grid_radiating_fin_ccx.skn)
  (radiation to the surroundings, `boundary(fin) - base`), and
  [grid_cooling_slab.skn](../examples/femahoy/grid_cooling_slab.skn)
  (`backward_euler` and `crank_nicolson`). Each passes its `expect`.
- Phase 2 closes G20. Spec [§9.1](skein_spec.md#91-mesh) states
  `crack_tip`, and the mesh of
  [fracture_centre_crack.skn](../examples/femahoy/fracture_centre_crack.skn)
  uses it.
- Phase 2 closes G24. Spec [§8.1](skein_spec.md#81-distance-functions)
  states `fill`, the region inside a `spline` and the line from its last
  point to its first. The airfoils of
  [naca0012_ladson.skn](../examples/femahoy/naca0012_ladson.skn) and
  [naca23012_flap.skn](../examples/femahoy/naca23012_flap.skn) and the
  hole of
  [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn) use
  it.
- Phase 2 closes G25. Spec [§20.2](skein_spec.md#202-study) gives
  `minimize` the arguments `xtol` and `step`, and the study of
  [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn)
  uses them.
- Phase 2 closes G19. Spec [§9.1](skein_spec.md#91-mesh) and
  [§9.2](skein_spec.md#92-mesh-generation) allow `hex` and `quad` on a
  union of boxes or rects, with a grid plane on each face of each box or
  rect of the model. The layers of
  [warpage_thick_copper.skn](../examples/femahoy/warpage_thick_copper.skn),
  [grid_layered_wall.skn](../examples/femahoy/grid_layered_wall.skn) and
  [grid_laminate_cell.skn](../examples/femahoy/grid_laminate_cell.skn),
  and the 2 boxes of
  [grid_heated_block_ccx.skn](../examples/femahoy/grid_heated_block_ccx.skn)
  and
  [grid_heated_block_stress_ccx.skn](../examples/femahoy/grid_heated_block_stress_ccx.skn),
  use it.
- Phase 2 closes G28. Spec [§3.1](skein_spec.md#31-order-of-items) and
  [§10.2](skein_spec.md#102-material) let a material property name an
  unknown of a physics that assigns the material. The conductivity of `T`
  of
  [grid_radiating_fin_ccx.skn](../examples/femahoy/grid_radiating_fin_ccx.skn)
  and
  [grid_conductivity_slab.skn](../examples/femahoy/grid_conductivity_slab.skn)
  is a material property.
- Phase 3 closes G36. Spec [§9.1](skein_spec.md#91-mesh) and
  [§9.2](skein_spec.md#92-mesh-generation) allow `quad` and `hex` on each
  region: a grid in r and θ on a quarter of a ring or a tube, Blossom-Quad
  on another `Surface`, a swept grid on a prism along z, and a mapped grid
  on another `Solid`. The ring of
  [submodel_kirsch_hole.skn](../examples/femahoy/submodel_kirsch_hole.skn)
  and the package of
  [warpage_2p5d_frame.skn](../examples/femahoy/warpage_2p5d_frame.skn) use
  the cells of femahoy. The ring of
  [hyperelastic_tube.skn](../examples/femahoy/hyperelastic_tube.skn), the
  rod of
  [taylor_bar_impact.skn](../examples/femahoy/taylor_bar_impact.skn), the
  bar of
  [dynamics_bar_step.skn](../examples/femahoy/dynamics_bar_step.skn) and
  the parts of
  [contact_cattaneo.skn](../examples/femahoy/contact_cattaneo.skn) use
  quadratic triangles and tetrahedra: the nonlinear solid and the contact
  of phase 5 (`src/elements/solid.goose`, `contact.goose`) run on
  simplices only.
- Phase 5 closes G32 and G35. The captain put them in phase 5 on
  2026-10-31. Spec [§8.4](skein_spec.md#84-region-algebra) gives `a & b`
  the lower dimension of `a` and `b`, so the 787 models
  [fuselage.skn](studies/b787/skein/fuselage.skn) and
  [door_corner.skn](studies/b787/skein/door_corner.skn) pass `skein check`
  with no mark. Spec [§8.1](skein_spec.md#81-distance-functions) states
  `copy(r)`, the index of the copy of a pattern at a point. `skein run`
  runs both models with the cuts, and
  [window_barrel.skn](../examples/shells/window_barrel.skn) runs `copy(r)`.
- Coupled physics on one mesh runs five more models with `skein run`:
  [grid_resistor_bar.skn](../examples/femahoy/grid_resistor_bar.skn),
  [grid_joule_bar.skn](../examples/femahoy/grid_joule_bar.skn),
  [grid_holm_spot.skn](../examples/femahoy/grid_holm_spot.skn),
  [line_coax.skn](../examples/femahoy/line_coax.skn) and
  [line_stripline.skn](../examples/femahoy/line_stripline.skn).
  The heat then stress of
  [grid_heated_block_stress_ccx.skn](../examples/femahoy/grid_heated_block_stress_ccx.skn)
  (hex cells on a union of boxes, laws by region, a film bc) and of
  [grid_bimetal_strip.skn](../examples/femahoy/grid_bimetal_strip.skn) runs
  in the elements lane. Each region has its own conductivity, source and
  thermal strain. The 3D mesher keeps the boundary of a region as an
  interface: it cuts each tetrahedron that the boundary crosses
  (`cut_interface` in `src/mesher/tet.goose`), so
  [grid_plate_hole_ccx.skn](../examples/femahoy/grid_plate_hole_ccx.skn) and
  [grid_plate_hole_stress_ccx.skn](../examples/femahoy/grid_plate_hole_stress_ccx.skn)
  pass.
- A written-out linear elastic law runs with `skein run` in the elements
  lane (`src/elements/gen_aniso.goose`): `let sigma` is any expression that
  is linear in `strain(u)`, with `outer`, `tr`, `dev`, `sym`, `dot`, `ddot`,
  `hooke`, the members `.xx` of a tensor, a thermal strain tensor and the
  lets and the material properties of the physics. The generator evaluates it
  at the 6 unit strains and at the zero strain, and checks that it is linear.
  Each material then has one 6 x 6 stiffness and one stress of the zero strain
  (`elements::anisotropic`, `anisotropic_in` for a region), which the cells of
  each element, the incompatible modes of the hexahedron and the stress output
  all use. The
  rotation of a ply is a `let` of the law, so a layup of plies is one material
  for each ply on a region. [grid_constrained_block.skn](../examples/femahoy/grid_constrained_block.skn)
  passes its `expect` (the orthotropic block under a heating of 100 K: the
  three stresses within 1e-6 MPa), and
  [two_ply_block_hex.skn](../examples/elements/two_ply_block_hex.skn) checks
  a region for each ply. A plane stress law is condensed to the strains xx, yy and xy.
  A transient or buckling step, a contact, the method `rkpm`, a submodel
  source and a plastic block with a written-out law are not in this generator,
  and `skein run` says so.
- The flow coupling runs two more models with `skein run`:
  [cavity_devahldavis.skn](../examples/femahoy/cavity_devahldavis.skn)
  (Boussinesq buoyancy, both Ra cases of the study) and
  [cht_heated_channel.skn](../examples/femahoy/cht_heated_channel.skn)
  (a solid and a fluid on one mesh, T within 0.006 K). `src/flow` has a T
  dof at each node, a conductivity for each cell, a body force that is
  affine in T, natural boundary fluxes, and the mean, the integral and the
  maximum of any field expression on a region. A banded LU solve
  (`src/solvers/band.goose`) replaces GMRES for the thin channel.
- [wedge_oblique_shock.skn](../examples/femahoy/wedge_oblique_shock.skn) runs
  with `skein run` on the conservation-law family `src/cons/`: any
  `ddt(U) + div(F) = S` system of an `fv` method with more than one unknown
  component, on triangles, with a MUSCL limiter, a Rusanov flux whose speed
  comes from the Jacobian of the model flux, and SSPRK3. All five expects pass
  within 0.25 % (the shock angle is 36.998 deg).
- Not run yet in the flow batch, each needs a new solver family:
  naca0012_ladson and naca23012_flap (Spalart-Allmaras with a nuTilde unknown,
  a wall-distance physics, a spline airfoil with a far field, boundary layer
  cells and a second physics), dam_break_vof and dam_break_martin_moyce
  (two-phase VOF with tables, surface tension and a Tait water; the flow
  generator rejects the Newtonian stress `let`).
- turek_hron_fsi1 runs in the flow lane (drag +0.02 %, lift -0.06 %, ux_a +0.1 %,
  uy_a -2.1 %). The flag is written out as St. Venant-Kirchhoff
  (`assume finite_strain`, `hooke` of `strain(u)`, `div(defgrad(u)*S) = 0`). It is solved by
  Newton in the total Lagrangian form (`src/flow/fsi.goose`). Deviations: the model uses
  `method flow: cg` with `order: { v: 2, p: 1 }` in place of `fv`, and tri cells for the flag
  in place of quads.
- slider_crank runs in the `ode` lane (`src/ode`): exact Taylor series in time and Newton on the
  history polynomial. `generalized_alpha(rho)` is mapped to 6-point backward differences, so the
  numerical dissipation of the alpha method is not reproduced.
- acoustics_pulsating_sphere runs in the wave lane: a Robin boundary term, complex sample
  expressions, and a graded mesh of the octant shell (`shell_tets` in `src/mesher/sweep.goose`).

- Structural models of phase 3 run with `skein run` and pass their `expect`:
  [buckling_column.skn](../examples/femahoy/buckling_column.skn),
  [birth_death_column.skn](../examples/femahoy/birth_death_column.skn),
  [submodel_kirsch_hole.skn](../examples/femahoy/submodel_kirsch_hole.skn)
  (with [submodel_kirsch_hole_global.skn](../examples/femahoy/submodel_kirsch_hole_global.skn))
  and [fracture_centre_crack.skn](../examples/femahoy/fracture_centre_crack.skn).
  Deviations: the quad cells of the submodel ring are meshed as P2
  triangles. The preload form `P = defgrad(u)*S` solves the static step to
  first order, and the buckling step takes its stress. A born node shows
  the change since its birth. Point location on a curved P2 cell uses
  Newton iteration. A `||` boundary with a round limit inside is refused
  (`tests/unsupported/boundary_or_limit.skn`).
- The harmonic step of a `Complex` vector unknown runs with `skein run`
  (`src/elements/harmonic.goose`): it is the complex lift of the real
  assembly of a physics, with the matrix `(1 + 2 im zeta) K - omega^2 M` and
  a complex skyline LDLT. The elastic physics uses it. The heat physics has
  the harmonic step in its natural form: `solves T: Complex<Temperature>`
  with `let q = -k*grad(T)` and `eq rho*cp*ddt(T) + div(q) = q_v`, where `ddt`
  is `im*omega`. It takes the route of the scalar wave physics with the
  diffusion `k`, `c = -im*omega*rho*cp/k` and the source `-q_v/k`, so a
  fixed temperature, a source and the direct sparse LU are the same.
  [harmonic_heat_bar_tri.skn](../examples/elements/harmonic_heat_bar_tri.skn)
  matches the closed form `cosh(m (L - x))/cosh(m L)` of a bar at two
  frequencies within 1e-6. [random_vibration_cantilever.skn](../examples/femahoy/random_vibration_cantilever.skn)
  passes its `expect`: the RMS tip acceleration is 402270 mm/s^2 (+0.34 %,
  Blevins 1979: 400925.4 mm/s^2). A report or an `expect` of a harmonic step
  takes a `Complex` expression of constants, `omega`, `mean(u.c, region)`,
  `u.c at p` and `rms(f, psd)`. The VTU output has the real part.
- The shape loop runs with `skein run --study shape`
  ([shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn)):
  `minimize(method: nelder_mead)` (`src/codegen/run.goose`) is a simplex
  search on the params in 0..1 of their ranges, and each case remeshes and
  solves the model. The mesher reader takes `fn` items, `sin`, `cos`, `tan`,
  `exp`, `ln`, `sqrt` and `abs` in a shape, and the `nondim` unit system runs.
  An `expect` of a param, a let or `max(von_mises(sigma), r)` runs. The 51
  cases take 23 min, because each builds its own program, so `./build.sh test`
  runs only the default case. The best case is ratio 2.0003, power 2.003 and a
  peak of 1.5008 (the reference: ratio 2, power 2, peak 1.5), and its three
  `expect` lines pass within 0.5 %.
- The wave-port batch runs three waveguide models with `skein run`:
  [grid_waveguide_te10.skn](../examples/femahoy/grid_waveguide_te10.skn)
  (phase of S21 within 0.1 deg),
  [grid_waveguide_step.skn](../examples/femahoy/grid_waveguide_step.skn)
  and [grid_waveguide_dielectric.skn](../examples/femahoy/grid_waveguide_dielectric.skn)
  (S11 and S21 within 0.035 and 0.012). The edge physics takes the stress
  tensor form `let Hs = -2*skew(grad(E))/mu`, `eq div(Hs) + eps*ddt2(E) = 0`
  (`ddt2` is `-omega^2` in a `harmonic` step) and the natural bc
  `Hs*n = z*(E - dot(E, n)*n) + g` with `g` a vector of x, y and z (a port
  with an incident wave). The bcs `E.y = 0` and `E.z = 0` on a plane of x
  hold the tangent of E. A `material a on r1, b on r2` gives each region its
  own mu and eps, and `integral(dot(E, g), s)` on a boundary projects E on a
  mode. A limit of a `&&` inside a `||` of a boundary cuts only the planes of
  that `&&` (`plane_cut` in `src/elements/api.goose`). Deviations: the mesh
  sizes are 1.6 mm (te10, for the phase), 2 mm (step) and 2 mm (dielectric),
  because the tetrahedra of the mesher need the finer mesh of te10 to put the
  phase of S21 within 0.1 deg. The step and the dielectric guide also run
  at the femahoy sizes (1.4 mm and 1.5 mm): the step takes 3 min 19 s and
  4.1 GB (142,468 dofs, estimate 4.05 GB), and the dielectric guide takes
  2 min 36 s and 2.7 GB (90,146 dofs, estimate 2.66 GB). Their four `expect`
  lines pass. They are too long for `./build.sh test`, so the logs are in
  `data/skein-fh-em2/` of the supervisor.
- The wave and edge systems take a direct solve (`src/solvers/slu.goose`): LU
  factors with no pivoting, on a nested dissection of the graph of the
  matrix (a level set of a search for each separator), with a pivot that is
  almost 0 replaced. The factor is also the preconditioner of GMRES. Reason:
  GMRES with ILU(0) did not converge on the step at 1.4 mm in 5 minutes, and
  its cap of 20,000 iterations was never reached. The plan counts the entries
  of L and the multiplies before the factors are sized, so the memory check
  and the time estimate come first. The run prints the estimate if it is
  over 20 s. A system with an estimate over 900 s (15 minutes) stops before
  the solve with a message that gives the estimate and says to use a coarser
  mesh, so a run never looks like a hang. The step at 1.4 mm (202 s estimate,
  186 s measured) and the dielectric at 1.5 mm (134 s estimate, 143 s
  measured) run within the limit.
- The batch of the eddy currents runs
  [line_wire_skin.skn](../examples/femahoy/line_wire_skin.skn) (R within
  0.1 % of the Bessel function answer) and
  [line_bars_getdp.skn](../examples/femahoy/line_bars_getdp.skn) (R and L
  within 0.3 % and 0.06 % of GetDP). The scalar wave physics takes the form
  `eq -div(Hr) = J` with `let Hr = grad(Az)/mu` and `let J = sigma_e*(E - ddt(Az))`:
  the generator takes c, s and 1/mu of each region from its material at the
  omega of the step (`ddt(u)` is `im*omega*u`), and `integral(J, r)` of a let
  over a region of cells is `a I + b A` with the material of the region. A
  `refine` of `boundary(shape)` or of `boundary(a) | boundary(b)` refines the
  surface of the shape.
- Not run yet in the electromagnetic batch, each for a large missing piece:
  grid_cavity_te101 needs the modal step of an edge physics and a wall
  impedance that depends on the eigenvalue, and its `expect` is within 3e-4
  of the frequency: the lowest order edge element needs a mesh of about 1
  million cells for that, so it needs the second order element first.
  team7_eddy_current needs the `stadium` shape of the mesher and a source in
  the edge physics. dipole_king_middleton needs order 2 edge elements and a
  PML in 3D. optics_mie_cylinder needs a source that varies in space in the
  scalar wave physics (`e_inc`), a PML on a layer that is not a box, order 2,
  the integral of a field on a curve inside the mesh, and the `study` of a
  param.
- Not run yet, each for a large missing piece:
  laminate_pagano and warpage_2p5d_frame do not run in this batch, and each
  needs more than one new piece. The written-out orthotropic law and the
  material of each ply on a region, which both need, run now
  (`two_ply_block_hex.skn`). laminate_pagano also needs: second order
  serendipity hexahedra with reduced integration (the hexahedra of the
  elements lane have order 1); a traction that varies with `x` and `y` on a
  face of the hexahedra (`sin(pi*x/a)`); a composite `reissner_mindlin`
  shell on quad cells of a flat plate (the shell lane takes only a part of a
  cylinder, and its `layup` table works there); and one model with two
  physics that run in two lanes (the elements lane for the 3D plate, the shell
  lane for the thin plate). warpage_2p5d_frame also needs: the second order
  serendipity hexahedra, and a hexahedron grid on boxes that a turned box cuts
  (the sloped fillet; the hex mesher takes a union of boxes).

The `// blocked-by:` mark of a model names the gaps of the model that stay
open after phase 1.

### Gaps that block a case

Skein cannot state the core physics of these cases, so they have no model.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G1 | Topology optimization: a density for each cell, SIMP, a sensitivity filter, an optimality criteria update, a volume limit. Then a solid from the density at 0.5. | `validation/cases/topology-mbb-beam/derivation.md:24`, `validation/cases/topology-mbb-beam-cad/derivation.md:23` | `minimize` (§20.2) varies only a few `param` values, and a `param` is a constant (§20.1). Geometry comes only from CSG of constants (§8). | A study procedure `topology(objective, region, volume, penalty, filter)`, and a geometry function `level(f, c)`: the region where a `Results` field `f` is `c` or more. | 5, closed by [§20.2](skein_spec.md#202-study): `topology(..., threshold: 0.5)` makes the solid of the cells at 0.5 or more, in place of `level(f, c)` |
| G2 | Random vibration: a modal damping ratio, a PSD input and the RMS response over a frequency band. | `validation/cases/random-vibration-cantilever/job.inp:14`, `:25`, `check.toml:9` | `harmonic` (§14.2) has no damping ratio for each mode, and no report function integrates over the frequencies of a step (§14.5). Also needs G8 and G9. | A `harmonic` argument `damping: Real`, and a report function `rms(f, psd)` over the frequencies of the step. | 4, closed by [§14.2](skein_spec.md#142-procedures) and [§14.5](skein_spec.md#145-output-and-report) |
| G3 | Shells (`S8R`, `S6`), also a composite shell section with plies. | `validation/cases/scordelis-lo-roof/job.inp:26`, `validation/cases/laminate-pagano/job.inp:16` | §4.2 names `reissner_mindlin(t)`, but gives no source, no types of `w` and `theta`, and no load argument. Its parameter `t` is also the built-in time. The thin plate of `laminate-pagano` is missing from its model. | Give the template source in §4.2, with a parameter `h: Length`, a load argument, typed unknowns and an optional `layup` table of plies (thickness, angle, material). | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) |
| G4 | A perfect electric conductor on a curved wall: only the tangential part of the field is zero. | `validation/cases/grid-coax-tem/coax/grid.toml:1` | The left side of an essential `bc` is an unknown or a component of it (§12.1). A component is along x, y or z, not along the wall. | Accept `cross(n, u)` as the left side of an essential condition. | 5, closed by [§12.1](skein_spec.md#121-form) |
| G5 | A least-squares curve fit to 40 points. | `validation/cases/fit-nist-pontius/pontius/fit.toml:4` | A table feeds only `interp` and the rows of an `expect` (§5.6, §15.3). No construct sums an expression over the rows. `minimize` has no gradients, so it cannot get to a relative error of 1e-10. | A function `sum(expr)` over the rows of a table. A `nondim` model can then state the normal equations. Or: state that a curve fit is not a Skein model. | 2, closed by [§5.6](skein_spec.md#56-tables) |

### Gaps that make a model fail `check`

The model states the physics, but the merged spec rejects the construct at
the `// gap` comment. The model fails `check` until the spec change of
phase 1 closes the gap.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G9 | The frequencies of a harmonic step: one frequency, many frequencies, or a band of N points. A model with one frequency fails `check`. A case with many frequencies or a band also loses the checks at the other frequencies. | `validation/cases/line-wire-skin/wire/grid.toml:27`, `validation/cases/line-bars-getdp/bars/grid.toml:35`, `validation/cases/random-vibration-cantilever/job.inp:24` | `harmonic(freq:)` needs a list (§14.2), but a procedure argument is an `expr` (§3), and `expr` has no list. | Let a procedure argument be a `list` (§3), as in a study. | 1, closed by [§3](skein_spec.md#3-grammar) |
| G30 | A power with an exponent that is not an integer constant: the Ogden terms with alpha = 1.3, 5.0 and -2.0, the points of a superellipse hole. | `validation/cases/hyperelastic-tube/job.inp:12`, `validation/cases/shape-hole-biaxial/L40/mesh.geo:15` | §6.1 needs an integer constant `k` in `a^k`. The models write `jac^(-al/3)` and `cos(k*dth)^(2/power)`. | In §6.1, let `k` be any constant `Real` when `a` is `Real`. | 1, closed by [§6.1](skein_spec.md#61-dimensions) |

### Gaps that remove a check

The model states the physics, but it cannot state some reference values of
`check.toml`. The model reports a related value, or it leaves a `// gap`
comment.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G6 | The frequency of mode k (`FREQ`) and the buckling factor of mode k (`BUCKLE`). | `validation/cases/frequency-cantilever/check.toml:9`, `validation/cases/buckling-column/check.toml:7`, `validation/cases/grid-free-block-modes-ccx/check.toml:13`, `validation/cases/grid-cavity-te101/check.toml:10` | `modal` and `buckling` give eigenpairs (§14.2), but no expression names an eigenvalue (§14.5, Appendix A). | A built-in `eigenvalue(k: Int)`, for `report` and `expect`: a `Frequency` after `modal`, a `Real` load factor after `buckling`. | 1, closed by [§14.5](skein_spec.md#145-output-and-report) |
| G7 | The real part, the imaginary part, the magnitude and the phase of a complex result (S-parameters, impedance, R and L of a line). | `validation/cases/grid-waveguide-te10/check.toml:10`, `validation/cases/line-bars-getdp/check.toml:11`, `validation/cases/acoustics-pulsating-sphere/check.toml:17`, `validation/cases/team7-eddy-current/50hz/case.pro:89` | §5.1 defines `Complex<D>` and `im`, but Appendix A has no function that takes a part of a complex value. | Add `re`, `imag`, `arg` and `conj` to Appendix A, and accept `Complex` in `abs`. | 1, closed by [Appendix A](skein_spec.md#appendix-a-built-in-names) |
| G8 | The angular frequency of the step in a load or a wall impedance: a base acceleration, a surface impedance of copper. | `validation/cases/random-vibration-cantilever/job.inp:27`, `validation/cases/grid-cavity-te101/cavity/grid.toml:20` | No built-in name gives the frequency of a `harmonic` or `modal` step (§14.2, Appendix A). | A built-in constant `omega`, the angular frequency of the current `harmonic` frequency or `modal` eigenvalue. | 4, closed by [§14.2](skein_spec.md#142-procedures) |
| G10 | A check after an earlier step, a check for each case of a parameter sweep, and a check of the optimum of a shape loop. | `validation/cases/hyperelastic-tube/check.toml:18`, `validation/cases/cavity-devahldavis/check.toml:40`, `validation/cases/naca0012-ladson/check.toml:19`, `validation/cases/shape-hole-biaxial/check.toml:23` | §15.3 runs each `expect` once, after the last step, with the default `param` values. | Let `expect` stand in a `step` (checked at the end of the step) and in a `study` block (checked for each case, and for the best case of `minimize`). | 2, closed by [§15.3](skein_spec.md#153-expect) |
| G11 | The time at which a value gets to a level: the front of a dam break. | `validation/cases/dam-break-vof/check.toml:19`, `validation/cases/dam-break-martin-moyce/check.toml:17` | Reports give values at the end of a step only (§14.5, §15.3). | A report function `when(expr, value)`: the first time at which `expr` gets to `value`. | 2, closed by [§14.5](skein_spec.md#145-output-and-report) |
| G12 | An angle from a ratio: the shock angle. | `validation/cases/wedge-oblique-shock/check.toml:19` | Appendix A has `tan` but no inverse. The model checks `tan(beta)`. | Add `atan` and `atan2` to Appendix A. | 1, closed by [Appendix A](skein_spec.md#appendix-a-built-in-names) |
| G13 | The directivity and the efficiency of an antenna from a near-to-far-field transform. | `validation/cases/dipole-king-middleton/omega12.5/openems.py:65` | Report functions act on regions of the mesh only (§14.5). | A report function `far_field(E, H, surface, direction)`. | 5, closed by [§14.5](skein_spec.md#145-output-and-report) |
| G14 | A stress intensity factor in MPa·mm^0.5. | `validation/cases/fracture-centre-crack/check.toml:7` | A dimension has integer exponents, and `sqrt` needs even exponents (§6.1). The model checks the shape factor `K_I/(S*sqrt(pi*a))`, with the same tolerance. | Allow exponents in halves (§6.1). | 2, closed by [§6.1](skein_spec.md#61-dimensions) |
| G15 | The effective Young's modulus of a periodic cell under a mean stress. | `validation/cases/grid-fibre-composite-cell/check.toml:31` | The mean strain must be a global unknown that a mean stress fixes, but `mean` is valid only in `output`, `report` and `expect` (§14.5). | Allow `mean(f, r)` and `integral(f, r)` in the `eq` of a global unknown (§10.5). | 2, closed by [§10.5](skein_spec.md#105-odes-and-global-unknowns) |
| G16 | A concentrated force on a node (`*CLOAD`). | `validation/cases/topology-mbb-beam/60x20/job.inp:2517`, `fixtures/results/cantilever-2step/cantilever-2step.inp` | §8.3 says that a `Point` can be the region of a point load, but §12.1 gives no left side for a force at a point: `n` has no value there. `cantilever_2step.skn` uses a traction with the same resultant. | On a `Point` region, read a natural condition `F*n = P` as the force `P`, with the dimension of `F` times area. | 1, closed by [§12.1](skein_spec.md#121-form) |
| G17 | A rigid constraint of a boundary (`*EQUATION`): the face stays straight and can turn. | `validation/cases/contact-cattaneo/mesh.inp:13915` | §12 has only essential and natural conditions. The model puts tractions with the same force and moment on the face. | A built-in `rigid(p)` for the right side of an essential `bc`. A `bc` on the `Point` `p` takes the force and the moment. | 2, closed by [§12.1](skein_spec.md#121-form) |

### Gaps that change the discretization or the physics

The model states the problem, but the result can differ from the femahoy
result, or the mesh or the solve is less efficient.

| Gap | Femahoy feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G18 | Element families: reduced integration (`C3D20R`, `C3D8R`, `CPE8R`, `S8R`), serendipity cells (`C3D20`, `CPE8`) and incompatible modes (`C3D8I`, the hex8 of the grid solver). | `validation/cases/frequency-cantilever/mesh.inp:2120`, `validation/cases/dynamics-bar-step/mesh.inp:2028`, `validation/cases/grid-cantilever/derivation.md:3` | The setting `quadrature` has no values, and `family` is `lagrange`, `nedelec` or `raviart_thomas` (§13). The models use `lagrange` of order 2. | Give `quadrature` the values `full` and `reduced`, and add `serendipity` and `incompatible` to `family`. | 2, closed by [§13](skein_spec.md#13-methods) |
| G19 | A structured hex grid through layers and blocks, with a grid plane on each material face. | `validation/cases/warpage-2p5d-frame/mesh.inp:5462`, `validation/cases/grid-heated-block-ccx/block/grid.toml:16`, `validation/cases/grid-layered-wall/wall/grid.toml:18` | `hex` and `quad` need one `box` or `rect` (or, for `hex`, a quarter of a tube) with one `size` (§9.1). §9.2 does not put a grid plane on an inner material face. The models use tetrahedra, or a size that puts the planes on the faces. | Allow `hex` and `quad` on a union of axis-aligned boxes or rects, with a grid plane on each face of each box. | 2, closed by [§9.2](skein_spec.md#92-mesh-generation) |
| G36 | A structured hex or quad grid on a region that is not a union of axis-aligned boxes or rects: the sloped fillet of a package, a ring about a hole, a tube, a disk on a block. | `validation/cases/warpage-2p5d-frame/mesh.inp:5462`, `validation/cases/contact-cattaneo/mesh.inp:13915` | §9.1 and §9.2 allow `quad` on another `Surface` (Blossom-Quad) and `hex` on a quarter of a tube, but the mesher makes `quad` and `hex` only on a union of boxes or rects. The spec has no `hex` grid on another `Solid` (the fillet, the rod). The models submodel-kirsch-hole, hyperelastic-tube, contact-cattaneo, turek-hron-fsi1, taylor-bar-impact and warpage-2p5d-frame use triangles or tetrahedra. | In the mesher, Blossom-Quad and the quarter-tube grid of §9.2. In §9.2, a mapped or swept `hex` grid on a `Solid` that is not a union of boxes. | 3, closed by [§9.2](skein_spec.md#92-mesh-generation) |
| G20 | Quarter-point cells at a crack tip. | `validation/cases/fracture-centre-crack/derivation.md:38` | No mesh setting moves the mid-side nodes (§9.1). The model uses a fine mesh. | A mesh entry `crack_tip r: size`. | 2, closed by [§9.1](skein_spec.md#91-mesh) |
| G21 | The stable increment of an explicit step, chosen by the solver. | `validation/cases/dynamics-bar-step/job.inp:14` | `transient` needs `dt` (§14.2). The model gives a `dt` below the stable limit. | Accept `dt: stable` with an explicit scheme. | 2, closed by [§14.2](skein_spec.md#142-procedures) |
| G22 | Finite-strain J2 plasticity with a multiplicative split (plastic strain up to 3.24). | `validation/cases/taylor-bar-impact/job.inp:12` | §10.6 states `plastic` with an additive plastic strain. Under `finite_strain` this is a different law at large strain. | One rule in §10.6: under `finite_strain`, the `strain` state is `Fp`, the elastic part is `defgrad(u)*inv(Fp)`, and the return map is the exponential map. | 5, closed by [§10.6](skein_spec.md#106-inelastic-materials) |
| G23 | The stick stiffness of penalty friction. It changes the reference traction by 0.34 % to 3.4 %. | `validation/cases/contact-cattaneo/job.inp:15` | `contact` has only `friction`, `enforce` and `stiffness` (§12.4). | A `contact` setting `stick_stiffness`. | 5, closed by [§12.4](skein_spec.md#124-contact) |
| G24 | Smooth closed shapes from points or parameters: the airfoil splines, a superellipse hole. | `validation/cases/naca0012-ladson/a4.04/system/blockMeshDict:67`, `validation/cases/shape-hole-biaxial/L40/mesh.geo:9` | `spline` gives a `Curve`, and no function makes a `Surface` from a closed curve (§8.1). The models use polygons. | A geometry function `fill(c: Curve)` for a closed curve. | 2, closed by [§8.1](skein_spec.md#81-distance-functions) |
| G25 | The stop rule and the first simplex of Nelder-Mead. The case tolerance comes from the stop rule. | `validation/cases/shape-hole-biaxial/derivation.md:50` | `minimize` has only `objective`, `method` and `max_cases` (§20.2). | `minimize` arguments `xtol` and `step`. | 2, closed by [§20.2](skein_spec.md#202-study) |
| G26 | An open boundary: a PML, or a map of a shell to infinity. | `validation/cases/team7-eddy-current/50hz/case.pro:50`, `validation/cases/dipole-king-middleton/omega12.5/openems.py:39`, `validation/cases/optics-mie-cylinder/res40/meep_case.py:38` | §8 and §13 have no absorbing layer and no infinite map. The models put a first-order absorbing `bc`, or a zero field, on a far boundary. | A method setting `pml: { region, thickness }`. | 4, closed by [§13](skein_spec.md#13-methods) |
| G27 | A tree gauge for edge elements in a region with no conductor. | `validation/cases/team7-eddy-current/50hz/case.pro:36` | `family: nedelec` has no gauge setting (§13). | A method setting `gauge: tree`. | 5, closed by [§13](skein_spec.md#13-methods) |
| G28 | A material property that is a function of the unknown of the physics that uses it: a conductivity of T. | `validation/cases/grid-radiating-fin-ccx/fin/grid.toml:9` | §10.2 allows a property of a temperature field, but §3.1 needs each name before its use, and the `material` comes before the `physics` that declares `T`. The model puts the property in a physics `let`. | In §3.1, let a material property name an unknown of a physics that assigns the material. | 2, closed by [§3.1](skein_spec.md#31-order-of-items) |
| G29 | A temperature below 0 °C, for example a sink at -5 °C. | `validation/cases/grid-layered-wall/wall/grid.toml:49` | `-` is a unary operator (§2.4) and `5 degC` is 278.15 K (§6.3), so `-5 degC` is -278.15 K. The model gives the values in K. | A `-` directly before a literal with `degC` is part of the literal. | 1, closed by [§6.3](skein_spec.md#63-units) |

### Gaps of the 787 study

The [787 fuselage study](studies/b787/skein/notes.md) is not a femahoy
case, but its two models use the same marks. Its evidence paths are
relative to the root of this repository. The composite skin of the models
uses the shell template of G3 (`docs/studies/b787/skein/fuselage.skn:150`, `docs/studies/b787/skein/door_corner.skn:158`).

| Gap | Study feature | Evidence | Why the spec cannot state it | Smallest spec change | Phase |
|---|---|---|---|---|---|
| G31 | Offset beams on a shell: the stringers and the frames, with a torsion constant, a section axis, an offset and the axial force. | `docs/studies/b787/skein/fuselage.skn:166`, `:184`, `:204`, `docs/studies/b787/skein/door_corner.skn:175` | §4.2 gives `timoshenko(A, I_y, I_z)` with no source, no torsion constant, no section axis, no offset and no axial force. No rule ties a beam to the shell points on its line. | `timoshenko` arguments `J`, `y_axis` and `offset`, and a template `let Nx` (axial force). In §12.1, let the right side of an essential `bc` be an unknown of another physics on the same points (a tie). | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) and [§12.1](skein_spec.md#121-form) |
| G32 | A shell or a beam line cut by a solid: `skin & door_box(pad_w)`, `stringers & left`. | `docs/studies/b787/skein/fuselage.skn:61`, `:70`, `:86`, `docs/studies/b787/skein/door_corner.skn:65`, `:75`, `:78`, `:91` | §5.3 and §8.4 do not give the dimension of `a & b` when `a` is a `Surface` or a `Curve` and `b` is a `Solid`. The checker gives the higher dimension, a `Solid`. | In §8.4: `a & b` has the lower dimension of `a` and `b`, and `a - b` has the dimension of `a`. | 5, closed by [§8.4](skein_spec.md#84-region-algebra) |
| G33 | Shell resultants: the force resultant `N`, the membrane strain `eps_m`, an edge load `N*n`, and a strain of the global model in a submodel. | `docs/studies/b787/skein/fuselage.skn:160`, `:201`, `docs/studies/b787/skein/door_corner.skn:168`, `:203`, `:208`, `:214`, `:217` | The shell template of G3 declares no resultant. A natural `bc` needs a flux of an `eq` (§12.1). §15.1 gives a `Results` value the unknowns of the other model, not its physics `let` fields. | The G3 template gives `let N` and `let eps_m`, and its `eq` has `div(N)`, so `N*n` is a natural condition. In §15.1, a `Results` value also gives the physics `let` fields. | 5, closed by [§4.3](skein_spec.md#43-shells-beams-and-laminates) and [§15.1](skein_spec.md#151-results) |
| G34 | A submodel for each case of a parameter sweep: the door corner takes the global case with the same `p_U`. | `docs/studies/b787/compare/edits/e4_sweep.skein.patch:23` | §15.1 gives `results(file, step:)` no study case, so a `Results` value refers to one state only. | A `case:` argument of `results`, which names a `param` value of the other model. | 3, closed by [§15.1](skein_spec.md#151-results) |
| G35 | A value for each copy of a `pattern`: the centre of each window, for the window load. | `docs/studies/b787/compare/edits/e3_windows.skein.patch:53` | §8.1 gives `pattern(a, v, n)` the type of `a`, with no index and no value for each copy. The model computes x - x_w with `atan2`. | `pattern` gives an index `k` of the copy, and the position of each copy, for use in a `bc`. | 5, closed by [§8.1](skein_spec.md#81-distance-functions): `copy(r)` |

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
| Spline airfoils, a superellipse hole | naca0012-ladson, naca23012-flap, shape-hole-biaxial | `fill(spline(...))` of type `Surface`, closed by a line | [§8.1](skein_spec.md#81-distance-functions) |
| A solid from a density field | topology-mbb-beam-cad | `study s: topology(..., threshold: 0.5)`: the cells at 0.5 or more | [§20.2](skein_spec.md#202-study) |

### Elements and meshes

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `C3D10`, and the 10-node tetrahedra of the grid solver | bar-pressure, buckling-column, contact-hertz, plasticity-sphere, heat-fin, cantilever-2step, grid-kirsch-hole, grid-lame-cylinder, grid-round-cantilever, grid-plate-hole-stress-ccx | `mesh m on r { cell: tet }`, default `cg` with `order: 2` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| Linear tetrahedra of the grid solver (heat, conduction) | grid-pipe-wall, grid-plate-hole-ccx, grid-coax-capacitor, grid-holm-spot, grid-fibre-composite-cell | `cell: tet`, `method p: cg { order: 1 }` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| `C3D8`, and the hex8 heat cells of the grid solver | birth-death-column, grid-cooling-slab, grid-heated-bar, grid-layered-wall, grid-radiating-fin-ccx, grid-radiation-slab, grid-resistor-bar, grid-joule-bar, grid-conductivity-slab | `cell: hex` on a `box`, `order: 1` | [§9.1](skein_spec.md#91-mesh), [§13](skein_spec.md#13-methods) |
| `C3D8I`, and the hex8 with incompatible modes of the grid solver | grid-cantilever, grid-bimetal-strip, grid-frequency-cantilever, grid-free-block-modes-ccx, turek-hron-fsi1 | `cell: hex`, `order: 1`, `family: incompatible` | [§13](skein_spec.md#13-methods) |
| `C3D20`, `C3D20R`, `C3D8R`, `CPE8R`, `S8R` (serendipity, reduced integration) | frequency-cantilever, laminate-pagano, taylor-bar-impact, submodel-kirsch-hole, hyperelastic-tube | `family: serendipity` with `order: 2`, `quadrature: reduced`, or both | [§13](skein_spec.md#13-methods) |
| A hex grid through layers or a union of blocks | warpage-thick-copper, grid-heated-block-ccx, grid-layered-wall, grid-laminate-cell | `cell: hex` on a union of `box` regions, with a grid plane on each face of each box | [§9.1](skein_spec.md#91-mesh), [§9.2](skein_spec.md#92-mesh-generation) |
| A hex or quad grid on a region that is not a union of boxes or rects | warpage-2p5d-frame, submodel-kirsch-hole, hyperelastic-tube, contact-cattaneo, turek-hron-fsi1, taylor-bar-impact | `cell: hex` or `cell: quad`: a grid in r and θ, Blossom-Quad, or a swept or mapped hex grid | [§9.1](skein_spec.md#91-mesh), [§9.2](skein_spec.md#92-mesh-generation) |
| Plane strain elements (`CPE4`, `CPE6`, `CPE8`, `CPE8R`) | submodel-kirsch-hole, fracture-centre-crack, contact-cattaneo, hyperelastic-tube, topology-mbb-beam | header `2d`, `cell: tri` or `cell: quad`, `assume plane_strain(th)` | [§3](skein_spec.md#3-grammar), [§9.1](skein_spec.md#91-mesh), [§11.1](skein_spec.md#111-assumption-table) |
| Plane stress elements (`CPS6`) | shape-hole-biaxial | `assume plane_stress(th)` | [§11.1](skein_spec.md#111-assumption-table), [§11.2](skein_spec.md#112-how-plane_stress-reduces-a-3d-law) |
| Quarter-point cells at a crack tip | fracture-centre-crack | `crack_tip tip: size` | [§9.1](skein_spec.md#91-mesh) |
| Shell elements (`S8R`, `S6`), a composite shell section | scordelis-lo-roof, laminate-pagano | the `reissner_mindlin` template on a `Surface`, with a `layup` table for a composite section | [§4.3](skein_spec.md#43-shells-beams-and-laminates), [§8.3](skein_spec.md#83-lower-dimension-entities) |
| First-order edge elements | grid-waveguide-te10, grid-waveguide-dielectric, grid-waveguide-step, grid-cavity-te101, grid-coax-tem, team7-eddy-current, dipole-king-middleton | `method p: cg { order: 1, family: nedelec }` | [§13](skein_spec.md#13-methods) |
| A tree gauge of edge elements | team7-eddy-current | `method p: cg { family: nedelec, gauge: tree }` | [§13](skein_spec.md#13-methods) |
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
| `*EXPANSION` with a reference temperature | warpage-2p5d-frame, warpage-thick-copper, grid-bimetal-strip, grid-constrained-block, grid-heated-block-stress-ccx, grid-plate-hole-stress-ccx | `alpha: Expansion` and `strain(u) - alpha*(T - T0)*I`. `skein run` solves it with constant temperatures on a union of boxes of hex cells, with one material on each region (`src/codegen/codegen.goose`): [warpage_thick_copper.skn](../examples/femahoy/warpage_thick_copper.skn) passes its `expect`. | [§10.4](skein_spec.md#104-coupling) |
| `*PLASTIC` (J2, linear hardening) | plasticity-sphere, taylor-bar-impact | `state`, `plastic { yield: von_mises(sigma) - ..., strain, hardening }` | [§10.6](skein_spec.md#106-inelastic-materials), gap [G22](#gaps-that-change-the-discretization-or-the-physics) |
| `*HYPERELASTIC` (Ogden) | hyperelastic-tube | a strain energy `W` and `S = 2*diff(W, C)` | [§10.6](skein_spec.md#106-inelastic-materials), [§6.1](skein_spec.md#61-dimensions) |
| St. Venant-Kirchhoff under `NLGEOM` | turek-hron-fsi1, buckling-column | `assume finite_strain` and `hooke` of the Green-Lagrange strain | [§11.1](skein_spec.md#111-assumption-table) |
| `*CONDUCTIVITY`, `*SPECIFIC HEAT`, grid heat materials | heat-fin, grid-cooling-slab and the other heat cases | `k: Conductivity`, `cp: SpecificHeat`, `q = fourier(k, T)` | [§6.2](skein_spec.md#62-named-dimension-types), [§10.2](skein_spec.md#102-material) |
| A conductivity table of the temperature | grid-radiating-fin-ccx, grid-conductivity-slab | `table` and `interp(..., T)` in a material property | [§5.6](skein_spec.md#56-tables), [§10.2](skein_spec.md#102-material), Appendix [A](skein_spec.md#appendix-a-built-in-names) |
| Electric conductivity, permittivity, permeability, loss | grid-resistor-bar, grid-coax-capacitor, grid-waveguide-dielectric, line-*, team7-eddy-current | properties with dimension expressions, for example `A/(V*m)` | [§6.2](skein_spec.md#62-named-dimension-types), [§6.3](skein_spec.md#63-units), [§10.2](skein_spec.md#102-material) |
| Viscosity, density, Prandtl number of a fluid | cavity-devahldavis, cht-heated-channel, naca0012-ladson, turek-hron-fsi1 | `mu: DynViscosity`, `rho: Density`, `k`, `cp` | [§10.2](skein_spec.md#102-material) |
| Boussinesq buoyancy | cavity-devahldavis | the body force `rho*(1 - beta*(T - T0))*g` in the momentum equation | [§10.3](skein_spec.md#103-physics), [§11.1](skein_spec.md#111-assumption-table) |
| Perfect gas | wedge-oblique-shock | `gamma: Real` and a gas constant `let` | [§10.2](skein_spec.md#102-material) |
| Two phases with a volume fraction (VOF) | dam-break-vof | mixture laws of `alpha` in physics `let` items | [§10.3](skein_spec.md#103-physics) |
| Tait equation of state | dam-break-martin-moyce | `p = B*((1/J)^7 - 1)` | [§10.3](skein_spec.md#103-physics) |
| A material on each element set or grid part | most multi-material cases | `material a on r1, b on r2` | [§10.3](skein_spec.md#103-physics) |
| A homogenized material of a periodic cell | grid-laminate-cell, grid-fibre-composite-cell | one physics for each load case, `symmetry s: bloch { k: 0 }`, `mean` in `expect`, a mean strain as a global unknown | [§12.6](skein_spec.md#126-symmetry-groups), [§14.5](skein_spec.md#145-output-and-report), [§10.5](skein_spec.md#105-odes-and-global-unknowns) |
| Rigid body mass and inertia | slider-crank | `let` constants in the kinetic energy | [§10.5](skein_spec.md#105-odes-and-global-unknowns) |

### Steps

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*STATIC` | most structural cases | `step s: static` | [§14.2](skein_spec.md#142-procedures) |
| `*STATIC` with increments and `NLGEOM` | hyperelastic-tube, plasticity-sphere, contact-cattaneo | `static(ramp: n)` and `assume finite_strain` | [§14.2](skein_spec.md#142-procedures), [§11.1](skein_spec.md#111-assumption-table) |
| Two or more steps, with the boundary conditions kept | cantilever-2step, contact-hertz, contact-cattaneo, hyperelastic-tube | steps in the order of the file, `release r` | [§14.1](skein_spec.md#141-form), [§12.2](skein_spec.md#122-where-bc-may-stand) |
| `*FREQUENCY`, a modes step of the grid solver | frequency-cantilever, grid-frequency-cantilever, grid-free-block-modes-ccx | `modal(count: n)` | [§14.2](skein_spec.md#142-procedures) |
| Modes near a frequency | grid-cavity-te101 | `modal(count: 1, shift: f)` | [§14.2](skein_spec.md#142-procedures) |
| `*BUCKLE` after a preload | buckling-column | `static` then `buckling(count: n)` | [§14.2](skein_spec.md#142-procedures) |
| `*DYNAMIC, EXPLICIT` | dynamics-bar-step, taylor-bar-impact | `transient(end, dt: stable, scheme: central_difference)` or `ssprk3` | [§14.2](skein_spec.md#142-procedures) |
| `*HEAT TRANSFER, STEADY STATE`, a steady grid heat step | heat-fin and the grid heat cases | `static` | [§14.2](skein_spec.md#142-procedures) |
| A transient heat step (backward Euler) | grid-cooling-slab | `transient(end, dt, scheme: backward_euler)` | [§14.2](skein_spec.md#142-procedures) |
| `*MODEL CHANGE, REMOVE` and `ADD=STRAIN FREE` | birth-death-column | `death r` and `birth r` | [§14.3](skein_spec.md#143-birth-and-death) |
| `*STEADY STATE DYNAMICS` with modal damping | random-vibration-cantilever | `harmonic(freq: [...], damping: 0.02)` | [§14.2](skein_spec.md#142-procedures) |
| One frequency of an electromagnetic or acoustic problem | grid-waveguide-*, acoustics-pulsating-sphere, dipole-king-middleton, optics-mie-cylinder, team7-eddy-current, line-wire-skin, line-bars-getdp | `harmonic(freq: [f])` with `Complex` unknowns | [§5.1](skein_spec.md#51-scalar-types), [§14.2](skein_spec.md#142-procedures) |
| Many frequencies in one step | line-wire-skin, line-bars-getdp, random-vibration-cantilever | `harmonic(freq: [f1, f2])` | [§14.2](skein_spec.md#142-procedures) |
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
| Case folders for each parameter value | cavity-devahldavis, naca0012-ladson, naca23012-flap, dipole-king-middleton, optics-mie-cylinder | `param`, and `study s: grid` with an `expect` for each case | [§20.1](skein_spec.md#201-param), [§20.2](skein_spec.md#202-study), [§15.3](skein_spec.md#153-expect) |
| Shape optimization (Nelder-Mead) | shape-hole-biaxial | `study s: minimize(objective, method: nelder_mead, max_cases, xtol, step)` | [§20.2](skein_spec.md#202-study) |
| Topology optimization (SIMP) | topology-mbb-beam, topology-mbb-beam-cad | `study s: topology(objective, region, volume, penalty, filter)` | [§20.2](skein_spec.md#202-study) |
| A least-squares curve fit | fit-nist-pontius | the normal equations: `sum(expr)` over the rows of a table, and a physics of global unknowns | [§5.6](skein_spec.md#56-tables), [§10.5](skein_spec.md#105-odes-and-global-unknowns) |

### Loads

| Femahoy feature | Cases | Skein construct | Spec |
|---|---|---|---|
| `*DSLOAD` and `*DLOAD` pressure | bar-pressure, buckling-column, plasticity-sphere, dynamics-bar-step, fracture-centre-crack | `bc r: sigma*n = -p*n` | [§12.1](skein_spec.md#121-form) |
| A traction vector on a face | grid-cantilever, grid-kirsch-hole, grid-round-cantilever | `bc r: sigma*n = (...)` | [§12.1](skein_spec.md#121-form) |
| A traction that varies along a face | laminate-pagano, contact-cattaneo | `sigma*n` with a function of `x`, `y` | [§7.6](skein_spec.md#76-coordinates-and-time), [§12.1](skein_spec.md#121-form) |
| `*CLOAD` at a node | cantilever-2step, topology-mbb-beam, contact-cattaneo | `bc p: sigma*n = P` on a `Point` | [§12.1](skein_spec.md#121-form) |
| `*DLOAD, GRAV` | birth-death-column, scordelis-lo-roof | the template argument `b: rho*(0, 0, -g)` | [§4.2](skein_spec.md#42-the-standard-library), [§4.3](skein_spec.md#43-shells-beams-and-laminates) |
| Gravity on a fluid | cavity-devahldavis, dam-break-vof, dam-break-martin-moyce | a body force term in `eq` | [§10.3](skein_spec.md#103-physics) |
| `*TEMPERATURE` and a uniform temperature of the grid solver | warpage-2p5d-frame, warpage-thick-copper, grid-bimetal-strip, grid-constrained-block | a constant `T` in the thermal strain, for example `alpha*(T_use - T_free)*I` with `degC` values. `skein run` solves it in the box lane (`src/codegen/codegen.goose`). | [§10.4](skein_spec.md#104-coupling) |
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
| `*BASE MOTION` with a PSD amplitude | random-vibration-cantilever | the load of the base acceleration, and the response `g - omega^2*u` | [§14.2](skein_spec.md#142-procedures) |
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
| `*EQUATION`: a face stays straight and can turn | contact-cattaneo | `bc top: u = rigid(c)`, and the force and the `moment` on `c` | [§12.1](skein_spec.md#121-form) |
| `*CONTACT PAIR`, `*SURFACE INTERACTION`, `*FRICTION` | contact-hertz, contact-cattaneo | `contact c: master, slave { friction, enforce: penalty, stiffness, stick_stiffness }` | [§12.4](skein_spec.md#124-contact), gap [G23](#gaps-that-change-the-discretization-or-the-physics) (closed) |
| A periodic cell | grid-laminate-cell, grid-fibre-composite-cell | `symmetry s: bloch { low, high, period, k: 0 }` | [§12.6](skein_spec.md#126-symmetry-groups) |
| A free body with no support | grid-free-block-modes-ccx | no `bc`, `allow rigid_body: "reason"` | [checks.md](checks.md) |
| No slip, inlet velocity, outlet pressure, slip wall, symmetry plane | cavity-devahldavis, cht-heated-channel, naca0012-ladson, wedge-oblique-shock, turek-hron-fsi1 | `bc r: v = ...`, `sigma*n = -p*n`, `dot(v, n) = 0` | [§12.1](skein_spec.md#121-form) |
| A pressure reference (`pRefCell`) | cavity-devahldavis | `bc p0: p = 0` on a `Point` | [§8.3](skein_spec.md#83-lower-dimension-entities), [§12.3](skein_spec.md#123-completeness) |
| The coupled temperature of 2 regions | cht-heated-channel | `bc r: Tf = Ts` and equal fluxes | [§12.1](skein_spec.md#121-form) |
| A perfect electric conductor on a flat wall | grid-waveguide-*, grid-cavity-te101 | the tangential components: `bc r: E.y = 0` | [§12.1](skein_spec.md#121-form) |
| A perfect electric conductor on a curved wall | grid-coax-tem | the tangential part: `bc r: cross(n, ef) = 0` | [§12.1](skein_spec.md#121-form) |
| A surface impedance wall | grid-cavity-te101 | a Robin `bc` with the impedance `(1 + im)*sqrt(omega*mu0/(2*sigma))` | [§12.1](skein_spec.md#121-form), [§14.2](skein_spec.md#142-procedures) |
| An open boundary (PML, a map to infinity) | team7-eddy-current, dipole-king-middleton, optics-mie-cylinder | the method setting `pml: { region: layer, thickness: h }` | [§13](skein_spec.md#13-methods) |
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
| `CSTR` contact pressure | contact-hertz, contact-cattaneo | `max(-sigma.zz, face)`, with σ at a node the mean of its cells | [§14.5](skein_spec.md#145-output-and-report) |
| `FREQ` and `BUCKLE` | frequency-cantilever, buckling-column, grid-frequency-cantilever, grid-free-block-modes-ccx, grid-cavity-te101 | `eigenvalue(k)` in an `expect` | [§14.5](skein_spec.md#145-output-and-report) |
| `A_RMS` | random-vibration-cantilever | `rms(f, psd)` in an `expect` | [§14.5](skein_spec.md#145-output-and-report) |
| Heat flow through a face (`energy.txt`) | grid-pipe-wall, grid-radiation-slab | `integral(dot(q, n), r)` | [§14.5](skein_spec.md#145-output-and-report) |
| Resistance, conductance, capacitance of terminals | grid-resistor-bar, grid-holm-spot, grid-coax-capacitor | a ratio of `integral(dot(J, n), r)` and the voltage | [§14.5](skein_spec.md#145-output-and-report) |
| S-parameters, impedance, R and L of a line | grid-waveguide-*, line-*, dipole-king-middleton | projection integrals, and `re`, `imag`, `abs` and `arg` in an `expect` | [§14.5](skein_spec.md#145-output-and-report), [§15.3](skein_spec.md#153-expect), [Appendix A](skein_spec.md#appendix-a-built-in-names) |
| Lift, drag and moment coefficients (`forceCoeffs`) | naca0012-ladson, naca23012-flap, turek-hron-fsi1 | `report cl = integral(dot(sigma*n, e_l), wall)/...` | [§14.5](skein_spec.md#145-output-and-report) |
| Nusselt number, peak velocities, mean temperatures | cavity-devahldavis, cht-heated-channel | `report` with `mean`, `max`, `integral` | [§14.5](skein_spec.md#145-output-and-report) |
| A curvature from a fit of the top face | grid-bimetal-strip | `integral` moments of `u.z` | [§14.5](skein_spec.md#145-output-and-report) |
| A stress intensity factor | fracture-centre-crack | `K_I` in `MPa*mm^0.5` | [§6.1](skein_spec.md#61-dimensions) |
| The time at which a front gets to a position | dam-break-vof, dam-break-martin-moyce | `when(front, value)` in an `expect` | [§14.5](skein_spec.md#145-output-and-report) |
| A shock angle | wedge-oblique-shock | `atan` of the slope | [Appendix A](skein_spec.md#appendix-a-built-in-names) |
| Directivity and efficiency of an antenna | dipole-king-middleton | the report function `far_field(E, H, s, d)`, and the power `integral(re(dot(cross(E, conj(H)), n))/2, s)` | [§14.5](skein_spec.md#145-output-and-report) |
| A check after an earlier step | hyperelastic-tube | `expect` in the step | [§15.3](skein_spec.md#153-expect) |
| The compliance of a design | topology-mbb-beam | `report compliance = -(u.y at p)*F`, the objective of `topology` | [§20.2](skein_spec.md#202-study) |
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
| `validation/cases/fit-nist-pontius` | [fit_nist_pontius.skn](../examples/femahoy/fit_nist_pontius.skn) |
| `validation/cases/fracture-centre-crack` | [fracture_centre_crack.skn](../examples/femahoy/fracture_centre_crack.skn) |
| `validation/cases/frequency-cantilever` | [frequency_cantilever.skn](../examples/femahoy/frequency_cantilever.skn) |
| `validation/cases/grid-bimetal-strip` | [grid_bimetal_strip.skn](../examples/femahoy/grid_bimetal_strip.skn) |
| `validation/cases/grid-cantilever` | [grid_cantilever.skn](../examples/femahoy/grid_cantilever.skn) |
| `validation/cases/grid-cavity-te101` | [grid_cavity_te101.skn](../examples/femahoy/grid_cavity_te101.skn) |
| `validation/cases/grid-coax-capacitor` | [grid_coax_capacitor.skn](../examples/femahoy/grid_coax_capacitor.skn) |
| `validation/cases/grid-coax-tem` | [grid_coax_tem.skn](../examples/femahoy/grid_coax_tem.skn) |
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
| `validation/cases/random-vibration-cantilever` | [random_vibration_cantilever.skn](../examples/femahoy/random_vibration_cantilever.skn) |
| `validation/cases/scordelis-lo-roof` | [scordelis_lo_roof.skn](../examples/femahoy/scordelis_lo_roof.skn) |
| `validation/cases/shape-hole-biaxial` | [shape_hole_biaxial.skn](../examples/femahoy/shape_hole_biaxial.skn) |
| `validation/cases/slider-crank` | [slider_crank.skn](../examples/femahoy/slider_crank.skn) |
| `validation/cases/submodel-kirsch-hole` | [submodel_kirsch_hole.skn](../examples/femahoy/submodel_kirsch_hole.skn), [submodel_kirsch_hole_global.skn](../examples/femahoy/submodel_kirsch_hole_global.skn) |
| `validation/cases/taylor-bar-impact` | [taylor_bar_impact.skn](../examples/femahoy/taylor_bar_impact.skn) |
| `validation/cases/team7-eddy-current` | [team7_eddy_current.skn](../examples/femahoy/team7_eddy_current.skn) |
| `validation/cases/topology-mbb-beam` | [topology_mbb_beam.skn](../examples/femahoy/topology_mbb_beam.skn), study `shape` |
| `validation/cases/topology-mbb-beam-cad` | [topology_mbb_beam.skn](../examples/femahoy/topology_mbb_beam.skn), study `solid` |
| `validation/cases/turek-hron-fsi1` | [turek_hron_fsi1.skn](../examples/femahoy/turek_hron_fsi1.skn) |
| `validation/cases/warpage-2p5d-frame` | [warpage_2p5d_frame.skn](../examples/femahoy/warpage_2p5d_frame.skn) |
| `validation/cases/warpage-thick-copper` | [warpage_thick_copper.skn](../examples/femahoy/warpage_thick_copper.skn) |
| `validation/cases/wedge-oblique-shock` | [wedge_oblique_shock.skn](../examples/femahoy/wedge_oblique_shock.skn) |
