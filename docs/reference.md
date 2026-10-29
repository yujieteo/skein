# The Skein Language Reference

<details>
<summary>Table of contents</summary>

- [1. Lexical structure](#1-lexical-structure)
- [2. The header](#2-the-header)
- [3. Declarations](#3-declarations)
- [4. Types](#4-types)
- [5. Units](#5-units)
- [6. Expressions](#6-expressions)
- [7. Geometry](#7-geometry)
- [8. Meshes and point clouds](#8-meshes-and-point-clouds)
- [9. Materials](#9-materials)
- [10. Fields and physics](#10-fields-and-physics)
- [11. Modelling assumptions](#11-modelling-assumptions)
- [12. Boundary conditions](#12-boundary-conditions)
- [13. Methods](#13-methods)
- [14. Steps](#14-steps)
- [15. Birth and death](#15-birth-and-death)
- [16. Restarts](#16-restarts)
- [17. Results and submodels](#17-results-and-submodels)
- [18. Tables](#18-tables)
- [19. Expectations](#19-expectations)
- [20. Turning off a check](#20-turning-off-a-check)
- [21. Inelastic materials](#21-inelastic-materials)
- [22. Coupled physics](#22-coupled-physics)
- [23. Contact](#23-contact)
- [24. Mesh adaptivity](#24-mesh-adaptivity)
- [25. ODEs](#25-odes)
- [26. Parameters and studies](#26-parameters-and-studies)
- [27. Scales and dimensionless groups](#27-scales-and-dimensionless-groups)
- [28. Manufacture and tolerances](#28-manufacture-and-tolerances)
- [29. Library files](#29-library-files)
- [30. Symmetry groups](#30-symmetry-groups)
- [Appendix A. Keywords](#appendix-a-keywords)

</details>

This is the reference for the Skein language: each construct, topic by
topic, with an example. Each section ends with a pointer to the rules in
[the specification](skein_spec.md). The examples are fragments. A complete
model needs the header and the other parts. The [samples](../samples/README.md)
are complete.

---

## 1. Lexical structure

A newline ends a statement. A line that ends with an operator, a comma or
an open bracket continues on the next line. `//` starts a comment. `///`
starts a documentation comment for the next declaration.

```skein
/// Young's modulus of the shaft, from the supplier's data sheet.
let E: Pressure = 205 GPa   // at 20 degC
let sigma: Tensor2<Pressure> = 2*mu*strain(u) +
    lam*tr(strain(u))*I
```

Spec: [§2](skein_spec.md#2-lexical-structure).

---

## 2. The header

Each model starts with one line: the name, the space and the unit system.

```skein
model bracket: 3d, mm_N_MPa_s
```

The spaces are `1d`, `2d`, `3d` and `axi` (axisymmetric). The unit systems
are `SI`, `mm_N_MPa_s`, `mm_kg_ms` and `nondim`, or a record:

```skein
model bracket: 3d, { length: mm, mass: t, time: s, temperature: K }
```

`use` lines for library files come after the header (§29).

Spec: [§3](skein_spec.md#3-grammar), [§6.5](skein_spec.md#65-unit-systems).

---

## 3. Declarations

`let` gives a name to a value. The type is required:

```skein
let L: Length = 200 mm
let modes: Int = 10
let plate: Surface = rect((0, 0) mm, (L, 50 mm))
```

`fn` declares a function of one expression. Its parameters and its result
have types:

```skein
fn lame_mu(E: Pressure, nu: Real) -> Pressure = E/(2*(1 + nu))
```

A name must be declared above its first use. A model has one namespace, so
two declarations with one name are an error. This also holds for the names
of the standard library, of `use` files and of the built-in functions.

Spec: [§3.1](skein_spec.md#31-order-of-items), [§4](skein_spec.md#4-names-and-scope).

---

## 4. Types

| Kind | Types |
|---|---|
| scalar quantities | `Length`, `Pressure`, … or a dimension expression such as `Force/Length` |
| other scalars | `Real`, `Angle`, `Int`, `Bool`, `Str` |
| number systems | `Complex<D>`, `Quaternion<D>`, `Dual<D>`. The list is closed. |
| tensors | `Vector<D>`, `Tensor2<D>`, `Tensor4<D>` |
| geometry | `Solid`, `Surface`, `Curve`, `Point` |
| model | `Mesh`, `Points`, `Results` |

```skein
let g: Vector<Acceleration> = (0, 0, -9.81) m/s^2
let C: Tensor4<Pressure> = elasticity(E, nu)
let k_spring: Force/Length = 2 kN/mm
```

A tensor has the space dimension: a `Vector` has 3 components in a `3d` model.
Components are `v.x`, `v.y`, `v.z` and `A.xx`, `A.xy`, ….

Spec: [§5](skein_spec.md#5-types).

---

## 5. Units

A literal is a number and a unit. A vector literal takes one unit after the
parenthesis:

```skein
let rho: Density = 1000 kg/m^3
let h: HeatTransfer = 25 W/(m^2*K)
let p0: Vector<Length> = (40, 0, 20) mm
let T0: Temperature = 20 degC
let T_sink: Temperature = -5 degC  // 268.15 K
```

Units have SI prefixes: `mm`, `kN`, `MPa`, `GPa`, `ms`. `degC` is valid only
alone on a temperature. A bare number with no `.` and no exponent is `Int`.
Each other bare number is `Real`. `Real` accepts `Int`. `0` is zero of any
type.

Spec: [§6](skein_spec.md#6-dimensions-and-units).

---

## 6. Expressions

The operators, from low to high precedence: `|| |`, `&& &`, comparisons,
`+ -`, `* /`, unary `- !`, `^`, then `.`, calls and `at`.

```skein
let eps: Tensor2<Real> = sym(grad(u))
let vm: Pressure = von_mises(sigma)
let dT: Temperature = T - T0
let w_tip: Length = u.z at (L, W/2, H/2)
```

`*` between two tensors is the contraction: `A*v` is a matrix-vector
product. `*` between two vectors is an error. Use `dot`, `cross` or `outer`.

The differential operators are `grad`, `div`, `curl`, `lap`, `ddt` and
`ddt2`. `x`, `y`, `z` and `t` are the coordinates and the time.

Spec: [§7](skein_spec.md#7-expressions-and-operators).

---

## 7. Geometry

Shapes: `box`, `cylinder`, `sphere`, `cone`, `torus`, `rect`, `circle`,
`polygon`, `extrude`, `revolve`, `sweep`, `line`, `arc`, `spline`, `point`.
Operators: `|` union, `&` intersection, `-` difference. Modifiers: `fillet`,
`chamfer`, `move`, `rotate`, `scale`.

```skein
let profile: Surface = polygon((0, 0) mm, (60, 0) mm, (60, 8) mm, (8, 8) mm, (8, 40) mm, (0, 40) mm)
let angle: Solid = fillet(extrude(profile, 100 mm), 2 mm)
let bolt: Solid = cylinder((30, 0, 50) mm, (0, 1, 0), 4 mm, 8 mm)
let part: Solid = angle - bolt
```

`boundary(g, s)` selects a part of the boundary by a condition or by a
shape:

```skein
let base: Surface = boundary(part, y == 0 mm)
let hole: Surface = boundary(part, bolt)
let all: Surface = boundary(part)
```

Spec: [§8](skein_spec.md#8-geometry-and-regions).

---

## 8. Meshes and point clouds

```skein
mesh m on part {
    size: 3 mm
    cell: tet
    grade: 1.3
    refine hole: 0.5 mm
}
```

For flow near a wall, add boundary layer cells:

```skein
mesh m on duct {
    size: 2 mm
    layers wall: { first: 0.05 mm, ratio: 1.2, count: 10 }
}
```

A `points` block makes a point cloud for a meshless method:

```skein
points c on plate {
    spacing: 2 mm
    refine hole: 0.4 mm
}
```

Spec: [§9](skein_spec.md#9-meshes-and-point-clouds).

---

## 9. Materials

A material is a record of typed properties. A property may use the
properties above it:

```skein
material water {
    rho: Density = 998 kg/m^3
    mu: DynViscosity = 1.0e-3 Pa*s
    nu_k: KinViscosity = mu/rho
}
```

A physics block assigns materials to regions:

```skein
material clay on ground, fill on lift1 | lift2
```

Spec: [§10.2](skein_spec.md#102-material).

---

## 10. Fields and physics

A `physics` block declares its unknowns in `solves` and states its
equations in strong form:

```skein
physics thermal on part solves T: Temperature = 293 K {
    material steel
    let q: Vector<HeatFlux> = -k*grad(T)
    eq rho*c*ddt(T) + div(q) = Q
}
```

The value after `=` is the initial value. One `eq` goes with each unknown
in `solves`, in the same order. The `let`
lines inside the block are visible in the rest of the model, for example in
a `bc` or an `output`.

A physics can take its statements from a standard-library template:

```skein
physics solid on beam: linear_elastic solves u: Vector<Length> {
    material steel
    assume plane_strain
}
```

The block may add `material`, `assume` and `let`, but not `eq`. `skein help
std` prints the source of each template.

Spec: [§10](skein_spec.md#10-fields-materials-and-physics), [§4.2](skein_spec.md#42-the-standard-library).

---

## 11. Modelling assumptions

```skein
assume small_strain, isotropic, plane_strain
assume finite_strain, isotropic
assume incompressible, newtonian, laminar
assume boussinesq(293 K, 2.1e-4/K)
```

| Assumption | Short meaning |
|---|---|
| `small_strain` / `finite_strain` | the meaning of `strain(u)` and of the reference shape |
| `isotropic` / `orthotropic(R)` | the symmetry of `hooke` and `fourier` |
| `plane_stress(t)` / `plane_strain(t)` | 2D reduction of a 3D law |
| `incompressible` / `compressible` | the continuity equation |
| `newtonian` | the stress is linear in the rate of strain |
| `laminar` / `turbulent(model)` | no turbulence model, or a closure |
| `boussinesq(T0, beta)` | buoyancy with constant density elsewhere |
| `isothermal` | no energy equation |

Spec: [§11](skein_spec.md#11-modelling-assumptions).

---

## 12. Boundary conditions

```skein
bc root: u = 0                        // essential
bc sym: u.y = 0                       // essential, one component
bc tip: sigma*n = (0, -2, 0) MPa      // natural: traction
bc wall: dot(q, n) = h*(T - T_inf)    // natural: convection (Robin)
bc inlet: v = ((2 m/s)*(1 - (y/R)^2), 0)
bc corner: p = 0                      // a point
bc load_pt: sigma*n = (0, 0, -1) kN   // a point force
```

A `bc` outside the steps applies in all steps. A `bc` inside a step applies
from that step on. A boundary with no `bc` is free.

Spec: [§12](skein_spec.md#12-boundary-conditions-and-loads).

---

## 13. Methods

```skein
method solid: cg { mesh: m, order: 2 }
method flow: cg { mesh: m, order: { v: 2, p: 1 } }
method flow: fv { mesh: m, gradient: least_squares, flux: linear_upwind, pv_coupling: rhie_chow }
method wave: dg { mesh: m, order: 3, flux: upwind }
method solid: rkpm { points: c, basis: linear, support: 2.5, integrate: scni, bc: nitsche }
method thermal: gfd { points: c, stencil: 12 }
```

A block may hold its settings on one line, with commas, or on separate
lines. `stabilize` adds residual-based terms for flow: `supg`, `supg_pspg`,
`supg_pspg_lsic` or `vms`. `limiter` (for `fv`) and `shock_capturing:
entropy_viscosity` (for `dg`) capture shocks.

A physics with no `method` item gets `cg` on its mesh, `order: 2`,
Taylor-Hood cells for incompressible flow, and `supg` when the equations
have convection.

Spec: [§13](skein_spec.md#13-methods), [§17](skein_spec.md#17-lowering-by-method).

---

## 14. Steps

```skein
step preload: static(ramp: 5) {
    bc bolt: sigma*n = -50 MPa*n
    solve { nonlinear: newton, rtol: 1e-8 }
    output u, sigma
}

step free: modal(count: 12) {
    output u
    report f1 = eigenvalue(1)
}

step heat_up: transient(end: 1 h, dt: 30 s, scheme: bdf2) {
    bc wall: T = 400 K
    output T every 4
    report max(T, part)
}
```

The procedures are `static`, `transient`, `modal`, `buckling` and
`harmonic(freq: [f])`, which takes a list and uses `Complex` unknowns.
`transient(adapt: { rtol: 1e-4, atol: 1e-8 })` with `scheme: bdf` changes
the time step and the order to keep the error below the limits. Steps run in the order of the file.

`re`, `imag`, `abs`, `arg` and `conj` take the parts of a `Complex` value.
`eigenvalue(k)` is the natural frequency f = ω/(2π) in Hz of mode `k` of a
`modal` step, or the load factor of a `buckling` step.

Spec: [§14](skein_spec.md#14-steps).

---

## 15. Birth and death

```skein
step deposit_2: transient(end: 20 s, dt: 0.5 s) {
    birth layer2
}

step burn_out: transient(end: 10 s, dt: 0.1 s) {
    death liner where T > 900 K
}
```

A born region starts strain-free on the current shape. A dead region leaves
the system. `where` selects the cells by a condition on the fields.

Spec: [§14.3](skein_spec.md#143-birth-and-death).

---

## 16. Restarts

Each run saves the state at the end of each step. To start again from a
step:

```
skein run embankment.skn --from lift_3
```

The steps above `lift_3` must not change since the saved run.

Spec: [§14.7](skein_spec.md#147-restarts).

---

## 17. Results and submodels

```skein
let global: Results = results("bracket.skn", step: "load")

step load: static {
    bc cut: u = global.u
}
```

`global.u` has the type of `u` in `bracket.skn`. The solver interpolates
the global field on the cut.

Spec: [§15](skein_spec.md#15-results-submodels-and-expectations).

---

## 18. Tables

```skein
table hardening{eps_p: Real, sy: Pressure in MPa}:
    0.00, 250
    0.02, 290
    0.10, 340
```

The rows are TOON rows. `in MPa` gives the unit of each bare number in the
column. `hardening.sy` is the column. A table also works as a piecewise
linear function: `interp(hardening.eps_p, hardening.sy, alpha)`.

Spec: [§5.6](skein_spec.md#56-tables).

---

## 19. Expectations

```skein
expect tip: u.z at (L, 0 mm, 0 mm) == -P*L^3/(3*E*I_z) within 2%
expect mms: T == T_exact within 1e-3 K
```

An expectation is a test on the results. A failed one gives exit status 1.

Spec: [§15.3](skein_spec.md#153-expect).

---

## 20. Turning off a check

```skein
allow volumetric_lock: "nu = 0.49 is for a sensitivity study; the result is not used"
```

`allow` takes the id of a lint rule or a verification check, and a reason.
Errors cannot be turned off.

Spec: [checks.md §1](checks.md#1-how-diagnostics-work).

---

## 21. Inelastic materials

```skein
physics solid on bar solves u: Vector<Length> {
    material steel
    assume small_strain, isotropic
    state eps_p: Tensor2<Real> = 0
    state alpha: Real = 0
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) - eps_p)
    plastic {
        yield: von_mises(sigma) - (sy + Hk*alpha)
        strain: eps_p
        hardening: alpha
    }
    eq div(sigma) = 0
}
```

The compiler makes the return map and the consistent tangent. For a
hyperelastic law, write the energy and use `diff` (spec §10.6).

Spec: [§10.6](skein_spec.md#106-inelastic-materials).

---

## 22. Coupled physics

Two physics blocks couple when one uses the field of the other:

```skein
physics solid on part solves u: Vector<Length> {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) - alpha*(T - T0)*I)
    eq div(sigma) = 0
}

method solid: cg { mesh: m, order: 2, couple: staggered }
```

For fluid-structure interaction, an `interface` joins the two boundaries:

```skein
interface fsi: wet_fluid, wet_solid {
    mesh_motion: harmonic
}
```

Spec: [§10.4](skein_spec.md#104-coupling), [§10.7](skein_spec.md#107-fsi-interfaces).

---

## 23. Contact

```skein
contact grip: jaw_face, part_face {
    friction: 0.2
    enforce: mortar
}
```

Spec: [§12.4](skein_spec.md#124-contact).

---

## 24. Mesh adaptivity

```skein
step load: static {
    adapt { estimator: zz, target: 1%, max_levels: 3 }
}
```

Spec: [§14.6](skein_spec.md#146-mesh-adaptivity).

---

## 25. ODEs

A physics with no region solves global unknowns. It is a system of ODEs:

```skein
physics oscillator solves q: Length = 1 mm {
    eq m*ddt2(q) + c*ddt(q) + k*q = 0
}

method oscillator: ode {}

step ring: transient(end: 2 s, dt: 1 ms, scheme: generalized_alpha(0.8)) {
    report q
}
```

Spec: [§10.5](skein_spec.md#105-odes-and-global-unknowns).

---

## 26. Parameters and studies

A `param` is a `let` that a study or `--set` can change. `in lo..hi` gives
its range:

```skein
param Re: Real = 20 in 1..200

study re_sweep: grid {
    Re: [10, 20, 40]
}
```

The procedures are `grid`, `zip`, `sample` and `minimize`. Each case runs
all the steps. Cases run at the same time on worker threads. A named
`report cd = …` gives a column of the results table.

```
skein run model.skn --study re_sweep --jobs 3
skein run model.skn --set Re=50
```

Spec: [§20](skein_spec.md#20-parameters-and-studies).

---

## 27. Scales and dimensionless groups

```skein
let mu: DynViscosity = rho*U*D/Re
scales D, U, rho
```

`scales` names the reference values of the scaled equations. Each model
with dimensions has one. The source stays dimensional. `skein check --groups` prints the dimensionless groups
of the model, with the known names such as `Re`, `Pr` and `Ma`.

Spec: [§21](skein_spec.md#21-dimensional-analysis-and-scaling).

---

## 28. Manufacture and tolerances

```skein
param t_wall: Length = 2 mm tol 0.1 mm in 1 mm..5 mm

manufacture housing: injection_mould {
    pull: (0, 0, 1)
    draft: 1 deg
    min_wall: 1.2 mm
}

tolerance fit: gap in 0.1 mm..0.7 mm { method: worst_case }
```

`skein verify --dfm` checks the part against its process: wall thickness,
draft, undercuts, inner radii, holes, overhangs and bends. A `tolerance`
item checks a stack-up of `tol` values, by `worst_case`, `rss` or
`monte_carlo`. `skein export` writes the part as STL or tessellated STEP.
`shell`, `loft` and `pattern` make the common features of a housing or a
bracket.

Spec: [§22](skein_spec.md#22-manufacture), [§8.1](skein_spec.md#81-distance-functions).

---

## 29. Library files

```skein
use "materials.skn"
```

A library file holds only `material`, `param`, `let` and `fn` items. Its
names enter the namespace of the model with no prefix. A library file
cannot `use` another file. Phase 2.

Spec: [§4.1](skein_spec.md#41-use-phase-2).

---

## 30. Symmetry groups

```skein
symmetry rotor: cyclic { n: 24, axis: (0, 0, 1), low: face_a, high: face_b }
```

A `symmetry` item solves one sector of a body that repeats. The groups are
`cyclic`, `dihedral` and `bloch` (periodic, with a wave vector `k`). The
list is closed.

Spec: [§12.6](skein_spec.md#126-symmetry-groups).

---

## Appendix A. Keywords

```
adapt       allow       assume      at          bc          birth
contact     death       eq          every       expect      false
fn          in          interface   layers      let         manufacture
material    mesh        method      model       on          output
param       physics     plastic     points      refine      release
report      scales      solve       solves      state       step
study       symmetry    table       template    tol         tolerance
true        use         where       within
```

A keyword can be a name only before the `:` of a setting, a named argument
or a record entry: `mesh: m`, `step: "load"`.
