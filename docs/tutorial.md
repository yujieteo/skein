# Skein, by example

<details>
<summary>Table of contents</summary>

- [1. The core idea](#1-the-core-idea)
- [2. Checking and running a model](#2-checking-and-running-a-model)
- [3. Declarations and types](#3-declarations-and-types)
- [4. Units](#4-units)
- [5. Geometry and regions](#5-geometry-and-regions)
- [6. The mesh](#6-the-mesh)
- [7. Physics in strong form](#7-physics-in-strong-form)
- [8. Modelling assumptions](#8-modelling-assumptions)
- [9. Boundary conditions](#9-boundary-conditions)
- [10. Methods](#10-methods)
- [11. Steps](#11-steps)
- [12. Verification by default](#12-verification-by-default)
- [13. Fluids and stabilization](#13-fluids-and-stabilization)
- [14. Birth and death](#14-birth-and-death)
- [15. Submodels](#15-submodels)
- [16. Without a mesh](#16-without-a-mesh)
- [17. Tables and expectations](#17-tables-and-expectations)
- [18. Where to go next](#18-where-to-go-next)

</details>

This tutorial teaches Skein with examples. It assumes that you know a little
about the finite element method and that you can read Rust or a similar
language. After it, read [the language reference](reference.md) for each
construct, and [the specification](skein_spec.md) for the exact rules. The
[samples](../samples/README.md) are complete models.

Skein is a language for numerical simulation. One file holds one model: the
geometry, the mesh, the physics, the modelling assumptions, the numerical
method, the analysis steps and the checks on the results. The compiler
writes [Goose](https://aardappel.github.io/goose/), and Goose writes C for a
native executable or for WebAssembly.

---

## 1. The core idea

This is a complete model of a steel bar under tension:

```skein
model bar
space 3d
units mm_N_MPa_s

let bar: Solid = box((0, 0, 0) mm, (100, 10, 10) mm)
let left: Surface = boundary(bar, x == 0 mm)
let right: Surface = boundary(bar, x == 100 mm)

material steel {
    E: Pressure = 210 GPa
    nu: Real = 0.3
}

mesh m on bar {
    size: 2 mm
}

field u: Vector<Length> on bar

physics solid on bar solves u {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}

method solid: cg {
    mesh: m
    order: 2
}

step pull: static {
    bc left: u = 0
    bc right: sigma*n = (100, 0, 0) MPa
    output u, sigma
}
```

Four things make Skein different from an Abaqus or APDL input deck:

- **The physics is an equation.** `div(sigma) = 0` is the equilibrium
  equation from a textbook. Skein has no element library with a fixed list
  of physics. The compiler derives the weak form from the equation.
- **Every quantity has a dimension.** `E: Pressure = 210 GPa` declares the
  dimension. The literal gives the unit. The checker finds `E + rho` before
  the solve.
- **There is no mesh file.** The model builds the geometry and meshes it.
- **Verification is on.** Each run checks the residual, the equilibrium of
  the reactions and the energy balance, and prints the result.

---

## 2. Checking and running a model

```
skein check bar.skn
skein run bar.skn
```

`check` reads the model and checks the names, the types, the units, the
regions, the assumptions and the boundary conditions. It needs no mesh, so
it is fast. Run it after each edit.

`run` checks the model, meshes it, runs the pre-solve checks, solves each
step and runs the post-solve checks. The output is TOON, a compact table
format:

```toon
verify[5]{stage,id,result,value,limit}:
  pre,rigid_body,pass,0,0
  pre,mesh_quality,pass,"14.1 deg","5 deg"
  post,residual,pass,2.2e-9,1e-7
  post,equilibrium,pass,"4.0e-6 N","1.0e-2 N"
  post,energy_balance,pass,1.1e-11,1e-6
```

`skein view bar.skn` shows the deformed shape and the contours in the
browser.

---

## 3. Declarations and types

Each declaration states its type. Skein does no type inference:

```skein
let L: Length = 200 mm
let n: Int = 10
let ratio: Real = 0.3
let bar: Solid = box((0, 0, 0) mm, (L, 10 mm, 10 mm))
```

The type comes after the name and a colon, as in Rust. An equation or an
expression has no type annotation, so the physics stays short:

```skein
let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
eq div(sigma) = 0
```

The checker computes the type of `div(sigma)`, a `Vector` of
`Pressure/Length`, and compares it with the other side. A `0` with no unit
is zero of any type.

A model has no variables that change. A `let` is a name for a value. A
model also has no loops and no functions. The compiler makes the loops, in
Goose.

---

## 4. Units

The type of a quantity is its dimension. The literal carries the unit:

```skein
let E: Pressure = 210 GPa
let rho: Density = 7.85 g/cm^3
let k: Force/Length = 5 kN/mm
```

Skein stores each value in SI, so `1 m + 1 mm` is valid. A sum of two
dimensions is not:

```skein
let bad: Pressure = E + rho      // error unit_mismatch: Pressure + Density
```

The `units` line in the header chooses the unit system of the solver and of
the output, for example `mm_N_MPa_s`. It does not change the meaning of a
literal, because each literal names its unit.

---

## 5. Geometry and regions

Geometry is built from simple shapes with `|` (union), `&` (intersection)
and `-` (difference):

```skein
let bore: Solid = cylinder((40, 0, 20) mm, (0, 1, 0), 6 mm, 20 mm)
let beam: Solid = box((0, 0, 0) mm, (200, 20, 40) mm) - bore
```

`boundary` selects a part of the boundary, by a condition on the
coordinates or by another shape:

```skein
let root: Surface = boundary(beam, x == 0 mm)
let hole: Surface = boundary(beam, bore)
```

Each shape is a signed distance function. The checker knows the dimension
of each value: `boundary` of a `Solid` is a `Surface`.

---

## 6. The mesh

```skein
mesh m on beam {
    size: 5 mm
    cell: tet
    refine hole: 1 mm
}
```

The built-in mesher makes the mesh from the geometry. It keeps the sharp
edges of the shapes, for example the edge where the hole meets the face. It
runs in the browser too. `refine` gives a smaller size on a region and near
it.

If a feature of the geometry is smaller than the local size, `skein lint
--mesh` says so and gives the size that resolves it.

---

## 7. Physics in strong form

A `physics` block states the equations, the materials and the unknowns:

```skein
field u: Vector<Length> on beam

physics solid on beam solves u {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}
```

- `field` declares an unknown and its region.
- `solves u` pairs the unknowns with the equations, one `eq` for each
  unknown.
- `material steel` gives `E` and `nu` to the equations.

You can write the law yourself. This is the same physics:

```skein
let eps: Tensor2<Real> = strain(u)
let sigma: Tensor2<Pressure> = E/(1 + nu)*(eps + nu/(1 - 2*nu)*tr(eps)*I)
```

Heat conduction has the same shape:

```skein
field T: Temperature on part = 293 K

physics heat on part solves T {
    material steel
    let q: Vector<HeatFlux> = -k*grad(T)
    eq rho*c*ddt(T) + div(q) = 0
}
```

---

## 8. Modelling assumptions

`assume` states the modelling assumptions. Each one has an exact meaning in
the specification:

```skein
assume small_strain, isotropic, plane_stress(2 mm)
```

`plane_stress` tells the compiler that σzz = 0. You still write the 3D law.
The compiler solves it for εzz and makes the 2D law. Thus one law serves
3D, plane strain and plane stress.

Two assumptions that conflict are an error:

```skein
assume small_strain, finite_strain   // error assumption_conflict
```

---

## 9. Boundary conditions

A boundary condition is an equation on a region:

```skein
bc root: u = 0
bc sym: u.x = 0
bc tip: sigma*n = (0, 0, -1) MPa
```

The checker sorts each condition by its left side. An unknown is an
essential (Dirichlet) condition. A flux times the normal `n` is a natural
(Neumann) condition. The flux is the expression in `div` in the equation.

A part of the boundary with no condition is free: zero traction, zero heat
flux. `check` lists each free part, so you see what you did not write:

```toon
free_boundary[1]{physics,field,region,condition}:
  solid,u,"boundary(beam) - root - tip","sigma*n = 0"
```

---

## 10. Methods

A `method` block chooses the numerical method of a physics:

```skein
method solid: cg {
    mesh: m
    order: 2
}
```

The physics does not change when the method changes. The kinds are `cg`
(continuous Galerkin), `dg` (discontinuous Galerkin), `fv` (finite volume),
`rkpm` and `gfd` (meshless).

---

## 11. Steps

The analysis is a sequence of steps, as in Abaqus. Each step has a
procedure:

```skein
step preload: static {
    bc root: u = 0
    bc tip: sigma*n = (0, 0, -1) MPa
}

step vibrate: modal(count: 10) {
    output u
}

step cool: transient(end: 600 s, dt: 5 s, scheme: bdf2) {
    bc wall: T = 293 K
    output T every 10
}
```

A condition stays in force in the steps after its step. `release` removes
it. A `bc` outside all steps applies in each step.

---

## 12. Verification by default

Each run does these checks, with no option:

| Before the solve | After the solve |
|---|---|
| units, types, assumptions, boundary conditions (all in `check`) | residual of the solver |
| no free rigid-body motion | equilibrium of the reactions and the loads |
| mesh quality | energy balance |
| physical material values | conservation of mass (flow) |

A failed check stops the run or gives a non-zero exit status. To turn a
check off, write the reason in the model:

```skein
allow rigid_body: "the model is free in space; the step is modal"
```

The reason stays in the file, so a reviewer sees it. The
[check catalogue](checks.md) lists all the checks.

---

## 13. Fluids and stabilization

The lid-driven cavity ([sample 02](../samples/02_cavity.skn)) has two
unknowns, the velocity and the pressure:

```skein
physics flow on cav solves v, p {
    material fluid
    assume incompressible, newtonian, laminar
    let sigma: Tensor2<Pressure> = -p*I + 2*mu*sym(grad(v))
    eq rho*(ddt(v) + grad(v)*v) - div(sigma) = 0
    eq div(v) = 0
}

method flow: cg {
    mesh: m
    order: { v: 1, p: 1 }
    stabilize: supg_pspg
}
```

Linear velocity with linear pressure is unstable in Galerkin form. Without
`stabilize`, `check` gives the error `inf_sup` and suggests `supg_pspg`.
The compiler has the strong form, so it makes the SUPG and PSPG terms from
the exact residual of the momentum equation.

---

## 14. Birth and death

`birth` adds cells to the model in a step, and `death` removes them:

```skein
step lift_2: static {
    birth lift2
}

step excavate: static {
    death berm
}
```

A born region starts with zero strain on the current deformed shape. A dead
region leaves the system. The mesh holds all regions from the start, so the
regions must be parts of the meshed region. [Sample 03](../samples/03_birth_death.skn)
builds an embankment in three lifts.

---

## 15. Submodels

A submodel takes the displacement of a global model on its cut boundary:

```skein
let global: Results = results("01_cantilever.skn", step: "load")
let cut: Surface = boundary(block, abs(x - 40 mm) == 15 mm || abs(z - 20 mm) == 15 mm)

step load: static {
    bc cut: u = global.u
}
```

The checker reads the global model, so `global.u` has a type. It checks that
the submodel lies inside the global body. After the solve, the check
`submodel_cut` compares the force on the cut with the global stress.

---

## 16. Without a mesh

A `points` block makes a point cloud instead of a mesh. The physics does not
change:

```skein
points c on quarter {
    spacing: 2 mm
    refine hole: 0.4 mm
}

method solid: rkpm {
    points: c
    basis: linear
    support: 2.5
    integrate: scni
    bc: nitsche
}
```

[Sample 05](../samples/05_meshless_plate.skn) solves a plate with a hole
this way.

---

## 17. Tables and expectations

A long list of data is a table. The header gives the type and the unit of
each column. The rows are TOON rows:

```skein
table ghia_u{y: Length in m, u: Velocity in m/s}:
    1.0000,  1.00000
    0.5000, -0.20581
    0.0000,  0.00000
```

An `expect` checks a result after the solve. A table column makes one test
for each row:

```skein
expect centreline: v.x at (a/2, ghia_u.y) == ghia_u.u within 0.02 m/s
```

---

## 18. Where to go next

- [The language reference](reference.md): each construct, with an example.
- [The specification](skein_spec.md): the exact rules, the grammar and the
  code generation.
- [The check catalogue](checks.md): every error, lint rule and
  verification check.
- [The samples](../samples/README.md): five complete models.
