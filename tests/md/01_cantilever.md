# The model `cantilever`

This document is the model in the file `samples/01_cantilever.skn`, written by `skein md`. The items are in the order of the file, and each item gives its line.

A steel cantilever with a lightening hole near the root, under an end load. Linear static, 3D, quadratic tetrahedra (the default method). Sample 04 is a submodel of the region around the hole.

## Model

The model `cantilever` (line 5) is three-dimensional. Its unit system is `mm_N_MPa_s`: length in mm, mass in t, time in s, force in N, stress in MPa, temperature in K. The unit system sets the units of the output. Each value in the source carries its own unit.

## Parameters and constants

| Name | Kind | Type | Value | Line |
|---|---|---|---|---|
| `L` | let | `Length` | $200\ \mathrm{mm}$ | 7 |
| `W` | let | `Length` | $20\ \mathrm{mm}$ | 8 |
| `H` | let | `Length` | $40\ \mathrm{mm}$ | 9 |
| `r` | let | `Length` | $6\ \mathrm{mm}$ | 10 |
| `P` | let | `Force` | $2\ \mathrm{kN}$ | 11 |

The reference values of the scales line (line 13) are $H$, $\frac{P}{W H}$. The solver divides each quantity by the product of powers of these values that has its dimension.

## Geometry and regions

- `bore` (line 15) is a Solid: the cylinder with its base centre at $(40, 0, 20)\ \mathrm{mm}$, axis $(0, 1, 0)$, radius $r$ and length $W$.
- `beam` (line 16) is a Solid: the box from $(0, 0, 0)\ \mathrm{mm}$ to $(L, W, H)$ minus `bore`.
- `root` (line 17) is a Surface: the part of the boundary of `beam` where $x = 0\ \mathrm{mm}$.
- `tip` (line 18) is a Surface: the part of the boundary of `beam` where $x = L$.
- `hole` (line 19) is a Surface: the part of the boundary of `beam` that touches `bore`.

## Materials

### Material `steel` (line 21)

| Property | Type | Value | Line |
|---|---|---|---|
| `E` | `Pressure` | $210\ \mathrm{GPa}$ | 22 |
| `nu` | `Real` | $0.3$ | 23 |

## Discretization

### Mesh `m` (line 26)

On `beam`.

- `size`: $5\ \mathrm{mm}$
- `cell`: `tet`
- refine `hole`: $1\ \mathrm{mm}$

## Physics

### Physics `solid` (line 32)

On `beam`.

Unknowns: $u$ of type `Vector<Length>`.

Material: `steel` (line 33).

The physics uses the template `linear_elastic` of the standard library, with $b = 0$. Its statements are:

Assumptions:

- `small_strain`
- `isotropic`

Let $\sigma$, of type `Tensor2<Pressure>`:

$$
\sigma = \operatorname{hooke}(E, \nu, \operatorname{strain}(u))
$$

Equation:

$$
\nabla\cdot \sigma + b = 0
$$

## Steps

### Step `load` (line 36)

Procedure: `static`.

| Region | Condition | Line |
|---|---|---|
| `root` | $u = 0$ | 37 |
| `tip` | $\sigma n = (0, 0, \frac{-P}{W H})$ | 38 |

Output (line 39): $u$, $\sigma$, $\operatorname{von}_{\mathrm{mises}}(\sigma)$.

Report (line 40): $\operatorname{reaction}(\mathit{root})$.
