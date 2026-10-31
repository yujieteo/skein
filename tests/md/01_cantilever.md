# The model `cantilever`

This document is the model in the file `samples/01_cantilever.skn`, written by `skein md`. The items are in the order of the file, and each item gives its line.

A steel cantilever with a lightening hole near the root, under an end load. Linear static, 3D, quadratic tetrahedra (the default method). Sample 04 is a submodel of the region around the hole. The study varies the length and the load: \`skein run 01\_cantilever.skn --study sizes\`, or one case with \`--set L=150mm --set P=3kN\`.

## Model

The model `cantilever` (line 7) is three-dimensional. Its unit system is `mm_N_MPa_s`: length in mm, mass in t, time in s, force in N, stress in MPa, temperature in K. The unit system sets the units of the output. Each value in the source carries its own unit.

## Parameters and constants

| Name | Kind | Type | Value | Line |
|---|---|---|---|---|
| `L` | param | `Length` | $200\ \mathrm{mm}$, in $100\ \mathrm{mm} \ldots 400\ \mathrm{mm}$ | 9 |
| `W` | let | `Length` | $20\ \mathrm{mm}$ | 10 |
| `H` | let | `Length` | $40\ \mathrm{mm}$ | 11 |
| `r` | let | `Length` | $6\ \mathrm{mm}$ | 12 |
| `P` | param | `Force` | $2\ \mathrm{kN}$, in $0.5\ \mathrm{kN} \ldots 10\ \mathrm{kN}$ | 13 |

The reference values of the scales line (line 15) are $H$, $\frac{P}{W H}$. The solver divides each quantity by the product of powers of these values that has its dimension.

## Geometry and regions

- `bore` (line 17) is a Solid: the cylinder with its base centre at $(40, 0, 20)\ \mathrm{mm}$, axis $(0, 1, 0)$, radius $r$ and length $W$.
- `beam` (line 18) is a Solid: the box from $(0, 0, 0)\ \mathrm{mm}$ to $(L, W, H)$ minus `bore`.
- `root` (line 19) is a Surface: the part of the boundary of `beam` where $x = 0\ \mathrm{mm}$.
- `tip` (line 20) is a Surface: the part of the boundary of `beam` where $x = L$.
- `hole` (line 21) is a Surface: the part of the boundary of `beam` that touches `bore`.

## Materials

### Material `steel` (line 23)

| Property | Type | Value | Line |
|---|---|---|---|
| `E` | `Pressure` | $210\ \mathrm{GPa}$ | 24 |
| `nu` | `Real` | $0.3$ | 25 |

## Discretization

### Mesh `m` (line 28)

On `beam`.

- `size`: $5\ \mathrm{mm}$
- `cell`: `tet`
- refine `hole`: $1\ \mathrm{mm}$

## Physics

### Physics `solid` (line 34)

On `beam`.

Unknowns: $u$ of type `Vector<Length>`.

Material: `steel` (line 35).

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

### Step `load` (line 38)

Procedure: `static`.

| Region | Condition | Line |
|---|---|---|
| `root` | $u = 0$ | 39 |
| `tip` | $\sigma n = (0, 0, \frac{-P}{W H})$ | 40 |

Output (line 41): $u$, $\sigma$, $\operatorname{von}_{\mathrm{mises}}(\sigma)$.

Report (line 42): $\operatorname{reaction}(\mathit{root})$.

## Other items

### Study `sizes` (line 45)

Procedure: `grid`.

- `L`: $[150, 175, 200, 225]\ \mathrm{mm}$
- `P`: $[1, 2, 3, 4]\ \mathrm{kN}$
