# The model `heat_fin`

This document is the model in the file `examples/femahoy/heat_fin.skn`, written by `skein md`. The items are in the order of the file, and each item gives its line.

femahoy: validation/cases/heat-fin/job.inp A steady aluminium pin fin, 100 mm long, 10 x 10 mm section. The base is at 100 degC. The sides and the tip give heat to the air at 20 degC. Reference (Incropera et al. 2007, Table 3.4, case A): the tip temperature 70.8757 degC.

## Model

The model `heat_fin` (line 6) is three-dimensional. Its unit system is `mm_N_MPa_s`: length in mm, mass in t, time in s, force in N, stress in MPa, temperature in K. The unit system sets the units of the output. Each value in the source carries its own unit.

## Parameters and constants

| Name | Kind | Type | Value | Line |
|---|---|---|---|---|
| `Lf` | let | `Length` | $100\ \mathrm{mm}$ | 8 |
| `k_al` | let | `Conductivity` | $200\ \mathrm{W}/(\mathrm{m}\cdot \mathrm{K})$ | 9 |
| `h_air` | let | `HeatTransfer` | $50\ \mathrm{W}/(\mathrm{m}^{2}\cdot \mathrm{K})$ | 10 |
| `T_air` | let | `Temperature` | $20\ {}^{\circ}\mathrm{C}$ | 11 |
| `dT` | let | `Temperature` | $80\ \mathrm{K}$ | 12 |

The reference values of the scales line (line 14) are $\mathit{Lf}$, $\mathit{dT}$, $k_{\mathrm{al}}$. The solver divides each quantity by the product of powers of these values that has its dimension.

## Geometry and regions

- `fin` (line 16) is a Solid: the box from $(0, 0, 0)\ \mathrm{mm}$ to $(\mathit{Lf}, 10\ \mathrm{mm}, 10\ \mathrm{mm})$.
- `base` (line 17) is a Surface: the part of the boundary of `fin` where $x = 0\ \mathrm{mm}$.
- `tip` (line 18) is a Surface: the part of the boundary of `fin` where $x = \mathit{Lf}$.
- `skin` (line 19) is a Surface: the part of the boundary of `fin` where $x > 0\ \mathrm{mm}$.

## Materials

### Material `alu` (line 21)

| Property | Type | Value | Line |
|---|---|---|---|
| `k` | `Conductivity` | $k_{\mathrm{al}}$ | 22 |

## Discretization

### Mesh `m` (line 25)

On `fin`.

- `size`: $5\ \mathrm{mm}$
- `cell`: `tet`

## Physics

### Physics `thermal` (line 30)

On `fin`.

Unknowns: $T$ of type `Temperature`, initial value $20\ {}^{\circ}\mathrm{C}$.

Material: `alu` (line 31).

Assumptions (line 32):

- `isotropic`

Let $q$, of type `Vector<HeatFlux>` (line 33):

$$
q = \operatorname{fourier}(k, T)
$$

Equation (line 34):

$$
\nabla\cdot q = 0
$$

## Boundary conditions

| Region | Condition | Line |
|---|---|---|
| `base` | $T = 100\ {}^{\circ}\mathrm{C}$ | 37 |
| `skin` | $q \cdot n = h_{\mathrm{air}} (T - T_{\mathrm{air}})$ | 38 |

## Steps

### Step `steady` (line 40)

Procedure: `static`.

Output (line 41): $T$.

## Checks

| Name | Quantity | Expected value | Limit | Line |
|---|---|---|---|---|
| `t_tip` | $\operatorname{mean}(T, \mathit{tip})$ | $70.8757\ {}^{\circ}\mathrm{C}$ | $0.1\ \mathrm{K}$ | 44 |
