# Notes on the Skein model

The Skein model of the [problem sheet](../problem.md) has three files:

| File | Contents |
|---|---|
| `barrel.skn` | the library file (§4.1) with the dimensions, the geometry and the materials that both models use |
| `fuselage.skn` | the global model (sheet 2 to 6.3, checks C1 to C5) |
| `door_corner.skn` | the door-corner submodel (sheet 6.4, check C6, outputs 5 to 7) |

The element order is 2 for the shells and the beams (the `method` lines).

## Check status

Both models pass `skein check` and `skein fmt --check`. `skein fmt` does
not read the library file `barrel.skn`, because it has no header.

## Run status

Both models run with `skein run` (`src/shells/gen.goose`): MITC9 shells
and 2-node offset beams of order 2, on a grid of the side of the cylinder.
The nodes near a cut move onto it, so the door cutout and the pad-up band
are cuts of the grid, with no triangles. The run on the desktop gives:

| Model | Cells and beams | Dofs | Solve | Peak memory | Checks |
|---|---|---|---|---|---|
| `fuselage.skn` | 45,079 | 714,834 | 34 s, nested dissection | 1.9 GB | C1 to C5 pass |
| `door_corner.skn` | 2,203 | 267,126 | 8 s | 0.7 GB | C6 passes |

| Value | Result |
|---|---|
| C1: hoop stress at P | 54.3 MPa, ratio 0.886 to p_U R/t |
| C2: axial stress at P | 29.4 MPa, ratio 0.960 to p_U R/(2t) |
| C3, C4, C5: reactions | -1272.96 kN, 0.03 kN and -3427.61 kN |
| C6: hoop strain on the cut, aft and top | 1621 and 289 µε, 1.5 % and 1.4 % from the global model |
| Outputs 5 to 7 | e_1 max 10,524 µε on the cutout edge, K = 13.1, MS = -0.62 |

The largest strain is on the upper aft corner arc, near the end of the
straight aft edge. The plug door load is a radial line load on the free
edge of the skin, so the skin bends there. 6 cm aft of the edge, the
hoop strain of the submodel is 2470 µε, and of the global model 2700 µε.

The submodel takes `global.w`, `global.theta` and `global.e_hoop` from the
state that the run of `fuselage.skn` saves. A run of `door_corner.skn`
runs `fuselage.skn` first if the state is not there.

## Gaps

No gap blocks the models. These gaps are closed, and `skein run` runs
them:

| Gap | Feature | Closed by |
|---|---|---|
| G3 | Composite shells: the skin and the pad-up, with a ply axis, the outer face and a load per area | spec §4.3: `reissner_mindlin(layup, mirror, axis, top, q)` |
| G31 (N1) | Offset beams on the skin and their axial force `Nx` | spec §4.3: `timoshenko(A, I_y, I_z, J, y_axis, offset, q)` and `let Nx`; spec §12.1: the tie |
| G32 | The skin, the pad-up, the stringers, the frames and the door edge cut by a `Solid` | spec §8.4: `a & b` has the lower dimension of `a` and `b` |
| G33 | The shell resultants `N` and `eps_m`, the edge loads `N*n`, and `global.e_hoop` in the submodel | spec §4.3: `let N` and `let eps_m`; spec §15.1: the `let` fields of a `Results` value |
| N2 | Outputs 1, 3 and 5: the point of a maximum, and the maximum of each beam | spec §14.5: `argmax(f, r)` gives the point, and `max(f, r)` on a `pattern` gives one value for each copy |
| N3 | Quadrilateral shells (sheet 6.3) | spec §9.1: `cell: quad` on each `Surface`, with triangles where the shape needs them. The models use `cell: quad`. |

## Nearest forms and choices

- The door corner radius is in the x-s plane (sheet 2.3). Skein has no
  shape in the x-s plane of a cylinder. The model fillets the radial edges
  of a box that turns to the door angle (`door_box`). The side faces cut
  the skin at the correct angles. On the curved skin, the corner is near a
  circle of radius 0.178 m, but it is not exact. The pad-up band uses the
  same function with 0.150 m more.
- The spec does not give the sign of `rotate`. The model uses the
  right-hand rule. Spec §8.1 gives the axis of `pattern_polar` through the
  origin, and the vector of `pattern` as the step between two copies.
- The beam ends get the symmetry conditions and the submodel cut values
  through the tie to the skin, so they have no `bc` of their own.
- The skin is one physics with two materials (`skin_lam`, `pad_lam`).
  Each material gives the ply values and the laminate thickness `h`.

## Parser and formatter findings

- A library file holds only `material`, `param`, `let` and `fn` items.
  Each model therefore repeats the ply table, mesh, physics and conditions.
- The parser accepts a list `[...]` only in a procedure argument, not in a
  template argument. The layup is a table for this reason.
- `skein fmt` breaks a physics header at a `|` in its region, and it
  puts the second unknown of `solves` on a new line. The model names the
  regions `skin` and `beams` to keep the `|` out of the header.
- `skein fmt` writes `table quasi{angle: ...}` with no space before `{`.

## Size

See the [study comparison](../readability.md#2-size) for the size of the
models and the count method.
