# Notes on the Skein model

The Skein model of the [problem sheet](../problem.md) has two files:

| File | Contents |
|---|---|
| `fuselage.skn` | the global model (sheet 2 to 6.3, checks C1 to C5) |
| `door_corner.skn` | the door-corner submodel (sheet 6.4, check C6, outputs 5 to 7) |

The element order is 2 for the shells and the beams (the `method` lines).

## Check status

`skein check` gives 15 errors in each file. Each error is on a line under
a `// gap` comment, and line 1 of each file names these gaps. Thus both
files pass the phase 1 gate of the
[femahoy coverage](../../../femahoy_coverage.md): a model passes `skein
check`, or it has a `// blocked-by:` mark. Both files also pass `skein fmt
--check`.

## Gaps

Each file names its open gaps on line 1 (`// blocked-by:`), and puts a
`// gap` comment where the gap applies. G31, G32 and G33 are in the
[femahoy coverage](../../../femahoy_coverage.md#gaps-of-the-787-study),
with the spec changes that close them. N2 and N3 do not make `check` fail,
so they are only on this page.

| Gap | Feature | Why the spec cannot state it | Smallest spec change |
|---|---|---|---|
| G3 | Composite shells: the skin and the pad-up | [femahoy coverage](../../../femahoy_coverage.md#gaps-that-block-a-case). The model also needs a ply axis, the outer face and a load per area. | The G3 change, with these template arguments: `layup` (a table column of ply angles), `mirror`, `axis`, `top` and `q`. The thickness `h` and the ply values come from the material. The template gives `N` (force per length, membrane and transverse shear) and `eps_m` (membrane strain). |
| G31 | Offset beams on the skin: the stringers and the frames, and their axial force `Nx` | [femahoy coverage](../../../femahoy_coverage.md#gaps-of-the-787-study) | The G31 change. |
| G32 | The skin, the pad-up, the stringers, the frames and the door edge cut by a `Solid` | [femahoy coverage](../../../femahoy_coverage.md#gaps-of-the-787-study) | The G32 change. |
| G33 | The shell resultants `N` and `eps_m`, the edge loads `N*n`, and `global.e_hoop` in the submodel | [femahoy coverage](../../../femahoy_coverage.md#gaps-of-the-787-study) | The G33 change. |
| N2 | Outputs 1, 3 and 5: the point of a maximum, and the maximum of each beam | `max(f, r)` gives one value on one region (§14.5). | A report function `argmax(f, r)` that gives the point. `max(f, r)` on a `pattern` gives one value for each copy. |
| N3 | Quadrilateral shells (sheet 6.3) | `cell: quad` needs a `box` or a `rect` (§9.1). | Let `cell: quad` apply to each `Surface`, with triangles where the shape needs them. The model uses `cell: tri`, the nearest valid form, so N3 is not in `blocked-by`. |

## Nearest forms and choices

- The door corner radius is in the x-s plane (sheet 2.3). Skein has no
  shape in the x-s plane of a cylinder. The model fillets the radial edges
  of a box that turns to the door angle (`door_box`). The side faces cut
  the skin at the correct angles. On the curved skin, the corner is near a
  circle of radius 0.178 m, but it is not exact. The pad-up band uses the
  same function with 0.150 m more.
- The spec does not give the sign of `rotate`, the axis point of
  `pattern_polar`, or the meaning of the vector of `pattern`. The model
  uses the right-hand rule, an axis through the origin, and the step
  between two copies. A comment states each choice.
- The beam ends get the symmetry conditions and the submodel cut values
  through the tie to the skin, so they have no `bc` of their own.
- The skin is one physics with two materials (`skin_lam`, `pad_lam`).
  Each material gives the ply values and the laminate thickness `h`.

## Parser and formatter findings

- A `use` library file (§4.1) has no header, and the parser rejects it
  with `header_order`. Thus each file holds its own copy of the dimensions,
  the materials, the geometry and the physics: about 165 shared lines.
  With `use`, one library file can hold them.
- The parser accepts a list `[...]` only in a procedure argument, not in a
  template argument. The layup is a table for this reason.
- `skein fmt` breaks a physics header at a `|` in its region, and it
  puts the second unknown of `solves` on a new line. The model names the
  regions `skin` and `beams` to keep the `|` out of the header.
- `skein fmt` writes `table quasi{angle: ...}` with no space before `{`.

## Size

Counted with tiktoken `o200k_base`:

| File | Lines | Tokens |
|---|---|---|
| `fuselage.skn` | 230 | 2914 |
| `door_corner.skn` | 218 | 2766 |
| Total | 448 | 5680 |
