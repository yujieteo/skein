# Notes on the Skein model

The Skein model of the [problem sheet](../problem.md) has two files:

| File | Contents |
|---|---|
| `fuselage.skn` | the global model (sheet 2 to 6.3, checks C1 to C5) |
| `door_corner.skn` | the door-corner submodel (sheet 6.4, check C6, outputs 5 to 7) |

The element order is 2 for the shells and the beams (the `method` lines).

## Check status

`skein check` is not built yet. It prints "check is not implemented yet"
for each file. The model uses `skein fmt --check` in its place. That
command parses each file and compares it with the standard layout. Both
files pass it.

## Gaps

Each file names its open gaps on line 2 (`// blocked-by:`), and puts a
`// gap` comment where the gap applies.

| Gap | Feature | Why the spec cannot state it | Smallest spec change |
|---|---|---|---|
| G3 | Composite shells: the skin and the pad-up | [femahoy coverage](../../../femahoy_coverage.md#gaps-that-block-a-case). The model also needs a ply axis, the outer face and a load per area. | The G3 change, with these template arguments: `layup` (a table column of ply angles), `mirror`, `axis`, `top` and `q`. The thickness `h` and the ply values come from the material. The template gives `N` (force per length, membrane and transverse shear) and `eps_m` (membrane strain). |
| N1 | Offset beams on the skin: the stringers and the frames | §4.2 gives `timoshenko(A, I_y, I_z)` with no source, no torsion constant, no section axis and no offset. No rule ties a beam to the shell points on its line. | `timoshenko` arguments `J`, `y_axis` and `offset`, and a let `Nx` (axial force). In §12.1, let the right side of an essential `bc` be an unknown of another physics on the same points (a tie). |
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

## Problem sheet issue

Check C6 compares the strain at the mid-point of each edge of the
submodel square. The mid-points of the forward edge (x = 0.0554 m) and
the lower edge (s = 3.3376 m) are in the door cutout, where there is no
skin. The submodel checks the other two mid-points (`c6_aft`, `c6_top`).

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
| `fuselage.skn` | 217 | 2662 |
| `door_corner.skn` | 202 | 2474 |
| Total | 419 | 5136 |
