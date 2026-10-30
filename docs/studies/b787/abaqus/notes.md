# Abaqus deck: notes

This folder holds the Abaqus input for the problem in
[`../problem.md`](../problem.md). The writer used only that sheet and
general Abaqus practice.

## Files

| File | Kind | Content |
|---|---|---|
| `global.inp` | hand-written | global quarter model, step, loads, checks |
| `door_corner.inp` | hand-written | door-corner submodel |
| `properties.inp` | hand-written | material, orientation, shell and beam sections; both decks include it |
| `mesh.py` | script | writes the four generated files below |
| `global_mesh.inp`, `door_corner_mesh.inp` | generated | nodes, elements, sets |
| `global_loads.inp`, `door_corner_loads.inp` | generated | nodal end closure load and door load |

The generated files are not in git. To make them and run:

```sh
python3 mesh.py
abaqus job=global
abaqus job=door_corner globalmodel=global
```

`mesh.py` gives the same output on each run. SHA-256 of the output:

```text
68eeecb81b3d5f49b26e095fbe294b3094e93b83836e4ce5ca24bcc5d92101c1  global_mesh.inp
e2a6ee110fff8ad6aea8537050e37dc8aabde82bb1f29c208696535e1565585d  global_loads.inp
dac853cc6bd8169bd6d09918de979e547f12a3b37f09296d8eab167a32b4d92c  door_corner_mesh.inp
f3f8c2aa5e86d66a8092e30ccaea9ec84d5ca203685048ac4396fa0ac37e1f35  door_corner_loads.inp
```

## Size

Lines count non-blank lines. Tokens use tiktoken `o200k_base`.

| Part | Lines | Lines without `**` comments | Tokens |
|---|---|---|---|
| Hand-written: `global.inp`, `door_corner.inp`, `properties.inp` | 237 | 164 | 3,580 |
| Mesh script: `mesh.py` | 344 | 344 | 5,572 |
| Generated: the four mesh and load files | 55,564 | 55,553 | 1,063,526 |

The 72 ply lines of the skin and pad-up sections are 72 of the 164
keyword lines.

## Decisions that the sheet leaves open

- All elements are first order: S4R, S3R and B31.
- The pad-up is the set of shells with a centroid within 0.150 m of the
  cutout edge. The band edge above the door falls on the stringer at 54°,
  so the band there is 0.142 m wide.
- G13 = G12 and G23 = 3.0 GPa. The sheet gives no value.
- The beam shear centre is at the centroid. The sheet gives no value.
- The end closure load and the door load are lumped nodal forces
  (`*CLOAD`). Thus the load directions are explicit.
- Abaqus has no assertion keyword. Each check is a print request
  (`*EL PRINT`, `*NODE PRINT`) with the expected value in a comment.

## Checks done

- **femahoy:** `fh run global.inp` stops with
  `*ERROR reading *BOUNDARY. Card image: SYM_X,XSYMM`. femahoy runs
  CalculiX, and CalculiX has no `XSYMM` label. The deck also uses
  `*BEAM GENERAL SECTION, SECTION=GENERAL` with `*CENTROID` and
  `*SUBMODEL, GLOBAL ELSET`, which are Abaqus only. Thus femahoy cannot
  solve this deck.
- **Normals:** all 21,393 global and 2,305 submodel shells have an outward
  normal. A script computed this from the generated nodes.
- **Mesh:** no inverted element. The smallest corner angle is 15.8° (global)
  and 17.0° (submodel). The S3R elements are only at the corner radii
  (4 global, 21 submodel).
- **Load resultants**, from the generated mesh (pressure on the flat
  element faces plus the nodal loads). `mesh.py` prints them.

| Model | Fx | Fy | Fz |
|---|---|---|---|
| Global | 1272.96 kN | 3408.10 kN | −5.09 kN |
| Submodel | 0 | 25.30 kN | 14.04 kN |

## Problems in the sheet

1. **C5 and C4 cannot pass.** The sheet puts the door load only on the
   straight part of the aft edge, so the door strips level with the two
   corner radii (0.178 m high) send no load to the skin. Thus the applied Fy is
   3408.1 kN, which is 0.55 % less than the 3427.1 kN in C5. The
   tolerance is 0.1 %. The applied Fz is −5.09 kN, which is 0.40 % of
   1273.0 kN. The C4 limit is 0.1 %. This deck applies the sheet loads as
   written, so the reactions will equal the table above, not the C4 and
   C5 values.
2. **C6 has two points only.** The mid-points of the forward edge
   (x = 0.0554) and the upper edge (s = 3.3376) of the submodel region are
   in the cutout. Only the lower edge and the aft edge have skin at the
   mid-point.

## Not checked without Abaqus

- A positive pressure on the `SNEG` face acts along the positive normal,
  so outward here.
- The `*CENTROID` offset uses I11 and I22 about the centroid.
- `*ORIENTATION` with `3, 0.` gives shell axis 1 along global X.
