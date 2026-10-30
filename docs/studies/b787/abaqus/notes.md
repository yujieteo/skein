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
76c0ef25277f7d8d76b031fba234e2b7a3007fb9ea2f2b033be9f7d5d9fb0729  global_mesh.inp
3345e72d7c9be74570f9340f911e03e3641e16e09e4fe51a8d352216bb01bec4  global_loads.inp
539656d15c1efc8c4821108b52abd71aad29e537c3df04b5e0be88ab2fe8e283  door_corner_mesh.inp
8babba889635509a3a575871c355a1920e70b9d2bbf1c5ce8be3e2282f9b24e9  door_corner_loads.inp
```

## Size

Lines count non-blank lines. Tokens use tiktoken `o200k_base`.

| Part | Lines | Lines without `**` comments | Tokens |
|---|---|---|---|
| Hand-written: `global.inp`, `door_corner.inp`, `properties.inp` | 240 | 164 | 3,639 |
| Mesh script: `mesh.py` | 352 | 352 | 5,690 |
| Generated: the four mesh and load files | 55,847 | 55,836 | 1,067,377 |

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
| Global | 1272.96 kN | 3427.08 kN | −0.6 N |
| Global, door load only | 0 | 89.05 kN | 23.86 kN |
| Submodel | 0 | 32.72 kN | 18.55 kN |

The global values agree with C3 and C5 within 0.01 %, and with C4 within
0.001 %. The door load agrees with the 89.12 kN and 23.88 kN of sheet
6.2 within 0.1 %. Thus the reactions must agree with C3 to C5 when the
solver balances the loads.

## Sheet corrections from Firstmate

The deck uses two corrections that Firstmate sent after the first sheet:

1. The door load is p_U × x_e(s) per unit s, on the straight aft edge and
   on the two aft corner arcs. The sheet copy in this branch has it.
2. C6 uses only the mid-points of the region edges that have skin: the
   top edge (s = 2.7373) and the aft edge (x = 0.6554). The mid-points of
   the forward and the lower edge are in the cutout. The deck names the
   sets `C6_TOP` and `C6_AFT`. The sheet copy in this branch has it.

## Not checked without Abaqus

- A positive pressure on the `SNEG` face acts along the positive normal,
  so outward here.
- The `*CENTROID` offset uses I11 and I22 about the centroid.
- `*ORIENTATION` with `3, 0.` gives shell axis 1 along global X.
