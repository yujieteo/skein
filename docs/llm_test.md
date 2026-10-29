# The LLM test

Skein must be easy for an LLM to write. This test measures it. It is the
gate of phase 1 in the [roadmap](roadmap.md): if Skein is not clearly better
than an Abaqus or APDL deck on both measures, the syntax changes.

**Status:** the protocol and the reference decks are ready. The test has not
run, because it needs `skein check`, which phase 1 builds. The results
table is empty.

---

## 1. Measures

For each of the five samples, and for each language (Skein, and an Abaqus
or APDL deck):

| Measure | Definition |
|---|---|
| first-try pass rate | the fraction of attempts whose first answer passes the checker with no error |
| tokens | the output tokens of the answer, counted with the tokenizer of the model under test. Only the model text counts, not the explanation around it. |

The checker for Skein is `skein check`. The checker for a deck is the
solver's own syntax check when a licence is available (`abaqus datacheck`,
or APDL in batch mode up to `SOLVE`). Without a licence, a person reviews
the deck against the reference deck in `llm_test/decks/` with this list:
the geometry, the mesh controls, the material, the boundary conditions,
the loads, the analysis sequence and the outputs. A deck passes when all
seven items are right.

**Clearly better** means both of these:

- the Skein first-try pass rate is at least 20 percentage points higher than
  the deck pass rate, summed over the five samples;
- the median Skein token count is at most 0.8 times the median deck token
  count.

Reason for the limits: smaller differences are inside the noise of 5
attempts for each sample.

---

## 2. Protocol

1. Use a fresh model session for each attempt, with no memory of earlier
   attempts.
2. For Skein, give the model the [tutorial](tutorial.md) and the
   [reference](reference.md) in the context. Do not give it the samples.
   For a deck, give no documentation, because the model knows APDL and
   Abaqus from its training.
3. Give the one-line request of §3, then "Write the complete Skein model."
   or "Write the complete APDL input deck." (Abaqus for the cavity).
4. Take the first answer. Do not correct it and do not ask again.
5. Run the checker on the answer, and count the tokens.
6. Do 5 attempts for each sample and each language: 50 attempts in total.
7. Record the model id and the date.

The context tokens of the Skein documentation are not in the token measure,
because they are the same for each answer. Record them in the results, so
that a reader can compare the total cost.

---

## 3. The requests

| # | Request |
|---|---|
| 01 | A 200 × 20 × 40 mm steel cantilever, fixed at x = 0, with a 6 mm radius hole through its 20 mm width at x = 40 mm, z = 20 mm, and a 2 kN downward load on the free end. Use quadratic tetrahedra, 5 mm cells, 1 mm at the hole. Report the root reaction and output the von Mises stress. |
| 02 | A 2D lid-driven cavity, 1 m square, lid speed 1 m/s, Re = 100, steady incompressible laminar flow. Use a 64 × 64 resolution and fix the pressure at one corner. Output velocity and pressure. |
| 03 | A 2D plane strain embankment: a 60 × 10 m clay foundation (E 40 MPa, ν 0.3, ρ 1900 kg/m³) and three 2 m fill lifts (E 25 MPa, ν 0.3, ρ 1800 kg/m³), base 24 m and top 12 m, plus a triangular toe berm. Apply gravity, add the lifts one per step, then remove the berm. Fix the base and the sides in x. |
| 04 | A submodel of the hole in model 01: a block x 25 to 55 mm, z 5 to 35 mm, 1.5 mm cells and 0.25 mm at the hole, with the cut-boundary displacements from model 01. Output the von Mises stress. |
| 05 | A quarter of a 200 × 200 × 2 mm aluminium plate (E 70 GPa, ν 0.33) with a central 10 mm radius hole, under 100 MPa tension in x, plane stress, with symmetry. Check that the peak stress at the hole is near 300 MPa. For Skein, use a meshless method (RKPM). |

Request 05 asks for a meshless method only in Skein, because APDL and
Abaqus have no RKPM for this problem. The deck uses the finite element
method for the same physical model.

---

## 4. Reference answers

| # | Skein | Deck |
|---|---|---|
| 01 | [samples/01_cantilever.skn](../samples/01_cantilever.skn) | [APDL](../llm_test/decks/01_cantilever.mac) |
| 02 | [samples/02_cavity.skn](../samples/02_cavity.skn) | [Abaqus/CFD](../llm_test/decks/02_cavity.inp) |
| 03 | [samples/03_birth_death.skn](../samples/03_birth_death.skn) | [APDL](../llm_test/decks/03_birth_death.mac) |
| 04 | [samples/04_submodel.skn](../samples/04_submodel.skn) | [APDL](../llm_test/decks/04_submodel.mac) |
| 05 | [samples/05_meshless_plate.skn](../samples/05_meshless_plate.skn) | [APDL](../llm_test/decks/05_plate.mac) |

The decks are written for this test. They copy no text from the Abaqus or
APDL manuals. They are not run, because no licence was available.

---

## 5. An early signal: the size of the reference answers

This is not the test, because a person wrote both answers. But it shows a
risk early. The byte counts of the reference answers, without comments and
blank lines:

| # | Skein (bytes) | Deck (bytes) | Skein / deck |
|---|---|---|---|
| 01 | 920 | 405 | 2.3 |
| 02 | 1236 | 914 | 1.4 |
| 03 | 1543 | 1049 | 1.5 |
| 04 | 921 | 505 | 1.8 |
| 05 | 1067 | 374 | 2.9 |

The Skein answers are 1.4 to 2.9 times longer in bytes. Bytes are not
tokens, and the APDL commands with many commas can use more tokens for
each byte. But the token measure will probably fail if the syntax stays as
it is. Some of the extra length is content that the decks do not have: the
equations, the units, and the checks (`expect`, `report`).

## 6. Syntax changes to try if the test fails

Try these in sequence, and run the test again after each one:

1. **Declare the field in `solves`.** `physics solid on beam solves u:
   Vector<Length> {` removes the `field` line.
2. **Physics templates in a standard library.** `physics solid: elastic(steel)
   on beam solves u` for the common physics. The strong form stays
   available for new physics.
3. **Defaults for `method`.** `cg` with order 2 for solids and Taylor-Hood
   for flow when the model has no `method` block.
4. **A shorter header.** `model cantilever 3d mm_N_MPa_s` on one line.

Each change removes text that an LLM must write, but each change also adds
a rule to learn. The test decides.

---

## 7. Results

| Date | Model | # | Language | Pass rate | Median tokens |
|---|---|---|---|---|---|
| | | | | | |
