# The LLM test

Skein must be easy for an LLM to write. This test measures it. It is the
gate of phase 1 in the [roadmap](roadmap.md): if Skein is not clearly better
than an Abaqus or APDL deck on both measures, the syntax changes. The test
also measures the size of the diff of small edits, because a person must
review each change that an LLM makes.

**Status:** this protocol has no recorded results in §7. `skein check`
is available. Separate experiments have results in the
[cantilever comparison](readability.md#41-results) and the
[B787 study](studies/b787/readability.md#3-llm-runs). Those experiments do
not complete this protocol.

---

## 1. Measures

For each of samples 01 to 05, and for each language (Skein, and an Abaqus
or APDL deck):

| Measure | Definition |
|---|---|
| pass rate | the fraction of attempts that pass the checker with no error within 5 rounds |
| total cost | the input and output tokens of all rounds of an attempt, until it passes or until round 5, counted with the tokenizer of the model under test. For Skein, the language card counts in each round. |
| repair rounds | the number of rounds after the first, until the answer passes |
| first-try pass rate | recorded, not part of the gate |

The checker for Skein is `skein check`. The checker for a deck is the
solver's own syntax check when a licence is available (`abaqus datacheck`,
or APDL in batch mode up to `SOLVE`). Without a licence, a person reviews
the deck against the reference deck in `llm_test/decks/` with this list:
the geometry, the mesh controls, the material, the boundary conditions,
the loads, the analysis sequence and the outputs. A deck passes when all
seven items are right.

**Clearly better** means both of these:

- the Skein pass rate is at least 20 percentage points higher than the deck
  pass rate, summed over samples 01 to 05;
- the median Skein total cost is at most 0.8 times the median deck total
  cost.

Reason for total cost and not the size of the first answer: a Skein answer
holds the equations, the units and the checks, so it is larger. But an
error in a deck often shows only in the solver, after many rounds. The cost
that a user pays is the cost of all the rounds. Reason for the limits:
smaller differences are inside the noise of 5 attempts for each sample.

---

## 2. Protocol

1. Use a fresh model session for each attempt, with no memory of earlier
   attempts.
2. For Skein, give the model the language card of `skein help --llm` in the
   context (about 2000 tokens: the grammar, the keywords, the standard
   library and one example). Do not give it the samples. For a deck, give
   no documentation, because the model knows APDL and Abaqus from its
   training.
3. Give the one-line request of §3, then "Write the complete Skein model."
   or "Write the complete APDL input deck." (Abaqus for the cavity).
4. Run the checker on the answer. If it fails, give the model the checker
   output and "Fix the model." This is one repair round. Stop when the
   answer passes or after round 5.
5. Count the tokens of all rounds.
6. Do 5 attempts for each sample and each language: 50 attempts in total.
7. Record the model id and the date.

For a deck without a licence, the feedback of a repair round is the list of
the items that the person found wrong.

### 2.1 Edit tasks

A person must read each change that an LLM makes. So the test also gives
three edit tasks on the reference answer of sample 01, in each language:

| Task | Request |
|---|---|
| E1 | change the end load from 2 kN to 3 kN |
| E2 | add a second static step that adds a 500 N load in y on the tip |
| E3 | make the hole radius a sweep parameter with the values 4, 6 and 8 mm |

Record the number of changed lines (added plus removed, from `git diff
--numstat`) for each answer that passes. The [layout rules](skein_spec.md#32-layout-for-small-diffs)
aim at 1 line for E1, at most 6 for E2 and at most 3 for E3. These numbers
are a record, not a gate, until a first run shows what a deck needs.

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
blank lines, before and after the syntax changes of §6:

| # | Skein before (bytes) | Skein after (bytes) | Deck (bytes) | Skein after / deck |
|---|---|---|---|---|
| 01 | 920 | 766 | 405 | 1.9 |
| 02 | 1236 | 1055 | 914 | 1.2 |
| 03 | 1543 | 1416 | 1049 | 1.3 |
| 04 | 921 | 766 | 505 | 1.5 |
| 05 | 1067 | 960 | 374 | 2.6 |

The changes cut the Skein answers by 8 to 17 percent. They are still 1.2 to
2.6 times longer than the decks in bytes. The required `scales` line adds
one line to each model. Thus the first answer of Skein stays larger, and
the gate compares the total cost (§1).
Bytes are not tokens, and the APDL commands with many commas can use more
tokens for each byte. Some of the extra length is content that the decks do
not have: the equations, the units, the scales and the checks (`expect`,
`report`).

## 6. Syntax changes made before the test

The captain chose to make these changes before phase 1 (grilling G6.1):

1. **The unknowns in `solves`.** `physics solid on beam solves u:
   Vector<Length> {` replaces the `field` line (spec §10.1).
2. **Physics templates in a standard library.** `physics solid on beam:
   linear_elastic solves u: Vector<Length>` for the common physics (spec §4.2).
   The strong form stays available for new physics.
3. **Defaults for `method`.** `cg` with order 2, Taylor-Hood for
   incompressible flow, and SUPG when the equations have convection (spec
   §13).
4. **A one-line header.** `model cantilever: 3d, mm_N_MPa_s` (spec §3).
5. **Functions of one expression** (spec §7.4). The samples do not need
   them, so they do not change the numbers of §5.

Each change removes text that an LLM must write, but each change also adds
a rule to learn. The test decides whether the syntax changes again.

---

## 7. Results

| Date | Model | # | Language | Pass rate | Median total cost | Median repair rounds | Diff lines E1, E2, E3 |
|---|---|---|---|---|---|---|---|
| | | | | | | | |
