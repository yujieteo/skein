# The LLM test

Skein must be easy for an LLM to write. This test measures it. It is the
gate of phase 1 in the [roadmap](roadmap.md): if Skein is not clearly better
than a CalculiX deck on both measures, the syntax changes. The test also
measures the size of the diff of small edits, because a person must review
each change that an LLM makes.

The baseline is a CalculiX deck that femahoy runs with `fh run`. Reason: no
Abaqus or APDL licence is available, so a check of those decks needs a
person. The driver `scripts/llm/workflow.py` runs both paths with no person.

**Status:** this protocol has no recorded results in §5. Separate
experiments have results in the
[cantilever comparison](readability.md#41-results) and the
[B787 study](studies/b787/readability.md#3-llm-runs). Those experiments do
not complete this protocol.

---

## 1. Measures

For each request of §3, and for each path of `scripts/llm/workflow.py`
(`skein` and `inp`):

| Measure | Definition |
|---|---|
| pass rate | the fraction of attempts that pass within 5 rounds (`passed` in the results file) |
| total cost | the input and output tokens of all rounds of an attempt, until it passes or until round 5, as the driver counts them. For Skein, the language card counts in each round. |
| repair rounds | the number of rounds after the first, until the answer passes |
| first-try pass rate | recorded, not part of the gate |

An attempt passes when the run passes. For `skein`, the run is `skein
check`, then `skein run`. For `inp`, the run is `fh run`, which gives the
deck to ccx.

**Clearly better** means both of these:

- the Skein pass rate is at least 20 percentage points higher than the deck
  pass rate, summed over the requests of §3;
- the median Skein total cost is at most 0.8 times the median deck total
  cost.

Reason for total cost and not the size of the first answer: a Skein answer
holds the equations, the units and the checks, and the card adds about
2700 input tokens. But an error in a deck often shows only in the solver,
after many rounds. The cost that a user pays is the cost of all the rounds.
Reason for the limits: smaller differences are inside the noise of 5
attempts for each request.

---

## 2. Protocol

For each request of §3 and each path:

```sh
./build.sh
scripts/llm/workflow.py skein --request cantilever --runs 5 --out /tmp/llm-test
scripts/llm/workflow.py inp --request cantilever --runs 5 --out /tmp/llm-test
```

The driver gives a fresh model session the request (for Skein, after the
card of `skein help --llm`). If the run fails, it gives the model the output
and "Fix the model." or "Fix the deck.", until the run passes or round 5.
[readability.md §4](readability.md#4-the-llm-workflow) gives the details.
Record the model id and the date from the results file.

### 2.1 Edit tasks

A person must read each change that an LLM makes. So the test also gives
three edit tasks on a passed answer of the cantilever request, in each
path. Give the answer and the task to a fresh session:

| Task | Request |
|---|---|
| E1 | change the tip load from 1 kN to 1.5 kN |
| E2 | add a second static step in which the tip load is 2 kN |
| E3 | sweep the beam height H over 100, 150 and 200 mm |

Record the number of changed lines (added plus removed, from `git diff
--numstat`) for each answer that passes. The [layout rules](skein_spec.md#32-layout-for-small-diffs)
aim at 1 line for E1, at most 6 for E2 and at most 3 for E3. These numbers
are a record, not a gate.
[readability.md §2](readability.md#2-three-edit-tasks) gives the same tasks,
written by a person, for Skein and CalculiX.

---

## 3. The requests

Each request is one short line, and its answer is a small file. Reason: a
long request or a large deck costs context and cache tokens in each round,
and adds no information to the test.

| Name | Request |
|---|---|
| `cantilever` | a steel cantilever 1 m x 50 mm x 100 mm, fixed at one end, 1 kN tip load, report the tip deflection |
| `plate` | an aluminium plate 100 mm x 50 mm x 2 mm, plane stress, one short edge fixed, 10 MPa tension on the other short edge, report the elongation |

The requests are in `REQUESTS` in `scripts/llm/workflow.py`. Euler-Bernoulli
gives 0.381 mm for the cantilever, and σL/E gives 0.0143 mm for the plate
with E = 70 GPa. The driver records the value, but the gate does not use it.

---

## 4. Syntax changes made before the test

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
   them.

Each change removes text that an LLM must write, but each change also adds
a rule to learn. The test decides whether the syntax changes again.

---

## 5. Results

| Date | Model | Request | Path | Pass rate | Median total cost | Median repair rounds | Diff lines E1, E2, E3 |
|---|---|---|---|---|---|---|---|
| | | | | | | | |
