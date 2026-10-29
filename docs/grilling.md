# Skein design grilling

This file holds the questions of the design review with the captain.
Firstmate put them to the captain, one round at a time. The captain
answered all rounds on 2026-10-30. The [answers](#answers) are at the end,
and the spec, the tradeoffs and the samples follow them. Each question has a
recommended answer and the reason for it. The captain can accept the
recommendation, choose another option, or give a new answer.

The questions come from the [tradeoffs](tradeoffs.md), the
[LLM test](llm_test.md) and the open items of the review of this spec. The
numbers (G0.1, G1.1, …) are stable, so an answer can name its question.

## How to use this file

1. Ask the questions of round 0 first, then round 1. Record each answer
   in the table at the end.
2. Do not start a round before the answers of the round before it are
   recorded, because later questions depend on them.
3. After the last round, the spec changes to match the answers, in a new PR.

## The design tree

```
Skein
├── G0 The ranking of the goals
│   ├── G0.1 simulation as the core
│   ├── G0.2 which goals are protected
│   ├── G0.3 the rank of LLM-native syntax
│   └── G0.4 which goals to cut first
├── G1 Scope and core language
│   ├── G1.1 simulation first, or CAD and simulation as equals
│   ├── G1.2 user functions
│   ├── G1.3 one namespace per model
│   └── G1.4 more than one file per model
├── G2 Dimensions
│   ├── G2.1 models with no units
│   ├── G2.2 default scales
│   └── G2.3 reuse of similar cases
├── G3 Parameters, studies and concurrency
│   ├── G3.1 only `param` can change
│   ├── G3.2 parallel cases before parallel assembly
│   ├── G3.2b a better method in the style of Go
│   ├── G3.3 adjoint gradients
│   └── G3.4 statistics in the study output
├── G4 Workload gaps
│   ├── G4.1 standard-library physics (shells, beams, turbulence)
│   ├── G4.2 radiation
│   ├── G4.3 shock capturing
│   ├── G4.4 element families (H(curl), H(div))
│   ├── G4.5 complex numbers
│   ├── G4.5b other number systems and symmetry groups
│   ├── G4.6 adaptive time steps
│   └── G4.7 processes in a step (reinitialization)
├── G5 Parts and manufacture
│   ├── G5.1 more geometry operations
│   ├── G5.2 exact STEP export
│   ├── G5.3 nonlinear tolerance stack-ups
│   └── G5.4 the list of processes
└── G6 LLM, reliability and the web
    ├── G6.1 the size of Skein text against APDL
    ├── G6.2 a language card for LLMs
    ├── G6.3 the Goose WebAssembly dependency
    └── G6.4 the exceptions to the Power of Ten
```

---

## Round 0: the ranking of the goals

Ask this round first. The [tradeoffs](tradeoffs.md#0-the-goals-in-priority-order)
rank nine goals. The answers here decide what the other rounds may cost.

**G0.1** Is simulation the core, with all other goals layered on it?

- A. Yes. Strong-form physics, methods, mesher and solvers come first, and no
  other goal may weaken them.
- B. No. Skein is a general engineering language, and simulation is one
  use of it.

Recommended: **A**. Reason: you said that Skein is "primarily for
simulation", and that there are too many goals. One core keeps the
tradeoffs decidable.

**G0.2** Which goals are protected (never weakened for a later goal)?

- A. Simulation core, units and static analysis with verification, and
  LLM- and AXI-native syntax (ranks 1 to 3).
- B. Only the simulation core and the units.
- C. All nine goals equally.

Recommended: **A**. Reason: units and verification protect each result, and
LLM-native syntax is the reason for a new language. C cannot be built,
because the goals conflict (for example, exact CAD against a small compiler).

**G0.3** When LLM-native syntax and strict units conflict, which wins?

- A. Units. A unit on each literal stays, even when it costs tokens.
- B. Tokens. Allow bare numbers in the unit system of the model.

Recommended: **A**. Reason: a bare number in an assumed unit is the fault
that Skein exists to remove. The token cost is a few per literal.

**G0.4** If the schedule is late, which goals go first?

- A. DFM (rank 9), then text-to-CAD export beyond STL (rank 8).
- B. Web output (rank 6).
- C. Dimensional analysis (rank 5).

Recommended: **A**. Reason: DFM and exact CAD export serve manufacture, which
is outside the simulation core. The web output and the dimensional analysis
serve each simulation user.

---

## Round 1: scope and core language

**G1.1** Is Skein a simulation language that also describes parts, or a CAD
language and a simulation language as equals?

- A. Simulation first. Parts are parametric geometry for the simulation, the
  studies and the DFM checks. No constraint sketches.
- B. Equals. Add constraint sketches, features and a B-rep kernel.

Recommended: **A**. Reason: you said that Skein is "primarily for
simulation". Option B needs a B-rep kernel, which is a project of its own
size. A still gives text-to-CAD for parts that CSG can describe.

**G1.2** Can a user declare functions?

- A. No functions (the spec today).
- B. Pure functions of one expression: `fn lame_mu(E: Pressure, nu: Real) -> Pressure = E/(2*(1 + nu))`.
  No recursion, no loops, no state.
- C. Functions with bodies.

Recommended: **B**. Reason: a function of one expression keeps rules 1 to 3
of the Power of Ten by construction, because it cannot loop or recurse. It
also cuts the size of models that repeat a law, which helps G6.1. C breaks
the declarative core.

**G1.3** Keep one namespace for each model?

- A. Yes. `sigma` from a physics is usable in a step with no prefix.
- B. A namespace for each physics: `solid.sigma`.

Recommended: **A**. Reason: models are short, and the prefix costs tokens
on each use. Two physics with a `sigma` each must name them `sigma_s` and
`sigma_f`, and the checker says so.

**G1.4** Can a model use more than one file?

- A. No. One file is one model (today). `results` reads another model.
- B. `use "materials.skn"` imports `material`, `param` and `let` items from
  a library file. No physics, no steps.

Recommended: **B**, in phase 2. Reason: a company keeps one materials file.
Copies in each model drift. Physics and steps stay in the model, so a
reviewer still reads one file for the analysis.

---

## Round 2: dimensions

**G2.1** Allow models with no units (`units nondim`)?

- A. Yes, declared in the header, and the check report says so.
- B. No. Each model has dimensions. A textbook PDE uses `Real` scales.

Recommended: **A**. Reason: general PDEs from mathematics have no units,
and you asked for general PDEs. The header makes the choice visible.

**G2.2** When a model has no `scales` line, scale by the base units of the
unit system?

- A. Yes (the spec today).
- B. No. Require `scales` in each model with a flow physics.

Recommended: **A**. Reason: a structural model in `mm_N_MPa_s` is already
well scaled. Make `scales` a lint suggestion for flow models, not a rule.

**G2.3** Reuse the solution of a similar case in a study (equal π groups)?

- A. Yes, on by default, marked `similar_to` in the table.
- B. No. Each case runs (the spec today).

Recommended: **B**. Reason: the similarity test must include each
dimensional constant of the SIR, also the literals in geometry, `bc` values
and step arguments. A missed constant gives a wrong result with no error.
For example, a fixed channel with a varied cylinder diameter changes the
blockage ratio. Add the reuse later only with that full test.

---

## Round 3: parameters, studies and concurrency

**G3.1** Only a `param` can change from a study or from `--set`?

- A. Yes.
- B. Any `let` can change.

Recommended: **A**. Reason: a reviewer finds each input in the `param`
lines, and the dependency analysis knows which stages each one reaches.

**G3.2** When a study has more cases than threads, run cases in parallel
before assembly in parallel?

- A. Yes (§20.5 of the spec).
- B. Always parallel assembly, one case at a time.

Recommended: **A**. Reason: independent cases scale almost perfectly. The
runner still gives all threads to one case when it is alone.

**G3.2b** (asked in round 4) The captain asked for a better method in the
style of Go. Use one shared work-stealing task pool?

- A. Yes, in phase 3. Cases first until then. Check whether Goose gives
  atomics, and use a lock on each queue if it does not.
- B. Keep the two levels of threads (cases, then assembly).

Recommended: **A**. Reason: the Go scheduler gives each thread a queue and
lets an idle thread take work from another. A study of a few large cases
then keeps all threads busy.

**G3.3** Add adjoint gradients for optimization?

- A. Yes, in phase 5, from the SIR. Gradient-free `minimize` before that.
- B. Not in v1.

Recommended: **A**. Reason: shape and topology optimization need gradients.
The SIR already holds the exact Jacobian, which is most of an adjoint.

**G3.4** Compute statistics (mean, deviation, quantiles) of a `sample` study
in the tool?

- A. No. The study writes a TOON table. Other tools compute statistics.
- B. Yes, for `sample` studies.

Recommended: **A**. Reason: each feature that nobody uses costs upkeep.
Add B only after a real study shows the need.

---

## Round 4: workload gaps

**G4.1** Add standard-library physics written in Skein (shells, beams,
turbulence models)?

- A. Yes. `physics shell: std.reissner_mindlin(steel) on skin solves w, theta`.
  The source of each one is Skein that a user can read and copy.
- B. No. Users write all physics in strong form.

Recommended: **A**. Reason: nobody writes a shell or a k-ω SST model from
memory. It also helps G6.1. It needs G1.2 B or a template form.

**G4.2** Add radiation between surfaces?

- A. A `radiation` item with view factors, phase 5.
- B. Out of scope.

Recommended: **A**. Reason: thermal models of electronics and furnaces need
it, and it is a non-local term that the strong form cannot express.

**G4.3** Add shock capturing for compressible flow?

- A. `limiter` for `fv` and `shock_capturing: entropy_viscosity` for `dg`.
- B. Compressible flow without shocks only.

Recommended: **A**. Reason: compressible CFD is in the scope of "any
engineering simulation", and shocks are the common case there.

**G4.4** Add element families for H(curl) and H(div)?

- A. `family: nedelec | raviart_thomas` in `cg`.
- B. No.

Recommended: **A**. Reason: electromagnetics and mixed Darcy flow are wrong
with nodal elements. The setting is small. The elements are phase 5.

**G4.5** Add complex numbers?

- A. A type `Complex<D>` and a `harmonic(freq:)` procedure.
- B. No. Use two real fields.

Recommended: **A**. Reason: acoustics and electromagnetics in the frequency
domain are complex. Two real fields double the size of each model.

**G4.5b** (asked in round 5) The captain asked for general representation
theory, "Complex or Quaternions or whatever". Which number systems?

- A. A closed list of the systems in practice: `Real`, `Complex`,
  `Quaternion`, `Dual`.
- B. User-defined number systems.

Recommended: **A**. Reason: each system needs its own solver arithmetic.
User systems would need user code in the solver. The captain added
symmetry groups, also as a closed list: cyclic sectors, dihedral, and Bloch
and Floquet periodic.

**G4.6** Add adaptive time steps with error control?

- A. `transient(adapt: { rtol, atol })` with BDF of variable order.
- B. Fixed `dt` only.

Recommended: **A**. Reason: stiff ODEs from chemistry and control need it,
and you asked for general ODEs.

**G4.7** How does a model state a process that is not an equation, such as
level-set reinitialization?

- A. A `process` statement in a step, from a fixed list in the spec.
- B. Out of scope for v1.

Recommended: **B** for v1, then A after a real model shows which processes
are necessary. Reason: each process in the list is an exception to the
strong form.

---

## Round 5: parts and manufacture

**G5.1** Add more geometry operations for text-to-CAD?

- A. `shell(s, t)` (constant-thickness offset), `loft`, `pattern` (linear
  and polar arrays).
- B. Only the operations in the spec today.

Recommended: **A**. Reason: housings and brackets, the common parts of
text-to-CAD, need shells and patterns. Each one is a distance function.

**G5.2** When does Skein export exact STEP?

- A. Phase 4 tessellated STEP, phase 5 exact B-rep for the analytic
  primitives.
- B. STL and tessellated STEP only.

Recommended: **A**. Reason: CAM tools work better with exact faces. The
analytic primitives are a closed set, so the B-rep is bounded work.

**G5.3** Support nonlinear tolerance stack-ups?

- A. `method: monte_carlo` that uses a `sample` study.
- B. Linear (worst case and RSS) only.

Recommended: **B** for v1. Reason: most stack-ups are linear, and the study
can already sample parameters with tolerances. Add A when a real part needs
it.

**G5.4** Is the list of processes right: `cnc_mill`, `injection_mould`,
`casting`, `sheet_metal`, `additive`, `turning`?

Recommended: **yes**, and add others only for a real part. Reason: these
cover most mechanical parts. Each process needs its own checks.

---

## Round 6: LLM, reliability and the web

**G6.1** The Skein reference answers are 1.4 to 2.9 times the size of the
APDL and Abaqus decks in bytes. The LLM test will probably fail on tokens.
What now?

- A. Keep the syntax until phase 1 measures it (the decision for this PR).
- B. Apply some candidate changes now: the field in `solves`, a one-line
  header, standard-library physics (G4.1), method defaults, functions (G1.2).
- C. Change the gate: compare total cost (tokens plus repair rounds), not
  only the first answer.

Recommended: **A and C**. Reason: the test needs `skein check`, so measure
first. But a Skein answer holds the equations, the units and the checks,
which a deck does not. A gate on tokens alone ignores the repairs that an
LLM makes to a deck that fails late in the solver.

**G6.2** Add a compact language card for LLMs (`skein help --llm`, about 2000
tokens)?

- A. Yes. It is the context of the LLM test and of each LLM session.
- B. No. LLMs read the tutorial and the reference.

Recommended: **A**. Reason: the tutorial and the reference are large. A card
with the grammar, the keywords and one example costs much less context.

**G6.3** Goose on WebAssembly is a proposal, not a feature. What if it does
not arrive?

- A. Ask the Goose author for the plan now. Build phase 1 native only.
- B. Also plan a fallback: the Skein runtime for the browser in another
  language.

Recommended: **A**. Reason: the front end needs only small fixed data
stacks, which the proposal already covers. A fallback doubles the runtime
work before the need is real.

**G6.4** Accept the two exceptions to the Power of Ten: generated kernels
longer than 60 lines, and one namespace per model?

Recommended: **yes**. Reason: the kernels are straight-line code with no
branches, from one template each. The namespace is G1.3.

---

## Answers

The captain answered on 2026-10-30. The words of the captain are exact.
Where an answer differs from the recommendation, the spec follows the
answer.

| Round | Words of the captain |
|---|---|
| 0 | "Q1. A. 1/3/ mix of 2 and 5 I realise dimensional analysis / units is very subtle. I want strong typing for real problems requiring dimensions, but yet I want dimensional analysis / parametric sweep when it counts. All other goals... Q3: A. Cut A." |
| 1 | "No constraint sketches, A. B. A. B." |
| 2 | "A. C (LLMs can easily write scales). B." |
| 3 | "A. For Q2 are there better methods I like Go's concurrency model, but A. A. A." |
| 4 | "A. A. A. A. A. For Q5 it should support general representation theory (Complex or Quaternions or whatever) but yes support complex numbers. Yes. B." |
| 5 | "A for those in practice. Symmetry groups as well. Closd lists. A. B. No exact B-rep. A. For Q4." |
| 6 | "Q0 correct. Make the candidate cuts now. Change the test and comparison should be total cost. LLM-native syntax is a good goal. Yes. A. No, we can convert the Goose to wasm similarly to our previous projects. Try to keep it less than 60 lines somehow tere should be some split up. I also want to minimize namespaces?" |
| after | "It should also maximize for human readable diffs, right now LLMs generate huge diffs that are hard to read, 1000s of lines long. OK I think it is fine now. The spec is fine." |

| Question | Answer | Where |
|---|---|---|
| G0.1 | A. Simulation is the core. | [tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order) |
| G0.2 | Protected: the simulation core, LLM-native syntax, and one merged goal of units, static analysis and dimensional analysis, with parameter sweeps where they count. | [tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order) |
| G0.3 | A. Units win over tokens. | [tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order) |
| G0.4 | A. Cut DFM first, then CAD export beyond STL. | [tradeoffs §0](tradeoffs.md#0-the-goals-in-priority-order) |
| G1.1 | A. Simulation first. No constraint sketches. | [tradeoffs §4](tradeoffs.md#4-parametric-parts-text-to-cad-and-export) |
| G1.2 | B. Pure functions of one expression. | spec [§7.4](skein_spec.md#74-functions) |
| G1.3 | A. One namespace for each model. | spec [§4](skein_spec.md#4-names-and-scope) |
| G1.4 | B. `use` reads materials, params, lets and functions, in phase 2. | spec [§4.1](skein_spec.md#41-use-phase-2) |
| G2.1 | A. `units nondim`, declared in the header. | spec [§6.5](skein_spec.md#65-unit-systems) |
| G2.2 | C (an answer that the question did not list). A `scales` line is required in each model with dimensions. | spec [§21.2](skein_spec.md#212-scales) |
| G2.3 | B. No reuse of similar cases. Each case runs. | [tradeoffs §8](tradeoffs.md#8-strict-dimensions-and-dimensionless-work) |
| G3.1 | A. Only a `param` can change. | spec [§20.1](skein_spec.md#201-param) |
| G3.2 | A, with G3.2b. | spec [§20.5](skein_spec.md#205-threads-for-a-case-and-threads-for-a-study) |
| G3.2b | A. One work-stealing task pool in phase 3. Check Goose atomics, and use a lock on each queue if there are none. | spec [§20.5](skein_spec.md#205-threads-for-a-case-and-threads-for-a-study) |
| G3.3 | A. Adjoint gradients in phase 5. Gradient-free `minimize` before. | [roadmap](roadmap.md) |
| G3.4 | A. No statistics in the tool. | [tradeoffs §2](tradeoffs.md#2-the-design-against-the-workloads) |
| G4.1 | A. Standard-library physics in Skein source. | spec [§4.2](skein_spec.md#42-the-standard-library) |
| G4.2 | A. Radiation with view factors, phase 5. | [roadmap](roadmap.md) |
| G4.3 | A. Limiters for `fv`, entropy viscosity for `dg`. | spec [§13](skein_spec.md#13-methods) |
| G4.4 | A. Nédélec and Raviart-Thomas, phase 5. | spec [§13](skein_spec.md#13-methods) |
| G4.5 | Complex numbers, in a design that generalizes (G4.5b). | spec [§5.1](skein_spec.md#51-scalar-types) |
| G4.5b | A. The closed list `Real`, `Complex`, `Quaternion`, `Dual`. Symmetry groups too, also a closed list: cyclic, dihedral, Bloch and Floquet. | spec [§5.1](skein_spec.md#51-scalar-types), [§12.6](skein_spec.md#126-symmetry-groups) |
| G4.6 | A. Adaptive time steps, variable-order BDF. | spec [§14.2](skein_spec.md#142-procedures) |
| G4.7 | B. No `process` statement in v1. | [tradeoffs §2](tradeoffs.md#2-the-design-against-the-workloads) |
| G5.1 | A. `shell`, `loft`, `pattern`. | spec [§8.1](skein_spec.md#81-distance-functions) |
| G5.2 | B. No exact B-rep. STL and tessellated STEP only. | spec [§22.5](skein_spec.md#225-export-for-manufacture) |
| G5.3 | A. Monte Carlo stack-ups through a `sample` study. | spec [§22.4](skein_spec.md#224-tolerances-and-stack-ups) |
| G5.4 | Yes. The list of processes stays. | spec [§22.2](skein_spec.md#222-manufacture) |
| G6.1 | B and C. Make the candidate cuts now, and compare the total cost. | [llm_test.md §6](llm_test.md#6-syntax-changes-made-before-the-test) |
| G6.2 | A. `skein help --llm`. | [command_line.md §1](command_line.md#1-commands) |
| G6.3 | No WebAssembly risk and no fallback: Goose writes C, and `emcc` compiles it, as in earlier Goose projects. | spec [§18.1](skein_spec.md#181-why-goose) |
| G6.4 | No exception for kernels: split each generated function to at most 60 lines. One namespace stays, and the standard library and `use` files enter it with no prefix. This is the one deviation from the Power of Ten (rule 6). | [tradeoffs §6](tradeoffs.md#6-reliability-the-power-of-ten) |
| after | Small human-readable diffs, under the LLM-native goal. | spec [§3.2](skein_spec.md#32-layout-for-small-diffs), [command_line.md §3.1](command_line.md#31-format-rules) |
