# The Skein compiler: implementation notes

<details>
<summary>Table of contents</summary>

- [1. Pipeline](#1-pipeline)
- [2. Front end](#2-front-end)
- [3. The checker](#3-the-checker)
- [4. From strong form to SIR](#4-from-strong-form-to-sir)
- [5. Methods and kernels](#5-methods-and-kernels)
- [6. The runtime](#6-the-runtime)
- [7. The mesher](#7-the-mesher)
- [8. The browser](#8-the-browser)
- [9. Results and the viewer](#9-results-and-the-viewer)
- [10. Testing and the verification suite](#10-testing-and-the-verification-suite)

</details>

These notes say how the compiler and the runtime will implement the
[specification](skein_spec.md). They are a plan: no code exists yet. The
[roadmap](roadmap.md) gives the sequence.

---

## 1. Pipeline

The compiler and the runtime are Goose programs. The compiler has these
passes:

| Pass | Input | Output | Phase |
|---|---|---|---|
| lex | `.skn` text | tokens with line and column | 1 |
| parse | tokens | AST, with error recovery | 1 |
| resolve | AST | names bound to declarations | 1 |
| check | AST | typed AST, dimensions, regions, diagnostics | 1 |
| fmt | AST and comments | text | 1 |
| lint | typed AST | diagnostics | 2 |
| lower | typed AST | SIR | 2 |
| discretize | SIR and method | discrete IR (kernels and plan) | 2 |
| emit | discrete IR | `.goose` source | 2 |

`skein` calls the Goose compiler as a separate program for `build` and
`run`. Reason: the Goose compiler is large, and a separate program keeps the
`check` path small and fast.

---

## 2. Front end

The lexer works on the bytes of the file and keeps the line and the column
of each token. Statement ends come from the rule of spec §2.2, so the parser
sees an explicit `NL` token.

The parser is recursive descent, one function for each production of spec
§3. After a syntax error, it skips to the next line that starts with a
keyword of the top level, and it continues. Thus one `check` reports the
errors of each item.

`fmt` works on the tokens, so that the comments and the parentheses stay.
The AST gives the role of each token. `fmt` applies the rules of
[command_line.md §3.1](command_line.md#31-format-rules), then parses its
output and compares the AST and the comments with those of the input
(`src/fmt/fmt.goose`). The test of `fmt` formats a model, formats the
result again, and compares the two.

---

## 3. The checker

### 3.1 Dimensions

A dimension is seven `i8` exponents. A quantity type is a dimension and a
rank. The checker computes the type of each expression bottom up and
compares it with the declared type at each declaration and on both sides of
each equation. The checker does no inference: an expression has one type,
from its parts.

A unit is a dimension and an `f64` scale. The checker converts each literal
to SI when it reads it. `degC` is the one affine unit (spec §6.3).

### 3.2 Regions

The checker keeps the CSG tree of each geometry value and computes its
topological dimension. For containment (spec §8.4), it first compares the
trees. When the trees do not decide, it samples the distance functions on a
grid of the bounding box and reports the test as "by sampling" in the
diagnostic.

### 3.3 Suggested fixes

Each diagnostic has a fix. The fixes come from a small set of strategies:

| Strategy | Used for |
|---|---|
| nearest name by edit distance | `unknown_name`, `unknown_unit`, `assumption_unknown` |
| the dimension that makes a sum valid | `unit_mismatch` |
| the template of a missing part | `method_missing`, `eq_count` |
| the line that conflicts | `bc_conflict`, `duplicate_name` |

---

## 4. From strong form to SIR

The lowering does these steps for each physics:

1. Apply the assumptions (spec §11): replace `strain`, reduce plane stress,
   add the axisymmetric terms.
2. Expand each `let` into the equations, so that each equation is a tree
   over the unknowns, their derivatives, the coordinates and time.
3. Put each equation in balance-law form (spec §16.1): collect the `div`
   terms into the flux, split the flux into advective and diffusive parts,
   and collect the time derivatives.
4. Differentiate each equation symbolically with respect to each unknown and
   its gradient. This gives the exact Jacobian for Newton's method.
5. Mark each equation `linear` and `symmetric`.

The symbolic algebra works on a small set of tensor operations. It
simplifies by constant folding, by zero and identity rules, and by common
subexpressions. It does not need a general computer algebra system.

---

## 5. Methods and kernels

Each method kind is a module of the compiler that maps SIR to a discrete IR
(spec §17.6). The kernels are straight-line code over fixed-size arrays:
the number of nodes of a cell and the number of quadrature points are
constants for each kernel. This lets the C compiler vectorize them, and
`simd fn` makes versions for the wider instruction sets.

A kernel longer than 60 lines is split into functions by blocks of the cell
matrix (spec §18.3). Each block function is straight-line code, and one
short function calls them in order. A test measures each generated
function.

The plan is the part with loops. It is short and the same for each model,
except the names of the kernels.

---

## 6. The runtime

`runtime` is a Goose package (spec §18.4). It uses the Goose data stacks for all arrays.
A mesh is a set of flat arrays (coordinates, connectivity, region tags).
A sparse matrix is CSR in three arrays. The solvers use no other storage.

| Part | Notes |
|---|---|
| CG | with Jacobi, ILU(0) or AMG as preconditioner |
| GMRES(m) | restarted, right preconditioned |
| Newton | with a backtracking line search and the load ramp of `static(ramp:)` |
| LOBPCG | for `modal` and `buckling` |
| AMG | smoothed aggregation |
| threads | assembly and the matrix-vector product run in Goose workers over a colouring |

Goose threads share no memory. A worker gets a range of cells as a flat
job, computes the cell matrices, and returns them in a flat result. The main
program adds them into the matrix. This costs one copy for each cell
matrix, and it needs no locks.

In phase 3, one work-stealing task pool runs the cases of a study and the
parallel parts of each case (spec §20.5). Each worker has a queue of
tasks. An idle worker takes a task from the queue of another worker.

---

## 7. The mesher

The mesher implements spec §9.2 from the published papers. It must not copy
code or text from Gmsh, Triangle, TetGen or CGAL.

The 3D path:

1. Build an octree over the bounding box. The cell size follows `size`,
   `refine` and `grade`.
2. Find the feature curves from the CSG tree, and place protecting balls on
   them.
3. Fill the octree with the body-centred cubic lattice of isosurface
   stuffing, and cut the lattice at φ = 0.
4. Smooth by optimal Delaunay triangulation, and refuse a move that makes a
   cell worse.
5. Measure each cell. Give the worst cells to the check `mesh_quality`.

For a cell order above 1, the mesher puts the extra nodes on the true
surface by a projection along ∇φ.

---

## 8. The browser

The browser runs the same Goose code, compiled to WebAssembly:

| Part | In the browser |
|---|---|
| `check`, `fmt`, `lint` | yes, in the editor, with no server |
| mesher | yes |
| solve | yes, for models that fit the memory of the browser |
| viewer | yes, WebGL |

The browser target compiles the C that Goose writes with `emcc`, as earlier
Goose projects do (spec §18.1). The front end needs little memory. The
solver needs much more, so the browser limit on memory sets the size of
the models that it can solve.

---

## 9. Results and the viewer

The run writes VTU files, one for each output time, which ParaView reads.
The WebGL viewer reads the same files. It shows:

- contours of any output, on the surface and on a cut plane;
- the deformed shape, with a scale factor;
- probes: the value at a point that the user clicks, and the line plot
  along a line;
- the active cells of each step, for birth and death.

---

## 10. Testing and the verification suite

The tests have three levels.

1. **Front end.** Each error, lint rule and fix in [checks.md](checks.md)
   has a test file that must give exactly that diagnostic at that line and
   column. Each sample must pass `check` and `fmt --check`. A `fmt` test
   changes one value in each sample and checks that `fmt` then changes no
   other line (spec §3.2).
2. **Methods.** A manufactured solution for each method and each element
   order, with `skein verify convergence`. The test fails if the observed
   order is lower than the expected order by more than 0.3.
3. **Benchmarks.** Published and analytic problems, not tied to one
   company:

| Benchmark | Physics | Reference |
|---|---|---|
| Patch test | linear elasticity | exact linear field |
| Cantilever tip deflection | linear elasticity | Euler-Bernoulli and Timoshenko beam theory |
| Plate with a hole | plane stress | Kirsch solution |
| Cook's membrane | near-incompressible, finite strain | Cook (1974), and the converged values in the literature |
| Thick cylinder under pressure | axisymmetric | Lamé solution |
| Lid-driven cavity, Re 100 and 1000 | Navier-Stokes | Ghia, Ghia and Shin (1982) |
| Flow past a cylinder, 2D-1 and 2D-2 | Navier-Stokes | Schäfer and Turek (1996) |
| Heat conduction in a slab | transient heat | the series solution |
| Euler buckling of a column | buckling | Euler load |
| Free vibration of a beam | modal | beam theory frequencies |
