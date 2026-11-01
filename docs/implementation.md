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

These notes describe the planned architecture for the
[specification](skein_spec.md). The repository already contains a front end,
mesh tools, native solvers and a browser solve for the simplex and hex models.
The SIR pipeline and the viewer remain planned work. See the [command reference](command_line.md) for the
current interface and the [roadmap](roadmap.md) for the planned phases.

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

A dimension is seven exponents, each stored as a count of halves (spec
§6.1). A quantity type is a dimension and a rank. The checker computes the
type of each expression bottom up and compares it with the declared type at
each declaration and on both sides of each equation. The checker does no
inference: an expression has one type, from its parts.

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
function. The emitter writes the kernels into `kernels.goose` of each
runtime package, and `build/test_codegen write` writes them again after a
change to the emitter.

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
| threads | assembly and the matrix-vector products of CG and of the flow GMRES(m) run in Goose workers (`src/solvers/pool.goose`) |

Goose threads share no memory. A worker gets a range of cells as a flat
job, computes the cell matrices, and returns them in a flat result. The main
program adds them into the matrix in cell order. This costs one copy for
each cell matrix, and it needs no locks and no colouring. For CG and for
the GMRES(m) of flow, each worker holds a block of the rows of the matrix,
and computes those rows of each product.

The package `flow` in `src/flow/` solves steady incompressible flow in 2D
(spec §17.2). Its notes:

- The unknowns are scaled (§21.2). `src/flow/gen.goose` gets the
  reference length, mass and time from the `scales` line, and divides each
  input by its reference value. It multiplies each report, `expect` and
  VTU value by its reference value again.
- The Jacobian of a cell is the forward difference of the cell residual.
  Reason: τ_M changes with the velocity, and the difference includes that
  derivative.
- The dofs of each node are together, in the reverse Cuthill-McKee
  sequence. The ILU(0) preconditioner adds the PSPG pressure term when the
  equations have no PSPG. Reason: Taylor-Hood pressure rows have a zero
  diagonal.
- `skein run` makes the mesh, and the generated program reads it from a
  file. Reason: the mesher has recursive functions, and Power of Ten rule
  1 forbids them in a generated program and its runtime.
- `k_omega_sst` (`src/flow/rans.goose`) runs on P1 cells with 5 dofs at
  each vertex. The solve is pseudo-transient continuation with the lumped
  mass. A step takes at most half of k or w at a vertex, the sources of k
  and w use the vertex rule, and the w rows are divided by the largest w
  of a wall. The comments of `rans.goose` and `solve.goose` give the run
  that showed the need for each of these.
- A `bloch` group with `k: 0` makes the flow periodic: the nodes of `high`
  share the dofs of the nodes of `low`. The mesh is half of the rect and
  its mirror. Reason: the mesher splits the boundary segments by the
  points inside, so the two ends of a rect do not get the same points.
- An `interface` between the flow and a `linear_elastic` solid in plane
  strain (`src/flow/fsi.goose`) is a partitioned steady solve: move the
  fluid mesh by the harmonic extension of u, solve the flow, load the
  solid with the reactions of the interface nodes, and relax u with the
  Aitken factor. The mesher meshes the flow and the solid as one shape
  with the solid as an interface, so the two meshes share the interface
  nodes. The solid stiffness is a band with a Cholesky factor. Reason:
  GMRES with ILU(0) did not reach 1e-10 on a cantilever in 5000
  iterations.

One task pool runs the cases of a study (spec §20.5). The runner starts
each case in case order when the memory ledger admits it, and the next
case when a case ends. A case is a separate program, so each case has its
own pool for its assembly and its matrix products.

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

The 2D mesher has a size field for `refine`: the size at a point is
`min(size, h + (grade - 1) d)`, with `d` the distance to the refined
segment or circle. The boundary points are at equal steps of the integral
of 1/size along each curve.

For a cell order above 1, the mesher puts the extra nodes on the true
surface by a projection along ∇φ.

---

## 8. The browser

The browser runs the same Goose code, compiled to WebAssembly:

| Part | In the browser |
|---|---|
| `check`, `fmt`, `lint` | planned, in the editor, with no server |
| mesher | yes, in the page of `skein build --target wasm` |
| solve | yes, for the simplex and hex models and the `fv` and `dg` grids whose memory fits in 15 GB. On the GPU with WebGPU: the explicit steps always, the implicit solves with the control of the page. |
| viewer | yes, WebGL 2: the page is the page of `skein report` (command_line.md §12), with its viewer. |

The browser target compiles the C that Goose writes with `emcc`, as earlier
Goose projects do (spec §18.1). `skein build --target wasm`
(`src/codegen/build.goose`) does these steps:

1. It generates the program of the model at the default values and at
   each end of the range of each param. The number literals that change
   become slots of the program of the page (`src/codegen/web.goose`).
   Reason: the program holds each value as a literal, and the page has no
   Goose compiler. The program of a union of boxes gives the boxes and
   the mesh size to `runtime::mesh_union`, which computes the grid lines.
   Reason: the number of grid lines changes with the size of a box, and
   a list of the lines would change the structure of the program.
2. It writes the program of the page. It holds the parser, the checker,
   the generator, the mesher and the solvers. For the values of the page,
   it generates the program again, takes the literals of the slots from
   it, and runs the program with them. A program that differs in another
   place than a slot stops the run with an error. Thus the run in the page
   gives the records and the VTU file of `skein run --set`, bit for bit.
3. It runs that program natively at the default values, at the corner
   of the ranges with the most dofs, and at the other end of each range
   with the other params at their defaults. Thus each end of each range
   has a run. It measures the peak of each Goose data stack and the peak
   resident memory. A run that fails stops the build with an error.
4. It compiles the C with `emcc` to 64-bit WebAssembly (memory64). Each
   data stack gets a fixed region: 1.5 times the largest peak of the
   runs, 1 MB at least. Reason: WebAssembly memory has no address-space reservation and
   no guard pages (Goose `docs/design/wasm.md`, fixed partitions). A
   canary gap behind each region replaces the guard gap: each end of the
   run checks it, and an overrun stops the run with exit status 4. The
   memory starts at 64 MB and grows to the regions plus 1.5 times the
   largest heap peak (the peak resident memory less the data stacks). A
   growth that the browser refuses stops the run, and the page shows an
   error with that limit. The page cannot catch the kernel OOM killer: a
   host with no free memory can still stop the tab.
5. It writes one HTML page with the module in it: the page of `skein
   report`, with the native run at the default values. Each run has a Web
   Worker and a module of its own. Thus a new value stops the run before
   it, and the memory of a run goes back when it ends.

The page uses WebGPU as goosestrike does: WGSL compute kernels do the
work (`src/gpu/kernels.wgsl`), and Goose sets up and controls the GPU
(`src/gpu/gpu.goose`):

1. The program of the page calls the solves of `gpu` in place of the
   solves of `elements`, `grid` and `meshless` (`GPU_CALLS` in
   `src/codegen/build.goose`). Each solve writes the device row of its
   step, and takes the CPU path when the page has no GPU device. Modal,
   buckling and the static steps of `fv`, `dg` and `rkpm` always use the
   CPU. Only the program of the page imports `gpu`, so `skein run` links
   no GPU code.
2. Goose converts the matrix and the vectors to f32 buffers, binds them to
   the kernels and records the dispatches. The GPU keeps the state
   between two read-backs. Each kernel binds at most 8 storage buffers,
   the default limit of WebGPU: `central_difference` and `dg` scale the
   rows of the matrix by 1/m or the stage weight on the CPU, so that no
   kernel needs them.
3. The C externs of `gpu` call the methods of `src/gpu/gpu.js` through
   `Module.skeinGpu`. The read-back waits for the GPU: JSPI suspends the
   module during the wait. The native profile run links stubs that report
   no GPU.
4. For `ssprk3`, the CPU gives the source and the boundary values of each
   stage (`grid::stage_terms`), and the GPU does the rest of the stage.
5. The implicit GPU path is CG with the Jacobi preconditioner in f32 on
   the GPU, in a defect correction in f64 on the CPU. Each pass solves
   K d = r to a relative residual of 1e-3, adds d to u and computes r
   again, until the residual of the CPU solve holds.

`tests/gpu/gpu.test.mjs` runs each page in Node with the `webgpu` package
(Dawn) and compares the GPU run with the CPU run.

The page runs the simplex models (`tri` or `tet` cells), the hex models
(a box, a sector or a union of boxes) and the `fv` and `dg` grids.
Sample 01 (the cantilever) and the Hisada interposer model
(`examples/package/interposer_warpage.skn`, `hex27` cells) are the
validated models of the page. The interposer at
the default values has 4839 cells and 139 449 dofs, and its memory limit
is 2.9 GB (native peak: 1.7 GB). In Chromium 156, a run meshes and solves
in 214 s on one thread (native: 263 s on all threads). For a value
between the study cases of each of `t_int`, `t_rdl`, `h_bump`, `f_bump`
and `x_off`, the warpage at point A in the page is that of `skein run
--set` to 2e-8. Native runs with a number of threads other than 1 add
their sums in another sequence, so they are not bit for bit.

The front end needs little memory. The solver needs much more, so the
browser limit on memory sets the size of the models that it can solve.
Sample 01 at its default values has 70 244 cells and 319 320 dofs. The
data stacks hold almost all of its memory, and the heap is near zero. Its
memory limit is 3.2 GB. In Chromium 156, a run at the default values
meshes and solves in 50 s on one thread, with a renderer peak of 1.8 GB
(native: 1.6 GB). At the corner L = 400 mm, P = 10 kN (87 524 cells,
401 160 dofs), it takes 87 s and 2.2 GB (native: 1.9 GB).

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
