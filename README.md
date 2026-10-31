# Skein

Skein is a typed, units-checked language for numerical simulation. It
compiles to [Goose](https://aardappel.github.io/goose/), which generates C
for native executables. The WebAssembly target remains planned work.

```skein
model cantilever: 3d, mm_N_MPa_s

scales 40 mm, 2.5 MPa

let beam: Solid = box((0, 0, 0) mm, (200, 20, 40) mm)
let root: Surface = boundary(beam, x == 0 mm)
let tip: Surface = boundary(beam, x == 200 mm)

material steel {
    E: Pressure = 210 GPa
    nu: Real = 0.3
}

mesh m on beam {
    size: 5 mm
}

physics solid on beam solves u: Vector<Length> {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}

step load: static {
    bc root: u = 0
    bc tip: sigma*n = (0, 0, -2.5) MPa
    output u, von_mises(sigma)
}
```

**Status:** the compiler has a front end and native solvers for a limited
set of models. `skein check` parses and type-checks a model. `skein lint` checks
it and reports lint diagnostics. With `--mesh`, it also applies the rules
that need a mesh. `skein verify` creates the mesh and checks its quality.
`skein md` writes a model as a Markdown document. `skein fmt` rewrites a
file in the standard layout. `skein export --format stl` writes the surface of a shape as STL.
`skein run --study` runs each case of a `grid`, `zip` or `sample`
study on parallel workers, and prints the table of the cases. The
[command reference](docs/command_line.md) states the current limits. The
[roadmap](docs/roadmap.md) gives the planned phases.

## Building

The compiler is written in [Goose](https://aardappel.github.io/goose/).
Put `goose` on `PATH`, or set `GOOSE` to its path. If `goose` cannot find
its standard library, set `GOOSE_STDLIB` to the `stdlib` directory of a
Goose source tree at the same commit as the compiler.

```
./build.sh          # writes build/skein
./build.sh test     # builds, then runs each src/<part>/test.goose and checks the generated C
```

The C that Goose generates from each model that `skein run` can run is in
the repository next to the model, for audit. For example,
`examples/femahoy/grid_cantilever.c` is the C of
`examples/femahoy/grid_cantilever.skn`. `./build.sh test` runs each of
these models and fails if the C differs. It needs the Goose version in
`goose_version` in `scripts/generated_c.sh`. After a change to the code
generator, the runtime or Goose, write the C again:

```
./build.sh && scripts/generated_c.sh write
```

The list of these models is in `scripts/generated_c.sh`.

The test of `src/md/` compares the documents of two models with the
golden files in `tests/md/`. After a change to `skein md`, write the
golden files again and read the diff:

```
./build.sh
./build/skein md samples/01_cantilever.skn > tests/md/01_cantilever.md
./build/skein md examples/femahoy/heat_fin.skn > tests/md/heat_fin.md
```

| Path | Purpose |
|---|---|
| `src/skein.goose` | command dispatch |
| `src/syntax/` | syntax tree, lexer, parser and library files |
| `src/check/` | type and unit checks, lint rules and parameter dependencies |
| `src/fmt/` | source layout |
| `src/cli/` | language card for `skein help --llm` |
| `src/md/` | model documents |
| `src/codegen/` | native code generation and study execution |
| `src/runtime/` | hexahedron meshes, solid elements, solvers and VTU output |
| `src/solvers/` | sparse matrices, ILU(0), GMRES and Newton |
| `src/elements/` | triangle and tetrahedron elements for elasticity and heat |
| `src/mesher/` | geometry meshes, mesh checks and STL export |
| `src/shells/` | shell, beam and laminate solver routines, plus code generation for cylindrical shells |
| `src/transient/` | transient heat solver routines |

## Why Skein

Abaqus, APDL and NASTRAN decks are hard to write and hard to check. Their
units are a convention that the user must remember. Their physics is a fixed
element library. Their errors appear late, in the solver.

The [command reference](docs/command_line.md) states which tools exist now.
The language design has these goals:

- **Physics in strong form.** You write `div(sigma) = 0` and the
  constitutive law, as in a textbook. The compiler derives the weak form,
  the finite-volume balance or the meshless equations. A standard library
  of physics templates, in Skein source, covers the common physics.
- **Named modelling assumptions.** `assume small_strain, plane_stress(2 mm)`
  has an exact meaning. The checker rejects assumptions that conflict.
- **Units in the type system.** `E: Pressure = 210 GPa`. The checker rejects
  `E + rho` before any solve.
- **One language for the whole model.** Geometry, mesh, physics, method,
  steps and checks are in one file. Skein does not import meshes. Its mesher
  runs in the browser too.
- **Many methods, one physics.** Continuous and discontinuous Galerkin,
  finite volume, and meshless methods (RKPM, generalized finite
  differences).
- **Advanced analysis from the start.** Element birth and death, submodels,
  CFD stabilization, coupled physics, modal and buckling steps.
- **Verification by default.** Each run checks the rigid-body modes, the
  mesh quality, the residual, the equilibrium, the energy balance and the
  conservation of mass, and prints the result.
- **Written for LLMs, reviewed by people.** Rust-like declarations, every
  type explicit, types only at declarations, compact TOON tables and
  output, and stable diagnostic ids with a fix. One item on each line and
  no column alignment, so a small change gives a small diff.
- **Dimensions done right.** Each model has a `scales` line, and the
  solver works in scaled variables. The planned `skein check --groups`
  command prints `Re`, `Pr` and the other groups of the model. A `param` and a `study` sweep a
  model, with the cases on parallel workers.

Simulation is the core. Parametric parts and export come later, and the
[tradeoffs](docs/tradeoffs.md) rank all the goals.

## Documentation

* [Documentation index](docs/index.md)
* [Tutorial](docs/tutorial.md)
* [Language reference](docs/reference.md)
* [Specification](docs/skein_spec.md)
* [Command line](docs/command_line.md)
* [Checks](docs/checks.md)
* [Samples](samples/README.md)
* [Roadmap](docs/roadmap.md)
* [Design decisions](docs/decisions.md)
* [Tradeoffs](docs/tradeoffs.md)
* [Grilling questions and answers](docs/grilling.md)
* [The LLM test](docs/llm_test.md)

## License

Apache 2.0. See [LICENSE](LICENSE).
