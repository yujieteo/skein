# Skein

Skein is a typed, units-checked language for numerical simulation. It
compiles to [Goose](https://aardappel.github.io/goose/), and Goose compiles
to C and WebAssembly.

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

**Status:** this repository holds the specification and the start of the
compiler. `skein check` parses and type-checks a model. `skein fmt` rewrites
a file in the standard layout. The [roadmap](docs/roadmap.md) gives the phases.
Phase 1 is the front end: the parser, the type and unit checker,
`skein fmt` and `skein check`.

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

| Directory | Contents |
|---|---|
| `src/skein.goose` | the command line |
| `src/syntax/` | the syntax tree (`tree.goose`), the lexer (`lexer.goose`) and the parser (`parser.goose`, `units.goose`) |
| `src/fmt/` | `skein fmt` (`fmt.goose`) |
| `src/runtime/` | the solver that generated code calls (`api.goose`): a 27-node hexahedron box mesher (`mesh.goose`), CG with Jacobi (`solve.goose`) and VTU output (`vtu.goose`) |
| `src/cli/` | the language card of `skein help --llm` (`card.goose`) |
| `src/codegen/` | `skein run` for one linear static box: the Goose code generator (`codegen.goose`) and the build and run (`run.goose`). It needs the runtime in `src/runtime/`. |
| `src/check/` | the type and unit checker and `skein check` (`check.goose`) |

## Why Skein

Abaqus, APDL and NASTRAN decks are hard to write and hard to check. Their
units are a convention that the user must remember. Their physics is a fixed
element library. Their errors appear late, in the solver.

Skein changes these things:

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
  CFD stabilization, restarts, coupled physics, modal and buckling steps.
- **Verification by default.** Each run checks the rigid-body modes, the
  mesh quality, the residual, the equilibrium, the energy balance and the
  conservation of mass, and prints the result.
- **Written for LLMs, reviewed by people.** Rust-like declarations, every
  type explicit, types only at declarations, compact TOON tables and
  output, and stable diagnostic ids with a fix. One item on each line and
  no column alignment, so a small change gives a small diff.
- **Dimensions done right.** Each model has a `scales` line, and the
  solver works in scaled variables. `skein check --groups` prints `Re`,
  `Pr` and the other groups of the model. A `param` and a `study` sweep a
  model, with the cases on parallel workers.

Simulation is the core. Parametric parts, export for manufacture and
design-for-manufacture checks come later, and the
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
