# Skein

Skein is a typed, units-checked language for numerical simulation. It
compiles to [Goose](https://aardappel.github.io/goose/), and Goose compiles
to C and WebAssembly.

```skein
model cantilever
space 3d
units mm_N_MPa_s

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

field u: Vector<Length> on beam

physics solid on beam solves u {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}

method solid: cg {
    mesh: m
    order: 2
}

step load: static {
    bc root: u = 0
    bc tip: sigma*n = (0, 0, -2.5) MPa
    output u, von_mises(sigma)
}
```

**Status:** this repository holds the specification only. No compiler
exists yet. The [roadmap](docs/roadmap.md) gives the phases. Phase 1 is the
front end: the parser, the type and unit checker, `skein fmt` and
`skein check`.

## Why Skein

Abaqus, APDL and NASTRAN decks are hard to write and hard to check. Their
units are a convention that the user must remember. Their physics is a fixed
element library. Their errors appear late, in the solver.

Skein changes these things:

- **Physics in strong form.** You write `div(sigma) = 0` and the
  constitutive law, as in a textbook. The compiler derives the weak form,
  the finite-volume balance or the meshless equations.
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
- **Written for LLMs.** Rust-like declarations, every type explicit, types
  only at declarations, compact TOON tables and output, and stable
  diagnostic ids with a fix.

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
* [The LLM test](docs/llm_test.md)

## License

Apache 2.0. See [LICENSE](LICENSE).
