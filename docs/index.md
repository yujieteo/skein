# Skein documentation

Skein is a typed language for numerical simulation. One file describes a
whole model: geometry, mesh, physics in strong form, modelling assumptions,
method, analysis steps and checks. The compiler writes Goose, and Goose
writes C for a native executable or for WebAssembly
([what it is and why](../README.md)).

**Status:** this is a specification. No compiler exists yet. The
[roadmap](roadmap.md) gives the phases.

## The language, in reading order

| Document | What it is | Read it when |
|---|---|---|
| [README](../README.md) | What Skein is and why, in one page. | You want to know whether to look further. |
| [Tutorial](tutorial.md) | The language by example, from a bar under tension to submodels and meshless methods. | You are new to Skein. |
| [Language reference](reference.md) | Each construct, topic by topic, with an example for each. | You want to look up a construct. |
| [Specification](skein_spec.md) | The exact rules: grammar, types, units, regions, assumptions, the intermediate form, the lowering by method and the code generation to Goose. | You want to know why a model did not pass `check`, or you write the compiler. |

## Tools and checks

* [The command line](command_line.md): each command of `skein`, the format
  rules of `skein fmt`, and the TOON output.
* [Checks](checks.md): the catalogue of errors, lint rules and verification
  checks, each with its id and its fix.
* [Samples](../samples/README.md): six complete models.

## The project

* [Implementation notes](implementation.md): how the compiler, the runtime
  and the mesher will work.
* [Roadmap](roadmap.md): the phases, their gates and the risks.
* [Design decisions](decisions.md): the 17 decisions and where the
  documentation implements each one.
* [Tradeoffs](tradeoffs.md): the goals in priority order, the design
  against diverse workloads, and the gaps.
* [Grilling](grilling.md): the design tree and the open questions, each
  with a recommended answer.
* [The LLM test](llm_test.md): how to measure whether an LLM writes Skein
  better than an Abaqus or APDL deck.
