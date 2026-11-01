# The Skein Language — Specification (v0 draft)

<details>
<summary>Table of contents</summary>

- [1. Core concepts](#1-core-concepts)
- [2. Lexical structure](#2-lexical-structure)
- [3. Grammar](#3-grammar)
- [4. Names and scope](#4-names-and-scope)
- [5. Types](#5-types)
- [6. Dimensions and units](#6-dimensions-and-units)
- [7. Expressions and operators](#7-expressions-and-operators)
- [8. Geometry and regions](#8-geometry-and-regions)
- [9. Meshes and point clouds](#9-meshes-and-point-clouds)
- [10. Fields, materials and physics](#10-fields-materials-and-physics)
- [11. Modelling assumptions](#11-modelling-assumptions)
- [12. Boundary conditions and loads](#12-boundary-conditions-and-loads)
- [13. Methods](#13-methods)
- [14. Steps](#14-steps)
- [15. Results, submodels and expectations](#15-results-submodels-and-expectations)
- [16. The intermediate form (SIR)](#16-the-intermediate-form-sir)
- [17. Lowering by method](#17-lowering-by-method)
- [18. Code generation to Goose](#18-code-generation-to-goose)
- [19. Static semantics checklist](#19-static-semantics-checklist)
- [20. Parameters and studies](#20-parameters-and-studies)
- [21. Dimensional analysis and scaling](#21-dimensional-analysis-and-scaling)
- [22. Export](#22-export)
- [Appendix A. Built-in names](#appendix-a-built-in-names)
- [Appendix B. References](#appendix-b-references)

</details>

This document gives the exact rules of Skein. It is for implementors, and
for users who want to know why a program did not pass `skein check`. The
[tutorial](tutorial.md) and the [reference](reference.md) teach the same
language with less precision. The [check catalogue](checks.md) lists every
diagnostic that this document names.

The words **must**, **must not** and **may** have their usual meaning for
specifications. A program that breaks a **must** rule does not pass
`skein check`. The diagnostic id that the checker reports is in brackets,
for example [`unit_mismatch`].

---

## 1. Core concepts

A Skein **program** is one source file that describes one **model**. The
model contains all of these parts, in one language:

| Part | Keyword | Section |
|---|---|---|
| Header: name, space dimension, unit system, on one line | `model` | §3 |
| Library files and the standard library | `use` | §4.1, §4.2 |
| Inputs that a study or the command line can change | `param` | §20 |
| Constants, geometry, regions, results of other models | `let` | §5, §8, §15 |
| Functions of one expression | `fn` | §7.4 |
| Tables of data | `table` | §5.6 |
| Materials | `material` | §10.2 |
| Discretizations of a region | `mesh`, `points` | §9 |
| Equations in strong form, with assumptions, and their unknowns | `physics … solves` | §10.1, §10.3, §11 |
| Numerical method for each physics, with defaults | `method` | §13 |
| Boundary conditions and loads | `bc` | §12 |
| Periodic symmetry | `symmetry` | §12.6 |
| The analysis sequence | `step` | §14 |
| Checks on the results | `expect` | §15.3 |
| Reference scales for the solver | `scales` | §21 |
| Parameter studies | `study` | §20 |
| Turned-off checks, each with a reason | `allow` | [checks.md](checks.md) |

Four rules shape the whole language.

1. **Every declaration has an explicit type.** The checker does no type
   inference. A type appears only at a declaration. An equation or an
   expression has no type annotation.
2. **Every quantity has a dimension.** The type of a quantity is its
   dimension, for example `Pressure`. The literal carries the unit, for
   example `210 GPa`. The checker rejects a sum of two different
   dimensions.
3. **Physics is strong form.** The user writes the partial differential
   equation as in a textbook. The compiler derives the weak form, the
   finite-volume balance, or the collocation equations from it. The user
   never writes a weak form.
4. **A program describes, it does not compute.** Skein has no loops, no
   pointers, no allocation and no mutable variables. A function is one
   expression (§7.4). The compiler makes all
   the loops, in Goose. Thus a Skein program cannot make a memory fault.

### 1.1 Pipeline

```
.skn source ──parse──> AST ──check──> typed AST ──lower──> SIR ──method──> discrete IR ──emit──> .goose
                                                                                               │
                                     goose (Goose compiler) ──> C ──> cc: native executable ◄──┤
                                                                     └> emcc: WebAssembly ◄────┘
```

- **SIR** (Skein intermediate representation) is method-neutral and has no
  mesh (§16).
- The **discrete IR** is per method: weak forms for CG and DG, cell
  balances for FV, point equations for meshless methods (§17).
- The generated Goose program calls the Skein runtime, a Goose package
  that holds the mesher, the solvers and the output writers (§18).

---

## 2. Lexical structure

### 2.1 Source text

A source file is UTF-8. The file extension is `.skn`. Identifiers use ASCII
only. Comments and strings may hold any UTF-8 text.

### 2.2 Lines and statements

A newline ends a statement. Skein has no semicolons. A statement continues
on the next line when the line ends with one of these tokens:

- a binary operator, `,`, `:` or `=`;
- an open bracket `(`, `[`, `{` or `<` whose close bracket is not on the
  line.

Reason: semicolons add one token to each statement and carry no
information for the parser.

### 2.3 Comments

`//` starts a comment that ends at the end of the line. Skein has no block
comments. `///` starts a documentation comment that attaches to the next
declaration. The language server shows it.

### 2.4 Tokens

```
ident      = letter { letter | digit | "_" }
letter     = "A".."Z" | "a".."z" | "_"
number     = digits [ "." digits ] [ ("e" | "E") [ "+" | "-" ] digits ]
digits     = digit { digit }
string     = '"' { any char except '"' and newline | '\"' | '\\' } '"'
percent    = number "%"
```

A number has no sign. A minus sign before a number is the unary operator,
except before a `degC` literal (§6.3).

### 2.5 Keywords

```
allow       assume      at          bc          birth       contact
crack_tip   death       eq          every       expect      false
fn          in          interface   layers      let         material
mesh        method      model       on          output      param
physics     plastic     points      refine      release     report
scales      solve       solves      state       step        study
symmetry    table       template    true        use         where
within
```

A keyword is not an identifier, except as the name before the `:` of a
setting, a named argument or a record entry (§3): `mesh: m`,
`step: "load"`. The names of types, units, built-in
functions and assumptions are not keywords. But a declaration with the name
of a built-in function or constant is an error [`duplicate_name`] (§4).

### 2.6 Unit names

The lexer reads a **unit expression** (§6.3) in these places: after a
number literal, after a vector literal, after the `]` of a list,
after the keyword `in` in a table column, and after the `:` of a field in
the header unit system. After each other `in`, an
expression follows. In the unit places, an identifier is a unit name, not
a variable. Thus `L` in `let L: Length = 2 L` is a variable in the first
place and the litre in the second place.

---

## 3. Grammar

The grammar is in EBNF. `{ x }` is zero or more, `[ x ]` is optional,
`( a | b )` is a choice. `NL` is a statement end (§2.2). Terminals are in
quotes.

```ebnf
program     = header { item } ;
header      = "model" ident ":" space "," unitsys NL ;
space       = "1d" | "2d" | "3d" | "axi" ;
unitsys     = ident | "nondim" | "{" unitfield { "," unitfield } "}" ;
unitfield   = ident ":" unitexpr ;

item        = usedecl | paramdecl | letdecl | fndecl | tabledecl | material
            | meshdecl | pointsdecl | physics | method | bcdecl | contact
            | interface | symmetry | step | expect | allow | scales | study ;

usedecl     = "use" string NL ;
fndecl      = "fn" ident "(" fnparam { "," fnparam } [ "," ] ")" "->" type "=" expr NL ;
fnparam     = ident ":" type ;

paramdecl   = "param" ident ":" type "=" expr [ "in" expr ".." expr ] NL ;
letdecl     = "let" ident ":" type "=" expr NL ;
scales      = "scales" expr { "," expr } NL ;
study       = "study" ident ":" procedure block NL ;  (* entry = ident ":" ( list | range ) | setting | check *)
list        = "[" expr { "," expr } "]" [ unitexpr ] ;
range       = expr ".." expr ;
tabledecl   = "table" ident "{" column { "," column } "}" ":" NL
              { INDENT row NL } ;
column      = ident ":" type [ "in" unitexpr ] ;
row         = expr { "," expr } ;

material    = "material" ident "{" NL { prop NL } "}" NL ;
prop        = ident ":" type "=" expr ;

block       = "{" [ NL ] [ entry { ( "," | NL ) entry } [ NL ] ] "}" ;
setting     = name ":" expr ;
name        = ident | keyword ;   (* a keyword only before ":" *)

meshdecl    = "mesh" ident "on" region block NL ;    (* entry = setting | refine | layers | cracktip *)
pointsdecl  = "points" ident "on" region block NL ;  (* entry = setting | refine *)
refine      = "refine" region ":" expr ;
layers      = "layers" region ":" expr ;
cracktip    = "crack_tip" region ":" expr ;

physics     = "physics" ident [ "on" region ] [ ":" ident [ "(" [ args ] ")" ] ]
              "solves" unknown { "," unknown } [ "{" NL { physstmt NL } "}" ] NL ;
unknown     = ident ":" type [ "=" expr ] ;
template    = "template" "physics" ident "(" [ tparam { "," tparam } ] ")"
              "solves" unknown { "," unknown } "{" NL { physstmt NL } "}" NL ;
tparam      = ident ":" type [ "=" expr ] ;      (* standard library only *)
physstmt    = matassign | assume | letdecl | eqn | statedecl | plastic ;
statedecl   = "state" ident ":" type "=" expr ;
plastic     = "plastic" block ;                      (* entry = setting *)
matassign   = "material" ident [ "on" region ] { "," ident "on" region } ;
assume      = "assume" assumption { "," assumption } ;
assumption  = ident [ "(" expr { "," expr } ")" ] ;
eqn         = "eq" expr "=" expr ;

method      = "method" ident ":" ident block NL ;    (* entry = setting *)

bcdecl      = "bc" region ":" bctarget "=" expr NL ;
bctarget    = expr ;
contact     = "contact" ident ":" region "," region block NL ;
interface   = "interface" ident ":" region "," region block NL ;
symmetry    = "symmetry" ident ":" procedure block NL ;  (* entry = setting *)

step        = "step" ident ":" procedure "{" NL { stepstmt NL } "}" NL ;
procedure   = ident [ "(" [ args ] ")" ] ;
stepstmt    = bcdecl | "release" region
            | "birth" region { "," region } [ "where" expr ]
            | "death" region { "," region } [ "where" expr ]
            | "solve" block
            | "output" expr { "," expr } [ "every" expr ]
            | "report" [ ident "=" ] expr { "," [ ident "=" ] expr }
            | check ;

expect      = check NL ;
check       = "expect" ident ":" addexpr "==" addexpr "within" expr ;  (* in a study, list | addexpr and list | expr *)
allow       = "allow" ident ":" string NL ;

region      = expr ;

type        = qtype | "Vector" "<" stype ">" | "Tensor2" "<" stype ">"
            | "Tensor4" "<" stype ">" | numsys | ident ;
stype       = qtype | numsys ;
numsys      = "Complex" "<" qtype ">" ;
qtype       = dimterm { ("*" | "/") dimterm } ;
dimterm     = ident [ "^" ["-"] number ] | "(" qtype ")" ;

expr        = orexpr ;
orexpr      = andexpr { ("||" | "|") andexpr } ;
andexpr     = cmpexpr { ("&&" | "&") cmpexpr } ;
cmpexpr     = addexpr [ ("==" | "!=" | "<" | "<=" | ">" | ">=") addexpr ] ;
addexpr     = mulexpr { ("+" | "-") mulexpr } ;
mulexpr     = unary { ("*" | "/") unary } ;
unary       = [ "-" | "!" ] power ;
power       = postfix [ "^" unary ] ;
postfix     = primary { "." ident | "(" [ args ] ")" | "at" primary } ;
primary     = number [ unitexpr ] | percent | string | "true" | "false"
            | ident | "(" expr ")" | "(" expr "," expr { "," expr } ")" [ unitexpr ]
            | "{" namedargs "}" ;
args        = arg { "," arg } [ "," ] ;
arg         = [ name ":" ] ( expr | list ) ;
namedargs   = name ":" expr { "," name ":" expr } ;

unitexpr    = unitfactor { ("*" | "/") unitfactor } ;
unitfactor  = unitname [ "^" ["-"] number ] | "(" unitexpr ")" ;
```

Notes on the grammar:

- A `block` holds entries on separate lines, or on one line with commas.
  `method solid: cg { mesh: m, order: 2 }` and the form with one entry on
  each line are the same.
- A `table` body ends at the first line that is not indented. This is the
  only place where indentation has a meaning. The rows are TOON rows
  ([TOON](https://github.com/toon-format/toon)): values that a comma
  separates. Skein omits the TOON row count `[N]`. Reason: the count changes
  with each new row, and the checker finds a short row without it.
- The header is one line: `model cantilever: 3d, mm_N_MPa_s`. `1d`, `2d`,
  `3d` and `axi` are tokens only in the header. Elsewhere `3d` is a syntax
  error.
- `axi` is the axisymmetric space (§11.3). It is not an assumption,
  because it changes the coordinates of every geometry function.
- The grammar admits more than the type rules accept. For example, a
  `region` is any `expr`. §5 and §8 give the rules.
- `|` and `&` are geometry union and intersection. `||` and `&&` are the
  Boolean operators. The two pairs have the same precedence.
- An argument may be a list, for example `harmonic(freq: [f0])`. A
  function takes no list literal `[...]` [`type_mismatch`]. A table
  column (§5.6) stays a valid function argument, for example
  `interp(t.col1, t.col2, x)`.
- An argument list and a study list may end with a comma. `skein fmt`
  writes a list that does not fit on one line with one item on each line
  and a comma after each item (§3.2).

### 3.1 Order of items

The header comes first. `use` lines come next. The items
after them may come in any order, but a name must be declared before its first
use [`use_before_decl`]. One exception: a material property may name an
unknown of a later physics that assigns the material (§10.2). Steps run in
the order of the file.

Reason: declaration before use lets the checker resolve names in one pass,
and it puts each error at the line of the use.

### 3.2 Layout for small diffs

A change to a model must give a small diff that a person can read. The
format rules of [command_line.md §3.1](command_line.md#31-format-rules) give
this layout: one statement on each line, no column alignment, one list item
on each line, and no reflow of a line that did not change.

Reason: LLMs write changes of thousands of lines that a person cannot
review. With these rules, a change to one load, one step or one study
parameter gives a diff of a few lines. The [LLM test](llm_test.md) measures
it.

---

## 4. Names and scope

A model has **one namespace**. These declarations add a name to it:
`param`, `let`, `fn`, `table`, `material`, `mesh`, `points`, `physics`, the
unknowns in `solves`, the `let` statements inside a `physics` block or a
physics template, `symmetry`, `step` and `expect`. The built-in functions
and constants (Appendix A), the standard library (§4.2) and the items of
each `use` file (§4.1) are in the same namespace, with no prefix. Two
declarations with one name are an error [`duplicate_name`], also when one
of them comes from a `use` file or from the standard library.

Reason: one namespace lets a step write `sigma*n` for the stress that a
physics block declares, with no prefix. A prefix such as `std.` or
`solid.` costs tokens on each use. The error on a clash keeps each name
unambiguous. This is the one deviation from the Power of Ten (rule 6,
smallest scope). The [tradeoffs](tradeoffs.md#6-reliability-the-power-of-ten)
give the reason.

These names have a smaller scope:

| Name | Visible in |
|---|---|
| A material property, for example `E` | the physics blocks that assign the material, the template arguments of those physics, and the `output`, `report` and `bc` expressions that use those physics |
| A column of a table, `t.col` | everywhere after the table |
| A field of a `Results` value, `r.u` | everywhere after the `let` |
| `x`, `y`, `z`, `t` (coordinates and time) | every expression that the solver evaluates at a point |
| `n` (outward unit normal) | `bc` expressions on a boundary, and the argument `f` of each `integral(f, r)` over a boundary region `r` in `output`, `report` and `expect` |
| `I` (identity tensor) | everywhere |

When two materials on the regions of one physics give a property, the
property must have the same type in both [`material_prop_type`]. When a
physics uses a property, each material that the physics assigns must give
it [`material_prop_missing`].

### 4.1 `use` (phase 2)

```skein
use "materials.skn"
```

A `use` line reads a library file. A library file holds only `material`,
`param`, `let` and `fn` items, with no header [`use_item`]. Its names enter
the namespace of the model with no prefix. A library file cannot `use`
another file. The path is relative to the model file.

Reason: a team keeps one file of materials. A copy in each model drifts.
Physics and steps stay in the model, so a reviewer reads one file for the
analysis. One level of `use` keeps the set of files of a model small and
free of cycles.

### 4.2 The standard library

The standard library holds physics templates and functions, written in
Skein source that a user can read (`skein help std`). A template is a
physics block with parameters:

```skein
template physics linear_elastic(b: Vector<Force/Volume> = 0) solves u: Vector<Length> {
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) + b = 0
}
```

A model uses it by name after the colon of a `physics` item:

```skein
physics solid on beam: linear_elastic solves u: Vector<Length> {
    material steel
}
```

The model names the unknowns, and their types must equal the types in the
template [`type_mismatch`]. An unknown may also be `Complex` of the type in
the template, for example `Vector<Complex<Length>>` for `linear_elastic`:
the template then solves its equations in complex numbers, for a
`harmonic` step (§14.2). The `let` names of the template, here `sigma`,
enter the namespace of the model. The block of the model may add `material`,
`assume` and `let` statements. It may not add `eq` [`template_eq`]. Only
the standard library declares templates in this version.

| Template | Unknowns |
|---|---|
| `linear_elastic(b)` | `u` |
| `heat(q_v)` | `T` |
| `stokes(f)`, `navier_stokes(f)` | `v`, `p` |
| `reissner_mindlin(layup, mirror, axis, top, q)` (shell), `timoshenko(A, I_y, I_z, J, y_axis, offset, q)` (beam), §4.3 | `w`, `theta` |
| `k_omega_sst(f)` (turbulence) | `v`, `p`, `k`, `w` (the specific dissipation rate ω: `omega` is the built-in of §14.2) |

The [roadmap](roadmap.md) gives the phase of each template.

Reason: nobody writes a shell or a k-ω SST model from memory, and a
template cuts the size of the common models. The source stays strong form,
so a user can read it and copy it into a model to change it.

### 4.3 Shells, beams and laminates

A shell is a physics on a `Surface`, and a beam is a physics on a `Curve`
(§8.3). Both solve a displacement `w: Vector<Length>` and a rotation vector
`theta: Vector<Angle>`. A physics that uses `reissner_mindlin` on a region
that is not a `Surface`, or `timoshenko` on a region that is not a `Curve`,
is an error [`geometry_dim`]. These are the two templates:

```skein
template physics reissner_mindlin(layup: Angle = 0 deg, mirror: Bool = false, axis: Vector<Real> = (1, 0, 0), top: Vector<Real>, q: Vector<Pressure> = 0) solves w: Vector<Length>, theta: Vector<Angle> {
    assume small_strain
    let eps_m: Tensor2<Real> = tangent(strain(w), top)
    let kappa: Tensor2<Real/Length> = tangent(sym(grad(cross(theta, top))), top)
    let gamma: Vector<Real> = tangent(grad(dot(w, top)) + cross(theta, top), top)
    let Q_ply: Tensor4<Pressure> = ply_stiffness(E1, E2, G12, nu12)
    let A_lam: Tensor4<Force/Length> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 0)
    let B_lam: Tensor4<Force> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 1)
    let D_lam: Tensor4<Force*Length> = laminate(Q_ply, layup, mirror, t_ply, h, axis, top, 2)
    let V_s: Vector<Force/Length> = 5/6*G12*h*gamma
    let N: Tensor2<Force/Length> = ddot(A_lam, eps_m) + ddot(B_lam, kappa) + outer(top, V_s)
    let M: Tensor2<Force> = ddot(B_lam, eps_m) + ddot(D_lam, kappa)
    eq div(N) + q = 0
    eq div(M) - V_s = 0
}

template physics timoshenko(A: Area, I_y: Length^4, I_z: Length^4, J: Length^4, y_axis: Vector<Real>, offset: Vector<Length> = 0, q: Vector<Force/Length> = 0) solves w: Vector<Length>, theta: Vector<Angle> {
    assume small_strain
    let e_yb: Vector<Real> = y_axis - dot(y_axis, tau)*tau
    let e_zb: Vector<Real> = cross(tau, e_yb)
    let G_b: Pressure = E/(2*(1 + nu))
    let gam_b: Vector<Real> = along(w + cross(theta, offset)) + cross(tau, theta)
    let kap_b: Vector<Real/Length> = along(theta)
    let Nx: Force = E*A*dot(tau, gam_b)
    let F_b: Vector<Force> = Nx*tau + 5/6*G_b*A*(dot(e_yb, gam_b)*e_yb + dot(e_zb, gam_b)*e_zb)
    let M_b: Vector<Force*Length> = G_b*J*dot(tau, kap_b)*tau + E*I_y*dot(e_yb, kap_b)*e_yb + E*I_z*dot(e_zb, kap_b)*e_zb
    eq along(F_b) + q = 0
    eq along(M_b) + cross(tau, F_b) = 0
}
```

**The shell.** The point at the distance z from the mid-surface, along
`top`, moves by `w + z*cross(theta, top)`. The material of the shell gives
the ply values `E1`, `E2`, `G12`, `nu12`, the ply thickness `t_ply` and
the laminate thickness `h`. An isotropic shell has `E1 = E2 = E`,
`G12 = E/(2*(1 + nu))`, `nu12 = nu` and `t_ply = h`.

| Argument | Meaning |
|---|---|
| `layup` | the ply angles from the `top` face, a table column (§5.6) or one angle |
| `mirror` | `true` adds the plies again in the reverse sequence: a symmetric laminate |
| `axis` | the 0 deg direction. Its projection on the surface is the axis 1 of the section. |
| `top` | the unit normal of the top face. A ply angle turns from `axis` by the right-hand rule about `top`. |
| `q` | the load per area, for example `p*n` for a pressure |

The plies of `layup` repeat to fill `h`, or to fill `h/2` with `mirror`.
Then `mirror` adds them again in the reverse sequence. The template gives `N`,
the force per length (the membrane force and the transverse shear force
`V_s` along `top`), the moment `M` and the membrane strain `eps_m`. The
membrane stress of a laminate is `N/h`. Because `div(N)` is in an `eq`,
`bc e: N*n = f` is an edge load per length (§12.1). The transverse shear
factor is 5/6, with G13 = G23 = G12.

**The beam.** `tau` is the unit tangent of the curve. The axis of the
section is at `offset` from the line of the physics: a rigid link. The
section y axis is `y_axis` made normal to `tau`, and z is `cross(tau, e_y)`.
`I_y` is `∫ z² dA` and `I_z` is `∫ y² dA`. The material gives `E` and
`nu`. The template gives the axial force `Nx` (positive in tension), the
force `F_b` and the moment `M_b`.

A `modal` step (§14.2) of beams takes the mass from the density `rho` of
each beam material: `rho*A` for the translation, `rho*(I_y + I_z)` for the
twist and `rho*I_y`, `rho*I_z` for the rotary inertia, on the axis of the
section. A `buckling` step of beams takes the geometric stiffness of the
axial force `Nx` of the `static` step just before it.
Reason: a stiffened panel needs the modes and the buckling load of its
stiffeners (implementation.md §10). The template has no time derivative,
so the density is the only input.

The built-in functions of the two templates:

| Function | Meaning |
|---|---|
| `tangent(f, top)` | `P*f` for a vector, `P*f*P` for a `Tensor2`, with `P = I - outer(top, top)` |
| `ply_stiffness(E1, E2, G12, nu12)` | the plane stress stiffness of a ply in its fibre axes, a `Tensor4` |
| `laminate(Q, layup, mirror, t_ply, h, axis, top, k)` | `∫ Q(z) z^k dz` over the plies, with each ply turned to its angle. `k` is 0, 1 or 2 [`type_mismatch`]: A, B and D of classical laminate theory. |
| `along(f)` | the derivative of `f` along the curve, `df/ds` |

Reason: a fuselage, a wing and a ship hull are shells with offset
stiffeners. A solid mesh of a 4.6 mm skin on a 3 m radius needs millions
of cells. The ABD matrix of classical laminate theory is the standard
section of a composite shell (Jones 1999, §4.2). The solver uses the MITC4
shell (Dvorkin and Bathe 1984) for `order: 1` and the MITC9 shell (Bucalem
and Bathe 1993) for `order: 2`, so the shell does not lock in shear.

---

## 5. Types

Skein has no type inference. Each `let`, `param`, `fn` parameter and
result, unknown in `solves`, material property and table column states its
type. The checker compares the declared type with
the type of the expression and reports [`type_mismatch`] or
[`unit_mismatch`] when they differ.

### 5.1 Scalar types

| Type | Values |
|---|---|
| a dimension, for example `Length` | a real number with that dimension (§6) |
| `Real` | a dimensionless real number |
| `Angle` | a dimensionless real number. Literals take `rad` or `deg`. |
| `Int` | an integer, for counts such as the number of modes |
| `Bool` | `true` or `false` |

A quantity can also be complex. `Complex<D>` holds `a + b*im`, with `a`
and `b` of dimension `D`. A tensor takes it in its argument:
`Vector<Complex<Velocity>>`. The built-in constant `im` is the imaginary
unit. Reason: the `harmonic` procedure (§14.2) of the acoustic and
electromagnetic models of femahoy needs complex unknowns. Skein has no
other number system, because no model needs one.

### 5.2 Tensor types

| Type | Rank | Components |
|---|---|---|
| `D` (a dimension) | 0 | 1 |
| `Vector<D>` | 1 | `d` |
| `Tensor2<D>` | 2 | `d × d` |
| `Tensor4<D>` | 4 | `d × d × d × d` |

`d` is the space dimension from the header (a `2d` model gives 2, an `axi` model
gives 2 for vectors and 3 for the hoop components of tensors, §11.3). All
components of a tensor have the one dimension `D`. Components are `.x`,
`.y`, `.z` for vectors and `.xx`, `.xy`, … for second-order tensors.

### 5.3 Geometry types

| Type | Topological dimension | Examples |
|---|---|---|
| `Solid` | 3 | a 3D body |
| `Surface` | 2 | a 2D domain, a face of a solid, a shell |
| `Curve` | 1 | a 2D boundary, an edge, a beam |
| `Point` | 0 | a vertex, a point load, a pressure reference |

The checker computes the topological dimension of each geometry expression.
A declaration with the wrong geometry type is an error [`geometry_dim`].

### 5.4 Model types

| Type | Made by |
|---|---|
| `Mesh` | a `mesh` block (§9) |
| `Points` | a `points` block (§9) |
| `Results` | `results(file, step:)` (§15.1) |

### 5.5 Fields

An unknown in `solves` (§10.1) has a tensor type and the region of its
physics. An unknown of a physics with no region is a global unknown, for
example the state of an ODE (§10.5).

### 5.6 Tables

A table has typed columns and TOON rows. `in unit` after a column type gives
the unit of each bare number in that column. A cell may also carry its own
unit. Each row must have one value for each column [`table_row_width`].

`t.col` is the column as a list. A list may stand where an expression of
the column type is expected in `expect`, `output`, `report` and in the
point of an `at`. The statement then applies to each row (§15.3).

`sum(expr)` is the sum of the scalar or tensor expression `expr` over the
rows of a table. It has the type of `expr`, and it is one value, not a
list. `expr` must use a column of a table [`sum_no_column`], and
all its columns must come from one table [`sum_mixed_tables`].

```skein
let sxy: Real = sum(obs.x*obs.y)
```

Reason: a least squares fit states its normal equations with sums over
the rows of the data (femahoy gap G5).

### 5.7 Assignability

A value of type `A` may initialize a declaration of type `B` only when `A`
and `B` are the same type. Skein has no implicit conversion. These are the
only exceptions:

1. The literal `0` (with no unit) has every dimension and every tensor rank.
   It means zero of the expected type.
2. A bare number with no `.` and no exponent has type `Int`. Each other
   bare number has type `Real`. `Real` and `Angle` accept `Int`. `Angle`
   accepts `Real`.
3. `Int` converts to `Real` in arithmetic.
4. `Complex<D>` accepts a value of `D`.

Reason for exception 1: `u = 0` is the most common boundary condition, and
`u = (0, 0, 0) mm` adds tokens and gives no safety.

---

## 6. Dimensions and units

### 6.1 Dimensions

A dimension is a vector of seven exponents over the SI base dimensions:
mass M, length L, time T, temperature Θ, current I, amount N and luminous
intensity J. Each exponent is a multiple of 1/2. Reason: a stress intensity
factor has the dimension `Pressure*Length^0.5` (femahoy gap G14). The
checker does arithmetic on these vectors:

- `a*b` adds the exponents, `a/b` subtracts them, `a^k` multiplies them by
  `k`. `k` must be a constant, and each exponent of the result must be a
  multiple of 1/2 [`unit_fractional`]. When `a` is `Real`, `k` may be any
  constant `Real`. A varying `k` is an error (§7.6) [`not_constant`].
  A function parameter is constant when each call gives it a constant.
- After `^` in a type or a unit, the exponent is a number that is a
  multiple of 0.5: `Pressure*Length^0.5`, `MPa*mm^0.5`, `Length^-1.5`.
- `a + b`, `a - b`, `a == b` and the other comparisons need equal
  dimensions [`unit_mismatch`].
- `sqrt(a)` halves the exponents. Each exponent of `a` must be a whole
  number [`unit_fractional`].
- `exp`, `ln`, `sin`, `cos`, `tan`, `tanh` take `Real` or `Angle`
  [`unit_transcendental`]. `atan` takes `Real`. `atan2(y, x)` takes two
  values of one dimension [`unit_mismatch`]. Both give an `Angle`.
- `grad`, `div` and `curl` divide by L. `ddt` divides by T. `ddt2` divides
  by T².

### 6.2 Named dimension types

These names are built in. For any other dimension, write a dimension
expression, for example `Force/Length` or `Energy/Volume`. Skein has no type
aliases (§6.6).

| Name | Dimension | | Name | Dimension |
|---|---|---|---|---|
| `Length` | L | | `Velocity` | L T⁻¹ |
| `Area` | L² | | `Acceleration` | L T⁻² |
| `Volume` | L³ | | `Frequency` | T⁻¹ |
| `Mass` | M | | `Density` | M L⁻³ |
| `Time` | T | | `Pressure` | M L⁻¹ T⁻² |
| `Temperature` | Θ | | `Stress` | M L⁻¹ T⁻² |
| `Force` | M L T⁻² | | `Energy` | M L² T⁻² |
| `Power` | M L² T⁻³ | | `DynViscosity` | M L⁻¹ T⁻¹ |
| `KinViscosity` | L² T⁻¹ | | `Conductivity` | M L T⁻³ Θ⁻¹ |
| `SpecificHeat` | L² T⁻² Θ⁻¹ | | `Expansion` | Θ⁻¹ |
| `HeatFlux` | M T⁻³ | | `HeatTransfer` | M T⁻³ Θ⁻¹ |
| `Current` | I | | `Amount` | N |

`Pressure` and `Stress` are the same type. The checker prints the name that
the declaration used.

### 6.3 Units

A unit is a dimension and a scale factor to SI. A unit expression combines
unit names with `*`, `/`, `^` and parentheses: `kg/m^3`, `W/(m^2*K)`.

Base and derived units: `m`, `g`, `s`, `K`, `A`, `mol`, `cd`, `N`, `Pa`,
`J`, `W`, `Hz`, `C`, `V`, `rad`, `deg`, `L` (litre), `t` (tonne), `min`,
`h` (hour), `bar`, `degC`.

The SI prefixes `n`, `u`, `m`, `c`, `k`, `M`, `G` apply to `m`, `g`, `s`,
`N`, `Pa`, `J`, `W`, `Hz`, `A`, `V` and `L`. A unit name is the longest
match: `mm` is millimetre and `ms` is millisecond. `min` is minute. `Mt` and
other prefixed tonnes are not units.

A unit expression after a literal ends at the first token that cannot
continue it. After `*` or `/`, it continues only when the next token is a
unit name, or `(` and then a unit name. Thus `2 m/s*(1 - y/R)` is
`2 m/s` times `(1 - y/R)`. A variable with the name of a unit, for example
`L`, directly after `*` in a unit expression is read as the unit. Put the
literal in parentheses to avoid this: `(2 kN)*L`.

`degC` is an affine unit. It is valid only alone, after a literal of
dimension `Temperature`: `20 degC` is 293.15 K. `degC` in a compound unit
is an error [`unit_affine`]. Write temperature differences and
coefficients in `K`. A unary `-` directly before a `degC` literal is part of
the literal: `-5 degC` is 268.15 K, not −278.15 K.

### 6.4 Literals

- `210 GPa` is a scalar. Its dimension comes from the unit.
- `(40, 0, 20) mm` is a vector. The unit after the parenthesis applies to
  each bare number in it.
- `(L, W, H)` is a vector of expressions. All components must have the same
  dimension. A `0` component takes that dimension.
- A bare number with no `.` and no exponent is `Int`. Each other bare
  number is `Real` (§5.7).

The value of a literal is stored in SI. The checker compares dimensions,
never units, so `1 m + 1 mm` is valid and is 1.001 m.

### 6.5 Unit systems

The `unitsys` part of the header names the unit system of the model:

| Name | Length | Mass | Time | Force | Stress | Temperature |
|---|---|---|---|---|---|---|
| `SI` | m | kg | s | N | Pa | K |
| `mm_N_MPa_s` | mm | t | s | N | MPa | K |
| `mm_kg_ms` | mm | kg | ms | kN | GPa | K |

A user system gives the base units in the header: `model m: 3d, { length:
mm, mass: t, time: s, temperature: K }`.

The unit system does not change the meaning of the source, because each
literal carries its unit. It sets the units of the output files and of the
numbers that the tools print. Reason: a user reads results in the units of
the drawing.

The `scales` line (§21.2), not the unit system, sets the size of the
numbers in the solver.

The header `model m: 2d, nondim` declares a model with no dimensions. Each
quantity has type `Real`, and a literal with a unit is an error
[`nondim_unit`]. Use it for
a PDE that has no physical units, for example a model problem from a
textbook in dimensionless form. The dimension checks then give no
protection, so the `check` report says `units: nondim` at its top.
Reason: the femahoy cases shape-hole-biaxial and fit-nist-pontius have no
units.

### 6.6 What Skein leaves out

- **No type aliases.** Reason: an alias puts one fact in two places, and an
  LLM must then find the alias to know the dimension.
- **No unit-generic functions.** Each parameter of a function has one
  fixed type (§7.4).

---

## 7. Expressions and operators

### 7.1 Precedence

From low to high: `||` `|`, `&&` `&`, comparisons, `+` `-`, `*` `/`, unary
`-` `!`, `^`, postfix (`.`, call, `at`). `^` is right-associative. All other
binary operators are left-associative. A comparison does not chain.

### 7.2 Arithmetic on tensors

| Expression | Rule |
|---|---|
| `a + b`, `a - b` | same rank, same dimension |
| `s*A`, `A*s`, `A/s` | `s` scalar, `A` any rank |
| `A*v` | `Tensor2` times `Vector` gives `Vector`, `(A*v)_i = A_ij v_j` |
| `A*B` | `Tensor2` times `Tensor2` gives `Tensor2` |
| `C*A` | `Tensor4` times `Tensor2` gives `Tensor2`, `(C*A)_ij = C_ijkl A_kl` |
| `v*w` | an error [`ambiguous_product`]. Use `dot`, `outer` or `cross`. |

### 7.3 Differential operators

| Operator | Argument | Result |
|---|---|---|
| `grad(f)` | rank r | rank r + 1. `grad(v)_ij = ∂v_i/∂x_j` |
| `div(F)` | rank r ≥ 1 | rank r − 1. `div(A)_i = ∂A_ij/∂x_j` |
| `curl(v)` | `Vector`, a `3d` model | `Vector` |
| `lap(f)` | rank r | rank r. Equal to `div(grad(f))` |
| `ddt(f)` | any | same rank. First time derivative |
| `ddt2(f)` | any | same rank. Second time derivative |

With these definitions `grad(v)*v` is the convective term `(v·∇)v`.

`ddt` and `ddt2` are zero in a `static` step. The `transient` procedure
replaces them with the time scheme of the step (§14.2).

### 7.4 Functions

```skein
fn lame_mu(E: Pressure, nu: Real) -> Pressure = E/(2*(1 + nu))
```

A function is one expression of its parameters and of the names declared
above it. It has no statements, no loops and no recursion: a function may
call only the functions declared above it [`use_before_decl`]. The result
type is required, and the checker compares it with the type of the
expression.

Reason: a law that a model uses in many places is written once, which cuts
the size of the model. One expression keeps rule 4 of §1 and rules 1 and 2
of the Power of Ten by construction.

### 7.5 `at`

`f at p` evaluates the field expression `f` at the point `p`. `p` is a
`Point` or a vector of `Length`. If `p` is outside the region of `f`, the
run fails with [`probe_outside`].

### 7.6 Coordinates and time

`x`, `y` and `z` are the coordinates of the point where the solver
evaluates the expression. `t` is the time. A `let` whose expression uses
them is a **varying** value. A setting that must be constant, for example a
mesh size or the end time of a step, must not use a varying value
[`not_constant`].

---

## 8. Geometry and regions

### 8.1 Distance functions

Each geometry value is a signed distance function φ: φ < 0 inside, φ = 0 on
the boundary, φ > 0 outside. The checker keeps the CSG tree, not a sampled
function, so that the mesher can find sharp edges (§9.2).

| Function | Type | Notes |
|---|---|---|
| `box(p0, p1)` | `Solid` | corners `p0` and `p1` |
| `cylinder(c, axis, r, h)` | `Solid` | base centre, axis direction, radius, height |
| `sphere(c, r)` | `Solid` | |
| `rect(p0, p1)` | `Surface` | 2D only |
| `circle(c, r)` | `Surface` | 2D only |
| `polygon(p1, p2, …)` | `Surface` | 2D only, counter-clockwise or clockwise |
| `fill(c)` | `Surface` | 2D only. The region inside the `Curve` `c` and the line from its last point to its first point |
| `extrude(s, h)` | `Solid` | `s` is a 2D `Surface` in the x-y plane |
| `revolve(s, angle)` | `Solid` | about the y axis |
| `sweep(s, path)` | `Solid` | `path` is a `Curve` |
| `line(p0, p1)`, `arc(c, r, a0, a1)`, `spline(p…)` | `Curve` | `spline` passes through the points in order |
| `nurbs(d, knots(k…), p1, w1, …, pn, wn)` | `Curve` | 2D only. The NURBS curve of degree `d` with the control points `pi` and their weights `wi`. `knots(k…)` is optional |
| `point(p)` | `Point` | |
| `shell(a, t)` | same as `a` | a wall of constant thickness `t` inside the surface of `a`: `max(φa, −φa − t)` |
| `loft(s1, s2, h)` | `Solid` | from 2D section `s1` at z = 0 to `s2` at z = `h` |
| `pattern(a, v, n)`, `pattern_polar(a, axis, n)` | same as `a` | `n` copies of `a`, along `v` or about `axis` through the origin. Copy k (0 to n − 1) is `move(a, k*v)`, or `a` turned by k·360°/n. |

`copy(r)` is a `Real` that varies in space: at a point, the index k of the
copy of the pattern `r` that is nearest to the point. `r` must come from a
`pattern` or a `pattern_polar` [`type_mismatch`]. Thus `x0 + copy(r)*dx` is
the x of the copy at the point, for example the centre of each window for a
load on the window edges:

```skein
let windows: Solid = pattern(window_0, (f_pitch, 0 m, 0 m), 8)
let window_edge: Curve = boundary(skin, windows)
let dx_w: Length = x - (x_w0 + copy(windows)*f_pitch)
```

Reason: gap G35. The load on the edge of a window depends on the distance
from the centre of that window, and the copies of a pattern have no other
name.

`spline(p1, …, pn)` is the centripetal Catmull-Rom spline through the
points (Yuksel et al. 2011). It has no cusp and no loop between two points.
At each end, it continues to a point that is the reflection of the second
point, for example `2*p1 - p2` before `p1`. A point equal to the point
before it adds nothing. `fill(spline(p1, …, pn))` closes the spline with a
line from `pn` to `p1`. Thus the ends make a corner, for example the sharp
trailing edge of an airfoil when `pn` equals `p1`, or the straight edges of
a quarter hole on 2 lines of symmetry:

```skein
let hole: Surface = fill(spline(rim(0), rim(1), rim(2), rim(3), rim(4))) | polygon((0, 0), rim(0), rim(4))
```

Reason: femahoy gap G24. The airfoils of naca0012-ladson and naca23012-flap
and the superellipse hole of shape-hole-biaxial are splines through
points. A polygon through the points has a different boundary.

`nurbs(d, knots(k1, …, km), p1, w1, …, pn, wn)` is the NURBS curve
(Piegl and Tiller 1997) of degree `d` (1, 2 or 3) with n ≥ d + 1 control
points `pi`, each with its weight `wi` > 0 (a `Real`). The m = n + d + 1
knots do not fall, the first and the last are each d + 1 times, and an
inner knot is d times or fewer [`type_mismatch`]. Without `knots`, the
knots are 0 (d + 1 times), 1, 2, …, n − d (d + 1 times). `fill(nurbs(…))` closes the curve with a line
from its end to its start, as for `spline`. A rational curve of degree 2
is a conic, so a circle is exact:

```skein
let c: Real = 0.7071067811865476
let bore: Surface = fill(nurbs(2, knots(0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4),
    (1, 0) mm, 1, (1, 1) mm, c, (0, 1) mm, 1, (-1, 1) mm, c, (-1, 0) mm, 1,
    (-1, -1) mm, c, (0, -1) mm, 1, (1, -1) mm, c, (1, 0) mm, 1))
```

The mesher cuts the curve at each knot into rational Bézier spans of
degree 3. Reason: a CAD curve is a NURBS curve, and a spline through
points of a circle is not the circle.

Operators and modifiers:

| Expression | Meaning | Distance function |
|---|---|---|
| `a | b` | union | `min(φa, φb)` |
| `a & b` | intersection | `max(φa, φb)` |
| `a - b` | difference | `max(φa, −φb)` |
| `fillet(a, r)` | round the convex and concave edges of `a` | smooth minimum with radius `r` |
| `move(a, v)`, `rotate(a, axis, angle)`, `scale(a, k)` | rigid motions and uniform scale | |

The `min` and `max` of distance functions give a bound on the distance, not
the exact distance. The mesher uses the CSG tree to correct this near the
surface (§9.2).

### 8.2 Regions

A **region** is a geometry value that a `physics`, `mesh`,
`points`, `bc` or `birth` statement names. The geometry types are the region
types.

`boundary(g, s)` selects a part of the boundary of `g`:

- `s` is a Boolean expression of `x`, `y`, `z`: the part where `s` is true.
  An equality holds within a tolerance of 10⁻⁶ times the size of the
  bounding box of `g`.
- `s` is a geometry value: the part of the boundary of `g` that lies on the
  boundary of `s`. This selects the surface of a hole that `g - s` cut.
- `boundary(g)` with no `s` is the whole boundary.

The result has dimension `dim(g) − 1`. An empty selection is an error
[`empty_region`].

### 8.3 Lower-dimension entities

A `Surface` in a `3d` model, or a `Curve` in a `2d` model or a `3d` model, can be
a region of its own physics, for example a shell or a beam. Such a region
must be the boundary of a solid region, or a geometry value that the
`mesh` block embeds. The mesher makes the cells of the entity conform to
the cells of the region around it. A `Point` can be the region of a `bc`,
for example a point load or a pressure reference.

### 8.4 Region algebra

Regions combine with `|`, `&` and `-`. `a | b` has the higher dimension of
`a` and `b`, `a & b` has the lower dimension, and `a - b` has the dimension
of `a`. Thus `skin & door_box` is the part of the `Surface` `skin` in the
`Solid` `door_box`, a `Surface`. Reason: gap G32. A shell or a beam line
is cut by a solid, for example a pad-up band or a door cutout.

The checker computes containment by
the CSG tree: `a` is inside `b` if `a & b` equals `a` in the tree, or if a
sampled test on the mesh shows it. A field, a physics, a `bc` or a
`birth`/`death` statement must name a region inside the region of the unknown
[`region_outside`].

---

## 9. Meshes and point clouds

Skein does not import meshes. Reason: decision 8 ([decisions.md](decisions.md)).
The built-in mesher makes every mesh from the geometry.

### 9.1 `mesh`

```skein
mesh m on beam {
    size: 5 mm
    cell: tet
    grade: 1.3
    feature_angle: 30 deg
    refine hole: 1 mm
    layers wall: { first: 0.1 mm, ratio: 1.2, count: 8 }
}
```

| Setting | Type | Default | Meaning |
|---|---|---|---|
| `size` | `Length` | required | target edge length |
| `cell` | `tri`, `quad`, `tet`, `hex` | `tri` in 2D, `tet` in 3D | `quad` and `hex` on a union of `rect` or `box` regions, or on a quarter of a ring or a tube, are a structured grid. `quad` on another `Surface` has triangles only where the shape needs them. `hex` on another `Solid` is a swept or a mapped grid (§9.2). |
| `grade` | `Real` | 1.3 | largest ratio of the size of two neighbour cells |
| `feature_angle` | `Angle` | 30 deg | an edge with a dihedral angle smaller than 180 deg minus this angle is sharp |
| `refine r` | `Length` | none | size on region `r` and near it |
| `layers r` | record | none | boundary layer cells on `r` (CFD) |
| `crack_tip p` | `Length` | none | quarter-point cells of this size at the crack tip `p` |

The cell order comes from the `method` (§13), not from the mesh. The mesher
makes the nodes for that order on the true geometry.

`crack_tip p: L` makes a fan of triangles with edges of length `L` around
the `Point` `p` of a 2d model. On each edge that touches `p`, the mesher
moves the mid-side node to `L/4` from `p`. These quarter-point cells give
the 1/sqrt(r) strain of a linear elastic crack tip (Barsoum 1976). The
crack opening at `r = L/4` and `r = L` then gives the stress intensity
factor (Shih et al. 1976). The error [`crack_tip_form`] (checks.md)
rejects `crack_tip` in a 3d model or in a `points` block, a region that is
not a `Point`, a size that is not a `Length`, and a displacement with an
order other than 2.
Reason: femahoy gap G20. The case fracture-centre-crack reads K_I from
the nodes of quarter-point cells, and no other mesh setting moves a node.

### 9.2 Mesh generation

The algorithms come from published papers only. The implementation must
not copy code or text from Gmsh, Triangle, TetGen or CGAL, because their
licences are not compatible with Apache 2.0.

| Case | Algorithm |
|---|---|
| 2D triangles | Delaunay refinement (Chew 1989, Ruppert 1995, Shewchuk 2002) on a boundary that the mesher samples from φ, with the sharp corners as fixed vertices |
| 3D tetrahedra | isosurface stuffing (Labelle and Shewchuk 2007) on an octree graded by `size`, `refine` and `grade` |
| 3D sharp edges | protecting balls on the feature curves (Cheng, Dey and Ramos 2007), then the stuffing |
| 3D tetrahedra on a prism along z with a `cylinder` or an `extrude`, and with each `refine` target in its full z range | the triangles of the 2D case on its section, swept through layers of equal height, each prism cut into 3 tetrahedra along the diagonal of each side face from its point of smallest index (Dompierre et al. 1999) |
| smoothing | optimal Delaunay triangulation smoothing (Chen and Xu 2004, Alliez et al. 2005), with no move that makes a cell worse |
| `quad`, `hex` | a structured grid on a union of `box` or `rect` regions, or in r, θ (and z) on a quarter of a ring or a tube: `(circle - circle) & rect` or `(cylinder - cylinder) & box`, with the corner of the rect or the box on the axis |
| `quad` on another `Surface` | the triangles of the 2D case, then the merge of pairs of neighbour triangles into quads (Blossom-Quad, Remacle et al. 2012): a matching of the largest size (Edmonds 1965) on the pairs whose quad has an acceptable quality. A triangle stays where it has no such pair. |
| `hex` on a prism along z | a `Solid` of `box`, `cylinder` along z and `extrude` regions that each hold its z range: the quads of Blossom-Quad on its section at twice the size, each cut into 4 quads and each remaining triangle into 3, swept through layers of equal height |
| `hex` on another `Solid` | the grid of the planes of the `box` regions in the bounds of the region, with the cells whose centre is in the region. Each node on the boundary of the cells moves to the nearest point of the boundary of the region (Schneiders 1996). |
| boundary layers | prisms or quads extruded along the normal of the surface, then the core mesh |

On a union of boxes or rects, each face of each `box` or `rect` of the
model is a grid plane: the faces of the region, and the faces of the boxes
that `let` items join with `|`, `&` and `-`, for example a material layer.
Between two neighbour planes, the grid has the fewest equal cells that are
not longer than `size`. A cell is in the mesh if its centre is in the
region. Thus no cell crosses a material face, and a thin layer has one
cell through it.

The sharp edges come from the CSG tree: each intersection of the surfaces
of two primitives with a dihedral angle that `feature_angle` marks is a
feature curve.

Reason for the swept tetrahedra: femahoy gap G4. The stuffing cuts the
sharp edge between the end plane and the curved wall of a coax line off
by about one cell. That chamfer is in neither region, so the port does not
touch the wall, and the port voltage floats.

Reason for the `quad` and `hex` cells on a region that is not a union of
boxes: femahoy gap G36. Femahoy meshes a ring about a hole, a tube, a rod
and a sloped fillet with quad or hex cells.

Isosurface stuffing alone guarantees dihedral angles between 10.7° and
164.8°. The feature protection and the boundary layers remove this
guarantee. Thus the mesher always measures the result, and the pre-solve
check [`mesh_quality`] (checks.md) reports the worst cells.

### 9.3 `points`

```skein
points c on quarter {
    spacing: 2 mm
    refine hole: 0.4 mm
}
```

A `points` block makes a point cloud for a meshless method. The interior
points come from Poisson disk sampling (Bridson 2007), graded like a mesh.
The boundary points lie on φ = 0, with the normal from ∇φ. Each sharp
corner gets a point. The pre-solve check [`points_quality`] measures the
fill distance and the separation ratio.

The current lowering (`src/meshless/`) takes the points of a 2D model in
one of two ways. On a `rect`, the points are a lattice of the `spacing`,
and each point inside moves by up to 0.15 of the spacing, so the cloud is
not a grid. On the plate with a hole of sample 05, the points are the
nodes of the P1 triangles of a mesh with the `spacing` as its size and
the `refine` at the hole. In both, the triangles between the points give
the nodal domains of `rkpm`, the boundary edges and the local spacing.
Poisson disk sampling and `points_quality` do not run yet.

---

## 10. Fields, materials and physics

### 10.1 Unknown fields

```skein
physics thermal on part solves T: Temperature = 293 K {
```

The `solves` list of a physics declares its unknowns. Each unknown has a
tensor type (§5.2). Its region is the region of the physics. The optional
`= expr` is the initial value (§14.1). It may be varying. Each unknown has
exactly one physics, the one that declares it.

Reason: a separate `field` line repeated the name and the region of each
unknown, and a model could declare a field that no physics solved.

### 10.2 `material`

```skein
material steel {
    E: Pressure = 210 GPa
    nu: Real = 0.3
    rho: Density = 7850 kg/m^3
}
```

A material is a record of typed properties. A property may use the
properties above it in the same material. A property may be varying (§7.6),
for example a function of a temperature field. A property may name an
unknown of each physics that assigns the material, also when the physics
comes after the material (§3.1):

```skein
material steel {
    k: Conductivity = interp(k_steel.temp, k_steel.k_val, T)
}

physics thermal on fin solves T: Temperature = 300 K {
    material steel
    eq div(-k*grad(T)) = 0
}
```

Reason: a material comes before the physics that assigns it, and a
conductivity of `T` is a property of the material, not of the physics.

### 10.3 `physics`

```skein
physics solid on beam solves u: Vector<Length> {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}
```

- `on r` is the region of the equations and of the unknowns. A physics
  with no `on` has only global unknowns (§10.5).
- `: template(args)` takes the statements from a standard-library template
  (§4.2). The block is then optional.
- `solves u: T, p: T, …` declares the unknowns. The block must have one `eq` for each
  unknown [`eq_count`]. Equation `k` pairs with unknown `k`: the test
  function of equation `k` is a variation of unknown `k` (§17.1).
- `material` assigns materials to parts of the region. Each point of the
  region must have exactly one material [`material_cover`].
- `assume` names modelling assumptions (§11).
- `let` declares a quantity for the equations, for example a stress.
- `eq lhs = rhs` is an equation. Both sides must have the rank and the
  dimension of the paired unknown's equation, and the same as each other
  [`unit_mismatch`], [`rank_mismatch`].

### 10.4 Coupling

Two physics blocks couple when one uses a field that the other solves:

```skein
physics thermal on part solves T: Temperature = 293 K { … }
physics solid on part solves u: Vector<Length> {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) - alpha*(T - T0)*I)
    eq div(sigma) = 0
}
```

The `method` of each physics gives the scheme of the coupled solve: one
system (`couple: monolithic`, the default) or one physics after the other
in each step (`couple: staggered`).

### 10.5 ODEs and global unknowns

A physics without `on` solves only global unknowns. It is a system of
ODEs:

```skein
physics oscillator solves q: Length {
    eq m*ddt2(q) + c*ddt(q) + k*q = 0
}
```

A global unknown may appear in a physics with a region. Then the system
couples an ODE and a PDE, for example a lumped spring on a solid.

In a physics, `mean(f, r)` and `integral(f, r)` (§14.5) are valid only in
the `eq` of a global unknown. In a physics with a region, the unknown of an
`eq` with `mean` or `integral` is a global unknown, but only if the `eq` has no
`grad`, `div`, `curl` or `strain` outside `mean` and `integral`. Outside
`mean` and `integral`, an `eq` with `mean` or `integral` uses no field: no
`x`, `y`, `z`, `n` or `tau`, no `grad`, `div`, `curl` or `strain`, and no
unknown with a region that is not global [`global_eq`]. A `let` in the `eq`
follows the same rule. Thus a field PDE with a mean term is not valid, for
example `eq -div(k*grad(T)) = Q - mean(Q, b)`.

```skein
physics cell on b solves u: Vector<Length>, ez: Real {
    material fibre on f, matrix on b - f
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) + ez*outer(e3, e3))
    eq div(sigma) = 0
    eq mean(sigma.zz, b) = 100 MPa
}
```

Reason: the effective modulus of a periodic cell needs a mean strain that
a mean stress sets (femahoy gap G15).

### 10.6 Inelastic materials

A physics can hold **state variables**: values at each integration point
that the solver updates at the end of each converged increment.

```skein
state eps_p: Tensor2<Real> = 0
state alpha: Real = 0
let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u) - eps_p)
plastic {
    yield: von_mises(sigma) - (sy + Hk*alpha)
    strain: eps_p
    hardening: alpha
}
```

`plastic` states a yield function `f(σ, α) ≤ 0` with associative flow:
`ddt(eps_p) = λ ∂f/∂σ`, `ddt(alpha) = λ √(2/3) |∂f/∂σ|`, `λ ≥ 0`, `λ f = 0`.
The compiler makes the closest-point return map and its consistent tangent
(Simo and Hughes 1998) by symbolic differentiation of `f`. `yield` must have
the dimension of the stress [`unit_mismatch`].

Under `finite_strain`, with the equation `div(defgrad(u)*S) = 0`, the same
block states the multiplicative split `F = Fe*Fp`. `strain(u) - eps_p` is
then the logarithmic elastic strain `ln(Ve)`, with `Fe = Ve*Re`, and
`hooke` gives the Kirchhoff stress from it (a Hencky law). The `strain`
state holds `inv(Cp)`, with `Cp = transpose(Fp)*Fp`, and the return map is
the exponential map of Simo (1992), J2 on the Kirchhoff stress. `S` is
`inv(F)*tau*inv(transpose(F))`. At small strain the law is the additive law
above. Reason: the Taylor bar of femahoy has plastic strains up to 3.24
(gap G22). There the additive split of the Green-Lagrange strain is a
different law, and its stress grows without limit at large strain.

A hyperelastic law needs no state. Write the strain energy and derive the stress with `diff`:

```skein
assume finite_strain, isotropic
let C: Tensor2<Real> = transpose(defgrad(u))*defgrad(u)
let W: Energy/Volume = mu/2*(tr(C) - 3) - mu*ln(det(defgrad(u))) + lam/2*ln(det(defgrad(u)))^2
let S: Tensor2<Pressure> = 2*diff(W, C)
eq div(defgrad(u)*S) = 0
```

### 10.7 FSI interfaces

```skein
interface fsi: wet_fluid, wet_solid {
    mesh_motion: harmonic
}
```

An `interface` joins a boundary of a fluid physics and a boundary of a solid
physics. It adds the kinematic condition (the fluid velocity equals
`ddt(u)` of the solid) and the dynamic condition (the tractions balance).
The fluid equations move to an arbitrary Lagrangian-Eulerian frame, and the
setting `mesh_motion` gives the equation of the mesh displacement.

---

## 11. Modelling assumptions

An assumption changes the meaning of the equations in a defined way. The
compiler applies it before the method sees the equations. A user can read
each effect in this table. The [reference](reference.md#11-modelling-assumptions)
gives examples.

### 11.1 Assumption table

| Assumption | Effect | Conflicts with |
|---|---|---|
| `small_strain` | `strain(u)` is `sym(grad(u))`. The equations are on the reference shape. | `finite_strain` |
| `finite_strain` | `strain(u)` is the Green-Lagrange strain `E = (FᵀF − I)/2`, `F = I + grad(u)`. `div` in the equation is on the reference shape, and the user writes the first Piola-Kirchhoff stress `F*S`. | `small_strain` |
| `isotropic` | material tensors from `hooke` and `fourier` are isotropic | `orthotropic` |
| `orthotropic(R)` | material axes are the columns of the rotation `R` | `isotropic` |
| `plane_stress(t)` | a `2d` model only. σzz = σxz = σyz = 0. The compiler solves the constitutive law for εzz and removes it (§11.2). Thickness `t` multiplies volume and surface integrals. | `plane_strain`, a `3d` model |
| `plane_strain(t)` | a `2d` model only. εzz = εxz = εyz = 0. With no t, integrals are per unit depth (§14.5). | `plane_stress` |
| `incompressible` | the physics must have an equation `div(v) = 0`, or a pressure unknown that enforces `tr(strain(u)) = 0` | `compressible` |
| `compressible` | density is a varying value or an unknown | `incompressible` |
| `newtonian` | the fluid stress must be linear in `sym(grad(v))` | |
| `laminar` | no turbulence model | `turbulent(model)` |
| `turbulent(model)` | adds the closure `model` | `laminar` |
| `boussinesq(T0, beta)` | density is constant except in the body force term, `rho*(1 − beta*(T − T0))` | `compressible` |
| `isothermal` | no energy equation is expected. A `T` in an expression is a constant. | |

Two assumptions in one "Conflicts with" pair are an error
[`assumption_conflict`]. An assumption whose condition fails is an error
[`assumption_invalid`], for example `newtonian` with a stress that is not
linear in the rate of strain.

### 11.2 How `plane_stress` reduces a 3D law

The user writes the 3D constitutive law. The compiler makes the 2D law:

1. It expands σ(ε) to the 3D components, with εzz as one more unknown at
   each point.
2. It solves σzz(εxx, εyy, εxy, εzz) = 0 for εzz. For a law linear in ε this
   is exact elimination. For a nonlinear law the generated code solves it
   by a local Newton iteration at each integration point.
3. It puts the result back into σxx, σyy, σxy.

For `hooke(E, nu, ε)` the result is the usual plane-stress law with E and
ν. Thus one law serves 3D, plane strain and plane stress.

### 11.3 Axisymmetry

an `axi` model uses the coordinates (r, z) and the names `x` for r and `y` for
z. `grad`, `div`, `strain` and the volume integral include the hoop terms
and the factor 2πr. A vector has two components. A `Tensor2` has the
components rr, rz, zr, zz and θθ.

---

## 12. Boundary conditions and loads

### 12.1 Form

```skein
bc root: u = 0
bc sym: u.x = 0
bc tip: sigma*n = (0, 0, -1) MPa
bc wall: dot(q, n) = h*(T - T_inf)
```

A `bc` names a region and states an equation on it. The checker sorts each
`bc` into one kind by the shape of its left side:

| Left side | Kind |
|---|---|
| an unknown, or a component of it | essential (Dirichlet) |
| `cross(n, u)` of a vector unknown `u` | essential, on the tangential part of `u` only |
| `F*n` or `dot(F, n)`, where `div(F)` is in an equation of the physics | natural (Neumann). The right side may use the unknowns (Robin). |
| `cross(n, F)`, where `curl(F)` is in an equation of the physics | natural, as `F*n` |
| anything else | an error [`bc_form`] |

A natural condition on the flux of the equation paired with unknown `u` is a
condition for `u`.

`cross(n, u) = cross(n, g)` holds the tangential part of `u` at the
tangential part of `g`, and leaves the normal part free. On a wall with
the normal `n`, `cross(n, ef) = 0` is a perfect electric conductor. A
component such as `ef.x` is along an axis, so it cannot state this on a
curved wall. Reason: the conductors of a coaxial line are curved
(femahoy gap G4).

```skein
bc inner_wall: cross(n, ef) = 0
bc port_out: cross(n, hc) = im*beta/mu*ef
```

The right side of an essential condition may be an unknown of another
physics, with the same type, on the same points: a **tie**.

```skein
bc beams: w_b = w
bc beams: theta_b = theta
```

The mesher puts the nodes of the `Curve` `beams` on nodes of the shell,
and the two fields share these dofs. Reason: a stringer and a frame share
the skin points on their lines. The offset of the beam (§4.3) is in the
beam, not in the tie.

On a `Point` region, `n` has no value. There the checker reads a natural
condition `F*n = P` as the force `P` at the point. `P` has the dimension of
`F` times area in a `3d` model, and of `F` times length in a `2d` or `axi`
model.

In a physics on a `Curve`, `along(F)` in an equation makes `F` a flux, as
`div(F)` does. There `F*n = P` stands only on a `Point`, and is the force
`P` at the point, with the dimension of `F`. A flux of `along` on a region
that is not a `Point` is an error [`bc_form`].

```skein
bc tip: F_b*n = (0, 0, -1)*P
bc tip: M_b*n = (0, M, 0)
```

Reason: a beam takes a force and a moment at a point, and the
implementation.md §10 cantilever puts its load on the tip.

The right side `rigid(p)` of an essential condition on a whole displacement
`u: Vector<Length>` makes the region a **rigid face**: the face stays
straight and can turn. Its displacement is the displacement of the `Point`
`p` plus a small rotation about `p`. Other uses of `rigid` are errors
[`rigid_form`].

```skein
let c: Point = point((0 mm, hp))
bc top: u = rigid(c)
bc c: u.x = 0
bc c: sigma*n = (Q, -P)
bc c: moment(sigma*n) = hp*Q
```

A `bc` on `p` gives the force and the moment on the face. An essential
condition on `u` fixes the displacement of `p`, and `F*n = P` is the force
`P` at `p`, as on each `Point`. `moment(F*n) = M` is the moment `M` about
`p`: a scalar about z in a `2d` model, and a vector in a `3d` model. `M`
has the dimension of the force times length. A `moment` on a region that
is not the point of an earlier `rigid` is an error [`rigid_form`]. With no
`moment`, the face turns freely about `p`.

Reason: a test puts its load on a face through a rigid plate, and the
force and the moment of the plate are in the test (femahoy gap G17,
CalculiX `*EQUATION`).

The right side `radiation(eps)` of a natural condition on a heat flux is
the net gray diffuse radiation out of the face, with the emissivity `eps`
(a constant `Real`). The faces of all such conditions of a step make one
enclosure, and they must close it. The view factors between the facets of
the mesh are from the crossed strings of Hottel, with the facets that hide
a pair in part. The temperature must be absolute (K). Other uses of
`radiation` are errors [`radiation_form`].

```skein
bc rod_face: dot(q, n) = radiation(0.8)
bc tube_face: dot(q, n) = radiation(0.5)
```

Reason: the radiation between surfaces is non-local, so the strong form
cannot state it (decision G4.2). A condition and not a separate item:
radiation is the flux of a face, like the other natural conditions.
`skein run` solves it in a `2d` model.

### 12.2 Where `bc` may stand

A `bc` at the top level applies in all steps. A `bc` in a step applies from
that step on, and replaces a `bc` on the same region and target from an
earlier step. `release r` in a step removes the conditions on `r`.

Reason: this is the propagation rule of Abaqus steps, which FEA users know.

### 12.3 Completeness

A part of the boundary with no `bc` for an unknown gets the natural
condition with a zero right side (traction-free, insulated, zero flux). The
`check` report lists these parts, so the user sees them.

These cases are errors:

- two essential conditions on one region and one component, with different
  right sides, in one step [`bc_conflict`];
- an essential and a natural condition on one region and one component
  [`bc_conflict`]. On the `Point` `p` of a `rigid(p)` only, a force with 0
  in the component of an essential condition does not fight it;
- a `bc` on a region that is not on the boundary of the region of the
  unknown, and is not a `Point` or a lower-dimension entity inside it
  [`bc_off_boundary`].

The pre-solve check [`rigid_body`] tests that the essential conditions of
each step remove all rigid-body motions. The check [`pressure_level`] tests
that an incompressible flow with essential velocity on the whole boundary
fixes the pressure at a point or by its mean.

### 12.4 Contact

```skein
contact grip: jaw_face, part_face {
    friction: 0.2
    enforce: mortar
}
```

A `contact` item names two boundary regions that can touch. The first is the
master, the second the slave. `enforce` is `penalty` (with `stiffness`),
`nitsche` or `mortar`. Contact makes the system nonlinear, so the step uses
Newton's method.

### 12.5 Values from other models

The right side of an essential condition may be a field of a `Results`
value (§15.1). The solver interpolates it at the nodes, or at the
integration points for a weak condition.

### 12.6 Symmetry groups

```skein
symmetry cell: bloch { period: (1, 0, 0) m, low: x0, high: x1, k: 0 }
```

A `symmetry` item solves one cell of a body that repeats. It ties the
unknowns on the boundary `high` to the unknowns on `low`. The one group is
`bloch`, with the settings `period` (a vector), `low`, `high` and `k` (a
wave vector). The value on `high` is the value on `low` multiplied by
`exp(im*dot(k, period))` (Bloch and Floquet). A `k` other than 0 needs
`Complex` unknowns (§5.1).

The checker tests that the shift by `period` maps `low` onto `high`
[`symmetry_map`]. Reason: a `period` that does not match the cell ties
the wrong points and gives a wrong result with no error.

Reason for `bloch`: the periodic cells of femahoy (grid-laminate-cell,
grid-fibre-composite-cell) use it. No femahoy case needs a rotational
group.

---

## 13. Methods

```skein
method solid: cg {
    mesh: m
    order: 2
}
```

`method p: kind { … }` gives the numerical method of physics `p`. A
physics with no `method` item takes the defaults:

- the kind `cg`, on the one mesh whose region holds the region of the
  physics. With no such mesh, or with more than one, the default does not
  apply [`method_missing`].
- `order: 2` for each unknown, except a pressure that pairs with a
  velocity under `incompressible`, which takes `order: 1` (Taylor-Hood);
- `stabilize: supg` when the equations have a convection term, and `none`
  when they do not.

A `cg`, `dg` or `fv` method with no `mesh` setting takes the same mesh as
the default.

Reason: these defaults are the right choice for most solid, heat and
laminar flow models, and they cut a block of 4 lines from each of them.
`skein check` prints the method that it chose.

The kinds:

| Kind | Settings | Notes |
|---|---|---|
| `cg` | `mesh`, `order` (Int or a record per field), `family` (`lagrange`, `serendipity`, `incompatible`, `nedelec`, `raviart_thomas`), `stabilize`, `quadrature` (`full`, `reduced`), `couple`, `pml` (`{ region, thickness }`), `gauge` (`tree`) | continuous Galerkin FEM. `nedelec` is H(curl) and `raviart_thomas` is H(div). |
| `dg` | `mesh`, `order`, `flux` (`upwind`, `lax_friedrichs`, `roe`), `penalty` | interior penalty for second-order terms |
| `fv` | `mesh`, `gradient` (`green_gauss`, `least_squares`), `flux` (`upwind`, `linear_upwind`, `central`), `limiter` (`none`, `minmod`, `van_leer`, `venkatakrishnan`), `pv_coupling` (`rhie_chow`) | cell-centred finite volume |
| `rkpm` | `points`, `basis` (`linear`, `quadratic`), `support` (Real, times the spacing), `integrate` (`scni`, `gauss`), `bc` (`nitsche`) | Galerkin meshless |
| `gfd` | `points`, `stencil` (Int), `weight` | strong-form collocation |
| `ode` | none | for a physics with only global unknowns (§10.5). The time scheme comes from the step. |

The setting `stabilize` takes `none`, `supg`, `supg_pspg`, `vms` or
`supg_pspg_lsic`. §17.2 gives the meaning.

The setting `family` of `cg` gives the cells of each order:

- `lagrange` (the default): the full tensor-product or simplex nodes of
  the order.
- `serendipity`: a quad or hex of order 2 with nodes only on the corners
  and the edges: 8 nodes on a quad (CPE8), 20 nodes on a hex (C3D20).
- `incompatible`: a quad or hex of order 1 with the incompatible modes of
  Wilson et al. (1973), with the correction of Taylor et al. (1976), so
  that the cell passes the patch test (CPE4I, C3D8I). The modes
  1 - ξ² along each cell axis are internal to each cell, and the solver
  condenses them out. They remove the shear locking of order 1 in bending.
- `nedelec`: the edge elements of Nédélec (1980), first kind, on
  triangles and tetrahedra. Order 1 has one dof on each edge, the
  integral of the tangential component of the vector unknown along it.
  The tangential part is continuous between cells and the normal part is
  not. A `bc` `cross(n, u) = cross(n, g)` holds these edge dofs.
- `raviart_thomas`: the face elements of Raviart and Thomas (1977) for a
  mixed physics with a vector `q` and a scalar `p`, such as
  `q + k*grad(p) = 0` and `div(q) = f`. Order 1 has one dof on each facet,
  the flux of `q` through it, and `p` is constant in each cell. The normal
  part of `q` is continuous between cells. In the mixed form, `p = v` is
  the natural condition and `dot(q, n) = g` is the essential one. A
  boundary with neither has no flux.

The setting `gauge: tree` of `nedelec` holds at zero the edges of a
spanning tree of the cells where the coefficient of `u` in the equation
is 0. There the curl-curl operator has the gradients in its null space,
and the tree makes the solution unique. The edges of the essential
conditions are in the tree first. Only `curl(u)` is unique there, not
`u`. `gauge` with another family is an error [`method_setting`].
Reason: an eddy current model has air and coils with no conductivity
(femahoy gap G27).

The setting `quadrature` takes `full` (the default), the Gauss rule that
integrates the stiffness of an undistorted cell exactly, or `reduced`,
one point fewer along each axis of a quad or hex and one order lower on a
simplex (C3D8R, C3D20R, CPE8R, S8R). With order 1 on a quad or hex,
`reduced` adds hourglass control (Flanagan and Belytschko 1981).

`serendipity` needs `order: 2` and `incompatible` needs `order: 1`, both
on a mesh with `cell: quad` or `cell: hex`. `incompatible` needs `full`
quadrature, because the gradients of its modes are 0 at the one centre
point. Another combination is an error [`family_cell`]. A word that is not a
value of `family` or `quadrature` is an error [`method_setting`].

The setting `pml: { region: layer, thickness: 1 um }` of `cg` makes
`region` a perfectly matched layer (PML) of an open boundary (Bérenger
1994). The layer is the part of the body within `thickness` of its outer
faces. In the layer, the solver stretches each coordinate axis along
which the point is past the rest of the body by the complex factor
s = 1 - im·σ, with σ = (3 ln(10⁶)/(2 k h))·(d/h)². Here d is the depth of
the point past the rest of the body along the axis, h is `thickness`, and
k is the wave number in the layer: the square root of the coefficient of
the unknown in the equation. An outgoing wave decays in the layer, and the
reflection of the layer at normal incidence is 10⁻⁶ before the
discretization. The outer faces of the layer keep their own conditions.
`region` and `thickness` are both required, and the physics needs a
`Complex` unknown [`method_setting`]. Reason: the antenna, the eddy current
and the optics models of femahoy need an open boundary (gap G26), and a
first-order absorbing condition reflects a wave that is not normal to it.

A setting that does not belong to the kind is an error
[`method_setting`]. The checker tests the method against the physics:

- an inf-sup unstable pair (equal order for velocity and pressure in `cg`)
  with no `pspg`, `vms` or `supg_pspg` is an error [`inf_sup`];
- a `cg` order for a field that the equations differentiate twice outside
  a `div` is an error [`order_too_low`];
- a `gfd` method with a physics that is not in a form that it can collocate
  is an error [`method_physics`].

---

## 14. Steps

### 14.1 Form

```skein
step heat_up: transient(end: 600 s, dt: 5 s, scheme: bdf2) {
    bc wall: T = 400 K
    output T every 10
}
```

The steps run in the order of the file. Each step starts from the state at
the end of the step before it. The first step starts from the initial value
of each unknown: the expression after `=` in `solves`, or zero.

```skein
physics thermal on part solves T: Temperature = 293 K { … }
```

### 14.2 Procedures

| Procedure | Arguments | Meaning |
|---|---|---|
| `static` | `ramp` (Int, default 1) | `ddt` = 0. Nonlinear problems take `ramp` load increments. |
| `transient` | `end`, `dt`, `scheme` | time stepping from the end time of the step before |
| `harmonic` | `freq` (a list of `Frequency`), `damping` (Real, default 0) | steady response to loads of the form `F·exp(im·ω·t)`, at each frequency of the list. The unknowns and loads are `Complex` (§5.1). |
| `modal` | `count`, `shift` | the `count` lowest eigenpairs of the linearized operator about the current state |
| `buckling` | `count` | eigenpairs of the stiffness and the geometric stiffness of the current state |

In a `harmonic` step, `ddt` is im·ω and `ddt2` is -ω². As in a `modal`
step, `linear_elastic` takes its inertia from the density `rho` of the
material. `damping: ζ` of `harmonic` is the damping ratio of each mode.
The solver multiplies the stiffness (the terms of the equations with no
`ddt`, `ddt2` or `omega`) by 1 + 2·im·ζ: structural damping, which gives
the damping ratio ζ at each resonance. Reason: the random vibration model of femahoy gives a
damping ratio for each mode (gap G2). On its cantilever, the RMS
response with this damping is 0.005 % from the one of viscous modal
damping.

The built-in constant `omega` is the angular frequency of the step: 2π
times the frequency that a `harmonic` step solves, and ω of the mode that
a `modal` step solves. In an `expect` after the steps, it is that of the
last frequency of the last `harmonic` step. In the other steps it is 0.
A load, an impedance or a coefficient that depends on the frequency uses
it, for example the surface impedance `(1 + im)*sqrt(omega*mu0/(2*sigma))`.
Reason: femahoy has a base acceleration and a wall impedance that depend
on the frequency (gap G8).

Time schemes:

| `scheme` | Order | Kind | Default for |
|---|---|---|---|
| `backward_euler` | 1 | implicit | |
| `bdf2` | 2 | implicit | first-order in time (`ddt`) |
| `crank_nicolson` | 2 | implicit | |
| `generalized_alpha(rho_inf)` | 2 | implicit | second-order in time (`ddt2`) |
| `ssprk3` | 3 | explicit | |
| `central_difference` | 2 | explicit, lumped mass | |

An explicit scheme with a `dt` above the stable limit gives the lint
[`dt_unstable`] (checks.md).

`dt: stable` tells the solver to choose the stable increment of an
explicit scheme. The solver estimates the limit from the highest frequency
of the mesh, and takes a step a little below it.

```skein
step pull: transient(end: 80 us, dt: stable, scheme: central_difference) { … }
```

`dt: stable` needs `scheme: ssprk3` or `scheme: central_difference`.
Otherwise the error is `procedure_args`.

### 14.3 Birth and death

```skein
step lift_2: static {
    birth lift2
}
step burn: transient(end: 10 s, dt: 0.1 s) {
    death fuel where T > 900 K
}
```

- `death r` removes the cells of region `r` from the active set at the
  start of the step. Their unknowns leave the system. A node that only dead
  cells hold leaves the system too. The removed cells keep their state for
  output.
- `birth r` adds the cells of `r` to the active set. The default is
  **strain-free birth**: the cells start with zero strain on the current
  deformed shape. The nodes of the new cells that touch active cells take the
  current values of those cells. The other new nodes take the values that
  the active nodes give by extrapolation along the boundary.
- The compiler makes strain-free birth as follows. Under `small_strain`,
  each born cell stores `eps_birth = strain(u)` at each integration point
  at the start of its birth step. In the cells of `r`, the compiler then
  replaces `strain(u)` with `strain(u) - eps_birth`. Under `finite_strain`,
  each born cell stores `F_birth = defgrad(u)` at its birth. The compiler
  then replaces `defgrad(u)` with `defgrad(u) * inv(F_birth)` in the cells
  of `r`. Thus the stress of a born cell is zero at its birth.
- `where cond` selects the cells of `r` where `cond` is true at the start of
  the step. `cond` may use the fields.
- The mesh of the whole model, alive or dead, is made once, before the first
  step. The region of a `birth` is inside the region of its unknown
  [`region_outside`] (§8.4), so the mesh holds it.

Reason for removal (not a small stiffness factor): the removed unknowns
give no ill conditioning to the solver.

### 14.4 `solve`

```skein
solve { nonlinear: newton, linear: gmres, precond: ilu0, rtol: 1e-8, max_iter: 50 }
```

Each setting is optional. The defaults:

| Setting | Default |
|---|---|
| `nonlinear` | `newton` if the system is nonlinear, else `none` |
| `linear` | `cg` if the system is symmetric positive definite, else `gmres` |
| `precond` | `amg` (phase 3) or `ilu0` (before phase 3) |
| `rtol` | 1e-8 |
| `max_iter` | 50 for Newton, 1000 for the linear solver |

All solvers are in Goose (§18.4).

A static step with `linear: cg` and `precond: amg` on a mesh of a lattice
(a box or a union of boxes of hexahedra, a box or a plate of simplices)
solves with a sparse direct factorization when it pays: its work is at
most 1e12 multiplies and adds and its factor fits in the free memory.
Else, or if the factorization fails (K is not positive definite), it solves
with CG and AMG. The program chooses, so the model needs no
setting for it. Reason: on `interposer_warpage` (139,449 dofs, thin layers)
CG with AMG takes 677 iterations and 268 s on 1 thread, and the direct
solve 11.6 s (in the page 43 s, and 3 s with the factor of the last run).

`linear: gmres` with `precond: ilu0` uses GMRES with a fixed restart of
300. It is for small models only, because each solve holds about 301 basis
vectors of n doubles for n unknowns. A finer mesh needs more iterations,
and the restart can then make GMRES stall. A solve that does not converge
writes the solver, the iterations and the restart to stderr.

### 14.5 `output` and `report`

`output e1, e2, …` writes the fields or field expressions to the result
file at the end of the step, or each `every k` time steps. `report e` prints
a scalar or a vector to the run log in TOON. Functions for `report`:
`reaction(r)`, `max(f, r)`, `min(f, r)`, `argmax(f, r)`, `integral(f, r)`,
`mean(f, r)`, `eigenvalue(k)`, `when(f, value)`, `rms(f, psd)`,
`far_field(E, H, s, d)`.

`argmax(f, r)` is the point of `r` where the scalar `f` has its largest
value, a `Vector<Length>` [`rank_mismatch`]. If `r` comes from a `pattern`
or a `pattern_polar` (§8.1), `max(f, r)` and `min(f, r)` give one value for
each copy, a list in the sequence of the copies. A region operation keeps
the copies of its left side. Reason: a report must give the largest axial
force of each stringer, not only of all of them.

`eigenvalue(k: Int)` is the eigenvalue of mode `k` of a `modal` step, or
of a `buckling` step. For a `modal` step, it is the natural frequency
f = ω/(2π) in Hz, a `Frequency`. For a `buckling` step, it is a `Real` load
factor. Mode 1 is the lowest. In a `report` or an `expect` of a step, the
step is that step. In another `expect`, the step is the last `modal` or
`buckling` step. A use in another place, or with no such step, is an error
[`eigenvalue_step`]. A `k` outside `1..count` of the step is also an error
[`eigenvalue_step`].

`when(f, value)` is the first time at which `f`, a scalar reduced over a
region, gets to `value`, a `Time` [`rank_mismatch`]. `value` has the
dimension of `f` [`unit_mismatch`], and it must be a constant
[`not_constant`]. The run finds the first time step at which `f` crosses
`value`, and interpolates linearly in that time step. If `f` never gets to
`value`, an `expect` of it fails. `when` is valid only in a `report` or an
`expect` of a `transient` step, or in another `expect` whose last step is
`transient` [`when_step`].

```skein
let front: Length = max(x + u.x, column)
expect z2: when(front, 2.11*a) - when(front, 1.44*a) == 0.2 s within 7%
```

Reason: a dam break test gives the times at which the front of the water
gets to positions on the floor (femahoy gap G11).

`rms(f, psd)` is the root mean square of a random response: the square
root of the integral of |f|²·psd over the frequencies of a `harmonic`
step, by the trapezoidal rule on the list `freq`. `f` is the scalar
response to a load of unit amplitude, at a point or reduced over a region,
and it can be `Complex`. `psd` is the power spectral density of the load,
a scalar that can use `omega`. Both must be scalars that do not vary in
space [`rank_mismatch`]. The result has the dimension of
f·sqrt(psd·Frequency): with a load of amplitude `g` and a PSD in g²/Hz,
`psd` is the PSD divided by g², per `Frequency`. `rms` is valid only in a
`report` or an `expect` of a `harmonic` step, or in another `expect`
whose last step is `harmonic` [`rms_step`].

```skein
expect a_rms: rms(g - omega^2*mean(u.z, tip), 0.04/(1 Hz)) == 400925.4 mm/s^2 within 0.5%
```

Reason: the random vibration model of femahoy checks the RMS acceleration
of the tip of a cantilever for the PSD of a base acceleration (gap G2).

`far_field(E, H, s, d)` is the radiation intensity, the power per unit
solid angle, in the direction `d` far from the closed `Surface` `s`. `E`
and `H` are the `Complex` electric and magnetic fields on `s`. The
transform takes the surface currents n × H and -n × E on `s` into the
medium outside `s` (Balanis 2016, §12.2), with the wave number and the
impedance of the material on `s`. The result has the dimension of
`E`·`H`·`Length^2`. The directivity is 4π times the intensity over the
radiated power, `integral(re(dot(cross(E, conj(H)), n))/2, s)`. `d` is a
constant vector, and only its direction counts. `far_field` is valid only
in a `3d` model, in a `report` or an `expect` of a `harmonic` step, or in
another `expect` whose last step is `harmonic` [`far_field_form`].

```skein
expect directivity: 4*pi*far_field(ef, hf, ntff, (1, 0, 0))/integral(re(dot(cross(ef, conj(hf)), n))/2, ntff) == 1.641 within 6%
```

Reason: the antenna model of femahoy checks the directivity and the
efficiency of a dipole (gap G13).

`report cd = expr` gives the value a name. A study (§20) uses the name as
the column of its results table. A report with no name gets its source
text as the column name.

`integral(f, r)` has the dimension of `f` times `Length^k`, where `k` is
the topological dimension of `r`. In a `2d` model with no thickness, a force
on a curve is thus `Force/Length`, the force per unit depth.

---

## 15. Results, submodels and expectations

### 15.1 `results`

```skein
let global: Results = results("01_cantilever.skn", step: "load")
```

`results(file, step:)` refers to the state of another model at the end of
the step. The checker reads the other file and gives `global.u` the type of
the field `u` in it. The `let` fields of a physics of the other model are
fields of the `Results` value too, for example `global.e_hoop`. The run uses the saved state, or runs the other model
first if no state exists.

```skein
param p_U: Pressure = 97.5 kPa
let global: Results = results("fuselage.skn", step: "pressurize", case: p_U)
```

`case: p` names a `param` of this model, declared before the `results`.
The other model must declare a `param` of the same name and unit
[`results_case`]. Each case of a study of this model then uses the state
of the other model run with the same value of `p`, as with
`skein run <file> --set p=<value>`. The cases that need the same state
run the other model once.

### 15.2 Submodels

A submodel is a model whose essential conditions on a **cut boundary** come
from a `Results` value:

```skein
bc cut: u = global.u
```

The checker tests:

- the region of the cut is inside the region of the field in the other
  model [`submodel_outside`];
- the field types are equal [`type_mismatch`];
- the unit systems may differ. The values are in SI.

The interpolation evaluates the global field at the submodel nodes with the
shape functions of the global method. For a global `rkpm` model, it uses the
reproducing kernel functions. The post-solve check [`submodel_cut`] compares
the resultant force on the cut with the global stress resultant on the same
surface.

### 15.3 `expect`

```skein
expect kirsch: sigma.xx at (0 mm, a) == 3*S within 5%
```

An expectation is a test on the results. The two sides are scalars or
vectors at points, or two field expressions on a region. For two field
expressions, the test compares the L2 norm of the difference with
`within`. This is how a model states a manufactured or analytic solution
for `skein verify convergence` ([command line](command_line.md)). `within`
takes an absolute value with the dimension of the sides, or a percent of
the right side. A table column in the expression makes one test for each
row. A failed expectation gives a non-zero exit status. The run reports
each expectation as a post-solve verification result. The place of the
`expect` sets when the run evaluates it:

| Place | The run evaluates it |
|---|---|
| an item | after the last step |
| a statement of a step | at the end of the step, after its reports |
| an entry of a `study` (§20.2) | after the last step of each case. In a `minimize` study, only in the best case. |

```skein
step inflate_25: static(ramp: 20) {
    bc inner: u = (15 mm)/A_in*(x, y)
    expect f_25: reaction(inner - sx).x == 5.74269 N within 0.2%
}

study ra_sweep: grid {
    Ra: [10000, 100000]
    expect nu: mean(-grad(T).x, hot)*a/dT == [2.243, 4.519] within 0.6%
}
```

In a `grid` or a `zip` study, the right side and `within` can each be a
list with one value for each case, in the order of the cases (§20.2). A
list in another place, or a list with a different number of values, is an
error [`expect_cases`]. Each case of a study also evaluates the `expect`
items and the `expect` statements of the steps, with the values of the
case.

Reason: a test can need a check after an earlier step, a check for each
case of a sweep with its own reference value, and a check of the optimum of
a shape study (femahoy gap G10).

---

## 16. The intermediate form (SIR)

The SIR is the contract between the front end and the methods. It is
method-neutral and has no mesh. It is a typed, unit-free (SI-scaled) list of
balance laws.

### 16.1 Balance-law form

The lowering writes each equation of each physics in this form:

```
M(U) ∂ₜ U  +  N(U) ∂ₜ² U  +  div F(U, ∇U)  +  S(U, ∇U)  =  0
```

- `F` is the flux: the sum of the terms in `div` in the source. The
  lowering also keeps a split `F = F_a(U) + F_d(U, ∇U)`: the advective part
  has no derivative of an unknown, the diffusive part has. Only the flux
  choices of `dg` and `fv` and the stabilization (§17.2) use this split.
  The natural conditions (§12.1) and the `cg` lowering (§17.1) use the full
  flux `F`.
- `S` holds all other terms, for example the body force and the convective
  term `grad(v)*v`.
- `M` and `N` are the coefficients of the time derivatives.

The SIR also keeps the **full strong residual** `R(U)` of each equation.
Residual-based stabilization (§17.2) and `gfd` (§17.5) need it.

### 16.2 SIR records

```
SirModel   { space, scale, fields[], equations[], conditions[], regions[], steps[], expects[] }
SirField   { name, rank, dim, region, global }
SirEq      { field, region, M, N, Fa, Fd, S, R, linear: Bool, symmetric: Bool }
SirCond    { region, field, component, kind: essential|natural, rhs }
SirRegion  { name, dim, csg }
SirStep    { name, procedure, args, conditions, active_set_changes, outputs, reports }
```

Each expression in the SIR is a typed tree of tensor operations over the
fields, their gradients, the coordinates and time. The compiler derives the
Jacobian of each equation by symbolic differentiation of this tree, so
Newton's method has the exact tangent.

The flags `linear` and `symmetric` come from the tree: an equation is
symmetric if its linearization is self-adjoint. The defaults of `solve`
(§14.4) use them.

---

## 17. Lowering by method

### 17.1 Continuous Galerkin (`cg`)

For equation `k`, paired with unknown `u_k` and the test function `w`:

1. Multiply the balance law by `w` and integrate over the region.
2. Integrate the full `div F` by parts when `F` holds a derivative of an
   unknown, or when a `bc` states a natural condition on `F`:
   `∫ w div F = −∫ grad(w) : F + ∮ w F·n`.
3. Keep each other `div F` as it is. Reason: such a flux has no derivative
   of an unknown and no natural condition, and integration by parts brings
   no gain for `cg`.
4. Replace the boundary term `∮ w F·n` with the natural conditions. Drop
   it where an essential condition holds, because `w` is zero there.

The rule in step 2 is the one rule for all physics. Thus `div(sigma)` with
`sigma` from `grad(u)` goes by parts, and `div(v)` in the continuity
equation does not. The full flux goes by parts, also its terms with no
derivative. Thus the natural condition `sigma*n` includes the thermal term
of `hooke(E, nu, strain(u) - alpha*(T - T0)*I)` and the pressure term `-p*I`
of a flow stress.

### 17.2 Stabilization

The SIR keeps the strong residual `R`, so the compiler builds the stabilized
forms exactly:

| `stabilize` | Added term (sum over cells) |
|---|---|
| `supg` | `τ_M (v·∇w, R_M)` |
| `supg_pspg` | `τ_M (v·∇w + ∇q/ρ, R_M)` |
| `supg_pspg_lsic` | the above plus `τ_C (div w, div v)` |
| `vms` | the residual-based VMS terms of Bazilevs et al. (2007), with the fine-scale velocity `−τ_M R_M` |

`τ_M = ((2/Δt)² + (2|v|/h)² + C_I (4ν/h²)²)^(−1/2)` and
`τ_C = 1/(τ_M · tr(G))`, after Tezduyar (1992), Shakib et al. (1991) and
Bazilevs et al. (2007). `G = (∂ξ/∂x)ᵀ(∂ξ/∂x)` is the metric tensor of the
cell, where `ξ` are the coordinates of the reference cell. `G` has the
dimension `1/Length^2`, so `τ_C` has the dimension `Length^2/Time`.
`R_M` is the residual of the momentum equation and `q` the pressure test
function. In a `static` step, the `2/Δt` term is zero.

### 17.3 Finite volume (`fv`)

The lowering integrates each balance law over each cell and applies the
divergence theorem to all `div` terms, advective and diffusive. The face
fluxes use the `flux` and `gradient` settings. Pressure-velocity coupling
uses Rhie-Chow interpolation (Rhie and Chow 1983) of the face velocity. The
system is solved fully coupled by Newton's method. The unknowns are cell
averages. Essential conditions act through the face flux.

The current lowering (`src/grid/`) takes one scalar unknown with
`ct ddt(u) + div(a*u) - d*div(grad(u)) = f` on a `rect` of `quad` cells,
with `ct`, `a` and `d` constant. The diffusive flux uses the two cell
values. `linear_upwind` adds the cell gradient to the upwind value, with
the face limiter `ψ(r)` of Darwish and Moukalled (2003) for `minmod` and
`van_leer`, and the cell limiter of Venkatakrishnan (1995) with `K = 1`.
An essential condition is a ghost cell. A `static` step is a defect
correction: each iteration solves the first order system with GMRES and
ILU(0). A `transient` step uses `ssprk3`.

### 17.4 Discontinuous Galerkin (`dg`)

The lowering applies §17.1 on each cell. The boundary terms on interior
faces use the numerical flux of the `flux` setting for `F_a` and the
symmetric interior penalty form (Arnold 1982, Arnold et al. 2002) for
`F_d`. Essential conditions are weak, by the same penalty.

The current lowering (`src/grid/`) takes the equation of §17.3 on the same
grid. The basis is the tensor Legendre basis `Q_p`, `order` 1 to 3, so the
mass matrix is diagonal. The penalty is `σ = penalty · d · (p+1)²/h`, with
the default `penalty: 2`. `roe` is `upwind` for one scalar.

### 17.5 Meshless

- `rkpm` uses the Galerkin form of §17.1 with reproducing kernel shape
  functions (Liu, Jun and Zhang 1995, Chen et al. 1996). The integration is
  stabilized conforming nodal integration (Chen et al. 2001). The shape
  functions do not interpolate, so essential conditions use Nitsche's method.
- `gfd` collocates the strong form at each point with generalized finite
  difference stencils (Benito, Ureña and Gavete 2001). It uses `R` from the
  SIR directly, with no weak form. The natural conditions are collocated at
  the boundary points.

The current lowering (`src/meshless/`) takes one scalar unknown with
`div(a*u) - d*div(grad(u)) = f` and static steps, and `rkpm` also takes 2D
linear elasticity in plane stress or plane strain.

- `rkpm` has `basis: linear`, the cubic B-spline kernel and the support
  `support` (default 2.5) times the spacing of each point. The spacing of
  a point is the largest length of the middle edge of its triangles: the
  lattice spacing of a lattice, and the larger spacing of a stretched
  triangle. The nodal
  domain of each point is cut into one quadrilateral in each of its
  triangles, and the stiffness takes the smoothed gradient of each
  quadrilateral (stabilized conforming nodal integration in the form of
  Puso, Chen, Zywicz and Elmer 2008). A fixed value is weak, by the
  symmetric method of Nitsche with the penalty 4 (λ + 2μ)/h or 4 d/h. The
  stress at a point is from the smoothed gradient there.
- `gfd` fits a Taylor polynomial of order 2 to the `stencil` (default 12)
  nearest points of each point, with the weight 1/r³. A point on a free
  edge collocates `(a·n) u - d ∂u/∂n = 0`.
- The solve is direct: a profile LU factorization in the reverse
  Cuthill-McKee order of the points.

### 17.6 Discrete IR

Each method writes a discrete IR: a list of **kernels** and a **plan**.

- A kernel computes the contribution of one cell, face or point to the
  residual and to the Jacobian. Its inputs are flat arrays of coordinates
  and unknowns. Its outputs are flat arrays of values.
- The plan gives the loops over cells, faces and points, the active sets of
  each step, the assembly into a sparse matrix, the solver calls, and the
  outputs.

---

## 18. Code generation to Goose

### 18.1 Why Goose

[Goose](https://aardappel.github.io/goose/) is memory-safe with no garbage
collector and no heap allocator, and it compiles to C. A Skein program
becomes a Goose program, so the generated code is memory-safe too. The
Goose compiler then writes C, and a C compiler makes a native executable or,
with Emscripten (`emcc`), a WebAssembly module.

The WebAssembly target uses the C that Goose writes: `emcc` compiles that C
to a WebAssembly module. Earlier projects of the author of Skein, beamdiag
and goosestrike, already run in the browser in this way. Thus the browser
target does not depend on a WebAssembly back end in Goose.

### 18.2 Files

`skein build model.skn` writes:

| File | Contents |
|---|---|
| `model.goose` | the plan of the model |
| `runtime` (package) | the runtime: mesher, sparse matrices, solvers, output writers, and the cell kernels that the emitter writes (§18.3). The same for each model. |

### 18.3 The shape of the generated code

The plan stores a mesh as flat arrays, which Goose keeps in data stacks with
no pointers. A cell kernel is a `simd fn`. Assembly runs in worker threads
on blocks of cells, and the thread of the case adds the cell matrices
(§20.5). Thus no two threads write one row, and the cells need no
colouring.

Each generated function has at most 60 lines (rule 4 of the Power of
Ten). When the kernel of a cell is longer, the emitter splits it into
functions by blocks of the cell matrix, one block for each pair of
unknowns or of node groups, and one function calls them in order. Reason:
a person can review a function of 60 lines, and a test checks the length
of each generated function.

A cell kernel depends only on the cell type and the physics. So the
emitter (`src/codegen/kernels.goose`) writes the kernels into the runtime
packages once, as `src/runtime/kernels.goose`, `src/elements/kernels.goose`
and `src/shells/kernels.goose`, and not into each `model.goose`. A test
fails while one of these files differs from what the emitter writes.

This is the kernel that the compiler writes for the linear elastic
tetrahedron of order 1. It is illustrative: the generated names and the
layout are not part of the specification.

```goose
import vec;

struct Tet4 { n: i32[4] }

inline fn comp(v: double3, i: i64) -> f64 {
    return if i == 0 { v.x } else if i == 1 { v.y } else { v.z };
}

// Stiffness of one linear tetrahedron, 12 x 12, row-major into k.
simd fn k_solid_tet4(x: double3[4], lam: f64, mu: f64, k: f64[144]&) {
    let gn = runtime::tet4_gradients(x);             // double3[4], and the volume
    for a in 4 {
        for b in 4 {
            let ga = gn.g[a];
            let gb = gn.g[b];
            for i in 3 {
                for j in 3 {
                    var v = lam * comp(ga, i) * comp(gb, j) + mu * comp(ga, j) * comp(gb, i);
                    if i == j { v += mu * dot(ga, gb); }
                    k[(3 * a + i) * 12 + 3 * b + j] = v * gn.vol;
                }
            }
        }
    }
}

thread_fn assemble_worker() {
    loop {
        let job = qget<solvers::CellJob>();
        guard job.lo >= 0 else { break; }
        qput(runtime::assemble_cells(job, k_solid_tet4));
    }
}
```

### 18.4 The runtime

`runtime` is a Goose package in `src/runtime/`. GMRES(m), ILU(0) and
Newton are in the package `solvers` in `src/solvers/`. `runtime` calls its
GMRES(m) and ILU(0). `runtime` holds no C, except the C that the Goose
standard library uses for files and the clock. The [roadmap](roadmap.md)
gives the phase of each part.

| Part | Contents |
|---|---|
| sparse | CSR matrices, assembly |
| linear | CG, GMRES(m), Jacobi, ILU(0) |
| nonlinear | Newton with line search, load ramps |
| eigen | LOBPCG (Knyazev 2001) for `modal` and `buckling` |
| multigrid | smoothed aggregation AMG (Vaněk, Mandel and Brezina 1996) |
| mesher | §9.2 |
| output | VTU writer, TOON reports, state files |

### 18.5 The compiler in Goose

The Skein compiler is a Goose program too. Thus `skein` builds for the
browser from the same source as the native tool, and the editor in the
browser runs `check`, `fmt` and `lint` with no server.

---

## 19. Static semantics checklist

`skein check` gives each error of [checks.md §2](checks.md#2-errors). That
catalogue is the one list of the checks, with the section of this document
that states each rule. An `allow` must name a lint rule or a verification
check, with a reason that is not empty [`allow_reason`].

---

## 20. Parameters and studies

### 20.1 `param`

```skein
param Re: Real = 100 in 1..1000
param D: Length = 0.1 m
```

A `param` is a `let` whose value can change from outside the file: from a
`study` (§20.2) or from `skein run --set Re=200`. It has these rules:

- The type is required, as for each declaration.
- The value must be a constant (§7.6) [`not_constant`].
- `in lo..hi` is an optional closed range. The default value, each `--set`
  value and each study value must be in it [`param_range`].
- A `--set` value is a literal of the parameter type, with its unit:
  `--set D=0.2m`. The checker type-checks it as source text.
- Only a `param` can change. A `let` is fixed. Thus a reviewer finds all the
  inputs of a model in its `param` lines.

Reason: a study must know which values may change, and the compiler must
know which parts of the model depend on them (§20.4).

### 20.2 `study`

```skein
study re_sweep: grid {
    Re: [10, 20, 40]
    D: [0.1, 0.2] m
}
```

A study runs the model many times with different parameter values. Each
run is a **case**. The procedure after the colon makes the cases:

| Procedure | Entries | Cases |
|---|---|---|
| `grid` | a list for each parameter | the Cartesian product of the lists |
| `zip` | lists of one length | case `i` takes value `i` of each list |
| `sample(n: Int, method: lhs \| sobol, seed: Int)` | a range `lo..hi` for each parameter | `n` points in the box of the ranges |
| `minimize(objective, method: nelder_mead \| bobyqa \| gradient, max_cases: Int, xtol: Real, step: Real)` | a range for each parameter, and the start value is the default | cases chosen one batch at a time to minimize a named `report` |
| `topology(objective, region: Surface, volume: Real, penalty: Real, filter: Length, max_iterations: Int, threshold: Real)` | no parameter, only `expect` entries | 1 case: the loop of topology optimization on the cells of `region` |

- Each entry names a `param` [`study_param`]. Each value has the type of the
  parameter [`type_mismatch`] and is in its range [`param_range`].
- A study does not change the steps. Each case runs all the steps.
- A model may hold many studies. `skein run model.skn` runs the model once,
  with the defaults. `skein run model.skn --study re_sweep` runs the study.
- `objective` names a scalar `report` of the last step
  [`study_objective`]. `nelder_mead` and `bobyqa` are gradient-free.
- `method: gradient` gets the derivative of the objective by each
  parameter from one adjoint solve in each case. It moves the values z in
  0..1 of the ranges by a projected gradient: the first step moves the
  parameter with the largest derivative by `step` of its range, each next
  step is that of Barzilai and Borwein, and a step that does not lower the
  objective (the Armijo rule) is halved. A parameter can change E, nu, the
  thickness, a load and the objective, but not the mesh, a region or a
  fixed value. `skein run` runs it for the models of `topology` below, with
  a range and a default in one unit. The defaults: `max_cases: 50`, `xtol:
  0.001`, `step: 0.1`. Reason: an objective of 4 parameters needs 8 more
  solves for each central difference, and the adjoint needs 1.
- `xtol` and `step` are fractions of the range of each parameter, above 0
  and at most 1 [`procedure_args`]. The first simplex of `nelder_mead` is
  the start and, for each parameter, the start moved by `step` of its
  range. The study stops when each point of the simplex is within `xtol`
  of the range from the best point, in each parameter, or after
  `max_cases` cases. Reason: femahoy gap G25. The case shape-hole-biaxial
  stops at 0.001 from a simplex of 0.1, and its tolerance comes from that
  stop rule.
- `topology` gives each cell of the mesh a density x in 0..1, and the
  modulus E (1e-9 + (1 - 1e-9) x^penalty) (SIMP). From x = `volume` in
  each cell, each iteration solves the step, gets the derivative of the
  objective by each density from one adjoint solve, filters it over the
  radius `filter`, and moves x by the optimality criteria update with the
  mean of x equal to `volume`. The loop stops when no density changes by
  more than 0.01, or after `max_iterations` (default 200). The default
  penalty is 3. With `threshold`, a cell at the threshold or more becomes
  solid (x = 1) and the others void (x = 0), and the step solves again: the
  reports and expects then give the solid. The objective is a `report` of
  the one static step that is the compliance of its one point load: a
  u_c(p) + b, with the load at p, the load component c only, and a of the
  sign of the load. `skein run` refuses another objective now. Reason: the
  optimality criteria update needs a derivative of the objective by each
  density of 0 or less.
  `objective`, `region`, `volume` and `filter` are necessary, `0 < volume
  <= 1`, `penalty >= 1`, `filter > 0`, `0 < threshold < 1`, and the study
  has no parameter entry [`procedure_args`]. `skein run` runs it for a 2d model with `cell: quad`
  on a `rect` and `linear_elastic`. The case records hold the table
  `topology` of each iteration, and the VTU holds `rho` in each cell.
  Reason: femahoy gap G1, the cases topology-mbb-beam and
  topology-mbb-beam-cad.
- An entry can be an `expect` (§15.3). The cases of a `grid` are in the
  order of the Cartesian product, with the last entry the fastest.

### 20.3 How cases run: Go-style concurrency

The study runner follows the model of Go: independent tasks that talk only
over typed channels.

```
             cases (typed queue)                 results (typed queue)
 scheduler ──────────────────────> worker 1 ─┐
           ──────────────────────> worker 2 ─┼──────────────────────> collector ──> table
           ──────────────────────> worker k ─┘
```

- The scheduler puts one **case record** for each case on a queue: the case
  number and the parameter values, as a flat value.
- Each worker takes a case, runs the whole model, and puts one **result
  record** on the result queue: the case number, the named reports, the
  `expect` results, the verification results, the run time and the status.
- The collector writes the records into one table, in case order, as they
  arrive.

On the native target, a worker is a Goose `thread_fn` started with
`thread_spawn_contained`, and the queues are Goose `qput` and `qget`.
Goose threads share no memory, so cases cannot interfere. A contained
worker that aborts ends only its own case. The study continues, and the
table gives that case the status `failed` and the abort message.

In the browser, a worker is a Web Worker that runs the WebAssembly module.
`postMessage` is the queue. The records are the same flat values.

`minimize` uses the same queues. Its scheduler puts a batch of cases, waits
for their results, and chooses the next batch.

### 20.4 Dependency analysis and reuse

A model has no mutable state, so the compiler knows the exact set of
declarations that each parameter reaches. It sorts the work of a case into
stages: geometry, mesh, assembly of the parts that do not depend on the
solution, and solve. A stage runs again only when a parameter that reaches
it changes.

| Parameter reaches | Each case does again |
|---|---|
| a material property or a load only | assembly and solve. The mesh is shared. |
| a size of the geometry or of the mesh | mesh, assembly and solve |
| only `expect` or `report` expressions | nothing but the post-processing |

A shared mesh is read-only. On the native target, each worker gets a copy
of it in its own data stacks. `skein check --study re_sweep` prints the
stages that each parameter reaches, so the user sees the cost before the
run.

### 20.5 Threads for a case and threads for a study

A machine has `H` hardware threads. A study runs at most `k` cases at a
time. The default `k` is `min(H, number of cases)`, and `--jobs k` sets
it. When a case starts and no case runs, it has `H` threads if no
second case can join it: there is no next case, `k` is 1, or the ledger
(below) does not admit the next case beside it. Such a case runs alone.
Each other case has 1 thread. Reason: on sample 01, 16 threads took 25%
more memory than 1 thread and the same time, so parallel cases scale
better than parallel assembly in one case.

**Memory.** Each solve checks its memory after the sparsity pattern of
its matrix and before the values (`src/solvers/memory.goose`). The
estimate is the resident memory at the check, plus a measured number of
bytes for each entry of the pattern and for each dof for its solution
path, plus the vectors of its solver. On more than 1 thread, it adds two
copies of the mesh for each thread: each worker of the assembly holds its
spawn argument and its globals. The pattern is what the program really builds, so
the estimate needs no table of cell types or physics. If the estimate
on the threads of the case does not fit but the estimate on 1 thread
fits, the solve goes on with 1 thread. If the estimate on 1 thread
does not fit in the budget, the program stops before the assembly with
the estimate on 1 thread and the free memory. Each case also sets its
`oom_score_adj` to 1000, so that the kernel kills a case first and not
another program.

The budget is 0.9 × min(MemAvailable, memory.max − memory.current of
each cgroup of the process), read once. Linux only: on another OS there
is no budget, and only `k` limits a study.

The study runner keeps a ledger: a case starts only if its estimate fits
in the budget less the estimates of the running cases, or if no case
runs. Reason: a case reaches its peak late, at the setup of its solver,
so the free memory at its start does not show what the running cases
take later. Before the first case, the runner probes the cases with the
fewest and the most dofs: each program prints the estimate of its first
solve on 1 thread and stops. The estimate of each case is the larger
bytes for each dof of the two probes times its dofs. Reason: the bytes
for each dof go down with the dofs on some paths (the RSS before the
solve, the GMRES restart of a wave solve) and up on others (the skyline
profile of a shell). When the model has more than one step, a later step can take
more memory than the first solve: then one case runs at a time until a
case gives its measured peak. Each case writes its measured peak
(`peak_mb`, from VmHWM, the `getrusage` peak of Linux) in its records,
and a measured peak for each dof that is larger than the estimate then
takes its place. When case 1 writes the stages that the other cases read, it runs
first and alone on `H` threads, and its peak is the first estimate.

**The task pool.** The thread of the runner starts the cases in case
order, each when the ledger admits it, and starts the next case when a
case ends. Thus a slow case does not stop the start of the other cases.

Goose threads share no memory, and Goose has no atomics or locks. Thus the
thread of the runner owns the queue and the ledger, and a worker gets its
case as the argument of its spawn. The table is in case order, so the
order of the runs does not change the results.

On the native target, a case is a separate program. Thus a task of the
pool is one whole case.

**The threads in a case.** Each case program has a second pool of
`SKEIN_THREADS` workers for the work in the case. Without
`SKEIN_THREADS`, it uses all hardware threads. Each worker has its own
copy of the data that it uses:

- Assembly: a task is a block of 16 cells. A free worker takes the next
  task from one queue, and gives the cell matrices of its cells. A worker
  copies the mesh when it starts.
- A CG solve, and a GMRES(m) solve of flow: each worker holds a block of
  the rows of the matrix, so the workers together hold one more copy of
  the matrix. For each product of the matrix and a vector, the thread of
  the case sends the vector to each worker, and each worker gives the
  product of its rows. The preconditioners and the other solvers run on
  the thread of the case.

The thread of the case adds the cell matrices in cell order. Each row of
a product is the same sum on each number of threads. Thus the results do
not change with the thread count or with the order of the messages, to
the last bit. A test compares a run on 1 thread with a run on 4 threads.

The product of the matrix and a vector reads all of the matrix, so the
memory bandwidth limits the gain of more threads in a solve. The shell
models (`src/shells`) use a direct skyline solve, which runs on one
thread.

### 20.6 The results table

The study writes `.skein/studies/<model>/<study>.toon` and prints the same
table:

```toon
study: re_sweep
jobs: 3
cases[3]{case,Re,status,verify,cd,time_s}:
  1,10,ok,pass,2.81,4.2
  2,20,ok,pass,2.05,4.6
  3,40,ok,pass,1.54,5.1
```

`jobs` is the most cases that ran at the same time (§20.5). One
column for each varied parameter, each named report and each `expect`.
`verify` is `pass` when each verification check of the case passes.
`--format json` gives the same records as JSON.

---

## 21. Dimensional analysis and scaling

Strict dimensions stay in force. This section adds the tools that engineers
use with dimensionless numbers: named groups, scaled equations and the
Buckingham π theorem. They work on the dimensional model and do not replace
it.

### 21.1 Dimensionless groups

A dimensionless group is an expression of type `Real`. A model states it as
a `param` or a `let`:

```skein
param Re: Real = 100
let mu: DynViscosity = rho*U*D/Re
```

This is the usual pattern for flow: the study varies `Re`, and the model
computes the dimensional value that gives it. The checker confirms that
`rho*U*D/Re` has the dimension `DynViscosity`.

### 21.2 `scales`

```skein
scales D, U, rho
```

Each model with dimensions has one `scales` line [`scales_missing`]. A
`nondim` model has none. Reason: an LLM writes a scales line as
easily as a unit system name, and the line makes the size of each unknown
explicit.

The `scales` line names reference values. Their dimensions must be
independent [`scales_dependent`], and together they must make each
dimension in the model [`scales_incomplete`]. The example spans length (D),
time (D/U) and mass (rho·D³), which is enough for incompressible flow
without heat.

With `scales`, the compiler writes each equation of the SIR in scaled
variables: each quantity `q` becomes `q* = q / q_ref`, where `q_ref` is the
product of powers of the scales with the dimension of `q`. The solver works
on the scaled system. Each input is dimensional in the source, and each
output, report and `expect` is dimensional again after the solve.

The scaled equations show the groups directly. For the cavity of sample 02
with `scales a, U, rho`, the momentum equation becomes
`ddt(v*) + grad(v*)*v* + grad(p*) − (1/Re)*lap(v*) = 0`, and the
`check --groups` report prints `Re = rho*U*a/mu`.

Reason for scales from the model, not from the unit system: the base units
are often far from the sizes of the problem, for example `m` and `s` for a
micro-channel, and scaled variables near 1 make the solver tolerances
meaningful.

### 21.3 The π groups of a model

`skein check --groups` applies the Buckingham π theorem to the model. Its
steps:

1. Take the constant inputs: each `param`, each constant `let` of a
   quantity, each material property, each mesh size.
2. Build the dimension matrix: one column for each input, one row for each
   base dimension.
3. Compute an integer basis of the null space. Each basis vector is a group:
   a product of powers of inputs with no dimension. Prefer the `scales` as
   the repeated variables.
4. Name the groups that match a known form.

```toon
groups[3]{name,expr,value}:
  Re,"rho*U*a/mu",100
  h_a,"m.size/a",0.0156
  h_lid,"m.refine(lid)/a",0.0078
```

The known forms:

| Name | Form | | Name | Form |
|---|---|---|---|---|
| `Re` | ρUL/μ | | `Pe` | UL/α |
| `Pr` | μc_p/k | | `Nu` | hL/k |
| `Ma` | U/c | | `Gr` | gβΔT L³/ν² |
| `Ra` | Gr·Pr | | `Fr` | U/√(gL) |
| `St` | fL/U | | `We` | ρU²L/σ |
| `Bi` | hL/k_solid | | `Fo` | αt/L² |
| `Da` | k_p/L² | | `Ec` | U²/(c_p ΔT) |

The groups drive two things: the lint rules that need a Reynolds or
Peclet number (`stabilization_missing`), and the report that a person
reads to judge the regime of the model.

### 21.4 What stays dimensional

- Source text: each quantity has its dimension and each literal its unit.
  `scales` does not change the type of any declaration.
- Boundary conditions, materials and outputs: dimensional.
- Only the SIR and the generated code use scaled variables.

Reason: dimensions find mistakes in the source, where the user writes. The
scaled form helps the solver, where the user does not write.

---

## 22. Export

```
skein export samples/01_cantilever.skn --geometry beam --format stl
```

| Format | Contents |
|---|---|
| `stl` | a triangulated surface from the mesher, at the size of `--size` |
| `step` | AP242 with tessellated geometry from the mesher |

Skein does not export exact surfaces (a B-rep). Reason: a B-rep kernel is a
project of its own size, and Skein is a simulation language first. See the
[tradeoffs](tradeoffs.md#4-parametric-parts-text-to-cad-and-export). Skein exports geometry. It still does not import meshes or CAD files
(decision 8).

---

## Appendix A. Built-in names

### Functions

| Function | Signature (informal) |
|---|---|
| `sym(A)`, `skew(A)`, `tr(A)`, `det(A)`, `inv(A)`, `transpose(A)`, `dev(A)` | `Tensor2` algebra |
| `dot(v, w)`, `cross(v, w)`, `outer(v, w)`, `norm(v)`, `ddot(A, B)` | products |
| `strain(u)` | strain under the strain assumption (§11.1) |
| `defgrad(u)` | deformation gradient `I + grad(u)` |
| `hooke(E, nu, eps)` | isotropic linear elastic stress |
| `elasticity(E, nu)` | the isotropic elasticity tensor, a `Tensor4` |
| `interp(xs, ys, x)` | piecewise linear function from two table columns |
| `sum(expr)` | the sum of `expr` over the rows of a table (§5.6) |
| `von_mises(sigma)`, `principal(sigma)` | stress measures |
| `fourier(k, T)` | heat flux `−k grad(T)` |
| `abs`, `sqrt`, `exp`, `ln`, `sin`, `cos`, `tan`, `tanh`, `atan`, `atan2(y, x)`, `min`, `max`, `pi` | scalar math |
| `re(c)`, `imag(c)`, `abs(c)`, `arg(c)`, `conj(c)` | parts of a `Complex<D>`: `re`, `imag` and `abs` give `D`, `arg` gives `Angle`, `conj` gives `Complex<D>` |
| `diff(W, A)` | symbolic derivative, for example `S = 2*diff(W, C)` from an energy |
| `tangent(f, top)`, `ply_stiffness(E1, E2, G12, nu12)`, `laminate(Q, layup, mirror, t_ply, h, axis, top, k)`, `along(f)` | shells and beams (§4.3) |
| `reaction(r)`, `integral(f, r)`, `mean(f, r)`, `argmax(f, r)`, `eigenvalue(k)`, `when(f, value)`, `rms(f, psd)`, `far_field(E, H, s, d)` | reports (§14.5), and `integral` and `mean` in the `eq` of a global unknown (§10.5) |
| `rigid(p)`, `moment(F*n)` | a rigid face and the moment on its point (§12.1) |
| `radiation(eps)` | the gray diffuse radiation of a face to the other faces of its enclosure, in a natural condition (§12.1) |
| `wall_distance()` | the distance to the nearest wall, a boundary with `v = 0`, in a turbulence model (§4.2) |
| `results(file, step:)` | §15.1 |
| `boundary(g, s)`, geometry functions, `copy(r)` | §8 |

### Constants

`I` (identity), `pi`, `x`, `y`, `z`, `t`, `n`, `tau` (the unit tangent of a `Curve`), `im` (imaginary unit), `omega` (the angular frequency of the step, §14.2).

### Standard library

The physics templates of §4.2 and §4.3. `skein help std` prints their source.

---

## Appendix B. References

- Alliez, P., Cohen-Steiner, D., Yvinec, M., Desbrun, M. (2005). Variational tetrahedral meshing. ACM Trans. Graph. 24(3).
- Arnold, D. N. (1982). An interior penalty finite element method with discontinuous elements. SIAM J. Numer. Anal. 19(4).
- Arnold, D. N., Brezzi, F., Cockburn, B., Marini, L. D. (2002). Unified analysis of discontinuous Galerkin methods for elliptic problems. SIAM J. Numer. Anal. 39(5).
- Barsoum, R. S. (1976). On the use of isoparametric finite elements in linear fracture mechanics. IJNME 10(1).
- Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. CMAME 197.
- Bérenger, J.-P. (1994). A perfectly matched layer for the absorption of electromagnetic waves. J. Comput. Phys. 114(2).
- Benito, J. J., Ureña, F., Gavete, L. (2001). Influence of several factors in the generalized finite difference method. Appl. Math. Model. 25.
- Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions. SIGGRAPH sketches.
- Chen, J. S., Pan, C., Wu, C. T., Liu, W. K. (1996). Reproducing kernel particle methods for large deformation analysis of non-linear structures. CMAME 139.
- Chen, J. S., Wu, C. T., Yoon, S., You, Y. (2001). A stabilized conforming nodal integration for Galerkin mesh-free methods. IJNME 50.
- Chen, L., Xu, J. (2004). Optimal Delaunay triangulations. J. Comput. Math. 22(2).
- Cheng, S.-W., Dey, T. K., Ramos, E. A. (2007). Delaunay refinement for piecewise smooth complexes. SODA.
- Chew, L. P. (1989). Guaranteed-quality triangular meshes. Cornell TR 89-983.
- Chung, J., Hulbert, G. M. (1993). A time integration algorithm for structural dynamics with improved numerical dissipation: the generalized-α method. J. Appl. Mech. 60.
- Dvorkin, E. N., Bathe, K. J. (1984). A continuum mechanics based four-node shell element for general non-linear analysis. Engineering Computations 1.
- Edmonds, J. (1965). Paths, trees, and flowers. Canadian J. Math. 17.
- Ghia, U., Ghia, K. N., Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. J. Comput. Phys. 48.
- Hughes, T. J. R., Franca, L. P., Balestra, M. (1986). A new finite element formulation for computational fluid dynamics: V. CMAME 59.
- Jones, R. M. (1999). Mechanics of Composite Materials, 2nd ed. Taylor & Francis.
- Knyazev, A. V. (2001). Toward the optimal preconditioned eigensolver: LOBPCG. SIAM J. Sci. Comput. 23(2).
- Labelle, F., Shewchuk, J. R. (2007). Isosurface stuffing: fast tetrahedral meshes with good dihedral angles. ACM Trans. Graph. 26(3).
- Liu, W. K., Jun, S., Zhang, Y. F. (1995). Reproducing kernel particle methods. IJNMF 20.
- Menter, F. R. (1994). Two-equation eddy-viscosity turbulence models for engineering applications. AIAA J. 32(8).
- Moser, R. D., Kim, J., Mansour, N. N. (1999). Direct numerical simulation of turbulent channel flow up to Re_tau = 590. Phys. Fluids 11(4).
- Remacle, J.-F., Lambrechts, J., Seny, B., Marchandise, E., Johnen, A., Geuzaine, C. (2012). Blossom-Quad: a non-uniform quadrilateral mesh generator using a minimum-cost perfect-matching algorithm. IJNME 89(9).
- Rhie, C. M., Chow, W. L. (1983). Numerical study of the turbulent flow past an airfoil with trailing edge separation. AIAA J. 21(11).
- Ruppert, J. (1995). A Delaunay refinement algorithm for quality 2-dimensional mesh generation. J. Algorithms 18(3).
- Schneiders, R. (1996). A grid-based algorithm for the generation of hexahedral element meshes. Engineering with Computers 12.
- Shakib, F., Hughes, T. J. R., Johan, Z. (1991). A new finite element formulation for computational fluid dynamics: X. CMAME 89.
- Shewchuk, J. R. (2002). Delaunay refinement algorithms for triangular mesh generation. Comput. Geom. 22.
- Shih, C. F., de Lorenzi, H. G., German, M. D. (1976). Crack extension modeling with singular quadratic isoparametric elements. Int. J. Fract. 12(4).
- Shu, C.-W., Osher, S. (1988). Efficient implementation of essentially non-oscillatory shock-capturing schemes. J. Comput. Phys. 77.
- Simo, J. C., Hughes, T. J. R. (1998). Computational Inelasticity. Springer.
- Tezduyar, T. E. (1992). Stabilized finite element formulations for incompressible flow computations. Adv. Appl. Mech. 28.
- Vaněk, P., Mandel, J., Brezina, M. (1996). Algebraic multigrid by smoothed aggregation for second and fourth order elliptic problems. Computing 56.
- Yuksel, C., Schaefer, S., Keyser, J. (2011). Parameterization and applications of Catmull-Rom curves. Computer-Aided Design 43(7).
