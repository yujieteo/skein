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
| Header: name, space dimension, unit system | `model`, `space`, `units` | §3 |
| Constants, geometry, regions, results of other models | `let` | §5, §8, §15 |
| Tables of data | `table` | §5.6 |
| Materials | `material` | §10.2 |
| Discretizations of a region | `mesh`, `points` | §9 |
| Unknown fields | `field` | §10.1 |
| Equations in strong form, with assumptions | `physics` | §10.3, §11 |
| Numerical method for each physics | `method` | §13 |
| Boundary conditions and loads | `bc` | §12 |
| The analysis sequence | `step` | §14 |
| Checks on the results | `expect` | §15.3 |
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
   pointers, no allocation and no mutable variables. The compiler makes all
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

A number has no sign. A minus sign before a number is the unary operator.

### 2.5 Keywords

```
adapt     allow     assume    at        bc        birth     contact   death
eq        every     expect    false     field     in        interface layers
let       material  mesh      method    model     on        output    physics
plastic   points    refine    release   report    solve     solves    space
state     step      table     true      units     where     within
```

A keyword is not an identifier. The names of types, units, built-in
functions and assumptions are not keywords. A declaration can hide them,
but the checker reports [`shadows_builtin`] as a lint.

### 2.6 Unit names

After a number literal, after a vector literal and after the keyword `in`,
the lexer reads a **unit expression** (§6.3). In those places an identifier
is a unit name, not a variable. Thus `L` in `let L: Length = 2 L` is a
variable in the first place and the litre in the second place.

---

## 3. Grammar

The grammar is in EBNF. `{ x }` is zero or more, `[ x ]` is optional,
`( a | b )` is a choice. `NL` is a statement end (§2.2). Terminals are in
quotes.

```ebnf
program     = header { item } ;
header      = "model" ident NL "space" ("1d" | "2d" | "3d" | "axi") NL
              "units" unitsys NL ;
unitsys     = ident | "{" unitfield { "," unitfield } "}" ;
unitfield   = ident ":" unitexpr ;

item        = letdecl | tabledecl | material | meshdecl | pointsdecl
            | fielddecl | physics | method | bcdecl | contact | interface
            | step | expect | allow ;

letdecl     = "let" ident ":" type "=" expr NL ;
tabledecl   = "table" ident "{" column { "," column } "}" ":" NL
              { INDENT row NL } ;
column      = ident ":" type [ "in" unitexpr ] ;
row         = expr { "," expr } ;

material    = "material" ident "{" NL { prop NL } "}" NL ;
prop        = ident ":" type "=" expr ;

block       = "{" [ NL ] [ entry { ( "," | NL ) entry } [ NL ] ] "}" ;
setting     = ident ":" expr ;

meshdecl    = "mesh" ident "on" region block NL ;    (* entry = setting | refine | layers *)
pointsdecl  = "points" ident "on" region block NL ;  (* entry = setting | refine *)
refine      = "refine" region ":" expr ;
layers      = "layers" region ":" expr ;

fielddecl   = "field" ident ":" type [ "on" region ] [ "=" expr ] NL ;

physics     = "physics" ident [ "on" region ] "solves" ident { "," ident }
              "{" NL { physstmt NL } "}" NL ;
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

step        = "step" ident ":" procedure "{" NL { stepstmt NL } "}" NL ;
procedure   = ident [ "(" [ namedargs ] ")" ] ;
stepstmt    = bcdecl | "release" region
            | "birth" region { "," region } [ "where" expr ]
            | "death" region { "," region } [ "where" expr ]
            | "solve" block | "adapt" block
            | "output" expr { "," expr } [ "every" expr ]
            | "report" expr { "," expr } ;

expect      = "expect" ident ":" addexpr "==" addexpr "within" expr NL ;
allow       = "allow" ident ":" string NL ;

region      = expr ;

type        = qtype | "Vector" "<" qtype ">" | "Tensor2" "<" qtype ">"
            | "Tensor4" "<" qtype ">" | ident ;
qtype       = dimterm { ("*" | "/") dimterm } ;
dimterm     = ident [ "^" ["-"] digits ] | "(" qtype ")" ;

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
args        = arg { "," arg } ;
arg         = expr | ident ":" expr ;
namedargs   = ident ":" expr { "," ident ":" expr } ;

unitexpr    = unitfactor { ("*" | "/") unitfactor } ;
unitfactor  = unitname [ "^" ["-"] digits ] | "(" unitexpr ")" ;
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
- `1d`, `2d`, `3d` and `axi` are tokens only after `space`. Elsewhere `3d`
  is a syntax error.
- `space axi` is the axisymmetric space (§11.3). It is not an assumption,
  because it changes the coordinates of every geometry function.
- The grammar admits more than the type rules accept. For example, a
  `region` is any `expr`. §5 and §8 give the rules.
- `|` and `&` are geometry union and intersection. `||` and `&&` are the
  Boolean operators. The two pairs have the same precedence.

### 3.1 Order of items

The header comes first, in the order `model`, `space`, `units`. The items
after it may come in any order, but a name must be declared before its first
use [`use_before_decl`]. Steps run in the order of the file.

Reason: declaration before use lets the checker resolve names in one pass,
and it puts each error at the line of the use.

---

## 4. Names and scope

A model has **one namespace**. These declarations add a name to it: `let`,
`table`, `material`, `mesh`, `points`, `field`, `physics`, the `let`
statements inside a `physics` block, `step` and `expect`. Two declarations
with one name are an error [`duplicate_name`].

Reason: one namespace lets a step write `sigma*n` for the stress that a
physics block declares, with no prefix.

These names have a smaller scope:

| Name | Visible in |
|---|---|
| A material property, for example `E` | the physics blocks that assign the material, and the `output`, `report` and `bc` expressions that use those physics |
| A column of a table, `t.col` | everywhere after the table |
| A field of a `Results` value, `r.u` | everywhere after the `let` |
| `x`, `y`, `z`, `t` (coordinates and time) | every expression that the solver evaluates at a point |
| `n` (outward unit normal) | `bc` expressions on a boundary |
| `I` (identity tensor) | everywhere |

When two materials on the regions of one physics give a property, the
property must have the same type in both [`material_prop_type`]. When a
physics uses a property, each material that the physics assigns must give
it [`material_prop_missing`].

---

## 5. Types

Skein has no type inference. Each `let`, `field`, material property and
table column states its type. The checker compares the declared type with
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
| `Str` | a string, for file paths |

### 5.2 Tensor types

| Type | Rank | Components |
|---|---|---|
| `D` (a dimension) | 0 | 1 |
| `Vector<D>` | 1 | `d` |
| `Tensor2<D>` | 2 | `d × d` |
| `Tensor4<D>` | 4 | `d × d × d × d` |

`d` is the space dimension from the header (`space 2d` gives 2, `space axi`
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

A `field` declaration gives the tensor type and the region of an unknown. A
field without `on` has no region: it is a global unknown, for example the
state of an ODE (§10.5).

### 5.6 Tables

A table has typed columns and TOON rows. `in unit` after a column type gives
the unit of each bare number in that column. A cell may also carry its own
unit. Each row must have one value for each column [`table_row_width`].

`t.col` is the column as a list. A list may stand where an expression of
the column type is expected in `expect`, `output`, `report` and in the
point of an `at`. The statement then applies to each row (§15.3).

### 5.7 Assignability

A value of type `A` may initialize a declaration of type `B` only when `A`
and `B` are the same type. Skein has no implicit conversion. These are the
only exceptions:

1. The literal `0` (with no unit) has every dimension and every tensor rank.
   It means zero of the expected type.
2. A bare number has type `Real`. `Angle` accepts `Real`.
3. `Int` converts to `Real` in arithmetic.

Reason for exception 1: `u = 0` is the most common boundary condition, and
`u = (0, 0, 0) mm` adds tokens and gives no safety.

---

## 6. Dimensions and units

### 6.1 Dimensions

A dimension is a vector of seven integer exponents over the SI base
dimensions: mass M, length L, time T, temperature Θ, current I, amount N and
luminous intensity J. The checker does arithmetic on these vectors:

- `a*b` adds the exponents, `a/b` subtracts them, `a^k` multiplies them by
  `k`. `k` must be an integer constant.
- `a + b`, `a - b`, `a == b` and the other comparisons need equal
  dimensions [`unit_mismatch`].
- `sqrt(a)` halves the exponents. Each exponent must be even
  [`unit_fractional`].
- `exp`, `ln`, `sin`, `cos`, `tan`, `tanh` take `Real` or `Angle`
  [`unit_transcendental`].
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
coefficients in `K`.

### 6.4 Literals

- `210 GPa` is a scalar. Its dimension comes from the unit.
- `(40, 0, 20) mm` is a vector. The unit after the parenthesis applies to
  each bare number in it.
- `(L, W, H)` is a vector of expressions. All components must have the same
  dimension. A `0` component takes that dimension.
- A bare number is `Real`.

The value of a literal is stored in SI. The checker compares dimensions,
never units, so `1 m + 1 mm` is valid and is 1.001 m.

### 6.5 Unit systems

The `units` line names the unit system of the model:

| Name | Length | Mass | Time | Force | Stress | Temperature |
|---|---|---|---|---|---|---|
| `SI` | m | kg | s | N | Pa | K |
| `mm_N_MPa_s` | mm | t | s | N | MPa | K |
| `mm_kg_ms` | mm | kg | ms | kN | GPa | K |

A user system gives the base units: `units { length: mm, mass: t, time: s,
temperature: K }`.

The unit system does not change the meaning of the source, because each
literal carries its unit. It sets two things:

1. the scale of the numbers in the generated code, so that the solver works
   with numbers near 1 for the model;
2. the units of the output files and of the numbers that the tools print.

Reason: a declared system keeps the numbers of the solver well scaled, and
it lets a user read results in the units of the drawing.

### 6.6 What Skein leaves out

- **No type aliases.** Reason: an alias puts one fact in two places, and an
  LLM must then find the alias to know the dimension.
- **No unit-generic functions.** Users do not declare functions (§7.4).

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
| `curl(v)` | `Vector`, `space 3d` | `Vector` |
| `lap(f)` | rank r | rank r. Equal to `div(grad(f))` |
| `ddt(f)` | any | same rank. First time derivative |
| `ddt2(f)` | any | same rank. Second time derivative |

With these definitions `grad(v)*v` is the convective term `(v·∇)v`.

`ddt` and `ddt2` are zero in a `static` step. The `transient` procedure
replaces them with the time scheme of the step (§14.2).

### 7.4 Functions

Users do not declare functions. Reason: physics needs closed-form
constitutive laws, and the built-in functions (Appendix A) and `let`
declarations cover them. A function body would need loops or recursion, and
these would break rule 4 of §1.

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
| `cone(c, axis, r0, r1, h)` | `Solid` | |
| `torus(c, axis, R, r)` | `Solid` | |
| `rect(p0, p1)` | `Surface` | 2D only |
| `circle(c, r)` | `Surface` | 2D only |
| `polygon(p1, p2, …)` | `Surface` | 2D only, counter-clockwise or clockwise |
| `extrude(s, h)` | `Solid` | `s` is a 2D `Surface` in the x-y plane |
| `revolve(s, angle)` | `Solid` | about the y axis |
| `sweep(s, path)` | `Solid` | `path` is a `Curve` |
| `line(p0, p1)`, `arc(c, r, a0, a1)`, `spline(p…)` | `Curve` | |
| `point(p)` | `Point` | |
| `nurbs(…)` | `Surface` or `Solid` | phase 5 (roadmap) |

Operators and modifiers:

| Expression | Meaning | Distance function |
|---|---|---|
| `a | b` | union | `min(φa, φb)` |
| `a & b` | intersection | `max(φa, φb)` |
| `a - b` | difference | `max(φa, −φb)` |
| `fillet(a, r)` | round the convex and concave edges of `a` | smooth minimum with radius `r` |
| `chamfer(a, d)` | bevel the edges of `a` | chamfer minimum with distance `d` |
| `move(a, v)`, `rotate(a, axis, angle)`, `scale(a, k)` | rigid motions and uniform scale | |

The `min` and `max` of distance functions give a bound on the distance, not
the exact distance. The mesher uses the CSG tree to correct this near the
surface (§9.2).

### 8.2 Regions

A **region** is a geometry value that a `field`, `physics`, `mesh`,
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

A `Surface` in `space 3d`, or a `Curve` in `space 2d` or `space 3d`, can be
a region of its own physics, for example a shell or a beam. Such a region
must be the boundary of a solid region, or a geometry value that the
`mesh` block embeds. The mesher makes the cells of the entity conform to
the cells of the region around it. A `Point` can be the region of a `bc`,
for example a point load or a pressure reference.

### 8.4 Region algebra

Regions combine with `|`, `&` and `-`. The checker computes containment by
the CSG tree: `a` is inside `b` if `a & b` equals `a` in the tree, or if a
sampled test on the mesh shows it. A field, a physics, a `bc` or a
`birth`/`death` statement must name a region inside the region of the field
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
| `cell` | `tri`, `quad`, `tet`, `hex` | `tri` in 2D, `tet` in 3D | `quad` and `hex` need a `box` or `rect` region |
| `grade` | `Real` | 1.3 | largest ratio of the size of two neighbour cells |
| `feature_angle` | `Angle` | 30 deg | an edge with a dihedral angle smaller than 180 deg minus this angle is sharp |
| `refine r` | `Length` | none | size on region `r` and near it |
| `layers r` | record | none | boundary layer cells on `r` (CFD) |

The cell order comes from the `method` (§13), not from the mesh. The mesher
makes the nodes for that order on the true geometry.

### 9.2 Mesh generation

The algorithms come from published papers only. The implementation must
not copy code or text from Gmsh, Triangle, TetGen or CGAL, because their
licences are not compatible with Apache 2.0.

| Case | Algorithm |
|---|---|
| 2D triangles | Delaunay refinement (Chew 1989, Ruppert 1995, Shewchuk 2002) on a boundary that the mesher samples from φ, with the sharp corners as fixed vertices |
| 3D tetrahedra | isosurface stuffing (Labelle and Shewchuk 2007) on an octree graded by `size`, `refine` and `grade` |
| 3D sharp edges | protecting balls on the feature curves (Cheng, Dey and Ramos 2007), then the stuffing |
| smoothing | optimal Delaunay triangulation smoothing (Chen and Xu 2004, Alliez et al. 2005), with no move that makes a cell worse |
| `quad`, `hex` | a structured grid on a `box` or `rect` |
| boundary layers | prisms or quads extruded along the normal of the surface, then the core mesh |

The sharp edges come from the CSG tree: each intersection of the surfaces
of two primitives with a dihedral angle that `feature_angle` marks is a
feature curve.

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

---

## 10. Fields, materials and physics

### 10.1 `field`

```skein
field u: Vector<Length> on beam
```

A field is an unknown of the model. Its type is a tensor type (§5.2). Its
region is the set of points where it exists. The optional `= expr` is the
initial value (§14.1). It may be varying. A field must have exactly one
physics that solves it [`field_unsolved`], [`field_solved_twice`].

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
for example a function of a temperature field.

### 10.3 `physics`

```skein
physics solid on beam solves u {
    material steel
    assume small_strain, isotropic
    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))
    eq div(sigma) = 0
}
```

- `on r` is the region of the equations. It defaults to the region of the
  first field in `solves`.
- `solves u, p, …` lists the unknowns. The block must have one `eq` for each
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
physics heat on part solves T { … }
physics solid on part solves u {
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

A field without a region is a global unknown. A physics without `on` that
solves only global unknowns is a system of ODEs:

```skein
field q: Length
physics oscillator solves q {
    eq m*ddt2(q) + c*ddt(q) + k*q = 0
}
```

A global unknown may appear in a physics with a region. Then the system
couples an ODE and a PDE, for example a lumped spring on a solid.

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
the dimension of the stress [`unit_mismatch`]. A hyperelastic law needs no
state. Write the strain energy and derive the stress with `diff`:

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
setting `mesh_motion` gives the equation of the mesh displacement. FSI is in
phase 5 of the [roadmap](roadmap.md).

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
| `plane_stress(t)` | `space 2d` only. σzz = σxz = σyz = 0. The compiler solves the constitutive law for εzz and removes it (§11.2). Thickness `t` multiplies volume and surface integrals. | `plane_strain`, `space 3d` |
| `plane_strain(t)` | `space 2d` only. εzz = εxz = εyz = 0. `t` defaults to 1 of the length unit. | `plane_stress` |
| `incompressible` | the physics must have an equation `div(v) = 0`, or a pressure unknown that enforces `tr(strain(u)) = 0` | `compressible` |
| `compressible` | density is a varying value or an unknown | `incompressible` |
| `newtonian` | the fluid stress must be linear in `sym(grad(v))` | |
| `laminar` | no turbulence model | `turbulent(model)` |
| `turbulent(model)` | adds the closure `model` (phase 5) | `laminar` |
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

`space axi` uses the coordinates (r, z) and the names `x` for r and `y` for
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
| `F*n` or `dot(F, n)`, where `div(F)` is in an equation of the physics | natural (Neumann). The right side may use the unknowns (Robin). |
| anything else | an error [`bc_form`] |

A natural condition on the flux of the equation paired with unknown `u` is a
condition for `u`.

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
  [`bc_conflict`];
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
Newton's method. Contact is in phase 5 of the [roadmap](roadmap.md).

### 12.5 Values from other models

The right side of an essential condition may be a field of a `Results`
value (§15.1). The solver interpolates it at the nodes, or at the
integration points for a weak condition.

---

## 13. Methods

```skein
method solid: cg {
    mesh: m
    order: 2
}
```

`method p: kind { … }` gives the numerical method of physics `p`. Each
physics must have one method [`method_missing`]. The kinds:

| Kind | Settings | Notes |
|---|---|---|
| `cg` | `mesh`, `order` (Int or a record per field), `stabilize`, `quadrature`, `couple` | continuous Galerkin FEM |
| `dg` | `mesh`, `order`, `flux` (`upwind`, `lax_friedrichs`, `roe`), `penalty` | interior penalty for second-order terms |
| `fv` | `mesh`, `gradient` (`green_gauss`, `least_squares`), `flux` (`upwind`, `linear_upwind`, `central`), `pv_coupling` (`rhie_chow`) | cell-centred finite volume |
| `rkpm` | `points`, `basis` (`linear`, `quadratic`), `support` (Real, times the spacing), `integrate` (`scni`, `gauss`), `bc` (`nitsche`) | Galerkin meshless |
| `gfd` | `points`, `stencil` (Int), `weight` | strong-form collocation |
| `ode` | none | for a physics with only global unknowns (§10.5). The time scheme comes from the step. |

The setting `stabilize` takes `none`, `supg`, `supg_pspg`, `vms` or
`supg_pspg_lsic`. §17.2 gives the meaning.

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
of each field: the expression after `=` in the `field` declaration, or zero.

```skein
field T: Temperature on part = 293 K
```

### 14.2 Procedures

| Procedure | Arguments | Meaning |
|---|---|---|
| `static` | `ramp` (Int, default 1) | `ddt` = 0. Nonlinear problems take `ramp` load increments. |
| `transient` | `end`, `dt`, `scheme` | time stepping from the end time of the step before |
| `modal` | `count`, `shift` | the `count` lowest eigenpairs of the linearized operator about the current state |
| `buckling` | `count` | eigenpairs of the stiffness and the geometric stiffness of the current state |

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
- `where cond` selects the cells of `r` where `cond` is true at the start of
  the step. `cond` may use the fields.
- The mesh of the whole model, alive or dead, is made once, before the first
  step. Thus each region that a `birth` names must be a part of the region
  of the mesh.

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

### 14.5 `output` and `report`

`output e1, e2, …` writes the fields or field expressions to the result
file at the end of the step, or each `every k` time steps. `report e` prints
a scalar or a vector to the run log in TOON. Functions for `report`:
`reaction(r)`, `max(f, r)`, `min(f, r)`, `integral(f, r)`, `mean(f, r)`.

### 14.6 Mesh adaptivity

```skein
adapt { estimator: zz, target: 1%, max_levels: 3 }
```

`adapt` in a step solves, estimates the error in each cell, refines the
cells with the largest error, and solves again. It stops when the estimated
relative error in the energy norm is below `target`, or after `max_levels`.
`estimator` is `zz` (Zienkiewicz and Zhu 1987) or `residual` (from `R` in
the SIR). The state of the step before moves to the new mesh by
interpolation.

### 14.7 Restarts

The run writes the state of the model at the end of each step to
`.skein/state/<model>/<step>.sks`. `skein run model.skn --from <step>`
starts at step `<step>` from the state at the end of the step before it. The
source of the steps before `<step>` must not change
[`restart_stale`]: the state file holds a hash of their text.

---

## 15. Results, submodels and expectations

### 15.1 `results`

```skein
let global: Results = results("01_cantilever.skn", step: "load")
```

`results(file, step:)` refers to the state of another model at the end of
the step. The checker reads the other file and gives `global.u` the type of
the field `u` in it. The run uses the saved state, or runs the other model
first if no state exists.

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
for `skein verify convergence` ([command line](command_line.md)). The run evaluates it after the
last step and reports it as a post-solve verification result. `within`
takes an absolute value with the dimension of the sides, or a percent of
the right side. A table column in the expression makes one test for each
row. A failed expectation gives a non-zero exit status.

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
  lowering keeps a split `F = F_a(U) + F_d(U, ∇U)`: the advective part has
  no derivative of an unknown, the diffusive part has.
- `S` holds all other terms, for example the body force and the convective
  term `grad(v)*v`.
- `M` and `N` are the coefficients of the time derivatives.

The SIR also keeps the **full strong residual** `R(U)` of each equation.
Residual-based stabilization (§17.2) and residual error estimates need it.

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
2. Integrate `div F_d` by parts: `∫ w div F_d = −∫ grad(w) : F_d + ∮ w F_d·n`.
3. Keep `div F_a` as it is. Reason: an advective flux has no derivative of
   an unknown, and integration by parts brings no gain for `cg`.
4. Replace the boundary term `∮ w F_d·n` with the natural conditions. Drop
   it where an essential condition holds, because `w` is zero there.

The rule in step 2 is the one rule for all physics: the compiler integrates
by parts the `div` of an expression that holds a derivative of an unknown.
Thus `div(sigma)` with `sigma` from `grad(u)` goes by parts, and `div(v)` in
the continuity equation does not.

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
`τ_C = h²/(τ_M · tr(G))`, after Tezduyar (1992) and Shakib et al. (1991).
`R_M` is the residual of the momentum equation and `q` the pressure test
function. In a `static` step, the `2/Δt` term is zero.

### 17.3 Finite volume (`fv`)

The lowering integrates each balance law over each cell and applies the
divergence theorem to all `div` terms, advective and diffusive. The face
fluxes use the `flux` and `gradient` settings. Pressure-velocity coupling
uses Rhie-Chow interpolation (Rhie and Chow 1983) of the face velocity. The
system is solved fully coupled by Newton's method. The unknowns are cell
averages. Essential conditions act through the face flux.

### 17.4 Discontinuous Galerkin (`dg`)

The lowering applies §17.1 on each cell. The boundary terms on interior
faces use the numerical flux of the `flux` setting for `F_a` and the
symmetric interior penalty form (Arnold 1982, Arnold et al. 2002) for
`F_d`. Essential conditions are weak, by the same penalty.

### 17.5 Meshless

- `rkpm` uses the Galerkin form of §17.1 with reproducing kernel shape
  functions (Liu, Jun and Zhang 1995, Chen et al. 1996). The integration is
  stabilized conforming nodal integration (Chen et al. 2001). The shape
  functions do not interpolate, so essential conditions use Nitsche's method.
- `gfd` collocates the strong form at each point with generalized finite
  difference stencils (Benito, Ureña and Gavete 2001). It uses `R` from the
  SIR directly, with no weak form. The natural conditions are collocated at
  the boundary points.

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

Note: Goose on WebAssembly is a proposal in the Goose repository, not a
feature (`docs/design/wasm.md` in Goose). Each Goose data stack needs
reserved address space, and wasm32 has no reservation. The Skein browser
target thus depends on the "fixed partitions" plan of that proposal. Phase
1 does not need it. See the risks in the [roadmap](roadmap.md#risks).

### 18.2 Files

`skein build model.skn` writes:

| File | Contents |
|---|---|
| `model.goose` | the kernels and the plan of the model |
| `skein_rt` (package) | the runtime: mesher, sparse matrices, solvers, output writers. The same for each model. |

### 18.3 The shape of the generated code

The plan stores a mesh as flat arrays, which Goose keeps in data stacks with
no pointers. A cell kernel is a `simd fn`. Assembly runs in worker threads
over a colouring of the cells, so no two workers write one row.

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
    let gn = skein_rt::tet4_gradients(x);          // double3[4], and the volume
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
        let job = qget<skein_rt::ColourJob>();
        guard job.first >= 0 else { break; }
        qput(skein_rt::assemble_range(job, k_solid_tet4));
    }
}
```

### 18.4 The runtime

`skein_rt` is a Goose package. It holds no C, except the C that the Goose
standard library uses for files and the clock.

| Part | Contents | Phase |
|---|---|---|
| sparse | CSR matrices, colouring, assembly | 2 |
| linear | CG, GMRES(m), Jacobi, ILU(0) | 2 |
| nonlinear | Newton with line search, load ramps | 2 |
| eigen | LOBPCG (Knyazev 2001) for `modal` and `buckling` | 3 |
| multigrid | smoothed aggregation AMG (Vaněk, Mandel and Brezina 1996) | 3 |
| mesher | §9.2 | 2 |
| output | VTU writer, TOON reports, state files | 2 |

### 18.5 The compiler in Goose

The Skein compiler is a Goose program too. Thus `skein` builds for the
browser from the same source as the native tool, and the editor in the
browser runs `check`, `fmt` and `lint` with no server.

---

## 19. Static semantics checklist

`skein check` applies these checks in this order. Each failed check gives
the diagnostic in brackets. The full list with the fixes is in
[checks.md](checks.md).

1. Lexical and syntax rules (§2, §3).
2. Names: declared before use, no duplicates (§3.1, §4).
3. Types and dimensions of each declaration and expression (§5, §6, §7).
4. Geometry dimensions and region containment (§8).
5. Fields, materials and physics: each unknown solved once, one equation
   per unknown, materials cover the region (§10).
6. Assumptions: no conflicts, conditions met (§11).
7. Boundary conditions: form, conflicts, location (§12).
8. Methods: present, settings valid, compatible with the physics (§13).
9. Steps: procedures valid, constant arguments, birth and death regions
   (§14).
10. Results and expectations: files exist, types match (§15).
11. `allow` names a check that exists, with a reason that is not empty
    [`allow_reason`].

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
| `von_mises(sigma)`, `tresca(sigma)`, `principal(sigma)` | stress measures |
| `fourier(k, T)` | heat flux `−k grad(T)` |
| `abs`, `sqrt`, `exp`, `ln`, `sin`, `cos`, `tan`, `tanh`, `min`, `max`, `pi` | scalar math |
| `diff(W, A)` | symbolic derivative, for example `S = 2*diff(W, C)` from an energy |
| `reaction(r)`, `integral(f, r)`, `mean(f, r)` | reports |
| `results(file, step:)` | §15.1 |
| `boundary(g, s)`, geometry functions | §8 |

### Constants

`I` (identity), `pi`, `x`, `y`, `z`, `t`, `n`.

---

## Appendix B. References

- Alliez, P., Cohen-Steiner, D., Yvinec, M., Desbrun, M. (2005). Variational tetrahedral meshing. ACM Trans. Graph. 24(3).
- Arnold, D. N. (1982). An interior penalty finite element method with discontinuous elements. SIAM J. Numer. Anal. 19(4).
- Arnold, D. N., Brezzi, F., Cockburn, B., Marini, L. D. (2002). Unified analysis of discontinuous Galerkin methods for elliptic problems. SIAM J. Numer. Anal. 39(5).
- Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. CMAME 197.
- Benito, J. J., Ureña, F., Gavete, L. (2001). Influence of several factors in the generalized finite difference method. Appl. Math. Model. 25.
- Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions. SIGGRAPH sketches.
- Chen, J. S., Pan, C., Wu, C. T., Liu, W. K. (1996). Reproducing kernel particle methods for large deformation analysis of non-linear structures. CMAME 139.
- Chen, J. S., Wu, C. T., Yoon, S., You, Y. (2001). A stabilized conforming nodal integration for Galerkin mesh-free methods. IJNME 50.
- Chen, L., Xu, J. (2004). Optimal Delaunay triangulations. J. Comput. Math. 22(2).
- Cheng, S.-W., Dey, T. K., Ramos, E. A. (2007). Delaunay refinement for piecewise smooth complexes. SODA.
- Chew, L. P. (1989). Guaranteed-quality triangular meshes. Cornell TR 89-983.
- Chung, J., Hulbert, G. M. (1993). A time integration algorithm for structural dynamics with improved numerical dissipation: the generalized-α method. J. Appl. Mech. 60.
- Ghia, U., Ghia, K. N., Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. J. Comput. Phys. 48.
- Hughes, T. J. R., Franca, L. P., Balestra, M. (1986). A new finite element formulation for computational fluid dynamics: V. CMAME 59.
- Knyazev, A. V. (2001). Toward the optimal preconditioned eigensolver: LOBPCG. SIAM J. Sci. Comput. 23(2).
- Labelle, F., Shewchuk, J. R. (2007). Isosurface stuffing: fast tetrahedral meshes with good dihedral angles. ACM Trans. Graph. 26(3).
- Liu, W. K., Jun, S., Zhang, Y. F. (1995). Reproducing kernel particle methods. IJNMF 20.
- Rhie, C. M., Chow, W. L. (1983). Numerical study of the turbulent flow past an airfoil with trailing edge separation. AIAA J. 21(11).
- Ruppert, J. (1995). A Delaunay refinement algorithm for quality 2-dimensional mesh generation. J. Algorithms 18(3).
- Shakib, F., Hughes, T. J. R., Johan, Z. (1991). A new finite element formulation for computational fluid dynamics: X. CMAME 89.
- Shewchuk, J. R. (2002). Delaunay refinement algorithms for triangular mesh generation. Comput. Geom. 22.
- Shu, C.-W., Osher, S. (1988). Efficient implementation of essentially non-oscillatory shock-capturing schemes. J. Comput. Phys. 77.
- Simo, J. C., Hughes, T. J. R. (1998). Computational Inelasticity. Springer.
- Tezduyar, T. E. (1992). Stabilized finite element formulations for incompressible flow computations. Adv. Appl. Mech. 28.
- Vaněk, P., Mandel, J., Brezina, M. (1996). Algebraic multigrid by smoothed aggregation for second and fourth order elliptic problems. Computing 56.
- Zienkiewicz, O. C., Zhu, J. Z. (1987). A simple error estimator and adaptive procedure for practical engineering analysis. IJNME 24.
