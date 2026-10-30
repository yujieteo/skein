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
- [22. Manufacture](#22-manufacture)
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
| Symmetry groups: cyclic, dihedral, periodic | `symmetry` | §12.6 |
| The analysis sequence | `step` | §14 |
| Checks on the results | `expect` | §15.3 |
| Reference scales for the solver | `scales` | §21 |
| Parameter studies | `study` | §20 |
| Manufacturing process and tolerances | `manufacture`, `tolerance` | §22 |
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
adapt       allow       assume      at          bc          birth
contact     death       eq          every       expect      false
fn          in          interface   layers      let         manufacture
material    mesh        method      model       on          output
param       physics     plastic     points      refine      release
report      scales      solve       solves      state       step
study       symmetry    table       template    tol         tolerance
true        use         where       within
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
            | interface | symmetry | step | expect | allow | scales | study
            | manufacture | tolerance ;

usedecl     = "use" string NL ;
fndecl      = "fn" ident "(" fnparam { "," fnparam } [ "," ] ")" "->" type "=" expr NL ;
fnparam     = ident ":" type ;

paramdecl   = "param" ident ":" type "=" expr [ "tol" expr ] [ "in" expr ".." expr ] NL ;
manufacture = "manufacture" region ":" procedure block NL ;  (* entry = setting *)
tolerance   = "tolerance" ident ":" expr "in" expr ".." expr [ block ] NL ;
letdecl     = "let" ident ":" type "=" expr NL ;
scales      = "scales" expr { "," expr } NL ;
study       = "study" ident ":" procedure block NL ;  (* entry = ident ":" ( list | range ) | setting *)
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

meshdecl    = "mesh" ident "on" region block NL ;    (* entry = setting | refine | layers *)
pointsdecl  = "points" ident "on" region block NL ;  (* entry = setting | refine *)
refine      = "refine" region ":" expr ;
layers      = "layers" region ":" expr ;

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
            | "solve" block | "adapt" block
            | "output" expr { "," expr } [ "every" expr ]
            | "report" [ ident "=" ] expr { "," [ ident "=" ] expr } ;

expect      = "expect" ident ":" addexpr "==" addexpr "within" expr NL ;
allow       = "allow" ident ":" string NL ;

region      = expr ;

type        = qtype | "Vector" "<" stype ">" | "Tensor2" "<" stype ">"
            | "Tensor4" "<" stype ">" | numsys | ident ;
stype       = qtype | numsys ;
numsys      = ( "Complex" | "Quaternion" | "Dual" ) "<" qtype ">" ;
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
args        = arg { "," arg } [ "," ] ;
arg         = [ name ":" ] ( expr | list ) ;
namedargs   = name ":" expr { "," name ":" expr } ;

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
use [`use_before_decl`]. Steps run in the order of the file.

Reason: declaration before use lets the checker resolve names in one pass,
and it puts each error at the line of the use.

### 3.2 Layout for small diffs

A change to a model must give a small diff that a person can read. Skein
text and `skein fmt` follow these rules (command_line.md §3.1):

- One item, one statement or one block entry on each line. A table row is
  one line.
- No column alignment. `fmt` never adds spaces to line up values, so a
  change to one value does not change the lines near it. This holds for
  TOON tables too.
- A list that does not fit on one line has one item on each line, with a
  comma after each item. A new item then adds one line.
- `fmt` keeps the order of items and never reflows a line that did not
  change.

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
template [`type_mismatch`]. The `let` names of the template, here `sigma`,
enter the namespace of the model. The block of the model may add `material`,
`assume` and `let` statements. It may not add `eq` [`template_eq`]. Only
the standard library declares templates in this version.

| Template | Unknowns | Phase |
|---|---|---|
| `linear_elastic(b)` | `u` | 2 |
| `heat(q_v)` | `T` | 2 |
| `stokes(f)`, `navier_stokes(f)` | `v`, `p` | 3 |
| `reissner_mindlin(t)` (shell), `timoshenko(A, I_y, I_z)` (beam) | `w`, `theta` | 5 |
| `k_omega_sst` (turbulence) | `k`, `omega` | 5 |

Reason: nobody writes a shell or a k-ω SST model from memory, and a
template cuts the size of the common models. The source stays strong form,
so a user can read it and copy it into a model to change it.

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
| `Str` | a string, for file paths |

Each quantity can also take a value in one of three other number systems.
The list is closed: Skein has no user number systems.

| Type | Values | Use |
|---|---|---|
| `Complex<D>` | `a + b*im` with `a`, `b` of dimension `D` | the `harmonic` procedure (§14.2): acoustics, electromagnetics, vibration |
| `Quaternion<D>` | four components of dimension `D` | rotations in rigid-body ODEs and beams with large rotations |
| `Dual<D>` | `a + b*dual_eps` with `dual_eps² = 0` | forward-mode derivatives, for example the sensitivity of a report to a `param` |

A tensor takes a number system in its argument: `Vector<Complex<Velocity>>`.
The built-in constants `im` and `dual_eps` are the units of the imaginary
part and of the dual part.
Reason: these four cover the number systems that engineering models use in
practice. Each one needs its own solver arithmetic, so the list stays
closed.

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
4. `Complex<D>`, `Quaternion<D>` and `Dual<D>` accept a value of `D`.

Reason for exception 1: `u = 0` is the most common boundary condition, and
`u = (0, 0, 0) mm` adds tokens and gives no safety.

---

## 6. Dimensions and units

### 6.1 Dimensions

A dimension is a vector of seven integer exponents over the SI base
dimensions: mass M, length L, time T, temperature Θ, current I, amount N and
luminous intensity J. The checker does arithmetic on these vectors:

- `a*b` adds the exponents, `a/b` subtracts them, `a^k` multiplies them by
  `k`. `k` must be an integer constant, except when `a` is `Real`. Then `k`
  may be any constant `Real`, not a varying value (§7.6) [`not_constant`].
  A function parameter is constant when each call gives it a constant.
- `a + b`, `a - b`, `a == b` and the other comparisons need equal
  dimensions [`unit_mismatch`].
- `sqrt(a)` halves the exponents. Each exponent must be even
  [`unit_fractional`].
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
numbers in the solver. Each model with dimensions has one
[`scales_missing`].

The header `model m: 2d, nondim` declares a model with no dimensions. Each
quantity has type `Real`, and a literal with a unit is an error
[`nondim_unit`]. Use it for
a PDE that has no physical units, for example a model problem from a
textbook in dimensionless form. The dimension checks then give no
protection, so the `check` report says `units: nondim` at its top.

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
| `shell(a, t)` | same as `a` | a wall of constant thickness `t` inside the surface of `a`: `max(φa, −φa − t)` |
| `loft(s1, s2, h)` | `Solid` | from 2D section `s1` at z = 0 to `s2` at z = `h` |
| `pattern(a, v, n)`, `pattern_polar(a, axis, n)` | same as `a` | `n` copies of `a`, along `v` or about `axis` |
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

Regions combine with `|`, `&` and `-`. The checker computes containment by
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
| `cell` | `tri`, `quad`, `tet`, `hex` | `tri` in 2D, `tet` in 3D | `quad` and `hex` need a `box` or `rect` region, and `hex` can also mesh a quarter of a tube |
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
| `quad`, `hex` | a structured grid on a `box` or `rect`, or in r, θ and z on a quarter of a tube: `(cylinder - cylinder) & box`, with the box corner on the axis |
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
for example a function of a temperature field.

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
| `plane_stress(t)` | a `2d` model only. σzz = σxz = σyz = 0. The compiler solves the constitutive law for εzz and removes it (§11.2). Thickness `t` multiplies volume and surface integrals. | `plane_strain`, a `3d` model |
| `plane_strain(t)` | a `2d` model only. εzz = εxz = εyz = 0. With no t, integrals are per unit depth (§14.5). | `plane_stress` |
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
| `F*n` or `dot(F, n)`, where `div(F)` is in an equation of the physics | natural (Neumann). The right side may use the unknowns (Robin). |
| anything else | an error [`bc_form`] |

A natural condition on the flux of the equation paired with unknown `u` is a
condition for `u`.

On a `Point` region, `n` has no value. There the checker reads a natural
condition `F*n = P` as the force `P` at the point. `P` has the dimension of
`F` times area in a `3d` model, and of `F` times length in a `2d` or `axi`
model.

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

### 12.6 Symmetry groups

```skein
symmetry rotor: cyclic { n: 24, axis: (0, 0, 1), low: face_a, high: face_b }
```

A `symmetry` item solves one sector of a body that repeats. It ties the
unknowns on the boundary `high` to the unknowns on `low`, through the
transformation of the group. The list of groups is closed:

| Group | Settings | Condition on `high` |
|---|---|---|
| `cyclic` | `n`, `axis`, `low`, `high`, `harmonic` (Int, default 0) | the value on `low`, rotated by 2π/n and multiplied by `exp(im*2*pi*harmonic/n)` |
| `dihedral` | `n`, `axis`, `low`, `high` | `cyclic` plus a mirror plane on each side of the sector |
| `bloch` | `period` (a vector), `low`, `high`, `k` (a wave vector) | the value on `low` multiplied by `exp(im*dot(k, period))` (Bloch and Floquet) |

The checker tests that the transformation maps `low` onto `high`
[`symmetry_map`]. A `harmonic` other than 0, or a `k` other than 0, needs
`Complex` unknowns (§5.1). A `modal` step over the harmonics gives the
modes of the full body.

Reason: rotors, gears and periodic structures are common, and a sector
model is many times smaller than the full body. A closed list of groups
keeps each one testable.

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

Reason: these defaults are the right choice for most solid, heat and
laminar flow models, and they cut a block of 4 lines from each of them.
`skein check` prints the method that it chose.

The kinds:

| Kind | Settings | Notes |
|---|---|---|
| `cg` | `mesh`, `order` (Int or a record per field), `family` (`lagrange`, `nedelec`, `raviart_thomas`), `stabilize`, `quadrature`, `couple` | continuous Galerkin FEM. `nedelec` (H(curl)) and `raviart_thomas` (H(div)) are phase 5. |
| `dg` | `mesh`, `order`, `flux` (`upwind`, `lax_friedrichs`, `roe`), `penalty`, `shock_capturing` (`none`, `entropy_viscosity`) | interior penalty for second-order terms |
| `fv` | `mesh`, `gradient` (`green_gauss`, `least_squares`), `flux` (`upwind`, `linear_upwind`, `central`), `limiter` (`none`, `minmod`, `van_leer`, `venkatakrishnan`), `pv_coupling` (`rhie_chow`) | cell-centred finite volume |
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
of each unknown: the expression after `=` in `solves`, or zero.

```skein
physics thermal on part solves T: Temperature = 293 K { … }
```

### 14.2 Procedures

| Procedure | Arguments | Meaning |
|---|---|---|
| `static` | `ramp` (Int, default 1) | `ddt` = 0. Nonlinear problems take `ramp` load increments. |
| `transient` | `end`, `dt`, `scheme`, `adapt` | time stepping from the end time of the step before. `adapt: { rtol, atol }` changes `dt` to keep the estimated local error below the limits. `dt` is then the first step. |
| `harmonic` | `freq` (a list of `Frequency`) | steady response to loads of the form `F·exp(im·ω·t)`. The unknowns and loads are `Complex` (§5.1). |
| `modal` | `count`, `shift` | the `count` lowest eigenpairs of the linearized operator about the current state |
| `buckling` | `count` | eigenpairs of the stiffness and the geometric stiffness of the current state |

Time schemes:

| `scheme` | Order | Kind | Default for |
|---|---|---|---|
| `backward_euler` | 1 | implicit | |
| `bdf2` | 2 | implicit | first-order in time (`ddt`) |
| `bdf` | 1 to 5, variable | implicit | `adapt` with `ddt` (stiff ODEs) |
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
`reaction(r)`, `max(f, r)`, `min(f, r)`, `integral(f, r)`, `mean(f, r)`,
`eigenvalue(k)`.

`eigenvalue(k: Int)` is the eigenvalue of mode `k` of a `modal` step, or
of a `buckling` step. For a `modal` step, it is the natural frequency
f = ω/(2π) in Hz, a `Frequency`. For a `buckling` step, it is a `Real` load
factor. Mode 1 is the lowest. In a `report`, the step is the step of the
report. In an `expect`, the step is the last `modal` or `buckling` step. A
use in another place, or with no such step, is an error
[`eigenvalue_step`]. A `k` outside `1..count` of the step is also an error
[`eigenvalue_step`].

`report cd = expr` gives the value a name. A study (§20) uses the name as
the column of its results table. A report with no name gets its source
text as the column name.

`integral(f, r)` has the dimension of `f` times `Length^k`, where `k` is
the topological dimension of `r`. In a `2d` model with no thickness, a force
on a curve is thus `Force/Length`, the force per unit depth.

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
  lowering also keeps a split `F = F_a(U) + F_d(U, ∇U)`: the advective part
  has no derivative of an unknown, the diffusive part has. Only the flux
  choices of `dg` and `fv` and the stabilization (§17.2) use this split.
  The natural conditions (§12.1) and the `cg` lowering (§17.1) use the full
  flux `F`.
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

The WebAssembly target uses the C that Goose writes: `emcc` compiles that C
to a WebAssembly module. Earlier projects of the author of Skein, beamdiag
and goosestrike, already run in the browser in this way. Thus the browser
target does not depend on a WebAssembly back end in Goose.

### 18.2 Files

`skein build model.skn` writes:

| File | Contents |
|---|---|
| `model.goose` | the kernels and the plan of the model |
| `runtime` (package) | the runtime: mesher, sparse matrices, solvers, output writers. The same for each model. |

### 18.3 The shape of the generated code

The plan stores a mesh as flat arrays, which Goose keeps in data stacks with
no pointers. A cell kernel is a `simd fn`. Assembly runs in worker threads
over a colouring of the cells, so no two workers write one row.

Each generated function has at most 60 lines (rule 4 of the Power of
Ten). When the kernel of a cell is longer, the emitter splits it into
functions by blocks of the cell matrix, one block for each pair of
unknowns or of node groups, and one function calls them in order. Reason:
a person can review a function of 60 lines, and a test checks the length
of each generated function.

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
        let job = qget<runtime::ColourJob>();
        guard job.first >= 0 else { break; }
        qput(runtime::assemble_range(job, k_solid_tet4));
    }
}
```

### 18.4 The runtime

`runtime` is a Goose package in `src/runtime/`. It holds no C, except the C
that the Goose standard library uses for files and the clock. The
[end-to-end slice](roadmap.md#the-end-to-end-slice) ships a box mesher, CG
with Jacobi and the VTU writer before phase 2.

| Part | Contents | Phase |
|---|---|---|
| sparse | CSR matrices, colouring, assembly | 2 |
| linear | CG, GMRES(m), Jacobi, ILU(0) | 2 |
| nonlinear | Newton with line search, load ramps | 2 |
| eigen | LOBPCG (Knyazev 2001) for `modal` and `buckling` | 3 |
| multigrid | smoothed aggregation AMG (Vaněk, Mandel and Brezina 1996) | 3 |
| mesher | §9.2 | 2 |
| output | VTU writer, TOON reports, state files | 2 |
| tasks | the work-stealing task pool of §20.5 | 3 |

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
11. Parameters, scales and studies (§20, §21), manufacture and tolerance
    items (§22).
12. `use` files, functions, templates and symmetry groups (§4.1, §4.2,
    §7.4, §12.6).
13. `allow` names a check that exists, with a reason that is not empty
    [`allow_reason`].

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
| `minimize(objective, method: nelder_mead \| bobyqa, max_cases: Int)` | a range for each parameter, and the start value is the default | cases chosen one batch at a time to minimize a named `report` |

- Each entry names a `param` [`study_param`]. Each value has the type of the
  parameter [`type_mismatch`] and is in its range [`param_range`].
- A study does not change the steps. Each case runs all the steps.
- A model may hold many studies. `skein run model.skn` runs the model once,
  with the defaults. `skein run model.skn --study re_sweep` runs the study.
- `minimize` is gradient-free. Gradients by adjoint are in the roadmap
  (phase 5). `objective` names a scalar `report` of the last step
  [`study_objective`].

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

A machine has `H` hardware threads. A study with `k` workers gives each case
`max(1, H / k)` threads for assembly and solves. The default `k` is
`min(H, number of cases)`. `--jobs k` sets it. Reason: cases are
independent, so parallel cases scale better than parallel assembly in one
case. But a study with fewer cases than threads still uses all threads.

The memory of one case times `k` must fit the target. The runner estimates
the memory of a case from the mesh size before it starts the workers, and
lowers `k` if necessary. It reports the value that it chose.

**Phase 3: one work-stealing task pool.** Phase 2 runs cases first, as
above. Phase 3 replaces the two levels of threads with one pool in the style
of the Go scheduler. Each hardware thread is a worker with its own queue of
tasks. A task is one case, or one part of the work in a case: a colour of
the assembly, a block of a solver. A worker that has no task takes one
from the queue of another worker. Thus a study with 3 large cases on 16
threads keeps each thread busy, and the user sets no `--jobs`.

The pool needs atomic operations for the queues. Before phase 3 starts,
check whether Goose gives atomics. If it does not, each queue uses a lock,
which costs more time for each task but gives the same results.

### 20.6 The results table

The study writes `.skein/studies/<model>/<study>.toon` and prints the same
table:

```toon
study: re_sweep
cases[3]{case,Re,status,verify,cd,time_s}:
  1,10,ok,pass,2.81,4.2
  2,20,ok,pass,2.05,4.6
  3,40,ok,pass,1.54,5.1
```

One column for each varied parameter, each named report and each `expect`.
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

## 22. Manufacture

Skein is a simulation language first. But the geometry of a model is a
parametric part, and the same parameters drive the part, the simulation
and the studies. This section adds three things that use that fact: the
manufacturing process of a part, its checks, and tolerance stack-ups. It
also says how the part leaves Skein for manufacture.

### 22.1 Parametric parts

A part is a geometry value built from `param` values (§8, §20):

```skein
param t_wall: Length = 2 mm tol 0.1 mm in 1 mm..5 mm
param r_fil: Length = 1.5 mm
param d_hole: Length = 6 mm tol 0.02 mm

let outer: Solid = box((0, 0, 0) mm, (80, 40, 30) mm)
let inner: Solid = move(box((0, 0, 0) mm, (80 mm - 2*t_wall, 40 mm - 2*t_wall, 30 mm)), (t_wall, t_wall, t_wall))
let housing: Solid = fillet(outer - inner, r_fil) - cylinder((40, 20, 0) mm, (0, 0, 1), d_hole/2, t_wall)
```

This is the form that text-to-CAD tools produce from a sentence: a short
script whose parameters have names. In Skein the script has types and
units, so `check` finds a wall thickness in `Pressure` or a hole that the
range of its parameter makes larger than the part.

### 22.2 `manufacture`

```skein
manufacture housing: injection_mould {
    pull: (0, 0, 1)
    draft: 1 deg
    min_wall: 1.2 mm
    max_wall: 4 mm
}
```

`manufacture r: process { … }` declares the process for region `r`. The
processes and their settings:

| Process | Settings |
|---|---|
| `cnc_mill` | `axes` (3 or 5), `tool_radius`, `max_depth_ratio`, `approach` (a list of directions for 3 axes) |
| `injection_mould` | `pull`, `draft`, `min_wall`, `max_wall` |
| `casting` | `pull`, `draft`, `min_wall` |
| `sheet_metal` | `thickness`, `min_bend_radius`, `min_flange` |
| `additive` | `build` (direction), `max_overhang`, `min_wall`, `min_hole` |
| `turning` | `axis`, `tool_radius` |

A setting that does not belong to the process is an error
[`process_setting`].

### 22.3 Design-for-manufacture checks

`skein verify --dfm` runs these checks on each `manufacture` region. They
are `verify` checks, so `allow` with a reason can turn one off. They need no
solve, only the geometry.

| Id | Test | Processes |
|---|---|---|
| `dfm_wall` | the local wall thickness is between `min_wall` and `max_wall` | all with `min_wall` |
| `dfm_draft` | each face that the pull direction drags along has at least the `draft` angle | `injection_mould`, `casting` |
| `dfm_undercut` | each surface is visible from the pull or approach directions | `injection_mould`, `casting`, `cnc_mill` |
| `dfm_inner_radius` | each concave edge has a radius of at least `tool_radius` | `cnc_mill`, `turning` |
| `dfm_hole` | each hole has a diameter of at least `min_hole`, and a depth of at most `max_depth_ratio` times the diameter | `additive`, `cnc_mill` |
| `dfm_overhang` | each downward face is within `max_overhang` of the build direction | `additive` |
| `dfm_bend` | each bend radius is at least `min_bend_radius`, each flange at least `min_flange` | `sheet_metal` |

The checks work on the distance function. The local wall thickness is twice
the radius of the largest ball inside the part at that point, which the
distance function gives directly. Draft and overhang use the normal ∇φ.
Undercuts and tool access use rays along the directions, on the distance
function. Each check reports the worst place as a point, so the viewer can
show it.

Each check runs on the default parameter values, and on each case of a
study. Thus a sweep of `t_wall` shows which cases cannot be made.

### 22.4 Tolerances and stack-ups

`tol` on a `param` gives a symmetric tolerance. A `tolerance` item states a
stack-up:

```skein
param L: Length = 80 mm tol 0.1 mm
param a: Length = 39.8 mm tol 0.05 mm
param b: Length = 39.8 mm tol 0.05 mm
let gap: Length = L - a - b

tolerance fit: gap in 0.1 mm..0.7 mm { method: worst_case }
```

The compiler differentiates the expression with respect to each parameter
that has a tolerance. It reports the nominal value and the range by the
worst case (the sum of |∂f/∂p|·tol) and by root sum square. The check
`dfm_stack` fails when the range of the chosen `method` (`worst_case`,
`rss` or `monte_carlo`) is not inside the limits.

`method: monte_carlo` is for a stack-up that is not linear, for example a
gap that depends on an angle. It runs a `sample` study (§20.2) with each
toleranced parameter uniform in its tolerance band, and takes the range of
the expression over the samples. `samples: n` sets the number (default
10000). The expression needs no solve, so the study is fast. The expression must be a constant: it can
use geometry and parameters, but not fields [`not_constant`].

### 22.5 Export for manufacture

```
skein export housing.skn --geometry housing --format stl
```

| Format | Contents | Phase |
|---|---|---|
| `stl` | a triangulated surface from the mesher, at the size of `--size` | 2 |
| `step` | AP242 with tessellated geometry from the mesher | 4 |

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
| `von_mises(sigma)`, `tresca(sigma)`, `principal(sigma)` | stress measures |
| `fourier(k, T)` | heat flux `−k grad(T)` |
| `abs`, `sqrt`, `exp`, `ln`, `sin`, `cos`, `tan`, `tanh`, `atan`, `atan2(y, x)`, `min`, `max`, `pi` | scalar math |
| `re(c)`, `imag(c)`, `abs(c)`, `arg(c)`, `conj(c)` | parts of a `Complex<D>`: `re`, `imag` and `abs` give `D`, `arg` gives `Angle`, `conj` gives `Complex<D>` |
| `diff(W, A)` | symbolic derivative, for example `S = 2*diff(W, C)` from an energy |
| `reaction(r)`, `integral(f, r)`, `mean(f, r)`, `eigenvalue(k)` | reports (§14.5) |
| `results(file, step:)` | §15.1 |
| `boundary(g, s)`, geometry functions | §8 |

### Constants

`I` (identity), `pi`, `x`, `y`, `z`, `t`, `n`, `im` (imaginary unit),
`dual_eps` (dual unit).

### Standard library

The physics templates of §4.2. `skein help std` prints their source.

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
