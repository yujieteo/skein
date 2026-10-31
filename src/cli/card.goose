// The language card that `skein help --llm` prints (docs/command_line.md,
// docs/llm_test.md §2). It is the only Skein text that an LLM gets before it
// writes a model, so it must stay true to the checker. Target: about 2000
// tokens.
namespace cli;

let CARD = """
    # Skein language card

    Skein is a typed language for physics simulation models. A file `name.skn`
    holds one model. `skein check name.skn` parses and type-checks it and prints
    errors in TOON: `class,id,file,line,col,msg,fix`. `skein run name.skn`
    meshes, solves and prints each `report` value.

    ## Layout

    - The first line is the header: `model NAME: 3d, mm_N_MPa_s` (space `1d`,
      `2d`, `3d` or `axi`; unit system `SI`, `mm_N_MPa_s` or `mm_kg_ms`).
    - One statement on each line. A newline ends it. No semicolons.
    - `//` starts a comment.
    - Declare a name before you use it. Steps run in file order.
    - Built-in names that a declaration may not use: `x`, `y`, `z`, `t`
      (time), `n`, `I`, `pi`, `im`, and the built-in functions.
    - Indent by 4 spaces. A block is `{ ... }` with one entry on each line, or
      `{ a: 1, b: 2 }` on one line.

    ## Values, types and units

    - `let NAME: Type = expr` declares a fixed value. The type is required.
    - `param NAME: Type = expr in lo..hi` declares a value that a `study` or
      `skein run --set NAME=value` can change. The range is optional.
    - Each literal carries its unit: `210 GPa`, `5 mm`, `2 kN`, `7850 kg/m^3`.
      A vector literal: `(0, 0, -1) MPa` or `(L, W, H)`. A bare `0` is zero in
      any dimension. The checker tests the dimension of each expression.
    - Scalar types: `Real`, `Int`, `Length`, `Area`, `Volume`, `Mass`, `Time`,
      `Force`, `Pressure` (= `Stress`), `Density`, `Temperature`, `Velocity`,
      `Acceleration`, `Energy`, `Power`, `Frequency`, `Angle`, or a product
      such as `Force/Length`.
    - Tensor types: `Vector<Length>`, `Tensor2<Pressure>`.
    - Geometry types: `Solid`, `Surface`, `Curve`, `Point`.
    - Units: `m`, `mm`, `cm`, `km`, `g`, `kg`, `t`, `s`, `ms`, `min`, `h`,
      `N`, `kN`, `MN`, `Pa`, `kPa`, `MPa`, `GPa`, `bar`, `J`, `W`, `Hz`, `K`,
      `degC`, `deg`, `rad`. Combine with `*`, `/`, `^`.
    - Operators: `+ - * / ^`, `== != < <= > >=`, `&& || !`. On geometry, `|`
      is union, `&` intersection and `-` difference.
    - `a.x`, `a.y`, `a.z` are the components of a vector. `sigma.xx` is a
      component of a tensor. `f at (x, y, z)` is the value of field `f` at a
      point.

    ## Required line: `scales`

    Each model with units has one `scales` line after its parameters. It names
    reference values whose dimensions together make each dimension of the
    model, and no name is a product of the others. For a static solid:
    `scales H, P/(W*H)` (a length and a stress).

    ## Geometry and regions

    - Solids: `box(p0, p1)` (opposite corners), `cylinder(base_centre, axis,
      r, h)`, `sphere(c, r)`. Combine with `|`, `&`, `-`.
    - `boundary(g, x == 0 mm)` is the part of the boundary of `g` where the
      condition holds. `boundary(g, other)` is the part on the surface of
      `other`. `boundary(g)` is the whole boundary.
    - `x`, `y`, `z` are the coordinates. `n` is the outward unit normal.

    ## Material, mesh, physics, method

    ```
    material steel {
        E: Pressure = 210 GPa
        nu: Real = 0.3
    }

    mesh m on part {
        size: 5 mm
        cell: hex
    }

    physics solid on part: linear_elastic solves u: Vector<Length> {
        material steel
    }
    ```

    - `mesh` settings: `size` (required), `cell` (`tet` by default, `hex` for a
      `box` or a quarter of a tube only), `refine region: size`.
    - `linear_elastic` is a standard template: small strain, isotropic. It
      needs the material properties `E` and `nu`. It defines the stress
      `sigma` for use in `bc`, `output` and `report`. Its body force is 0
      unless the model gives `b`, for example `linear_elastic(b: rho*(0, 0, -g))`
      on `tet` cells. `skein help std` prints the source of the templates.
    - The default method is continuous Galerkin of order 2 on the one mesh of
      the region. Write it only to change it:
      `method solid: cg { mesh: m, order: 2 }`.

    ## Shells

    `skein run` solves a shell on a part of the side of a cylinder on the x
    axis (MITC4 quads, order 1). theta turns from +z to +y.

    ```
    let tube: Solid = cylinder((0, 0, 0) m, (1, 0, 0), R, L)
    let wall: Surface = boundary(tube, y^2 + z^2 == R^2 && y >= 0 m)
    let x0: Curve = boundary(wall, x == 0 m)
    let e_t: Vector<Real> = (0, z/R, -y/R)
    table plies{angle: Angle in deg}:
        45
        0
        -45
        90
    material ply {
        E1: Pressure = 148 GPa
        E2: Pressure = 10.1 GPa
        G12: Pressure = 3.94 GPa
        nu12: Real = 0.32
        t_ply: Length = 0.191 mm
        h: Length = 24*0.191 mm
    }
    mesh m on wall {
        size: 50 mm
        cell: quad
    }
    physics skin on wall: reissner_mindlin(layup: plies.angle, mirror: true, top: n, q: p*n) solves w: Vector<Length>,
        theta: Vector<Angle> {
        material ply
        let s_hoop: Pressure = dot(e_t, N*e_t)/h
    }
    method skin: cg { mesh: m, order: 1 }
    ```

    - `layup` is one ply angle or a table column of angles from the top face.
      The plies repeat to fill `h` (`h/2` with `mirror: true`).
    - `N` is the force per length and `eps_m` the membrane strain (tensors).
    - A `Curve` is a side of the shell: `boundary(wall, x == 0 m)` or
      `boundary(wall, y == 0 m)`. A `Point` is `point((x, y, z))`.
    - `bc x0: w.x = 0`, `bc x0: theta.y = 0`, `bc end: N*n = (f, 0, 0)` (a
      force per length on a side).
    - Report `w.z at p` or a physics `let` at p (p a node), or
      `reaction(x0)`.

    ## Boundary conditions

    - `bc region: u = 0` fixes the displacement. `bc region: u.x = 0` fixes one
      component.
    - `bc region: sigma*n = (0, 0, -p)` is a traction (a force per area) on
      a surface. A force `P` on a face of area `A` is the traction `P/A`.
    - A `bc` at the top level applies in all steps. A `bc` in a step applies
      from that step on, and replaces a `bc` on the same region and target.
    - Each part of the boundary with no `bc` is traction free. `skein check`
      lists those parts.

    ## Steps, output and report

    ```
    step NAME: static {
        bc tip: sigma*n = (0, 0, -s)
        output u, sigma
        report reaction(root)
        report tip_uz = u.z at (L, W/2, H/2)
    }
    ```

    - Procedures: `static`, `transient(end: 1 s, dt: 0.01 s)`,
      `modal(count: 6)`.
    - `output` writes fields to the result file (VTU).
    - `report` prints a value in TOON. `report name = expr` names it.
      Functions: `reaction(region)`, `max(f, region)`, `min(f, region)`,
      `mean(f, region)`, `integral(f, region)`, `von_mises(sigma)`.
    - A second step starts from the state at the end of the first.

    ## Studies

    ```
    param W: Length = 20 mm in 5 mm..100 mm

    study width: grid {
        W: [10, 20, 40] mm
    }
    ```

    A study runs the model once for each value (`grid`: all combinations).
    Each entry names a `param`. `skein run file.skn --study width` runs it.

    ## Expectations

    `expect name: u.z at (L, 0 mm, 0 mm) == -0.5 mm within 2%` tests a result
    after the last step.

    ## Example

    ```
    // An aluminium plate 200 x 100 x 5 mm with a 10 mm hole, pulled in x.

    model plate: 3d, mm_N_MPa_s

    let a: Length = 200 mm
    let w: Length = 100 mm
    let th: Length = 5 mm
    let r: Length = 10 mm
    param S: Pressure = 50 MPa in 1 MPa..200 MPa

    scales w, S

    let hole: Solid = cylinder((a/2, w/2, 0 mm), (0, 0, 1), r, th)
    let plate: Solid = box((0, 0, 0) mm, (a, w, th)) - hole
    let left: Surface = boundary(plate, x == 0 mm)
    let right: Surface = boundary(plate, x == a)

    material aluminium {
        E: Pressure = 70 GPa
        nu: Real = 0.33
    }

    mesh m on plate {
        size: 5 mm
        refine boundary(plate, hole): 1 mm
    }

    physics solid on plate: linear_elastic solves u: Vector<Length> {
        material aluminium
    }

    bc left: u = 0

    step pull: static {
        bc right: sigma*n = (S, 0, 0)
        output u, sigma, von_mises(sigma)
        report reaction(left)
        report peak = max(von_mises(sigma), plate)
    }

    study load: grid {
        S: [25, 50, 100] MPa
    }
    ```
    """;
