// Tests of the code generator. `./build.sh test` runs this program from the
// root of the repository, and a non-zero exit status is a failure.
import std;
import os;
import syntax from "../syntax";
import .codegen;

var failures = 0;

// Reads `text` and gives the generated program, or the first error.
fn generate(text: u8[:]) -> u8[>..] {
    syntax::reset();
    let f = syntax::add_file("t.skn", text);
    syntax::lex(f);
    let prog = syntax::parse(f);
    var out: u8[>..] = [];
    if syntax::errors.len > 0 || !codegen::read(prog) {
        let e = syntax::errors[0];
        let t = syntax::tokens[e.tok];
        format(out, "error ", e.id, " ", t.line, ":", t.col, " ", e.msg);
        return out;
    }
    codegen::generate("RUNTIME", "RESULTS", out);
    return out;
}

// The generated program for `text` holds each of `wants`.
fn expect_lines(what: u8[:], text: u8[:], wants: const u8[:][:]) {
    let got = generate(text);
    for w in wants {
        if find(got, w) < 0 {
            failures++;
            print("FAIL ", what, ": no '", w, "' in:\n", got);
            return;
        }
    }
}

fn expect_error(what: u8[:], text: u8[:], want: u8[:]) {
    let got = generate(text);
    if find(got, want) != 0 {
        failures++;
        print("FAIL ", what, "\n--- got:\n", got, "\n--- want:\n", want);
    }
}

// The number of times `w` is in `s`.
fn occurrences(s: u8[:], w: u8[:]) -> i64 {
    var n = 0;
    var at = 0;
    while true {
        let i = find(s[at..], w);
        if i < 0 { return n; }
        n++;
        at += i + w.len;
    }
    return n;
}

// `s` with each `from` changed to `to`.
fn replace(s: u8[:], from: u8[:], to: u8[:]) -> u8[>..] {
    var out: u8[>..] = [];
    for part, i in split(s, from) {
        if i > 0 { out.append(to); }
        out.append(part);
    }
    return out;
}

fn model(header: u8[:], body: u8[:]) -> u8[>..] {
    return str(header, """

        let L: Length = 100 mm
        let beam: Solid = box((0, 0, 0) mm, (L, 10 mm, 10 mm))
        let root: Surface = boundary(beam, x == 0 mm)
        let tip: Surface = boundary(beam, x == L)

        material steel {
            E: Pressure = 210 GPa
            nu: Real = 0.3
        }

        """, body);
}

let MESH = """
    mesh m on beam {
        size: 2.5 mm
        cell: hex
    }

    physics solid on beam: linear_elastic solves u: Vector<Length> {
        material steel
    }

    """;

let STEP = """
    step load: static {
        bc root: u = 0
        bc tip: sigma*n = (0, 0, -1) MPa
        output u
        report reaction(root)
    }
    """;

let MM = "model t: 3d, mm_N_MPa_s";

fn test_grid_cantilever() {
    var text: u8[>..] = [];
    read_file("examples/femahoy/grid_cantilever.skn", text);
    expect_lines("grid_cantilever", text, [
        "import runtime from \"RUNTIME\";",
        "runtime::mesh_box(0.0, 0.0, 0.0, 100.0, 10.0, 10.0,",
        "runtime::cells_for(100.0, 2.5),",
        "runtime::material(210000.0, 0.3);",
        "face(runtime::fix(\"x0\"), \"root\");",
        "face(runtime::traction(\"x1\", 0.0, 0.0, -1.0), \"tip\");",
        "largest = max(largest, 1.0 * runtime::face_area(\"x1\"));",
        "equilibrium(rx, ry, rz, largest, \"N\");",
        "solve(\"load\", \"cg\", \"jacobi\");",
        "let x, y, z = runtime::reaction(\"x0\");",
        "report(\"load\", \"reaction(root)\", vector(x, y, z, \"N\"));",
        "output(\"RESULTS/load.vtu\");",
        "let ux, uy, uz = runtime::displacement_at(100.0, 5.0, 5.0);",
        "expect(\"tip_uz\", uz, -0.191933, 0.015, true, \"mm\", \"-0.191933 mm within 1.5%\");",
        "let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(51.25, 5.0, 10.0);",
        "expect(\"top_sxx\", sxx, 29.25, 0.015, true, \"MPa\", \"29.25 MPa within 1.5%\");",
    ]);
}

fn test_closed_cylinder() {
    var text: u8[>..] = [];
    read_file("examples/closed_cylinder.skn", text);
    expect_lines("closed_cylinder", text, [
        "runtime::mesh_annulus(0.0, 0.0, 100.0, 105.0, 0.0, 1.5707963267948966, 0.0, 40.0,",
        "runtime::cells_for(1.5707963267948966 * 105.0, 2.5),",
        "face(runtime::fix_component(\"t1\", runtime::X), \"sym_x\");",
        "face(runtime::fix_component(\"t0\", runtime::Y), \"sym_y\");",
        "face(runtime::fix_component(\"z0\", runtime::Z), \"bottom\");",
        "face(runtime::normal_traction(\"r0\", -1.0), \"inner\");",
        "face(runtime::traction(\"z1\", 0.0, 0.0, 9.75609756097561), \"cap\");",
        "let rx, ry, rz = runtime::total_reaction();",
        "largest = max(largest, 1.0 * runtime::face_area(\"r0\"));",
        "let x, y, z = runtime::reaction(\"z0\");",
        "report(\"press\", \"reaction(bottom)\", vector(0.0, 0.0, z, \"N\"));",
        "let sxx, syy, szz, syz, sxz, sxy = runtime::stress_at(102.5, 0.0, 20.0);",
        "expect(\"hoop_lame\", syy, 19.993906066365838, 0.001, true, \"MPa\", \"p*a^2/(b^2 - a^2)*(1 + b^2/r_m^2) within 0.1%\");",
        "expect(\"axial_force\", z, -7853.981633974483, 1e-06, true, \"N\", \"-p*pi*a^2/4 within 0.0001%\");",
    ]);
    expect_lines("a label that starts and ends with a group",
                 replace(text, "== p*a^2/(b^2 - a^2)*(1 + b^2/r_m^2)", "== (p*a^2/(b^2 - a^2))*(1 + b^2/r_m^2)"),
                 ["\"(p*a^2/(b^2 - a^2))*(1 + b^2/r_m^2) within 0.1%\""]);
    let tube = "let wall: Solid = (barrel - bore) & quarter";
    expect_error("a full tube", replace(text, tube, "let wall: Solid = barrel - bore"),
                 "error unsupported 28:31 the solid 'wall' is not a box or a quarter of a tube");
    expect_error("a cut off the centre", replace(text, "box((0, 0, 0) mm, (b, b, L))", "box((1, 0, 0) mm, (b, b, L))"),
                 "error unsupported 27:35 this slice cuts only the quarter x >= cx, y >= cy of a tube");
    expect_lines("an edge with the same fixed component", replace(text, "bc bottom: u.z = 0", "bc bottom: u = 0"), [
        "let x, y, z = runtime::reaction(\"z0\");\n        var alone = true;\n        if shared(\"bottom\", \"z0\", \"sym_x\", \"t1\") { alone = false; }\n        if shared(\"bottom\", \"z0\", \"sym_y\", \"t0\") { alone = false; }\n        if alone { report(",
        "let x, y, z = runtime::reaction(\"z0\");\n        var alone = true;\n        if alone { expect(\"axial_force\"",
    ]);
    expect_error("a reaction with no support", replace(text, "reaction(bottom).z", "reaction(bottom).x"),
                 "error unsupported 76:8 expect axial_force: the face bottom does not fix u.x");
    expect_error("a transposed sigma component", replace(text, "sigma.zz at (r_m, 0 mm, L/2) == p*a/(2", "sigma.zy at (r_m, 0 mm, L/2) == p*a/(2"),
                 "error unsupported 63:26 sigma has no component 'zy'");
}

fn test_units() {
    // 1 kN on 50 mm x 100 mm is 0.2 MPa, and in SI 200000 Pa.
    let load = """
        step load: static {
            bc root: u = 0
            bc tip: sigma*n = (0, 0, -1 kN/(50 mm*100 mm^0.5*mm^0.5))
            report tip_uz = u.z at (1, 0.025, 0.05) m
        }
        """;
    expect_lines("units in mm_N_MPa_s", model(MM, str(MESH, load)), [
        "traction(\"x1\", 0.0, 0.0, -0.2000000000",
        "displacement_at(1000.0, 25.0, 50.0)",
        "quantity(uz, \"mm\")",
    ]);
    expect_lines("units in SI", model("model t: 3d, SI", str(MESH, load)), [
        "runtime::mesh_box(0.0, 0.0, 0.0, 0.1, 0.01, 0.01,",
        "runtime::material(210000000000.0, 0.3);",
        "traction(\"x1\", 0.0, 0.0, -200000.0000",
        "displacement_at(1.0, 0.025, 0.05)",
    ]);
}

fn test_two_steps() {
    let steps = """
        bc root: u = 0

        step one: static {
            bc tip: sigma*n = (0, 0, -1) MPa
        }

        step two: static {
            bc tip: sigma*n = (0, 0, -2) MPa
        }
        """;
    expect_lines("two steps", model(MM, str(MESH, steps)), [
        "fn step_one() {",
        "fn step_two() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -2.0) MPa\n    face(runtime::traction(\"x1\", 0.0, 0.0, -2.0), \"tip\");",
        "    step_one();\n    step_two();\n",
    ]);
}

// Spec §12.2: a top-level bc holds in each step wherever it is in the
// file, and a step bc holds from its step on.
fn test_bc_scope() {
    let after = """
        step load: static {
            bc root: u = 0
        }

        bc tip: sigma*n = (0, 0, -1) MPa
        """;
    expect_lines("a top-level bc after the step", model(MM, str(MESH, after)), [
        "fn step_load() {\n    // bc tip: sigma*n = (0.0, 0.0, -1.0) MPa\n    face(runtime::traction(\"x1\", 0.0, 0.0, -1.0), \"tip\");",
    ]);
    let carry = """
        let bottom: Surface = boundary(beam, z == 0 mm)

        step one: static {
            bc root: u = 0
            bc tip: sigma*n = (0, 0, -1) MPa
        }

        step two: static {
            bc bottom: sigma*n = (0, 0, 1) MPa
        }

        step three: static {
            bc tip: sigma*n = (0, 0, -3) MPa
        }
        """;
    expect_lines("step bcs carry on", model(MM, str(MESH, carry)), [
        "fn step_two() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -1.0) MPa\n    face(runtime::traction(\"x1\", 0.0, 0.0, -1.0), \"tip\");\n    // bc bottom: sigma*n = (0.0, 0.0, 1.0) MPa\n    face(runtime::traction(\"z0\", 0.0, 0.0, 1.0), \"bottom\");\n    solve(",
        "fn step_three() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -3.0) MPa\n    face(runtime::traction(\"x1\", 0.0, 0.0, -3.0), \"tip\");\n    // bc bottom:",
    ]);
}

// The equilibrium check sums the reaction of each fixed node once.
fn test_reaction() {
    let twice = """
        bc root: u = 0

        step load: static {
            bc root: u = 0
            bc tip: sigma*n = (0, 0, -1) MPa
        }
        """;
    let got = generate(model(MM, str(MESH, twice)));
    if occurrences(got, "runtime::total_reaction()") != 1 || occurrences(got, "runtime::fix(") != 1 {
        failures++;
        print("FAIL a face fixed twice: not one fix and one reaction in:\n", got);
    }
    // Loads that cancel: the limit is 1e-6 of the largest load, 1 MPa on
    // 100 mm x 10 mm, not of the resultant 0.
    expect_lines("loads that cancel", model(MM, str(MESH, """
        let bottom: Surface = boundary(beam, z == 0 mm)
        let top: Surface = boundary(beam, z == 10 mm)

        step load: static {
            bc root: u = 0
            bc bottom: sigma*n = (0, 0, 1) MPa
            bc top: sigma*n = (0, 0, -1) MPa
        }
        """)), [
        "largest = max(largest, 1.0 * runtime::face_area(\"z0\"));\n    largest = max(largest, 1.0 * runtime::face_area(\"z1\"));\n    equilibrium(rx, ry, rz, largest, \"N\");",
    ]);
    expect_lines("fixed faces with an edge", model(MM, str(MESH, """
        let bottom: Surface = boundary(beam, z == 0 mm)

        step load: static {
            bc root: u = 0
            bc bottom: u = 0
        }
        """)), ["face(runtime::fix(\"x0\"), \"root\");\n    // bc bottom: u = 0\n    face(runtime::fix(\"z0\"), \"bottom\");",
                "let rx, ry, rz = runtime::total_reaction();"]);
    expect_lines("the reaction of a face with a shared edge", model(MM, str(MESH, """
        let bottom: Surface = boundary(beam, z == 0 mm)

        step load: static {
            bc root: u = 0
            bc bottom: u.x = 0
            bc tip: u = 0
            report reaction(root)
        }
        """)), ["let x, y, z = runtime::reaction(\"x0\");\n        var alone = true;\n        if shared(\"root\", \"x0\", \"bottom\", \"z0\") { alone = false; }\n        if shared(\"root\", \"x0\", \"tip\", \"x1\") { alone = false; }\n        if alone { report(\"load\", \"reaction(root)\", vector(x, y, z, \"N\")); }"]);
    expect_error("a later step fixes a new face", model(MM, str(MESH, """
        step one: static {
            bc root: u = 0
        }

        step two: static {
            bc tip: u = 0
        }
        """)), "error unsupported 24:5 step two fixes a face that step one does not fix");
    expect_error("the reaction of a loaded face", replace(model(MM, str(MESH, STEP)), "reaction(root)", "reaction(tip)"),
                 "error unsupported 23:21 the face tip is not fixed in step load");
}

// The dispatch tables: a template or a procedure of a lane that has no
// generator yet gives an error at its item.
fn test_dispatch() {
    expect_error("the heat template", replace(model(MM, str(MESH, STEP)), "linear_elastic", "heat"),
                 "error unsupported 16:9 the template heat has no generator in this build");
    expect_error("a transient step", replace(model(MM, str(MESH, STEP)), "load: static", "load: transient"),
                 "error unsupported 19:12 the procedure transient has no generator in this build");
    // The linear solver of a step: cg with jacobi, or a pair of the solvers lane.
    let solve = "    solve { linear: cg, precond: jacobi }\n";
    expect_lines("cg with jacobi", model(MM, str(MESH, replace(STEP, "    output u\n", solve))),
                 ["solve(\"load\", \"cg\", \"jacobi\");"]);
    expect_error("gmres", model(MM, str(MESH, replace(STEP, "    output u\n", replace(solve, "cg", "gmres")))),
                 "error unsupported 22:5 the linear solver gmres with the preconditioner jacobi is not in this build");
    expect_error("a solve setting", model(MM, str(MESH, replace(STEP, "    output u\n", replace(solve, "precond", "rtol")))),
                 "error unsupported 22:25 the solve setting is not in this slice");
}

fn test_unsupported() {
    expect_error("a cylinder", str(MM, """

        let bore: Solid = cylinder((40, 0, 20) mm, (0, 1, 0), 6 mm, 20 mm)
        """), "error unsupported 2:27 this slice has cylinders on the axis (0, 0, 1) only");
    expect_error("tet cells", replace(model(MM, str(MESH, STEP)), "cell: hex", "cell: tet"),
                 "error unsupported 13:11 this slice meshes with hex cells only");
    expect_error("order 1", model(MM, str(MESH, "method solid: cg { order: 1 }\n", STEP)),
                 "error unsupported 19:20 this slice has only order 2");
    expect_error("a face off the box", replace(model(MM, str(MESH, STEP)), "x == L", "x == 50 mm"),
                 "error unsupported 5:37 the surface is not on a face of beam");
    expect_error("no support", replace(model(MM, str(MESH, STEP)), "bc root: u = 0\n", ""),
                 "error unsupported 16:9 no face is fixed");
    expect_error("a unit outside the slice", replace(model(MM, str(MESH, STEP)), "-1) MPa", "-1) J"),
                 "error unsupported 21:34 the unit 'J' is not in this slice");
    expect_error("2d", model("model t: 2d, mm_N_MPa_s", str(MESH, STEP)), "error unsupported 1:10 this slice is 3d only");
}

fn main() {
    test_grid_cantilever();
    test_closed_cylinder();
    test_units();
    test_two_steps();
    test_bc_scope();
    test_reaction();
    test_dispatch();
    test_unsupported();
    if failures > 0 {
        print(failures, " codegen tests failed");
        exit(1);
    }
    print("codegen: all tests pass");
}
