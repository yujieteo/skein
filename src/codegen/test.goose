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
        "face(runtime::fix_face(runtime::X, 0.0), \"root\");",
        "face(runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -1.0), \"tip\");",
        "equilibrium(rx, ry, rz, 100.0, \"N\");",
        "report(\"load\", \"reaction(root)\", vector(x, y, z, \"N\"));",
        "output(\"RESULTS/load.vtu\");",
        "let ux, uy, uz = runtime::displacement_at(100.0, 5.0, 5.0);",
        "expect(\"tip_uz\", uz, -0.191933, 0.015, true, \"mm\", \"-0.191933 mm within 1.5%\");",
        "not_run(\"top_sxx: sigma.xx at (51.25, 5, 10) mm is not in this slice\");",
    ]);
}

fn test_units() {
    // 1 kN on 50 mm x 100 mm is 0.2 MPa, and in SI 200000 Pa.
    let load = """
        step load: static {
            bc root: u = 0
            bc tip: sigma*n = (0, 0, -1 kN/(50 mm*100 mm))
            report tip_uz = u.z at (1, 0.025, 0.05) m
        }
        """;
    expect_lines("units in mm_N_MPa_s", model(MM, str(MESH, load)), [
        "traction_face(runtime::X, 100.0, 0.0, 0.0, -0.2)",
        "displacement_at(1000.0, 25.0, 50.0)",
        "quantity(uz, \"mm\")",
    ]);
    expect_lines("units in SI", model("model t: 3d, SI", str(MESH, load)), [
        "runtime::mesh_box(0.0, 0.0, 0.0, 0.1, 0.01, 0.01,",
        "runtime::material(210000000000.0, 0.3);",
        "traction_face(runtime::X, 0.1, 0.0, 0.0, -199999.9999",
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
        "fn step_two() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -2.0) MPa\n    face(runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -2.0), \"tip\");",
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
        "fn step_load() {\n    // bc tip: sigma*n = (0.0, 0.0, -1.0) MPa\n    face(runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -1.0), \"tip\");",
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
        "fn step_two() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -1.0) MPa\n    face(runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -1.0), \"tip\");\n    // bc bottom: sigma*n = (0.0, 0.0, 1.0) MPa\n    face(runtime::traction_face(runtime::Z, 0.0, 0.0, 0.0, 1.0), \"bottom\");\n    solve(",
        "fn step_three() {\n    // The supports of the first step stay.\n    runtime::clear_loads();\n    // bc tip: sigma*n = (0.0, 0.0, -3.0) MPa\n    face(runtime::traction_face(runtime::X, 100.0, 0.0, 0.0, -3.0), \"tip\");\n    // bc bottom:",
    ]);
}

// The equilibrium check sums the reaction of each fixed face once.
fn test_reaction() {
    let twice = """
        bc root: u = 0

        step load: static {
            bc root: u = 0
            bc tip: sigma*n = (0, 0, -1) MPa
        }
        """;
    let got = generate(model(MM, str(MESH, twice)));
    if occurrences(got, "runtime::reaction_face(") != 1 || occurrences(got, "fix_face(") != 1 {
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
        """)), ["equilibrium(rx, ry, rz, 1000.0, \"N\");"]);
    expect_error("fixed faces with an edge", model(MM, str(MESH, """
        let bottom: Surface = boundary(beam, z == 0 mm)

        step load: static {
            bc root: u = 0
            bc bottom: u = 0
        }
        """)), "error unsupported 23:5 the fixed faces root and bottom share an edge");
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

fn test_unsupported() {
    expect_error("a cylinder", str(MM, """

        let bore: Solid = cylinder((40, 0, 20) mm, (0, 1, 0), 6 mm, 20 mm)
        """), "error unsupported 2:27 the solid 'bore' is not one box");
    expect_error("tet cells", replace(model(MM, str(MESH, STEP)), "cell: hex", "cell: tet"),
                 "error unsupported 13:11 this slice meshes a box with hex cells only");
    expect_error("order 1", model(MM, str(MESH, "method solid: cg { order: 1 }\n", STEP)),
                 "error unsupported 19:20 this slice has only order 2");
    expect_error("a face off the box", replace(model(MM, str(MESH, STEP)), "x == L", "x == 50 mm"),
                 "error unsupported 5:37 the surface is not on a face of the box");
    expect_error("no support", replace(model(MM, str(MESH, STEP)), "bc root: u = 0\n", ""),
                 "error unsupported 16:9 no face is fixed");
    expect_error("a unit outside the slice", replace(model(MM, str(MESH, STEP)), "-1) MPa", "-1) J"),
                 "error unsupported 21:34 the unit 'J' is not in this slice");
    expect_error("2d", model("model t: 2d, mm_N_MPa_s", str(MESH, STEP)), "error unsupported 1:10 this slice is 3d only");
}

// The dispatch tables: a template or a procedure of a lane that has no
// generator yet gives an error at its item.
fn test_dispatch() {
    expect_error("the heat template", replace(model(MM, str(MESH, STEP)), "linear_elastic", "heat"),
                 "error unsupported 16:9 the template heat has no generator in this build");
    expect_error("a transient step", replace(model(MM, str(MESH, STEP)), "load: static", "load: transient"),
                 "error unsupported 19:12 the procedure transient has no generator in this build");
}

fn main() {
    test_grid_cantilever();
    test_dispatch();
    test_units();
    test_two_steps();
    test_bc_scope();
    test_reaction();
    test_unsupported();
    if failures > 0 {
        print(failures, " codegen tests failed");
        exit(1);
    }
    print("codegen: all tests pass");
}
