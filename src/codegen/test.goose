// Tests of the code generator. `./build.sh test` runs this program from the
// root of the repository, and a non-zero exit status is a failure.
import std;
import os;
import syntax from "../syntax";
import .codegen;
import .study;
import .split;
import .status;
import check from "../check";   // study.goose writes TOON with it

extern "libc.so.6" {
    fn system(command: const u8&) -> i32;
}

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
    codegen::generate("RUNTIME", "RESULTS", codegen::no_share(), out);
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
        "runtime::mesh_box_incompatible(0.0, 0.0, 0.0, 100.0, 10.0, 10.0,",
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
    expect_error("incompatible of order 2", model(MM, str(MESH, "method solid: cg { order: 2, family: incompatible }\n", STEP)),
                 "error unsupported 19:38 this slice has the incompatible modes only for linear_elastic on hex cells of order 1");
    expect_error("serendipity", model(MM, str(MESH, "method solid: cg { family: serendipity }\n", STEP)),
                 "error unsupported 19:28 this slice has only the families lagrange and incompatible");
    expect_error("a face off the box", replace(model(MM, str(MESH, STEP)), "x == L", "x == 50 mm"),
                 "error unsupported 5:37 the surface is not on a face of beam");
    expect_error("no support", replace(model(MM, str(MESH, STEP)), "bc root: u = 0\n", ""),
                 "error unsupported 16:9 no face is fixed");
    expect_error("a unit outside the slice", replace(model(MM, str(MESH, STEP)), "-1) MPa", "-1) J"),
                 "error unsupported 21:34 the unit 'J' is not in this slice");
    expect_error("2d", model("model t: 2d, mm_N_MPa_s", str(MESH, STEP)), "error unsupported 1:10 this slice is 3d only");
}

// --- studies (spec §20) ---------------------------------------------------------

let STUDY_MODEL = """
    model t: 3d, SI

    param a: Length = 1 m in 0.5 m..4 m
    param b: Real = 2
    let c: Length = a

    step s: static {
        report r = reaction(c)
    }

    expect e: a == 1 m within 1%
    """;

// Reads the study `study` of `STUDY_MODEL` with the study item `item`. It
// gives the cases as lines of the varied values, or the first error.
fn cases(item: u8[:]) -> u8[>..] {
    let text = str(STUDY_MODEL, "\n", item);
    let prog = syntax::parse_model("t.skn", text);
    var out: u8[>..] = [];
    if prog >= 0 && codegen::study_cases(prog, codegen::find_study(prog, "s")) {
        for c in 0..codegen::case_count {
            for p in 0..codegen::varied.len { format(out, if p > 0 { "; " } else { "" }, codegen::case_value(c, p)); }
            out.push('\n');
        }
        return out;
    }
    let e = syntax::errors[0];
    let t = syntax::tokens[e.tok];
    format(out, "error ", e.id, " ", t.line, ":", t.col, " ", e.msg);
    return out;
}

fn expect_text(what: u8[:], got: u8[:], want: u8[:]) {
    if got != want {
        failures++;
        print("FAIL ", what, "\n--- got:\n", got, "\n--- want:\n", want);
    }
}

let SAMPLE_RANGES = ["1 m..2 m", "-2 m..-1 m"];
let SAMPLE_LOWS = [1.0, -2.0];

fn test_study_cases() {
    expect_text("grid", cases("study s: grid {\n    a: [1, 2] m\n    b: [3, 4, 5]\n}\n"),
                "1 m; 3\n1 m; 4\n1 m; 5\n2 m; 3\n2 m; 4\n2 m; 5\n");
    expect_text("zip", cases("study s: zip {\n    a: [1 m, 2*a]\n    b: [3, 4]\n}\n"), "1 m; 3\n2*a; 4\n");
    expect_text("zip lengths", cases("study s: zip {\n    a: [1, 2] m\n    b: [3]\n}\n"),
                "error procedure_args 14:5 the lists of a zip have one length: 'b' has 1 values, not 2");
    expect_text("a grid range", cases("study s: grid {\n    a: 1 m..2 m\n}\n"),
                "error unsupported 13:11 a grid or a zip takes a list of values, not a range");
    expect_text("minimize", cases("study s: minimize(objective: r, method: nelder_mead, max_cases: 9) {\n    a: 1 m..2 m\n}\n"),
                "error unsupported 12:10 the study procedure minimize is not in this build");
    expect_text("a sample list", cases("study s: sample(n: 2, method: lhs, seed: 1) {\n    a: [1, 2] m\n}\n"),
                "error unsupported 13:5 a sample takes a range for each param");
    expect_text("a sample method", cases("study s: sample(n: 2, method: grid, seed: 1) {\n    a: 1 m..2 m\n}\n"),
                "error procedure_args 12:10 sample needs n, method and seed");
    expect_text("a sample n", cases("study s: sample(n: 0, method: lhs, seed: 1) {\n    a: 1 m..2 m\n}\n"),
                "error procedure_args 12:20 n of sample is a whole number of 1 or more");
    expect_text("a sample unit", cases("study s: sample(n: 2, method: lhs, seed: 1) {\n    a: 500 mm..2 m\n}\n"),
                "error unsupported 13:5 a sample takes a range of two numbers in one unit");
    expect_text("an empty list", cases("study s: grid {\n    a: [] m\n}\n"),
                "error procedure_args 13:5 the list of 'a' has no values");
    expect_text("a param twice", cases("study s: grid {\n    a: [1] m\n    a: [2] m\n}\n"),
                "error study_param 14:5 the study varies 'a' twice");
    // A bare number in a list takes the unit of the list, also after a `-`.
    expect_text("a negative value", cases("study s: grid {\n    a: [-1, 2] m\n}\n"), "-1 m\n2 m\n");
    expect_text("a unit with a power", cases("study s: grid {\n    a: [1, 2] m^-1*m^2\n}\n"), "1 m^-1*m^2\n2 m^-1*m^2\n");
    // A case of a sample is the value lo + t*(hi - lo), with t in [0, 1), in
    // the unit of the range.
    for i in 0..SAMPLE_RANGES.len {
        let lo = SAMPLE_LOWS[i];
        let item = str("study s: sample(n: 4, method: lhs, seed: 7) {\n    a: ", SAMPLE_RANGES[i], "\n}\n");
        let got = cases(item);
        for line in split(got, "\n") {
            if line.len == 0 { continue; }
            let v, good = parse_flt(line[..max(line.len - 2, 0)]);
            if !(good && v >= lo && v < lo + 1.0 && ends_with(line, " m")) {
                failures++;
                print("FAIL a sample case of ", item, ": '", line, "'");
            }
        }
    }
    // lhs and sobol: one t in each quarter.
    for method in ["lhs", "sobol"] {
        var seen = [false, false, false, false];
        for c in 0..4 {
            let t = codegen::sample_points(4, 2, method == "sobol", 7)[c * 2 + 1];
            if t >= 0.0 && t < 1.0 { seen[(t * 4.0) as! i64] = true; }
        }
        if !(seen[0] && seen[1] && seen[2] && seen[3]) {
            failures++;
            print("FAIL ", method, ": the 4 points do not fill the 4 strata");
        }
    }
    // The text of a case: each varied default in parentheses.
    let text = str(STUDY_MODEL, "\nstudy s: grid {\n    b: [3]\n    a: [2] m\n}\n");
    let prog = syntax::parse_model("t.skn", text);
    codegen::study_cases(prog, codegen::find_study(prog, "s"));
    let case_text = codegen::case_source(text, 0);
    expect_text("the case text", case_text[..find(case_text, "let c")], "model t: 3d, SI\n\nparam a: Length = (2 m) in 0.5 m..4 m\nparam b: Real = (3)\n");
}

// The table of a study: a case that ran, a case whose expect fails, a case
// that did not build, and a case that printed its records and stopped with
// the exit status 1, for example when a solve does not converge.
fn test_study_table() {
    let text = str(STUDY_MODEL, "\nstudy s: grid {\n    a: [1, 2, 3, 4] m\n}\n");
    let prog = syntax::parse_model("t.skn", text);
    codegen::study_cases(prog, codegen::find_study(prog, "s"));
    let ran = """
        report[1]{step,expr,value}:
          s,r = reaction(c),"(1.0, 0.0, 0.0) N"
        verify[2]{stage,id,result,value,limit}:
          post,equilibrium,pass,0 N,1e-6 N
          post,expect,pass,"e: 1.0 m","1 m within 1%"
        """;
    codegen::outs.clear();
    codegen::outs.append(ran);
    codegen::outs.append(replace(ran, "post,expect,pass", "post,expect,fail"));
    var ends: codegen::CaseEnd[>..] = [];
    for c in 0..4 {
        var e = codegen::CaseEnd { ended: true, built: c != 2, status: if c == 3 { 1 } else { 0 }, seconds: 1.234, err: [], out_lo: 0, out_hi: 0 };
        let o = if c == 1 { 1 } else { 0 };
        if c != 2 {
            e.out_lo = o * ran.len;
            e.out_hi = (o + 1) * ran.len;
        } else {
            e.err.append("main.goose:3: error\n");
        }
        ends.push(e);
    }
    var out: u8[>..] = [];
    var toon: u8[>..] = [];
    let status = codegen::study_table("s", codegen::tok(codegen::find_study(prog, "s")), 2, prog, ends, false, out, toon);
    expect_text("the study table", out, str("""
        study: s
        jobs: 2
        cases[4]{case,a,status,verify,r,e,time_s}:
          1,1 m,ok,pass,"(1.0, 0.0, 0.0) N",pass,1.23
          2,2 m,ok,fail,"(1.0, 0.0, 0.0) N",fail,1.23
          3,3 m,failed,fail,"","",1.23
          4,4 m,ok,fail,"(1.0, 0.0, 0.0) N",pass,1.23
        diagnostics[1]{class,id,file,line,col,msg,fix}:
          error,case_failed,t.skn,12,7,"case 3 of the study s stopped","goose cannot build the program of the case: main.goose:3: error"
        summary: errors=1 lints=0
        """, "\n"));
    if status != 1 {
        failures++;
        print("FAIL the study table: exit status ", status, ", not 1");
    }
}

// The workers of a study (spec §20.5): k cases that fit in the memory, and
// max(1, H / k) threads for each.
fn test_study_workers() {
    let runs = [
        [4.0, 16.0, 1.0e6, 1.0e9, 4.0, 4.0],    // the memory holds 1000 cases
        [4.0, 16.0, 3.0e8, 1.0e9, 3.0, 5.0],    // the memory holds 3 cases
        [4.0, 16.0, 2.0e9, 1.0e9, 1.0, 16.0],   // one case is larger than the memory
        [8.0, 4.0, 0.0, 1.0e9, 8.0, 1.0],       // --jobs 8 on 4 threads
    ];
    for r in runs {
        let k, threads = codegen::study_workers(r[0] as! i64, r[1] as! i64, r[2], r[3]);
        if k != r[4] as! i64 || threads != r[5] as! i64 {
            failures++;
            print("FAIL study_workers(", r[0], ", ", r[1], ", ", r[2], ", ", r[3], "): ", k, " and ", threads, ", not ", r[4], " and ", r[5]);
        }
    }
}

// `skein run --study` of an example, and the command line errors. `./build.sh
// test` builds build/skein first. The runs write .skein under build/.
fn test_study_run() {
    let runs = [
        ["--study lengths --jobs 2", "0"],
        ["--study nothing", "2"],
        ["--jobs 2", "2"],
        ["--study lengths --jobs 0", "2"],
        ["--study sources --jobs 2", "1"],
        ["--study bands --jobs 2", "1"],
    ];
    for r in runs {
        var c: u8[>..] = [];
        format(c, "cd build && ./skein run ../examples/elements/heated_bar_study_tet.skn ", r[0], " >/dev/null 2>&1; test $? = ", r[1]);
        c.push(0);
        if system(c[0]) != 0 {
            failures++;
            print("FAIL skein run --study: '", r[0], "' does not exit with ", r[1]);
        }
    }
    var table: u8[>..] = [];
    read_file("build/.skein/studies/heated_bar_study_tet/lengths.toon", table);
    for want in ["jobs: 2\n", "cases[3]{case,lx,status,verify,heat_out,centre,time_s}:\n", "  1,0.005 m,ok,pass,-20.0", "  3,0.01 m,ok,pass,-40.0"] {
        if find(table, want) < 0 {
            failures++;
            print("FAIL skein run --study: no '", want, "' in:\n", table);
        }
    }
    // lx reaches the mesh, so each case meshes. qv reaches only the physics:
    // case 2 reads the mesh of case 1, and gives twice its heat.
    let studies = "build/.skein/studies/heated_bar_study_tet";
    if is_dir(str(studies, "/lengths/shared")) || !is_dir(str(studies, "/sources/shared")) {
        failures++;
        print("FAIL skein run --study: the study lengths shares stages, or the study sources shares none");
    }
    table.clear();
    read_file(str(studies, "/sources.toon"), table);
    for want in ["  1,1e7 W/m^3,ok,pass,-20.0", "  2,2e7 W/m^3,ok,fail,-40.0"] {
        if find(table, want) < 0 {
            failures++;
            print("FAIL skein run --study sources: no '", want, "' in:\n", table);
        }
    }
    // band reaches only the expect: case 2 reads the solution of case 1,
    // so the study solves once, and the tolerance of case 2 is too small.
    table.clear();
    read_file(str(studies, "/bands.toon"), table);
    for want in ["  1,0.0147 K,ok,pass,-20.0", "  2,1e-12 K,ok,fail,-20.0"] {
        if find(table, want) < 0 {
            failures++;
            print("FAIL skein run --study bands: no '", want, "' in:\n", table);
        }
    }
    // The solve record of a case (command_line.md §6) has one row for each solve.
    for c in ["1", "2"] {
        var record: u8[>..] = [];
        read_file(str(studies, "/bands/", c, "/run.toon"), record);
        let want = if c == "1" { "\nsolve[1]{" } else { "\nsolve[0]{" };
        if find(record, want) < 0 {
            failures++;
            print("FAIL skein run --study bands: case ", c, " has no '", want[1..], "' in:\n", record);
        }
    }
}

// The number of lines of the longest function of `program`, from its `fn`
// line to its closing brace.
fn longest_fn(program: u8[:]) -> i64 {
    var longest = 0;
    var start = -1;
    for line, i in split(program, '\n') {
        if starts_with(line, "fn ") { start = i; }
        if line == "}" && start >= 0 {
            longest = max(longest, i - start + 1);
            start = -1;
        }
    }
    return longest;
}

// `fn step_s() {` with the body lines `body`, after the splitter.
fn split_step(body: u8[:]) -> u8[>..] {
    var program = str("import std;\nfn step_s() {\n", body, "}\n");
    codegen::split_long_functions(program);
    return program;
}

// `count` statements `    call_<i>();`.
fn calls(from: i64, count: i64) -> u8[>..] {
    var out: u8[>..] = [];
    for i in from..from + count { format(out, "    call_", i, "();\n"); }
    return out;
}

fn expect_split(what: u8[:], got: u8[:], wants: const u8[:][:]) {
    if longest_fn(got) > 60 {
        failures++;
        print("FAIL split ", what, ": a function has ", longest_fn(got), " lines:\n", got);
    }
    for w in wants {
        if find(got, w) < 0 {
            failures++;
            print("FAIL split ", what, ": no '", w, "' in:\n", got);
            return;
        }
    }
}

fn expect_unsplit(what: u8[:], body: u8[:]) {
    let got = split_step(body);
    if find(got, "part1_") >= 0 {
        failures++;
        print("FAIL split ", what, ": it splits:\n", got);
    }
}

fn test_split() {
    let short = calls(0, 58);
    if split_step(short) != str("import std;\nfn step_s() {\n", short, "}\n") {
        failures++;
        print("FAIL split: a function of 60 lines changes");
    }
    // 130 statements: parts of 58, 58 and 14, called in order.
    expect_split("130 calls", split_step(calls(0, 130)), [
        "fn part1_step_s() {\n    call_0();\n",
        "    call_57();\n}\nfn part2_step_s() {\n    call_58();\n",
        "    call_115();\n}\nfn part3_step_s() {\n    call_116();\n",
        "    call_129();\n}\nfn step_s() {\n    part1_step_s();\n    part2_step_s();\n    part3_step_s();\n}\n",
    ]);
    // No cut between a local and its last use. A comment stays with the
    // statement below it, and a block is one statement.
    let locals = str(calls(0, 50), "    let rx, ry = reaction();\n", calls(50, 10),
                     "    // the check\n    {\n        check(rx);\n    }\n", calls(60, 40));
    expect_split("a local", split_step(locals), [
        "    call_49();\n}\nfn part2_step_s() {\n    let rx, ry = reaction();\n",
        "    // the check\n    {\n        check(rx);\n    }\n",
    ]);
    expect_unsplit("a local used at the end", str("    var a = 0;\n", calls(0, 70), "    use(a);\n"));
    expect_unsplit("a statement of 60 lines", str(calls(0, 10), "    {\n", calls(10, 58), "    }\n", calls(68, 10)));
    expect_unsplit("a return", str(calls(0, 40), "    if done { return; }\n", calls(40, 40)));
}

fn test_unchecked_status() {
    let checked = "fn solve() {\n    converged = runtime::solve_static();\n    let ok = shells::solve();\n    if !elements::solve() { failed = true; }\n}\n";
    if codegen::unchecked_status(checked) != 0 {
        failures++;
        print("FAIL unchecked status: a checked solver call is line ", codegen::unchecked_status(checked));
    }
    for call in ["runtime::solve_static();", "elements::solve();", "shells::solve();"] {
        let program = str("fn solve() {\n    report();\n    ", call, "\n}\n");
        if codegen::unchecked_status(program) != 3 {
            failures++;
            print("FAIL unchecked status: '", call, "' is line ", codegen::unchecked_status(program), ", not 3");
        }
    }
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
    test_study_cases();
    test_study_table();
    test_study_workers();
    test_study_run();
    test_split();
    test_unchecked_status();
    if failures > 0 {
        print(failures, " codegen tests failed");
        exit(1);
    }
    print("codegen: all tests pass");
}
