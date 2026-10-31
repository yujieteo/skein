// Tests of `skein check`. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
import std;
import os;
import syntax from "../syntax";
import .state;
import .check;

var failures = 0;

let H = "model m: 3d, SI\nscales 1 m, 1 kg, 1 s, 1 K\n";
let H2 = "model m: 2d, SI\nscales 1 m, 1 kg, 1 s, 1 K\n";
let SOLID = """
    let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)
    let f: Surface = boundary(b, x == 0 m)
    material steel {
        E: Pressure = 200 GPa
        nu: Real = 0.3
        rho: Density = 7800 kg/m^3
    }
    mesh m1 on b { size: 0.1 m }
    """;
let FLOW = """
    let c: Surface = rect((0, 0) m, (1, 1) m)
    let lid: Curve = boundary(c, y == 1 m)
    material fluid {
        rho: Density = 1 kg/m^3
        mu: DynViscosity = 0.01 Pa*s
    }
    mesh m1 on c { size: 0.1 m }
    physics flow on c: navier_stokes solves v: Vector<Velocity>, p: Pressure {
        material fluid
    }
    """;

// A shell with ribs: the templates of G3 and N1 (spec §4.3), a tie, an edge
// load N*n, argmax and max on a pattern (N2), and quad cells on a Surface (N3).
let SHELL = """

    let barrel: Solid = cylinder((0, 0, 0) m, (1, 0, 0), 1 m, 2 m)
    let wall: Surface = boundary(barrel, y^2 + z^2 == 1 m^2)
    let root: Curve = boundary(wall, x == 0 m)
    let tip: Curve = boundary(wall, x == 2 m)
    let s0: Curve = line((0, 0, 1) m, (2, 0, 1) m)
    let ribs: Curve = pattern_polar(s0, (1, 0, 0), 8)

    table plies{angle: Angle in deg}:
        0
        90

    material lam {
        E1: Pressure = 148 GPa
        E2: Pressure = 10 GPa
        G12: Pressure = 4 GPa
        nu12: Real = 0.3
        t_ply: Length = 1 mm
        h: Length = 4 mm
    }

    material sec {
        E: Pressure = 70 GPa
        nu: Real = 0.3
        A: Area = 1e-4 m^2
        I_y: Length^4 = 1e-8 m^4
        I_z: Length^4 = 2e-8 m^4
        J: Length^4 = 1e-9 m^4
    }

    mesh m1 on wall | ribs {
        size: 0.1 m
        cell: quad
    }

    physics skin on wall: reissner_mindlin(layup: plies.angle, mirror: true, top: n, q: 1 kPa*n) solves w: Vector<Length>,
        theta: Vector<Angle> {
        material lam
    }

    physics rib on ribs: timoshenko(A: A, I_y: I_y, I_z: I_z, J: J, y_axis: n, offset: -10 mm*n) solves w_b: Vector<Length>,
        theta_b: Vector<Angle> {
        material sec
    }

    method skin: cg { mesh: m1, order: 1 }

    method rib: cg { mesh: m1, order: 1 }

    bc ribs: w_b = w
    bc ribs: theta_b = theta
    bc root: w = 0
    bc root: theta = 0

    step go: static {
        bc tip: N*n = (1, 0, 0) kN/m
        report peak = argmax(norm(w), wall)
        report rib_max = max(Nx, ribs)
    }
    """;

// The ids of the diagnostics of `src`, one per line.
fn ids(src: u8[:]) -> u8[>..] {
    check::check_text("t.skn", src);
    var out: u8[>..] = [];
    for d in check::diags {
        out.append(d.id);
        out.push('\n');
    }
    return out;
}

// `src` gives the error `id` (and perhaps others).
fn expect_id(id: u8[:], src: u8[:]) {
    let got = ids(src);
    for line in split(got, '\n') {
        if line == id { return; }
    }
    failures++;
    print("FAIL ", id, ": got [", got, "]");
}

// The free boundary report of `src`, one "field,region,condition" per line.
fn free_rows(src: u8[:]) -> u8[>..] {
    check::check_text("t.skn", src);
    var out: u8[>..] = [];
    for f in check::free {
        out.append(f.field);
        out.push(',');
        out.append(f.region);
        out.push(',');
        out.append(f.condition);
        out.push('\n');
    }
    return out;
}

fn expect_free(what: u8[:], src: u8[:], want: u8[:]) {
    let got = free_rows(src);
    if got == want && check::diags.len == 0 { return; }
    failures++;
    print("FAIL ", what, ": got [", got, "] and ", check::diags.len, " diagnostics");
}

// `src` gives exactly the ids `want`, one per line.
fn expect_ids(what: u8[:], src: u8[:], want: u8[:]) {
    let got = ids(src);
    if got == want { return; }
    failures++;
    print("FAIL ", what, ": got [", got, "]");
}

// `src` gives no diagnostic.
fn expect_clean(what: u8[:], src: u8[:]) {
    let got = ids(src);
    if got.len == 0 { return; }
    failures++;
    print("FAIL ", what, ": got [", got, "]");
}

fn test_names() {
    expect_id("syntax", "model m: 3d, SI\nlet a: Real = (1, 2\n");
    expect_id("header_order", "let a: Real = 1\n");
    expect_id("unknown_name", str(H, "let a: Real = bb\n"));
    expect_id("use_before_decl", str(H, "let a: Real = bb\nlet bb: Real = 1\n"));
    expect_id("duplicate_name", str(H, "let a: Real = 1\nlet a: Real = 2\n"));
    expect_id("duplicate_name", str(H, "let t: Time = 1 s\n"));
    expect_id("unknown_unit", str(H, "let a: Length = 5 mms\n"));
    expect_id("table_row_width", str(H, "table tb{a: Length in m, b: Real}:\n    1, 2, 3\n"));
}

// `use` library files (spec §4.1). The test writes them into build/.
fn test_use() {
    write_file("build/use_ok.skn", """
        // The steel of the tests.
        let E_s: Pressure = 200 GPa
        param rho_s: Density = 7800 kg/m^3
        fn twice(a: Length) -> Length = 2*a
        material steel {
            E: Pressure = E_s
            nu: Real = 0.3
            rho: Density = rho_s
        }
        """);
    write_file("build/use_header.skn", "model lib: 3d, SI\nlet a: Real = 1\n");
    write_file("build/use_table.skn", "let a: Real = 1\ntable tb{a: Real}:\n    1\n");
    write_file("build/use_use.skn", "use \"use_ok.skn\"\n");
    write_file("build/use_bad.skn", "// 1\n// 2\n// 3\nlet a: Real = bb\n");
    expect_clean("use", str(H, "use \"build/use_ok.skn\"\nlet L: Length = twice(1 m)\nlet E2: Pressure = 2*E_s\n"));
    expect_id("use_item", str(H, "use \"build/no_such_file.skn\"\n"));
    expect_id("use_item", str(H, "use \"build/use_header.skn\"\n"));
    expect_id("use_item", str(H, "use \"build/use_table.skn\"\n"));
    expect_id("use_item", str(H, "use \"build/use_use.skn\"\n"));
    expect_id("duplicate_name", str(H, "use \"build/use_ok.skn\"\nlet E_s: Pressure = 1 GPa\n"));
    expect_id("duplicate_name", str(H, "use \"build/use_ok.skn\"\nuse \"build/use_ok.skn\"\n"));
    // An error in the library and one in the model at the same line and column.
    expect_ids("a library and a model error at one place",
               str(H, "use \"build/use_bad.skn\"\nlet c: Real = dd\n"), "unknown_name\nunknown_name\n");
    if check::diags.len != 2 || check::diags[0].file != "t.skn" || check::diags[1].file != "build/use_bad.skn" {
        failures++;
        print("FAIL the model diagnostics come before the library diagnostics");
    }
}

fn test_types() {
    expect_id("type_mismatch", str(H, "let count: Int = 2.5\n"));
    expect_id("unit_mismatch", str(H, "let a: Pressure = 1 GPa + 1 kg/m^3\n"));
    expect_id("rank_mismatch", str(H, "let a: Length = (1, 2, 3) m + 1 m\n"));
    expect_id("ambiguous_product", str(H, "let a: Vector<Length> = (1, 2, 3) m\nlet c: Real = a*a\n"));
    expect_id("unit_transcendental", str(H, "let a: Temperature = 300 K\nlet c: Real = exp(a)\n"));
    expect_id("unit_affine", str(H, "let k: ThermalConductivity = 1.2 W/(m*degC)\n"));
    expect_id("not_constant", str(H, "param a: Length = x\n"));
    // sum over the rows of a table (spec §5.6).
    let tb = "table tb{a: Length in m, b: Real}:\n    1, 2\n    3, 4\n";
    expect_clean("sum of a table column", str(H, tb, "let s: Area = sum(tb.a^2*tb.b)\n"));
    expect_id("unit_mismatch", str(H, tb, "let s: Real = sum(tb.a)\n"));
    expect_id("sum_no_column", str(H, "let s: Length = sum(2 m)\n"));
    let tc = "table tc{c: Real}:\n    1\n";
    expect_id("sum_mixed_tables", str(H, tb, tc, "let s: Length = sum(tb.a*tc.c)\n"));
    expect_id("sum_mixed_tables", str(H, tb, tc, "let l: Length = tb.a\nlet s: Length = sum(l*tc.c)\n"));
    expect_clean("sum of interp at the rows of another table", str(H, tb, tc, "let s: Length = sum(tc.c*interp(tb.b, tb.a, tc.c))\n"));
    expect_id("sum_no_column", str(H, tb, "let s: Length = sum(interp(tb.b, tb.a, 2))\n"));
    expect_clean("sum of one table through a let", str(H, tb, "let l: Length = tb.a\nlet s: Area = sum(l*tb.a)\n"));
    expect_id("geometry_dim", str(H2, "let s: Solid = rect((0, 0) m, (1, 1) m)\n"));
    expect_id("geometry_space", str(H, "let s: Surface = circle((0, 0) m, 1 m)\n"));
    expect_clean("a 2D profile of extrude and revolve in a 3d model", str(H, "let s: Solid = extrude(rect((0, 0) m, (1, 1) m), 2 m) | revolve(circle((2, 0) m, 1 m), 1 rad)\n"));
    expect_id("geometry_space", str(H, "let s: Solid = extrude(rect((0, 0) m, (1, 1) m), circle((0, 0) m, 1 m))\n"));
    expect_clean("loft, pattern and pattern_polar", str(H, "let s: Solid = loft(circle((0, 0) m, 2 m), rect((-1, -1) m, (1, 1) m), 3 m) | pattern(box((0, 0, 0) m, (1, 1, 1) m), (2, 0, 0) m, 4) | pattern_polar(sphere((5, 0, 0) m, 1 m), (0, 0, 1), 6)\n"));
    expect_id("type_mismatch", str(H, "let s: Solid = loft(circle((0, 0) m, 2 m), 3 m)\n"));
    expect_id("type_mismatch", str(H, "let s: Solid = loft(circle((0, 0) m, 2 m), 1 m, 3 m)\n"));
    expect_id("type_mismatch", str(H, "let s: Solid = pattern(box((0, 0, 0) m, (1, 1, 1) m), (2, 0, 0) m)\n"));
    expect_id("nondim_unit", "model m: 2d, nondim\nlet a: Real = 1 m\n");
    // Gaps G7, G12, G29, G30: these are valid.
    expect_clean("G7 re, imag, arg, conj", str(H, "let w: Complex<Real> = 1 + 2*im\nlet a: Real = re(w) + imag(w) + abs(w)\nlet b: Angle = arg(conj(w))\n"));
    expect_clean("G12 atan, atan2", str(H, "let a: Angle = atan(0.5) + atan2(1 m, 2 m)\n"));
    expect_clean("G29 a negative degC", str(H, "let a: Temperature = -5 degC\n"));
    expect_clean("a vector with a unit, a number and a name", str(H, "let h: Length = 1 m\nlet p: Vector<Length> = (0.5, h, 0) m\n"));
    expect_clean("G30 a real exponent", str(H, "param q: Real = 3.0\nlet a: Real = 0.5^(2/q)\n"));
    // Gap G14: exponents in halves (spec §6.1).
    expect_clean("G14 exponents in halves", str(H, "let a: Length = 25 mm\nlet k: Pressure*Length^0.5 = 1 MPa*sqrt(a)\n",
                                                 "let j: Pressure*Length^0.5 = 1051.27 MPa*mm^0.5\nlet c: Length^-1.5 = 1/a^1.5\n"));
    expect_id("unit_fractional", str(H, "let a: Length = 1 m\nlet c: Real = sqrt(sqrt(a))\n"));
    expect_id("unit_fractional", str(H, "let a: Length = 1 m\nlet c: Real = a^0.25\n"));
    expect_id("unit_fractional", str(H, "let a: Length = 1 m\nlet c: Real = sqrt(a)^0.5\n"));
}

fn test_physics() {
    // skein help std prints the templates, with no model line.
    let std = check::std_text();
    let want = "template physics linear_elastic(b: Vector<Force/Volume> = 0) solves u: Vector<Length> {";
    if std.len < want.len || std[..want.len] != want {
        failures++;
        print("FAIL help std: the text does not start with the linear_elastic template");
    }
    let elastic = str(H, SOLID, "\n");
    expect_id("template_eq", str(elastic, "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    eq div(u) = 0\n}\n"));
    expect_id("eq_count", str(elastic, "physics s on b solves u: Vector<Length>, q: Pressure {\n    material steel\n    eq div(hooke(E, nu, strain(u))) = 0\n}\n"));
    expect_id("material_cover", str(elastic, "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel, steel\n}\n"));
    expect_id("material_prop_missing", str(elastic, "material soft {\n    E: Pressure = 1 GPa\n    nu: Real = 0.3\n}\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel on f, soft on b - f\n    let w: Density = rho\n}\n"));
    expect_id("material_prop_missing", str(elastic, "material soft {\n    E: Pressure = 1 GPa\n    nu: Real = 0.3\n}\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material soft\n    let w: Density = rho\n}\n"));
    expect_id("unknown_name", str(elastic, "let w: Density = rho\n"));
    expect_id("material_prop_type", str(elastic, "material soft {\n    E: Pressure = 1 GPa\n    nu: Length = 1 m\n}\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel on f, soft on b - f\n}\n"));
    expect_id("assumption_unknown", str(elastic, "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plain_stress\n}\n"));
    expect_id("assumption_conflict", str(elastic, "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume finite_strain\n}\n"));
    expect_id("assumption_space", str(elastic, "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plane_stress(1 mm)\n}\n"));
    expect_id("assumption_invalid", str(elastic, "physics s on b solves u: Vector<Length> {\n    material steel\n    assume incompressible\n    eq div(hooke(E, nu, strain(u))) = 0\n}\n"));
}

fn test_conditions() {
    let elastic = str(H, SOLID, "\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n");
    expect_id("bc_form", str(elastic, "bc f: sigma = 0\n"));
    expect_ids("a bc with a property and no flux", str(elastic, "bc f: E*n = 0\n"), "bc_form\n");
    expect_id("bc_conflict", str(elastic, "bc f: u = 0\nbc f: sigma*n = (0, 0, 1) MPa\n"));
    expect_id("bc_off_boundary", str(elastic, "bc b: u = 0\n"));
    expect_clean("G16 a point force", str(elastic, "let tip: Point = point((1, 1, 1) m)\nbc tip: sigma*n = (0, 0, -1) kN\n"));
    expect_clean("a bc on a face", str(elastic, "bc f: u = 0\nbc f: u.x = 0 m\n"));
    // Gap G17: a rigid face, with the force and the moment on its point.
    let tip = "let c: Point = point((1, 1, 1) m)\nbc f: u = rigid(c)\n";
    expect_clean("G17 a rigid face with a force and a moment", str(elastic, tip, "bc c: sigma*n = (0, 0, -1) kN\nbc c: moment(sigma*n) = (0, 1, 0) kN*m\n"));
    expect_clean("G17 a fixed component and a force with 0 in it", str(elastic, tip, "bc c: u.x = 0\nbc c: sigma*n = (0, 0, -1) kN\n"));
    expect_clean("G17 a fixed point with a moment", str(elastic, tip, "bc c: u = 0\nbc c: moment(sigma*n) = (0, 1, 0) kN*m\n"));
    expect_id("bc_conflict", str(elastic, tip, "bc c: u.x = 0\nbc c: sigma*n = (1, 0, -1) kN\n"));
    expect_id("bc_conflict", str(elastic, "let c: Point = point((1, 1, 1) m)\nbc c: u.x = 0\nbc c: sigma*n = (0, 0, -1) kN\n"));
    expect_id("unit_mismatch", str(elastic, tip, "bc c: moment(sigma*n) = (0, 1, 0) kN\n"));
    expect_id("rigid_form", str(elastic, "let c: Point = point((1, 1, 1) m)\nbc f: u.x = rigid(c)\n"));
    expect_id("rigid_form", str(elastic, "bc f: u = rigid(f)\n"));
    expect_id("rigid_form", str(elastic, "let c: Point = point((1, 1, 1) m)\nbc f: u = 2*rigid(c)\n"));
    expect_id("rigid_form", str(elastic, "let c: Point = point((1, 1, 1) m)\nbc c: moment(sigma*n) = (0, 1, 0) kN*m\n"));
    expect_id("rigid_form", str(elastic, tip, "bc f: moment(sigma*n) = (0, 1, 0) kN*m\n"));
    expect_ids("G17 a 2D moment is a scalar", "model m: 2d, mm_N_MPa_s\nlet b: Surface = rect((0, 0) mm, (10, 10) mm)\nlet f: Curve = boundary(b, y == 10 mm)\nlet g: Curve = boundary(b, y == 0 mm)\nlet c: Point = point((5, 10) mm)\nscales 10 mm, 210 GPa\nmaterial steel {\n    E: Pressure = 210 GPa\n    nu: Real = 0.3\n}\nmesh m1 on b { size: 1 mm }\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plane_strain(1 mm)\n}\nbc g: u = 0\nbc f: u = rigid(c)\nbc c: sigma*n = (1, -2) N\nbc c: moment(sigma*n) = 5 N*mm\n", "");
    expect_clean("a vector and a component bc that agree", str(elastic, "bc f: u = (0, 1, 0) mm\nbc f: u.x = 0 m\n"));
    expect_clean("a named vector and a component bc", str(elastic, "let d: Vector<Length> = (1, 0, 0) mm\nbc f: u = d\nbc f: u.y = 0 m\n"));
    expect_id("bc_conflict", str(elastic, "bc f: u.x = 0 m\nbc f: sigma*n = (0, 0, 1) MPa\n"));
    expect_id("bc_conflict", str(elastic, "bc f: u = (1, 0, 0) mm\nbc f: u.x = 0 m\n"));
    expect_id("bc_conflict", str(elastic, "bc boundary(b, x == 1 m): u = 0\nbc boundary(b,x==1 m): sigma*n = (0, 0, 1) MPa\n"));
    let heat = str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nlet f: Surface = boundary(b, x == 0 m)\nmaterial cu {\n    k: Conductivity = 400 W/(m*K)\n}\nmesh m1 on b { size: 0.1 m }\n");
    expect_clean("a flux with other spaces", str(heat, "physics h on b solves T: Temperature {\n    material cu\n    eq div(k * grad(T) + k * grad(T) + k * grad(T) + k * grad(T) + k * grad(T)) = 0 W/m^3\n}\nbc f: dot(k*grad(T) + k*grad(T) + k*grad(T) + k*grad(T) + k*grad(T), n) = 0 W/m^2\n"));
    expect_free("a sum flux", str(heat, "physics h on b solves T: Temperature {\n    material cu\n    eq div(k * grad(T) + k * grad(T)) = 0 W/m^3\n}\nbc f: T = 300 K\n"),
                "T,boundary(b) - f,(k * grad(T) + k * grad(T))*n = 0\n");
    expect_free("a bc on a compound region", str(heat, "physics h on b solves T: Temperature {\n    material cu\n    eq div(k * grad(T)) = 0 W/m^3\n}\nbc boundary(b)-f: T = 300 K\n"),
                "T,boundary(b) - (boundary(b)-f),(k * grad(T))*n = 0\n");
    expect_free("a symmetry plane", str(elastic, "let plane: Surface = boundary(b, y == 0 m)\nbc f: u = 0\nbc plane: u.y = 0 m\n"),
                "u,boundary(b) - f - plane,sigma*n = 0\nu.x,plane,(sigma*n).x = 0\nu.z,plane,(sigma*n).z = 0\n");
    expect_id("method_missing", str(H, SOLID, "\nmesh m2 on b { size: 0.2 m }\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n"));
    expect_id("method_setting", str(elastic, "method s: cg {\n    points: m1\n}\n"));
    expect_id("method_physics", str(H, SOLID, "\npoints cloud on b { spacing: 0.1 m }\nphysics s on b solves u: Vector<Length> {\n    material steel\n    state a: Real = 0\n    let sigma: Tensor2<Pressure> = hooke(E, nu, strain(u))\n    plastic {\n        yield: von_mises(sigma) - 200 MPa\n        strain: a\n    }\n    eq div(sigma) = 0\n}\nmethod s: gfd {\n    points: cloud\n}\n"));
    expect_id("inf_sup", str(H2, FLOW, "\nmethod flow: cg {\n    mesh: m1\n    order: { v: 1, p: 1 }\n    stabilize: none\n}\n"));
    expect_clean("a stable pair", str(H2, FLOW, "\nmethod flow: cg {\n    mesh: m1\n    order: { v: 1, p: 1 }\n    stabilize: supg_pspg\n}\n"));
    expect_id("order_too_low", str(H2, "let c: Surface = rect((0, 0) m, (1, 1) m)\nmesh m1 on c { size: 0.1 m }\nphysics plate on c solves w: Length {\n    eq lap(lap(w)) = 0 m^-3\n}\nmethod plate: cg {\n    mesh: m1\n    order: 1\n}\n"));
}

// Gap G20: quarter-point cells at a crack tip (spec §9.1).
fn test_crack_tip() {
    let plate = "let b: Surface = rect((0, 0) mm, (10, 10) mm)\nlet f: Curve = boundary(b, y == 0 mm)\nlet c: Point = point((5, 0) mm)\nscales 10 mm, 210 GPa\nmaterial steel {\n    E: Pressure = 210 GPa\n    nu: Real = 0.3\n}\n";
    let solid = "physics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plane_strain(1 mm)\n}\nbc f: u = 0\n";
    let model = "model m: 2d, mm_N_MPa_s\n";
    expect_clean("G20 a crack tip", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2 mm\n}\n", solid));
    expect_id("crack_tip_form", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip f: 0.2 mm\n}\n", solid));
    expect_id("crack_tip_form", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2\n}\n", solid));
    expect_id("crack_tip_form", str(model, plate, "points p1 on b {\n    spacing: 1 mm\n    crack_tip c: 0.2 mm\n}\nmesh m1 on b { size: 1 mm }\n", solid));
    expect_id("crack_tip_form", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2 mm\n}\n", solid, "method s: cg { mesh: m1, order: 1 }\n"));
    expect_id("crack_tip_form", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2 mm\n}\n", solid, "method s: cg { order: 1 }\n"));
    expect_clean("G20 a meshless method beside a crack tip mesh", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2 mm\n}\npoints cloud on b { spacing: 1 mm }\n", solid, "method s: rkpm { points: cloud }\n"));
    expect_id("crack_tip_form", str(H, SOLID, "\nlet c: Point = point((1, 1, 1) m)\nmesh m2 on b {\n    size: 0.1 m\n    crack_tip c: 0.01 m\n}\n"));
    expect_id("not_constant", str(model, plate, "mesh m1 on b {\n    size: 1 mm\n    crack_tip c: 0.2 mm*(1 + x/(1 mm))\n}\n", solid));
}

fn test_steps() {
    let elastic = str(H, SOLID, "\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n");
    expect_id("procedure_args", str(elastic, "step go: transient(end: 1 s) {\n    bc f: u = 0\n}\n"));
    expect_id("procedure_args", str(elastic, "step go: static(rate: 2) {\n    bc f: u = 0\n}\n"));
    expect_id("not_constant", str(elastic, "step go: transient(end: t, dt: 1 s) {\n    bc f: u = 0\n}\n"));
    expect_clean("G9 a list of frequencies", str(elastic, "step go: harmonic(freq: [10, 20] Hz) {\n    bc f: u = 0\n}\n"));
    expect_clean("G6 an eigenvalue", str(elastic, "step go: modal(count: 4) {\n    bc f: u = 0\n    report f1 = eigenvalue(1)\n}\nexpect f1: eigenvalue(1) == 10 Hz within 1%\n"));
    expect_clean("an eigenvalue with no units", "model m: 3d, nondim\nlet b: Solid = box((0, 0, 0), (1, 1, 1))\nlet f: Surface = boundary(b, x == 0)\nmaterial steel {\n    E: Pressure = 1000\n    nu: Real = 0.3\n    rho: Density = 1\n}\nmesh m1 on b { size: 0.1 }\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\nstep go: modal(count: 4) {\n    bc f: u = 0\n    report f1 = eigenvalue(1)\n}\nexpect f1: eigenvalue(1) == 2.5 within 1%\n");
    expect_id("type_mismatch", str(elastic, "step go: static {\n    bc f: u = 0\n    report eigenvalue(1)\n}\n"));
    // when: the first time at which a value gets to a level (spec §14.5).
    let shake = "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report t1 = when(max(norm(u), b), 1 mm)\n}\n";
    expect_clean("when in a report of a transient step", str(elastic, shake));
    expect_clean("when in an expect after a transient step", str(elastic, shake, "expect t1: when(max(norm(u), b), 2 mm) - when(max(norm(u), b), 1 mm) == 0.1 s within 5%\n"));
    expect_clean("when at each row of a table", str(elastic, "table lv{d: Length in mm, t: Time in s}:\n    1, 0.1\n    2, 0.2\n", shake, "expect t2: when(max(norm(u), b), lv.d) == lv.t within 5%\n"));
    expect_id("when_step", str(elastic, "step go: static {\n    bc f: u = 0\n    report when(max(norm(u), b), 1 mm)\n}\n"));
    expect_id("when_step", str(elastic, "step go: static {\n    bc f: u = 0\n}\nexpect t1: when(max(norm(u), b), 1 mm) == 0.1 s within 5%\n"));
    expect_id("when_step", str(elastic, shake, "let t2: Time = when(max(norm(u), b), 1 mm)\n"));
    expect_id("unit_mismatch", str(elastic, "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report when(max(norm(u), b), 1 s)\n}\n"));
    expect_id("rank_mismatch", str(elastic, "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report when(max(u, b), 1 mm)\n}\n"));
    expect_id("rank_mismatch", str(elastic, "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report when(norm(u), 1 mm)\n}\n"));
    expect_clean("when of a scalar reduced over a region", str(elastic, "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report when(min(norm(u), b), 1 mm)\n}\n"));
    expect_id("not_constant", str(elastic, "step go: transient(end: 1 s, dt: 1 ms) {\n    bc f: u = 0\n    report when(max(norm(u), b), x)\n}\n"));
    expect_id("allow_reason", str(H, "allow rigid_body: \"\"\n"));
    expect_id("allow_reason", str(H, "allow rigid_bodies: \"a reason\"\n"));
    expect_id("results_missing", str(H, "let g: Results = results(\"no_such_file.skn\", step: \"load\")\n"));
    expect_id("results_missing", str(H, "let g: Results = results(\"samples/01_cantilever.skn\", step: \"unload\")\n"));
    expect_id("expect_form", str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nexpect e: b == 1 m within 1%\n"));
    expect_id("unit_mismatch", str(H, "let a: Length = 1 m\nexpect e: a == 1 m within 1 s\n"));
}

fn test_studies() {
    expect_id("param_range", str(H, "param Re: Real = 2000 in 1..1000\n"));
    expect_id("param_range", str(H, "param Re: Real = 100 in 1..1000\nstudy s: grid {\n    Re: [10, 2000]\n}\n"));
    expect_id("study_param", str(H, "let D: Length = 1 m\nstudy s: grid {\n    D: [0.1, 0.2] m\n}\n"));
    expect_clean("a study", str(H, "param D: Length = 1 m in 0.1 m..2 m\nstudy s: grid {\n    D: [0.1, 0.2] m\n}\n"));
    expect_id("study_objective", str(H, "param D: Length = 1 m\nstudy s: minimize(drag) {\n    D: 0.1 m..0.2 m\n}\n"));
    expect_id("scales_missing", "model m: 3d, SI\nlet a: Length = 1 m\n");
    expect_id("scales_dependent", "model m: 3d, SI\nlet a: Length = 1 m\nscales a, 2*a\n");
    expect_id("scales_incomplete", "model m: 3d, SI\nlet a: Length = 1 m\nlet r: Density = 1 kg/m^3\nscales a\n");
    expect_id("process_setting", str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nmanufacture b: cnc_mill {\n    draft: 1 deg\n}\n"));
}

// The JSON report holds no raw control byte, also for a tab in the model.
fn test_output() {
    check::check_text("t.skn", str(H, SOLID, "\nbc f: sigma *\t2 = 0\n"));
    var out: u8[>..] = [];
    check::write_json(out);
    for c in out[..out.len - 1] {
        if c < 0x20 {
            failures++;
            print("FAIL a tab in JSON: got [", out, "]");
            return;
        }
    }
    if find(out[..], "\\t") < 0 {
        failures++;
        print("FAIL a tab in JSON: got [", out, "]");
    }
}

// The default methods of `src` (spec §13), one "physics,mesh,order,stabilize" per line.
fn expect_chosen(what: u8[:], src: u8[:], want: u8[:]) {
    check::check_text("t.skn", src);
    var got: u8[>..] = [];
    for c in check::chosen {
        format(got, c.physics, ",", c.mesh, ",", c.order, ",", c.stabilize, "\n");
    }
    if got == want && check::diags.len == 0 { return; }
    failures++;
    print("FAIL ", what, ": got [", got, "] and ", check::diags.len, " diagnostics");
}

fn test_method_defaults() {
    let elastic = str(H, SOLID, "\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n");
    expect_chosen("a solid", elastic, "s,m1,2,none\n");
    expect_chosen("the mesh on the region of the physics", str(elastic, "mesh m2 on f { size: 0.2 m }\n"), "s,m1,2,none\n");
    expect_chosen("a method item", str(elastic, "method s: cg { order: 1 }\n"), "");
    expect_chosen("Taylor-Hood and SUPG", str(H2, FLOW), "flow,m1,{ v: 2, p: 1 },supg\n");
    let heat = str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nmaterial cu {\n    k: Conductivity = 400 W/(m*K)\n    rho: Density = 8900 kg/m^3\n    c_p: SpecificHeat = 385 J/(kg*K)\n}\nmesh m1 on b { size: 0.1 m }\n");
    expect_chosen("conduction", str(heat, "physics h on b: heat solves T: Temperature {\n    material cu\n}\n"), "h,m1,2,none\n");
    expect_chosen("convection of a heat", str(heat, "let w: Vector<Velocity> = (1, 0, 0) m/s\nphysics h on b solves T: Temperature {\n    material cu\n    eq rho*c_p*dot(w, grad(T)) - div(k*grad(T)) = 0 W/m^3\n}\n"),
                  "h,m1,2,supg\n");
    expect_chosen("convection in a let", str(heat, "let w: Vector<Velocity> = (1, 0, 0) m/s\nphysics h on b solves T: Temperature {\n    material cu\n    let adv: Power/Volume = rho*c_p*dot(w, grad(T))\n    let net: Power/Volume = adv\n    eq net - div(k*grad(T)) = 0 W/m^3\n}\n"),
                  "h,m1,2,supg\n");
    // The report: TOON and JSON.
    check::check_text("t.skn", str(H2, FLOW));
    var toon: u8[>..] = [];
    check::write_toon(toon);
    var json: u8[>..] = [];
    check::write_json(json);
    if find(toon[..], "method[1]{physics,kind,mesh,order,stabilize}:\n  flow,cg,m1,\"{ v: 2, p: 1 }\",supg\n") < 0 {
        failures++;
        print("FAIL method TOON: got [", toon, "]");
    }
    if find(json[..], "\"method\":[{\"physics\":\"flow\",\"kind\":\"cg\",\"mesh\":\"m1\",\"order\":\"{ v: 2, p: 1 }\",\"stabilize\":\"supg\"}]") < 0 {
        failures++;
        print("FAIL method JSON: got [", json, "]");
    }
}

// The ids of the diagnostics of `skein lint` on `src`, one per line.
fn lint_ids(src: u8[:]) -> u8[>..] {
    check::lint_text("t.skn", src);
    var out: u8[>..] = [];
    for d in check::diags {
        out.append(d.id);
        out.push('\n');
    }
    return out;
}

// `skein lint` on `src` gives exactly the ids `want`, one per line.
fn expect_lints(what: u8[:], src: u8[:], want: u8[:]) {
    let got = lint_ids(src);
    if got == want { return; }
    failures++;
    print("FAIL lint ", what, ": got [", got, "]");
}

fn test_lint() {
    let elastic = str(H, SOLID, "\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n");
    let fixed = str(elastic, "step go: static {\n    bc f: u = 0\n    output u\n}\n");
    expect_lints("a clean model", fixed, "");
    expect_lints("unused_let", str(fixed, "let a: Length = 1 m\n"), "unused_let\n");
    expect_lints("unused_let in a physics", str(H, SOLID, "\nphysics s on b solves T: Temperature {\n    material steel\n    let q: Vector<HeatFlux> = fourier(1 W/(m*K), T)\n    eq div(-fourier(1 W/(m*K), T)) = 0\n}\nstep go: static {\n    bc f: T = 300 K\n    output T\n}\n"), "unused_let\n");
    expect_lints("allow", str(fixed, "let a: Length = 1 m\nallow unused_let: \"a is for the reader\"\n"), "");
    expect_lints("unused_material", str(H, SOLID), "unused_let\nunused_material\n");
    expect_lints("unused_field", str(elastic, "step go: static {\n    bc f: u = 0\n}\n"), "unused_field\n");
    expect_lints("a field in a report", str(elastic, "step go: static {\n    bc f: u = 0\n    report max(von_mises(sigma), b)\n}\n"), "");
    // A load on a region before its birth, and after its death.
    let staged = str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nlet top: Solid = box((0, 0, 1) m, (1, 1, 2) m)\nlet all: Solid = b | top\nlet f: Surface = boundary(b, z == 0 m)\nlet lid: Surface = boundary(top, z == 2 m)\nmaterial steel {\n    E: Pressure = 200 GPa\n    nu: Real = 0.3\n}\nmesh m1 on all { size: 0.1 m }\nphysics s on all: linear_elastic solves u: Vector<Length> {\n    material steel\n}\nbc f: u = 0\n");
    expect_lints("load_no_cells before the birth", str(staged, "step one: static {\n    bc lid: sigma*n = (0, 0, -1) MPa\n    output u\n}\nstep two: static {\n    birth top\n}\n"), "load_no_cells\n");
    expect_lints("a load after the birth", str(staged, "step one: static {\n    output u\n}\nstep two: static {\n    birth top\n    bc lid: sigma*n = (0, 0, -1) MPa\n}\n"), "");
    expect_lints("load_no_cells after the death", str(staged, "bc lid: sigma*n = (0, 0, -1) MPa\nstep one: static {\n    output u\n}\nstep two: static {\n    death top\n}\n"), "load_no_cells\n");
    let soft = str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nlet f: Surface = boundary(b, x == 0 m)\nmaterial rubber {\n    E: Pressure = 10 MPa\n    nu: Real = 0.499\n}\nmesh m1 on b { size: 0.1 m }\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material rubber\n}\nbc f: u = 0\nstep go: static {\n    output u\n}\n");
    expect_lints("volumetric_lock", str(soft, "method s: cg { mesh: m1, order: 1 }\n"), "volumetric_lock\n");
    expect_lints("order 2 does not lock", str(soft, "method s: cg { mesh: m1, order: 2 }\n"), "");
    let impact = str(elastic, "bc f: u = 0\nstep go: transient(end: 1 s, dt: ");
    expect_lints("dt_unstable", str(impact, "1 ms, scheme: ssprk3) {\n    output u\n}\n"), "dt_unstable\n");
    expect_lints("a stable dt", str(impact, "1 us, scheme: ssprk3) {\n    output u\n}\n"), "");
    expect_lints("an implicit scheme", str(impact, "1 ms, scheme: bdf2) {\n    output u\n}\n"), "");
    expect_lints("dt_unstable at a crack tip", str(H2, "let c: Surface = rect((0, 0) m, (1, 1) m)\nlet f: Curve = boundary(c, x == 0 m)\nlet tip: Point = point((0.5, 0) m)\nmaterial steel {\n    E: Pressure = 200 GPa\n    nu: Real = 0.3\n    rho: Density = 7800 kg/m^3\n}\nmesh m1 on c {\n    size: 0.1 m\n    crack_tip tip: 1 mm\n}\nphysics s on c: linear_elastic solves u: Vector<Length> {\n    material steel\n    assume plane_strain(1 m)\n}\nbc f: u = 0\nstep go: transient(end: 1 s, dt: 1 us, scheme: ssprk3) {\n    output u\n}\n"), "dt_unstable\n");
    let flow = str("model m: 2d, SI\nscales 1 m, 1 m/s, 1 kg/m^3\n", FLOW, "\nmethod flow: cg {\n    mesh: m1\n    order: { v: 2, p: 1 }\n    stabilize: ");
    let steady = "\n}\nstep go: static {\n    bc lid: v = (1, 0) m/s\n    output v\n}\n";
    expect_lints("stabilization_missing", str(flow, "none", steady), "stabilization_missing\n");
    expect_lints("supg", str(flow, "supg", steady), "");
    // An error stops the lint rules.
    expect_lints("an error first", str(H, "let a: Length = 1 s\n"), "unit_mismatch\n");
}

// The lint output: the class `lint`, and the count in the summary (command_line.md §6.1).
fn test_lint_output() {
    let src = str(H, "let a: Length = 1 m\n");
    check::lint_text("t.skn", src);
    var toon: u8[>..] = [];
    check::write_toon(toon);
    var json: u8[>..] = [];
    check::write_json(json);
    if find(toon[..], "  lint,unused_let,t.skn,3,5,") < 0 || find(toon[..], "summary: errors=0 lints=1\n") < 0 {
        failures++;
        print("FAIL lint TOON: got [", toon, "]");
    }
    if find(json[..], "{\"class\":\"lint\",\"id\":\"unused_let\"") < 0 || find(json[..], "\"summary\":{\"errors\":0,\"lints\":1}") < 0 {
        failures++;
        print("FAIL lint JSON: got [", json, "]");
    }
}

// Each model passes `check` and `lint`, except a femahoy model with a `// blocked-by:` mark.
fn test_corpus(dir: u8[:]) {
    var names: u8[>..] = [];
    guard list_dir(dir, names) else { return; }
    for name in split(names, '\n') {
        guard name.len > 4 && name[name.len - 4..] == ".skn" else { continue; }
        let path = str(dir, "/", name);
        var src: u8[>..] = [];
        guard read_file(path, src) else {
            failures++;
            print("FAIL ", path, ": cannot read the file");
            continue;
        }
        if find(src[..], "// blocked-by:") >= 0 { continue; }
        check::lint_text(path, src);
        for d in check::diags {
            failures++;
            print("FAIL ", path, ":", d.line, ":", d.col, ": ", d.id, ": ", d.msg);
        }
    }
}

fn test_shells() {
    expect_clean("a shell with offset beams", str(H, SHELL));
    let mat = "material m2 {\n    E1: Pressure = 148 GPa\n    E2: Pressure = 10 GPa\n    G12: Pressure = 4 GPa\n    nu12: Real = 0.3\n    t_ply: Length = 1 mm\n    h: Length = 4 mm\n    E: Pressure = 70 GPa\n    nu: Real = 0.3\n    A: Area = 1e-4 m^2\n    I_y: Length^4 = 1e-8 m^4\n    I_z: Length^4 = 2e-8 m^4\n    J: Length^4 = 1e-9 m^4\n}\n";
    expect_id("geometry_dim", str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\n", mat,
        "physics p on b: reissner_mindlin(top: (0, 0, 1)) solves w: Vector<Length>, theta: Vector<Angle> {\n    material m2\n}\n"));
    expect_id("geometry_dim", str(H, "let b: Solid = box((0, 0, 0) m, (1, 1, 1) m)\nlet s: Surface = boundary(b, x == 0 m)\n", mat,
        "physics p on s: timoshenko(A: A, I_y: I_y, I_z: I_z, J: J, y_axis: (0, 1, 0)) solves w: Vector<Length>, theta: Vector<Angle> {\n    material m2\n}\n"));
    expect_id("type_mismatch", str(H, "let D3: Real = laminate(ply_stiffness(1 Pa, 1 Pa, 1 Pa, 0.3), 0 deg, false, 1 mm, 1 mm, (1, 0, 0), (0, 0, 1), 3)\n", SHELL));
    expect_id("rank_mismatch", str(H, "let g: Real = tangent(1, (0, 0, 1))\n", SHELL));
    expect_id("rank_mismatch", str(H, SHELL, "\nstep more: static {\n    report argmax(w, wall)\n}\n"));
}

fn main() {
    test_names();
    test_use();
    test_types();
    test_physics();
    test_conditions();
    test_crack_tip();
    test_method_defaults();
    test_steps();
    test_shells();
    test_studies();
    test_output();
    test_lint();
    test_lint_output();
    test_corpus("samples");
    test_corpus("examples");
    test_corpus("examples/femahoy");
    if failures > 0 {
        print(failures, " failed");
        exit(1);
    }
    print("check: ok");
}
