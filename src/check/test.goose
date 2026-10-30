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
    expect_id("use_item", str(H, "use \"lib.skn\"\n"));
    expect_id("unknown_unit", str(H, "let a: Length = 5 mms\n"));
    expect_id("table_row_width", str(H, "table tb{a: Length in m, b: Real}:\n    1, 2, 3\n"));
}

fn test_types() {
    expect_id("type_mismatch", str(H, "let count: Int = 2.5\n"));
    expect_id("unit_mismatch", str(H, "let a: Pressure = 1 GPa + 1 kg/m^3\n"));
    expect_id("rank_mismatch", str(H, "let a: Length = (1, 2, 3) m + 1 m\n"));
    expect_id("ambiguous_product", str(H, "let a: Vector<Length> = (1, 2, 3) m\nlet c: Real = a*a\n"));
    expect_id("unit_fractional", str(H, "let a: Pressure = 1 Pa\nlet c: Real = sqrt(a)\n"));
    expect_id("unit_transcendental", str(H, "let a: Temperature = 300 K\nlet c: Real = exp(a)\n"));
    expect_id("unit_affine", str(H, "let k: ThermalConductivity = 1.2 W/(m*degC)\n"));
    expect_id("not_constant", str(H, "param a: Length = x\n"));
    expect_id("geometry_dim", str(H2, "let s: Solid = rect((0, 0) m, (1, 1) m)\n"));
    expect_id("geometry_space", str(H, "let s: Surface = circle((0, 0) m, 1 m)\n"));
    expect_id("nondim_unit", "model m: 2d, nondim\nlet a: Real = 1 m\n");
    // Gaps G7, G12, G29, G30: these are valid.
    expect_clean("G7 re, imag, arg, conj", str(H, "let w: Complex<Real> = 1 + 2*im\nlet a: Real = re(w) + imag(w) + abs(w)\nlet b: Angle = arg(conj(w))\n"));
    expect_clean("G12 atan, atan2", str(H, "let a: Angle = atan(0.5) + atan2(1 m, 2 m)\n"));
    expect_clean("G29 a negative degC", str(H, "let a: Temperature = -5 degC\n"));
    expect_clean("a vector with a unit, a number and a name", str(H, "let h: Length = 1 m\nlet p: Vector<Length> = (0.5, h, 0) m\n"));
    expect_clean("G30 a real exponent", str(H, "param q: Real = 3.0\nlet a: Real = 0.5^(2/q)\n"));
}

fn test_physics() {
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

fn test_steps() {
    let elastic = str(H, SOLID, "\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\n");
    expect_id("procedure_args", str(elastic, "step go: transient(end: 1 s) {\n    bc f: u = 0\n}\n"));
    expect_id("procedure_args", str(elastic, "step go: static(rate: 2) {\n    bc f: u = 0\n}\n"));
    expect_id("not_constant", str(elastic, "step go: transient(end: t, dt: 1 s) {\n    bc f: u = 0\n}\n"));
    expect_clean("G9 a list of frequencies", str(elastic, "step go: harmonic(freq: [10, 20] Hz) {\n    bc f: u = 0\n}\n"));
    expect_clean("G6 an eigenvalue", str(elastic, "step go: modal(count: 4) {\n    bc f: u = 0\n    report f1 = eigenvalue(1)\n}\nexpect f1: eigenvalue(1) == 10 Hz within 1%\n"));
    expect_clean("an eigenvalue with no units", "model m: 3d, nondim\nlet b: Solid = box((0, 0, 0), (1, 1, 1))\nlet f: Surface = boundary(b, x == 0)\nmaterial steel {\n    E: Pressure = 1000\n    nu: Real = 0.3\n    rho: Density = 1\n}\nmesh m1 on b { size: 0.1 }\nphysics s on b: linear_elastic solves u: Vector<Length> {\n    material steel\n}\nstep go: modal(count: 4) {\n    bc f: u = 0\n    report f1 = eigenvalue(1)\n}\nexpect f1: eigenvalue(1) == 2.5 within 1%\n");
    expect_id("type_mismatch", str(elastic, "step go: static {\n    bc f: u = 0\n    report eigenvalue(1)\n}\n"));
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

// Each model passes `check`, except a femahoy model with a `// blocked-by:` mark.
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
        check::check_text(path, src);
        for d in check::diags {
            failures++;
            print("FAIL ", path, ":", d.line, ":", d.col, ": ", d.id, ": ", d.msg);
        }
    }
}

fn main() {
    test_names();
    test_types();
    test_physics();
    test_conditions();
    test_steps();
    test_studies();
    test_output();
    test_corpus("samples");
    test_corpus("examples/femahoy");
    if failures > 0 {
        print(failures, " failed");
        exit(1);
    }
    print("check: ok");
}
