// Tests of `skein report`. `./build.sh test` runs this program from the root
// of the repository, and a non-zero exit status is a failure.
import std;
import os;
import syntax from "../syntax";
import check from "../check";
import .report;

var failures = 0;

fn expect_text(what: u8[:], got: u8[:], want: u8[:]) {
    if got == want { return; }
    failures++;
    print("FAIL ", what, ": got [", got, "]");
}

fn expect_has(what: u8[:], got: u8[:], want: u8[:]) {
    if find(got, want) >= 0 { return; }
    failures++;
    print("FAIL ", what, ": no [", want, "] in [", got, "]");
}

// A source text cannot end the script element of the page.
fn test_json() {
    var out: u8[>..] = [];
    report::json_string(out, "a \"b\"\n</script>\t");
    expect_text("json string", out, "\"a \\\"b\\\"\\n<\\/script>\\u0009\"");
}

// One tet10 becomes one tet4 with its 4 corner nodes and their point data.
// The 6 mid-side nodes and their values go. The largest value 3 becomes
// 32767 with the scale 3/32767.
fn test_small_vtu() {
    var text: u8[>..] = [];
    text.append("""
        <VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="10" NumberOfCells="1">
        <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1  .5 0 0  .5 .5 0  0 .5 0  0 0 .5  .5 0 .5  0 .5 .5
        </DataArray></Points>
        <Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 4 5 6 7 8 9</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">10</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">24</DataArray></Cells>
        <PointData><DataArray type="Float64" Name="T" NumberOfComponents="1" format="ascii">0 1 2 3 4 5 6 7 8 9</DataArray></PointData>
        </Piece></UnstructuredGrid></VTKFile>
        """);
    var mesh: u8[>..] = [];
    var fields: u8[>..] = [];
    expect_text("small vtu", if report::small_vtu(text, mesh, fields) { "yes" } else { "no" }, "yes");
    expect_has("small vtu size", mesh, "NumberOfPoints=\"4\" NumberOfCells=\"1\"");
    expect_has("small vtu field", fields, str("Name=\"T\" NumberOfComponents=\"1\" Scale=\"", 3.0 / 32767.0, "\""));
    let n, scale = report::shorts([0.0, -1.5, 3.0]);
    expect_text("shorts", str(n[0], " ", n[1], " ", n[2]), "0 -16384 32767");
    expect_text("no cells", if report::small_vtu("<VTKFile/>", mesh, fields) { "yes" } else { "no" }, "no");
}

fn test_case_count() {
    expect_text("case count", str(report::case_count("study:\ncases[16]{case,L,P}:\n")), "16");
    expect_text("no study", str(report::case_count("solve[1]{step}:\n")), "0");
}

// With `--no-run`, the report uses the table of a study run only if the run
// wrote the hash of the current source.
fn test_last_run() {
    let dir = ".skein/studies/report_test";
    report::make_dirs(dir);
    write_file(str(dir, "/s.toon"), "study:\ncases[1]{case,status}:\n  1,ok\n");
    write_file(str(dir, "/s.hash"), "");
    let toon = ".skein/studies/report_test/run.toon";
    syntax::parse_model("m.skn", "model m\n");
    expect_text("no hash", str(report::last_run("report_test", "s", toon)), "2");
    write_file(str(dir, "/s.hash"), syntax::source_hash());
    expect_text("same source", str(report::last_run("report_test", "s", toon)), "0");
    expect_text("table copy", report::text_of(toon), "study:\ncases[1]{case,status}:\n  1,ok\n");
    syntax::parse_model("m.skn", "model m\nparam L = 2 m\n");
    expect_text("edited source", str(report::last_run("report_test", "s", toon)), "2");
}

fn main() {
    test_json();
    test_small_vtu();
    test_case_count();
    test_last_run();
    if failures > 0 {
        print(failures, " report tests fail");
        exit(1);
    }
    print("report: ok");
}
