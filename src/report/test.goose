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

fn main() {
    test_json();
    test_small_vtu();
    test_case_count();
    if failures > 0 {
        print(failures, " report tests fail");
        exit(1);
    }
    print("report: ok");
}
