// Tests of `skein md`. `./build.sh test` runs this program from the root of
// the repository, and a non-zero exit status is a failure.
//
// Each model must give its golden document in tests/md/ exactly. README.md
// gives the command that writes the golden files again.
import std;
import os;
import syntax from "../syntax";
import check from "../check";
import .md;

var failures = 0;

// The models and their golden documents.
let GOLDEN = [
    "samples/01_cantilever.skn", "tests/md/01_cantilever.md",
    "examples/femahoy/heat_fin.skn", "tests/md/heat_fin.md",
];

// The first line where `got` and `want` differ, 1-based.
fn first_difference(got: u8[:], want: u8[:]) -> i64 {
    var line = 1;
    for i in min(got.len, want.len) {
        if got[i] != want[i] { return line; }
        if got[i] == '\n' { line++; }
    }
    return line;
}

// The document of `model`, or an empty one after a failure.
fn document(model: u8[:]) -> u8[>..] {
    var text: u8[>..] = [];
    var result: u8[>..] = [];
    if !read_file(model, text) {
        failures++;
        print("FAIL ", model, ": cannot read it");
        return result;
    }
    if !md::md_text(model, text, result) {
        failures++;
        print("FAIL ", model, ": the model has ", check::diags.len, " errors");
    }
    return result;
}

fn test_golden() {
    var i = 0;
    while i < GOLDEN.len {
        // Copies: the checker empties the syntax tables, and the compiler
        // cannot see that the names do not point into them.
        let model = str(GOLDEN[i]);
        let golden = str(GOLDEN[i + 1]);
        let got = document(model);
        var want: u8[>..] = [];
        if !read_file(golden, want) {
            failures++;
            print("FAIL ", golden, ": cannot read it");
        } else if got != want {
            failures++;
            print("FAIL ", model, ": the document differs from ", golden, " at line ",
                  first_difference(got, want));
        }
        i += 2;
    }
}

// The same model gives the same document again: no date and no state from
// the run before.
fn test_deterministic() {
    let a = document(str(GOLDEN[0]));
    let b = document(str(GOLDEN[2]));
    let c = document(str(GOLDEN[0]));
    if a != c || a == b {
        failures++;
        print("FAIL a second document of ", GOLDEN[0], " is not the same as the first");
    }
}

// Text outside the math has the Markdown escapes. The math has no backslash
// before a punctuation character.
fn test_escapes() {
    let text = """
        // Costs 5% & *more*: $x_1 [a] #2 <b> `c`
        model esc_1: 3d, SI

        fn twice_of(y: Real) -> Real = 2*y
        let a: Real = twice_of(3)
        let r: Real = 5%
        let k: Pressure/Length^0.5 = 2 MPa*mm^-0.5

        scales 1 m, 1 kg, 1 s, 1 K
        """;
    var result: u8[>..] = [];
    if !md::md_text("t.skn", text, result) {
        failures++;
        print("FAIL escapes: the model has ", check::diags.len, " errors");
        return;
    }
    let wants = [
        "Costs 5% & \\*more\\*: \\$x\\_1 \\[a\\] \\#2 \\<b\\> \\`c\\`",
        "`esc_1`",
        "$\\operatorname{twice}_{\\mathrm{of}}(3)$",
        "\\mathit{twice}_{\\mathrm{of}}(y) = 2 y",
        "| `r` | let | `Real` | $5$ % |",
        "`Pressure/Length^0.5`",
        "\\mathrm{mm}^{-0.5}",
    ];
    for w in wants {
        if find(result, w) < 0 {
            failures++;
            print("FAIL escapes: no '", w, "' in:\n", result);
        }
    }
}

fn main() {
    test_golden();
    test_deterministic();
    test_escapes();
    if failures > 0 {
        print("md: ", failures, " failures");
        exit(1);
    }
    print("md: ok");
}
