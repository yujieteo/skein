// `use` lines (spec §4.1). The items of each library file go into the
// program of the model, after the `use` line, so each later pass reads one
// list of items in one namespace.
import std;
import os;
namespace syntax;
import .tree;
import .lexer;
import .parser;

// Parses the model `text` of file `path` with the items of its `use` files.
// It gives the N_PROGRAM node, or -1 after an error, which is in `errors`.
fn parse_model(path: u8[:], text: u8[:]) -> i64 {
    reset();
    let f = add_file(path, text);
    lex(f);
    let prog = parse(f);
    read_uses(prog);
    if errors.len > 0 { return -1; }
    return prog;
}

// Reads the file of each `use` line of program `prog`, and puts its items
// after that line.
fn read_uses(prog: i64) {
    var item = nodes[prog].first;
    while item >= 0 {
        let next = nodes[item].next;
        if nodes[item].kind == N_USE { read_use(prog, item); }
        item = next;
    }
}

fn read_use(prog: i64, u: i64) {
    let tok = nodes[u].tok;
    let name = text(tok);
    guard name.len >= 2 else { return; };
    // The path is relative to the directory of the model file.
    let model = files[tokens[tok].file].path;
    var cut = -1;
    for c, i in model {
        if c == '/' { cut = i; }
    }
    let path = str(model[..cut + 1], name[1..name.len - 1]);
    var src: u8[>..] = [];
    guard read_file(path, src) else {
        error_at(tok, "use_item", str("the library file '", path, "' does not exist"),
                 "give the path of the file, relative to this model");
        return;
    };
    let f = add_file(path, src);
    lex(f);
    let lib = parse_library(f);
    let first = nodes[lib].first;
    guard first >= 0 else { return; };
    let last = nodes[lib].last;
    nodes[last].next = nodes[u].next;
    nodes[u].next = first;
    if nodes[prog].last == u { nodes[prog].last = last; }
}
