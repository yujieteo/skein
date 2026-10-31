// `results(file, step:, case:)` (spec §15.1): the checker reads the other
// model and gives each of its unknowns a type.
import std;
import os;
namespace check;
import syntax from "../syntax";
import .types;
import .state;
import .expr;
import .typeexpr;

// The directory of the file of token `tok`, with its final '/', or "".
fn dir_of(tok: i64) -> u8[>..] {
    let path = syntax::files[syntax::tokens[tok].file].path;
    var cut = -1;
    for c, i in path {
        if c == '/' { cut = i; }
    }
    var d: u8[>..] = [];
    if cut >= 0 { d.append(path[..cut + 1]); }
    return d;
}

// The string literal at token `tok`, without its quotes.
fn string_value(tok: i64) -> const u8[:] {
    let s = syntax::text(tok);
    if s.len < 2 { return s; }
    return s[1..s.len - 1];
}

// The type of `results(file, step: name, case: p)`. `file`, `step` and
// `case` are the argument nodes (N_ARG children were already unwrapped).
// `step` and `case` may be -1.
fn results_ty(tok: i64, file: i64, step: i64, case: i64) -> Ty {
    guard syntax::nodes[file].kind == syntax::N_STRING else {
        error(first_tok(file), "type_mismatch", "the file of results must be a string", "write the path in quotes");
        return any_ty();
    }
    let path = str(dir_of(tok), string_value(syntax::nodes[file].tok));
    var text: u8[>..] = [];
    guard read_file(path, text) else {
        error(first_tok(file), "results_missing", str("the file '", path, "' does not exist"),
              "give the path of a model, relative to this file");
        return any_ty();
    }
    let res = results_count;
    results_count++;
    let f = syntax::add_file(path, text);
    // The syntax errors of the other file stay in syntax::errors: the
    // checker copied the errors of the model before it came here.
    syntax::lex(f);
    let prog = syntax::parse(f);
    var steps: u8[>..] = [];
    var found = step < 0;
    let want = if step >= 0 && syntax::nodes[step].kind == syntax::N_STRING { string_value(syntax::nodes[step].tok) } else { "" };
    // `case: p` names a param of this model. The other model has a param
    // of the same name and type, and runs with the value of p (spec §15.1).
    let case_name = if case >= 0 && syntax::nodes[case].kind == syntax::N_NAME { tok_text(syntax::nodes[case].tok) } else { "" };
    var case_sym = -1;
    if case >= 0 {
        let s = visible.get(key(case_name));
        if s && syms[s].kind == S_PARAM { case_sym = s; }
        if case_sym < 0 {
            error(first_tok(case), "results_case", "the case of results is not a param of this model",
                  "write case: p, with p a param of this model and of the other model");
        }
    }
    var case_found = case_sym < 0;
    var item = syntax::nodes[prog].first;
    while item >= 0 {
        let k = syntax::nodes[item];
        if k.kind == syntax::N_PARAM && case_sym >= 0 && syntax::text(k.tok) == case_name {
            case_found = true;
            assign(first_tok(case), type_of(k.first, false), syms[case_sym].ty, str("the param '", case_name, "' of '", path, "'"));
        }
        if k.kind == syntax::N_STEP {
            if steps.len > 0 { steps.append(", "); }
            steps.append(syntax::text(k.tok));
            if syntax::text(k.tok) == want { found = true; }
        }
        if k.kind == syntax::N_PHYSICS {
            var c = k.first;
            while c >= 0 {
                // The unknowns and the `let` fields of the physics (spec §15.1).
                if syntax::nodes[c].kind == syntax::N_UNKNOWN || syntax::nodes[c].kind == syntax::N_LET {
                    var t = type_of(syntax::nodes[c].first, false);
                    t.varying = true;
                    result_fields.push(ResultField { res: res, name: key(syntax::text(syntax::nodes[c].tok)), ty: t });
                }
                c = syntax::nodes[c].next;
            }
        }
        item = k.next;
    }
    if !case_found {
        error(first_tok(case), "results_case", str("the model '", path, "' has no param '", case_name, "'"),
              "give the case a param that both models declare");
    }
    if !found {
        error(first_tok(step), "results_missing", str("the model '", path, "' has no step '", want, "'"),
              str("the steps in that file: ", steps));
    }
    var t = kind_ty(K_RESULTS);
    t.aux = res;
    return t;
}
