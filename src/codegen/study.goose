// The cases and the table of `skein run --study` (spec §20, docs/command_line.md
// §6.5). run.goose runs the cases. Each case is the model text with the case
// values in place of the defaults of the varied `param` items.
import std;
import os;
import math;
namespace codegen;
import syntax from "../syntax";
import check from "../check";
import .codegen;

// --- the cases (spec §20.2) ----------------------------------------------------

// A varied parameter: its name and the offsets of its default value in the
// model text.
struct Varied {
    name: u8[..64],
    lo: i64,
    hi: i64,
}

var varied: Varied[>..];
// The value text of parameter `p` in case `c` is `case_values[c * varied.len + p]`.
var case_values: u8[..256][>..];
var case_count = 0;
// The `expect` lines of the study (spec §15.3). Case c appends
// `case_expects[expect_at[c]..expect_at[c + 1]]` to its model text.
var case_expects: u8[>..];
var expect_at: i64[>..];

fn case_value(c: i64, p: i64) -> const u8[:] { return case_values[c * varied.len + p]; }

// The `study` item named `name` of program `prog`, or -1.
fn find_study(prog: i64, name: u8[:]) -> i64 {
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        if kind(c) == syntax::N_STUDY && codegen::name(c) == name { return c; }
        c = syntax::nodes[c].next;
    }
    return -1;
}

// The `param` item named `word` of program `prog`, or -1.
fn find_param(prog: i64, word: u8[:]) -> i64 {
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        if kind(c) == syntax::N_PARAM && name(c) == word { return c; }
        c = syntax::nodes[c].next;
    }
    return -1;
}

// The named argument `word` of procedure `pr`, or -1.
fn proc_arg(pr: i64, word: u8[:]) -> i64 {
    var a = syntax::nodes[pr].first;
    while a >= 0 {
        if syntax::nodes[a].aux >= 0 && syntax::text(syntax::nodes[a].aux) == word { return syntax::nodes[a].first; }
        a = syntax::nodes[a].next;
    }
    return -1;
}

// The texts of the values of study entry `v`: each item of a list, with the
// unit of the list after a bare number.
fn list_texts(v: i64, out: u8[..256][>..]&) {
    if kind(v) != syntax::N_LIST {
        var t: u8[..256] = [];
        set256(t, source_of(v));
        out.push(t);
        return;
    }
    let unit = syntax::child_of_kind(v, syntax::N_UNIT);
    var c = syntax::nodes[v].first;
    while c >= 0 {
        if kind(c) != syntax::N_UNIT {
            var t: u8[..256] = [];
            if unit >= 0 && check::is_bare(c) { set256(t, str(source_of(c), " ", source_of(unit))); } else { set256(t, source_of(c)); }
            out.push(t);
        }
        c = syntax::nodes[c].next;
    }
}

// Reads the cases of study `s` of program `prog` into `varied` and
// `case_values`. False on an error, and `syntax::errors` holds it.
fn study_cases(prog: i64, s: i64) -> bool {
    varied.clear();
    case_values.clear();
    case_count = 0;
    let pr = syntax::nodes[s].first;
    let procedure = name(pr);
    let body = syntax::nodes[pr].next;
    var entries: i64[>..] = [];
    var c = syntax::nodes[body].first;
    while c >= 0 {
        if kind(c) == syntax::N_SETTING { entries.push(c); }
        c = syntax::nodes[c].next;
    }
    guard entries.len > 0 else {
        unsupported(tok(s), "the study varies no param", "add an entry, for example Re: [10, 20]");
        return false;
    };
    for e in entries {
        for w in varied {
            guard w.name != name(e) else {
                syntax::error_at(tok(e), "study_param", str("the study varies '", name(e), "' twice"), "give each param one entry");
                return false;
            };
        }
        let p = find_param(prog, name(e));
        guard p >= 0 else {
            unsupported(tok(e), str("the param '", name(e), "' is not in the model file"), "declare the param in the model file, not in a use file");
            return false;
        };
        let lo, hi = source_span(syntax::nodes[syntax::nodes[p].first].next);
        guard hi <= syntax::files[0].len else {
            unsupported(tok(e), str("the param '", name(e), "' is not in the model file"), "declare the param in the model file, not in a use file");
            return false;
        };
        var w = Varied { name: [], lo: lo, hi: hi };
        set64(w.name, name(e));
        varied.push(w);
    }
    var good = false;
    if procedure == "grid" || procedure == "zip" {
        good = list_cases(procedure == "grid", entries);
    } else if procedure == "sample" {
        good = sample_cases(pr, entries);
    } else {
        unsupported(tok(pr), str("the study procedure ", procedure, " is not in this build"),
                    "use grid, zip or sample");
    }
    if good { study_expects(body); }
    return good;
}

// The text in case `c` of the expected value or the `within` `v` of an
// expect of a study: value c of a list, else `v`.
fn case_text(v: i64, c: i64) -> u8[>..] {
    guard kind(v) == syntax::N_LIST else { return str(source_of(v)); };
    var texts: u8[..256][>..] = [];
    list_texts(v, texts);
    return str(texts[c]);
}

// The `expect` lines of each case of the study with the block `body`. The
// checker gives each list of an expect one value for each case.
fn study_expects(body: i64) {
    case_expects.clear();
    expect_at.clear();
    for c in 0..case_count {
        expect_at.push(case_expects.len);
        var e = syntax::nodes[body].first;
        while e >= 0 {
            if kind(e) == syntax::N_EXPECT {
                format(case_expects, "\nexpect ", name(e), ": ", source_of(syntax::child(e, 0)), " == ",
                       case_text(syntax::child(e, 1), c), " within ", case_text(syntax::child(e, 2), c), "\n");
            }
            e = syntax::nodes[e].next;
        }
    }
    expect_at.push(case_expects.len);
}

// The cases of `grid` (the Cartesian product of the lists, with the last
// entry the fastest) or of `zip` (value i of each list).
fn list_cases(grid: bool, entries: i64[:]) -> bool {
    // The values of entry e are texts[first[e]..first[e] + count[e]].
    var texts: u8[..256][>..] = [];
    var first: i64[>..] = [];
    var count: i64[>..] = [];
    for e in entries {
        let v = syntax::nodes[e].first;
        guard kind(v) != syntax::N_RANGE else {
            unsupported(tok(v), "a grid or a zip takes a list of values, not a range", "write a list, for example [10, 20, 40]");
            return false;
        };
        first.push(texts.len);
        list_texts(v, texts);
        count.push(texts.len - first[first.len - 1]);
        guard count[count.len - 1] > 0 else {
            syntax::error_at(tok(e), "procedure_args", str("the list of '", name(e), "' has no values"), "give the list one value or more");
            return false;
        };
    }
    if grid {
        case_count = 1;
        for n in count { case_count *= n; }
    } else {
        case_count = count[0];
        for n, i in count {
            guard n == case_count else {
                syntax::error_at(tok(entries[i]), "procedure_args",
                                 str("the lists of a zip have one length: '", name(entries[i]), "' has ", n, " values, not ", case_count),
                                 str("give each list ", case_count, " values"));
                return false;
            };
        }
    }
    for c in 0..case_count {
        for p in 0..entries.len { case_values.push(texts[0]); }
        var rest = c;
        var p = entries.len;
        while p > 0 {
            p--;
            let i = if grid { rest % count[p] } else { c };
            if grid { rest /= count[p]; }
            case_values[c * entries.len + p] = texts[first[p] + i];
        }
    }
    return true;
}

// The Sobol direction numbers of dimensions 2 to 8 (Joe and Kuo, 2008):
// the degree s, the polynomial a and the first s numbers m.
let SOBOL_S = [1, 2, 3, 3, 4, 4, 5];
let SOBOL_A = [0, 1, 1, 2, 1, 4, 2];
let SOBOL_M = [[1, 0, 0, 0, 0], [1, 3, 0, 0, 0], [1, 3, 1, 0, 0], [1, 1, 1, 0, 0], [1, 1, 3, 3, 0],
               [1, 3, 5, 13, 0], [1, 1, 5, 5, 17]];

// The 32 direction numbers of Sobol dimension `d`, from 0.
fn sobol_directions(d: i64) -> i64[32] {
    var v = default<i64[32]>();
    for k in 0..32 { v[k] = 1 << (31 - k); }
    if d == 0 { return v; }
    let s = SOBOL_S[d - 1];
    let a = SOBOL_A[d - 1];
    for k in 0..s { v[k] = SOBOL_M[d - 1][k] << (31 - k); }
    for k in s..32 {
        var x = v[k - s] ^ (v[k - s] >> s);
        for j in 1..s {
            if ((a >> (s - 1 - j)) & 1) == 1 { x ^= v[k - j]; }
        }
        v[k] = x;
    }
    return v;
}

// The points in [0, 1) of `sample(n, method, seed)`: `n` rows of `dims`
// values. `lhs` takes one point in each of `n` strata of each dimension, in
// a random order. `sobol` takes points 0 to n - 1 of the Sobol sequence,
// with a random digital shift of each dimension.
fn sample_points(n: i64, dims: i64, sobol: bool, seed: i64) -> f64[>..] {
    var r = rng { s: seed as! u64 };
    var t: f64[>..] = [];
    for i in 0..n * dims { t.push(0.0); }
    for d in 0..dims {
        if sobol {
            let v = sobol_directions(d);
            let shift = (rand_u64(r) >> 32) as! i64;
            var x = 0;
            for i in 0..n {
                t[i * dims + d] = (x ^ shift) / 4294967296.0;
                // The next point, in Gray code order, differs in the direction of the lowest zero bit of i.
                var c = 0;
                while ((i >> c) & 1) == 1 { c++; }
                x ^= v[c];
            }
        } else {
            var strata: i64[>..] = [];
            for i in 0..n { strata.push(i); }
            shuffle(strata, r);
            for i in 0..n { t[i * dims + d] = (strata[i] + rand_flt(r)) / n; }
        }
    }
    return t;
}

// The value and the unit text of the literal `n`: a number, or `-` and a
// number, with an optional unit. False for another expression.
recursive fn literal(n: i64) -> f64, const u8[:], bool {
    if kind(n) == syntax::N_UNARY && syntax::text(tok(n)) == "-" {
        let v, unit, ok = literal(syntax::nodes[n].first);
        return -v, unit, ok;
    }
    guard kind(n) == syntax::N_NUMBER else { return 0.0, "", false; };
    let v = if syntax::nodes[n].aux >= 0 { -number(tok(n)) } else { number(tok(n)) };
    let u = syntax::child_of_kind(n, syntax::N_UNIT);
    return v, if u >= 0 { source_of(u) } else { "" }, true;
}

// The cases of `sample(n: Int, method: lhs | sobol, seed: Int)`: each value
// is `lo + t*(hi - lo)` for a point t in [0, 1), in the unit of the range.
fn sample_cases(pr: i64, entries: i64[:]) -> bool {
    let fix = "write sample(n: 8, method: lhs, seed: 1), with method lhs or sobol";
    let n_arg = proc_arg(pr, "n");
    let method = proc_arg(pr, "method");
    let seed_arg = proc_arg(pr, "seed");
    guard n_arg >= 0 && method >= 0 && seed_arg >= 0 && kind(n_arg) == syntax::N_NUMBER && kind(seed_arg) == syntax::N_NUMBER &&
          kind(method) == syntax::N_NAME && (name(method) == "lhs" || name(method) == "sobol") else {
        syntax::error_at(tok(pr), "procedure_args", "sample needs n, method and seed", fix);
        return false;
    };
    let n = number(tok(n_arg)) as! i64;
    guard n >= 1 && n == number(tok(n_arg)) else {
        syntax::error_at(tok(n_arg), "procedure_args", "n of sample is a whole number of 1 or more", fix);
        return false;
    };
    let sobol = name(method) == "sobol";
    guard !sobol || entries.len <= SOBOL_S.len + 1 else {
        syntax::error_at(tok(method), "procedure_args", str("sobol takes at most ", SOBOL_S.len + 1, " params"), "use method: lhs");
        return false;
    };
    // The range of entry d is lows[d] to highs[d], in the unit units[d].
    var lows: f64[>..] = [];
    var highs: f64[>..] = [];
    var units: u8[..64][>..] = [];
    for e in entries {
        let range = syntax::nodes[e].first;
        guard kind(range) == syntax::N_RANGE else {
            unsupported(tok(e), "a sample takes a range for each param", "write a range, for example Re: 10..40");
            return false;
        };
        let lo, lo_unit, lo_ok = literal(syntax::nodes[range].first);
        let hi, hi_unit, hi_ok = literal(syntax::nodes[syntax::nodes[range].first].next);
        guard lo_ok && hi_ok && lo_unit == hi_unit else {
            unsupported(tok(e), "a sample takes a range of two numbers in one unit", "write a range, for example D: 0.1 m..0.2 m");
            return false;
        };
        lows.push(lo);
        highs.push(hi);
        var w: u8[..64] = [];
        set64(w, lo_unit);
        units.push(w);
    }
    let t = sample_points(n, entries.len, sobol, number(tok(seed_arg)) as! i64);
    case_count = n;
    for c in 0..n {
        for d in 0..entries.len {
            var v: u8[..256] = [];
            set256(v, str(lows[d] + t[c * entries.len + d] * (highs[d] - lows[d]), if units[d].len > 0 { " " } else { "" }, units[d]));
            case_values.push(v);
        }
    }
    return true;
}

// The model `text` of case `c`: the case values in place of the defaults,
// and the `expect` lines of the study at the end.
fn case_source(text: u8[:], c: i64) -> u8[>..] {
    var order: i64[>..] = [];
    for p in 0..varied.len { order.push(p); }
    sort(order) { a, b => varied[a].lo < varied[b].lo }
    var out: u8[>..] = [];
    var at = 0;
    for p in order {
        out.append(text[at..varied[p].lo]);
        format(out, "(", case_value(c, p), ")");
        at = varied[p].hi;
    }
    out.append(text[at..]);
    out.append(case_expects[expect_at[c]..expect_at[c + 1]]);
    return out;
}

// The model `text` of `path` with the value of each `--set name=value` of
// `sets` in place of the default value of its `param` (spec §20.1). The
// checker then checks the value as source text, with its type and its
// range. False if a word is not name=value or names no param of the model
// file: the message is on stderr.
fn set_source(path: u8[:], text: u8[:], sets: const u8[:][:], out: u8[>..]&) -> bool {
    out.clear();
    let prog = syntax::parse_model(path, text);
    guard prog >= 0 else {
        out.append(text);
        return true;
    };
    var los: i64[>..] = [];
    var his: i64[>..] = [];
    var values: i64[>..] = [];
    for s, i in sets {
        var eq = -1;
        for c, j in s { if c == '=' && eq < 0 { eq = j; } }
        guard eq > 0 && eq + 1 < s.len else {
            write_stderr(str("skein run: --set takes name=value, not '", s, "'\n"));
            return false;
        };
        let p = find_param(prog, s[..eq]);
        var lo = -1;
        var hi = -1;
        if p >= 0 { lo, hi = source_span(syntax::nodes[syntax::nodes[p].first].next); }
        guard p >= 0 && hi <= syntax::files[0].len else {
            write_stderr(str("skein run: the model file has no param '", s[..eq], "' for --set\n"));
            return false;
        };
        for l in los {
            guard l != lo else {
                write_stderr(str("skein run: --set gives '", s[..eq], "' two values\n"));
                return false;
            };
        }
        los.push(lo);
        his.push(hi);
        values.push(i);
    }
    var order: i64[>..] = [];
    for k in los.len { order.push(k); }
    sort(order) { a, b => los[a] < los[b] }
    var at = 0;
    for k in order {
        let s = sets[values[k]];
        var eq = 0;
        while s[eq] != '=' { eq++; }
        out.append(text[at..los[k]]);
        format(out, "(", s[eq + 1..], ")");
        at = his[k];
    }
    out.append(text[at..]);
    return true;
}

// --- the memory ledger (spec §20.5) -------------------------------------------

// The ledger of the study runner. The runner reads the budget once, at the
// start of the study. A case starts only if its estimate fits in the
// budget less the estimates of the cases that run, or if no case runs:
// then the case checks itself before its assembly, and stops with a
// message if it does not fit (solvers::memory_check). Reason: a case
// reaches its peak late, at the setup of its solver, so the free memory at
// its start does not show what the running cases take later.
struct Ledger {
    budget: f64,    // the bytes for all cases, or 0 if they are not known
    per_dof: f64,   // the bytes of a case for each dof of its mesh, or 0 if they are not known
    source: i64,    // the case of per_dof, or -1
    measured: bool, // per_dof is a measured peak, not an estimate
    used: f64,      // the sum of the estimates of the cases that run
    most: i64,      // the most cases that ran at a time
}

fn new_ledger(budget: f64) -> Ledger {
    return Ledger { budget: budget, per_dof: 0.0, source: -1, measured: false, used: 0.0, most: 0 };
}

// The estimate in bytes of a case with `dofs` dofs: the whole budget while
// the ledger knows no case.
fn case_need(l: Ledger, dofs: i64) -> f64 {
    return if l.per_dof > 0.0 { l.per_dof * dofs } else { l.budget };
}

// Whether a case of `need` bytes may start while `running` cases run.
fn admits(l: Ledger, need: f64, running: i64) -> bool {
    return running == 0 || l.budget <= 0.0 || l.used + need <= l.budget;
}

// Learns from case c of `dofs` dofs, which ended with the measured peak
// `peak`. The first measured peak takes the place of the estimate, and a
// larger peak for each dof takes the place of a smaller one.
fn learn(l: Ledger&, c: i64, dofs: i64, peak: f64) {
    guard peak > 0.0 && dofs > 0 else { return; };
    let per = peak / dofs;
    if l.measured && per <= l.per_dof { return; }
    l.per_dof = per;
    l.source = c;
    l.measured = true;
}

// The value of the line "key: value" of the TOON `text`, or -1.
fn toon_number(text: u8[:], key: u8[:]) -> f64 {
    for line in split(text, '\n') {
        if !starts_with(line, key) { continue; }
        let v, good = parse_flt(trim(line[key.len..]));
        return if good { v } else { -1.0 };
    }
    return -1.0;
}

fn set1024(into: u8[..1024]&, s: u8[:]) {
    into = [];
    into.append(s[..min(s.len, 1024)]);
}

// The end of a case: its result record, with its standard output in
// `outs[out_lo..out_hi]`.
struct CaseEnd {
    ended: bool,
    built: bool,
    status: i64,
    seconds: f64,
    err: u8[..1024],
    out_lo: i64,
    out_hi: i64,
}

var outs: u8[>..];

// The fields of the TOON row `line`, after its indent, with the quotes
// and escapes taken off.
fn toon_fields(line: u8[:]) -> u8[..512][>..] {
    var fields: u8[..512][>..] = [];
    var f: u8[..512] = [];
    var quoted = false;
    var i = 0;
    while i < line.len && line[i] == ' ' { i++; }
    while i < line.len {
        let c = line[i];
        if quoted && c == '\\' && i + 1 < line.len {
            i++;
            let e = line[i];
            if f.len < 512 { f.push(if e == 'n' { '\n' } else if e == 't' { '\t' } else { e }); }
        } else if c == '"' {
            quoted = !quoted;
        } else if c == ',' && !quoted {
            fields.push(f);
            f = [];
        } else if f.len < 512 {
            f.push(c);
        }
        i++;
    }
    fields.push(f);
    return fields;
}

// The columns of the table after `verify`: the named reports, then the
// expects of the items and the steps, then the expects of the study.
// `cells[c * columns.len + j]` is the value of column j in case c.
var columns: u8[..64][>..];
var report_columns = 0;
var cells: u8[..256][>..];

fn add_column(t: i64) {
    var w: u8[..64] = [];
    set64(w, syntax::text(t));
    columns.push(w);
}

// The named reports and the expects of program `prog` and of its study `s`.
fn result_columns(prog: i64, s: i64) {
    columns.clear();
    cells.clear();
    var c = syntax::nodes[prog].first;
    while c >= 0 {
        if kind(c) == syntax::N_STEP {
            var st = syntax::nodes[syntax::nodes[c].first].next;
            while st >= 0 {
                if kind(st) == syntax::N_REPORT {
                    var r = syntax::nodes[st].first;
                    while r >= 0 {
                        let aux = syntax::nodes[r].aux;
                        if aux >= 0 && position(columns, syntax::text(aux)) < 0 { add_column(aux); }
                        r = syntax::nodes[r].next;
                    }
                }
                st = syntax::nodes[st].next;
            }
        }
        c = syntax::nodes[c].next;
    }
    report_columns = columns.len;
    c = syntax::nodes[prog].first;
    while c >= 0 {
        if kind(c) == syntax::N_EXPECT { add_column(tok(c)); }
        if kind(c) == syntax::N_STEP {
            var st = syntax::nodes[c].first;
            while st >= 0 {
                if kind(st) == syntax::N_EXPECT { add_column(tok(st)); }
                st = syntax::nodes[st].next;
            }
        }
        c = syntax::nodes[c].next;
    }
    c = syntax::nodes[syntax::nodes[syntax::nodes[s].first].next].first;
    while c >= 0 {
        if kind(c) == syntax::N_EXPECT { add_column(tok(c)); }
        c = syntax::nodes[c].next;
    }
}

// Reads the output of case `e`: it adds the value of each named report and
// the result of each expect to `cells`. It gives the message of a failed
// case, empty for a case that ran, and true when no verification record of
// the case fails.
fn read_case(e: CaseEnd) -> u8[>..], bool {
    var pass = true;
    var section: u8[..16] = [];
    var printed = false;
    let row = cells.len;
    for j in 0..columns.len { cells.push([]); }
    for line in split(outs[e.out_lo..e.out_hi], "\n") {
        if line.len > 0 && line[0] != ' ' {
            section = [];
            if starts_with(line, "report[") { section.append("report"); }
            if starts_with(line, "verify[") { section.append("verify"); printed = true; }
            continue;
        }
        let f = toon_fields(line);
        for col, j in columns {
            if section == "report" && f.len == 3 && j < report_columns && starts_with(f[1], str(col, " =")) {
                set256(cells[row + j], f[2]);
            }
            if section == "verify" && f.len == 5 && j >= report_columns && f[1] == "expect" && starts_with(f[3], str(col, ":")) {
                set256(cells[row + j], f[2]);
            }
        }
        if section == "verify" && f.len == 5 && f[2] == "fail" { pass = false; }
    }
    // A program that printed its records stops with a non-zero status when
    // a solve does not converge or a result file cannot be written.
    if e.ended && e.built && e.status != 0 { pass = false; }
    var msg: u8[>..] = [];
    if !e.ended {
        msg.append(e.err);
    } else if !e.built {
        format(msg, "goose cannot build the program of the case: ", e.err);
    } else if !printed {
        format(msg, "the program stopped with the exit status ", e.status, ": ", e.err);
    }
    while msg.len > 0 && (msg[msg.len - 1] == '\n' || msg[msg.len - 1] == ' ') { msg.pop(); }
    return msg, pass;
}

// The table of study `study` (spec §20.6, command_line.md §6.5): one row
// for each case, then the diagnostic `case_failed` for each case that
// stopped. `at` is the token of the study, and `prog` the model. It puts
// the TOON in `toon`, and what to print in `out`. The result is the exit
// status (command_line.md §10).
fn study_table(study: u8[:], at: i64, jobs: i64, prog: i64, ends: CaseEnd[:], json: bool, out: u8[>..]&, toon: u8[>..]&) -> i64 {
    result_columns(prog, find_study(prog, study));
    var verify: bool[>..] = [];
    var msgs: u8[..1024][>..] = [];
    var failed: i64[>..] = [];
    var exit = 0;
    for e, c in ends {
        let msg, pass = read_case(e);
        var m: u8[..1024] = [];
        set1024(m, msg);
        msgs.push(m);
        verify.push(msg.len == 0 && pass);
        if msg.len > 0 { failed.push(c); }
        if !verify[c] { exit = 1; }
    }
    let t = syntax::tokens[at];
    let file = syntax::files[t.file].path;
    format(toon, "study: ", study, "\njobs: ", jobs, "\ncases[", ends.len, "]{case");
    for p in varied { format(toon, ",", p.name); }
    toon.append(",status,verify");
    for col in columns { format(toon, ",", col); }
    toon.append(",time_s}:\n");
    for e, c in ends {
        format(toon, "  ", c + 1);
        for p in 0..varied.len {
            toon.push(',');
            check::toon_value(toon, case_value(c, p));
        }
        format(toon, ",", if msgs[c].len == 0 { "ok" } else { "failed" }, ",", if verify[c] { "pass" } else { "fail" });
        for j in 0..columns.len {
            toon.push(',');
            check::toon_value(toon, cells[c * columns.len + j]);
        }
        format(toon, ",", round(e.seconds * 100.0) / 100.0, "\n");
    }
    if failed.len > 0 {
        format(toon, "diagnostics[", failed.len, "]{class,id,file,line,col,msg,fix}:\n");
        for c in failed {
            toon.append("  error,case_failed,");
            check::toon_value(toon, file);
            format(toon, ",", t.line, ",", t.col, ",");
            check::quoted(toon, str("case ", c + 1, " of the study ", study, " stopped"));
            toon.push(',');
            check::quoted(toon, msgs[c]);
            toon.push('\n');
        }
        format(toon, "summary: errors=", failed.len, " lints=0\n");
    }
    if !json {
        out.append(toon);
        return exit;
    }
    out.append("{\"study\":");
    check::quoted(out, study);
    format(out, ",\"jobs\":", jobs, ",\"cases\":[");
    for e, c in ends {
        if c > 0 { out.push(','); }
        format(out, "{\"case\":", c + 1);
        for p, i in varied {
            out.push(',');
            check::quoted(out, p.name);
            out.push(':');
            check::quoted(out, case_value(c, i));
        }
        format(out, ",\"status\":\"", if msgs[c].len == 0 { "ok" } else { "failed" }, "\",\"verify\":\"", if verify[c] { "pass" } else { "fail" }, "\"");
        for col, j in columns {
            out.push(',');
            check::quoted(out, col);
            out.push(':');
            check::quoted(out, cells[c * columns.len + j]);
        }
        format(out, ",\"time_s\":", round(e.seconds * 100.0) / 100.0, "}");
    }
    out.append("],\"diagnostics\":[");
    for c, i in failed {
        if i > 0 { out.push(','); }
        out.append("{\"class\":\"error\",\"id\":\"case_failed\",\"file\":");
        check::quoted(out, file);
        format(out, ",\"line\":", t.line, ",\"col\":", t.col, ",\"msg\":");
        check::quoted(out, str("case ", c + 1, " of the study ", study, " stopped"));
        out.append(",\"fix\":");
        check::quoted(out, msgs[c]);
        out.push('}');
    }
    out.append("]}\n");
    return exit;
}
