// The state that a page keeps from one solve of a model to the next, so
// that a solve after a slider move starts from the last solution and the
// last AMG levels. The page sets SKEIN_WARM to a directory of its file
// system, and copies the files of that directory from each run to the
// next (src/codegen/page.goose). A program without SKEIN_WARM keeps
// nothing: each native run solves from 0.
//
// Each solve of a run has a key: its number in the run. The files of solve
// k are <dir>/<k>.x (the solution) and <dir>/<k>.amg.* (the levels).
//
// The levels of the last run are a preconditioner of the new matrix when
// the matrix has the same dofs and pattern and the same fixed dofs: the
// fine level takes the new values, and the coarse levels stay. CG then
// converges to the same tolerance; the levels change only its iterations.
import std;
import os;
namespace solvers;
import .amg;

var warm_solves = 0;

// The directory of SKEIN_WARM, or "" if the program keeps nothing.
fn warm_dir() -> u8[>..] {
    var dir: u8[>..] = [];
    if !env("SKEIN_WARM", dir) { dir.clear(); }
    return dir;
}

// The key of the next solve of the run.
fn warm_key() -> i64 {
    warm_solves++;
    return warm_solves;
}

fn warm_path(key: i64, part: u8[:]) -> u8[>..] { return str(warm_dir(), "/", key, ".", part); }

fn load_f64s(path: u8[:]) -> f64[>..], bool {
    var bytes: u8[>..] = [];
    guard read_file(path, bytes) else { return [], false; }
    return from_bytes<f64[>..]>(bytes);
}

fn load_i64s(path: u8[:]) -> i64[>..], bool {
    var bytes: u8[>..] = [];
    guard read_file(path, bytes) else { return [], false; }
    return from_bytes<i64[>..]>(bytes);
}

fn load_i32s(path: u8[:]) -> i32[>..], bool {
    var bytes: u8[>..] = [];
    guard read_file(path, bytes) else { return [], false; }
    return from_bytes<i32[>..]>(bytes);
}

// The free dofs of x from the solution of solve `key` of the last run, if
// it has the size of x. The fixed dofs keep their values.
fn warm_start(key: i64, x: f64[:], fixed: u8[:]) -> bool {
    guard x.len == fixed.len else { abort("warm_start: x and fixed of one size"); }
    guard warm_dir().len > 0 else { return false; }
    let last, ok = load_f64s(warm_path(key, "x"));
    guard ok && last.len == x.len else { return false; }
    for d in x.len { if fixed[d] == 0 { x[d] = last[d]; } }
    return true;
}

// Keeps x as the solution of solve `key`.
fn warm_keep(key: i64, x: f64[:]) {
    guard warm_dir().len > 0 else { return; }
    var copy: f64[>..] = [];
    copy.append(x);
    write_file(warm_path(key, "x"), to_bytes(copy));
}

// The signature of the dofs and the pattern of a matrix and of its fixed
// dofs: the levels of one matrix fit another of the same signature.
fn pattern_signature(row_start: i64[:], col: i32[:], fixed: u8[:]) -> i64[4] {
    var s: i64[4] = [row_start.len - 1, col.len, 0, 0];
    for r in row_start.len - 1 {
        let w = r % 1009 + 1;
        for e in row_start[r]..row_start[r + 1] { s[2] = (s[2] + w * (col[e] + 1)) % 1000000007; }
    }
    for d in fixed.len { if fixed[d] != 0 { s[3] = (s[3] + d + 1) % 1000000007; } }
    return s;
}

// The AMG levels of solve `key` of the last run, with the matrix
// (row_start, col, val) as the fine level, if they fit it. False if the
// page keeps nothing or the levels do not fit.
fn amg_reuse(key: i64, row_start: i64[:], col: i32[:], val: f64[:], fixed: u8[:]) -> bool {
    guard warm_dir().len > 0 else { return false; }
    let sig, ok = load_i64s(warm_path(key, "amg.sig"));
    let now = pattern_signature(row_start, col, fixed);
    guard ok && sig.len == 4 && sig[0] == now[0] && sig[1] == now[1] && sig[2] == now[2] && sig[3] == now[3] else { return false; }
    let n, ok_n = load_i64s(warm_path(key, "amg.n"));
    let rs, ok_rs = load_i64s(warm_path(key, "amg.rs"));
    let cl, ok_col = load_i32s(warm_path(key, "amg.col"));
    let vl, ok_val = load_f64s(warm_path(key, "amg.val"));
    let prs, ok_prs = load_i64s(warm_path(key, "amg.prs"));
    let pcol, ok_pcol = load_i32s(warm_path(key, "amg.pcol"));
    let pval, ok_pval = load_f64s(warm_path(key, "amg.pval"));
    let chol, ok_chol = load_f64s(warm_path(key, "amg.chol"));
    guard ok_n && ok_rs && ok_col && ok_val && ok_prs && ok_pcol && ok_pval && ok_chol && n.len >= 1 &&
          n[0] == row_start.len - 1 else { return false; }
    return amg_restore(row_start, col, val, n, rs, cl, vl, prs, pcol, pval, chol);
}

// Keeps the AMG levels of the last amg_setup as the levels of solve `key`.
fn amg_keep(key: i64, row_start: i64[:], col: i32[:], fixed: u8[:]) {
    guard warm_dir().len > 0 else { return; }
    var sig: i64[>..] = [];
    sig.append(pattern_signature(row_start, col, fixed));
    var n: i64[>..] = [];
    var rs: i64[>..] = [];
    var cl: i32[>..] = [];
    var vl: f64[>..] = [];
    var prs: i64[>..] = [];
    var pcol: i32[>..] = [];
    var pval: f64[>..] = [];
    var chol: f64[>..] = [];
    amg_coarse(n, rs, cl, vl, prs, pcol, pval, chol);
    write_file(warm_path(key, "amg.n"), to_bytes(n));
    write_file(warm_path(key, "amg.rs"), to_bytes(rs));
    write_file(warm_path(key, "amg.col"), to_bytes(cl));
    write_file(warm_path(key, "amg.val"), to_bytes(vl));
    write_file(warm_path(key, "amg.prs"), to_bytes(prs));
    write_file(warm_path(key, "amg.pcol"), to_bytes(pcol));
    write_file(warm_path(key, "amg.pval"), to_bytes(pval));
    write_file(warm_path(key, "amg.chol"), to_bytes(chol));
    // The signature last: a run stopped before it leaves no levels to use.
    write_file(warm_path(key, "amg.sig"), to_bytes(sig));
}
