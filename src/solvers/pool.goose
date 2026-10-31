// The pool of the threads in a case (spec §20.5): the assembly on blocks
// of cells, and the product of the matrix and a vector on blocks of rows.
// Goose threads share no memory and have no atomics or locks, so the
// workers get their data as copies, and they give their results in
// messages to the thread of the case. Each cell matrix and each row is the
// same sum in the same order on each number of threads, and the thread of
// the case adds the cell matrices in cell order. Thus the results do not
// change with the thread count or with the order of the messages.
import std;
import os;
namespace solvers;
import .sparse;

// The threads of the case, or 0 before the first call of case_threads.
var threads = 0;

// SKEIN_THREADS, which the study runner sets for each case, or all
// hardware threads if it is not set (docs/command_line.md, `--threads`).
fn case_threads() -> i64 {
    if threads == 0 {
        threads = max(1, hardware_threads());
        var s: u8[>..] = [];
        if env("SKEIN_THREADS", s) {
            let n, good = parse_int(s);
            if good && n >= 1 { threads = n; }
        }
    }
    return threads;
}

// Sets the threads of the case to n, in place of SKEIN_THREADS.
fn use_threads(n: i64) {
    guard n >= 1 else { abort("use_threads: n >= 1"); }
    threads = n;
}

// --- the assembly on blocks of cells --------------------------------------------

// The cells of a task of the assembly. Reason: one message for each cell
// costs more than the cell matrix of a tet4.
let CELL_BLOCK = 16;

// A task of the assembly: the cells lo <= c < hi. A task with lo < 0 stops
// the worker.
struct CellJob { lo: i64, hi: i64 }

// The data of the cells of a task, one cell after the other, and the first
// cell that has no element routine, or -1.
struct CellData { lo: i64, bad: i64, data: f64[] }

// Gives the tasks of `cells` cells to the `workers` workers that the
// caller started, and calls add(lo, hi, data) for each task in cell order.
// A worker takes a CellJob, and puts the CellData of its cells, at most
// `size` values for each cell. At most 2 workers tasks are out at a time,
// and task b waits in slot b % (2 workers) until the tasks before it are
// added. Then the stop tasks end the workers. It gives the first cell
// that has no element routine, or -1, and then it adds no task.
fn gather_cells<ADD>(cells: i64, workers: i64, size: i64) -> i64 {
    guard cells >= 1 && workers >= 1 && size >= 1 else { abort("gather_cells: cells, workers and size >= 1"); }
    let tasks = (cells + CELL_BLOCK - 1) / CELL_BLOCK;
    let slots = 2 * workers;
    let stride = CELL_BLOCK * size;
    var store: f64[>..] = [];
    zeros(store, slots * stride);
    var count: i64[>..] = [];   // the values in each slot, or -1 if it is empty
    for s in slots { count.push(-1); }
    var next = 0;
    while next < min(slots, tasks) { put_cells(next, cells); next++; }
    var bad = -1;
    var added = 0;
    while added < tasks {
        let d = qget<CellData>();
        let s = (d.lo / CELL_BLOCK) % slots;
        guard d.data.len <= stride else { abort("gather_cells: at most size values for each cell"); }
        for v, i in d.data { store[s * stride + i] = v; }
        count[s] = d.data.len;
        if d.bad >= 0 && (bad < 0 || d.bad < bad) { bad = d.bad; }
        while added < tasks && count[added % slots] >= 0 {
            let at = (added % slots) * stride;
            let lo = added * CELL_BLOCK;
            if bad < 0 { ADD(lo, min(cells, lo + CELL_BLOCK), store[at..at + count[added % slots]]); }
            count[added % slots] = -1;
            added++;
            if next < tasks { put_cells(next, cells); next++; }
        }
    }
    for w in workers { qput(CellJob { lo: -1, hi: -1 }); }
    return bad;
}

fn put_cells(task: i64, cells: i64) {
    qput(CellJob { lo: task * CELL_BLOCK, hi: min(cells, (task + 1) * CELL_BLOCK) });
}

// --- the product on blocks of rows -------------------------------------------------

// The rows lo <= r < lo + row_start.len - 1 of a matrix, as the argument of
// the spawn of a worker. row_start[0] is the start of row lo in the matrix.
struct Rows { lo: i64, row_start: i64[], col: i32[], val: f64[] }

// The vector x of a product. Two products in a row use the two types in
// turn, and a worker takes one of each in turn. The next product starts
// only after each worker gave its rows, so each worker takes exactly one
// x of each product.
struct XEven { stop: bool, x: f64[] }
struct XOdd { stop: bool, x: f64[] }

// The rows lo <= r < lo + y.len of a product.
struct RowsY { lo: i64, y: f64[] }

var rows_ids: i64[>..] = [];
var rows_products = 0;

thread_fn rows_worker(b: Rows) {
    var y: f64[>..] = [];
    zeros(y, b.row_start.len - 1);
    loop {
        let e = qget<XEven>();
        guard !e.stop else { break; };
        multiply_rows(b.row_start, b.col, b.val, e.x, y);
        qput(RowsY { lo: b.lo, y: y[..] });
        let o = qget<XOdd>();
        guard !o.stop else { break; };
        multiply_rows(b.row_start, b.col, b.val, o.x, y);
        qput(RowsY { lo: b.lo, y: y[..] });
    }
}

// Starts case_threads() workers, each with a block of the rows of the
// matrix. Together they hold one copy of the matrix. With one thread, it
// starts none.
fn rows_start(row_start: i64[:], col: i32[:], val: f64[:]) {
    guard rows_ids.len == 0 && row_start.len >= 1 && row_start[row_start.len - 1] == col.len && col.len == val.len else { abort("rows_start: no workers yet, and a matrix"); }
    let rows = row_start.len - 1;
    let workers = min(case_threads(), rows);
    if workers < 2 { return; }
    for w in workers {
        let lo = rows * w / workers;
        let hi = rows * (w + 1) / workers;
        let a = row_start[lo];
        let b = row_start[hi];
        rows_ids.push(thread_spawn(rows_worker, Rows { lo: lo, row_start: row_start[lo..hi + 1], col: col[a..b], val: val[a..b] }));
    }
    rows_products = 0;
}

// y = A x, on the workers of rows_start, or on this thread if it started
// none. A is the matrix of rows_start.
fn rows_multiply(row_start: i64[:], col: i32[:], val: f64[:], x: f64[:], y: f64[:]) {
    if rows_ids.len == 0 {
        multiply(row_start, col, val, x, y);
        return;
    }
    guard x.len == y.len && row_start.len == y.len + 1 else { abort("rows_multiply: x and y have one entry for each row"); }
    for w in rows_ids.len { put_x(false, x); }
    for w in rows_ids.len {
        let p = qget<RowsY>();
        for v, i in p.y { y[p.lo + i] = v; }
    }
    rows_products++;
}

// Stops the workers of rows_start.
fn rows_stop() {
    for w in rows_ids.len { put_x(true, []); }
    for id in rows_ids { thread_wait(id); }
    rows_ids.clear();
}

fn put_x(stop: bool, x: f64[:]) {
    if rows_products % 2 == 0 { qput(XEven { stop: stop, x: x }); } else { qput(XOdd { stop: stop, x: x }); }
}
