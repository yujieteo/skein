// Tests of the finite volume method for systems of conservation laws.
// `./build.sh test` runs this program from the root of the repository, and
// a non-zero exit status is a failure.
//
// The mesh is the unit square in n × n squares of 2 triangles, written as
// skein run writes a mesh. The tests check that the method keeps a uniform
// state exactly, that the spectral radius of the Jacobian of a flux is the
// speed of its fastest wave, and that a bump of the advection equation
// moves with an error that falls when the mesh is finer.
import std;
import math;
import os;
import cons from "../cons";
import solvers from "../solvers";   // cons imports it

var failures = 0;

fn fail(what: u8[:]) {
    failures++;
    print("FAIL ", what);
}

// Writes the mesh of n × n squares, with the mask 1 on each side of the
// boundary, in the directory dir.
fn write_square(dir: u8[:], n: i64) -> bool {
    var xy: f64[>..] = [];
    var tri: i64[>..] = [];
    for j in n + 1 {
        for i in n + 1 {
            xy.push(i / n as f64);
            xy.push(j / n as f64);
        }
    }
    for j in n {
        for i in n {
            let a = j * (n + 1) + i;
            let b = a + 1;
            let c = a + n + 2;
            let d = a + n + 1;
            tri.push(a); tri.push(b); tri.push(c);
            tri.push(a); tri.push(c); tri.push(d);
        }
    }
    var mask: i64[>..] = [];
    for t in tri { mask.push(1); }
    var flags: u8[>..] = [];
    for t in tri.len / 3 { flags.push(1); }
    // make_dir makes one level, and a fresh clone has no .skein yet.
    for c, i in dir {
        if c == '/' && i > 0 { make_dir(dir[..i]); }
    }
    return make_dir(dir) && write_file(str(dir, "/cons.xy"), to_bytes(xy[..])) && write_file(str(dir, "/cons.tri"), to_bytes(tri[..])) &&
           write_file(str(dir, "/cons.sides"), to_bytes(mask[..])) && write_file(str(dir, "/cons.cells"), to_bytes(flags[..]));
}

// The advection of one component with the velocity (1, 0.5).
fn advect_flux(u: f64[:], nx: f64, ny: f64, f: f64[:]) { f[0] = (nx + 0.5 * ny) * u[0]; }
fn no_source(u: f64[:], x: f64, y: f64, t: f64, s: f64[:]) { s[0] = 0.0; }
var inflow_bump = false;
// The state outside: the bump at its place at time t, or 1.
fn exact(x: f64, y: f64, t: f64) -> f64 {
    if !inflow_bump { return 1.0; }
    let dx = x - t - 0.3;
    let dy = y - 0.5 * t - 0.3;
    return 1.0 + exp(-80.0 * (dx * dx + dy * dy));
}
fn advect_ghost(mask: i64, u: f64[:], x: f64, y: f64, t: f64, nx: f64, ny: f64, g: f64[:]) { g[0] = exact(x, y, t); }
fn advect_wall(mask: i64, u: f64[:], x: f64, y: f64, t: f64, nx: f64, ny: f64, f: f64[:]) {}
fn advect_start(x: f64, y: f64, k: i64) -> f64 { return exact(x, y, 0.0); }
fn one(u: f64[:], x: f64, y: f64) -> f64 { return u[0]; }

// The wave system (p, m) with ddt(p) + div(c² m) = 0, ddt(m) + grad(p) = 0
// in x: the waves go at ±c.
let SPEED = 3.0;
fn wave_flux(u: f64[:], nx: f64, ny: f64, f: f64[:]) {
    f[0] = SPEED * SPEED * u[1] * nx;
    f[1] = u[0] * nx;
}

// The L1 error of the bump at time 0.2 on n × n squares.
fn bump_error(dir: u8[:], n: i64, limiter: i64) -> f64 {
    guard write_square(dir, n) && cons::read_mesh(dir, 1) else {
        fail("the mesh of the square cannot be written or read");
        return 1.0;
    };
    inflow_bump = true;
    cons::setup(1, limiter);
    cons::start(advect_start);
    if !cons::solve_transient(0.2, 0.2 / (4 * n), advect_flux, no_source, advect_ghost, advect_wall) { fail("the bump does not finish"); }
    var e = 0.0;
    for c in cons::cell_count() {
        e += abs(cons::sol[c] - exact(cons::centre_x[c], cons::centre_y[c], 0.2)) * cons::area[c];
    }
    return e;
}

fn main() {
    let dir = ".skein/results/test_cons";
    // A uniform state stays.
    guard write_square(dir, 8) && cons::read_mesh(dir, 1) else {
        fail("the mesh of the square cannot be written or read");
        exit(1);
    };
    inflow_bump = false;
    cons::setup(1, cons::VAN_LEER);
    cons::start(advect_start);
    if !cons::solve_transient(0.5, 0.02, advect_flux, no_source, advect_ghost, advect_wall) { fail("the uniform state does not finish"); }
    var worst = 0.0;
    for c in cons::cell_count() { worst = max(worst, abs(cons::sol[c] - 1.0)); }
    if worst > 1.0e-12 { fail(str("a uniform state changes by ", worst)); }
    // Each cell has 3 neighbours or boundary sides.
    var sides = 0;
    for c in cons::cell_count() { for s in 3 { if cons::neighbour[3 * c + s] >= 0 { sides++; } } }
    if sides != 2 * (3 * 8 * 8 - 2 * 8) { fail(str("the mesh has ", sides, " sides with a neighbour")); }
    // The spectral radius of the waves.
    cons::setup(2, cons::VAN_LEER);
    cons::start(advect_start);
    cons::jacobians(cons::sol, wave_flux);
    let q = cons::spectral(0, 1.0, 0.0);
    if q < SPEED || q > 1.05 * SPEED { fail(str("the spectral radius of the waves is ", q, ", not ", SPEED)); }
    let q45 = cons::spectral(0, 0.6, 0.8);
    if abs(q45 - 0.6 * SPEED) > 0.05 * SPEED { fail(str("the radius on the normal (0.6, 0.8) is ", q45, ", not ", 0.6 * SPEED)); }
    // The bump moves, and the error falls with the mesh.
    let e1 = bump_error(dir, 16, cons::VAN_LEER);
    let e2 = bump_error(dir, 32, cons::VAN_LEER);
    if !(e2 < 0.75 * e1) { fail(str("the error of the bump falls from ", e1, " to ", e2, " only")); }
    if e2 > 0.01 { fail(str("the error of the bump is ", e2)); }
    if failures > 0 { exit(1); }
    print("cons: ok");
}
