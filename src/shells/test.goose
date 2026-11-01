// Tests of the shells package: the gates of roadmap lane 5. `./build.sh
// test` runs this program from the root of the repository, and a non-zero
// exit status is a failure. Each gate prints its value, the reference and
// the error.
import std;
import math;
import os;
import shells from "../shells";
import solvers from "../solvers";   // the shells import it

var failures = 0;

fn expect_near(what: u8[:], got: f64, want: f64, rel: f64) {
    let err = abs(got - want) / abs(want);
    print("  ", what, ": ", got, " (reference ", want, ", error ", 100.0 * err, " %, tolerance ",
        100.0 * rel, " %)");
    if err > rel {
        failures++;
        print("FAIL ", what);
    }
}

fn expect_solve(what: u8[:]) {
    if !shells::solve() {
        failures++;
        print("FAIL ", what, ": the matrix is singular at dof ", shells::singular_dof);
    }
    print(what, ": ", shells::node_count(), " nodes, ", shells::shells.len, " shells, ",
        shells::beams.len, " beams, ", shells::dof_count(), " dofs, ", shells::last_seconds, " s");
}

// The sum of the reactions and the loads must be 0: the equilibrium check of
// docs/checks.md, 1e-6 of the largest load.
fn expect_equilibrium(what: u8[:], scale: f64) {
    let r = shells::reactions();
    var s: f64[3] = [0.0, 0.0, 0.0];
    for d in r.len { if d % 6 < 3 { s[d % 6] += r[d] + shells::load[d]; } }
    if abs(s[0]) + abs(s[1]) + abs(s[2]) > 1.0e-6 * scale {
        failures++;
        print("FAIL ", what, " equilibrium: reactions + loads = ", s[0], " ", s[1], " ", s[2]);
    }
}

fn fix_all(n: i64, dofs: i64[:]) { for d in dofs { shells::fix(n, d); } }

// n + 1 equal lines from a to b.
fn lines(a: f64, b: f64, n: i64) -> f64[>..] {
    var v: f64[>..] = [];
    for i in n + 1 { v.push(a + (b - a) * i / n); }
    return v;
}

// The section s for each of n cells.
fn same(s: i64, n: i64) -> i64[>..] {
    var v: i64[>..] = [];
    for i in n { v.push(s); }
    return v;
}

// --- the laminate --------------------------------------------------------------

// T800S/3900 (docs/studies/b787/problem.md §1.2), SI units.
let PLY_E1 = 148.0e9;
let PLY_E2 = 10.1e9;
let PLY_G12 = 3.94e9;
let PLY_NU12 = 0.32;
let PLY_T = 0.191e-3;

// [45/0/-45/90]ns from the outer face: the skin is n = 3 (24 plies).
fn quasi_isotropic(n: i64) -> i64 {
    shells::laminate_clear();
    let angles = [45.0, 0.0, -45.0, 90.0];
    for k in n {
        for a in angles {
            shells::laminate_ply(radians(a), PLY_T, PLY_E1, PLY_E2, PLY_G12, PLY_NU12, PLY_G12, PLY_G12);
        }
    }
    shells::laminate_mirror();
    return shells::laminate_section();
}

// A quasi-isotropic laminate has A = h (U1, U4, 0; U4, U1, 0; 0, 0, U5) and
// B = 0, with the invariants of Jones §2.9. The sheet (problem.md §3.3)
// gives its in-plane values E = 55.85 GPa and nu = 0.326.
fn test_laminate() {
    print("laminate [45/0/-45/90]3s of T800S/3900");
    shells::clear();
    let k = quasi_isotropic(3);
    let s = shells::sections[k];
    let nu21 = PLY_NU12 * PLY_E2 / PLY_E1;
    let q11 = PLY_E1 / (1.0 - PLY_NU12 * nu21);
    let q22 = PLY_E2 / (1.0 - PLY_NU12 * nu21);
    let q12 = PLY_NU12 * q22;
    let u1 = (3.0 * q11 + 3.0 * q22 + 2.0 * q12 + 4.0 * PLY_G12) / 8.0;
    let u4 = (q11 + q22 + 6.0 * q12 - 4.0 * PLY_G12) / 8.0;
    let h = 24.0 * PLY_T;
    expect_near("h", s.h, h, 1.0e-12);
    expect_near("A11 against h U1", s.a[0], h * u1, 1.0e-12);
    expect_near("A22 against h U1", s.a[4], h * u1, 1.0e-12);
    expect_near("A12 against h U4", s.a[1], h * u4, 1.0e-12);
    expect_near("A66 against h (U1 - U4)/2", s.a[8], h * (u1 - u4) / 2.0, 1.0e-12);
    var off = abs(s.a[2]) + abs(s.a[5]);
    for i in 9 { off += abs(s.b[i]); }
    if off > 1.0e-9 * s.a[0] {
        failures++;
        print("FAIL laminate: A16, A26 or B is not 0: ", off);
    }
    expect_near("E in-plane against the sheet", (s.a[0] * s.a[0] - s.a[1] * s.a[1]) / (s.a[0] * h), 55.85e9, 0.001);
    expect_near("nu in-plane against the sheet", s.a[1] / s.a[0], 0.326, 0.001);
}

// --- the Scordelis-Lo roof -------------------------------------------------------

// MacNeal and Harder (1985), "A proposed standard set of problems to test
// finite element accuracy", Finite Elements in Analysis and Design 1, 3-20:
// R = 25 ft, L = 50 ft, 40 deg on each side of the crown, t = 0.25 ft,
// E = 4.32e8 lb/ft², nu = 0, a weight of 90 lb/ft², diaphragms at the ends.
// The vertical displacement at the midside of the free edge is 0.3024 ft.
// The quarter x in [0, 25], θ in [0, 40 deg], with n x n cells of order 1
// (MITC4) or 2 (MITC9).
fn scordelis_lo(n: i64, order: i64) -> f64 {
    shells::clear();
    let s = shells::isotropic_section(4.32e8, 0.0, 0.25);
    let g = shells::mesh_cylinder_lines(25.0, lines(0.0, 25.0, n), lines(0.0, radians(40.0), n), order, same(s, n * n), 1.0, 0.0, 0.0);
    let m = order * n;
    for k in m + 1 {
        fix_all(shells::grid_node(g, 0, k), [shells::UX, shells::RY, shells::RZ]);
        fix_all(shells::grid_node(g, m, k), [shells::UY, shells::UZ]);
        fix_all(shells::grid_node(g, k, 0), [shells::UY, shells::RX, shells::RZ]);
    }
    for c in shells::shells.len { shells::shell_load(c, 0.0, 0.0, -90.0, false); }
    expect_solve(str("Scordelis-Lo roof ", n, " x ", n, " ", element(order)));
    expect_equilibrium("Scordelis-Lo", 90.0 * 25.0 * 25.0);
    return -shells::displacement(shells::grid_node(g, 0, m), shells::UZ);
}

fn test_scordelis_lo() {
    for n in [8, 16] { print("  ", n, " x ", n, ": ", scordelis_lo(n, 1), " ft"); }
    // MITC4 converges to the reference from below. The error of a 32 x 32
    // quarter is about 0.6 %, so 1 % finds a wrong term but not mesh error.
    expect_near("Scordelis-Lo uz at the free edge, 32 x 32 MITC4", scordelis_lo(32, 1), 0.3024, 0.01);
    // MITC9 gives 0.3001, 0.3009 and 0.3013 ft. MITC4 32 x 32 and MITC9
    // 16 x 16 agree to 0.2 %: each converges to about 0.3015 ft.
    var last = 1.0;
    var uz = 0.0;
    for n in [4, 8, 16] {
        uz = scordelis_lo(n, 2);
        let err = abs(uz - 0.3024) / 0.3024;
        print("  MITC9 ", n, " x ", n, ": ", uz, " ft, error ", 100.0 * err, " %");
        if err >= last {
            failures++;
            print("FAIL Scordelis-Lo MITC9: the error does not fall with the cell size");
        }
        last = err;
    }
    expect_near("Scordelis-Lo uz at the free edge, 16 x 16 MITC9", uz, 0.3024, 0.006);
}

// --- the pinched cylinder --------------------------------------------------------

// MacNeal and Harder (1985): R = 300 in, L = 600 in, t = 3 in, E = 3e6 psi,
// nu = 0.3, rigid diaphragms at the ends, two forces P = 1 lb at mid-length.
// The radial displacement under a force is 1.8248e-5 in. The octant
// x in [0, 300] (x = 0 is a diaphragm), θ in [0, 90 deg].
fn pinched_cylinder(n: i64, order: i64) -> f64 {
    shells::clear();
    let s = shells::isotropic_section(3.0e6, 0.3, 3.0);
    let g = shells::mesh_cylinder_lines(300.0, lines(0.0, 300.0, n), lines(0.0, radians(90.0), n), order, same(s, n * n), 1.0, 0.0, 0.0);
    let m = order * n;
    for k in m + 1 {
        fix_all(shells::grid_node(g, 0, k), [shells::UY, shells::UZ, shells::RX]);
        fix_all(shells::grid_node(g, m, k), [shells::UX, shells::RY, shells::RZ]);
        fix_all(shells::grid_node(g, k, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(shells::grid_node(g, k, m), [shells::UZ, shells::RX, shells::RY]);
    }
    shells::force(shells::grid_node(g, m, 0), shells::UZ, -0.25);
    expect_solve(str("pinched cylinder ", n, " x ", n, " ", element(order)));
    expect_equilibrium("pinched cylinder", 0.25);
    return -shells::displacement(shells::grid_node(g, m, 0), shells::UZ);
}

fn test_pinched_cylinder() {
    for n in [8, 16] { print("  ", n, " x ", n, ": ", pinched_cylinder(n, 1), " in"); }
    // The point load converges slowly. The error of a 32 x 32 octant is about
    // 1.1 %, so 2 % finds a wrong term but not mesh error.
    expect_near("pinched cylinder radial displacement, 32 x 32 MITC4", pinched_cylinder(32, 1), 1.8248e-5, 0.02);
    for n in [4, 8] { print("  MITC9 ", n, " x ", n, ": ", pinched_cylinder(n, 2), " in"); }
    expect_near("pinched cylinder radial displacement, 16 x 16 MITC9", pinched_cylinder(16, 2), 1.8248e-5, 0.02);
}

// --- the thin cylinder under pressure ----------------------------------------------

// A closed thin cylinder with the 787 skin (problem.md §1, §3): R = 2.883 m,
// p = 97.5 kPa, the [45/0/-45/90]3s laminate. Membrane theory gives the
// hoop stress pR/t, the axial stress pR/(2t), and the hoop strain
// a21 Nx + a22 Nθ with a = A⁻¹. The quarter θ in [0, 90 deg] of a bay
// x in [0, 1.2192 m], with symmetry at x = 0 and the closure load pR/2 per
// length at x = 1.2192 m.
fn test_thin_cylinder(order: i64) {
    let r = 2.883;
    let p = 97.5e3;
    let len = 1.2192;
    let nx = 8 / order;
    let nt = 24 / order;
    shells::clear();
    let s = quasi_isotropic(3);
    let h = shells::sections[s].h;
    let g = shells::mesh_cylinder_lines(r, lines(0.0, len, nx), lines(0.0, radians(90.0), nt), order, same(s, nx * nt), 1.0, 0.0, 0.0);
    let la = order * nx;
    let lb = order * nt;
    for j in lb + 1 { fix_all(shells::grid_node(g, 0, j), [shells::UX, shells::RY, shells::RZ]); }
    for i in la + 1 {
        fix_all(shells::grid_node(g, i, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(shells::grid_node(g, i, lb), [shells::UZ, shells::RX, shells::RY]);
    }
    for c in shells::shells.len { shells::shell_load(c, p, 0.0, 0.0, true); }
    for j in nt { closure_load(g, la, order, j, p * r / 2.0); }
    expect_solve(str("thin cylinder ", nx, " x ", nt, " ", element(order), ", 24-ply skin"));
    expect_equilibrium("thin cylinder", p * r * len);
    // The cell at x = 0.61 m, θ = 45 deg.
    let c = shells::grid_cell(g, nx / 2, nt / 2);
    let n, eps = shells::shell_resultants(c);
    let t = radians(45.0 + 45.0 / nt);
    let et: f64[3] = [0.0, cos(t), -sin(t)];
    var hoop = 0.0;
    var hoop_eps = 0.0;
    for i in 3 {
        for j in 3 {
            hoop += et[i] * n[3 * i + j] * et[j];
            hoop_eps += et[i] * eps[3 * i + j] * et[j];
        }
    }
    let axial = n[0];
    expect_near("hoop stress against pR/t", hoop / h, p * r / h, 0.005);
    expect_near("axial stress against pR/(2t)", axial / h, p * r / (2.0 * h), 0.005);
    let a = shells::membrane_compliance(s);
    expect_near("hoop strain against A⁻¹ (pR/2, pR)", hoop_eps, a[3] * p * r / 2.0 + a[4] * p * r, 0.005);
}

// The force per length f along x on the side of cell row j of grid g at
// the x lattice line la: a 2-node edge for order 1, a 3-node edge for 2.
fn closure_load(g: i64, la: i64, order: i64, j: i64, f: f64) {
    let b = order * j;
    if order == 1 {
        shells::edge_load(shells::grid_node(g, la, b), shells::grid_node(g, la, b + 1), f, 0.0, 0.0);
    } else {
        shells::edge_load3(shells::grid_node(g, la, b), shells::grid_node(g, la, b + 1), shells::grid_node(g, la, b + 2), f, 0.0, 0.0);
    }
}

fn element(order: i64) -> u8[>..] { return if order == 1 { str("MITC4") } else { str("MITC9") }; }

// --- the offset beam --------------------------------------------------------------

// A cantilever of 4 Timoshenko beams on the x axis, L = 1 m, with its axis
// at z = -e below the nodes. A tip force P along z bends it in x-z:
// w = PL³/(3EIy) + PL/(kGA). A tip force P along x at the node gives the
// axial force P and the moment eP: w = -eP L²/(2EIy).
fn offset_cantilever(fx: f64, fz: f64) -> f64, f64 {
    shells::clear();
    let b = shells::beam_section(70.0e9, 0.3, 2.0e-4, 4.0e-8, 1.0e-8, 1.0e-9);
    for i in 5 { shells::add_node(0.25 * i, 0.0, 0.0); }
    for i in 4 { shells::add_beam(i, i + 1, b, [0.0, 1.0, 0.0], [0.0, 0.0, -0.02], [0.0, 0.0, -0.02]); }
    for d in 6 { shells::fix(0, d); }
    shells::force(4, shells::UX, fx);
    shells::force(4, shells::UZ, fz);
    expect_solve("offset cantilever, 4 Timoshenko beams");
    expect_equilibrium("offset cantilever", abs(fx) + abs(fz));
    return shells::displacement(4, shells::UZ), shells::beam_axial_force(3);
}

fn test_offset_beam() {
    let e = 70.0e9;
    let g = e / 2.6;
    let p = 1000.0;
    let w1, n1 = offset_cantilever(0.0, p);
    expect_near("tip w against PL³/(3EIy) + PL/(kGA)", w1, p / (3.0 * e * 4.0e-8) + p / (5.0 / 6.0 * g * 2.0e-4), 1.0e-9);
    let w2, n2 = offset_cantilever(p, 0.0);
    expect_near("tip w against -eP L²/(2EIy)", w2, -0.02 * p / (2.0 * e * 4.0e-8), 1.0e-9);
    expect_near("Nx against P", n2, p, 1.0e-9);
}

// --- the thin Pagano plate --------------------------------------------------------

// Pagano (1970), "Exact solutions for rectangular bidirectional composites
// and sandwich plates", J. Composite Materials 4, 20-34: a square [0/90/0]
// plate, a = 100 mm, h = 1 mm (a/h = 100), simply supported, under
// q0 sin(pi x/a) sin(pi y/a) with q0 = 1e-3 MPa. The ply has E_L = 175000 MPa,
// E_T = 7000 MPa, G_LT = 3500 MPa, G_TT = 1400 MPa, nu = 0.25. The exact
// 3D solution (femahoy validation/cases/laminate-pagano, derivation.md)
// gives w̄ = 100 E_T h³ w/(q0 a⁴) = 0.4347 at the centre and
// σ̄_xx = σ_xx/(q0 (a/h)²) = 0.5393 at the centre of the bottom face.
// The quarter 0 <= x, y <= 50 mm.
fn test_pagano_thin() {
    let n = 20;
    let a = 100.0;
    let h = 1.0;
    let q0 = 1.0e-3;
    let el = 175000.0;
    let et = 7000.0;
    shells::clear();
    shells::laminate_clear();
    for ang in [0.0, 90.0, 0.0] { shells::laminate_ply(radians(ang), h / 3.0, el, et, 3500.0, 0.25, 3500.0, 1400.0); }
    let s = shells::laminate_section();
    for j in n + 1 {
        for i in n + 1 {
            let node = shells::add_node(0.5 * a * i / n, 0.5 * a * j / n, 0.0);
            shells::set_director(node, 0.0, 0.0, 1.0);
        }
    }
    for j in n {
        for i in n {
            let k = i + (n + 1) * j;
            let c = shells::add_shell(k, k + 1, k + n + 2, k + n + 1, s, 1.0, 0.0, 0.0);
            let q = q0 * sin(PI * (i + 0.5) / (2.0 * n)) * sin(PI * (j + 0.5) / (2.0 * n));
            shells::shell_load(c, 0.0, 0.0, -q, false);
        }
    }
    for k in n + 1 {
        fix_all(k * (n + 1), [shells::UY, shells::UZ]);
        fix_all(k, [shells::UX, shells::UZ]);
        fix_all(k * (n + 1) + n, [shells::UX, shells::RY, shells::RZ]);
        fix_all(n * (n + 1) + k, [shells::UY, shells::RX, shells::RZ]);
    }
    expect_solve(str("Pagano [0/90/0] plate a/h = 100, ", n, " x ", n, " MITC4"));
    let centre = n * (n + 1) + n;
    let w = shells::displacement(centre, shells::UZ);
    expect_near("Pagano w̄ at the centre", -100.0 * et * h * h * h * w / (q0 * a * a * a * a), 0.4347, 0.01);
    // The bottom ply is at 0 deg: σ_xx = Q11 (eps_11 - h/2 kappa_11) + Q12 (eps_22 - h/2 kappa_22),
    // at the corner of the centre cell.
    let eps, kap = shells::shell_strains(n * n - 1, 1.0, 1.0);
    let nu21 = 0.25 * et / el;
    let q11 = el / (1.0 - 0.25 * nu21);
    let q12 = 0.25 * et / (1.0 - 0.25 * nu21);
    let sx = q11 * (eps[0] - 0.5 * h * kap[0]) + q12 * (eps[1] - 0.5 * h * kap[1]);
    expect_near("Pagano σ̄_xx at the bottom centre", sx / (q0 * (a / h) * (a / h)), 0.5393, 0.01);
}

// --- the 787 barrel -------------------------------------------------------------

// The global model of docs/studies/b787/problem.md on a structured grid, as
// docs/studies/b787/skein/fuselage.skn states it: the quarter barrel
// 0 <= x <= 6.096 m, 0 <= θ <= 180 deg, composite skin and pad-up shells,
// offset beams for 45 stringers and 10 frames, ultimate pressure, the end
// closure load and the plug door load. The grid has the stringers on its
// θ lines and the frames on its x lines. The door cutout is the cells with
// their centre in the cutout, so its edge is a staircase on the grid; the
// pad-up is the cells with their centre within 0.150 m of the cutout.
let B_R = 2.883;
let B_XEND = 6.096;
let B_PU = 97.5e3;
let B_WD = 1.067;
let B_HD = 1.829;
let B_RC = 0.178;
let B_THD = 75.0;       // deg
let B_PAD = 0.150;

// True if (x, s) is in the cutout grown by g in the x-s plane.
fn in_door(x: f64, s: f64, g: f64) -> bool {
    let a = B_WD / 2.0 + g;
    let b = B_HD / 2.0 + g;
    let r = B_RC + g;
    let dx = abs(x) - a;
    let ds = abs(s - B_R * radians(B_THD)) - b;
    if dx > 0.0 || ds > 0.0 { return false; }
    let cx = dx + r;
    let cs = ds + r;
    return cx <= 0.0 || cs <= 0.0 || cx * cx + cs * cs <= r * r;
}

fn e_r(t: f64) -> f64[3] { return [0.0, sin(t), cos(t)]; }

fn scaled(v: f64[3], k: f64) -> f64[3] { return [k * v[0], k * v[1], k * v[2]]; }

var b_stringers: i64[>..] = [];
var b_frames: i64[>..] = [];
var b_hole = 0;

// The barrel with nx x nt cells of order 1 or 2 (mesh_cylinder_lines): nx
// a multiple of 20 (frames) and nt a multiple of 90 (stringers). Returns
// the grid. A beam line on a grid of order 2 has 2 beams on each cell side.
fn b787_model(nx: i64, nt: i64, order: i64) -> i64 {
    shells::clear();
    let skin = quasi_isotropic(3);
    let pad_sec = quasi_isotropic(6);
    let dx = B_XEND / nx;
    let dt = PI / nt;
    var secs: i64[>..] = [];
    b_hole = 0;
    for j in nt {
        for i in nx {
            let xc = dx * (i + 0.5);
            let sc = B_R * dt * (j + 0.5);
            let cut = in_door(xc, sc, 0.0);
            if cut { b_hole++; }
            secs.push(if cut { -1 } else if in_door(xc, sc, B_PAD) { pad_sec } else { skin });
        }
    }
    let g = shells::mesh_cylinder_lines(B_R, lines(0.0, B_XEND, nx), lines(0.0, PI, nt), order, secs, 1.0, 0.0, 0.0);
    for c in shells::shells.len { shells::shell_load(c, B_PU, 0.0, 0.0, true); }
    b787_door(g, nx, nt, order);
    b787_beams(g, nx, nt, order);
    // Symmetry at x = 0 and y = 0, one point at the bottom along z, and
    // the end closure load p_U R/2 per length at x = 6.096 m.
    let la = order * nx;
    let lb = order * nt;
    for j in lb + 1 { fix_all(shells::grid_node(g, 0, j), [shells::UX, shells::RY, shells::RZ]); }
    for i in la + 1 {
        fix_all(shells::grid_node(g, i, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(shells::grid_node(g, i, lb), [shells::UY, shells::RX, shells::RZ]);
    }
    shells::fix(shells::grid_node(g, 0, lb), shells::UZ);
    for j in nt { closure_load(g, la, order, j, B_PU * B_R / 2.0); }
    return g;
}

// The door load: the pressure of the cutout cells of row j, on the aft
// edge of the cutout in the row (problem.md §6.2), with the weights of the
// edge (1/2, 1/2 or 1/6, 4/6, 1/6).
fn b787_door(g: i64, nx: i64, nt: i64, order: i64) {
    let dx = B_XEND / nx;
    let dt = PI / nt;
    let w = if order == 1 { [0.5, 0.5, 0.0] } else { [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0] };
    for j in nt {
        var aft = -1;
        var door: f64[3] = [0.0, 0.0, 0.0];
        for i in nx {
            if !in_door(dx * (i + 0.5), B_R * dt * (j + 0.5), 0.0) { continue; }
            aft = i + 1;
            let t0 = dt * j;
            let t1 = dt * (j + 1);
            for k in 3 { door[k] += B_PU * dx * B_R * [0.0, cos(t0) - cos(t1), sin(t1) - sin(t0)][k]; }
        }
        if aft < 0 { continue; }
        for e in order + 1 {
            for k in 3 { shells::force(shells::grid_node(g, order * aft, order * j + e), k, w[e] * door[k]); }
        }
    }
}

// Stringers at θ = 2 deg + 4 deg k, frames at x = 0.3048 m + 0.6096 m k,
// offset inward along the radius, except in the cutout: one beam between
// each two next nodes of the lattice on the line.
fn b787_beams(g: i64, nx: i64, nt: i64, order: i64) {
    let stringer = shells::beam_section(55.85e9, 0.326, 2.016e-4, 4.908e-8, 3.639e-8, 1.569e-10);
    let frame = shells::beam_section(55.85e9, 0.326, 7.148e-4, 8.305e-8, 2.237e-6, 2.225e-9);
    let dx = B_XEND / nx;
    let dt = PI / nt;
    b_stringers.clear();
    b_frames.clear();
    for k in 45 {
        let j = (nt / 90) * (1 + 2 * k);
        let o = scaled(e_r(dt * j), -0.01750);
        for i in nx {
            if in_door(dx * (i + 0.5), B_R * dt * j, 0.0) { continue; }
            for e in order {
                let a = order * i + e;
                b_stringers.push(shells::add_beam(shells::grid_node(g, a, order * j), shells::grid_node(g, a + 1, order * j), stringer, e_r(dt * j), o, o));
            }
        }
    }
    for k in 10 {
        let i = (nx / 20) * (1 + 2 * k);
        for j in nt {
            if in_door(dx * i, B_R * dt * (j + 0.5), 0.0) { continue; }
            for e in order {
                let b = order * j + e;
                let t0 = dt * b / order;
                let t1 = dt * (b + 1) / order;
                b_frames.push(shells::add_beam(shells::grid_node(g, order * i, b), shells::grid_node(g, order * i, b + 1), frame,
                                               e_r(0.5 * (t0 + t1)), scaled(e_r(t0), -0.07415), scaled(e_r(t1), -0.07415)));
            }
        }
    }
}

fn test_b787(nx: i64, nt: i64, order: i64) {
    let g = b787_model(nx, nt, order);
    let t_skin = shells::sections[0].h;
    expect_solve(str("787 barrel ", nx, " x ", nt, " ", element(order), ", ", b_hole, " cutout cells, ",
        b_stringers.len, " stringer and ", b_frames.len, " frame beams"));
    // The reactions by axis: only x = 0 takes x, only y = 0 takes y, only
    // the bottom point takes z.
    let r = shells::reactions();
    var sum: f64[3] = [0.0, 0.0, 0.0];
    for d in r.len { if d % 6 < 3 { sum[d % 6] += r[d]; } }
    // The membrane stress at P (x = 4.8768 m, θ = 140 deg): the mean of the
    // 4 cells at the node.
    let ip = floor(4.8768 / (B_XEND / nx) + 0.5) as i64;
    let jp = floor(140.0 / 180.0 * nt + 0.5) as i64;
    let tp = radians(140.0);
    let et: f64[3] = [0.0, cos(tp), -sin(tp)];
    var hoop = 0.0;
    var axial = 0.0;
    for a in 2 {
        for b in 2 {
            let n, eps = shells::shell_resultants(shells::grid_cell(g, ip - 1 + a, jp - 1 + b));
            for i in 3 { for j in 3 { hoop += 0.25 * et[i] * n[3 * i + j] * et[j]; } }
            axial += 0.25 * n[0];
        }
    }
    let wr_max, s_max, f_max = b787_maxima();
    let s_hoop = hoop / t_skin;
    let s_axial = axial / t_skin;
    let hi = B_PU * B_R / t_skin;
    let lo = B_PU * B_R / (t_skin + 7.148e-4 / 0.6096);
    let f_end = B_PU * PI * B_R * B_R / 2.0;
    print("  w_r_max ", wr_max, " m, stringer |Nx| max ", s_max, " N, frame |Nx| max ", f_max, " N");
    print("  c1_ratio ", s_hoop / hi, ", c2_ratio ", s_axial / (B_PU * B_R / (2.0 * t_skin)));
    print("  C1 hoop stress at P: ", s_hoop, " Pa, range ", lo, " to ", hi, " Pa");
    if s_hoop < lo || s_hoop > hi {
        failures++;
        print("FAIL C1: the hoop stress at P is out of the range");
    }
    // C2: the sheet value (problem.md §7.3): the same axial strain in the
    // skin and the stringers, with the Poisson term of the skin at the
    // sheet hoop stress p_U R/t and nu = 0.326.
    let smear = 2.016e-4 / (2.0 * PI * B_R / 90.0);
    let sheet_c2 = (B_PU * B_R / 2.0 + 0.326 * hi * smear) / (t_skin + smear);
    print("  C2 axial stress at P: ", s_axial, " Pa against the sheet ", sheet_c2, " Pa, ratio ", s_axial / sheet_c2);
    expect_near("C2 axial stress at P against the sheet value", s_axial, sheet_c2, 0.05);
    expect_near("C3 sum of x reactions against -p pi R²/2", sum[0], -f_end, 0.001);
    print("  C4 sum of z reactions: ", sum[2], " N (limit ", f_end / 1000.0, " N)");
    if abs(sum[2]) > f_end / 1000.0 {
        failures++;
        print("FAIL C4: the z reaction is not near 0");
    }
    expect_near("C5 sum of y reactions against -p 2R x_end", sum[1], -B_PU * 2.0 * B_R * B_XEND, 0.001);
    if order == 2 { expect_same_as_skyline(); }
}

// The largest radial displacement, and the largest |Nx| of the stringers
// and of the frames.
fn b787_maxima() -> f64, f64, f64 {
    var wr_max = 0.0;
    for n in shells::node_count() {
        let p = shells::node_point(n);
        let t = atan2(p[1], p[2]);
        wr_max = max(wr_max, shells::displacement(n, shells::UY) * sin(t) + shells::displacement(n, shells::UZ) * cos(t));
    }
    var s_max = 0.0;
    for b in b_stringers { s_max = max(s_max, abs(shells::beam_axial_force(b))); }
    var f_max = 0.0;
    for b in b_frames { f_max = max(f_max, abs(shells::beam_axial_force(b))); }
    return wr_max, s_max, f_max;
}

// --- the nested dissection ---------------------------------------------------------

// The solve of the model now, which must be the nested dissection
// (dissect.goose), against the skyline: a residual |K u - f| / |f| below
// 1e-10, and the same displacements to 1e-9 of the largest.
fn expect_same_as_skyline() {
    guard shells::last_method == shells::SOLVE_DISSECTION else { abort("expect_same_as_skyline: the last solve is the dissection"); }
    let res = shells::relative_residual();
    print("  nested dissection: ", shells::front_entries, " entries, relative residual ", res);
    if !(res < 1.0e-10) {
        failures++;
        print("FAIL dissection: the relative residual is not below 1e-10");
    }
    var u: f64[>..] = [];
    for v in shells::disp { u.push(v); }
    shells::solve_method = shells::SOLVE_SKYLINE;
    expect_solve("  the same model with the skyline");
    shells::solve_method = shells::SOLVE_AUTO;
    var diff = 0.0;
    var size = 0.0;
    for d in u.len {
        diff = max(diff, abs(shells::disp[d] - u[d]));
        size = max(size, abs(u[d]));
    }
    print("  skyline: ", shells::profile_size(), " entries; largest difference of the displacements: ", diff / size, " of the largest");
    if diff > 1.0e-9 * size {
        failures++;
        print("FAIL dissection: the displacements differ from the skyline");
    }
}

// --- modes and buckling of the beams ------------------------------------------------

// Steel, SI units, and a rectangular section ty × tz.
let ST_E = 210.0e9;
let ST_NU = 0.3;
let ST_RHO = 7850.0;

fn steel_section(ty: f64, tz: f64) -> i64 {
    let iy = ty * tz * tz * tz / 12.0;
    let iz = tz * ty * ty * ty / 12.0;
    let s = shells::beam_section(ST_E, ST_NU, ty * tz, iy, iz, iy + iz);
    shells::beam_density(s, ST_RHO);
    return s;
}

// n beams of section s from the origin to len dir, with y axis y and the
// offset o at both ends of each beam.
fn beam_line(n: i64, len: f64, dir: f64[3], y: f64[3], o: f64[3], s: i64) {
    for i in n + 1 {
        let x = len * i / n;
        shells::add_node(x * dir[0], x * dir[1], x * dir[2]);
    }
    for i in n { shells::add_beam(i, i + 1, s, y, o, o); }
}

fn expect_eigen(what: u8[:], ok: bool) {
    if !ok {
        failures++;
        print("FAIL ", what, ": ", shells::eigen_error);
    }
    print(what, ": ", shells::dof_count(), " dofs, ", shells::eigen_iterations, " iterations, ", shells::last_seconds, " s");
}

// The first root b of the frequency equation of a Timoshenko cantilever
// (T. C. Huang, J. Appl. Mech. 28, 1961, 579-584; S. M. Han, H. Benaroya and
// T. Wei, J. Sound Vib. 225, 1999, 935-988, eq. (75)), with b² = ρAL⁴ω²/(EI),
// r2 = I/(AL²) and s2 = EI/(kGAL²): bisection on 0.5 .. 3.6.
fn timoshenko_cantilever_b(r2: f64, s2: f64) -> f64 {
    var lo = 0.5;
    var hi = 3.6;
    let flo = huang(lo, r2, s2);
    guard flo * huang(hi, r2, s2) < 0.0 else { abort("timoshenko_cantilever_b: no root in 0.5 .. 3.6"); }
    for it in 200 {
        let mid = 0.5 * (lo + hi);
        if huang(mid, r2, s2) * flo > 0.0 { lo = mid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

fn huang(b: f64, r2: f64, s2: f64) -> f64 {
    let root = sqrt((r2 - s2) * (r2 - s2) + 4.0 / (b * b));
    let a = sqrt(0.5 * (root - r2 - s2));
    let c = sqrt(0.5 * (root + r2 + s2));
    return 2.0 + (b * b * (r2 - s2) * (r2 - s2) + 2.0) * 0.5 * (exp(b * a) + exp(-b * a)) * cos(b * c)
         - b * (r2 + s2) / sqrt(1.0 - b * b * r2 * s2) * 0.5 * (exp(b * a) - exp(-b * a)) * sin(b * c);
}

// The first frequency of a Timoshenko cantilever of section ty × tz, length
// len, in the x-z plane (I = Iy).
fn timoshenko_cantilever_f(ty: f64, tz: f64, len: f64) -> f64 {
    let a = ty * tz;
    let i = ty * tz * tz * tz / 12.0;
    let g = ST_E / (2.0 * (1.0 + ST_NU));
    let b = timoshenko_cantilever_b(i / (a * len * len), ST_E * i / (5.0 / 6.0 * g * a * len * len));
    return b * sqrt(ST_E * i / (ST_RHO * a * len * len * len * len)) / (2.0 * PI);
}

// The modes of a steel cantilever of n beams on the skew line dir, L = 1 m,
// 15 × 10 mm, with ey normal to dir. Returns the first 3 frequencies, and
// the parts along ez and ey of the tip translation of modes 1 and 2.
fn skew_cantilever(n: i64) -> f64[3], f64, f64 {
    shells::clear();
    let dir: f64[3] = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let ey: f64[3] = [-2.0 / sqrt(5.0), 1.0 / sqrt(5.0), 0.0];
    let ez = shells::cross(dir, ey);
    beam_line(n, 1.0, dir, ey, [0.0, 0.0, 0.0], steel_section(0.015, 0.01));
    for d in 6 { shells::fix(0, d); }
    let ok = shells::modes(3);
    expect_eigen(str("skew steel cantilever, ", n, " Timoshenko beams, modal"), ok);
    var part: f64[2] = [0.0, 0.0];
    for k in 2 {
        shells::show_mode(k + 1);
        let t: f64[3] = [shells::displacement(n, 0), shells::displacement(n, 1), shells::displacement(n, 2)];
        part[k] = abs(shells::dot3(t, if k == 0 { ez } else { ey })) / shells::norm(t);
    }
    return [shells::eigenvalue(1), shells::eigenvalue(2), shells::eigenvalue(3)], part[0], part[1];
}

// (a) and (d): bending in x-z (Iy) gives mode 1 and in x-y (Iz = 2.25 Iy)
// mode 2, against the Timoshenko root (Huang) and Euler-Bernoulli,
// 1.8751² √(EI/(ρAL⁴))/(2π). Mode 3 is the second x-z mode, 4.6941².
fn test_beam_modes() {
    let eb = 1.8751040687 * 1.8751040687 / (2.0 * PI) * sqrt(ST_E * 0.015 * 1.0e-6 / 12.0 / (ST_RHO * 1.5e-4));
    let ref1 = timoshenko_cantilever_f(0.015, 0.01, 1.0);
    let ref2 = timoshenko_cantilever_f(0.01, 0.015, 1.0);
    print("  Timoshenko f1 ", ref1, " Hz, Euler-Bernoulli ", eb, " Hz");
    var last = 1.0;
    for n in [1, 2, 4, 8, 16] {
        let f, along_z, along_y = skew_cantilever(n);
        let err = abs(f[0] - ref1) / ref1;
        print("  n = ", n, ": f1 ", f[0], " Hz, error ", 100.0 * err, " %");
        if !(err < last) {
            failures++;
            print("FAIL the error of f1 does not fall with n");
        }
        last = err;
        if n == 16 {
            expect_near("cantilever f1 (Iy) against Timoshenko", f[0], ref1, 1.0e-4);
            expect_near("cantilever f1 (Iy) against Euler-Bernoulli", f[0], eb, 1.0e-3);
            expect_near("cantilever f2 (Iz) against Timoshenko", f[1], ref2, 1.0e-4);
            expect_near("cantilever f2 (Iz) against 1.5 Euler-Bernoulli", f[1], 1.5 * eb, 1.0e-3);
            expect_near("cantilever f3 against 4.6941² / 1.8751² Euler-Bernoulli", f[2], 4.6940911330 * 4.6940911330 / (1.8751040687 * 1.8751040687) * eb, 2.0e-3);
            expect_near("mode 1 tip translation along ez", along_z, 1.0, 1.0e-6);
            expect_near("mode 2 tip translation along ey", along_y, 1.0, 1.0e-6);
        }
    }
    // A short beam, L/h = 5: the shear and the rotary inertia lower f1.
    shells::clear();
    beam_line(32, 0.1, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], steel_section(0.03, 0.02));
    for d in 6 { shells::fix(0, d); }
    let ok = shells::modes(1);
    expect_eigen("short steel cantilever L/h = 5, 32 Timoshenko beams, modal", ok);
    let short_eb = 1.8751040687 * 1.8751040687 / (2.0 * PI) * sqrt(ST_E * 0.02 * 0.02 / 12.0 / ST_RHO) / 0.01;
    print("  Euler-Bernoulli ", short_eb, " Hz");
    expect_near("short cantilever f1 against Timoshenko", shells::eigenvalue(1), timoshenko_cantilever_f(0.03, 0.02, 0.1), 1.0e-4);
}

// (c) A uniform offset e of a clamped cantilever moves all of it, so the
// frequencies stay. A simply supported beam held in x at both nodes gets
// an axial restraint at e from its axis, so f1 goes up.
fn offset_modes(e: f64, clamped: bool) -> f64[3] {
    shells::clear();
    beam_line(8, 1.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -e], steel_section(0.015, 0.01));
    if clamped {
        for d in 6 { shells::fix(0, d); }
    } else {
        fix_all(0, [shells::UX, shells::UY, shells::UZ, shells::RX]);
        fix_all(8, [shells::UX, shells::UY, shells::UZ]);
    }
    let ok = shells::modes(3);
    expect_eigen(str(if clamped { "clamped" } else { "simply supported" }, " beam, offset ", e, " m, modal"), ok);
    return [shells::eigenvalue(1), shells::eigenvalue(2), shells::eigenvalue(3)];
}

fn test_offset_modes() {
    let c0 = offset_modes(0.0, true);
    let c1 = offset_modes(0.05, true);
    for k in 3 { expect_near(str("clamped f", k + 1, " with offset against no offset"), c1[k], c0[k], 1.0e-9); }
    let s0 = offset_modes(0.0, false);
    let s1 = offset_modes(0.05, false);
    let eb = PI / 2.0 * sqrt(ST_E * 0.015 * 1.0e-6 / 12.0 / (ST_RHO * 1.5e-4));
    expect_near("simply supported f1 against π/2 √(EIy/(ρAL⁴))", s0[0], eb, 1.0e-3);
    print("  simply supported f1 with offset 0.05 m: ", s1[0], " Hz, without: ", s0[0], " Hz");
    if !(s1[0] > 1.01 * s0[0]) {
        failures++;
        print("FAIL the offset does not raise f1 of the simply supported beam");
    }
}

// (b) Euler buckling of a steel column on x, L = 1 m, 15 × 10 mm, under
// P0 = 1000 N: π² EIy/L² pinned-pinned, π² EIy/(4L²) fixed-free.
fn column_buckling(n: i64, pinned: bool) -> f64 {
    shells::clear();
    beam_line(n, 1.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], steel_section(0.015, 0.01));
    if pinned {
        fix_all(0, [shells::UX, shells::UY, shells::UZ, shells::RX]);
        fix_all(n, [shells::UY, shells::UZ]);
    } else {
        for d in 6 { shells::fix(0, d); }
    }
    shells::force(n, shells::UX, -1000.0);
    let what = str(if pinned { "pinned-pinned" } else { "fixed-free" }, " column, ", n, " Timoshenko beams");
    expect_solve(what);
    let ok = shells::buckling(2);
    expect_eigen(str(what, ", buckling"), ok);
    return shells::eigenvalue(1);
}

fn test_buckling() {
    let pe = PI * PI * ST_E * 0.015 * 1.0e-6 / 12.0 / 1000.0;
    for pinned in [true, false] {
        let want = if pinned { pe } else { pe / 4.0 };
        // The shear of the Timoshenko beam lowers the load (Engesser):
        // P = P_E/(1 + P_E/(kGA)). The mesh converges to it.
        let kga = 5.0 / 6.0 * ST_E / (2.0 * (1.0 + ST_NU)) * 1.5e-4;
        let engesser = want / (1.0 + 1000.0 * want / kga);
        print("  Euler ", want, ", Engesser ", engesser);
        var last = 1.0;
        for n in [1, 2, 4, 8, 16] {
            if pinned && n == 1 { continue; }
            let got = column_buckling(n, pinned);
            let err = abs(got - engesser) / engesser;
            print("  n = ", n, ": load factor ", got, ", error against Engesser ", 100.0 * err, " %");
            if !(err < last) {
                failures++;
                print("FAIL the error of the load factor does not fall with n");
            }
            last = err;
            if n == 8 { expect_near(if pinned { "pinned-pinned load factor against π² EI/(L² P0)" } else { "fixed-free load factor against π² EI/(4 L² P0)" }, got, want, 0.01); }
            if n == 16 { expect_near("load factor against Engesser", got, engesser, 1.0e-5); }
        }
    }
}

fn main() {
    test_laminate();
    test_offset_beam();
    test_beam_modes();
    test_offset_modes();
    test_buckling();
    test_thin_cylinder(1);
    test_thin_cylinder(2);
    test_scordelis_lo();
    test_pinched_cylinder();
    test_pagano_thin();
    test_b787(40, 90, 1);
    test_b787(20, 90, 2);
    if failures > 0 { exit(1); }
    print("shells: all tests pass");
}
