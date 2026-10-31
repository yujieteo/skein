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
// The quarter x in [0, 25], θ in [0, 40 deg].
fn scordelis_lo(n: i64) -> f64 {
    shells::clear();
    let s = shells::isotropic_section(4.32e8, 0.0, 0.25);
    let g = shells::mesh_cylinder(25.0, 0.0, 25.0, 0.0, radians(40.0), n, n, s, 1.0, 0.0, 0.0);
    for k in n + 1 {
        fix_all(shells::grid_node(g, 0, k), [shells::UX, shells::RY, shells::RZ]);
        fix_all(shells::grid_node(g, n, k), [shells::UY, shells::UZ]);
        fix_all(shells::grid_node(g, k, 0), [shells::UY, shells::RX, shells::RZ]);
    }
    for c in shells::shells.len { shells::shell_load(c, 0.0, 0.0, -90.0, false); }
    expect_solve(str("Scordelis-Lo roof ", n, " x ", n, " MITC4"));
    expect_equilibrium("Scordelis-Lo", 90.0 * 25.0 * 25.0);
    return -shells::displacement(shells::grid_node(g, 0, n), shells::UZ);
}

fn test_scordelis_lo() {
    for n in [8, 16] { print("  ", n, " x ", n, ": ", scordelis_lo(n), " ft"); }
    // MITC4 converges to the reference from below. The error of a 32 x 32
    // quarter is about 0.6 %, so 1 % finds a wrong term but not mesh error.
    expect_near("Scordelis-Lo uz at the free edge, 32 x 32", scordelis_lo(32), 0.3024, 0.01);
}

// --- the pinched cylinder --------------------------------------------------------

// MacNeal and Harder (1985): R = 300 in, L = 600 in, t = 3 in, E = 3e6 psi,
// nu = 0.3, rigid diaphragms at the ends, two forces P = 1 lb at mid-length.
// The radial displacement under a force is 1.8248e-5 in. The octant
// x in [0, 300] (x = 0 is a diaphragm), θ in [0, 90 deg].
fn pinched_cylinder(n: i64) -> f64 {
    shells::clear();
    let s = shells::isotropic_section(3.0e6, 0.3, 3.0);
    let g = shells::mesh_cylinder(300.0, 0.0, 300.0, 0.0, radians(90.0), n, n, s, 1.0, 0.0, 0.0);
    for k in n + 1 {
        fix_all(shells::grid_node(g, 0, k), [shells::UY, shells::UZ, shells::RX]);
        fix_all(shells::grid_node(g, n, k), [shells::UX, shells::RY, shells::RZ]);
        fix_all(shells::grid_node(g, k, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(shells::grid_node(g, k, n), [shells::UZ, shells::RX, shells::RY]);
    }
    shells::force(shells::grid_node(g, n, 0), shells::UZ, -0.25);
    expect_solve(str("pinched cylinder ", n, " x ", n, " MITC4"));
    expect_equilibrium("pinched cylinder", 0.25);
    return -shells::displacement(shells::grid_node(g, n, 0), shells::UZ);
}

fn test_pinched_cylinder() {
    for n in [8, 16] { print("  ", n, " x ", n, ": ", pinched_cylinder(n), " in"); }
    // The point load converges slowly. The error of a 32 x 32 octant is about
    // 1.1 %, so 2 % finds a wrong term but not mesh error.
    expect_near("pinched cylinder radial displacement, 32 x 32", pinched_cylinder(32), 1.8248e-5, 0.02);
}

// --- the thin cylinder under pressure ----------------------------------------------

// A closed thin cylinder with the 787 skin (problem.md §1, §3): R = 2.883 m,
// p = 97.5 kPa, the [45/0/-45/90]3s laminate. Membrane theory gives the
// hoop stress pR/t, the axial stress pR/(2t), and the hoop strain
// a21 Nx + a22 Nθ with a = A⁻¹. The quarter θ in [0, 90 deg] of a bay
// x in [0, 1.2192 m], with symmetry at x = 0 and the closure load pR/2 per
// length at x = 1.2192 m.
fn test_thin_cylinder() {
    let r = 2.883;
    let p = 97.5e3;
    let len = 1.2192;
    let nx = 8;
    let nt = 24;
    shells::clear();
    let s = quasi_isotropic(3);
    let h = shells::sections[s].h;
    let g = shells::mesh_cylinder(r, 0.0, len, 0.0, radians(90.0), nx, nt, s, 1.0, 0.0, 0.0);
    for j in nt + 1 { fix_all(shells::grid_node(g, 0, j), [shells::UX, shells::RY, shells::RZ]); }
    for i in nx + 1 {
        fix_all(shells::grid_node(g, i, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(shells::grid_node(g, i, nt), [shells::UZ, shells::RX, shells::RY]);
    }
    for c in shells::shells.len { shells::shell_load(c, p, 0.0, 0.0, true); }
    for j in nt { shells::edge_load(shells::grid_node(g, nx, j), shells::grid_node(g, nx, j + 1), p * r / 2.0, 0.0, 0.0); }
    expect_solve(str("thin cylinder ", nx, " x ", nt, " MITC4, 24-ply skin"));
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

fn barrel_node(nx: i64, i: i64, j: i64) -> i64 { return i + (nx + 1) * j; }

fn e_r(t: f64) -> f64[3] { return [0.0, sin(t), cos(t)]; }

fn scaled(v: f64[3], k: f64) -> f64[3] { return [k * v[0], k * v[1], k * v[2]]; }

// nx and nt are the cells along x and θ: nx a multiple of 20 (frames) and
// nt a multiple of 90 (stringers).
fn test_b787(nx: i64, nt: i64) {
    shells::clear();
    let skin = quasi_isotropic(3);
    let pad_sec = quasi_isotropic(6);
    let t_skin = shells::sections[skin].h;
    let stringer = shells::beam_section(55.85e9, 0.326, 2.016e-4, 4.908e-8, 3.639e-8, 1.569e-10);
    let frame = shells::beam_section(55.85e9, 0.326, 7.148e-4, 8.305e-8, 2.237e-6, 2.225e-9);
    let dx = B_XEND / nx;
    let dt = PI / nt;
    // The nodes, x first: the profile is about 6 (nx + 1) wide.
    for j in nt + 1 {
        for i in nx + 1 {
            let t = dt * j;
            let n = shells::add_node(dx * i, B_R * sin(t), B_R * cos(t));
            shells::set_director(n, 0.0, sin(t), cos(t));
        }
    }
    var cell: i64[>..] = [];
    var hole = 0;
    for j in nt {
        // The door load: the pressure of the cutout cells of row j, on the
        // aft edge of the cutout in the row (problem.md §6.2).
        var aft = -1;
        var door: f64[3] = [0.0, 0.0, 0.0];
        for i in nx {
            let xc = dx * (i + 0.5);
            let sc = B_R * dt * (j + 0.5);
            if in_door(xc, sc, 0.0) {
                cell.push(-1);
                hole++;
                aft = i + 1;
                let t0 = dt * j;
                let t1 = dt * (j + 1);
                for k in 3 { door[k] += B_PU * dx * B_R * [0.0, cos(t0) - cos(t1), sin(t1) - sin(t0)][k]; }
            } else {
                let sec = if in_door(xc, sc, B_PAD) { pad_sec } else { skin };
                let c = shells::add_shell(barrel_node(nx, i, j), barrel_node(nx, i + 1, j), barrel_node(nx, i + 1, j + 1), barrel_node(nx, i, j + 1), sec, 1.0, 0.0, 0.0);
                shells::shell_load(c, B_PU, 0.0, 0.0, true);
                cell.push(c);
            }
        }
        if aft >= 0 {
            for k in 3 {
                shells::force(barrel_node(nx, aft, j), k, 0.5 * door[k]);
                shells::force(barrel_node(nx, aft, j + 1), k, 0.5 * door[k]);
            }
        }
    }
    // Stringers at θ = 2 deg + 4 deg k, frames at x = 0.3048 m + 0.6096 m k,
    // offset inward along the radius, except in the cutout.
    var stringers: i64[>..] = [];
    var frames: i64[>..] = [];
    for k in 45 {
        let j = (nt / 90) * (1 + 2 * k);
        let t = dt * j;
        for i in nx {
            if in_door(dx * (i + 0.5), B_R * t, 0.0) { continue; }
            let o = scaled(e_r(t), -0.01750);
            stringers.push(shells::add_beam(barrel_node(nx, i, j), barrel_node(nx, i + 1, j), stringer, e_r(t), o, o));
        }
    }
    for k in 10 {
        let i = (nx / 20) * (1 + 2 * k);
        for j in nt {
            let t0 = dt * j;
            let t1 = dt * (j + 1);
            if in_door(dx * i, B_R * dt * (j + 0.5), 0.0) { continue; }
            frames.push(shells::add_beam(barrel_node(nx, i, j), barrel_node(nx, i, j + 1), frame, e_r(0.5 * (t0 + t1)),
                                         scaled(e_r(t0), -0.07415), scaled(e_r(t1), -0.07415)));
        }
    }
    // Symmetry at x = 0 and y = 0, one point at the bottom along z, and
    // the end closure load p_U R/2 per length at x = 6.096 m.
    for j in nt + 1 { fix_all(barrel_node(nx, 0, j), [shells::UX, shells::RY, shells::RZ]); }
    for i in nx + 1 {
        fix_all(barrel_node(nx, i, 0), [shells::UY, shells::RX, shells::RZ]);
        fix_all(barrel_node(nx, i, nt), [shells::UY, shells::RX, shells::RZ]);
    }
    shells::fix(barrel_node(nx, 0, nt), shells::UZ);
    for j in nt { shells::edge_load(barrel_node(nx, nx, j), barrel_node(nx, nx, j + 1), B_PU * B_R / 2.0, 0.0, 0.0); }
    expect_solve(str("787 barrel ", nx, " x ", nt, " MITC4, ", hole, " cutout cells, ",
        stringers.len, " stringer and ", frames.len, " frame beams"));
    // The reactions by axis: only x = 0 takes x, only y = 0 takes y, only
    // the bottom point takes z.
    let r = shells::reactions();
    var sum: f64[3] = [0.0, 0.0, 0.0];
    for d in r.len { if d % 6 < 3 { sum[d % 6] += r[d]; } }
    // The membrane stress at P (x = 4.8768 m, θ = 140 deg): the mean of the
    // 4 cells at the node.
    let ip = floor(4.8768 / dx + 0.5) as i64;
    let jp = floor(140.0 / 180.0 * nt + 0.5) as i64;
    let tp = radians(140.0);
    let et: f64[3] = [0.0, cos(tp), -sin(tp)];
    var hoop = 0.0;
    var axial = 0.0;
    for a in 2 {
        for b in 2 {
            let n, eps = shells::shell_resultants(cell[(ip - 1 + a) + nx * (jp - 1 + b)]);
            for i in 3 { for j in 3 { hoop += 0.25 * et[i] * n[3 * i + j] * et[j]; } }
            axial += 0.25 * n[0];
        }
    }
    var wr_max = 0.0;
    for j in nt + 1 {
        for i in nx + 1 {
            let t = dt * j;
            let w = shells::displacement(barrel_node(nx, i, j), shells::UY) * sin(t) + shells::displacement(barrel_node(nx, i, j), shells::UZ) * cos(t);
            wr_max = max(wr_max, w);
        }
    }
    var s_max = 0.0;
    for b in stringers { s_max = max(s_max, abs(shells::beam_axial_force(b))); }
    var f_max = 0.0;
    for b in frames { f_max = max(f_max, abs(shells::beam_axial_force(b))); }
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
}

fn main() {
    test_laminate();
    test_offset_beam();
    test_thin_cylinder();
    test_scordelis_lo();
    test_pinched_cylinder();
    test_pagano_thin();
    test_b787(40, 90);
    if failures > 0 { exit(1); }
    print("shells: all tests pass");
}
