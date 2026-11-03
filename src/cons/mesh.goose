// The mesh and the state of the conservation laws: the triangles that
// skein run meshes with src/mesher/, the sides of each cell, and the
// state of N unknown components in each cell.
//
// Cell c has the vertices tri[3c], tri[3c+1] and tri[3c+2]. Side k of the
// cell is the edge from vertex k to vertex (k + 1) % 3. Its outward
// normal times its length is (side_nx, side_ny). A side on the boundary
// has no neighbour (-1), and the generator gives it a mask of the
// conditions that hold on it (side_mask, bit b for the bc line b).
import std;
import math;
import os;
namespace cons;

let FV = 0;
let UPWIND = 0;      // the flux of the system: the Rusanov flux with the spectral radius of the Jacobian
let NONE = 0;
let MINMOD = 1;
let VAN_LEER = 2;
let FREE = 0;

var cell_total = 0;
var vertex_total = 0;
var vx: f64[>..] = [];
var vy: f64[>..] = [];
var tri: i64[>..] = [];
var centre_x: f64[>..] = [];
var centre_y: f64[>..] = [];
var area: f64[>..] = [];
var neighbour: i64[>..] = [];     // 3 for each cell
var side_nx: f64[>..] = [];
var side_ny: f64[>..] = [];
var side_mx: f64[>..] = [];       // the midpoint of each side
var side_my: f64[>..] = [];
var side_mask: i64[>..] = [];
var cell_region: u8[>..] = [];    // region_count flags for each cell
var region_count = 0;

var comps = 0;                    // N: the components of the unknowns
var sol: f64[>..] = [];           // comps for each cell
var time = 0.0;
var limiter = VAN_LEER;

fn cell_count() -> i64 { return cell_total; }
fn dof_count() -> i64 { return cell_total * comps; }

// Reads the mesh that skein run wrote in dir: the vertices, the triangles,
// the masks of the sides and the flags of the cell regions. False if a file
// is missing or does not fit the others.
fn read_mesh(dir: const u8[:], regions: i64) -> bool {
    var a: u8[>..] = [];
    var b: u8[>..] = [];
    var c: u8[>..] = [];
    var d: u8[>..] = [];
    guard read_file(str(dir, "/cons.xy"), a) && read_file(str(dir, "/cons.tri"), b) && read_file(str(dir, "/cons.sides"), c) &&
          read_file(str(dir, "/cons.cells"), d) else { return false; }
    let xy, xy_ok = from_bytes<f64[>..]>(a);
    let tr, tr_ok = from_bytes<i64[>..]>(b);
    let mask, mask_ok = from_bytes<i64[>..]>(c);
    let flags, flags_ok = from_bytes<u8[>..]>(d);
    guard xy_ok && tr_ok && mask_ok && flags_ok && tr.len > 0 && tr.len % 3 == 0 else { return false; }
    guard mask.len == tr.len && flags.len == regions * (tr.len / 3) else { return false; }
    vx.clear();
    vy.clear();
    tri.clear();
    side_mask.clear();
    cell_region.clear();
    for i in xy.len / 2 {
        vx.push(xy[2 * i]);
        vy.push(xy[2 * i + 1]);
    }
    for t in tr { tri.push(t); }
    for m in mask { side_mask.push(m); }
    for f in flags { cell_region.push(f); }
    region_count = regions;
    vertex_total = vx.len;
    cell_total = tri.len / 3;
    for t in tri { guard t >= 0 && t < vertex_total else { return false; }; }
    build_sides();
    return true;
}

// The centre, the area, the normals and the neighbours of each cell.
fn build_sides() {
    centre_x.clear();
    centre_y.clear();
    area.clear();
    neighbour.clear();
    side_nx.clear();
    side_ny.clear();
    side_mx.clear();
    side_my.clear();
    var keys: i64[>..] = [];
    for c in cell_total {
        var sx = 0.0;
        var sy = 0.0;
        for k in 3 {
            let a = tri[3 * c + k];
            let b = tri[3 * c + (k + 1) % 3];
            sx += vx[a];
            sy += vy[a];
            side_nx.push(vy[b] - vy[a]);
            side_ny.push(vx[a] - vx[b]);
            side_mx.push(0.5 * (vx[a] + vx[b]));
            side_my.push(0.5 * (vy[a] + vy[b]));
            neighbour.push(-1);
            keys.push((min(a, b) * vertex_total + max(a, b)) * (3 * cell_total) + 3 * c + k);
        }
        centre_x.push(sx / 3.0);
        centre_y.push(sy / 3.0);
        let a = tri[3 * c];
        let b = tri[3 * c + 1];
        let d = tri[3 * c + 2];
        area.push(0.5 * abs((vx[b] - vx[a]) * (vy[d] - vy[a]) - (vx[d] - vx[a]) * (vy[b] - vy[a])));
    }
    sort(keys);
    let n3 = 3 * cell_total;
    for i in 1..keys.len {
        if keys[i] / n3 == keys[i - 1] / n3 {
            let p = keys[i - 1] % n3;
            let q = keys[i] % n3;
            neighbour[p] = q / 3;
            neighbour[q] = p / 3;
        }
    }
}

// Sets the number of components of the unknowns and the limiter, and sizes
// the state.
fn setup(n: i64, face_limiter: i64) {
    guard cell_total > 0 && n >= 1 && face_limiter >= NONE && face_limiter <= VAN_LEER else { abort("setup: a mesh, comps >= 1 and a limiter"); }
    comps = n;
    limiter = face_limiter;
    sol.clear();
    for d in cell_total * n { sol.push(0.0); }
    time = 0.0;
    size_solver();
}

// Sets the state to the field initial(x, y, k) of component k at the centre
// of each cell.
fn start<INIT>() {
    guard comps >= 1 && sol.len == dof_count() else { abort("start: setup"); }
    time = 0.0;
    for c in cell_total {
        for k in comps { sol[c * comps + k] = INIT(centre_x[c], centre_y[c], k); }
    }
}
