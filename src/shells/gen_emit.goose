// The program of the shell generator: the prelude, the mesh, the beams, the
// steps, the reports and the expects, as Goose code that calls src/shells/.
import std;
import math;
namespace shells_gen;
import syntax from "../syntax";
import codegen.codegen;

fn lit(out: u8[>..]&, x: f64) { codegen::lit(out, x); }

fn lits(out: u8[>..]&, v: f64[3]) {
    lit(out, v[0]);
    out.append(", ");
    lit(out, v[1]);
    out.append(", ");
    lit(out, v[2]);
}

fn quoted(out: u8[>..]&, s: u8[:]) { codegen::quoted(out, s); }

fn f64s(out: u8[>..]&, v: f64[:]) {
    out.push('[');
    for x, i in v {
        if i > 0 { out.append(", "); }
        lit(out, x);
    }
    out.push(']');
}

fn i64s(out: u8[>..]&, v: i64[:]) {
    out.push('[');
    for x, i in v {
        if i > 0 { out.append(", "); }
        format(out, x);
    }
    out.push(']');
}

// The part of the generated program that is the same for each model.
let PRELUDE = """
    extern "libc.so.6" {
        fn system(command: const u8&) -> i32;
    }

    // The reports and the verification records, in TOON (command_line.md §6).
    var reports: u8[>..] = [];
    var report_count = 0;
    var checks: u8[>..] = [];
    var check_count = 0;
    var solves: u8[>..] = [];
    var solve_count = 0;
    var results: u8[>..] = [];
    var result_count = 0;
    var failed = false;
    var grid = -1;
    var secs: i64[>..] = [];        // the section of each material of the shell
    var run_first: i64[>..] = [];   // the beams of each run: a range of shells::beams
    var run_count: i64[>..] = [];
    var pres_dofs: i64[>..] = [];   // the dofs with a prescribed value, and the values
    var pres_vals: f64[>..] = [];

    // Appends `s` as a TOON value, in quotes when it needs them.
    fn toon(into: u8[>..]&, s: u8[:]) {
        var needs = s.len == 0;
        for c in s {
            if c == ',' || c == '"' || c == ':' || c == '\\' || c == ' ' { needs = true; }
        }
        if !needs { into.append(s); return; }
        into.push('"');
        for c in s {
            if c == '"' || c == '\\' { into.push('\\'); }
            into.push(c);
        }
        into.push('"');
    }

    fn quantity(x: f64, unit: u8[:]) -> u8[>..] { if unit.len == 0 { return str(x); } return str(x, " ", unit); }

    fn vector(v: f64[3], unit: u8[:]) -> u8[>..] { str("(", v[0], ", ", v[1], ", ", v[2], ") ", unit) }

    fn report(step: u8[:], expr: u8[:], value: u8[:]) {
        reports.append("  ");
        toon(reports, step);
        reports.push(',');
        toon(reports, expr);
        reports.push(',');
        toon(reports, value);
        reports.push('\n');
        report_count++;
    }

    fn check(id: u8[:], pass: bool, value: u8[:], limit: u8[:]) {
        checks.append("  post,");
        toon(checks, id);
        checks.append(if pass { ",pass," } else { ",fail," });
        toon(checks, value);
        checks.push(',');
        toon(checks, limit);
        checks.push('\n');
        check_count++;
        if !pass { failed = true; }
    }

    fn fault(label: u8[:], what: u8[:]) {
        write_stderr(str("skein run: the point of '", label, "' ", what, "\n"));
        failed = true;
    }

    // --- the values of the expressions ------------------------------------------

    fn pack1(x: f64) -> f64[9] {
        var r: f64[9] = [0.0; 9];
        r[0] = x;
        return r;
    }

    fn pack3(a: f64[3]) -> f64[9] {
        var r: f64[9] = [0.0; 9];
        for i in 3 { r[i] = a[i]; }
        return r;
    }

    fn pack9(a: f64[9]) -> f64[9] { return a; }
    fn at1(a: f64[9]) -> f64 { return a[0]; }
    fn at3(a: f64[9]) -> f64[3] { return [a[0], a[1], a[2]]; }
    fn comp(a: f64[3], i: i64) -> f64 { return a[i]; }
    fn neg1(x: f64) -> f64 { return -x; }
    fn neg3(a: f64[3]) -> f64[3] { return [-a[0], -a[1], -a[2]]; }
    fn add3(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
    fn sub3(a: f64[3], b: f64[3]) -> f64[3] { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
    fn scale3(s: f64, a: f64[3]) -> f64[3] { return [s * a[0], s * a[1], s * a[2]]; }
    fn dot3(a: f64[3], b: f64[3]) -> f64 { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
    fn cross3(a: f64[3], b: f64[3]) -> f64[3] { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
    fn norm3(a: f64[3]) -> f64 { return sqrt(dot3(a, a)); }

    fn neg9(a: f64[9]) -> f64[9] {
        var r = a;
        for i in 9 { r[i] = -a[i]; }
        return r;
    }

    fn add9(a: f64[9], b: f64[9]) -> f64[9] {
        var r = a;
        for i in 9 { r[i] += b[i]; }
        return r;
    }

    fn sub9(a: f64[9], b: f64[9]) -> f64[9] {
        var r = a;
        for i in 9 { r[i] -= b[i]; }
        return r;
    }

    fn scale9(s: f64, a: f64[9]) -> f64[9] {
        var r = a;
        for i in 9 { r[i] *= s; }
        return r;
    }

    // t v, for a tensor t (row-major) and a vector v.
    fn tdotv(t: f64[9], v: f64[3]) -> f64[3] {
        return [t[0] * v[0] + t[1] * v[1] + t[2] * v[2], t[3] * v[0] + t[4] * v[1] + t[5] * v[2],
                t[6] * v[0] + t[7] * v[1] + t[8] * v[2]];
    }

    // The value of material m: a property of the materials of the shell.
    fn pick(m: i64, values: f64[:]) -> f64 { return values[clamp(m, 0, values.len - 1)]; }

    // The outward unit normal of the cylinder at p.
    fn radial(p: f64[3]) -> f64[3] {
        let r = sqrt(p[1] * p[1] + p[2] * p[2]);
        return [0.0, p[1] / r, p[2] / r];
    }

    // --- the views ---------------------------------------------------------------

    // A view holds the solution of a model at its nodes, for the reports
    // and for a submodel: view 0 is this model, view 1 the model of the
    // results value. It is [order, nx, nt, r, nodes, the x lines, the θ
    // lines, 1 or 0 for each cell (i + nx j) that is kept], then VIEW_SIZE values for each node: x, y and z, the 6 dofs,
    // N, eps_m and the material. The nodes of a grid are its lattice
    // points, a + (order nx + 1) b. A view of beams has order 0 and the
    // nodes of the model.
    let VIEW_SIZE = 28;
    var own: f64[>..] = [];
    var other: f64[>..] = [];
    var node_view: i64[>..] = [];   // the view node of each node of the model

    fn view(v: i64) -> const f64[:] {
        if v == 0 { return own[..]; }
        return other[..];
    }

    fn v_order(v: i64) -> i64 { return view(v)[0] as! i64; }
    fn v_nx(v: i64) -> i64 { return view(v)[1] as! i64; }
    fn v_nt(v: i64) -> i64 { return view(v)[2] as! i64; }
    fn v_count(v: i64) -> i64 { return view(v)[4] as! i64; }
    fn v_base(v: i64) -> i64 { return if v_order(v) == 0 { 5 } else { 7 + v_nx(v) + v_nt(v) + v_nx(v) * v_nt(v) }; }
    fn v_kept(v: i64, i: i64, j: i64) -> bool { return view(v)[7 + v_nx(v) + v_nt(v) + i + v_nx(v) * j] > 0.0; }
    fn v_val(v: i64, k: i64, c: i64) -> f64 { return view(v)[v_base(v) + VIEW_SIZE * k + c]; }
    fn v_point(v: i64, k: i64) -> f64[3] { return [v_val(v, k, 0), v_val(v, k, 1), v_val(v, k, 2)]; }
    fn v_w(v: i64, k: i64) -> f64[3] { return [v_val(v, k, 3), v_val(v, k, 4), v_val(v, k, 5)]; }
    fn v_th(v: i64, k: i64) -> f64[3] { return [v_val(v, k, 6), v_val(v, k, 7), v_val(v, k, 8)]; }
    fn v_mat(v: i64, k: i64) -> i64 { return v_val(v, k, 27) as! i64; }

    fn v_tn(v: i64, k: i64) -> f64[9] {
        var r: f64[9] = [0.0; 9];
        for c in 9 { r[c] = v_val(v, k, 9 + c); }
        return r;
    }

    fn v_te(v: i64, k: i64) -> f64[9] {
        var r: f64[9] = [0.0; 9];
        for c in 9 { r[c] = v_val(v, k, 18 + c); }
        return r;
    }

    // The node of lattice point k of the grid.
    fn lnode(k: i64) -> i64 {
        let la = shells::grids[grid].order * shells::grids[grid].nx + 1;
        return shells::grid_node(grid, k % la, k / la);
    }

    fn lk(a: i64, b: i64) -> i64 { return a + (shells::grids[grid].order * shells::grids[grid].nx + 1) * b; }

    fn to_node(k: i64) -> i64 { return if grid >= 0 { lnode(k) } else { k }; }

    // The (xi, eta) of the nodes of a cell: the corners, the mid-sides and
    // the centre (shells::Shell).
    let CELL_XI = [-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
    let CELL_ETA = [-1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, 0.0];

    // View 0 from the solution: N and eps_m at a node are the mean over
    // the shell cells at the node of the values of each cell at the node.

    fn make_own() {
        own.clear();
        var count = shells::node_count();
        if grid >= 0 {
            let g = shells::grids[grid];
            count = (g.order * g.nx + 1) * (g.order * g.nt + 1);
            own.append([1.0 * g.order, 1.0 * g.nx, 1.0 * g.nt, g.r, 1.0 * count]);
            for i in g.nx + 1 { own.push(shells::grid_line(grid, 0, i)); }
            for j in g.nt + 1 { own.push(shells::grid_line(grid, 1, j)); }
            for j in g.nt { for i in g.nx { own.push(if shells::grid_cell(grid, i, j) >= 0 { 1.0 } else { 0.0 }); } }
        } else {
            own.append([0.0, 0.0, 0.0, 0.0, 1.0 * count]);
        }
        let base = own.len;
        node_view.clear();
        for n in shells::node_count() { node_view.push(-1); }
        var cells_at: i64[>..] = [];
        for k in count {
            let n = to_node(k);
            node_view[n] = k;
            let p = shells::node_point(n);
            for c in 3 { own.push(p[c]); }
            for d in 6 { own.push(shells::displacement(n, d)); }
            for c in 19 { own.push(0.0); }
            cells_at.push(0);
        }
        // The resultants at the 2 x 2 Gauss points of each cell, taken
        // bilinearly out to its nodes.
        let g = 1.0 / sqrt(3.0);
        let gauss = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        var at: f64[72] = [0.0; 72];
        for c, ci in shells::shells {
            for q in 4 {
                let tn, te = shells::shell_resultants_at(ci, g * gauss[q][0], g * gauss[q][1]);
                for e in 9 {
                    at[18 * q + e] = tn[e];
                    at[18 * q + 9 + e] = te[e];
                }
            }
            let m = position(secs[..], c.section);
            for j in c.count {
                let k = node_view[c.n[j]];
                let r = sqrt(3.0) * CELL_XI[j];
                let s = sqrt(3.0) * CELL_ETA[j];
                for q in 4 {
                    let w = 0.25 * (1.0 + r * gauss[q][0]) * (1.0 + s * gauss[q][1]);
                    for e in 18 { own[base + VIEW_SIZE * k + 9 + e] += w * at[18 * q + e]; }
                }
                own[base + VIEW_SIZE * k + 27] = 1.0 * m;
                cells_at[k]++;
            }
        }
        for k in count {
            if cells_at[k] < 2 { continue; }
            for e in 18 { own[base + VIEW_SIZE * k + 9 + e] /= cells_at[k]; }
        }
    }

    fn save_view(dir: u8[:], step: u8[:]) {
        if !write_file(str(dir, "/", step, ".shell_view"), to_bytes(own[..])) {
            write_stderr(str("skein run: cannot write the state of step '", step, "' in '", dir, "'\n"));
            failed = true;
        }
    }

    fn load_view(dir: u8[:], step: u8[:]) -> bool {
        var b: u8[>..] = [];
        guard read_file(str(dir, "/", step, ".shell_view"), b) else { return false; };
        let a, good = from_bytes<f64[>..]>(b);
        other.clear();
        other.append(a);
        return good && other.len > 5;
    }

    // View 1: the state of step `step` of the model `file`. If the state is
    // not there, `command` runs that model first.
    fn results(dir: u8[:], step: u8[:], file: u8[:], command: u8[:]) {
        if load_view(dir, step) { return; }
        var c: u8[>..] = [];
        c.append(command);
        c.push(0);
        if system(c[0]) != 0 { write_stderr(str("skein run: the run of '", file, "' for its results failed\n")); }
        if !load_view(dir, step) {
            write_stderr(str("skein run: no state of step '", step, "' of '", file, "' in '", dir, "'\n"));
            exit(1);
        }
    }

    fn line_index(v: const f64[:], x: f64) -> i64 {
        var lo = 0;
        var hi = v.len - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if v[mid] <= x { lo = mid; } else { hi = mid; }
        }
        return lo;
    }

    fn quad_weights(s: f64) -> f64[3] { return [0.5 * s * (s - 1.0), 1.0 - s * s, 0.5 * s * (s + 1.0)]; }

    // The view nodes and the weights of point p in view v: the nodes of the
    // cell that holds p, or the nearest node of a view of beams. The count
    // is 0 if p is farther than tol from the view.
    fn locate(v: i64, p: f64[3], tol: f64, label: u8[:]) -> i64[9], f64[9], i64 {
        var nodes: i64[9] = [0; 9];
        var w: f64[9] = [0.0; 9];
        let o = v_order(v);
        if o == 0 {
            var best = -1;
            var d = 0.0;
            for k in v_count(v) {
                let e = norm3(sub3(v_point(v, k), p));
                if best < 0 || e < d {
                    best = k;
                    d = e;
                }
            }
            if best < 0 || d > tol {
                fault(label, "is not a node of the mesh");
                return nodes, w, 0;
            }
            nodes[0] = best;
            w[0] = 1.0;
            return nodes, w, 1;
        }
        let s = view(v);
        let nx = v_nx(v);
        let nt = v_nt(v);
        let r = s[3];
        let xl = s[5..6 + nx];
        let tl = s[6 + nx..7 + nx + nt];
        let span = tl[nt] - tl[0];
        let d = arc_of(p, tl[0], span);
        let off_x = max(max(xl[0] - p[0], p[0] - xl[nx]), 0.0);
        let off_t = r * max(max(-d, d - span), 0.0);
        if off_x > tol || off_t > tol || abs(sqrt(p[1] * p[1] + p[2] * p[2]) - r) > tol {
            fault(label, "is not on the shell");
            return nodes, w, 0;
        }
        // The nodes move onto the cuts, so the cell of the lines is only a
        // start: the kept cell near it whose map holds p, or comes nearest.
        let i0 = line_index(xl, p[0]);
        let j0 = line_index(tl, tl[0] + d);
        let best_i, best_j, bxi, beta, best = nearest_cell(v, i0, j0, p[0], r * d, tl[0], span);
        if best_i < 0 || best > 0.5 * (xl[min(i0 + 1, nx)] - xl[i0] + r * (tl[min(j0 + 1, nt)] - tl[j0])) {
            fault(label, "is not on the mesh of the shell");
            return nodes, w, 0;
        }
        let la = o * nx + 1;
        let lx, dx = shape1(o, bxi);
        let lt, dt = shape1(o, beta);
        var c = 0;
        for b in o + 1 {
            for a in o + 1 {
                nodes[c] = o * best_i + a + la * (o * best_j + b);
                w[c] = lx[a] * lt[b];
                c++;
            }
        }
        return nodes, w, c;
    }

    // The kept cell within 2 cells of cell (i0, j0) of view v that holds
    // the point (x, s), or comes nearest to it: the cell, its (xi, eta) and
    // the distance. The cell is -1 if no cell is kept there.
    fn nearest_cell(v: i64, i0: i64, j0: i64, x: f64, s: f64, t0: f64, span: f64) -> i64, i64, f64, f64, f64 {
        var best_i = -1;
        var best_j = -1;
        var best = 1.0e300;
        var bxi = 0.0;
        var beta = 0.0;
        for dj in 5 {
            for di in 5 {
                let i = i0 + di - 2;
                let j = j0 + dj - 2;
                if i < 0 || j < 0 || i >= v_nx(v) || j >= v_nt(v) || !v_kept(v, i, j) { continue; }
                let xi, eta, gap = cell_coords(v, i, j, x, s, t0, span);
                if gap < best {
                    best = gap;
                    best_i = i;
                    best_j = j;
                    bxi = xi;
                    beta = eta;
                }
            }
        }
        return best_i, best_j, bxi, beta, best;
    }

    // The θ of p from t0, in [0, 2π), or a little below 0 near t0.
    fn arc_of(p: f64[3], t0: f64, span: f64) -> f64 {
        var d = atan2(p[1], p[2]) - t0;
        d -= 2.0 * PI * floor(d / (2.0 * PI));
        if d > span && 2.0 * PI - d < d - span { d -= 2.0 * PI; }
        return d;
    }

    // The shape functions of order o (1 or 2) on [-1, 1] at u, and their
    // derivatives.
    fn shape1(o: i64, u: f64) -> f64[3], f64[3] {
        if o == 1 { return [0.5 * (1.0 - u), 0.5 * (1.0 + u), 0.0], [-0.5, 0.5, 0.0]; }
        return quad_weights(u), [u - 0.5, -2.0 * u, u + 0.5];
    }

    // The (xi, eta) of the point (x, s) in cell (i, j) of view v, from the
    // nodes of the cell, and the distance in x and s from the cell (0 in
    // it). s is r times the θ from t0.
    fn cell_coords(v: i64, i: i64, j: i64, x: f64, s: f64, t0: f64, span: f64) -> f64, f64, f64 {
        let o = v_order(v);
        let la = o * v_nx(v) + 1;
        let r = view(v)[3];
        var px: f64[9] = [0.0; 9];
        var ps: f64[9] = [0.0; 9];
        var c = 0;
        for b in o + 1 {
            for a in o + 1 {
                let q = v_point(v, o * i + a + la * (o * j + b));
                px[c] = q[0];
                ps[c] = r * arc_of(q, t0, span);
                c++;
            }
        }
        var xi = 0.0;
        var eta = 0.0;
        for it in 30 {
            let fx, fs, j = cell_map(o, px, ps, xi, eta);
            let det = j[0] * j[3] - j[1] * j[2];
            if abs(det) <= 1.0e-300 { break; }
            let sx = (j[3] * (fx - x) - j[1] * (fs - s)) / det;
            let se = (j[0] * (fs - s) - j[2] * (fx - x)) / det;
            xi = clamp(xi - sx, -3.0, 3.0);
            eta = clamp(eta - se, -3.0, 3.0);
            if abs(sx) + abs(se) <= 1.0e-12 { break; }
        }
        // The distance from the nearest point of the cell.
        let cx = clamp(xi, -1.0, 1.0);
        let ce = clamp(eta, -1.0, 1.0);
        let qx, qs, j = cell_map(o, px, ps, cx, ce);
        return cx, ce, sqrt((qx - x) * (qx - x) + (qs - s) * (qs - s));
    }

    // The (x, s) of a cell of order o with nodes (px, ps) at (xi, eta), and
    // its Jacobian [dx/dxi, dx/deta, ds/dxi, ds/deta].
    fn cell_map(o: i64, px: f64[9], ps: f64[9], xi: f64, eta: f64) -> f64, f64, f64[4] {
        let lx, dx = shape1(o, xi);
        let lt, dt = shape1(o, eta);
        var x = 0.0;
        var s = 0.0;
        var j: f64[4] = [0.0; 4];
        var c = 0;
        for b in o + 1 {
            for a in o + 1 {
                x += lx[a] * lt[b] * px[c];
                s += lx[a] * lt[b] * ps[c];
                j[0] += dx[a] * lt[b] * px[c];
                j[1] += lx[a] * dt[b] * px[c];
                j[2] += dx[a] * lt[b] * ps[c];
                j[3] += lx[a] * dt[b] * ps[c];
                c++;
            }
        }
        return x, s, j;
    }

    // Field f at point p: the shape functions of its view.
    fn field_at(f: i64, p: f64[3], label: u8[:]) -> f64[9] {
        let v = field_view(f);
        let nodes, w, count = locate(v, p, POINT_TOL, label);
        var r: f64[9] = [0.0; 9];
        for c in count {
            let a = field(f, v, nodes[c]);
            for e in 9 { r[e] += w[c] * a[e]; }
        }
        return r;
    }

    // The largest (big) or the smallest value of scalar field f at the view
    // nodes `nodes`.
    fn node_extreme(f: i64, nodes: i64[:]) -> f64 {
        var best = 0.0;
        for k, i in nodes {
            let x = field(f, 0, k)[0];
            if i == 0 || x > best { best = x; }
        }
        return best;
    }

    // The view nodes of the shell cells of the grid.
    fn surface_nodes() -> i64[>..] {
        let g = shells::grids[grid];
        let la = g.order * g.nx + 1;
        var mark: u8[>..] = [];
        for k in la * (g.order * g.nt + 1) { mark.push(0); }
        for j in g.nt {
            for i in g.nx {
                if shells::grid_cell(grid, i, j) < 0 { continue; }
                for b in g.order + 1 { for a in g.order + 1 { mark[g.order * i + a + la * (g.order * j + b)] = 1; } }
            }
        }
        var out: i64[>..] = [];
        for m, k in mark { if m != 0 { out.push(k); } }
        return out;
    }

    // The sum of the reactions at the view nodes `nodes`.
    fn reaction_of(nodes: i64[:]) -> f64[3] {
        let r = shells::reactions();
        var ns: i64[>..] = [];
        for k in nodes { ns.push(to_node(k)); }
        let x, y, z = shells::reaction_sum(r, ns);
        return [x, y, z];
    }

    // The largest (big) or the smallest axial force of the beams of the runs.
    fn beam_extreme(runs: i64[:]) -> f64 {
        var best = 0.0;
        var seen = false;
        for r in runs {
            for b in run_first[r]..run_first[r] + run_count[r] {
                let x = shells::beam_axial_force(b);
                if !seen || x > best { best = x; }
                seen = true;
            }
        }
        return best;
    }

    // The node at (x, y, z), for a point of a model of beams.
    fn node_at(x: f64, y: f64, z: f64, tol: f64, expr: u8[:]) -> i64 {
        let n, d = shells::nearest_node(x, y, z);
        if d > tol { fault(expr, "is not a node of the mesh"); }
        return n;
    }

    // The axial force of the beam with its middle nearest to (x, y, z).
    fn axial_at(x: f64, y: f64, z: f64, label: u8[:]) -> f64 {
        var best = -1;
        var dist = 0.0;
        for b, i in shells::beams {
            let p = shells::node_point(b.n[0]);
            let e = sub3(shells::node_point(b.n[1]), p);
            let u = sub3([x, y, z], p);
            let t = max(0.0, min(1.0, dot3(u, e) / dot3(e, e)));
            let d = norm3([u[0] - t * e[0], u[1] - t * e[1], u[2] - t * e[2]]);
            if best < 0 || d < dist {
                best = i;
                dist = d;
            }
        }
        if best < 0 || dist > POINT_TOL {
            fault(label, "is not on a beam");
            return 0.0;
        }
        return shells::beam_axial_force(best);
    }

    // --- the mesh, the supports and the loads --------------------------------------

    // Moves lattice point k to (x, θ) on the cylinder of the grid.
    fn move_node(k: i64, x: f64, t: f64) {
        let n = lnode(k);
        let r = shells::grids[grid].r;
        shells::coords[3 * n] = x;
        shells::coords[3 * n + 1] = r * sin(t);
        shells::coords[3 * n + 2] = r * cos(t);
        shells::set_director(n, 0.0, sin(t), cos(t));
    }

    // A beam from node a to node b. With `y_radial`, y is the normal at the
    // middle; with `o_radial`, the offset at each node is o[0] times the
    // normal there.
    fn beam(a: i64, b: i64, section: i64, y: f64[3], y_radial: bool, o: f64[3], o_radial: bool) {
        let pa = shells::node_point(a);
        let pb = shells::node_point(b);
        var ya = y;
        if y_radial { ya = radial([0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]); }
        var oa = o;
        var ob = o;
        if o_radial {
            oa = scale3(o[0], radial(pa));
            ob = scale3(o[0], radial(pb));
        }
        shells::add_beam(a, b, section, ya, oa, ob);
    }

    // A run of beams on the lattice points from (a0, b0) to (a1, b1), a line
    // of x or of θ.
    fn beam_run(a0: i64, b0: i64, a1: i64, b1: i64, section: i64, y: f64[3], y_radial: bool, o: f64[3], o_radial: bool) {
        run_first.push(shells::beams.len);
        let along_x = b0 == b1;
        for k in (if along_x { a1 - a0 } else { b1 - b0 }) {
            let a = if along_x { lnode(lk(a0 + k, b0)) } else { lnode(lk(a0, b0 + k)) };
            let b = if along_x { lnode(lk(a0 + k + 1, b0)) } else { lnode(lk(a0, b0 + k + 1)) };
            beam(a, b, section, y, y_radial, o, o_radial);
        }
        run_count.push(shells::beams.len - run_first[run_first.len - 1]);
    }

    // The node at (x, y, z) if one is within tol, else a new node.
    fn node_or_new(p: f64[3], tol: f64) -> i64 {
        if shells::node_count() > 0 {
            let n, d = shells::nearest_node(p[0], p[1], p[2]);
            if d <= tol { return n; }
        }
        return shells::add_node(p[0], p[1], p[2]);
    }

    // A run of `count` beams on the line from p to q. An end at a node
    // shares it.
    fn line_beams(p: f64[3], q: f64[3], count: i64, tol: f64, section: i64, y: f64[3], o: f64[3]) {
        run_first.push(shells::beams.len);
        var a = node_or_new(p, tol);
        for k in count {
            let f = 1.0 * (k + 1) / count;
            let at: f64[3] = [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])];
            let b = if k == count - 1 { node_or_new(q, tol) } else { shells::add_node(at[0], at[1], at[2]) };
            shells::add_beam(a, b, section, y, o, o);
            a = b;
        }
        run_count.push(count);
    }

    // A force per length on the edge of lattice points a, m (the middle, or
    // -1) and b.
    fn edge_force(a: i64, m: i64, b: i64, fx: f64, fy: f64, fz: f64) {
        if m < 0 {
            shells::edge_load(lnode(a), lnode(b), fx, fy, fz);
        } else {
            shells::edge_load3(lnode(a), lnode(m), lnode(b), fx, fy, fz);
        }
    }

    fn fix_nodes(nodes: i64[:], dof: i64) {
        for k in nodes { shells::fix(to_node(k), dof); }
    }

    // Dofs dof..dof + count of the view nodes `nodes` take the values of
    // field f at their points.
    fn prescribe(f: i64, nodes: i64[:], dof: i64, count: i64, label: u8[:]) {
        for k in nodes {
            let n = to_node(k);
            let v = field_at(f, shells::node_point(n), label);
            for j in count {
                shells::fix(n, dof + j);
                pres_dofs.push(6 * n + dof + j);
                pres_vals.push(v[if count == 1 { dof % 3 } else { j }]);
            }
        }
    }

    // K u = f with the prescribed values: u = u_p + d, with K d = f - K u_p
    // on the free dofs and d = 0 on the fixed dofs.
    fn solve_prescribed() -> bool {
        var f: f64[>..] = [];
        f.append(shells::load);
        for &x in shells::disp { x = 0.0; }
        for d, i in pres_dofs { shells::disp[d] = pres_vals[i]; }
        let r = shells::out_of_balance();
        var up: f64[>..] = [];
        up.append(shells::disp);
        for d in shells::dof_count() { shells::load[d] = if shells::fixed[d] != 0 { 0.0 } else { -r[d] }; }
        let converged = shells::solve();
        for d in shells::dof_count() {
            shells::disp[d] = if shells::fixed[d] != 0 { up[d] } else { shells::disp[d] + up[d] };
            shells::load[d] = f[d];
        }
        return converged;
    }

    // The direct solve has no iterations: cg_iterations is 0.
    fn solve(step: u8[:]) {
        let converged = if pres_dofs.len > 0 { solve_prescribed() } else { shells::solve() };
        format(solves, "  ", step, ",", shells::shells.len + shells::beams.len, ",", shells::dof_count(), ",0,",
               shells::last_seconds, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The lowest `count` eigenpairs of a modal or a buckling step.
    fn eigen(step: u8[:], buckling: bool, count: i64) {
        let t = clock();
        let converged = if buckling { shells::buckling(count) } else { shells::modes(count) };
        if !converged { write_stderr(str("skein run: step ", step, ": ", shells::eigen_error, "\n")); }
        format(solves, "  ", step, ",", shells::shells.len + shells::beams.len, ",", shells::dof_count(), ",0,",
               clock() - t, ",", if converged { "yes" } else { "no" }, "\n");
        solve_count++;
        if !converged { failed = true; }
    }

    // The sum of the reactions and the applied loads is 0 (checks.md,
    // `equilibrium`), to 1e-6 of the sum of the sizes of the loads and the
    // reactions.
    fn equilibrium(unit: u8[:]) {
        let r = shells::reactions();
        var nodes: i64[>..] = [];
        for n in shells::node_count() { nodes.push(n); }
        let rx, ry, rz = shells::reaction_sum(r, nodes);
        let fx, fy, fz = shells::load_sum();
        let e = sqrt((rx + fx) * (rx + fx) + (ry + fy) * (ry + fy) + (rz + fz) * (rz + fz));
        let limit = 1.0e-6 * shells::load_scale(r);
        check("equilibrium", e <= limit, quantity(e, unit), quantity(limit, unit));
    }

    fn expect(name: u8[:], got: f64, want: f64, within: f64, relative: bool, unit: u8[:], limit: u8[:]) {
        let tol = if relative { within * abs(want) } else { within };
        check("expect", abs(got - want) <= tol, str(name, ": ", quantity(got, unit)), limit);
    }

    // Appends the point array of field f (of view 0) to xml.
    fn point_field(xml: u8[>..]&, name: u8[:], f: i64, rank: i64) {
        format(xml, "<DataArray type=\"Float64\" Name=\"", name, "\" NumberOfComponents=\"", rank, "\" format=\"ascii\">\n");
        for n in shells::node_count() {
            let v = if node_view[n] >= 0 { field(f, 0, node_view[n]) } else { pack1(0.0) };
            for c in rank { format(xml, if c > 0 { " " } else { "" }, v[c]); }
            xml.push('\n');
        }
        xml.append("</DataArray>\n");
    }

    // Appends the cell array of the axial force of the beams to xml: 0 on
    // the quads of the shells.
    fn axial_cells(xml: u8[>..]&, name: u8[:]) {
        format(xml, "<DataArray type=\"Float64\" Name=\"", name, "\" format=\"ascii\">\n");
        for c in shells::shells { for q in (if c.count == 9 { 4 } else { 1 }) { xml.append("0\n"); } }
        for b in shells::beams.len { format(xml, shells::beam_axial_force(b), "\n"); }
        xml.append("</DataArray>\n");
    }

    fn output(path: u8[:], w: u8[:], theta: u8[:], points: u8[:], cells: u8[:]) {
        if !shells::write_vtu(path, w, theta, points, cells) {
            write_stderr(str("skein run: cannot write '", path, "'\n"));
            failed = true;
            return;
        }
        format(results, "  ", path, "\n");
        result_count++;
    }

    fn print_records() {
        var s: u8[>..] = [];
        format(s, "report[", report_count, "]{step,expr,value}:\n", reports);
        format(s, "verify[", check_count, "]{stage,id,result,value,limit}:\n", checks);
        format(s, "solve[", solve_count, "]{step,cells,dofs,cg_iterations,seconds,converged}:\n", solves);
        format(s, "results[", result_count, "]:\n", results);
        s.append(solvers::memory_records());
        write_stdout(s);
    }

    """;

let DOF_NAMES = ["shells::UX", "shells::UY", "shells::UZ", "shells::RX", "shells::RY", "shells::RZ"];

// --- the regions -------------------------------------------------------------------

// The view nodes of region r, as code: lattice points, or the node at a
// point of a model of beams.
fn region_code(r: i64) -> u8[>..] {
    let reg = regions[r];
    var out: u8[>..] = [];
    if reg.point && physics_tok < 0 {
        out.append("[node_at(");
        lits(out, reg.at);
        out.append(", ");
        lit(out, 1.0e-3 * mesh_size);
        out.append(", ");
        quoted(out, reg.name);
        out.append(")]");
        return out;
    }
    i64s(out, region_lattice(r));
    return out;
}

// The lattice points of region r: the point nearest to a Point, or the
// points of the edges on a curve of the shell.
fn region_lattice(r: i64) -> i64[>..] {
    let reg = regions[r];
    var out: i64[>..] = [];
    if reg.point {
        var best = -1;
        var d = 0.0;
        for k in px_lat.len {
            let p = lat_point(k);
            let e = sqrt(sq(p[0] - reg.at[0]) + sq(p[1] - reg.at[1]) + sq(p[2] - reg.at[2]));
            if best < 0 || e < d {
                best = k;
                d = e;
            }
        }
        if best < 0 || d > 1.0e-3 * mesh_size {
            unsupported(shell_tok, str("the point ", reg.name, " is not a node of the mesh"), "put the point on a line of the grid, such as a side of the shell");
            return out;
        }
        out.push(best);
        return out;
    }
    if reg.node < 0 { return out; }
    var mark: u8[>..] = [];
    for k in px_lat.len { mark.push(0); }
    for e in edges {
        if !edge_on(reg.node, e) { continue; }
        mark[e.a] = 1;
        mark[e.b] = 1;
        if e.m >= 0 { mark[e.m] = 1; }
    }
    for m, k in mark { if m != 0 { out.push(k); } }
    if out.len == 0 { unsupported(shell_tok, str("the curve ", reg.name, " has no edge of the mesh"), "check that the curve is on the boundary of the meshed surface"); }
    return out;
}

// --- the beams on the grid ---------------------------------------------------------

// A run of beams on lattice line points (a0, b0) to (a1, b1), of piece p
// with beam material mat.
struct Run {
    a0: i64,
    b0: i64,
    a1: i64,
    b1: i64,
    piece: i64,
    mat: i64,
}

var runs: Run[>..];

// Whether lattice segment k0 to k1 is on piece p: its middle is.
fn segment_on(p: Piece, k0: i64, k1: i64) -> bool {
    let tol = 1.0e-6 * mesh_size;
    if p.ring {
        let t = 0.5 * (pt_lat[k0] + pt_lat[k1]);
        return t >= min(p.ta, p.tb) - 1.0e-9 && t <= max(p.ta, p.tb) + 1.0e-9;
    }
    let x = 0.5 * (px_lat[k0] + px_lat[k1]);
    return x >= min(p.p0[0], p.p1[0]) - tol && x <= max(p.p0[0], p.p1[0]) + tol;
}

// The runs of beams of each piece: the lattice segments of its line that
// are on the piece and on a side of a shell cell.
fn make_runs() {
    runs.clear();
    let ps = beam_pieces();
    let tol = geo_tol();
    for p, pi in ps {
        let mat = piece_mat(p);
        if mat < 0 {
            unsupported(beam_tok, "a beam has no material", "give each beam a material, or a material with no on");
            return;
        }
        if physics_tok < 0 {
            runs.push(Run { a0: 0, b0: 0, a1: 0, b1: 0, piece: pi, mat: mat });
            continue;
        }
        // The line of the piece: a θ line, or an x line for a ring.
        var line = -1;
        let lines = if p.ring { xs[..] } else { ts[..] };
        let c = if p.ring { p.p0[0] } else { theta_of(p.p0) };
        for v, i in lines { if abs(v - c) <= (if p.ring { tol } else { tol / radius }) { line = i; } }
        if line < 0 {
            unsupported(beam_tok, "a beam is not on a line of the mesh of the shell", "put each beam on a ring x = c or a line theta = c");
            return;
        }
        let first_run = runs.len;
        let count = if p.ring { lat_nt() } else { lat_nx() };
        var start = -1;
        for s in count + 1 {
            var on = s < count;
            if on {
                let k0 = if p.ring { lat(order * line, s) } else { lat(s, order * line) };
                let k1 = if p.ring { lat(order * line, s + 1) } else { lat(s + 1, order * line) };
                let cell = s / order;
                let side = if p.ring { kept(line - 1, cell) || kept(line, cell) } else { kept(cell, line - 1) || kept(cell, line) };
                on = segment_on(p, k0, k1);
                if on && !side {
                    unsupported(beam_tok, "a beam crosses a removed cell of the shell", "cut the beam curve with the same Solid, such as - door_box(0 m)");
                    return;
                }
            }
            if on && start < 0 { start = s; }
            if !on && start >= 0 {
                if p.ring {
                    runs.push(Run { a0: order * line, b0: start, a1: order * line, b1: s, piece: pi, mat: mat });
                } else {
                    runs.push(Run { a0: start, b0: order * line, a1: s, b1: order * line, piece: pi, mat: mat });
                }
                start = -1;
            }
        }
        if runs.len == first_run {
            unsupported(beam_tok, "a beam is not on a kept cell of the shell", "keep the beam off the cuts of the shell");
            return;
        }
    }
}

// The runs of the pieces of region r, as code.
fn region_runs_code(r: i64) -> u8[>..] {
    if runs.len == 0 { make_runs(); }
    var ids: i64[>..] = [];
    let ps = beam_pieces();
    for u, i in runs { if piece_in(r, ps[u.piece]) { ids.push(i); } }
    var out: u8[>..] = [];
    i64s(out, ids);
    return out;
}

// --- the program -------------------------------------------------------------------

// Writes the Goose program of the model that `read` read. The program
// imports the shells from `shells_path` and writes its VTU files and its
// state under `results_dir`.
fn generate(shells_path: u8[:], results_dir: u8[:], out: u8[>..]&) {
    let sys = codegen::system;
    fields.clear();
    field_ranks.clear();
    runs.clear();
    if beam_tok >= 0 { make_runs(); }
    format(out, "// Generated by skein run from the model ", codegen::model_name, ". Do not edit it: edit the model.\n");
    format(out, "// All values are in the units of the model: length in ", sys.length, ", force in ", sys.force,
           ", angle in rad.\n");
    out.append("import std;\nimport os;\nimport math;\nimport shells from ");
    quoted(out, shells_path);
    out.append(";\nimport solvers from ");
    quoted(out, str(shells_path, "/../solvers"));
    out.append(";   // the shells import it\n");
    out.append(PRELUDE);
    out.append("let POINT_TOL = ");
    lit(out, 1.0e-3 * mesh_size);
    out.append(";\n\n");
    var body: u8[>..] = [];
    for si in 0..steps.len { emit_step(body, si, results_dir); }
    body.append("fn main() {\n    shells::clear();\n");
    if source_name.len > 0 { emit_source(body); }
    if physics_tok >= 0 { emit_shell(body); }
    if beam_tok >= 0 { emit_beams(body); }
    for s in steps { format(body, "    step_", s.name, "();\n"); }
    for x in expects { if x.step < 0 { emit_expect(body, x); } }
    body.append("    print_records();\n");
    body.append("    if failed { exit(1); }\n");
    body.append("}\n");
    // The fields that the steps use.
    for k in fields.len { emit_field_fn(out, k); }
    out.append("fn field(f: i64, v: i64, k: i64) -> f64[9] {\n");
    for k in fields.len { format(out, "    if f == ", k, " { return field_", k, "(v, k); }\n"); }
    out.append("    return pack1(0.0);\n}\n\nfn field_view(f: i64) -> i64 {\n");
    for k in fields.len { if source_member(fields[k]) { format(out, "    if f == ", k, " { return 1; }\n"); } }
    out.append("    return 0;\n}\n\n");
    out.append(body);
}

// Reads the state of the results value, and runs its model first if the
// state is not there. The run of the other model is locked, so the cases
// of a study that need the same state run it once.
fn emit_source(out: u8[>..]&) {
    var sets: const u8[:][>..] = [];
    let dir = codegen::results_path(source_model, sets);
    var inner: u8[>..] = [];
    format(inner, "test -e ", codegen::shell_quoted(str(dir, "/", source_step, ".shell_view")), " || exec ",
           codegen::shell_quoted(str(resource_dir(), "skein")), " run ", codegen::shell_quoted(source_path), " 1>&2");
    var command: u8[>..] = [];
    format(command, "mkdir -p ", codegen::shell_quoted(dir), " && flock ", codegen::shell_quoted(str(dir, "/.lock")), " sh -c ",
           codegen::shell_quoted(inner));
    format(out, "    // ", source_name, ": the state of step ", source_step, " of ", source_path, "\n    results(");
    quoted(out, dir);
    out.append(", ");
    quoted(out, source_step);
    out.append(", ");
    quoted(out, source_path);
    out.append(", ");
    quoted(out, command);
    out.append(");\n");
}

// The laminates, the grid and the moves of its nodes onto the cuts.
fn emit_shell(out: u8[>..]&) {
    let sys = codegen::system;
    for l, k in lams {
        format(out, "    // material ", material_names[shell_mats[k]], ": ", l.count, " plies from the top face\n");
        out.append("    shells::laminate_clear();\n");
        for a in angles[l.first..l.first + l.count] {
            out.append("    shells::laminate_ply(");
            lit(out, a);
            for x in [l.t_ply, l.e1, l.e2, l.g12, l.nu12, l.g12, l.g12] {
                out.append(", ");
                lit(out, x);
            }
            out.append(");\n");
        }
        out.append("    secs.push(shells::laminate_section());\n");
    }
    // The section of each cell, as runs of equal cells.
    var rle: i64[>..] = [];
    for m, i in cell_mat {
        if i > 0 && m == cell_mat[i - 1] {
            rle[rle.len - 2]++;
        } else {
            rle.push(1);
            rle.push(m);
        }
    }
    format(out, "    // mesh on ", physics_name, ": ", nx_cells(), " x ", nt_cells(), if order == 2 { " MITC9" } else { " MITC4" },
           " cells of at most ", mesh_size, " ", sys.length, "\n    {\n        let cells: i64[>..] = ");
    i64s(out, rle);
    out.append(";\n        var cs: i64[>..] = [];\n");
    out.append("        for i in cells.len / 2 { for k in cells[2 * i] { cs.push(if cells[2 * i + 1] >= 0 { secs[cells[2 * i + 1]] } else { -1 }); } }\n");
    out.append("        let xs: f64[>..] = ");
    f64s(out, xs);
    out.append(";\n        let ts: f64[>..] = ");
    f64s(out, ts);
    out.append(";\n        grid = shells::mesh_cylinder_lines(");
    lit(out, radius);
    format(out, ", xs, ts, ", order, ", cs, ");
    lits(out, axis);
    out.append(");\n    }\n");
    // The nodes that the cuts and the material boundaries moved.
    var mk: i64[>..] = [];
    var mx: f64[>..] = [];
    var mt: f64[>..] = [];
    for b in lat_nt() + 1 {
        for a in lat_nx() + 1 {
            let k = lat(a, b);
            let x, t = default_place(a, b);
            if px_lat[k] == x && pt_lat[k] == t { continue; }
            mk.push(k);
            mx.push(px_lat[k]);
            mt.push(pt_lat[k]);
        }
    }
    if mk.len > 0 {
        format(out, "    // ", mk.len, " nodes on the cuts\n    {\n        let ks: i64[>..] = ");
        i64s(out, mk);
        out.append(";\n        let xs: f64[>..] = ");
        f64s(out, mx);
        out.append(";\n        let ts: f64[>..] = ");
        f64s(out, mt);
        out.append(";\n        for k, i in ks { move_node(k, xs[i], ts[i]); }\n    }\n");
    }
}

// The beam sections and the beams of each run: on the lattice lines of the
// shell, or on lines of their own.
fn emit_beams(out: u8[>..]&) {
    for b, k in bmats {
        format(out, "    // material ", material_names[beam_mats[k]], ": a beam section\n    let beam_section_", k, " = shells::beam_section(");
        for v, i in b.section {
            if i > 0 { out.append(", "); }
            lit(out, v);
        }
        out.append(");\n");
        if b.rho > 0.0 {
            format(out, "    shells::beam_density(beam_section_", k, ", ");
            lit(out, b.rho);
            out.append(");\n");
        }
    }
    let ps = beam_pieces();
    format(out, "    // physics ", beam_name, " on ", regions[beam_region].name, ": ", runs.len, " runs of beams\n");
    for u in runs {
        let b = bmats[u.mat];
        if physics_tok < 0 {
            let p = ps[u.piece];
            out.append("    line_beams([");
            lits(out, p.p0);
            out.append("], [");
            lits(out, p.p1);
            format(out, "], ", line_cells(p), ", ");
            lit(out, 1.0e-3 * mesh_size);
            format(out, ", beam_section_", u.mat, ", [");
            lits(out, b.y);
            out.append("], [");
            lits(out, b.off);
            out.append("]);\n");
            continue;
        }
        format(out, "    beam_run(", u.a0, ", ", u.b0, ", ", u.a1, ", ", u.b1, ", beam_section_", u.mat, ", [");
        lits(out, b.y);
        format(out, "], ", if b.y_radial { "true" } else { "false" }, ", [");
        lits(out, b.off);
        format(out, "], ", if b.off_radial { "true" } else { "false" }, ");\n");
    }
}

fn emit_bc(out: u8[>..]&, b: Bc) {
    let r = regions[b.region];
    format(out, "    // bc ", r.name, "\n");
    if b.kind == BC_FIX && b.value < 0 {
        let nodes = region_code(b.region);
        for j in b.count { format(out, "    fix_nodes(", nodes, ", ", DOF_NAMES[b.dof + j], ");\n"); }
    } else if b.kind == BC_FIX {
        format(out, "    prescribe(", field_index(b.value), ", ", region_code(b.region), ", ", b.dof, ", ", b.count, ", ");
        quoted(out, source_of(b.value));
        out.append(");\n");
    } else if b.kind == BC_FORCE {
        let load = vec3(value(b.value));
        format(out, "    {\n        let n = to_node(", region_code(b.region), "[0]);\n");
        for j in 3 {
            if load[j] == 0.0 { continue; }
            format(out, "        shells::force(n, ", DOF_NAMES[b.dof + j], ", ");
            lit(out, load[j]);
            out.append(");\n");
        }
        out.append("    }\n");
    } else {
        // The load of each edge at its middle, with n the outward normal of
        // the edge on the surface.
        edge_mode = true;
        var count = 0;
        for e in edges {
            if !edge_on(r.node, e) { continue; }
            count++;
            edge_n = edge_normal(e);
            let f = value_at(b.value, edge_mid(e));
            if !ok() { break; }
            if f.n != 3 {
                unsupported(tok(b.value), "the force per length is a vector", BC_FIX_TEXT);
                break;
            }
            format(out, "    edge_force(", e.a, ", ", e.m, ", ", e.b, ", ");
            lits(out, vec3(f));
            out.append(");\n");
        }
        edge_mode = false;
        if count == 0 { unsupported(tok(b.value), str("the curve ", r.name, " has no edge on the boundary of the mesh"), "put the load on a side or a cut of the shell"); }
    }
}

fn emit_step(out: u8[>..]&, si: i64, results_dir: u8[:]) {
    let s = steps[si];
    eigen_modal = s.kind == STEP_MODAL;
    format(out, "// step ", s.name, ": ", if s.kind == STEP_STATIC { "static" } else if s.kind == STEP_MODAL { "modal" } else { "buckling" }, "\n");
    format(out, "fn step_", s.name, "() {\n");
    if s.kind != STEP_BUCKLING {
        out.append("    shells::clear_supports();\n    shells::clear_loads();\n    pres_dofs.clear();\n    pres_vals.clear();\n");
        var own: i64[>..] = [];
        for i in global_bcs { own.push(i); }
        for i in s.first_own..s.first_own + s.own_count { own.push(i); }
        for i in own { if bcs[i].kind == BC_FIX { emit_bc(out, bcs[i]); } }
        if s.kind == STEP_STATIC {
            if q_normal || q[0] != 0.0 || q[1] != 0.0 || q[2] != 0.0 {
                format(out, "    // physics ", physics_name, ": q\n    for c in shells::shells.len { shells::shell_load(c, ");
                lits(out, q);
                format(out, ", ", if q_normal { "true" } else { "false" }, "); }\n");
            }
            for i in own { if bcs[i].kind != BC_FIX { emit_bc(out, bcs[i]); } }
        }
    }
    if s.kind == STEP_STATIC {
        out.append("    solve(");
        quoted(out, s.name);
        format(out, ");\n    equilibrium(\"", codegen::system.force, "\");\n    make_own();\n    save_view(");
        quoted(out, results_dir);
        out.append(", ");
        quoted(out, s.name);
        out.append(");\n");
    } else {
        out.append("    eigen(");
        quoted(out, s.name);
        format(out, ", ", if s.kind == STEP_BUCKLING { "true" } else { "false" }, ", ", s.count, ");\n");
    }
    for r in reports[s.first_report..s.first_report + s.report_count] { emit_report(out, s.name, r); }
    for x in expects { if x.step == si { emit_expect(out, x); } }
    if s.output_count > 0 {
        if s.kind == STEP_STATIC {
            emit_output(out, s, str(results_dir, "/", s.name, ".vtu"));
        } else {
            format(out, "    for k in ", s.count, " {\n        shells::show_mode(k + 1);\n        make_own();\n");
            emit_output(out, s, str(results_dir, "/", s.name, "_mode\", k + 1, \".vtu"));
            out.append("    }\n");
        }
    }
    out.append("}\n\n");
}

// The output of step s to the file `path` (code in a str call).
fn emit_output(out: u8[>..]&, s: Step, path: u8[:]) {
    out.append("    {\n        var points: u8[>..] = [];\n        var cells: u8[>..] = [];\n");
    for o in outputs[s.first_output..s.first_output + s.output_count] {
        if is_output_word(o) {
            if name(o) == "Nx" {
                out.append("        axial_cells(cells, \"Nx\");\n");
            }
            continue;
        }
        let k = field_index(o);
        out.append("        point_field(points, ");
        quoted(out, source_of(o));
        format(out, ", ", k, ", ", if field_ranks[k] == 1 { 1 } else if field_ranks[k] == 3 { 3 } else { 9 }, ");\n");
    }
    format(out, "        output(str(\"", path, "\"), ");
    quoted(out, if physics_tok >= 0 { w_name[..] } else { wb_name[..] });
    out.append(", ");
    quoted(out, if physics_tok >= 0 { theta_name[..] } else { thb_name[..] });
    out.append(", points, cells);\n    }\n");
}

fn emit_report(out: u8[>..]&, step: u8[:], r: Report) {
    format(out, "    // report ", r.label, "\n    {\n        let v = ");
    let rank, d = emit_global_expr(r.node, out);
    format(out, ";\n        report(");
    quoted(out, step);
    out.append(", ");
    quoted(out, r.label);
    format(out, ", ", if rank == 1 { "quantity" } else { "vector" }, "(v, ");
    quoted(out, unit_text(d));
    out.append("));\n    }\n");
}

fn emit_expect(out: u8[>..]&, x: Expect) {
    format(out, "    // expect ", x.name, "\n    {\n        let got = ");
    let r, d = emit_global_expr(x.node, out);
    out.append(";\n        let want = ");
    let rw, dw = emit_global_expr(x.want, out);
    out.append(";\n        let within = ");
    if x.relative {
        lit(out, scalar_value(x.within));
    } else {
        let ri, di = emit_global_expr(x.within, out);
    }
    out.append(";\n        expect(");
    quoted(out, x.name);
    format(out, ", got, want, within, ", if x.relative { "true" } else { "false" }, ", ");
    quoted(out, unit_text(d));
    out.append(", ");
    quoted(out, x.label);
    out.append(");\n    }\n");
}
