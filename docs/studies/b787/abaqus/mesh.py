#!/usr/bin/env python3
"""Write the mesh includes for global.inp and door_corner.inp.

Usage: python3 mesh.py
Output, next to this file: global_mesh.inp and door_corner_mesh.inp (model
data), global_loads.inp and door_corner_loads.inp (nodal loads for the step).

The skin is a structured grid in the developed (x, s) plane, s = R * theta.
Grid lines lie on each frame, stringer, door edge and pad-up edge. Cells
cut by a door corner radius are clipped to the arc, which gives the only
triangles. Beam elements are the skin element edges on each beam line.
The script also prints the applied load resultants for checks C3 to C5.
"""

import math
import os

R = 2.883
P_U = 97.5e3
L_HALF = 6.096
PITCH_FRAME = 0.6096
FRAMES = [0.3048 + PITCH_FRAME * k for k in range(10)]
STRINGERS_DEG = [2.0 + 4.0 * k for k in range(45)]  # 0 to 180 deg only

DOOR_A = 0.5334                    # half width along x
DOOR_B = 1.829 / 2.0               # half height along s
DOOR_R = 0.178                     # corner radius
DOOR_S0 = R * math.radians(75.0)   # centre
S_UPPER = DOOR_S0 - DOOR_B
S_LOWER = DOOR_S0 + DOOR_B
PAD = 0.150

ORIENT_NODE_BASE = 900000  # beam orientation nodes, not part of the mesh
BEAM_BASE = 100000


def door_distance(x, s):
    """Signed distance to the cutout edge: > 0 in the skin, < 0 in the cutout."""
    qx = abs(x) - (DOOR_A - DOOR_R)
    qs = abs(s - DOOR_S0) - (DOOR_B - DOOR_R)
    out = math.hypot(max(qx, 0.0), max(qs, 0.0))
    return out + min(max(qx, qs), 0.0) - DOOR_R


def corner_centre(x, s):
    """Centre of the corner arc nearest to (x, s), or None off the corner zones."""
    cx = math.copysign(DOOR_A - DOOR_R, x)
    cs = DOOR_S0 + math.copysign(DOOR_B - DOOR_R, s - DOOR_S0)
    if abs(x) > abs(cx) and abs(s - DOOR_S0) > DOOR_B - DOOR_R:
        return cx, cs
    return None


def merge(points, keep, tol=0.02):
    """Sorted breakpoints; a point within tol of a kept point is dropped."""
    out = sorted(set(keep))
    for p in sorted(points):
        if all(abs(p - q) > tol for q in out):
            out.append(p)
    return sorted(out)


def subdivide(breaks, size):
    """Equal spacing near size inside each interval; size(a, b) -> list."""
    pts = [breaks[0]]
    for a, b in zip(breaks, breaks[1:]):
        pts += size(a, b)[1:]
    return pts


def uniform(h):
    def f(a, b):
        n = max(1, round((b - a) / h))
        return [a + (b - a) * i / n for i in range(n + 1)]
    return f


def graded(h_fine, h_coarse, fine_at, growth=1.25):
    """Geometric spacing from h_fine at the end fine_at ('a' or 'b')."""
    def f(a, b):
        length = b - a
        steps, h = [], h_fine
        while sum(steps) < length:
            steps.append(h)
            h = min(h * growth, h_coarse)
        scale = length / sum(steps)
        steps = [d * scale for d in steps]
        if fine_at == "b":
            steps.reverse()
        pts = [a]
        for d in steps:
            pts.append(pts[-1] + d)
        pts[-1] = b
        return pts
    return f


class Mesh:
    def __init__(self, xs, ss, beam_x, beam_s):
        self.xy = {}            # node id -> (x, s)
        self.shells = []        # (eid, [nids])
        self.beams = []         # (eid, n1, n2, orientation node, set name)
        self.beam_x, self.beam_s = beam_x, beam_s
        self._build(xs, ss)

    def _build(self, xs, ss):
        nx = len(xs)
        grid = {}
        for j, s in enumerate(ss):
            for i, x in enumerate(xs):
                nid = j * nx + i + 1
                grid[i, j] = nid
                self.xy[nid] = (x, s)
        self._snap(xs, ss, grid)
        d = {n: door_distance(*p) for n, p in self.xy.items()}
        for n in d:
            if abs(d[n]) < 1e-12:
                d[n] = 0.0
        cut = {}
        next_id = max(self.xy) + 1
        eid = 0
        for j in range(len(ss) - 1):
            for i in range(nx - 1):
                ring = [grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]]
                if all(d[n] >= 0 for n in ring):
                    mid = [sum(self.xy[n][k] for n in ring) / 4 for k in (0, 1)]
                    assert door_distance(*mid) > 0, f"cell {i},{j} hides the cutout"
                    poly = ring
                elif all(d[n] <= 0 for n in ring):
                    continue
                else:
                    poly = []
                    for a, b in zip(ring, ring[1:] + ring[:1]):
                        if d[a] >= 0:
                            poly.append(a)
                        if d[a] * d[b] < 0:
                            key = (min(a, b), max(a, b))
                            if key not in cut:
                                cut[key] = next_id
                                self.xy[next_id] = self._arc_point(a, b)
                                d[next_id] = 0.0
                                next_id += 1
                            poly.append(cut[key])
                    poly = [n for k, n in enumerate(poly) if n != poly[k - 1]]
                    if len(poly) < 3 or self._area(poly) < 1e-4 * self._area(ring):
                        continue
                for part in self._split(poly):
                    eid += 1
                    self.shells.append((eid, part))
        used = {n for _, c in self.shells for n in c}
        self.xy = {n: p for n, p in self.xy.items() if n in used}
        self._beams()

    def _snap(self, xs, ss, grid):
        """Move grid nodes near a corner arc onto the arc, so no cut is a sliver."""
        hx = {i: min(b - a for a, b in zip(xs[max(i - 1, 0):i + 2], xs[max(i - 1, 0) + 1:i + 2])) for i in range(len(xs))}
        hs = {j: min(b - a for a, b in zip(ss[max(j - 1, 0):j + 2], ss[max(j - 1, 0) + 1:j + 2])) for j in range(len(ss))}
        for (i, j), nid in grid.items():
            x, s = self.xy[nid]
            c = corner_centre(x, s)
            if c is None:
                continue
            dist = math.hypot(x - c[0], s - c[1]) - DOOR_R
            if abs(dist) >= 0.3 * min(hx[i], hs[j]):
                continue
            on_s = any(abs(s - b) < 1e-9 for b in self.beam_s)
            on_x = any(abs(x - b) < 1e-9 for b in self.beam_x)
            if on_s:   # keep the node on its stringer: move along x
                x = c[0] + math.copysign(math.sqrt(DOOR_R ** 2 - (s - c[1]) ** 2), x - c[0])
            elif on_x:  # keep the node on its frame: move along s
                s = c[1] + math.copysign(math.sqrt(DOOR_R ** 2 - (x - c[0]) ** 2), s - c[1])
            else:
                k = DOOR_R / (dist + DOOR_R)
                x, s = c[0] + (x - c[0]) * k, c[1] + (s - c[1]) * k
            self.xy[nid] = (x, s)

    def _arc_point(self, a, b):
        (xa, sa), (xb, sb) = self.xy[a], self.xy[b]
        lo, hi = 0.0, 1.0
        da = door_distance(xa, sa)
        for _ in range(60):
            t = 0.5 * (lo + hi)
            if (door_distance(xa + t * (xb - xa), sa + t * (sb - sa)) >= 0) == (da >= 0):
                lo = t
            else:
                hi = t
        return xa + lo * (xb - xa), sa + lo * (sb - sa)

    def _area(self, poly):
        pts = [self.xy[n] for n in poly]
        return 0.5 * sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))

    def _split(self, poly):
        if len(poly) <= 4:
            return [poly]
        assert len(poly) == 5, poly
        best = max(range(5), key=lambda k: min(
            self._area([poly[k], poly[(k + 1) % 5], poly[(k + 2) % 5]]),
            self._area([poly[k], poly[(k + 2) % 5], poly[(k + 3) % 5], poly[(k + 4) % 5]])))
        k = best
        return [[poly[(k + m) % 5] for m in (0, 2, 3, 4)], [poly[k], poly[(k + 1) % 5], poly[(k + 2) % 5]]]

    def _beams(self):
        edges = set()
        for _, c in self.shells:
            for a, b in zip(c, c[1:] + c[:1]):
                edges.add((min(a, b), max(a, b)))
        eid = BEAM_BASE
        for a, b in sorted(edges):
            (xa, sa), (xb, sb) = self.xy[a], self.xy[b]
            if abs(sa - sb) < 1e-9 and any(abs(sa - v) < 1e-9 for v in self.beam_s):
                n1, n2 = (a, b) if xa < xb else (b, a)
                eid += 1
                self.beams.append((eid, n1, n2, ORIENT_NODE_BASE, "STRINGERS"))
            elif abs(xa - xb) < 1e-9 and any(abs(xa - v) < 1e-9 for v in self.beam_x):
                n1, n2 = (a, b) if sa < sb else (b, a)
                k = min(range(len(self.beam_x)), key=lambda m: abs(self.beam_x[m] - xa))
                eid += 1
                self.beams.append((eid, n1, n2, ORIENT_NODE_BASE + 1 + k, "FRAMES"))

    # -- sets ---------------------------------------------------------------

    def nodes_where(self, test):
        return sorted(n for n, (x, s) in self.xy.items() if test(x, s))

    def elements_at(self, x0, s0):
        """Shell elements that touch the point (x0, s0), nearest centroid first."""
        def cen(c):
            return [sum(self.xy[n][k] for n in c) / len(c) for k in (0, 1)]
        ranked = sorted(self.shells, key=lambda e: math.dist(cen(e[1]), (x0, s0)))
        near = [e for e in ranked if all(
            min(p[k] for p in map(self.xy.get, e[1])) - 1e-9 <= v <= max(p[k] for p in map(self.xy.get, e[1])) + 1e-9
            for k, v in ((0, x0), (1, s0)))]
        return [e[0] for e in near] or [ranked[0][0]]

    def edge_loads(self, chain, along, q, direction):
        """Lumped nodal forces of a line load q(x, s), per unit of coordinate
        along (0 = x, 1 = s), on a chain of nodes (trapezoidal rule)."""
        f = {n: 0.0 for n in chain}
        for a, b in zip(chain, chain[1:]):
            half = 0.25 * (q(*self.xy[a]) + q(*self.xy[b])) * (self.xy[b][along] - self.xy[a][along])
            f[a] += half
            f[b] += half
        lines = []
        for n in chain:
            for dof, comp in direction(*self.xy[n]):
                lines.append((n, dof, f[n] * comp))
        return lines

    # -- output -------------------------------------------------------------

    def write(self, path, header, extra):
        out = [f"** {header}", "** Generated by mesh.py. Do not edit by hand.", "*NODE, NSET=NALL"]
        for n in sorted(self.xy):
            x, s = self.xy[n]
            out.append(f"{n}, {x:.6f}, {R * math.sin(s / R):.6f}, {R * math.cos(s / R):.6f}")
        out.append("*NODE, NSET=BEAM_ORIENT")
        out.append(f"{ORIENT_NODE_BASE}, 0., 0., 0.")
        for k, x in enumerate(self.beam_x):
            out.append(f"{ORIENT_NODE_BASE + 1 + k}, {x:.6f}, 0., 0.")
        pad = set(self.pad_elements())
        for etype, nn in (("S4R", 4), ("S3R", 3)):
            for region in ("SKIN", "PADUP"):
                rows = [(e, c) for e, c in self.shells if len(c) == nn and (e in pad) == (region == "PADUP")]
                if rows:
                    out.append(f"*ELEMENT, TYPE={etype}, ELSET={region}")
                    out += [f"{e}, " + ", ".join(map(str, c)) for e, c in rows]
        for name in ("STRINGERS", "FRAMES"):
            rows = [b for b in self.beams if b[4] == name]
            if rows:
                out.append(f"*ELEMENT, TYPE=B31, ELSET={name}")
                out += [f"{e}, {a}, {b}, {o}" for e, a, b, o, _ in rows]
        out.append("*ELSET, ELSET=SHELLS")
        out.append("SKIN, PADUP")
        loads = [f"** {header}, nodal loads (history data, include inside the step)",
                 "** Generated by mesh.py. Do not edit by hand."]
        for kind, name, items in extra:
            if kind == "CLOAD":
                loads += [f"** {name}", "*CLOAD"]
                loads += [f"{n}, {dof}, {v:.4f}" for n, dof, v in items]
            else:
                out.append(f"*{kind}, {kind}={name}")
                for k in range(0, len(items), 16):
                    out.append(", ".join(map(str, items[k:k + 16])))
        for name, text in ((path + "_mesh.inp", out), (path + "_loads.inp", loads)):
            with open(name, "w") as fh:
                fh.write("\n".join(text) + "\n")

    def pad_elements(self):
        return [e for e, c in self.shells
                if door_distance(*[sum(self.xy[n][k] for n in c) / len(c) for k in (0, 1)]) <= PAD]


DOOR_LOAD = ("Door load on the aft half of the cutout edge, p_U * x_e(s) per unit s, "
             "along the radius (%.1f N/m on the straight edge)" % (P_U * DOOR_A))


def door_edge(mesh):
    """Cutout edge nodes in the order of s. The door strip from x = 0 to the
    edge at each s loads that edge point, so the straight upper and lower
    edges (one s each) get no load. Ties sort so that the chain goes from
    the upper edge into the upper arc, and from the lower arc into the
    lower edge, at the node nearest to the arc."""
    nodes = mesh.nodes_where(lambda x, s: abs(door_distance(x, s)) < 1e-6)
    return sorted(nodes, key=lambda n: (mesh.xy[n][1], mesh.xy[n][0] if mesh.xy[n][1] < DOOR_S0 else -mesh.xy[n][0]))


def radial(x, s):
    return [(2, math.sin(s / R)), (3, math.cos(s / R))]


def axial(x, s):
    return [(1, 1.0)]


def resultant(mesh, cloads):
    """Pressure (1-point per element, flat facets) plus nodal loads, in N."""
    f = [0.0, 0.0, 0.0]
    for _, c in mesh.shells:
        p = [(x, R * math.sin(s / R), R * math.cos(s / R)) for x, s in map(mesh.xy.get, c)]
        tris = [(p[0], p[k], p[k + 1]) for k in range(1, len(p) - 1)]
        for a, b, cc in tris:
            u = [b[k] - a[k] for k in range(3)]
            v = [cc[k] - a[k] for k in range(3)]
            n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            for k in range(3):
                f[k] += 0.5 * P_U * n[k]
    for _, _, items in cloads:
        for _, dof, v in items:
            f[dof - 1] += v
    return f


def global_model(here):
    stringer_s = [R * math.radians(t) for t in STRINGERS_DEG]
    xs_break = merge([DOOR_A + PAD], [0.0, DOOR_A, L_HALF] + FRAMES)
    ss_break = merge([S_UPPER - PAD, S_LOWER + PAD], [0.0, S_UPPER, S_LOWER, R * math.pi] + stringer_s)
    xs = subdivide(xs_break, uniform(0.05))
    ss = subdivide(ss_break, uniform(0.05))
    m = Mesh(xs, ss, FRAMES, stringer_s)

    tol = 1e-6
    p_x, p_s = 4.8768, R * math.radians(140.0)
    end = m.nodes_where(lambda x, s: abs(x - L_HALF) < tol)
    door = door_edge(m)
    cloads = [
        ("CLOAD", "End closure load, p_U R / 2 = %.1f N/m along +x" % (P_U * R / 2),
         m.edge_loads(sorted(end, key=lambda n: m.xy[n][1]), 1, lambda x, s: P_U * R / 2, axial)),
        ("CLOAD", DOOR_LOAD, m.edge_loads(door, 1, lambda x, s: P_U * x, radial)),
    ]
    sub_x0, sub_s0 = DOOR_A - DOOR_R - 0.3, S_UPPER + DOOR_R - 0.3
    sets = [
        ("NSET", "SYM_X", m.nodes_where(lambda x, s: abs(x) < tol)),
        ("NSET", "SYM_Y", m.nodes_where(lambda x, s: s < tol or abs(s - R * math.pi) < tol)),
        ("NSET", "RIGID_Z", m.nodes_where(lambda x, s: abs(x) < tol and abs(s - R * math.pi) < tol)),
        ("NSET", "END", end),
        ("NSET", "DOOR_EDGE", door),
        ("NSET", "CHECK_P", m.nodes_where(lambda x, s: abs(x - p_x) < tol and abs(s - p_s) < tol)),
        ("ELSET", "CHECK_P", m.elements_at(p_x, p_s)[:4]),
        ("ELSET", "C6_TOP", m.elements_at(sub_x0 + 0.3, sub_s0)[:2]),
        ("ELSET", "C6_AFT", m.elements_at(sub_x0 + 0.6, sub_s0 + 0.3)[:2]),
    ]
    m.write(os.path.join(here, "global"), "Global model mesh, quarter barrel", sets + cloads)
    report("global", m, cloads)
    assert len(sets[5][2]) == 1, "check point P is not a node"


def submodel(here):
    x0, x1 = DOOR_A - DOOR_R - 0.3, DOOR_A - DOOR_R + 0.3
    s_c = S_UPPER + DOOR_R
    s0, s1 = s_c - 0.3, s_c + 0.3
    beam_x = [x for x in FRAMES if x0 < x < x1]
    beam_s = [R * math.radians(t) for t in STRINGERS_DEG if s0 < R * math.radians(t) < s1]
    xs = subdivide([x0] + beam_x + [DOOR_A], uniform(0.005)) + graded(0.005, 0.02, "a")(DOOR_A, x1)[1:]
    ss = graded(0.005, 0.02, "b")(s0, S_UPPER) + subdivide([S_UPPER] + beam_s + [s1], uniform(0.005))[1:]
    m = Mesh(xs, ss, beam_x, beam_s)

    tol = 1e-6
    door = door_edge(m)
    cloads = [("CLOAD", DOOR_LOAD, m.edge_loads(door, 1, lambda x, s: P_U * x, radial))]
    sets = [
        ("NSET", "DRIVEN", m.nodes_where(lambda x, s: min(abs(x - x0), abs(x - x1), abs(s - s0), abs(s - s1)) < tol)),
        ("NSET", "DOOR_EDGE", door),
        ("NSET", "CUTOUT_EDGE", m.nodes_where(lambda x, s: abs(door_distance(x, s)) < 1e-6)),
        ("ELSET", "C6_TOP", m.elements_at(x0 + 0.3, s0)[:2]),
        ("ELSET", "C6_AFT", m.elements_at(x1, s0 + 0.3)[:2]),
    ]
    edge = set(sets[2][2])
    sets.append(("ELSET", "CUTOUT_EDGE", [e for e, c in m.shells if edge.intersection(c)]))
    m.write(os.path.join(here, "door_corner"), "Door-corner submodel mesh", sets + cloads)
    report("door corner", m, cloads)


def report(name, m, cloads):
    tri = sum(1 for _, c in m.shells if len(c) == 3)
    print(f"{name}: {len(m.xy)} nodes, {len(m.shells)} shells ({tri} S3R), {len(m.beams)} beams, "
          f"{len(m.pad_elements())} pad-up shells")
    f = resultant(m, cloads)
    print(f"  applied load resultant (N): Fx {f[0]:.1f}, Fy {f[1]:.1f}, Fz {f[2]:.1f}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    global_model(here)
    submodel(here)
