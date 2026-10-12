"""Reusable medial-axis reconstruction: ink polygon -> strokes (centreline, width profile, caps, primitives)."""
import math, time
import numpy as np
import shapely
from shapely.geometry import Polygon, MultiPoint, LineString, Point
from shapely.ops import voronoi_diagram, unary_union
from scipy.spatial import cKDTree
from geom import *



def extend_prim(kind, data, ext, at_start):
    """continue a primitive by `ext` units beyond its start or end"""
    if kind == 'line':
        a, b = np.asarray(data[0], float), np.asarray(data[1], float)
        u = (a - b) / (np.linalg.norm(a - b) + 1e-12)
        return ('line', [a + u * ext, b]) if at_start else ('line', [a, b - u * ext])
    if kind == 'arc':
        cx, cy, r, a0, a1 = data
        sgn = 1 if a1 >= a0 else -1
        da = ext / max(r, 1e-6)
        return ('arc', (cx, cy, r, a0 - sgn * da, a1)) if at_start else ('arc', (cx, cy, r, a0, a1 + sgn * da))
    q = primitive_points(kind, data, 0.03)
    t = end_tangent(q, at_start, 0.4)
    if at_start:
        return ('multi', [('line', [q[0] - t * ext, q[0]]), (kind, data)])
    return ('multi', [(kind, data), ('line', [q[-1], q[-1] - t * ext])])


class MAT:
    def __init__(self, ink, step=0.03, log=print, min_residual=0.0):
        self.ink = ink
        self.step = step
        self.log = log
        self.boundary = ink.boundary
        t0 = time.time()
        rings = []
        for poly in getattr(ink, 'geoms', [ink]):
            for ring in [poly.exterior, *poly.interiors]:
                rings.append(resample(np.asarray(ring.coords), step))
        self.rings = rings
        self.islands = [np.asarray(p.exterior.coords) for p in getattr(ink, 'geoms', [ink]) if p.area < 15 and len(p.interiors) == 0]
        B, BR, BP = [], [], []
        for ri, pts in enumerate(rings):
            B.append(pts); BR.append(np.full(len(pts), ri)); BP.append(np.arange(len(pts)))
        self.B = np.vstack(B); self.BR = np.concatenate(BR); self.BP = np.concatenate(BP)
        self.ring_len = {ri: len(p) for ri, p in enumerate(rings)}
        self.tree = cKDTree(self.B)
        # Voronoi -> medial axis edges
        vd = voronoi_diagram(MultiPoint(self.B), edges=True)
        edges = []

        def walk(g):
            if g.geom_type == 'LineString':
                c = np.asarray(g.coords)
                for a, b in zip(c[:-1], c[1:]):
                    edges.append((a, b))
            else:
                for h in getattr(g, 'geoms', []):
                    walk(h)
        walk(vd)
        E = np.array(edges)
        mid = E.mean(axis=1)
        inside = (shapely.contains_xy(ink, mid[:, 0], mid[:, 1]) & shapely.contains_xy(ink, E[:, 0, 0], E[:, 0, 1])
                  & shapely.contains_xy(ink, E[:, 1, 0], E[:, 1, 1]))
        E, mid = E[inside], mid[inside]
        _, idx2 = self.tree.query(mid, k=2)
        same = self.BR[idx2[:, 0]] == self.BR[idx2[:, 1]]
        gap = np.abs(self.BP[idx2[:, 0]] - self.BP[idx2[:, 1]])
        rl = np.array([self.ring_len[r] for r in self.BR[idx2[:, 0]]])
        gap = np.minimum(gap, rl - gap)
        # Ogniewicz residual: boundary arc between the two generators minus their chord.
        # Noise spurs (anti-aliased raster bumps) have residual ~0; real axis branches grow past it.
        chord = np.linalg.norm(self.B[idx2[:, 0]] - self.B[idx2[:, 1]], axis=1)
        resid = gap * step - chord
        E = E[~(same & ((gap <= 3) | (resid < min_residual)))]
        self.adj = {}
        for a, b in E:
            ka, kb = self.key(a), self.key(b)
            if ka == kb:
                continue
            self.adj.setdefault(ka, set()).add(kb)
            self.adj.setdefault(kb, set()).add(ka)
        self.raw_branches = self.extract_branches()
        self.raw_nodes = [np.array(k) for k, nb in self.adj.items() if len(nb) != 2]
        log(f'medial axis: {len(self.B)} boundary samples, {len(E)} edges, {len(self.raw_branches)} branches ({time.time() - t0:.0f}s)')

    @staticmethod
    def key(p):
        return (round(float(p[0]), 6), round(float(p[1]), 6))

    def radius(self, p):
        return float(self.tree.query(np.asarray(p, float))[0])

    def deg(self, k):
        return len(self.adj.get(k, ()))

    def extract_branches(self):
        seen = set(); branches = []
        nodes = [k for k, nb in self.adj.items() if len(nb) != 2]
        for n0 in nodes:
            for nb in list(self.adj[n0]):
                if (n0, nb) in seen:
                    continue
                path = [n0, nb]; seen.add((n0, nb)); seen.add((nb, n0))
                prev, cur = n0, nb
                while len(self.adj[cur]) == 2:
                    nxt = [x for x in self.adj[cur] if x != prev][0]
                    if (cur, nxt) in seen:
                        break
                    seen.add((cur, nxt)); seen.add((nxt, cur))
                    path.append(nxt); prev, cur = cur, nxt
                    if cur == n0:
                        break
                branches.append(np.array(path, float))
        return branches

    def remove_branch(self, pts):
        ks = [self.key(p) for p in pts]
        for a, b in zip(ks[:-1], ks[1:]):
            self.adj.get(a, set()).discard(b); self.adj.get(b, set()).discard(a)
        for k in ks:
            if k in self.adj and not self.adj[k]:
                del self.adj[k]

    def prune(self, spur_len=1.6, end_radius=0.25):
        for _ in range(4):
            changed = False
            for pts in self.extract_branches():
                k0, k1 = self.key(pts[0]), self.key(pts[-1])
                d0, d1 = self.deg(k0), self.deg(k1)
                L = length(pts)
                if L < spur_len and ((d0 == 1 and self.radius(pts[0]) < end_radius and d1 >= 3) or (d1 == 1 and self.radius(pts[-1]) < end_radius and d0 >= 3)):
                    self.remove_branch(pts); changed = True
                elif L < 0.5 and d0 == 1 and d1 == 1:
                    self.remove_branch(pts); changed = True
            if not changed:
                break

    def merge(self, straight_deg=30):
        branches = self.extract_branches()
        ends = {i: (self.key(p[0]), self.key(p[-1])) for i, p in enumerate(branches)}
        bmap = {i: p for i, p in enumerate(branches)}
        through = set()

        def node_branches(k):
            out = []
            for i, (a, b) in ends.items():
                if a == k: out.append((i, 0))
                if b == k: out.append((i, 1))
            return out

        def merge_pass():
            nodes = set()
            for a, b in ends.values():
                nodes.add(a); nodes.add(b)
            for n in sorted(nodes):
                inc = node_branches(n)
                if len(inc) < 2:
                    continue
                best = None
                for i in range(len(inc)):
                    for j in range(i + 1, len(inc)):
                        (b1, e1), (b2, e2) = inc[i], inc[j]
                        if b1 == b2:
                            continue
                        d1 = end_tangent(bmap[b1], e1 == 0, 0.6); d2 = end_tangent(bmap[b2], e2 == 0, 0.6)
                        ang = math.degrees(math.acos(float(np.clip(-np.dot(d1, d2), -1, 1))))
                        lim = 181 if len(inc) == 2 else straight_deg
                        if ang < lim and (best is None or ang < best[0]):
                            best = (ang, b1, e1, b2, e2)
                if best is None:
                    continue
                _, b1, e1, b2, e2 = best
                p1 = bmap[b1] if e1 == 1 else bmap[b1][::-1]
                p2 = bmap[b2] if e2 == 0 else bmap[b2][::-1]
                new_ends = (ends[b1][0] if e1 == 1 else ends[b1][1], ends[b2][1] if e2 == 0 else ends[b2][0])
                del bmap[b2], ends[b2]
                bmap[b1] = np.vstack([p1, p2[1:]]); ends[b1] = new_ends
                if len(inc) > 2:
                    through.add(n)
                return True
            return False

        while merge_pass():
            pass
        self.bmap, self.ends, self.through = bmap, ends, through
        self.log(f'strokes after merging: {len(bmap)}')

    def ray_hit(self, p, dvec, maxu=3.0):
        ls = LineString([p, p + dvec * maxu])
        g = ls.intersection(self.boundary)
        if g.is_empty:
            return maxu
        pts = [np.asarray(q.coords[0]) for q in getattr(g, 'geoms', [g]) if q.geom_type == 'Point']
        return float(min(np.linalg.norm(q - p) for q in pts)) if pts else maxu

    def strokes(self, tol=0.02, exclude=None, width_snap=None, cap_snap_deg=None, smooth=0, corner_deg=35,
                snap_tol=None, min_junction_len=0.0, taper_ratio=0.0, rms=None):
        """build strokes with widths (optionally snapped to a set of weights), caps and primitives.
        smooth: half window (in 0.05 u samples) of a moving average applied to the axis (raster noise)
        snap_tol: only snap a width to a weight when it is this close, else keep it rounded to 0.05
        min_junction_len: drop junction-to-junction fragments shorter than this (covered by their neighbours)
        taper_ratio: a free-ended stroke whose end radius falls below this fraction of its base is a tapered spike"""
        out = []
        for bid, raw in self.bmap.items():
            if exclude is not None and np.mean(shapely.contains_xy(exclude, raw[:, 0], raw[:, 1])) > 0.6:
                continue
            k0, k1 = self.ends[bid]
            free0, free1 = self.deg(k0) == 1, self.deg(k1) == 1
            pts = resample(raw, 0.05)
            if len(pts) < 3:
                continue
            L0 = length(pts)
            if not free0 and not free1 and L0 < min_junction_len:
                continue
            rad = np.array([self.radius(p) for p in pts])
            n = len(pts)
            s_along = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
            if smooth and n > 2 * smooth + 2:
                k = int(smooth)
                pts = np.stack([np.convolve(np.pad(pts[:, i], (k, k), mode='reflect', reflect_type='odd'),
                                            np.ones(2 * k + 1) / (2 * k + 1), mode='valid') for i in (0, 1)], 1)
            win = 20
            prof = np.array([2 * np.percentile(rad[max(0, i - win):i + win + 1], 20) for i in range(n)])
            axis0 = pts.copy()
            # width from the radii away from both ends (junctions inflate, caps shrink)
            lo = 0.7 if free0 else 1.0; hi = 0.7 if free1 else 1.0
            keep = (s_along >= lo) & (s_along <= L0 - hi)
            core = rad[keep]
            w_raw = 2 * float(np.median(core)) if len(core) >= 5 else 2 * float(np.percentile(rad, 40))
            # tapered spike: one free end, the radius falls towards it
            taper = None
            if taper_ratio and (free0 != free1) and L0 < 6:
                q = max(3, n // 4)
                r_end = rad[-q:] if free1 else rad[:q]
                r_base = rad[q:2 * q + 1] if free0 else rad[-2 * q - 1:-q]
                if len(r_base) and np.median(r_end) < taper_ratio * np.median(r_base):
                    wb = 2 * float(np.percentile(r_base, 70))
                    taper = (round(wb, 2), 0.1)
                    w_raw = wb
            w = w_raw
            if width_snap and taper is None:
                cand = min(width_snap, key=lambda s: abs(s - w_raw))
                w = cand if (snap_tol is None or abs(cand - w_raw) <= snap_tol) else round(w_raw / 0.05) * 0.05
            pr_core = prof[int(n * 0.1):max(int(n * 0.9), int(n * 0.1) + 1)]
            var_w = (pr_core.max() - pr_core.min()) > 0.12 if len(pr_core) > 10 else False
            if taper is not None:
                var_w = False
            hw = w / 2

            # junction ends: the axis bends towards the node there, so trim it and later extend the fitted
            # primitive straight on (a line stays a line, an arc stays an arc)
            cut = hw * 1.25
            ext0 = ext1 = 0.0
            short = L0 < 2 * cut + 0.8          # too short to trim: keep the whole axis, just reach 0.15 past each node
            if not free0:
                if short:
                    ext0 = 0.15
                else:
                    d = np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1)); k = int(np.searchsorted(d, cut)) + 1
                    if k < len(pts) - 3:
                        ext0 = float(d[k - 1]) + 0.15; pts = pts[k:]
            if not free1:
                if short:
                    ext1 = 0.15
                else:
                    d = np.cumsum(np.linalg.norm(np.diff(pts[::-1], axis=0), axis=1)); k = int(np.searchsorted(d, cut)) + 1
                    if k < len(pts) - 3:
                        ext1 = float(d[k - 1]) + 0.15; pts = pts[:len(pts) - k]
            if short and (not free0 or not free1) and len(pts) >= 2:
                pts = np.array([pts[0], pts[-1]]) if length(pts) < 1.0 else pts   # a stub between two nodes is a straight connector
            caps = [None, None]
            for end in (0, 1):
                if not (free0 if end == 0 else free1):
                    continue
                p = pts[0] if end == 0 else pts[-1]
                t = -end_tangent(pts, end == 0, 0.5)
                nn = np.array([-t[1], t[0]])
                dc = min(self.ray_hit(p, t), 2.9)
                rl = min(hw, self.radius(p)) - 0.1
                pl, pr = p + nn * rl, p - nn * rl
                dl = self.ray_hit(pl, t) if self.ink.contains(Point(pl)) else dc
                dr = self.ray_hit(pr, t) if self.ink.contains(Point(pr)) else dc
                dl = dc if (dl < 0.02 or dl >= 2.9) else dl
                dr = dc if (dr < 0.02 or dr >= 2.9) else dr
                if taper is not None:
                    dl = dr = dc
                if cap_snap_deg is not None:
                    ang = math.degrees(math.atan2(abs(dl - dr), w))
                    if ang < cap_snap_deg:
                        dl = dr = (dl + dr) / 2
                caps[end] = (float(dl), float(dr))
                ext = max(dl, dr)
                pts = np.vstack([p + t * ext, pts]) if end == 0 else np.vstack([pts, p + t * ext])
            cuts = [0] + corner_indices(pts, ang_deg=corner_deg, min_seg=0.7) + [len(pts) - 1]
            prims = []

            def flat(k, dd):
                if k == 'multi':
                    for k2, d2 in dd:
                        yield from flat(k2, d2)
                else:
                    yield (k, dd)
            for a, b in zip(cuts, cuts[1:]):
                piece = pts[a:b + 1]
                if len(piece) >= 2:
                    prims.extend(flat(*classify(piece, tol, rms=rms)))
            if prims and ext0 > 0:
                prims = list(flat(*extend_prim(*prims[0], ext0, True))) + prims[1:]
            if prims and ext1 > 0:
                prims = prims[:-1] + list(flat(*extend_prim(*prims[-1], ext1, False)))
            out.append(dict(id=int(bid), width=round(w, 3), width_raw=round(w_raw, 3), len=round(length(pts), 2),
                            free=[bool(free0), bool(free1)], caps=caps, prims=prims, var_w=bool(var_w), axis=axis0, prof=prof,
                            taper=taper))
        return out

    def render_clean(self, strokes, fills=(), join='mitre'):
        geoms = []
        for s in strokes:
            pts = []
            for kind, data in s['prims']:
                q = primitive_points(kind, data, 0.02)
                pts.extend(q if not pts else q[1:])
            pts = np.array(pts)
            if len(pts) < 2:
                continue
            hw = s['width'] / 2
            if s.get('taper'):
                wb, wt = s['taper']
                base_first = s['free'][1]            # tip is the free end
                sa = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
                u = sa / max(sa[-1], 1e-9)
                wq = wb + (wt - wb) * (u if base_first else 1 - u)
                tan = np.gradient(pts, axis=0); tan /= (np.linalg.norm(tan, axis=1)[:, None] + 1e-12)
                nor = np.stack([-tan[:, 1], tan[:, 0]], 1)
                g = Polygon(np.vstack([pts + nor * (wq / 2)[:, None], (pts - nor * (wq / 2)[:, None])[::-1]])).buffer(0)
                geoms.append(g)
                continue
            if s['var_w']:
                tr = cKDTree(s['axis']); wq = s['prof'][tr.query(pts)[1]]
                tan = np.gradient(pts, axis=0); tan /= (np.linalg.norm(tan, axis=1)[:, None] + 1e-12)
                nor = np.stack([-tan[:, 1], tan[:, 0]], 1)
                g = Polygon(np.vstack([pts + nor * (wq / 2)[:, None], (pts - nor * (wq / 2)[:, None])[::-1]])).buffer(0)
                hw0, hw1 = float(wq[0] / 2), float(wq[-1] / 2)
            else:
                g = LineString(pts).buffer(hw, cap_style='flat', join_style=join, mitre_limit=6, quad_segs=16)
                hw0 = hw1 = hw
            for end in (0, 1):
                c = s['caps'][end]
                if c is None or abs(c[0] - c[1]) < 0.01:
                    continue
                dl, dr = c
                p = pts[0] if end == 0 else pts[-1]
                t = -end_tangent(pts, end == 0, 0.5); nn = np.array([-t[1], t[0]])
                m = max(dl, dr); h_ = hw0 if end == 0 else hw1
                cl = p - t * (m - dl) + nn * h_; cr = p - t * (m - dr) - nn * h_
                back = -t
                g = g.intersection(Polygon([cl, cr, cr + back * 8, cl + back * 8]).buffer(0).union(
                    Polygon([cl + nn * 3, cr - nn * 3, cr - nn * 3 + back * 8, cl + nn * 3 + back * 8]).buffer(0)))
            geoms.append(g)
        return unary_union(geoms + list(fills))

    def render_exact(self, fills=()):
        """union of discs along every axis branch with exact radii: reproduces the ink"""
        pts = []
        for b in self.raw_branches:
            pts.extend(resample(b, 0.03))
        pts.extend(self.raw_nodes)
        pts = np.array(pts)
        r = self.tree.query(pts)[0]
        discs = [Point(p).buffer(rr, quad_segs=12) for p, rr in zip(pts, r) if rr > 0.004]
        return unary_union(discs + list(fills))


def compare(a, b):
    inter = a.intersection(b).area; uni = a.union(b).area
    return dict(iou=inter / uni, a_only=a.difference(b).area, b_only=b.difference(a).area)
