"""Vector re-draw of the original: exact medial axis from the outline, primitives per line, exact widths and cuts."""
import json, math, os, re, sys, subprocess, time
import numpy as np
import shapely
from shapely.geometry import Polygon, MultiPoint, LineString, Point, MultiLineString, GeometryCollection
from shapely.ops import voronoi_diagram, unary_union
from scipy.spatial import cKDTree
from geom import *

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'art'))
from lib import svgpath, fmt

STEP = 0.03          # outline sampling step (units)
TOL = 0.02           # primitive fitting tolerance
t_start = time.time()

# ---------------- original as an exact polygon (nonzero winding, like the SVG renderer) ----------------
svg = open(os.path.join(HERE, '..', 'logo.svg')).read()
d = re.search(r' d="([^"]+)"', svg).group(1)
subs = sample_path(d, STEP)
raw_rings = [pts for pts, closed in subs]


def winding(pt, ring):
    """winding number of pt w.r.t. a closed polyline (crossing-number with direction)"""
    x, y = pt
    xs, ys = ring[:, 0], ring[:, 1]
    x2, y2 = np.roll(xs, -1), np.roll(ys, -1)
    up = (ys <= y) & (y2 > y)
    dn = (ys > y) & (y2 <= y)
    # cross product sign: is pt left of the edge?
    cross = (x2 - xs) * (y - ys) - (x - xs) * (y2 - ys)
    return int(np.sum(up & (cross > 0)) - np.sum(dn & (cross < 0)))


noded = unary_union([LineString(np.vstack([r, r[:1]])) for r in raw_rings])
from shapely.ops import polygonize as _polygonize
faces = list(_polygonize(noded))
keep = []
for f in faces:
    c = f.representative_point()
    wn = sum(winding((c.x, c.y), r) for r in raw_rings)
    if wn != 0:
        keep.append(f)
ink = unary_union(keep)
print(f'original polygon (nonzero winding): area {ink.area:.2f} u², faces {len(faces)} kept {len(keep)}, parts {len(getattr(ink, "geoms", [ink]))}')
# boundary rings of the final polygon, resampled uniformly
rings = []
for poly in getattr(ink, 'geoms', [ink]):
    for ring in [poly.exterior, *poly.interiors]:
        pts = resample(np.asarray(ring.coords), STEP)
        rings.append((pts, True))
islands = [np.asarray(p.exterior.coords) for p in getattr(ink, 'geoms', [ink]) if p.area < 15 and len(p.interiors) == 0]
print('rings', len(rings), 'islands', len(islands))

# boundary samples with ring ids and positions
bpts, bring, bpos = [], [], []
for ri, (pts, cw) in enumerate(rings):
    bpts.append(pts); bring.append(np.full(len(pts), ri)); bpos.append(np.arange(len(pts)))
B = np.vstack(bpts); BR = np.concatenate(bring); BP = np.concatenate(bpos)
ring_len = {ri: len(pts) for ri, (pts, cw) in enumerate(rings)}
tree = cKDTree(B)
print('boundary samples', len(B))

# ---------------- medial axis via Voronoi ----------------
vd = voronoi_diagram(MultiPoint(B), edges=True)
edges = []
def _walk(g):
    if g.geom_type == 'LineString':
        c = np.asarray(g.coords)
        for a, b in zip(c[:-1], c[1:]):
            edges.append((a, b))
    else:
        for h in getattr(g, 'geoms', []):
            _walk(h)
_walk(vd)
E = np.array(edges)                      # (n, 2, 2)
mid = E.mean(axis=1)
inside = shapely.contains_xy(ink, mid[:, 0], mid[:, 1]) & shapely.contains_xy(ink, E[:, 0, 0], E[:, 0, 1]) & shapely.contains_xy(ink, E[:, 1, 0], E[:, 1, 1])
E, mid = E[inside], mid[inside]
dist2, idx2 = tree.query(mid, k=2)
same_ring = BR[idx2[:, 0]] == BR[idx2[:, 1]]
gap = np.abs(BP[idx2[:, 0]] - BP[idx2[:, 1]])
rl = np.array([ring_len[r] for r in BR[idx2[:, 0]]])
gap = np.minimum(gap, rl - gap)
spur = same_ring & (gap <= 3)
E = E[~spur]
print(f'voronoi edges inside the ink: {inside.sum()}, after dropping neighbour spurs: {len(E)}  ({time.time() - t_start:.0f}s)')

# ---------------- graph ----------------
def key(p):
    return (round(float(p[0]), 6), round(float(p[1]), 6))

adj = {}
for a, b in E:
    ka, kb = key(a), key(b)
    if ka == kb:
        continue
    adj.setdefault(ka, set()).add(kb)
    adj.setdefault(kb, set()).add(ka)


def radius(p):
    return float(tree.query(np.asarray(p, float))[0])


def extract_branches():
    seen = set()
    branches = []
    nodes = [k for k, nb in adj.items() if len(nb) != 2]
    for n0 in nodes:
        for nb in list(adj[n0]):
            if (n0, nb) in seen:
                continue
            path = [n0, nb]
            seen.add((n0, nb)); seen.add((nb, n0))
            prev, cur = n0, nb
            while len(adj[cur]) == 2:
                nxt = [x for x in adj[cur] if x != prev][0]
                if (cur, nxt) in seen:
                    break
                seen.add((cur, nxt)); seen.add((nxt, cur))
                path.append(nxt); prev, cur = cur, nxt
                if cur == n0:
                    break
            branches.append(np.array(path, float))
    # pure cycles (no junction nodes) are rare here; ignore
    return branches


def length(pts):
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))


branches = extract_branches()
raw_branches = [b.copy() for b in branches]
raw_nodes = [np.array(k) for k, nb in adj.items() if len(nb) != 2]
print('raw branches', len(branches))

# prune cap/corner spurs: a free end whose radius ~0 (runs into a convex corner), short
def deg(k):
    return len(adj.get(k, ()))


def remove_branch(pts):
    ks = [key(p) for p in pts]
    for a, b in zip(ks[:-1], ks[1:]):
        adj.get(a, set()).discard(b); adj.get(b, set()).discard(a)
    for k in ks:
        if k in adj and not adj[k]:
            del adj[k]


for _ in range(4):
    changed = False
    for pts in extract_branches():
        k0, k1 = key(pts[0]), key(pts[-1])
        d0, d1 = deg(k0), deg(k1)
        L = length(pts)
        if L < 1.6 and ((d0 == 1 and radius(pts[0]) < 0.25 and d1 >= 3) or (d1 == 1 and radius(pts[-1]) < 0.25 and d0 >= 3)):
            remove_branch(pts); changed = True
        elif L < 0.5 and d0 == 1 and d1 == 1:
            remove_branch(pts); changed = True
    if not changed:
        break
branches = extract_branches()
print('branches after pruning spurs', len(branches))

# merge: through degree-2 nodes always; through junctions when nearly straight
ends = {i: (key(p[0]), key(p[-1])) for i, p in enumerate(branches)}
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
    for n in nodes:
        inc = node_branches(n)
        if len(inc) < 2:
            continue
        best = None
        for i in range(len(inc)):
            for j in range(i + 1, len(inc)):
                (b1, e1), (b2, e2) = inc[i], inc[j]
                if b1 == b2:
                    continue
                d1 = end_tangent(bmap[b1], e1 == 0, 0.6)
                d2 = end_tangent(bmap[b2], e2 == 0, 0.6)
                ang = math.degrees(math.acos(float(np.clip(-np.dot(d1, d2), -1, 1))))
                lim = 181 if len(inc) == 2 else 30
                if ang < lim and (best is None or ang < best[0]):
                    best = (ang, b1, e1, b2, e2)
        if best is None:
            continue
        _, b1, e1, b2, e2 = best
        p1 = bmap[b1] if e1 == 1 else bmap[b1][::-1]
        p2 = bmap[b2] if e2 == 0 else bmap[b2][::-1]
        merged = np.vstack([p1, p2[1:]])
        new_ends = (ends[b1][0] if e1 == 1 else ends[b1][1], ends[b2][1] if e2 == 0 else ends[b2][0])
        del bmap[b2], ends[b2]
        bmap[b1] = merged; ends[b1] = new_ends
        if len(inc) > 2:
            through.add(n)
        return True
    return False


while merge_pass():
    pass
print('strokes after merging', len(bmap), f'({time.time() - t_start:.0f}s)')

# ---------------- per stroke: width, ends, primitives ----------------
boundary = ink.boundary
prep = shapely.prepared.prep(ink)


def ray_hit(p, dvec, maxu=3.0):
    """distance from p along dvec to the first boundary crossing"""
    ls = LineString([p, p + dvec * maxu])
    g = ls.intersection(boundary)
    if g.is_empty:
        return maxu
    pts = [np.asarray(q.coords[0]) for q in getattr(g, 'geoms', [g]) if q.geom_type == 'Point']
    if not pts:
        return maxu
    return float(min(np.linalg.norm(q - p) for q in pts))


def ser(v):
    if isinstance(v, np.ndarray):
        return [round(float(x), 4) for x in v.ravel()] if v.ndim == 1 else [[round(float(x), 4) for x in r] for r in v]
    if isinstance(v, (list, tuple)):
        return [ser(x) for x in v]
    if isinstance(v, (np.floating, float)):
        return round(float(v), 4)
    return v

ARROW = Polygon([[0, 16.04], [9.225, 10.155], [7.923, 13.166], [8.322, 14.543], [10.917, 16.543]]).buffer(0.05)
ISLANDS = unary_union([Polygon(p).buffer(0.05) for p in islands])
def split_at_step(raw):
    """if the width profile has one clear step, split the polyline there -> list of (pts, forced_junction_flags)"""
    pts = resample(raw, 0.05)
    if len(pts) < 40:
        return [pts]
    r = np.array([radius(p) for p in pts])
    sm = np.array([np.median(r[max(0, i - 6):i + 7]) for i in range(len(r))])
    n = len(sm)
    best = None
    for k in range(12, n - 12):
        a, b = sm[:k], sm[k:]
        cost = a.var() * len(a) + b.var() * len(b)
        if best is None or cost < best[0]:
            best = (cost, k)
    k = best[1]
    a, b = sm[:k], sm[k:]
    step = abs(np.median(a) - np.median(b)) * 2
    if step > 0.09 and (np.percentile(a, 80) - np.percentile(a, 20)) * 2 < 0.08 and (np.percentile(b, 80) - np.percentile(b, 20)) * 2 < 0.08 \
            and length(pts[:k + 1]) > 1.5 and length(pts[k:]) > 1.5:
        return [pts[:k + 1], pts[k:]]
    return [pts]


strokes = []
work = []
for bid, raw in bmap.items():
    if np.mean(shapely.contains_xy(ARROW, raw[:, 0], raw[:, 1])) > 0.6 or np.mean(shapely.contains_xy(ISLANDS, raw[:, 0], raw[:, 1])) > 0.6:
        continue
    k0, k1 = ends[bid]
    f0, f1 = deg(k0) == 1, deg(k1) == 1
    parts = split_at_step(raw)
    for j, part in enumerate(parts):
        work.append((bid * 10 + j, part, f0 if j == 0 else False, f1 if j == len(parts) - 1 else False, k0 if j == 0 else None, k1 if j == len(parts) - 1 else None))
print('strokes split at width steps:', sum(1 for b in bmap if len(split_at_step(bmap[b])) > 1))
for bid, pts, free0, free1, k0, k1 in work:
    if len(pts) < 3:
        continue
    # width profile from the exact distance to the boundary
    rad = np.array([radius(p) for p in pts])
    n = len(pts)
    win = 20                                   # samples (0.05 u each) -> +-1.0 u window
    prof = np.array([2 * np.percentile(rad[max(0, i - win):i + win + 1], 20) for i in range(n)])
    axis0 = pts.copy()
    core = rad[int(n * 0.2):max(int(n * 0.8), int(n * 0.2) + 1)]
    w = 2 * float(np.percentile(core, 30)) if len(core) else 2 * float(np.median(rad))
    spread = float(np.percentile(core, 90) - np.percentile(core, 10)) * 2 if len(core) > 4 else 0.0
    pct = [round(2 * float(np.percentile(core, q)), 2) for q in (10, 30, 50, 70, 90)] if len(core) > 4 else []
    # junction ends: replace the last hw*1.2 of the axis (which bends into the junction) by a straight extension
    hw = w / 2
    split_end0, split_end1 = (k0 is None), (k1 is None)
    def straighten(pts, at_start):
        cut = hw * 1.25
        if at_start:
            d = np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1)); k = int(np.searchsorted(d, cut)) + 1
            if k >= len(pts) - 3:
                return pts
            t = end_tangent(pts[k:], True, 0.6)
            return np.vstack([pts[k] - t * (cut + 0.15), pts[k:]])
        d = np.cumsum(np.linalg.norm(np.diff(pts[::-1], axis=0), axis=1)); k = int(np.searchsorted(d, cut)) + 1
        if k >= len(pts) - 3:
            return pts
        t = end_tangent(pts[:len(pts) - k], False, 0.6)
        return np.vstack([pts[:len(pts) - k], [pts[len(pts) - 1 - k] - t * (cut + 0.15)]])
    if not free0 and not split_end0:
        pts = straighten(pts, True)
    if not free1 and not split_end1:
        pts = straighten(pts, False)
    # free ends: exact cap corners from the vector boundary
    caps = [None, None]
    for end in (0, 1):
        if not (free0 if end == 0 else free1):
            continue
        p = pts[0] if end == 0 else pts[-1]
        t = -end_tangent(pts, end == 0, 0.5)          # outward
        nn = np.array([-t[1], t[0]])
        dc = ray_hit(p, t)
        rl = min(hw, radius(p)) - 0.1
        pl, pr = p + nn * rl, p - nn * rl
        dl = ray_hit(pl, t) if ink.contains(Point(pl)) else dc
        dr = ray_hit(pr, t) if ink.contains(Point(pr)) else dc
        dl = dc if (dl < 0.02 or dl >= 2.9) else dl
        dr = dc if (dr < 0.02 or dr >= 2.9) else dr
        dc = min(dc, 2.9)
        caps[end] = (float(dl), float(dr))
        ext = max(dl, dr)
        if end == 0:
            pts = np.vstack([p + t * ext, pts])
        else:
            pts = np.vstack([pts, p + t * ext])
    # primitives: split at corners, classify each piece
    cuts = [0] + corner_indices(pts) + [len(pts) - 1]
    prims = []
    for a, b in zip(cuts, cuts[1:]):
        piece = pts[a:b + 1]
        if len(piece) < 2:
            continue
        kind, data = classify(piece, TOL)
        def flat(k, dd):
            if k == 'multi':
                for k2, d2 in dd:
                    yield from flat(k2, d2)
            else:
                yield (k, dd)
        prims.extend(flat(kind, data))
    # variable width when the (junction-free) profile really changes along the stroke
    pr_core = prof[int(n * 0.1):max(int(n * 0.9), int(n * 0.1) + 1)]
    var_w = (pr_core.max() - pr_core.min()) > 0.1 if len(pr_core) > 10 else False
    strokes.append(dict(id=int(bid), width=round(w, 3), spread=round(spread, 3), pct=pct, len=round(length(pts), 2),
                        free=[bool(free0), bool(free1)], caps=caps, prims=prims, var_w=bool(var_w),
                        axis=axis0, prof=prof.copy()))
    if np.linalg.norm(axis0 - np.array([23.5, 21.0]), axis=1).min() < 1.0:
        json.dump(dict(axis=axis0.tolist(), pts=pts.tolist(), prims=[[k, ser(dd)] for k, dd in prims], caps=caps, w=w, free=[bool(free0), bool(free1)]),
                  open(os.path.join(HERE, 'debug_shoulder.json'), 'w'))
print('strokes', len(strokes), 'primitives:', sum(len(s['prims']) for s in strokes),
      'lines', sum(1 for s in strokes for k, _ in s['prims'] if k == 'line'),
      'arcs', sum(1 for s in strokes for k, _ in s['prims'] if k == 'arc'),
      'beziers', sum(len(dd) for s in strokes for k, dd in s['prims'] if k == 'bezier'),
      'polys', sum(1 for s in strokes for k, _ in s['prims'] if k == 'poly'), f'({time.time() - t_start:.0f}s)')

# ---------------- render: exact buffers with mitred joins, cut caps ----------------
geoms = []
for s in strokes:
    pts = []
    for kind, data in s['prims']:
        q = primitive_points(kind, data, 0.02)
        pts.extend(q if not pts else q[1:])
    pts = np.array(pts)
    hw = s['width'] / 2
    if s['var_w']:
        ax = s['axis']; pf = s['prof']
        # width at each rendered point = profile at the nearest axis sample
        tr = cKDTree(ax)
        wq = pf[tr.query(pts)[1]]
        tan = np.gradient(pts, axis=0); tan /= (np.linalg.norm(tan, axis=1)[:, None] + 1e-12)
        nor = np.stack([-tan[:, 1], tan[:, 0]], 1)
        left = pts + nor * (wq / 2)[:, None]; right = pts - nor * (wq / 2)[:, None]
        g = Polygon(np.vstack([left, right[::-1]])).buffer(0)
        hw = float(wq[0] / 2) if True else hw
        hw_end = float(wq[-1] / 2)
    else:
        g = LineString(pts).buffer(hw, cap_style='flat', join_style='mitre', mitre_limit=6, quad_segs=16)
        hw_end = hw
    for end in (0, 1):
        c = s['caps'][end]
        if c is None or abs(c[0] - c[1]) < 0.01:
            continue
        dl, dr = c
        p = pts[0] if end == 0 else pts[-1]
        t = -end_tangent(pts, end == 0, 0.5)
        nn = np.array([-t[1], t[0]])
        m = max(dl, dr)
        h_ = hw if end == 0 else hw_end
        cl = p - t * (m - dl) + nn * h_
        cr = p - t * (m - dr) - nn * h_
        # keep the half-plane behind the cut line (towards the stroke)
        u = cr - cl; u /= np.linalg.norm(u)
        back = -t
        big = Polygon([cl + u * 5 + back * 0.0, cr - u * 5, cr - u * 5 + back * 6, cl + u * 5 + back * 6])
        g = g.intersection(big.union(Polygon([cl, cr, cr + back * 6, cl + back * 6]).buffer(0)))
    geoms.append(g)
fills = [Polygon([[0, 16.04], [9.225, 10.155], [7.923, 13.166], [8.322, 14.543], [10.917, 16.543]])]
for p in islands:
    fills.append(Polygon(p).buffer(0))
mine_lines = unary_union(geoms + fills)
t_e = time.time()
exact = []
for b in raw_branches:
    pts = resample(b, 0.02)
    if len(pts) < 2:
        continue
    r = np.array([radius(p) for p in pts])
    tan = np.gradient(pts, axis=0); tan /= (np.linalg.norm(tan, axis=1)[:, None] + 1e-12)
    nor = np.stack([-tan[:, 1], tan[:, 0]], 1)
    left = pts + nor * r[:, None]; right = pts - nor * r[:, None]
    poly = Polygon(np.vstack([left, right[::-1]])).buffer(0)
    if not poly.is_empty:
        exact.append(poly)
for nd in raw_nodes:
    rr = radius(nd)
    if rr > 0.01:
        exact.append(Point(nd).buffer(rr, quad_segs=24))
mine = unary_union(exact + fills)
# diagnostics: are the big node radii real?
for nd in sorted(raw_nodes, key=lambda p: -radius(p))[:6]:
    print(f'  node ({nd[0]:.2f}, {nd[1]:.2f}) kd-radius {radius(nd):.3f}  true boundary distance {boundary.distance(Point(nd)):.3f}  inside {ink.contains(Point(nd))}')
ex_rank = sorted(((g.difference(ink).area, i) for i, g in enumerate(exact)), key=lambda x: -x[0])[:8]
for a, i in ex_rank:
    g = exact[i]; c = g.centroid
    print(f'  overshoot {a:5.2f} u² part {i} ({"node disc" if i >= len(raw_branches) else "branch"}) near ({c.x:.1f}, {c.y:.1f}) bbox {[round(v, 1) for v in g.bounds]}')
print(f'exact MAT render: original {ink.area:.2f}, mine {mine.area:.2f}, original-only {ink.difference(mine).area:.2f}, '
      f'mine-only {mine.difference(ink).area:.2f}, IoU {mine.intersection(ink).area / mine.union(ink).area:.4f}  ({time.time() - t_e:.0f}s)')
t_d = time.time()
disc_pts = []
for b in raw_branches:
    pts = resample(b, 0.03)
    disc_pts.extend(pts)
disc_pts.extend(raw_nodes)
disc_pts = np.array(disc_pts)
# true distance to the boundary for every disc centre (vectorised GEOS distance against the boundary rings)
rings_geom = [LineString(p) for p, cw in rings]
from shapely.strtree import STRtree
ring_tree = STRtree(rings_geom)
pts_geom = shapely.points(disc_pts[:, 0], disc_pts[:, 1])
cand = ring_tree.query_nearest(pts_geom, all_matches=False)
true_r = np.array([rings_geom[j].distance(pts_geom[i]) for i, j in zip(cand[0], cand[1])])
kd_r = tree.query(disc_pts)[0]
print(f'disc centres {len(disc_pts)}: kd radius vs true radius: mean diff {np.mean(kd_r - true_r):+.4f}, max diff {np.max(kd_r - true_r):.4f}  ({time.time() - t_d:.0f}s)')
discs = [Point(p).buffer(r, quad_segs=12) for p, r in zip(disc_pts, true_r) if r > 0.004]
mine_disc = unary_union(discs + fills)
print(f'disc-union render: original {ink.area:.2f}, mine {mine_disc.area:.2f}, original-only {ink.difference(mine_disc).area:.2f}, '
      f'mine-only {mine_disc.difference(ink).area:.2f}, IoU {mine_disc.intersection(ink).area / mine_disc.union(ink).area:.4f}  ({time.time() - t_d:.0f}s)')
mine_cls = mine_lines
print(f'classified-lines render: IoU {mine_cls.intersection(ink).area / mine_cls.union(ink).area:.4f}')
mine = mine_disc
# construction view: classified centerlines coloured by weight class over the exact ink
col = {'light': '#2a9df4', 'standard': '#1c1b19', 'heavy': '#d6473b'}
def wclass(w):
    return 'light' if w < 1.27 else ('standard' if w < 1.5 else 'heavy')
cparts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 54 64" width="540" height="640">',
          f'<path d="{svgpath(mine, prec=3, tol=0.004)}" fill="#dedbd3" fill-rule="evenodd"/>']
Ls = {'light': 0.0, 'standard': 0.0, 'heavy': 0.0}; Ns = {'light': 0, 'standard': 0, 'heavy': 0}
for s in strokes:
    pts = np.vstack([primitive_points(k, dd, 0.05) for k, dd in s['prims']])
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in pts[::2])
    cls = wclass(s['width']); Ls[cls] += s['len']; Ns[cls] += 1
    cparts.append(f'<path d="{pd}" fill="none" stroke="{col[cls]}" stroke-width="0.28" stroke-linecap="butt"/>')
for p in islands + [np.array(fills[0].exterior.coords)]:
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in p[::3]) + 'Z'
    cparts.append(f'<path d="{pd}" fill="none" stroke="#7a3fb8" stroke-width="0.22"/>')
cparts.append('</svg>')
open(os.path.join(HERE, 'construction_v.svg'), 'w').write(''.join(cparts))
tot = sum(Ls.values())
for c in ('light', 'standard', 'heavy'):
    ws = [s['width'] for s in strokes if wclass(s['width']) == c]
    print(f'  {c:9s}: {Ns[c]:3d} lines, {Ls[c]:6.1f} u ({Ls[c] / tot:.0%}), median width {np.median(ws) if ws else 0:.2f}')
caps_sq = sum(1 for s in strokes for c in s['caps'] if c and abs(c[0] - c[1]) < 0.06)
caps_an = sum(1 for s in strokes for c in s['caps'] if c and abs(c[0] - c[1]) >= 0.06)
print(f'  free ends: {caps_sq + caps_an}, square {caps_sq}, angled {caps_an}; primitives: lines {sum(1 for s in strokes for k, _ in s["prims"] if k == "line")}, '
      f'arcs {sum(1 for s in strokes for k, _ in s["prims"] if k == "arc")}, curves {sum(1 for s in strokes for k, _ in s["prims"] if k in ("bezier", "poly"))}')
# diagnostics: largest regions of the original not covered
miss = ink.difference(mine)
parts = sorted(getattr(miss, 'geoms', [miss]), key=lambda g: -g.area)[:8]
for g in parts:
    c = g.centroid
    near = sorted(strokes, key=lambda s: Point(primitive_points(*s['prims'][0])[0]).distance(c))[:1]
    print(f"  missing {g.area:5.2f} u² at ({c.x:.1f}, {c.y:.1f}) bbox {[round(v, 1) for v in g.bounds]}")
# diagnostics: radius profiles of the longest strokes (percentiles of 2*radius over the core)
for s in sorted(strokes, key=lambda s: -s['len'])[:10]:
    print(f"  profile stroke {s['id']:3d} len {s['len']:6.2f}: w chosen {s['width']:.2f}  p10/30/50/70/90 = {s['pct']}")
# diagnostics: which strokes add ink the original does not have
excess = sorted(((g.difference(ink).area, s) for g, s in zip(geoms, strokes)), key=lambda x: -x[0])
for a, s in excess[:4]:
    print(f"  excess {a:5.2f} u²  stroke {s['id']:3d} len {s['len']:6.2f} w {s['width']:.2f} spread {s['spread']:.2f} free {s['free']} caps {s['caps']} prims {[k for k, _ in s['prims']]}")
sym = mine.symmetric_difference(ink)
print(f'vector compare: original {ink.area:.2f} u², mine {mine.area:.2f} u², original-only {ink.difference(mine).area:.2f} u², '
      f'mine-only {mine.difference(ink).area:.2f} u², IoU {mine.intersection(ink).area / mine.union(ink).area:.4f}  ({time.time() - t_start:.0f}s)')

# save drawing + svgs
json.dump(dict(strokes=[dict(id=s['id'], width=s['width'], spread=s['spread'], len=s['len'], free=s['free'], caps=ser(s['caps']), var_w=s['var_w'],
                              prims=[[k, ser(dd)] for k, dd in s['prims']]) for s in strokes]),
          open(os.path.join(HERE, 'drawing_v.json'), 'w'))
path = svgpath(mine, prec=4, tol=0.002)
open(os.path.join(HERE, 'redraw_v.svg'), 'w').write(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 54 64" width="1296" height="1536"><rect x="0" y="0" width="54" height="64" fill="#fff"/><path d="{path}" fill="#000" fill-rule="evenodd"/></svg>')
open(os.path.join(HERE, 'redraw_v-clean.svg'), 'w').write(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 54 64" width="540" height="640"><path d="{path}" fill="#1c1b19" fill-rule="evenodd"/></svg>')
# diff svg: original red, mine blue, overlap grey (exact vector overlay)
both = mine.intersection(ink)
open(os.path.join(HERE, 'overlay_v.svg'), 'w').write(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 54 64" width="1296" height="1536"><rect x="0" y="0" width="54" height="64" fill="#fff"/>'
    f'<path d="{svgpath(ink.difference(mine), prec=4, tol=0.002)}" fill="#d6473b" fill-rule="evenodd"/>'
    f'<path d="{svgpath(mine.difference(ink), prec=4, tol=0.002)}" fill="#2a9df4" fill-rule="evenodd"/>'
    f'<path d="{svgpath(both, prec=4, tol=0.002)}" fill="#c6c3bc" fill-rule="evenodd"/></svg>')
print('done', f'({time.time() - t_start:.0f}s)')

# pickle the result for the Step 2 composition
import pickle as _pk
_pk.dump(dict(strokes=strokes, geoms=geoms, mine=mine, ink=ink, fills=[Polygon(p).buffer(0) for p in islands] if 'islands' in dir() else []),
         open(os.path.join(HERE, 'step1.pkl'), 'wb'))
print('pickled step1.pkl')
