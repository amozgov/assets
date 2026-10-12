"""Step 2 final assembly, v5: exact handovers.
Body = the original's exact ink minus the removed strokes. Every kept line that the removed ink touched is protected by its
own exact band (discs of the medial axis' true radius, carried straight across the removed crossings), so nothing of the
old ink survives and nothing of the kept ink is nibbled. Every join is cut on a line across the kept stroke and the new
band starts exactly at the ink's own edge points on that line. The arm pose comes from the raster: arms, hands, thorns,
the fold (chest line + light lower leg) and the armpit bend."""
import os, sys, pickle, math, numpy as np, shapely
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import unary_union, substring
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from lib import svgpath, fmt
from geom import *
from geom import _line_resid, reparam
from mat import MAT, compare
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = pickle.load(open(os.path.join(HERE, 'step1.pkl'), 'rb'))
C2 = pickle.load(open(os.path.join(HERE, 'compose2.pkl'), 'rb'))
C4 = pickle.load(open(os.path.join(HERE, 'compose4.pkl'), 'rb'))
P_ = pickle.load(open(os.path.join(HERE, 'prep.pkl'), 'rb'))
A = pickle.load(open(os.path.join(HERE, 'align2.pkl'), 'rb'))
ink1, fills1, geoms1, strokes1 = S1['ink'], S1['fills'], S1['geoms'], S1['strokes']
MAT_ = A['mat']; K = A['k']
ink2v, boltv = C4['ink2v'], C4['boltv']
tris = C4['spikes']
pieces = {p['id']: p for p in C4['pieces']}
old93_line = LineString(pieces['r93']['prims'][0][1].copy())
from scipy.interpolate import splprep, splev
def smooth_open(q, dev=0.05):
    """a C2 smoothing spline through a traced polyline: no kinks between primitives, rms deviation about `dev`"""
    q = np.asarray(q, float); q = resample(q, 0.05)
    tck, _ = splprep([q[:, 0], q[:, 1]], s=len(q) * dev ** 2, k=3)
    x, y = splev(np.linspace(0, 1, max(int(length(q) / 0.03), 12)), tck); return np.stack([x, y], 1)
for _pid in ('r11', 'r13', 'r93', 'r103'):
    pieces[_pid]['prims'] = [('poly', smooth_open(pieces[_pid]['prims'][0][1]))]
HAND_RING = {h: MAT.render_clean(None, [pieces[h + '_oval']], []) for h in ('handL', 'handR')}
HAND_FILL = {h: Polygon(g.exterior) for h, g in HAND_RING.items()}
HAND_HOLE = {h: Polygon(max(g.interiors, key=lambda r: Polygon(r).area)) for h, g in HAND_RING.items()}
def to_loop(q, hand, run=1.0):
    q = np.asarray(q, float); t = unit(q[-1] - q[max(len(q) - 6, 0)]); return np.vstack([q, [q[-1] + t * run]])
R = {st['id']: st for st in P_['strokes2']}
O = {st['id']: st for st in strokes1}; G = {st['id']: g for st, g in zip(strokes1, geoms1)}
FORCE_DROP = {1430, 320, 540, 880}     # the original's single armpit thorn axis, the fold's old stub, the bun junction's 0.6 u stump, the neck's junk continuation
LIGHT, STD = 1.15, 1.38
notes = []


def perp(u): return np.array([-u[1], u[0]])
def cross2(a, b): return float(a[0] * b[1] - a[1] * b[0])
def unit(v): return v / (np.linalg.norm(v) + 1e-12)
def buf(pts, h): return LineString(pts).buffer(h, cap_style='flat', join_style='mitre', mitre_limit=8, quad_segs=16)
def halfplane(p0, nrm, size=40.0):
    nrm = unit(nrm); t = perp(nrm)
    return Polygon([p0 + t * size, p0 - t * size, p0 - t * size + nrm * size, p0 + t * size + nrm * size])
def ray_hit(p, u, q, v):
    M = np.array([[u[0], -v[0]], [u[1], -v[1]]])
    if abs(np.linalg.det(M)) < 1e-6: return None
    s, t = np.linalg.solve(M, q - p); return float(s), float(t)
def arc_index(pts, L):
    d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]); return int(np.searchsorted(d, L))
def cum(pts): return np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
def at_len(pts, L):
    """point of the polyline at arc length L (interpolated)"""
    d = cum(pts); L = min(max(L, 0.0), d[-1]); i = int(np.searchsorted(d, L)); i = min(max(i, 1), len(pts) - 1)
    f = (L - d[i - 1]) / max(d[i] - d[i - 1], 1e-12); return pts[i - 1] + (pts[i] - pts[i - 1]) * f
def fit_dir(pts, L0, L1):
    i0, i1 = arc_index(pts, L0), arc_index(pts, L1); q = pts[i0:max(i1, i0 + 2)]
    m = q.mean(axis=0); _, _, vt = np.linalg.svd(q - m); d = vt[0]
    if np.dot(d, q[-1] - q[0]) < 0: d = -d
    return unit(d), m
def spts(st, xf=None):
    pts = []
    for k, dd in st['prims']:
        q = primitive_points(k, dd, 0.02); pts.extend(q if not pts else q[1:])
    pts = np.array(pts); return xf(pts) if xf else pts
def flat(k, dd):
    if k == 'multi':
        for k2, d2 in dd: yield from flat(k2, d2)
    else: yield (k, dd)
def seg_near(geom, p, maxlen=2.6, maxdist=0.3):
    """the piece of a line-intersection result nearest p"""
    segs = [s for s in getattr(geom, 'geoms', [geom]) if s.geom_type == 'LineString' and s.length < maxlen and s.distance(Point(p)) < maxdist]
    return min(segs, key=lambda s: s.distance(Point(p))) if segs else None
def mitre_straight(c, uA, uB, hA, hB, LA=3.0, LB=3.0):
    """exact mitre of two straight legs meeting at c: A arrives along uA (unit, pointing into c), B leaves along uB.
    Returns the polygon of both legs (lengths LA, LB) with the mitred corner, and the outer / inner corner points."""
    nA, nB = perp(uA), perp(uB)
    turn = 1.0 if cross2(uA, uB) > 0 else -1.0
    pts = []
    for side in (-turn, turn):
        pa = c + nA * side * hA; pb = c + nB * side * hB
        st = ray_hit(pa, uA, pb, uB)
        P = pa + uA * st[0] if st is not None else (pa + pb) / 2
        pts.append(P)
    Po, Pi = pts
    poly = Polygon([c - uA * LA + nA * (-turn) * hA, Po, c + uB * LB + nB * (-turn) * hB,
                    c + uB * LB + nB * turn * hB, Pi, c - uA * LA + nA * turn * hA]).buffer(0)
    return poly, Po, Pi
def mitre_join(Apts, Bpts, hA, hB, ext=2.5):
    """exact mitred join of leg A (ending near the corner) and leg B (starting near it), each at its own half-width."""
    uA = unit(Apts[-1] - Apts[max(len(Apts) - 8, 0)]); uB = unit(Bpts[min(7, len(Bpts) - 1)] - Bpts[0])
    A2 = np.vstack([Apts, Apts[-1] + uA * ext]); B2 = np.vstack([Bpts[0] - uB * ext, Bpts])
    la, lb = LineString(A2), LineString(B2)
    X = la.intersection(lb)
    cands = [np.asarray(g.coords[0]) for g in getattr(X, 'geoms', [X]) if g.geom_type == 'Point']
    if not cands:
        return buf(np.vstack([Apts, Bpts]), max(hA, hB)), None
    c = min(cands, key=lambda p: np.linalg.norm(p - Apts[-1]))
    if abs(cross2(uA, uB)) < 0.03:
        return buf(np.vstack([A2[:-1], [c], B2[1:]]), max(hA, hB)), c
    turn = 1.0 if cross2(uA, uB) > 0 else -1.0
    corners = {}
    for side in (-1.0, 1.0):
        oa = la.offset_curve(side * hA, join_style='mitre', mitre_limit=8); ob = lb.offset_curve(side * hB, join_style='mitre', mitre_limit=8)
        Y = oa.intersection(ob)
        pts = [np.asarray(g.coords[0]) for g in getattr(Y, 'geoms', [Y]) if g.geom_type == 'Point']
        corners[side] = min(pts, key=lambda p: np.linalg.norm(p - c)) if pts else None
    Po, Pi = corners[-turn], corners[turn]
    if Po is None or Pi is None:
        return buf(np.vstack([A2[:-1], [c], B2[1:]]), max(hA, hB)), c
    nsplit = perp(unit(Pi - Po))
    hpA = halfplane(Po, nsplit if np.dot(nsplit, Apts[0] - Po) > 0 else -nsplit)
    hpB = halfplane(Po, -nsplit if np.dot(nsplit, Apts[0] - Po) > 0 else nsplit)
    poly = buf(A2, hA).intersection(hpA).union(buf(B2, hB).intersection(hpB)).buffer(0)
    return poly, c
def end_dir(sid, at_end=True):
    """unit tangent of stroke sid at its end (pointing out of the stroke) or at its start (pointing into it)"""
    q = spts(O[sid]); kind, d = O[sid]['prims'][-1 if at_end else 0]
    if kind == 'line':
        a, b = np.asarray(d[0]), np.asarray(d[1]); return unit(b - a)
    if kind == 'arc':
        cx, cy, r, a0, a1 = d; a = a1 if at_end else a0; s = 1.0 if a1 > a0 else -1.0
        return unit(np.array([-math.sin(a), math.cos(a)]) * s)
    return unit(q[-1] - q[arc_index(q, length(q) - 0.6)]) if at_end else unit(q[arc_index(q, 0.6)] - q[0])
def local_dir(q, L, span=0.25):
    """direction of polyline q around arc length L"""
    return unit(at_len(q, L + span) - at_len(q, max(L - span, 0.0)))


# ---------------- kept / dropped ----------------
kept = [r for r in C2['kept_recs'] if r['id'] not in FORCE_DROP]
kept.append(dict(id=560, st=O[560], g=G[560], pts=spts(O[560]), kind='whole', cut=[False, False]))   # hair strand under the old spear: Step 1 keeps it
kept_ids = {r['id'] for r in kept}
dropped_ids = [i for i in O if i not in kept_ids]
BND = ink1.boundary
RT = {}
def r_true(sid):
    if sid not in RT:
        ax = np.asarray(O[sid]['axis'], float); RT[sid] = shapely.distance(BND, shapely.points(ax))
    return RT[sid]
def discs(ax, r, qs=16):
    return unary_union([Point(p).buffer(rr, quad_segs=qs) for p, rr in zip(ax, r) if rr > 0.004])
def axis_tan(sid, P, side=0, span=0.7):
    """unit tangent of the stroke's medial axis at the axis point nearest P, along the axis (increasing arc length):
    a local quadratic fit over `span` u on the clean side only (side -1: before P, +1: after P, 0: both sides)"""
    ax = np.asarray(O[sid]['axis'], float); d = cum(ax); k = int(np.argmin(np.linalg.norm(ax - P, axis=1)))
    lo, hi = (d[k] - span, d[k] + 0.02) if side < 0 else ((d[k] - 0.02, d[k] + span) if side > 0 else (d[k] - span, d[k] + span))
    m = (d >= lo) & (d <= hi)
    if m.sum() < 5: return unit(ax[min(k + 1, len(ax) - 1)] - ax[max(k - 1, 0)])
    s = d[m] - d[k]; cx = np.polyfit(s, ax[m, 0], 2); cy = np.polyfit(s, ax[m, 1], 2)
    return unit(np.array([cx[1], cy[1]]))
def raw_band(sid): return discs(np.asarray(O[sid]['axis'], float), r_true(sid))
# the removed strokes' ink: their exact discs, plus the classified polygon near them (junction ends and genuine caps; not the long bogus cap extensions)
DZ_core = unary_union([raw_band(j).buffer(0.0).union(G[j].intersection(raw_band(j).buffer(0.9))) for j in dropped_ids if len(O[j]['axis']) > 1] + [G[j] for j in dropped_ids if len(O[j]['axis']) <= 1])
KNOTS = unary_union([Point(p).buffer(r_) for p, r_ in [((25.08, 9.57), 1.6), ((23.8, 24.41), 1.0), ((23.2, 25.6), 1.2), (tuple(spts(O[1180])[-1]), 1.6), ((36.3, 6.5), 1.3), ((34.1, 17.0), 1.3), (tuple(spts(O[890])[0]), 1.5), (tuple(spts(O[330])[-1]), 1.2)]])
zone_all = DZ_core
EDG = {}
def edges(sid):
    """per axis sample: the original ink's own left/right edge points on the (smoothed) normal, their signed offsets,
    and the normals; cached"""
    if sid in EDG: return EDG[sid]
    st = O[sid]; ax0 = np.asarray(st['axis'], float); n = len(ax0); r0 = r_true(sid)
    T_ = np.array([unit(ax0[min(k + 5, n - 1)] - ax0[max(k - 5, 0)]) for k in range(n)]); N_ = np.stack([-T_[:, 1], T_[:, 0]], 1)
    loc = ink1.intersection(box(*(np.array(unary_union([Point(p).buffer(0.01) for p in ax0[::4]]).bounds) + np.array([-3, -3, 3, 3]))))
    L = np.zeros((n, 2)); Rr = np.zeros((n, 2)); oL = np.zeros(n); oR = np.zeros(n)
    for k in range(n):
        p = ax0[k]; nn = N_[k]; reach = r0[k] + 0.8
        X = loc.intersection(LineString([p - nn * reach, p + nn * reach]))
        segs = [s for s in getattr(X, 'geoms', [X]) if s.geom_type == 'LineString' and s.distance(Point(p)) < 0.02]
        if segs:
            s = max(segs, key=lambda g_: g_.length); a_, b_ = np.asarray(s.coords[0]), np.asarray(s.coords[-1])
            o1, o2 = float(np.dot(a_ - p, nn)), float(np.dot(b_ - p, nn))
            if o1 >= o2: L[k], Rr[k], oL[k], oR[k] = a_, b_, o1, o2
            else: L[k], Rr[k], oL[k], oR[k] = b_, a_, o2, o1
        else:
            L[k], Rr[k], oL[k], oR[k] = p + nn * r0[k], p - nn * r0[k], r0[k], -r0[k]
    EDG[sid] = (L, Rr, oL, oR, N_, T_); return EDG[sid]
BI = {}
def band_info(sid, tol=0.03):
    """(exact band, bulge, whole) for a kept stroke. The band is a ribbon whose two edges are the original ink's own
    edge points at every axis sample. Where a run of samples meets removed ink and an edge is inflated there (its offset
    above the robust profile), that edge alone is carried across the run as the straight / Hermite continuation of its
    own clean course either side (a start or end inside the removed ink takes the one clean side); the other edge keeps
    its real points. bulge = the inflated cross-sections (to be removed). whole = all junction (a fragment)."""
    if sid in BI: return BI[sid]
    st = O[sid]; ax0 = np.asarray(st['axis'], float); n = len(ax0); r0 = r_true(sid); pr = np.asarray(st['prof'], float) / 2
    bulge = Polygon(); whole = False
    if n < 4:
        E = discs(ax0, r0); BI[sid] = (E, Polygon(), zone_all.distance(Point(ax0[0])) < r0[0] + 0.15); return BI[sid]
    d = cum(ax0); L0, R0, oL, oR, N_, T_ = edges(sid); L = L0.copy(); Rr = R0.copy()
    inflL = oL > pr + tol; inflR = -oR > pr + tol; infl = inflL | inflR
    dz = shapely.distance(zone_all, shapely.points(ax0))
    runs = []; i = 0; bad = np.zeros(n, bool)
    while i < n:
        if not infl[i]: i += 1; continue
        j = i
        while j < n and infl[j]: j += 1
        if dz[i:j].min() < 0.1:                                                        # only runs that meet the removed ink; then two samples' margin
            i2, j2 = max(i - 2, 0), min(j + 2, n)
            if runs and i2 <= runs[-1][1]: runs[-1] = (runs[-1][0], j2)
            else: runs.append((i2, j2))
            bad[i2:j2] = True
        i = j
    runs = [(i, j) for i, j in runs]
    if any(j - i >= n - 2 for i, j in runs): whole = True; runs = []
    def edge_tan(pts, k, side):
        m = (d >= d[k] - 0.6) & (d <= d[k] + 0.02) if side < 0 else (d >= d[k] - 0.02) & (d <= d[k] + 0.6)
        if m.sum() < 5: return unit(pts[min(k + 1, n - 1)] - pts[max(k - 1, 0)])
        s = d[m] - d[k]; cx = np.polyfit(s, pts[m, 0], 2); cy = np.polyfit(s, pts[m, 1], 2); return unit(np.array([cx[1], cy[1]]))
    slabs = []
    for i, j in runs:
        a, b = i - 1, j
        if a < 0 and b >= n: whole = True; continue
        slabs.append(halfplane(ax0[a], T_[a]) if b >= n else (halfplane(ax0[b], -T_[b]) if a < 0 else halfplane(ax0[a], T_[a]).intersection(halfplane(ax0[b], -T_[b]))))
        for pts, fl in ((L, inflL), (Rr, inflR)):
            if not fl[i:j].any(): continue                                          # this edge is clean through the run
            if a < 0:
                t = edge_tan(pts, b, +1)
                for k in range(b): pts[k] = pts[b] - t * (d[b] - d[k])
            elif b >= n:
                t = edge_tan(pts, a, -1)
                for k in range(a + 1, n): pts[k] = pts[a] + t * (d[k] - d[a])
            else:
                Lr = d[b] - d[a]; ta = edge_tan(pts, a, -1) * Lr; tb = edge_tan(pts, b, +1) * Lr; pa, pb = pts[a].copy(), pts[b].copy()
                for k in range(a + 1, b):
                    s = (d[k] - d[a]) / Lr; h00 = 2 * s ** 3 - 3 * s ** 2 + 1; h10 = s ** 3 - 2 * s ** 2 + s; h01 = -2 * s ** 3 + 3 * s ** 2; h11 = s ** 3 - s ** 2
                    pts[k] = h00 * pa + h10 * ta + h01 * pb + h11 * tb
    if runs and d[-1] < 1.0: whole = True
    if os.environ.get('BI_DEBUG'): print(f'  band_info {sid}: n {n} runs {[(round(d[i], 2), round(d[min(j, n - 1)], 2)) for i, j in runs]} whole {whole}')
    if runs and not whole:
        bulge = unary_union([Point(ax0[k]).buffer(max(oL[k], -oR[k], r0[k]) + 0.02, quad_segs=16) for k in range(n) if bad[k]]).intersection(unary_union(slabs))
    E = unary_union([shapely.MultiPoint([L[k], L[k + 1], Rr[k + 1], Rr[k]]).convex_hull for k in range(n - 1)]).buffer(0)
    if E.geom_type != 'Polygon': E = max(E.geoms, key=lambda g_: g_.area) if hasattr(E, 'geoms') else discs(ax0, r0)
    if whole: bulge = discs(ax0, r0)
    BI[sid] = (E, bulge, whole); return BI[sid]
def clean_len(sid, from_end, margin=0.25, tol=0.04):
    """arc length along the stroke's axis of the first sample that is the stroke's own clean band, counted from the start
    or the end, moved `margin` u further into the stroke. A sample is clean when neither edge is inflated by a junction
    (offset above the robust profile and above the stroke 0.3-1.5 u further in) and the band is not tapering out."""
    st = O[sid]; ax = np.asarray(st['axis'], float); r = r_true(sid); pr = np.asarray(st['prof'], float) / 2; d = cum(ax)
    _, _, oL, oR, _, _ = edges(sid); w = np.maximum(oL, -oR)
    idx = range(len(ax) - 1, -1, -1) if from_end else range(len(ax))
    for k in idx:
        m = ((d >= d[k] - 1.5) & (d <= d[k] - 0.3)) if from_end else ((d >= d[k] + 0.3) & (d <= d[k] + 1.5))
        if m.sum() < 3: continue
        inflated = (w[k] - pr[k] > 0.03) and (w[k] > w[m].min() + tol)
        tapered = r[k] < r[m].max() - 0.08 or abs(r[k] - float(np.median(r[m]))) > 0.08
        if not inflated and not tapered:
            return (d[k] - margin) if from_end else (d[k] + margin), ax
    return (d[-1] - margin) if from_end else margin, ax
keep = {}                                  # sid -> region of the stroke that stays (the rest is removed and rebuilt)
def ink_cut(sid, P0, t, nref=None, label=''):
    """the body's band of stroke sid across the cut line through P0 (normal perp(t)): (mid, t, half-width, edge1, edge2),
    the edges on the line itself, taken from the original ink inside the stroke's exact band (what the body keeps there).
    edge1 is on the nref side when nref is given."""
    n = perp(t); E = band_info(sid)[0]
    src = ink1.intersection(E.buffer(0.012))
    s = seg_near(src.intersection(LineString([P0 - n * 3.0, P0 + n * 3.0])), P0)
    if s is None:
        s = seg_near(E.intersection(LineString([P0 - n * 3.0, P0 + n * 3.0])), P0); print(f'  ink_cut {label or sid}: no ink on the cut line; band used')
    e1, e2 = np.asarray(s.coords[0]), np.asarray(s.coords[-1])
    if nref is not None and np.dot(e1 - e2, nref) < 0: e1, e2 = e2, e1
    w = float(np.linalg.norm(e1 - e2)); st = O[sid]
    if abs(w - st['width']) > 0.2: print(f'  ink_cut {label or sid}: width at the cut {w:.2f} vs the stroke\'s {st["width"]:.2f}')
    return (e1 + e2) / 2, t, w / 2, e1, e2
def nearest_axis_pt(sid, P):
    ax = np.asarray(O[sid]['axis'], float); return ax[int(np.argmin(np.linalg.norm(ax - P, axis=1)))]


# ---------------- head: the bun outline through the old spear junction: one arc from 520 to 570 ----------------
# 520 comes down the outline, meets the shaft, the strand 560 and 570 at one knot; 540 is a 0.6 u stump of that knot.
# The arc starts where 520's ink is clean again above the knot and lands tangentially on 570 where 570 is clean below it.
LB_, ax570 = clean_len(570, False, 0.3); PB0 = at_len(ax570, LB_); uB = axis_tan(570, PB0, +1)
B_, _, hB_b, eB1, eB2 = ink_cut(570, PB0, uB, label='570')
LA_, ax520 = clean_len(520, True, 0.3)
_ax520c = ax520[cum(ax520) <= LA_]                                                   # 520's clean course above the knot
_turnB = 1.0 if cross2(axis_tan(520, _ax520c[-1], -1), uB) > 0 else -1.0; nB = perp(uB) * _turnB
def _f(R_):                                                                           # circle tangent to 570's line at B_, radius R_: gap to 520's course
    c_ = B_ + nB * R_; return float(np.linalg.norm(_ax520c - c_, axis=1).min()) - R_
_lo, _hi = 0.3, 30.0
print(f'  bun: B_ ({B_[0]:.2f},{B_[1]:.2f}) uB {uB.round(3)} nB {nB.round(3)} 520 clean end {_ax520c[-1].round(2)} tan {axis_tan(520, _ax520c[-1], -1).round(3)} turnB {_turnB}; f: ' + ' '.join(f'{R_:.1f}:{_f(R_):+.2f}' for R_ in (0.3, 1, 2, 3, 5, 10, 30)))
assert _f(_lo) > 0 and _f(_hi) < 0, (_f(_lo), _f(_hi))
for _ in range(60):
    _m = (_lo + _hi) / 2
    if _f(_m) > 0: _lo = _m
    else: _hi = _m
R_bun = (_lo + _hi) / 2; cen = B_ + nB * R_bun; A_ = _ax520c[int(np.argmin(np.linalg.norm(_ax520c - cen, axis=1)))]
for _it in range(3):                                                                  # now exactly tangent at 520's own mid A_; 570's tangent point moves a hair along its line
    uA = axis_tan(520, A_, -1); A_, _, hA_b, eA1, eA2 = ink_cut(520, A_, uA, label='520')
    turn = 1.0 if cross2(uA, uB) > 0 else -1.0; nA = perp(uA) * turn; nB = perp(uB) * turn
    R_bun = float(np.dot(A_ - B_, nB)) / (1.0 - float(np.dot(nA, nB))); cen = A_ + nA * R_bun; B_new = cen - nB * R_bun
    B_, _, hB_b, eB1, eB2 = ink_cut(570, B_new, uB, label='570')
    print(f'  bun arc it{_it}: R {R_bun:.3f}, B moved {np.linalg.norm(B_new - B_):.4f} lateral, A ({A_[0]:.3f},{A_[1]:.3f}) B ({B_[0]:.3f},{B_[1]:.3f})')
cen = A_ + nA * R_bun
print(f'  bun arc: tangent point residual {np.linalg.norm((B_ + nB * R_bun) - cen):.4f}')
a0, a1 = math.atan2(A_[1] - cen[1], A_[0] - cen[0]), math.atan2(B_[1] - cen[1], B_[0] - cen[0])
if turn > 0 and a1 < a0: a1 += 2 * math.pi
if turn < 0 and a1 > a0: a1 -= 2 * math.pi
_ang = np.linspace(a0, a1, 96)
bun_pts = np.stack([cen[0] + R_bun * np.cos(_ang), cen[1] + R_bun * np.sin(_ang)], 1)
_rad = np.stack([np.cos(_ang), np.sin(_ang)], 1)
_hw = np.linspace(hA_b, hB_b, len(bun_pts))[:, None]
bun_poly = Polygon(np.vstack([bun_pts + _rad * _hw, (bun_pts - _rad * _hw)[::-1]])).buffer(0)
keep[520] = halfplane(A_, -uA); keep[570] = halfplane(B_, uB)
theta = math.acos(np.clip(np.dot(uA, uB), -1, 1)); cj = (A_ + B_) / 2
ax560 = np.asarray(O[560]['axis'], float); r560 = r_true(560); d560 = cum(ax560)
_s_exit = None
for k in range(len(ax560)):
    tk = axis_tan(560, ax560[k], 0, 0.3); nk = perp(tk)
    if not bun_poly.contains(LineString([ax560[k] + nk * r560[k], ax560[k] - nk * r560[k]])):
        _s_exit = d560[k]; break
_s_cut = max((_s_exit if _s_exit is not None else 0.0) - 0.1, 0.0); P560 = at_len(ax560, _s_cut); u560 = axis_tan(560, P560, +1)
keep[560] = halfplane(P560, u560)
notes.append(f'head: Step 1 hair; the spear shaft, blade and the junction knot removed; bun outline closed with one arc R {R_bun:.2f}, tangent to 520 at ({A_[0]:.2f},{A_[1]:.2f}) and to 570 at ({B_[0]:.2f},{B_[1]:.2f}), bend {math.degrees(theta):.0f} deg; width {2 * hA_b:.2f} to {2 * hB_b:.2f}, both the ink\'s own at the cuts')
print(notes[-1])

# ---------------- forehead: the hair outline (790) turns into the face line (1290); the old arm (1450) used to cover the corner ----------------
q790 = spts(O[790]); q1290 = spts(O[1290])
_A = q790[arc_index(q790, length(q790) - 2.5):]; _B = q1290[arc_index(q1290, length(q1290) - 2.5):][::-1]
_, c_fh0 = mitre_join(_A, _B, 0.7, 0.6)
L790 = length(q790) - float(np.linalg.norm(q790[-1] - c_fh0)) - 0.45; PA_fh = at_len(q790, L790); u790 = axis_tan(790, PA_fh, -1)
L1290 = length(q1290) - float(np.linalg.norm(q1290[-1] - c_fh0)) - 0.45; PB_fh = at_len(q1290, L1290); u1290 = -axis_tan(1290, PB_fh, -1)   # leaving the corner down the face line
if np.dot(u790, c_fh0 - PA_fh) < 0: u790 = -u790
if np.dot(u1290, PB_fh - c_fh0) < 0: u1290 = -u1290
mA_fh, _, hA_fh, _, _ = ink_cut(790, PA_fh, u790, label='790')
mB_fh, _, hB_fh, _, _ = ink_cut(1290, PB_fh, u1290, label='1290')
st_ = ray_hit(mA_fh, u790, mB_fh, u1290); c_fh = mA_fh + u790 * st_[0]
fh_poly, _, _ = mitre_straight(c_fh, u790, u1290, hA_fh, hB_fh, LA=float(np.linalg.norm(c_fh - mA_fh)), LB=float(np.linalg.norm(c_fh - mB_fh)))
keep[790] = halfplane(mA_fh, -u790); keep[1290] = halfplane(mB_fh, u1290)
notes.append(f'forehead: hair outline and face line mitred at ({c_fh[0]:.2f},{c_fh[1]:.2f}), bend {math.degrees(math.acos(np.clip(np.dot(u790, u1290), -1, 1))):.0f} deg')

# ---------------- J5 fold: the contour (160) turns up a short upright at its own weight and into the light chest line (330); both corners mitred ----------------
q160 = spts(O[160]); q330 = spts(O[330])
ax320 = np.asarray(O[320]['axis'], float); d320 = cum(ax320)
u_up = unit(at_len(ax320, d320[-1] - 1.1) - at_len(ax320, d320[-1] - 0.25))         # the old upright's own course (its lower, straight part), pointing up
L160c, ax160 = clean_len(160, True, 0.35); P160 = at_len(ax160, L160c); u160 = axis_tan(160, P160, -1)
if np.dot(u160, ax160[-1] - P160) < 0: u160 = -u160
m160, _, h160, _, _ = ink_cut(160, P160, u160, label='160')
L330c, ax330 = clean_len(330, False, 0.35); P330s = at_len(ax330, L330c); u330s = axis_tan(330, P330s, +1)
if np.dot(u330s, ax330[-1] - P330s) < 0: u330s = -u330s
m330s, _, h330s, _, _ = ink_cut(330, P330s, u330s, label='330 start')
_st = ray_hit(m160, u160, ax320[-1], u_up); c_k = m160 + u160 * _st[0]
_st = ray_hit(c_k, u_up, m330s, u330s); c5 = c_k + u_up * _st[0]
L_up = float(np.linalg.norm(c5 - c_k))
knot_poly, _, _ = mitre_straight(c_k, u160, u_up, h160, h160, LA=float(np.linalg.norm(c_k - m160)), LB=L_up)
fold_poly, _, _ = mitre_straight(c5, u_up, u330s, h160, h330s, LA=L_up, LB=float(np.linalg.norm(c5 - m330s)))
keep[160] = halfplane(m160, -u160); keep[330] = halfplane(m330s, u330s)
_knotU = unary_union([raw_band(s_) for s_ in (310, 440, 480)])                        # the breast point: Step 1's knot of three stubs, one round dot here
dot_c = np.asarray(_knotU.centroid.coords[0]); dot_r = math.sqrt(_knotU.area / math.pi)
dot_poly = Point(dot_c).buffer(dot_r, quad_segs=32)
notes.append(f'breast point: one dot of radius {dot_r:.2f} at ({dot_c[0]:.2f},{dot_c[1]:.2f}) in place of the three junction stubs')
notes.append(f'J5 fold: contour corner at ({c_k[0]:.2f},{c_k[1]:.2f}) bend {math.degrees(math.acos(np.clip(np.dot(u160, u_up), -1, 1))):.0f} deg, upright {L_up:.2f} u at the contour\'s weight ({2 * h160:.2f}), fold corner at ({c5[0]:.2f},{c5[1]:.2f}) bend {math.degrees(math.acos(np.clip(np.dot(u_up, u330s), -1, 1))):.0f} deg, chest line {2 * h330s:.2f} from {np.linalg.norm(c5 - m330s):.2f} u up')

# ---------------- J2 left armpit: the chest line (330) runs on into the new arm line, whose end face is flush with the chest line's outer edge ----------------
r11 = pieces['r11']; q11 = r11['prims'][0][1]; u11 = unit(q11[-1] - q11[arc_index(q11, length(q11) - 1.5)])
q11 = np.vstack([q11[:arc_index(q11, length(q11) - 0.5)], [q11[-1] + u11 * 3.0]])
L330e, _ = clean_len(330, True, 0.3); P330e = at_len(ax330, L330e); u330 = axis_tan(330, P330e, -1)
if np.dot(u330, ax330[-1] - P330e) < 0: u330 = -u330
m330e, _, h330e, _, _ = ink_cut(330, P330e, u330, label='330 end')
X = LineString([m330e - u330 * 0.5, m330e + u330 * 8.0]).intersection(LineString(q11))
P2 = np.asarray(min(getattr(X, 'geoms', [X]), key=lambda g_: g_.distance(Point(m330e))).coords[0])
n330 = perp(u330) * (1.0 if np.dot(perp(u330), u11) > 0 else -1.0)
r11['prims'] = [('poly', q11)]; r11['cuts'] = [halfplane(m330e + n330 * h330e, n330)]
chest_ext = buf(np.array([m330e, P2]), h330e)
keep[330] = keep[330].intersection(halfplane(m330e, -u330))
notes.append(f'J2: chest line 330 runs on {np.linalg.norm(P2 - m330e):.2f} u from ({m330e[0]:.2f},{m330e[1]:.2f}) into the arm line at ({P2[0]:.2f},{P2[1]:.2f}); the arm line ends flush with the chest line\'s outer edge')

# ---------------- J4 right armpit: the contour's own two edges carried straight on to the corner; the arm's lower line from the corner; one bush of three thorns ----------------
q1180 = spts(O[1180]); L1180 = length(q1180)
L1180c, ax1180 = clean_len(1180, True, 0.4); P1180 = at_len(ax1180, L1180c); t_base = axis_tan(1180, P1180, -1)
if np.dot(t_base, ax1180[-1] - P1180) < 0: t_base = -t_base
n_out = perp(t_base) * (-1.0 if perp(t_base)[0] > 0 else 1.0)
base1180, _, h1180, cut_o, cut_i = ink_cut(1180, P1180, t_base, nref=n_out, label='1180')
do_, di_ = t_base, t_base; mo, mi = cut_o, cut_i
w1180 = 2 * h1180
def on_line(m, d, p0, n0):
    st = ray_hit(m, d, p0, perp(n0)); return m + d * st[0]
keep[1180] = halfplane(base1180, -t_base)
r93 = pieces['r93']; q93 = r93['prims'][0][1]; h93 = r93['width'] / 2
u93, m93 = fit_dir(q93, 1.5, 4.5)
u_mid = t_base; m_mid = base1180
st = ray_hit(m_mid, u_mid, m93, u93); C4 = m_mid + u_mid * st[0]; u1180 = u_mid
arm_raw = np.vstack([[C4 + u93 * k for k in np.arange(0.0, 4.01, 0.5)], q93[(q93 - C4) @ u93 > 4.4]])
arm_sm = smooth_open(arm_raw, 0.05)
_k1 = int(np.argmax((arm_sm - C4) @ u93 > 2.8)); P1 = arm_sm[_k1]; u93 = unit(arm_sm[_k1 + 1] - arm_sm[_k1 - 1])
st = ray_hit(m_mid, u_mid, P1, u93); C4 = m_mid + u_mid * st[0]
arm_sm = np.vstack([[C4], arm_sm[_k1:]])
arm_pts = to_loop(arm_sm, 'handR')
n93o = perp(u93) * (1.0 if np.dot(perp(u93), n_out) > 0 else -1.0)
Po4 = on_line(mo, do_, C4 + n93o * h93, perp(u93)); Pi4 = on_line(mi, di_, C4 - n93o * h93, perp(u93))
corner4 = Polygon([cut_o, Po4, C4 + u93 * 0.45 + n93o * h93, C4 + u93 * 0.45 - n93o * h93, Pi4, cut_i]).buffer(0)
arm93 = buf(arm_pts, h93).intersection(halfplane(C4 + u93 * 0.3, u93)).difference(HAND_HOLE['handR'])
contour_ext = Polygon()
r93['prims'] = [('poly', arm_pts)]; bend = arm_pts; del pieces['r93']
n_in = -n93o
s_pi = float(np.dot(Pi4 - C4, u93))
def bush_pt(s_, d_):                                                                 # s along the arm line from the corner, d below its lower edge
    p_ = at_len(arm_sm, s_); t_ = unit(at_len(arm_sm, s_ + 0.2) - at_len(arm_sm, max(s_ - 0.2, 0.0))); return p_ + perp(t_) * np.sign(np.dot(perp(t_), n_in)) * (h93 + d_)
BUSH = [(s_pi, 0.0), (s_pi + 0.6, 3.4), (s_pi + 1.3, 1.0), (s_pi + 2.3, 2.8), (s_pi + 2.7, 0.95), (s_pi + 3.6, 2.3), (s_pi + 4.0, 0.0), (s_pi + 4.0, -0.5), (s_pi, -0.5)]
tris = [Polygon([bush_pt(s_, d_) for s_, d_ in BUSH]).buffer(0)]
assert np.linalg.norm(Po4 - C4) < 2.7, 'outer corner point beyond the straight start of the arm'
notes.append(f'J4: contour edges carried straight to the corner at ({C4[0]:.2f},{C4[1]:.2f}), bend {math.degrees(math.acos(np.clip(np.dot(u1180, u93), -1, 1))):.0f} deg; one bush rising from the crease corner, thorns 3.4 / 2.8 / 2.3 u below the arm line')

# ---------------- J3 throat L ----------------
thr = [r for r in kept if r['id'] == 1300][0]['pts']
r103 = pieces['r103']; q103 = r103['prims'][0][1]
_, c3_0 = mitre_join(thr, q103, 0.7, r103['width'] / 2)
Pc3 = nearest_axis_pt(1300, c3_0 - unit(thr[-1] - thr[arc_index(thr, length(thr) - 0.5)]) * 0.7)
u_thr = axis_tan(1300, Pc3, -1)
if np.dot(u_thr, c3_0 - Pc3) < 0: u_thr = -u_thr
m_thr, _, h_thr, _, _ = ink_cut(1300, Pc3, u_thr, label='1300')
legA = np.array([m_thr, m_thr + u_thr * (float(np.linalg.norm(c3_0 - m_thr)) + 0.5)])
corner3, c3 = mitre_join(legA, q103, h_thr, r103['width'] / 2)
corner3 = corner3.intersection(Point(c3).buffer(3.0))
keep[1300] = halfplane(m_thr, -u_thr)
notes.append(f'J3 throat/shoulder mitre at ({c3[0]:.2f},{c3[1]:.2f}), throat {2 * h_thr:.2f} wide at the cut')

# ---------------- J1 left shoulder: 890 carried straight on past the old join as one band; the arm's top line stops flush inside it ----------------
q890 = spts(O[890])
L890c, ax890 = clean_len(890, False, 0.3); P890 = at_len(ax890, L890c); u_down = -axis_tan(890, P890, +1)     # 890's axis runs upward from its lower end
if u_down[1] < 0: u_down = -u_down
m_neck, _, h_neck, eN1, eN2 = ink_cut(890, P890, u_down, label='890')
P_nbot = m_neck + u_down * ((18.82 - m_neck[1]) / u_down[1])
neck_band = buf(np.array([m_neck, P_nbot]), h_neck)
keep[890] = halfplane(m_neck, -u_down)
r13 = pieces['r13']; q13 = r13['prims'][0][1]; u13 = unit(q13[-1] - q13[arc_index(q13, length(q13) - 1.5)])
q13 = np.vstack([q13[:arc_index(q13, length(q13) - 1.0)], [q13[-1] + u13 * 1.8]])
side = 1.0 if np.dot(perp(u_down), u13) > 0 else -1.0
e_far = eN1 if np.dot(eN1 - m_neck, perp(u_down) * side) > 0 else eN2
r13['prims'] = [('poly', q13)]; r13['cuts'] = [halfplane(e_far, perp(u_down) * side)]
notes.append(f'J1: back-of-neck line carried straight down from ({m_neck[0]:.2f},{m_neck[1]:.2f}) as one {2 * h_neck:.2f} band; the shoulder line runs into it, cut flush with its far edge')

# ---------------- new ink ----------------
def fit_cubic_free(q, iters=6):
    """one cubic through q[0] and q[-1]; the two inner control points by least squares (chord parameters, Newton-refined)"""
    q = np.asarray(q, float); p0, p3 = q[0], q[-1]
    d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]); u = d / d[-1]
    c = None
    for _ in range(iters):
        B1 = 3 * (1 - u) ** 2 * u; B2 = 3 * (1 - u) * u ** 2
        rhs = q - np.outer((1 - u) ** 3, p0) - np.outer(u ** 3, p3)
        sol, *_ = np.linalg.lstsq(np.stack([B1, B2], 1), rhs, rcond=None)
        c = (p0, sol[0], sol[1], p3); u = np.clip(reparam(q, u, c), 0, 1)
    return c
def single_fit(q, tol=0.1):
    """one primitive for a short polyline, end points kept: line if nearly straight, else arc if round, else one cubic"""
    q = np.asarray(q, float)
    res = _line_resid(q)
    if res.max() <= tol:
        return primitive_points('line', [q[0], q[-1]], 0.03), 'line'
    cx, cy, r, dev = fit_circle(q)
    if dev <= tol and r < 20:
        a0, a1 = arc_from_pts(q, cx, cy, r)
        return arc_points(cx, cy, r, a0, a1, 0.03), 'arc'
    c = fit_cubic_free(q)
    return primitive_points('bezier', [c], 0.03), 'cubic'
for pid, p in pieces.items():
    p['var_w'] = False
# the forearm lines end on the loop centreline: run them on and cut them with the loop's hole
for pid, hand in (('r11', 'handL'), ('r13', 'handL'), ('r103', 'handR')):
    p = pieces[pid]; q = p['prims'][0][1]; ring = HAND_RING[hand]
    if ring.distance(Point(q[0])) < 0.3 and ring.distance(Point(q[-1])) > 0.3:
        q = q[::-1]
    if ring.distance(Point(q[-1])) < 0.3:
        p['prims'] = [('poly', to_loop(q, hand))]; p.setdefault('cuts', []).append(HAND_HOLE[hand])
for hand in ('handL', 'handR'):
    arms_raw = [pieces[k]['prims'][0][1] for k in pieces if k.startswith(hand + '_in')]
    oval = LineString(pieces[hand + '_oval']['prims'][0][1])
    # the junction: the end of each arm that is not on the loop; the longest arm passes through it (split it there)
    inner_ends = []
    for q in arms_raw:
        for e in (q[0], q[-1]):
            if oval.distance(Point(e)) > 0.6:
                inner_ends.append(e)
    Jpt = np.mean(inner_ends, axis=0)
    arms = []
    for q in arms_raw:
        d = np.linalg.norm(q - Jpt, axis=1); i = int(np.argmin(d))
        if 0.3 < length(q[:i + 1]) and 0.3 < length(q[i:]) and oval.distance(Point(q[0])) < 0.6 and oval.distance(Point(q[-1])) < 0.6:
            arms += [q[:i + 1][::-1], q[i:]]                                    # passes through the junction: two arms
        else:
            arms.append(q if d[0] < d[-1] else q[::-1])                            # junction end first
    for k in [k for k in pieces if k.startswith(hand + '_in')]:
        del pieces[k]
    dirs = []
    HF = pickle.load(open(os.path.join(HERE, 'handfit.pkl'), 'rb'))[hand] if os.path.exists(os.path.join(HERE, 'handfit.pkl')) else None
    if HF is not None:                                                             # junction + four single cubics fitted to the raster by IoU (handfit.py)
        Jpt = HF['J']; fitted = [(np.asarray(c), 'cubic (raster fit)') for c in HF['curves']]
        notes.append(f'{hand}: one junction at ({Jpt[0]:.2f},{Jpt[1]:.2f}), four single cubics fitted to the raster (inner IoU {HF["iou"]:.2f})')
    else:
        fitted = []
        for q in arms:
            if np.linalg.norm(q[0] - Jpt) < 1.2:                                  # the skeleton bends into the junction: start the arm at the junction itself
                q = np.vstack([[Jpt], q[arc_index(q, 0.5):]])
            pts_, kind = single_fit(q); fitted.append((to_loop(pts_, hand), kind))
    for n_, (pts_, kind) in enumerate(fitted):
        dirs.append(math.degrees(math.atan2(*unit(pts_[min(6, len(pts_) - 1)] - pts_[0])[::-1])) % 360)
        pieces[f'{hand}_arm{n_}'] = dict(id=f'{hand}_arm{n_}', prims=[('poly', pts_)], width=0.9, caps=[None, None], free=[False, False], var_w=False, taper=None, len=round(length(pts_), 2), clip=HAND_FILL[hand])
        notes.append(f'{hand} arm {n_}: {kind}, {length(pts_):.1f} u')
    dirs = sorted(dirs); gaps = [(dirs[(i + 1) % len(dirs)] - dirs[i]) % 360 for i in range(len(dirs))]
    print(f'{hand}: junction ({Jpt[0]:.2f},{Jpt[1]:.2f}), arm directions {np.round(dirs, 0)}, gaps {np.round(gaps, 0)} (all < 180 -> the square caps bury each other)')
def render_piece(p):
    g = MAT.render_clean(None, [p], [])
    for hp in p.get('cuts', []): g = g.difference(hp)
    if p.get('clip') is not None: g = g.intersection(p['clip'])
    return g
piece_poly = {pid: render_piece(p) for pid, p in pieces.items()}


new_ink = unary_union(list(piece_poly.values()) + list(tris) + [contour_ext, arm93, corner4, corner3, bun_poly, chest_ext, knot_poly, fold_poly, fh_poly, neck_band, dot_poly]).buffer(0)

# ---------------- body: Step 1 ink minus the removed strokes; kept strokes protected by their exact bands ----------------
import time; _t0 = time.time()
Zcore = {sid: unary_union([G[sid], raw_band(sid)]).difference(kp) for sid, kp in keep.items()}       # the parts of cut strokes beyond their cuts
D = DZ_core.buffer(0.06)
Z = {sid: Zcore[sid].buffer(0.15).difference(keep[sid]) for sid in keep}
near = unary_union([D, KNOTS, unary_union(list(Z.values()))]).buffer(0.8)
prot = {}; bulges = []; fragments = []
for sid in kept_ids:
    if not G[sid].intersects(near) and not raw_band(sid).intersects(near): continue
    E, bulge, whole = band_info(sid)
    if whole: fragments.append(sid); bulges.append(bulge); continue
    if not bulge.is_empty: bulges.append(bulge)
    g = E.buffer(0.012)
    st = O[sid]; hw = st['width'] / 2; q_ = spts(st)
    for end, pt in ((0, q_[0]), (1, q_[-1])):                                               # genuine flat caps: the classified polygon's end
        c = st['caps'][end]
        if c is not None and max(c) <= hw + 0.35: g = g.union(G[sid].intersection(Point(pt).buffer(hw * 1.2)))
    if sid in keep: g = g.intersection(keep[sid])
    prot[sid] = g
P_all = unary_union(list(prot.values()) + list(fills1))
T = unary_union(bulges) if bulges else Polygon()
print(f'protection: {len(prot)} exact bands, {len(bulges)} bulges trimmed, junction fragments {sorted(fragments)}  ({time.time() - _t0:.0f}s)')
body = ink1.difference(unary_union([D, T]).difference(P_all))
for sid in keep:
    others = unary_union([prot[s] for s in prot if s != sid and prot[s].intersects(Z[sid])] + list(fills1))
    body = body.difference(Z[sid].difference(others))
body = body.difference(box(19.0, 20.3, 23.5, 26.3).difference(P_all)).difference(C2['arrow'])                 # the spear head is a filled shape of its own                     # the old arm's underside at the left armpit: ink no kept stroke claims
cand = body.intersection(unary_union([zone_all, T]).buffer(0.45)).difference(unary_union([P_all, new_ink.buffer(0.02)]))
bumps = [q for q in getattr(cand, 'geoms', [cand]) if q.area > 0.0005]
if bumps: body = body.difference(unary_union(bumps).buffer(0.001))
parts = list(getattr(body, 'geoms', [body])); good = [q for q in parts if q.area >= 0.3]
print(f'body: {len(parts)} parts, {len(parts) - len(good)} crumbs removed ({sum(q.area for q in parts if q not in good):.2f} u²); {len(bumps)} residue bits shaved ({sum(b.area for b in bumps):.3f} u²)  ({time.time() - _t0:.0f}s)')
body = unary_union(good)
comp = unary_union([body, new_ink]).buffer(0)
for n_ in notes: print('  ' + n_)
b_ = box(0, 34, 60, 70); print('tail+fin vs Step 1 IoU %.4f' % compare(ink1.intersection(b_), comp.intersection(b_))['iou'])
b_ = box(19, 3, 33, 15); print('head vs Step 1 (spear removed) IoU %.4f' % compare(ink1.difference(G[130]).difference(G[150]).difference(G[510]).intersection(b_), comp.intersection(b_))['iou'])
def fill_slivers(g, amax=0.02):
    out = []; n = 0
    for p in getattr(g, 'geoms', [g]):
        keep_ = [h for h in p.interiors if Polygon(h).area >= amax]; n += len(p.interiors) - len(keep_)
        out.append(Polygon(p.exterior, keep_))
    return unary_union(out), n
comp = comp.buffer(0.0005, join_style='mitre').buffer(-0.0005, join_style='mitre')
comp, n_sliv = fill_slivers(comp, 0.05); print(f'slivers filled: {n_sliv}')
pickle.dump(dict(comp=comp, body=body, new_ink=new_ink, pieces=list(pieces.values()), tris=tris, contour_ext=contour_ext, arm93=arm93, corner4=corner4, corner3=corner3, bun_poly=bun_poly, chest_ext=chest_ext,
                 bend=bend, c3=c3, C4=C4, cj=cj, P2=P2, c5=c5, c_k=c_k, knot_poly=knot_poly, fold_poly=fold_poly, neck_band=neck_band, dot_poly=dot_poly, fh_poly=fh_poly, ink2v=ink2v, boltv=boltv, notes=notes, kept_ids=kept_ids,
                 keep=keep, A_=A_, B_=B_, bun_pts=bun_pts, P_all=P_all, T=T, D=D, KNOTS=KNOTS, prot=prot, Z=Z, fragments=fragments, Zcore=Zcore), open(os.path.join(HERE, 'final2.pkl'), 'wb'))

b = comp.bounds
vb = (b[0] - 2, b[1] - 2, b[2] - b[0] + 4, b[3] - b[1] + 4); vbs = ' '.join(fmt(v, 2) for v in vb)
W, H = 540, int(540 * vb[3] / vb[2])
def svg(body_): return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vbs}" width="{W}" height="{H}">{body_}</svg>'
INK = '#1c1b19'
open(os.path.join(HERE, 'step2-clean.svg'), 'w').write(svg(f'<path d="{svgpath(comp, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/>'))
parts_ = [f'<path d="{svgpath(comp, prec=3, tol=0.004)}" fill="#dedbd3" fill-rule="evenodd"/>']
for r in kept:
    q = r['pts']
    parts_.append(f'<path d="M{" L".join(f"{fmt(x, 2)} {fmt(y, 2)}" for x, y in q[::3])}" fill="none" stroke="#1c1b19" stroke-width="0.28"/>')
for q in (np.array([m330e, P2]), bun_pts, np.array([base1180, C4]), np.array([c_k - u160 * 2.0, c_k, c5, c5 + u330s * 1.5]), np.array([m_neck, P_nbot]), np.array([mA_fh, c_fh, mB_fh]), np.array([m_thr, c3])) + tuple(p['prims'][0][1] for p in pieces.values()):
    parts_.append(f'<path d="M{" L".join(f"{fmt(x, 2)} {fmt(y, 2)}" for x, y in np.asarray(q)[::3])}" fill="none" stroke="#e3761b" stroke-width="0.28"/>')
for t in tris:
    parts_.append(f'<path d="{svgpath(t, prec=3, tol=0.003)}" fill="none" stroke="#7a3fb8" stroke-width="0.22"/>')
for c_ in (cj, P2, C4, c3, m_neck, c5, c_fh):
    parts_.append(f'<circle cx="{fmt(c_[0], 2)}" cy="{fmt(c_[1], 2)}" r="2.2" fill="none" stroke="#1f8f80" stroke-width="0.18" stroke-dasharray="0.4 0.3"/>')
open(os.path.join(HERE, 'step2-construction.svg'), 'w').write(svg(''.join(parts_)))
print(f'written ({time.time() - _t0:.0f}s)')
