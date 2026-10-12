"""handfit.py: the inner strokes of each hand as ONE junction + four single cubics, fitted to the raster by IoU.
Keeps the construction clean (uniform 0.9 pen, square caps burying each other at the junction, ends run into the loop)
and lets the raster decide the junction point, the start tangents and the end points. -> handfit.pkl"""
import os, sys, math, pickle, numpy as np
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import unary_union
from scipy.optimize import minimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from geom import *
from mat import MAT
HERE = os.path.dirname(os.path.abspath(__file__))
C4 = pickle.load(open(os.path.join(HERE, 'compose4.pkl'), 'rb'))
pieces = {p['id']: p for p in C4['pieces']}; ink2v = C4['ink2v']
HW = 0.45
def unit(v): return v / (np.linalg.norm(v) + 1e-12)
def arc_index(pts, L):
    d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]); return int(np.searchsorted(d, L))
def strip(pts, hw=HW): return LineString(pts).buffer(hw, cap_style='flat', join_style='mitre', mitre_limit=6, quad_segs=16)
out = {}
for hand in ('handL', 'handR'):
    oval_pts = pieces[hand + '_oval']['prims'][0][1]; oval = LineString(oval_pts)
    ring = MAT.render_clean(None, [pieces[hand + '_oval']], [])
    fill = Polygon(ring.exterior); hole = Polygon(max(ring.interiors, key=lambda r: Polygon(r).area))
    target = ink2v.intersection(hole)
    arms_raw = [pieces[k]['prims'][0][1] for k in pieces if k.startswith(hand + '_in')]
    inner_ends = [e for q in arms_raw for e in (q[0], q[-1]) if oval.distance(Point(e)) > 0.6]
    Jpt = np.mean(inner_ends, axis=0)
    arms = []
    for q in arms_raw:
        d = np.linalg.norm(q - Jpt, axis=1); i = int(np.argmin(d))
        if 0.3 < length(q[:i + 1]) and 0.3 < length(q[i:]) and oval.distance(Point(q[0])) < 0.6 and oval.distance(Point(q[-1])) < 0.6:
            arms += [q[:i + 1][::-1], q[i:]]
        else:
            arms.append(q if d[0] < d[-1] else q[::-1])
    # initial parameters per arm: start direction (from the junction to the polyline 1 u out), end point on the loop, end tangent
    L_oval = oval.length
    init = []; ends = []
    for q in arms:
        k = arc_index(q, 1.0); th0 = math.atan2(*(q[min(k, len(q) - 1)] - Jpt)[::-1])
        E = q[-1]; s_E = oval.project(Point(E))
        t1 = -end_tangent(q, False, 0.7); th1 = math.atan2(t1[1], t1[0])
        Lc = np.linalg.norm(E - Jpt)
        init += [th0, Lc / 3, Lc / 3, 0.0, th1]
        ends.append(s_E)
    x0 = np.array([0.0, 0.0] + init)
    prev = pickle.load(open(os.path.join(HERE, 'handfit.pkl'), 'rb')).get(hand) if os.path.exists(os.path.join(HERE, 'handfit.pkl')) else None
    if prev is not None and len(prev['x']) == len(x0): x0 = np.asarray(prev['x'], float)
    def build(x, run=1.0):
        J = Jpt + x[:2]; polys = []; curves = []
        for a in range(len(arms)):
            th0, a1, a2, ds, th1 = x[2 + 5 * a: 7 + 5 * a]
            E = np.asarray(oval.interpolate((ends[a] + ds) % L_oval).coords[0])
            t0 = np.array([math.cos(th0), math.sin(th0)]); t1 = np.array([math.cos(th1), math.sin(th1)])
            c = (J, J + t0 * a1, E - t1 * a2, E)
            pts = primitive_points('bezier', [c], 0.03)
            pts = np.vstack([pts, [E + t1 * run]])
            curves.append(pts); polys.append(strip(pts))
        return J, curves, polys
    ras_holes = [g for g in getattr(hole.difference(target), 'geoms', [hole.difference(target)]) if g.area > 0.15]   # the white inside the loop, per region
    def hole_score(ours):
        our_holes = list(getattr(hole.difference(ours), 'geoms', [hole.difference(ours)]))
        tot = 0.0
        for rh in ras_holes:
            best = 0.0
            for oh in our_holes:
                i = rh.intersection(oh).area
                if i > 0: best = max(best, i / rh.union(oh).area)
            tot += 1 - best
        return tot / max(len(ras_holes), 1)
    def score(x):
        J, curves, polys = build(x)
        ours = unary_union(polys).intersection(hole)
        inter = ours.intersection(target).area; uni = ours.union(target).area
        hs = hole_score(ours)
        pen = 0.0
        for a in range(len(arms)):                                             # keep the cubics sane: no loops, tangents positive
            th0, a1, a2, ds, th1 = x[2 + 5 * a: 7 + 5 * a]
            Lc = np.linalg.norm(np.asarray(oval.interpolate((ends[a] + ds) % L_oval).coords[0]) - (Jpt + x[:2]))
            pen += max(0.0, a1 + a2 - 1.1 * Lc) ** 2 + max(0.0, 0.15 - a1) ** 2 + max(0.0, 0.15 - a2) ** 2
        return 1 - inter / uni + 0.5 * hs + 2.0 * pen
    s0 = score(x0)
    bounds = [(-1.0, 1.0), (-1.0, 1.0)]
    for a in range(len(arms)):
        th0, a1, a2, ds, th1 = init[5 * a: 5 * a + 5]
        bounds += [(th0 - 1.1, th0 + 1.1), (0.15, 4.0), (0.15, 4.0), (-1.2, 1.2), (th1 - 0.8, th1 + 0.8)]
    res = minimize(score, x0, method='Powell', bounds=bounds, options=dict(maxiter=40, xtol=1e-3, ftol=1e-5, disp=False))
    res = minimize(score, res.x, method='Powell', bounds=bounds, options=dict(maxiter=40, xtol=1e-3, ftol=1e-6, disp=False))
    J, curves, polys = build(res.x)
    ours = unary_union(polys).intersection(hole)
    iou = ours.intersection(target).area / ours.union(target).area
    dirs = sorted(math.degrees(math.atan2(*unit(c[6] - c[0])[::-1])) % 360 for c in curves)
    gaps = [(dirs[(i + 1) % len(dirs)] - dirs[i]) % 360 for i in range(len(dirs))]
    print(f'{hand}: junction ({Jpt[0]:.2f},{Jpt[1]:.2f}) -> ({J[0]:.2f},{J[1]:.2f}); score {s0:.3f} -> {res.fun:.3f}; inner IoU {iou:.3f}; holes {1 - hole_score(ours):.3f}; arm dirs {np.round(dirs)} gaps {np.round(gaps)}')
    for a, c in enumerate(curves):
        th0, a1, a2, ds, th1 = res.x[2 + 5 * a: 7 + 5 * a]
        print(f'   arm {a}: start {math.degrees(th0):.0f} deg, tangents {a1:.2f}/{a2:.2f}, end slid {ds:+.2f} u, length {length(c):.2f}')
    out[hand] = dict(J=J, curves=curves, iou=iou, iou0=1 - s0, x=res.x, fill=fill, hole=hole)
pickle.dump(out, open(os.path.join(HERE, 'handfit.pkl'), 'wb')); print('written handfit.pkl')
