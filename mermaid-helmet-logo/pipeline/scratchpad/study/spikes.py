"""Armpit spikes: the only tapered marks in the drawing. Rebuilt as clean filled triangles hanging from the arm line,
sized from the raster (base width, length, direction). Replaces the traced spike pieces in the composite."""
import os, sys, pickle, numpy as np, shapely
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import unary_union
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from lib import svgpath, fmt
from geom import *
from mat import MAT, compare
HERE = os.path.dirname(os.path.abspath(__file__))
C = pickle.load(open(os.path.join(HERE, 'compose3.pkl'), 'rb'))
pieces, kept_ink, kept_fills, ink2v, boltv = C['pieces'], C['kept_ink'], C['kept_fills'], C['ink2v'], C['boltv']
BOX = box(36.5, 18, 45, 25)
ARM_ID = 'r93'                                   # the lower line of the raised right arm; the spikes hang from it
arm = next(p for p in pieces if p['id'] == ARM_ID)
arm_line = LineString(arm['prims'][0][1])
arm_ink = arm_line.buffer(arm['width'] / 2 + 0.12)
other = unary_union([LineString(p['prims'][0][1]).buffer(p['width'] / 2 + 0.12) for p in pieces if p['id'] != ARM_ID and p['len'] > 6])
reg = ink2v.intersection(BOX).difference(arm_ink).difference(other).difference(kept_ink.buffer(0.1))
blobs = [g for g in getattr(reg, 'geoms', [reg]) if g.area > 0.5]
spikes = []
for g in blobs:
    er = g.buffer(-0.22)
    parts_ = [q for q in getattr(er, 'geoms', [er]) if not q.is_empty and q.area > 0.08]
    gpts = np.asarray(g.exterior.coords)
    owner = np.array([int(np.argmin([q.distance(Point(x)) for x in [pt] for q in parts_])) for pt in gpts]) if parts_ else np.array([])
    print(f'blob area {g.area:.2f} bounds {np.round(g.bounds, 1)} -> {len(parts_)} parts')
    for pi, part in enumerate(parts_):
        own = part.buffer(0.45).intersection(g)
        pts = gpts[owner == pi]
        if len(pts) < 3:
            continue
        d = np.array([arm_line.distance(Point(q)) for q in pts])
        apex = pts[np.argmax(d)]; L = float(d.max())
        contact = own.intersection(arm_ink.boundary.buffer(0.12))   # where this spike meets the arm line
        if contact.is_empty or L - arm['width'] / 2 < 1.0:
            continue
        foot = np.array([contact.centroid.x, contact.centroid.y])
        if apex[1] < foot[1]:                                        # must hang below the arm
            continue
        axis = apex - foot; axis /= np.linalg.norm(axis)
        nrm = np.array([-axis[1], axis[0]])
        band = LineString([foot + axis * 0.45 - nrm * 3, foot + axis * 0.45 + nrm * 3]).intersection(own)
        bw = float(band.length) if not band.is_empty else 1.0
        spikes.append(dict(apex=apex, foot=foot, axis=axis, nrm=nrm, base_w=bw, length=round(float(np.linalg.norm(apex - foot)), 1), area=own.area))
spikes.sort(key=lambda s: s['foot'][0])
for s in spikes:
    print(f"spike: foot ({s['foot'][0]:.2f},{s['foot'][1]:.2f}) apex ({s['apex'][0]:.2f},{s['apex'][1]:.2f}) length {s['length']:.1f} base {s['base_w']:.2f} angle {np.degrees(np.arctan2(s['axis'][1], s['axis'][0])):.0f}°")
# regularise: same base width for all three, tips sharp, base set 0.25 inside the arm line for a clean join
BW = 1.0                                                       # measured 0.9-1.2 at the base; one width for all three
BULGE = 0.12                                                   # the artist's spikes are fat: sides bowed out a little
tris = []
for s in spikes:
    base_c = s['foot'] - s['axis'] * 0.55                          # base tucked into the arm line
    tip = s['foot'] + s['axis'] * s['length']
    a, b_ = base_c + s['nrm'] * BW / 2, base_c - s['nrm'] * BW / 2
    def side(p0, p1, out):
        m = (p0 + p1) / 2 + out * BULGE * np.linalg.norm(p1 - p0)
        t = np.linspace(0, 1, 24)[:, None]
        return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * m + t ** 2 * p1
    tris.append(Polygon(np.vstack([side(a, tip, s['nrm']), side(tip, b_, -s['nrm'])[1:]])).buffer(0))

new_pieces = [p for p in pieces if not (BOX.contains(Point(p['prims'][0][1].mean(axis=0))) and p['len'] < 6)]

# the left arm's lower line meets the body contour in a flush corner: its end face continues the contour's outer edge
def flush_corner(pid):
    p = next(q for q in new_pieces if q['id'] == pid)
    q = p['prims'][0][1]
    end = 1 if kept_ink.distance(Point(q[-1])) < kept_ink.distance(Point(q[0])) else 0
    pt = q[0] if end == 0 else q[-1]
    t = -end_tangent(q, end == 0, 0.6)
    cands = [(np.linalg.norm(sub[i] - pt), sub, i, w) for (kind, sub), w in zip(C['kept_lines'], C['kept_w']) for i in (0, -1)]
    d, sub, i, w = min(cands, key=lambda c: c[0])
    if d > 1.5:
        print(f'flush corner {pid}: no body line end near the arm end ({d:.2f}); left as is')
        return MAT.render_clean(None, [p], [])
    db = end_tangent(sub, i == 0, 1.8)                 # body contour direction (over its top 1.8 u), into the body
    nb = np.array([-db[1], db[0]])
    if np.dot(nb, t) < 0:
        nb = -nb                                       # outward: towards the arm's end
    edge_pt = sub[i] + nb * (w / 2)                    # a point on the contour's outer edge line
    # the arm's lower edge (body side) meets that edge line at X; the arm is cut square through X
    n_arm = np.array([-t[1], t[0]])
    if np.dot(n_arm, db) < 0:
        n_arm = -n_arm                                 # towards the body
    hw = p['width'] / 2
    a0 = pt + n_arm * hw                               # a point on the arm's lower edge
    M = np.array([[t[0], -db[0]], [t[1], -db[1]]])
    lam = np.linalg.solve(M, edge_pt - a0)
    X = a0 + t * lam[0]
    q = np.vstack([pt + t * 3.0, q]) if end == 0 else np.vstack([q, pt + t * 3.0])
    p['prims'] = [('poly', q)]; p['len'] = round(length(q), 2)
    half = Polygon([X + n_arm * 20, X - n_arm * 20, X - n_arm * 20 - t * 20, X + n_arm * 20 - t * 20])
    g = MAT.render_clean(None, [p], []).intersection(half)
    print(f'flush corner {pid}: contour end {np.round(sub[i], 2)} dir {np.round(db, 2)}, corner X {np.round(X, 2)}, arm end moved {np.linalg.norm(X - pt):.2f}')
    return g
new_ink = MAT.render_clean(None, new_pieces, tris)
comp = unary_union([kept_ink, new_ink])
c = compare(ink2v, comp)
print(f'composite v4 vs aligned raster: IoU {c["iou"]:.4f}, raster-only {c["a_only"]:.1f}, composite-only {c["b_only"]:.1f}')
cb = compare(ink2v.intersection(BOX), comp.intersection(BOX)); print(f'  spikes box IoU {cb["iou"]:.3f}')
pickle.dump(dict(pieces=new_pieces, spikes=tris, kept=C['kept'], kept_w=C.get('kept_w'), kept_fills=kept_fills, kept_lines=C['kept_lines'], comp=comp, ink2v=ink2v, boltv=boltv, kept_ink=kept_ink, new_ink=new_ink),
            open(os.path.join(HERE, 'compose4.pkl'), 'wb'))
b = comp.union(ink2v).bounds
vb = (b[0] - 2, b[1] - 2, b[2] - b[0] + 4, b[3] - b[1] + 4); vbs = ' '.join(fmt(v, 2) for v in vb)
W, Hh = 540, int(540 * vb[3] / vb[2])
def svg(body): return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vbs}" width="{W}" height="{Hh}">{body}</svg>'
INK = '#1c1b19'
open(os.path.join(HERE, 'step2-clean.svg'), 'w').write(svg(f'<path d="{svgpath(comp, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/>'))
both = comp.intersection(ink2v)
open(os.path.join(HERE, 'step2-overlay.svg'), 'w').write(svg(
    f'<rect x="{fmt(vb[0])}" y="{fmt(vb[1])}" width="{fmt(vb[2])}" height="{fmt(vb[3])}" fill="#fff"/>'
    f'<path d="{svgpath(ink2v.difference(comp), prec=3, tol=0.003)}" fill="#d6473b" fill-rule="evenodd"/>'
    f'<path d="{svgpath(comp.difference(ink2v), prec=3, tol=0.003)}" fill="#2a9df4" fill-rule="evenodd"/>'
    f'<path d="{svgpath(both, prec=3, tol=0.003)}" fill="#c6c3bc" fill-rule="evenodd"/>'))
parts = [f'<path d="{svgpath(comp, prec=3, tol=0.004)}" fill="#dedbd3" fill-rule="evenodd"/>']
for kind, q in C['kept_lines']:
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    parts.append(f'<path d="{pd}" fill="none" stroke="#1c1b19" stroke-width="0.28"/>')
for p in new_pieces:
    q = p['prims'][0][1]
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    parts.append(f'<path d="{pd}" fill="none" stroke="#e3761b" stroke-width="0.28"/>')
for t in tris:
    parts.append(f'<path d="{svgpath(t, prec=3, tol=0.002)}" fill="none" stroke="#7a3fb8" stroke-width="0.22"/>')
open(os.path.join(HERE, 'step2-construction.svg'), 'w').write(svg(''.join(parts)))
print('written')
