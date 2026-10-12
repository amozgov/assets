"""Step 2 re-pose, v2: the Step-1 strokes (per-stroke clean polygons from vector.py) are kept wherever the lightning
raster confirms them; the raster's new strokes (arms, hands, armpit spikes) are transplanted in the vector frame."""
import os, sys, pickle, numpy as np, shapely
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import unary_union
from shapely import affinity
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from lib import svgpath, fmt
from geom import *
from mat import MAT, compare
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = pickle.load(open(os.path.join(HERE, 'step1.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'prep.pkl'), 'rb'))
A = pickle.load(open(os.path.join(HERE, 'align2.pkl'), 'rb'))
strokes1, geoms1, fills1, ink1 = S1['strokes'], S1['geoms'], S1['fills'], S1['ink']
strokes2, ink2, bolt = P['strokes2'], P['ink2'], P['bolt']
MAT_ = A['mat']; K = A['k']
TOL = float(os.environ.get('ZONE_TOL', 0.5)); EXT = float(os.environ.get('CUT_EXT', 0.4))

def r2v(pts):
    pts = np.asarray(pts, float)
    a, b, c, d, xo, yo = MAT_
    return np.stack([a * pts[:, 0] + b * pts[:, 1] + xo, c * pts[:, 0] + d * pts[:, 1] + yo], 1)
ink2v = affinity.affine_transform(ink2, MAT_)
boltv = affinity.affine_transform(bolt, MAT_)

def stroke_pts(st, xf=None):
    pts = []
    for k, dd in st['prims']:
        q = primitive_points(k, dd, 0.02)
        pts.extend(q if not pts else q[1:])
    pts = np.array(pts)
    return xf(pts) if xf else pts

def runs(pts, zone, keep_inside, ext=EXT, min_len=0.6):
    inside = shapely.contains_xy(zone, pts[:, 0], pts[:, 1])
    want = inside if keep_inside else ~inside
    n = len(pts); out = []; i = 0
    while i < n:
        if not want[i]:
            i += 1; continue
        j = i
        while j + 1 < n and want[j + 1]:
            j += 1
        sub = pts[i:j + 1]
        if len(sub) >= 2 and length(sub) >= min_len:
            c0, c1 = i > 0, j < n - 1
            if c0:
                t = end_tangent(sub, True, 0.5); sub = np.vstack([sub[0] - t * ext, sub])
            if c1:
                t = end_tangent(sub, False, 0.5); sub = np.vstack([sub, sub[-1] - t * ext])
            out.append((sub, c0, c1, want[i:j + 1].sum() / max(want.sum(), 1)))
        i = j + 1
    return out, float(want.mean())

# pass 1: vector strokes confirmed by the raster. A stroke is old (dropped) when a run of it >= 2 u lies where the
# raster has no ink; a stroke that is partly shared keeps only its long confirmed runs (the shoulder lines).
MIS = float(os.environ.get('MISALIGN', 0.3))   # strict: a centreline is confirmed where it runs inside the other drawing's line, allowing this much misalignment
NEAR = float(os.environ.get('NEAR', 0.6))      # lenient: ink of the other drawing within this distance counts as "the same line, slightly moved"
_zones = {}
def zone_for(poly, width, key):
    k = (key, round(width, 2))
    if k not in _zones:
        _zones[k] = poly.buffer(-max(width / 2 - MIS, 0.05))
    return _zones[k]
near_zone = ink2v.buffer(NEAR)
kept, kept_lines, dropped = [], [], []
dropped_polys, kept_w, kept_recs = [], [], []
short_items = []
for st, g in zip(strokes1, geoms1):
    pts = stroke_pts(st)
    L = length(pts)
    rr, frac = runs(pts, zone_for(ink2v, st['width'], 'r'), True, ext=0.1)      # frac = strictly confirmed share
    if L < 3.0:
        short_items.append((st, g, pts, frac)); continue               # decided after the long strokes
    old_runs, _ = runs(pts, near_zone, False, ext=0.0, min_len=0.0)
    if max([length(u[0]) for u in old_runs], default=0.0) < 1.5 and frac >= 0.3:
        # an unconfirmed run at either end (the old shoulder stub) is cut off; short gaps in the middle are tolerated
        inside = shapely.contains_xy(near_zone, pts[:, 0], pts[:, 1])
        i0 = int(np.argmax(inside)); i1 = len(pts) - 1 - int(np.argmax(inside[::-1]))
        head, tail = length(pts[:i0 + 1]) if i0 > 0 else 0.0, length(pts[i1:]) if i1 < len(pts) - 1 else 0.0
        if head > 0.4 or tail > 0.4:
            sub = pts[i0:i1 + 1]
            if head > 0.4:
                t = end_tangent(sub, True, 0.5); sub = np.vstack([sub[0] - t * 0.25, sub])
            if tail > 0.4:
                t = end_tangent(sub, False, 0.5); sub = np.vstack([sub, sub[-1] - t * 0.25])
            region = LineString(sub).buffer(st['width'] / 2 + 0.1, cap_style='flat')
            kept.append(g.intersection(region)); kept_lines.append(('part', sub)); kept_w.append(st['width']); dropped_polys.append(g.difference(region))
            kept_recs.append(dict(id=st['id'], st=st, g=g.intersection(region), pts=sub, kind='part', cut=[head > 0.4, tail > 0.4]))
        else:
            kept.append(g); kept_lines.append(('whole', pts)); kept_w.append(st['width'])
            kept_recs.append(dict(id=st['id'], st=st, g=g, pts=pts, kind='whole', cut=[False, False]))
        continue
    good = [r for r in rr if length(r[0]) >= 3.0 or length(r[0]) / L >= 0.5]
    if not good:
        dropped.append(st['id']); dropped_polys.append(g)
    for sub, c0, c1, _ in good:
        region = LineString(sub).buffer(st['width'] / 2 + 0.1, cap_style='flat')
        kept.append(g.intersection(region)); kept_lines.append(('part', sub)); kept_w.append(st['width']); dropped_polys.append(g.difference(region))
        kept_recs.append(dict(id=st['id'], st=st, g=g.intersection(region), pts=sub, kind='part', cut=[bool(c0), bool(c1)]))
# short strokes: junction connectors of the body touch kept long strokes; leftovers of the old hands do not
long_ink = unary_union(kept).buffer(0.5)
# provisional new ink: raster strokes that run mostly outside the original's ink (the new arms, hands, thorns)
_prov = []
_orig06 = ink1.buffer(0.6)
for st2 in strokes2:
    q2 = stroke_pts(st2, r2v)
    outs, _ = runs(q2, _orig06, False, ext=0.0, min_len=0.0)
    if max([length(u[0]) for u in outs], default=0.0) >= 3.0:        # a raster stroke with a long stretch the original never had: new
        _prov.append(LineString(q2).buffer(st2['width'] * K / 2 + 0.1, cap_style='flat'))
new_prov = unary_union(_prov) if _prov else Polygon()
for st, g, pts, frac in short_items:
    core = LineString(pts).buffer(st['width'] / 2, cap_style='flat')
    rest = core.difference(new_prov)                                  # the part of this stroke not hidden under a new line
    hidden = rest.area < 0.2 * core.area
    cover = rest.intersection(ink2v).area / max(rest.area, 1e-9)      # does the raster draw this stroke where it shows?
    body_region = pts[:, 1].min() > 30.0                              # torso, tail and fin: identical in both drawings, never dropped
    keep_it = body_region or hidden or cover >= 0.55
    print(f"  short {st['id']:5d} at ({pts[:, 0].mean():5.1f},{pts[:, 1].mean():5.1f}) cover {cover:.2f}{' hidden' if hidden else ''} -> {'keep' if keep_it else 'DROP'}")
    if keep_it:
        kept.append(g); kept_lines.append(('whole', pts)); kept_w.append(st['width'])
        kept_recs.append(dict(id=st['id'], st=st, g=g, pts=pts, kind='whole', cut=[False, False]))
    else:
        dropped.append(st['id']); dropped_polys.append(g)
# The kept body is the ORIGINAL outline itself (exact), minus the dropped strokes; the notches that subtracting a
# dropped stroke leaves in a kept one at their junction are patched with the kept stroke's clean polygon there.
U_k = unary_union(kept)
U_d = unary_union(dropped_polys) if dropped_polys else Polygon()
arrow = ink1.intersection(Polygon([[0, 16.04], [9.225, 10.155], [7.923, 13.166], [8.322, 14.543], [10.917, 16.543]]).buffer(0.3))   # spear head (a fill)
U_d02 = U_d.buffer(0.02)
fillets = ink1.intersection(U_d.buffer(0.5)).difference(U_k.buffer(0.1)).difference(U_d02)   # junction fillets of old strokes
base = ink1.difference(U_d02).difference(fillets).difference(arrow)
patch = U_k.intersection(U_d.buffer(1.1))     # also refills the knockout gaps the old arms cut into body lines (gap = 1/2 to 3/4 line)
kept_ink = unary_union([base, patch]).buffer(0)
parts_ = [q for q in getattr(kept_ink, 'geoms', [kept_ink]) if q.area > 0.6 and q.intersection(ink2v).area / q.area > 0.3]
print(f'kept ink: {len(parts_)} parts (dropped {len(getattr(kept_ink, "geoms", [kept_ink])) - len(parts_)} stray bits), area {sum(q.area for q in parts_):.1f} of original {ink1.area:.1f}')
kept_ink = unary_union(parts_)
kept_fills = []
print('dropped ids/len:', [(st['id'], st['len']) for st in strokes1 if st['id'] in dropped])
print('dropped in the body (y>35):', [(st['id'], st['len'], np.round(stroke_pts(st).mean(axis=0), 1).tolist()) for st in strokes1 if st['id'] in dropped and stroke_pts(st)[:, 1].mean() > 35])
print(f'vector: {sum(1 for k, _ in kept_lines if k == "whole")} strokes kept whole, {sum(1 for k, _ in kept_lines if k == "part")} clipped, {len(dropped)} dropped; fills kept {len(kept_fills)}/{len(fills1)}')

# pass 2: raster strokes where the kept vector has no ink
pieces = []
far_zone = kept_ink.buffer(0.8)
reach_zone = kept_ink.buffer(-0.6)   # ~ the centreline of a standard line
ARM_REGION = float(os.environ.get('ARM_Y', 30.0))          # new content lives above this y (vector units)
for st in strokes2:
    pts = stroke_pts(st, r2v)
    if np.mean(pts[:, 1] < ARM_REGION) < 0.5:
        continue
    new_runs, _ = runs(pts, far_zone, False, ext=0.0, min_len=0.0)
    if max([length(u[0]) for u in new_runs], default=0.0) < 1.5:
        continue                                   # near-duplicate of kept ink (or a junction bit): not new content
    covered = zone_for(kept_ink, st['width'] * K, 'k')
    rr, frac = runs(pts, covered, False)
    for sub, c0, c1, _ in rr:
        t0 = c0 or (not st['free'][0] and np.allclose(sub[0], pts[0]) and kept_ink.distance(Point(sub[0])) < 0.9)
        t1 = c1 or (not st['free'][1] and np.allclose(sub[-1], pts[-1]) and kept_ink.distance(Point(sub[-1])) < 0.9)
        bridge = t0 and t1 and length(sub) < 3.5
        if not bridge and not any(not far_zone.contains(Point(q)) for q in sub[::5]):
            continue                               # this run never leaves the neighbourhood of kept ink
        # a cut end is pushed on along its tangent until it sits in the middle of the kept line it joins
        for end in (0, 1):
            if (st['free'][end] and not (c0 if end == 0 else c1)):
                continue                           # a real free end keeps its cap; junction and cut ends land on the kept line's centre
            if end == 1:
                sub = sub[::-1]
            # trim: if the end already runs into a kept line, stop at the first point that reaches its centre zone
            ins = shapely.contains_xy(kept_ink, sub[:, 0], sub[:, 1])
            k = 0
            while k < len(sub) - 2 and ins[k]:
                k += 1
            if k > 0:
                j = next((i for i in range(k - 1, -1, -1) if reach_zone.contains(Point(sub[i]))), None)
                if j is not None and j > 0:
                    sub = sub[j:]
            t = -end_tangent(sub, True, 0.6); p = sub[0]
            hit = next((e for e in np.arange(0.05, 1.55, 0.05) if reach_zone.contains(Point(p + t * e))), None)
            if hit is not None and not reach_zone.contains(Point(p)):
                sub = np.vstack([p + t * hit, sub])
            if end == 1:
                sub = sub[::-1]
        d = dict(id=f"r{st['id']}", prims=[('poly', sub)], width=round(st['width'] * K, 3),
                 caps=[None if c0 else st['caps'][0], None if c1 else st['caps'][1]],
                 free=[(not c0) and st['free'][0], (not c1) and st['free'][1]], var_w=bool(st.get('var_w')) and not (c0 or c1),
                 taper=st.get('taper'), len=round(length(sub), 2), src_len=st['len'])
        if d['var_w']:
            d['axis'] = r2v(st['axis']); d['prof'] = st['prof'] * K
        if d['taper']:
            d['taper'] = (st['taper'][0] * K, st['taper'][1] * K)
        pieces.append(d)
new_ink = MAT.render_clean(None, pieces, [])
print(f'raster: {len(pieces)} pieces transplanted, {sum(p["len"] for p in pieces):.1f} u of line')
comp = unary_union([kept_ink, new_ink])
c = compare(ink2v, comp)
print(f'composite vs aligned raster: IoU {c["iou"]:.4f}, raster-only {c["a_only"]:.1f}, composite-only {c["b_only"]:.1f}')
for name, bb in {'left arm': (-2, 0, 24, 30), 'right arm': (34, -2, 58, 26), 'shoulders': (20, 18, 40, 30), 'body+tail': (10, 26, 56, 66)}.items():
    b_ = box(*bb); a1, a2 = ink2v.intersection(b_), comp.intersection(b_)
    print(f'  {name:10s} IoU {a1.intersection(a2).area / max(a1.union(a2).area, 1e-9):.3f}')
pickle.dump(dict(kept=kept, kept_w=kept_w, kept_recs=kept_recs, dropped=dropped, dropped_polys=dropped_polys, arrow=arrow, kept_fills=kept_fills, kept_lines=kept_lines, pieces=pieces, comp=comp, ink2v=ink2v, boltv=boltv, kept_ink=kept_ink, new_ink=new_ink),
            open(os.path.join(HERE, 'compose2.pkl'), 'wb'))

b = comp.union(ink2v).bounds
vb = (b[0] - 2, b[1] - 2, b[2] - b[0] + 4, b[3] - b[1] + 4); vbs = ' '.join(fmt(v, 2) for v in vb)
W, H = 540, int(540 * vb[3] / vb[2])
def svg(body): return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vbs}" width="{W}" height="{H}">{body}</svg>'
INK = '#1c1b19'
open(os.path.join(HERE, 'step2-clean.svg'), 'w').write(svg(f'<path d="{svgpath(comp, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/>'))
open(os.path.join(HERE, 'step2-raster.svg'), 'w').write(svg(f'<path d="{svgpath(boltv, prec=3, tol=0.003)}" fill="#ad9673"/><path d="{svgpath(ink2v, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/>'))
both = comp.intersection(ink2v)
open(os.path.join(HERE, 'step2-overlay.svg'), 'w').write(svg(
    f'<rect x="{fmt(vb[0])}" y="{fmt(vb[1])}" width="{fmt(vb[2])}" height="{fmt(vb[3])}" fill="#fff"/>'
    f'<path d="{svgpath(ink2v.difference(comp), prec=3, tol=0.003)}" fill="#d6473b" fill-rule="evenodd"/>'
    f'<path d="{svgpath(comp.difference(ink2v), prec=3, tol=0.003)}" fill="#2a9df4" fill-rule="evenodd"/>'
    f'<path d="{svgpath(both, prec=3, tol=0.003)}" fill="#c6c3bc" fill-rule="evenodd"/>'))
parts = [f'<path d="{svgpath(comp, prec=3, tol=0.004)}" fill="#dedbd3" fill-rule="evenodd"/>']
for kind, q in kept_lines:
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    parts.append(f'<path d="{pd}" fill="none" stroke="{"#1c1b19" if kind == "whole" else "#2a9df4"}" stroke-width="0.28"/>')
for p in pieces:
    q = p['prims'][0][1]
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    parts.append(f'<path d="{pd}" fill="none" stroke="#e3761b" stroke-width="0.28"/>')
open(os.path.join(HERE, 'step2-construction.svg'), 'w').write(svg(''.join(parts)))
dparts = [f'<path d="{svgpath(ink2v, prec=3, tol=0.004)}" fill="#f1d9c0" fill-rule="evenodd"/>']
for st in strokes1:
    q = stroke_pts(st); pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    dparts.append(f'<path d="{pd}" fill="none" stroke="{"#d6473b" if st["id"] in dropped else "#1c1b19"}" stroke-width="0.3"/>')
open(os.path.join(HERE, 'diag_status.svg'), 'w').write(svg(''.join(dparts)))
print('viewBox', vbs)
