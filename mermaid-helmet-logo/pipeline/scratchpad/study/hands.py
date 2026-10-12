"""Rebuild the two hands as clean curves: an egg-shaped loop (palm + fingers) fitted to the raster's outer silhouette,
and the inner thumb hook from the traced pieces, both smoothed. Replaces the traced hand pieces in the composite."""
import os, sys, pickle, numpy as np, shapely
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import unary_union
from scipy.interpolate import splprep, splev
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from lib import svgpath, fmt
from geom import *
from mat import MAT, compare
HERE = os.path.dirname(os.path.abspath(__file__))
C = pickle.load(open(os.path.join(HERE, 'compose2.pkl'), 'rb'))
pieces, kept_ink, kept_fills, ink2v, boltv = C['pieces'], C['kept_ink'], C['kept_fills'], C['ink2v'], C['boltv']
HW = 0.45          # hand pen 0.9
HANDS = {'L': dict(box=(3, -4, 15, 9.5), arms=('r11', 'r13')), 'R': dict(box=(44, -8, 56, 3.5), arms=('r93', 'r103'))}


def ppts(p):
    return p['prims'][0][1]


def smooth_closed(pts, s=0.6, n=240):
    pts = resample(np.vstack([pts, pts[:1]]), 0.1)[:-1]
    tck, _ = splprep([pts[:, 0], pts[:, 1]], s=s * len(pts) * 0.01, per=True, k=3)
    u = np.linspace(0, 1, n, endpoint=False)
    x, y = splev(u, tck)
    return np.stack([x, y], 1)


def smooth_open(pts, s=0.6, n=120):
    pts = resample(pts, 0.1)
    if len(pts) < 6:
        return pts
    tck, _ = splprep([pts[:, 0], pts[:, 1]], s=s * len(pts) * 0.01, k=3)
    x, y = splev(np.linspace(0, 1, n), tck)
    return np.stack([x, y], 1)


def chain(polys, tol=0.7):
    """join polylines end-to-end by proximity into chains (longest path first)"""
    polys = [np.asarray(p) for p in polys]
    chains = []
    while polys:
        cur = polys.pop(0)
        changed = True
        while changed:
            changed = False
            for i, q in enumerate(polys):
                for flip in (False, True):
                    qq = q[::-1] if flip else q
                    if np.linalg.norm(cur[-1] - qq[0]) < tol:
                        cur = np.vstack([cur, qq[1:]]); polys.pop(i); changed = True; break
                    if np.linalg.norm(cur[0] - qq[-1]) < tol:
                        cur = np.vstack([qq[:-1], cur]); polys.pop(i); changed = True; break
                if changed:
                    break
        chains.append(cur)
    return sorted(chains, key=lambda c: -length(c))


def skeleton_path(poly, res=25):
    """longest skeleton path of a polygon (units), via a raster skeleton at `res` px per unit"""
    from PIL import Image, ImageDraw
    from skimage.morphology import skeletonize
    from skan import Skeleton
    b = poly.bounds
    x0, y0 = b[0] - 0.5, b[1] - 0.5
    W, Hh = int((b[2] - b[0] + 1) * res), int((b[3] - b[1] + 1) * res)
    img = Image.new('L', (W, Hh), 0); dr = ImageDraw.Draw(img)
    for g in getattr(poly, 'geoms', [poly]):
        dr.polygon([((x - x0) * res, (y - y0) * res) for x, y in g.exterior.coords], fill=1)
        for h in g.interiors:
            dr.polygon([((x - x0) * res, (y - y0) * res) for x, y in h.coords], fill=0)
    sk = skeletonize(np.array(img) > 0)
    S = Skeleton(sk)
    # longest geodesic path: concatenate branches greedily from the longest one
    paths = [S.path_coordinates(i) for i in range(S.n_paths)]
    paths = sorted(paths, key=lambda q: -len(q))
    pts = [np.stack([q[:, 1], q[:, 0]], 1) / res + [x0, y0] for q in paths]
    return chain(pts, tol=0.3)[0] if pts else None


new_pieces = [p for p in pieces]
hand_strokes = {}
for key, H in HANDS.items():
    bb = box(*H['box'])
    arm = sorted([p for p in pieces if p['id'] in H['arms']], key=lambda p: -p['len'])[:2]
    hand_ink = ink2v.intersection(bb)
    hand_ink = max(getattr(hand_ink, 'geoms', [hand_ink]), key=lambda g: g.area)
    # wrist line: through the two forearm line ends; the loop lives on the hand side of it
    wrist_pts = []
    for p in arm:
        q = ppts(p)
        wrist_pts.append(q[0] if bb.contains(Point(q[0])) else q[-1])
    w0, w1 = wrist_pts
    wmid = (w0 + w1) / 2
    adir = np.mean([end_tangent(ppts(p), bb.contains(Point(ppts(p)[0])), 3.0) for p in arm], axis=0)
    adir /= np.linalg.norm(adir)                                 # points from the hand down the forearm
    outer = Polygon(hand_ink.exterior)
    shr = outer.buffer(-HW); shr = max(getattr(shr, 'geoms', [shr]), key=lambda g: g.area)
    ring = resample(np.asarray(shr.exterior.coords), 0.1)
    keep = ((ring - wmid) @ adir) < -0.2                         # drop the forearm side
    # longest kept arc (the ring is a loop: rotate so the arc is contiguous)
    idx = np.where(keep)[0]
    if len(idx) and not keep.all():
        start = next(i for i in range(len(ring)) if keep[i] and not keep[i - 1])
        order = [(start + k) % len(ring) for k in range(len(ring))]
        arc = ring[[i for i in order if keep[i]]]
    else:
        arc = ring
    arc = arc[np.r_[True, np.linalg.norm(np.diff(arc, axis=0), axis=1) > 1e-6]]
    arc_c = np.vstack([arc, arc[:1]])
    tck, _ = splprep([arc_c[:, 0], arc_c[:, 1]], s=1.0 * len(arc_c) * 0.01, per=True, k=3)
    x, y = splev(np.linspace(0, 1, 240, endpoint=False), tck)
    oval = np.stack([x, y], 1)
    oval_poly = Polygon(oval)
    # inner ink: ink strictly inside the loop, away from the wrist mass
    inner_ink = hand_ink.intersection(oval_poly.buffer(-HW - 0.12))
    cutter = Polygon([wmid - adir * 1.3 + np.array([-adir[1], adir[0]]) * 6, wmid - adir * 1.3 - np.array([-adir[1], adir[0]]) * 6,
                      wmid + adir * 9 - np.array([-adir[1], adir[0]]) * 6, wmid + adir * 9 + np.array([-adir[1], adir[0]]) * 6])
    print(f'  [{key}] hand ink {hand_ink.area:.1f}, oval poly {oval_poly.area:.1f} ({len(arc)} arc pts of {len(ring)}), inner before cutter {inner_ink.area:.2f}, wrist {np.round(w0,1)} {np.round(w1,1)} adir {np.round(adir,2)}')
    dbg = [f'<path d="{svgpath(hand_ink, prec=3, tol=0.003)}" fill="#f1d9c0" fill-rule="evenodd"/>',
           f'<path d="{svgpath(cutter, prec=3, tol=0.003)}" fill="#c6e4f7" fill-opacity="0.6"/>',
           'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in oval) + 'Z']
    dbg[2] = f'<path d="{dbg[2]}" fill="none" stroke="#d6473b" stroke-width="0.15"/>'
    dbg.append('<path d="M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in arc) + f'" fill="none" stroke="#2a9df4" stroke-width="0.1"/>')
    b_ = bb.bounds
    open(os.path.join(HERE, f'dbg_hand{key}.svg'), 'w').write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{b_[0]} {b_[1]} {b_[2]-b_[0]} {b_[3]-b_[1]}" width="600" height="{int(600*(b_[3]-b_[1])/(b_[2]-b_[0]))}">{"".join(dbg)}</svg>')
    inner_ink = inner_ink.difference(cutter).buffer(-0.06).buffer(0.06)
    print(f'  [{key}] inner after cutter {inner_ink.area:.2f}')
    inner_ink = max(getattr(inner_ink, 'geoms', [inner_ink]), key=lambda g: g.area)
    # inner strokes: medial axis of the inner ink, one stroke per branch (no chaining across junctions)
    mi = MAT(inner_ink, step=0.03, log=lambda *a: None, min_residual=0.15)
    mi.prune(spur_len=0.8, end_radius=0.25); mi.merge(straight_deg=35)
    hooks = []
    for bid, raw in mi.bmap.items():
        if length(raw) < 0.6:
            continue
        h = smooth_open(raw, s=0.8)
        k0, k1 = mi.ends[bid]
        free = [mi.deg(k0) == 1, mi.deg(k1) == 1]
        for end in (0, 1):
            p = h[0] if end == 0 else h[-1]
            t = -end_tangent(h, end == 0, 0.5)
            if oval_poly.exterior.distance(Point(p)) < 1.1:        # ends at the loop: run on to its centreline
                e = 0.0
                while e < 2.0 and oval_poly.contains(Point(p + t * e)):
                    e += 0.05
                h = np.vstack([p + t * e + t * 0.1, h]) if end == 0 else np.vstack([h, p + t * e + t * 0.1])
                free[end] = False
            elif not free[end]:                                     # inner junction: reach 0.15 past the node
                h = np.vstack([p + t * 0.15, h]) if end == 0 else np.vstack([h, p + t * 0.15])
        hooks.append((h, free))
    print(f'  [{key}] inner strokes: {[round(length(h), 2) for h, _ in hooks]}')
    hook = max(hooks, key=lambda x: length(x[0]))[0]
    strokes = [dict(id=f'hand{key}_oval', prims=[('poly', np.vstack([oval, oval[:1]]))], width=0.9, caps=[None, None], free=[False, False], var_w=False, taper=None, len=round(length(oval), 2), src_len=0, closed=True)]
    for i, (h, free) in enumerate(hooks):
        strokes.append(dict(id=f'hand{key}_in{i}', prims=[('poly', h)], width=0.9, caps=[None, None], free=free, var_w=False, taper=None, len=round(length(h), 2), src_len=0))
    hand_strokes[key] = strokes
    # arm lines: their wrist ends land on the oval centreline
    for p in arm:
        q = ppts(p)
        for end in (0, 1):
            pt = q[0] if end == 0 else q[-1]
            if not bb.contains(Point(pt)):
                continue
            # trim while inside the oval
            k = 0
            while k < len(q) - 3 and oval_poly.contains(Point(q[k] if end == 0 else q[len(q) - 1 - k])):
                k += 1
            if k:
                q = q[k:] if end == 0 else q[:len(q) - k]
            pt = q[0] if end == 0 else q[-1]
            t = -end_tangent(q, end == 0, 0.6)
            e = 0.0
            while e < 2.5 and not oval_poly.contains(Point(pt + t * e)):
                e += 0.05
            q = np.vstack([pt + t * e, q]) if end == 0 else np.vstack([q, pt + t * e])
        p['prims'] = [('poly', q)]; p['len'] = round(length(q), 2)
    new_pieces = [p for p in new_pieces if p['id'] in H['arms'] or not bb.contains(Point(ppts(p).mean(axis=0)))] + strokes
    g_hand = MAT.render_clean(None, strokes + arm, [])
    c = compare(ink2v.intersection(bb), g_hand.intersection(bb))
    print(f'hand {key}: oval {length(oval):.1f} u, hook {length(hook):.1f} u, inner ink {inner_ink.area:.1f}; IoU in box {c["iou"]:.3f}')

new_ink = MAT.render_clean(None, new_pieces, [])
comp = unary_union([kept_ink, new_ink])
c = compare(ink2v, comp)
print(f'composite v3 vs aligned raster: IoU {c["iou"]:.4f}, raster-only {c["a_only"]:.1f}, composite-only {c["b_only"]:.1f}')
pickle.dump(dict(pieces=new_pieces, kept=C['kept'], kept_w=C.get('kept_w'), kept_fills=kept_fills, kept_lines=C['kept_lines'], comp=comp, ink2v=ink2v, boltv=boltv, kept_ink=kept_ink, new_ink=new_ink),
            open(os.path.join(HERE, 'compose3.pkl'), 'wb'))
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
    q = ppts(p)
    pd = 'M' + ' L'.join(f'{fmt(x, 2)} {fmt(y, 2)}' for x, y in q[::3])
    parts.append(f'<path d="{pd}" fill="none" stroke="#e3761b" stroke-width="0.28"/>')
open(os.path.join(HERE, 'step2-construction.svg'), 'w').write(svg(''.join(parts)))
print('written')
