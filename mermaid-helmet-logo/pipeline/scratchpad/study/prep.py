"""Build stroke sets for both drawings with the MAT pipeline and pickle them (vector = Step 1, raster = lightning)."""
import os, sys, re, pickle, time, numpy as np
from PIL import Image
from skimage.measure import find_contours
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union, polygonize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from geom import *
from mat import MAT, compare
HERE = os.path.dirname(os.path.abspath(__file__))

# vector (Step 1)
svg = open(os.path.join(HERE, '..', 'logo.svg')).read()
d = re.search(r' d="([^"]+)"', svg).group(1)
rings = [pts for pts, closed in sample_path(d, 0.03)]
def winding(pt, ring):
    x, y = pt; xs, ys = ring[:, 0], ring[:, 1]; xe, ye = np.roll(xs, -1), np.roll(ys, -1)
    cross = (xe - xs) * (y - ys) - (x - xs) * (ye - ys)
    return int(np.sum((ys <= y) & (ye > y) & (cross > 0)) - np.sum((ys > y) & (ye <= y) & (cross < 0)))
faces = list(polygonize(unary_union([LineString(np.vstack([r, r[:1]])) for r in rings])))
ink1 = unary_union([f for f in faces if sum(winding((f.representative_point().x, f.representative_point().y), r) for r in rings) != 0])
t0 = time.time()
m1 = MAT(ink1, step=0.03)
m1.prune(1.6, 0.25); m1.merge(30)
ex1 = unary_union([Polygon(p).buffer(0.05) for p in m1.islands])
strokes1 = m1.strokes(tol=0.02, exclude=ex1)
fills1 = [Polygon(p).buffer(0) for p in m1.islands]
clean1 = m1.render_clean(strokes1, fills1)
print('vector: %d strokes, clean IoU %.4f (%.0fs)' % (len(strokes1), compare(ink1, clean1)['iou'], time.time() - t0))

# raster (lightning)
im = np.array(Image.open(os.path.join(HERE, '..', '..', 'images', '3.png')).convert('RGBA')).astype(float)
a = im[:, :, 3] / 255; rgb = im[:, :, :3] * a[:, :, None] + 255 * (1 - a[:, :, None]); mx = rgb.max(axis=2) / 255
PX = 8.816
def contours_to_polygon(field, level, smooth=0):
    polys = []
    for c in find_contours(field, level):
        if len(c) < 8: continue
        pts = np.stack([c[:, 1], c[:, 0]], 1) / PX
        if smooth:
            k = int(smooth)
            pts = np.stack([np.convolve(np.pad(pts[:, i], (k, k), mode='wrap'), np.ones(2 * k + 1) / (2 * k + 1), mode='valid') for i in (0, 1)], 1)
        p = Polygon(pts).buffer(0)
        if not p.is_empty and p.area > 0.02: polys.append(p)
    polys.sort(key=lambda p: -p.area); result = None
    for p in polys:
        depth = sum(1 for q in polys if q is not p and q.area > p.area and q.contains(p.representative_point()))
        result = p if result is None else (result.difference(p) if depth % 2 else result.union(p))
    return result.buffer(0)
ink2 = contours_to_polygon(mx, 0.5, smooth=1)
tan = (rgb[:, :, 0] - rgb[:, :, 2]) / 255
bolt = contours_to_polygon(np.where(mx > 0.4, tan, 0.0), 0.12)
bolt = max(getattr(bolt, 'geoms', [bolt]), key=lambda g: g.area)
t0 = time.time()
m2 = MAT(ink2, step=0.03, min_residual=0.3)
m2.prune(1.6, 0.25); m2.merge(30)
ex2 = unary_union([Polygon(p).buffer(0.05) for p in m2.islands]) if m2.islands else None
strokes2 = m2.strokes(tol=0.05, exclude=ex2, width_snap=(0.9, 1.15, 1.38, 1.60), cap_snap_deg=12, smooth=3, corner_deg=50, snap_tol=0.1, min_junction_len=1.0, taper_ratio=0.5, rms=0.03)
fills2 = [Polygon(p).buffer(0) for p in m2.islands]
clean2 = m2.render_clean(strokes2, fills2)
print('raster: %d strokes, clean IoU %.4f (%.0fs)' % (len(strokes2), compare(ink2, clean2)['iou'], time.time() - t0))
align = pickle.load(open(os.path.join(HERE, 'align.pkl'), 'rb'))
pickle.dump(dict(ink1=ink1, strokes1=strokes1, fills1=fills1, ink2=ink2, strokes2=strokes2, fills2=fills2, bolt=bolt, align=align),
            open(os.path.join(HERE, 'prep.pkl'), 'wb'))
print('pickled')
