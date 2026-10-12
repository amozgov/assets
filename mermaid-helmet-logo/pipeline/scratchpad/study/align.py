"""Align the Step-1 vector (logo.svg) to the lightning-mermaid raster by scale + translation; overlay to see what differs."""
import os, sys, re, numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from skimage.measure import find_contours
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union, polygonize
from shapely import affinity
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from geom import sample_path
HERE = os.path.dirname(os.path.abspath(__file__))

# --- raster ink polygon (units: standard line = 1.38) ---
im = np.array(Image.open(os.path.join(HERE, '..', '..', 'images', '3.png')).convert('RGBA')).astype(float)
a = im[:, :, 3] / 255; rgb = im[:, :, :3] * a[:, :, None] + 255 * (1 - a[:, :, None]); mx = rgb.max(axis=2) / 255
PX = 8.816
def contours_to_polygon(field, level):
    polys = []
    for c in find_contours(field, level):
        if len(c) < 8: continue
        pts = np.stack([c[:, 1], c[:, 0]], 1) / PX
        p = Polygon(pts).buffer(0)
        if not p.is_empty and p.area > 0.02: polys.append(p)
    polys.sort(key=lambda p: -p.area); result = None
    for p in polys:
        depth = sum(1 for q in polys if q is not p and q.area > p.area and q.contains(p.representative_point()))
        result = p if result is None else (result.difference(p) if depth % 2 else result.union(p))
    return result.buffer(0)
ink2 = contours_to_polygon(mx, 0.5)
tan = (rgb[:, :, 0] - rgb[:, :, 2]) / 255
bolt = contours_to_polygon(np.where(mx > 0.4, tan, 0.0), 0.12)
bolt = max(getattr(bolt, 'geoms', [bolt]), key=lambda g: g.area)

# --- Step-1 vector polygon (nonzero winding) ---
svg = open(os.path.join(HERE, '..', 'logo.svg')).read()
d = re.search(r' d="([^"]+)"', svg).group(1)
rings = [pts for pts, closed in sample_path(d, 0.03)]
def winding(pt, ring):
    x, y = pt; xs, ys = ring[:, 0], ring[:, 1]; xe, ye = np.roll(xs, -1), np.roll(ys, -1)
    cross = (xe - xs) * (y - ys) - (x - xs) * (ye - ys)
    return int(np.sum((ys <= y) & (ye > y) & (cross > 0)) - np.sum((ys > y) & (ye <= y) & (cross < 0)))
faces = list(polygonize(unary_union([LineString(np.vstack([r, r[:1]])) for r in rings])))
ink1 = unary_union([f for f in faces if sum(winding((f.representative_point().x, f.representative_point().y), r) for r in rings) != 0])
print('vector area %.1f, raster area %.1f' % (ink1.area, ink2.area))

# --- bitmaps ---
R = 4; W, H = 100 * R, 110 * R   # canvas 0..100 x 0..110 units
def rasterize(poly, R, W, H, sx=1.0, tx=0.0, ty=0.0):
    img = Image.new('L', (W, H), 0); dr = ImageDraw.Draw(img)
    parts = sorted(getattr(poly, 'geoms', [poly]), key=lambda p: -p.area)
    for p in parts:
        dr.polygon([((x * sx + tx) * R, (y * sx + ty) * R) for x, y in p.exterior.coords], fill=1)
        for h in p.interiors:
            dr.polygon([((x * sx + tx) * R, (y * sx + ty) * R) for x, y in h.coords], fill=0)
    return np.array(img, dtype=np.float32)
B2 = rasterize(ink2, R, W, H)
best = None
for s in np.arange(0.80, 1.30, 0.02):
    B1 = rasterize(ink1, R, W, H, s, 20, 20)
    corr = fftconvolve(B2, B1[::-1, ::-1], mode='full')
    iy, ix = np.unravel_index(np.argmax(corr), corr.shape)
    dy, dx = iy - (B1.shape[0] - 1), ix - (B1.shape[1] - 1)
    B1s = np.roll(np.roll(B1, dy, 0), dx, 1)
    iou = (B1s * B2).sum() / np.maximum(B1s, B2).sum()
    if best is None or iou > best[0]:
        best = (iou, s, 20 + dx / R, 20 + dy / R)
print('coarse best: IoU %.3f at scale %.3f, shift (%.2f, %.2f)' % best)
iou0, s0, tx0, ty0 = best
for s in np.arange(s0 - 0.02, s0 + 0.021, 0.004):
    for tx in np.arange(tx0 - 0.6, tx0 + 0.61, 0.2):
        for ty in np.arange(ty0 - 0.6, ty0 + 0.61, 0.2):
            B1s = rasterize(ink1, R, W, H, s, tx, ty)
            iou = (B1s * B2).sum() / np.maximum(B1s, B2).sum()
            if iou > best[0]: best = (iou, s, tx, ty)
print('fine best:   IoU %.3f at scale %.3f, shift (%.2f, %.2f)' % best)
iou, s, tx, ty = best
v = affinity.affine_transform(ink1, [s, 0, 0, s, tx, ty])
both = v.intersection(ink2); ronly = ink2.difference(v); vonly = v.difference(ink2)
print('units: vector->raster scale %.3f (vector standard line 1.38 -> %.2f raster units)' % (s, 1.38 * s))
print('area both %.1f, raster-only %.1f, vector-only %.1f' % (both.area, ronly.area, vonly.area))
# per-region: boxes in raster units
boxes = {'head': (31, 23, 49, 40), 'torso': (28, 38, 50, 62), 'tail+scales': (26, 60, 50, 90), 'fin': (46, 52, 66, 72), 'arms (both)': (13, 13, 66, 40)}
from shapely.geometry import box
for k, (x0, y0, x1, y1) in boxes.items():
    bb = box(x0, y0, x1, y1); a_ = ink2.intersection(bb); b_ = v.intersection(bb)
    print('  %-12s IoU %.3f  raster %.1f vector %.1f' % (k, a_.intersection(b_).area / max(a_.union(b_).area, 1e-9), a_.area, b_.area))
# overlay image at 8 px/u
R2 = 8; W2, H2 = 90 * R2, 100 * R2
def mask(poly): return rasterize(poly, R2, W2, H2) > 0
mb, mr, mv = mask(both), mask(ronly), mask(vonly)
out = np.full((H2, W2, 3), 255, np.uint8)
out[mb] = (198, 195, 188); out[mr] = (214, 71, 59); out[mv] = (42, 157, 244)
bm = mask(bolt); out[bm & ~(mb | mr | mv)] = (230, 215, 190)
Image.fromarray(out).crop((10 * R2, 10 * R2, 70 * R2, 92 * R2)).save(os.path.join(HERE, 'align_overlay.png'))
import pickle; pickle.dump(dict(scale=s, tx=tx, ty=ty, iou=iou), open(os.path.join(HERE, 'align.pkl'), 'wb'))
print('wrote align_overlay.png')
