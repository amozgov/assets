"""spikes_check.py: list sharp boundary features of the composite (interior angle < AMAX deg) and tiny holes/parts"""
import pickle, sys, numpy as np
from shapely.geometry import Polygon
F = pickle.load(open('final2.pkl', 'rb')); comp = F['comp']
AMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
def ring_feats(coords, hole):
    pts = np.array(coords[:-1]); n = len(pts); out = []
    for i in range(n):
        p = pts[i]; a = pts[i - 1]; b = pts[(i + 1) % n]
        # walk back/forward ~0.12 u for a stable direction
        j = i; d = 0.0
        while d < 0.12:
            j2 = (j - 1) % n; d += np.linalg.norm(pts[j2] - pts[j]); j = j2
            if j == i: break
        a = pts[j]; k = i; d = 0.0
        while d < 0.12:
            k2 = (k + 1) % n; d += np.linalg.norm(pts[k2] - pts[k]); k = k2
            if k == i: break
        b = pts[k]
        u, v = a - p, b - p
        cosang = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-12)
        ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
        if ang < AMAX:
            cross = u[0] * v[1] - u[1] * v[0]
            # exterior rings are CCW in shapely? orientation varies; report sign only
            out.append((round(float(p[0]), 2), round(float(p[1]), 2), round(float(ang), 1), 'hole' if hole else 'outer', 'L' if cross > 0 else 'R'))
    return out
feats = []
for g in getattr(comp, 'geoms', [comp]):
    feats += ring_feats(list(g.exterior.coords), False)
    for h in g.interiors: feats += ring_feats(list(h.coords), True)
feats.sort()
print(len(feats), 'sharp features (<', AMAX, 'deg)')
for f in feats: print(f)
print('parts:', len(getattr(comp, 'geoms', [comp])), 'smallest parts:', sorted(round(g.area, 3) for g in getattr(comp, 'geoms', [comp]))[:5])
holes = [Polygon(h).area for g in getattr(comp, 'geoms', [comp]) for h in g.interiors]
print('holes:', len(holes), 'smallest:', sorted(round(a, 4) for a in holes)[:8])
