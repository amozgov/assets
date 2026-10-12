"""notches.py [AMIN DMAX CMIN]: turn-based outline inspector of final2.pkl's composite, with Step 1's ink as baseline.
Lists (a) jogs/steps: two opposite turns of >= AMIN deg within DMAX u of path, (b) sharp corners >= CMIN deg,
leaving out those that Step 1's own outline has at the same place (within 0.08 u)."""
import pickle, sys, math, numpy as np
from shapely.geometry import Point
from shapely.strtree import STRtree
F = pickle.load(open('final2.pkl', 'rb')); comp = F['comp']
S1 = pickle.load(open('step1.pkl', 'rb')); mine = S1['ink']
AMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
DMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 0.12
CMIN = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
def ring_turns(coords):
    P = np.array(coords[:-1]); n = len(P)
    seg = np.roll(P, -1, axis=0) - P; L = np.linalg.norm(seg, axis=1); D = seg / (L[:, None] + 1e-12)
    Dp = np.roll(D, 1, axis=0)
    cr = Dp[:, 0] * D[:, 1] - Dp[:, 1] * D[:, 0]; dt = (Dp * D).sum(1)
    return P, L, np.degrees(np.arctan2(cr, dt))
def inspect(geom):
    jogs, corners = [], []
    for g in getattr(geom, 'geoms', [geom]):
        gs = g.simplify(0.002)
        for r in [gs.exterior] + list(gs.interiors):
            P, L, ang = ring_turns(list(r.coords)); n = len(P)
            for i in range(n):
                if abs(ang[i]) >= CMIN: corners.append((round(abs(ang[i]), 1), round(float(P[i][0]), 2), round(float(P[i][1]), 2)))
                if abs(ang[i]) < AMIN: continue
                d = 0.0; j = i
                while True:
                    d += L[j % n]; j += 1
                    if d > DMAX: break
                    if abs(ang[j % n]) >= AMIN and np.sign(ang[j % n]) != np.sign(ang[i]):
                        jogs.append((round(min(abs(ang[i]), abs(ang[j % n])), 1), round(d, 3), round(float(P[i][0]), 2), round(float(P[i][1]), 2))); break
    return jogs, corners
j2, c2 = inspect(comp); j1, c1 = inspect(mine)
def novel(items, base, k):
    pts = [Point(b[k], b[k + 1]) for b in base]; tree = STRtree(pts) if pts else None
    out = []
    for it in items:
        p = Point(it[k], it[k + 1])
        if tree is None or not any(pts[i].distance(p) < 0.08 for i in tree.query(p.buffer(0.08))): out.append(it)
    return out
jn, cn = novel(j2, j1, 2), novel(c2, c1, 1)
jn.sort(reverse=True); cn.sort(reverse=True)
print(len(j2), 'opposite-turn pairs in Step 2,', len(jn), f'not in Step 1  (turn, gap, x, y); AMIN {AMIN} DMAX {DMAX}')
for j in jn[:80]: print(tuple(float(v) for v in j))
print(len(c2), 'corners >=', CMIN, 'deg in Step 2,', len(cn), 'not in Step 1')
for c in cn[:60]: print(tuple(float(v) for v in c))
