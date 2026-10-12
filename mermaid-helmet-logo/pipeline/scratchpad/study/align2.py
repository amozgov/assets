"""Refine the raster->vector alignment with a similarity transform (scale, rotation, shift) by ICP on the shared body."""
import os, pickle, numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import box
from shapely import affinity
from geom import resample
HERE = os.path.dirname(os.path.abspath(__file__))
P = pickle.load(open(os.path.join(HERE, 'prep.pkl'), 'rb'))
ink1, ink2 = P['ink1'], P['ink2']
S, TX, TY = P['align']['scale'], P['align']['tx'], P['align']['ty']
# dense boundary samples of the vector
def boundary_pts(poly, step=0.05):
    out = []
    for p in getattr(poly, 'geoms', [poly]):
        for r in [p.exterior, *p.interiors]:
            out.append(resample(np.asarray(r.coords), step))
    return np.vstack(out)
V = boundary_pts(ink1); tree = cKDTree(V)
# raster boundary samples inside the shared body (raster units): head+torso and tail+fin, arms excluded
M = box(30, 23, 50, 62).union(box(24, 52, 66, 90))
R = boundary_pts(ink2)
import shapely
R = R[shapely.contains_xy(M, R[:, 0], R[:, 1])]
# transform raster -> vector: p_v = (Rot(th) * p_r) * k + t ; start from the coarse fit (inverse of vector->raster)
k, th, t = 1 / S, 0.0, np.array([-TX / S, -TY / S])
def apply(pts, k, th, t):
    c, s_ = np.cos(th), np.sin(th)
    return (pts @ np.array([[c, -s_], [s_, c]]).T) * k + t
for it in range(25):
    Q = apply(R, k, th, t)
    d, idx = tree.query(Q)
    keep = d < (0.6 if it < 5 else 0.35)
    A, B = R[keep], V[idx[keep]]
    # Umeyama similarity fit A -> B
    ma, mb = A.mean(0), B.mean(0)
    H = (A - ma).T @ (B - mb) / len(A)
    U, D, Vt = np.linalg.svd(H)
    Rm = Vt.T @ U.T
    if np.linalg.det(Rm) < 0:
        Vt[-1] *= -1; Rm = Vt.T @ U.T
    var_a = ((A - ma) ** 2).sum() / len(A)
    k = float(np.trace(np.diag(D)) / var_a)
    th = float(np.arctan2(Rm[1, 0], Rm[0, 0]))
    t = mb - apply(ma[None], k, th, np.zeros(2))[0]
    if it in (0, 4, 9, 24):
        print(f'iter {it}: scale {k:.4f} rot {np.degrees(th):.3f} deg shift ({t[0]:.3f}, {t[1]:.3f}); matched {keep.mean():.0%}, median dist {np.median(d[keep]):.3f}, 90% {np.percentile(d[keep], 90):.3f}')
Q = apply(R, k, th, t); d, _ = tree.query(Q)
print('final residual over the shared body: median %.3f u, 90%% %.3f u, >0.3u: %.1f%%' % (np.median(d), np.percentile(d, 90), 100 * (d > 0.3).mean()))
# per-part residuals
for name, bb in {'head': (31, 23, 49, 40), 'torso': (28, 38, 50, 62), 'tail+scales': (26, 60, 50, 90), 'fin': (46, 52, 66, 72)}.items():
    m_ = shapely.contains_xy(box(*bb), R[:, 0], R[:, 1])
    print(f'  {name:12s} median {np.median(d[m_]):.3f}  90% {np.percentile(d[m_], 90):.3f}')
c, s_ = np.cos(th), np.sin(th)
# shapely affine: x' = a x + b y + xoff ; y' = d x + e y + yoff
mat = [k * c, -k * s_, k * s_, k * c, t[0], t[1]]
ink2v = affinity.affine_transform(ink2, mat)
print('IoU vector vs transformed raster (whole):', round(ink1.intersection(ink2v).area / ink1.union(ink2v).area, 4))
pickle.dump(dict(k=k, th=th, t=t, mat=mat), open(os.path.join(HERE, 'align2.pkl'), 'wb'))
