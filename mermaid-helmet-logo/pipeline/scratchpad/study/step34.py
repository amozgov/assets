"""Step 3: mirror the finished Step 2 figure (tail to the left). Step 4: a trident in the hand that is on the left after
the flip (the pre-flip right hand), built by the rules: a standard straight shaft, a three-prong fork head with filled
tips, the shaft knocked out where it passes behind the hand's loop. Writes step3.svg (mirror only) and step4.svg."""
import os, sys, pickle, math, json, numpy as np, shapely
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import unary_union
from shapely import affinity
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art'))
from lib import svgpath, fmt
HERE = os.path.dirname(os.path.abspath(__file__))
F = pickle.load(open(os.path.join(HERE, 'final2.pkl'), 'rb'))
comp = F['comp']
b = comp.bounds; cx = (b[0] + b[2]) / 2
def mirror(g): return affinity.scale(g, xfact=-1, yfact=1, origin=(cx, 0))
mir = mirror(comp)
INK = '#1c1b19'
def write_svg(name, g, extra=''):
    bb = g.bounds; vb = (bb[0] - 2, bb[1] - 2, bb[2] - bb[0] + 4, bb[3] - bb[1] + 4)
    W, H = 540, int(540 * vb[3] / vb[2])
    open(os.path.join(HERE, name), 'w').write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{" ".join(fmt(v, 2) for v in vb)}" width="{W}" height="{H}"><path d="{svgpath(g, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/>{extra}</svg>')
    print(name, 'bounds', [round(v, 2) for v in bb])
write_svg('step3.svg', mir)
pickle.dump(dict(mir=mir, cx=cx), open(os.path.join(HERE, 'step3.pkl'), 'wb'))

# ---------------- Step 4: the trident, built in the pre-flip frame and mirrored with the figure ----------------
from mat import MAT
C4 = pickle.load(open(os.path.join(HERE, 'compose4.pkl'), 'rb'))
pieces = {p['id']: p for p in C4['pieces']}
HAND_RING = MAT.render_clean(None, [pieces['handR_oval']], [])
HAND_FILL = Polygon(HAND_RING.exterior)                                            # the loop filled to its outer edge: the palm
F_ = np.asarray(HAND_FILL.centroid.coords[0])
STD = 1.38
u = np.array([0.0, -1.0]); n = np.array([1.0, 0.0])                                 # shaft straight up through the fist
def perp(v): return np.array([-v[1], v[0]])
def unit(v): return v / np.linalg.norm(v)
def band(pts, w=STD): return LineString(pts).buffer(w / 2, cap_style='flat', join_style='mitre', mitre_limit=8, quad_segs=24)
L_DOWN, L_UP, R_U, L_P, TIP_H, TIP_B = 34.0, 6.0, 2.6, 3.0, 2.8, 2.3
P_butt = F_ - u * L_DOWN; P_ubot = F_ + u * L_UP; C_u = P_ubot + u * R_U
th = np.linspace(0, math.pi, 64)
U_pts = np.array([C_u - n * R_U * math.cos(t) - u * R_U * math.sin(t) for t in th])  # the U: from the left prong base down round the bottom to the right one
prong_top = C_u + u * L_P
tips = []
def tip(apex_base_centre):
    """a filled isosceles tip: base 2.0 across the prong, 2.4 high, pointing along the shaft"""
    c = apex_base_centre; return Polygon([c - n * TIP_B / 2, c + n * TIP_B / 2, c + u * TIP_H])
parts = []
parts.append(band([P_butt, prong_top + u * 0.8 - u * 0.3]))                           # shaft + middle prong, buried 0.3 in its tip
parts.append(band(U_pts))
for s in (-1, 1):
    base = C_u + n * s * R_U
    parts.append(band([base - u * 0.6, base + u * (L_P - 0.3)]))                        # side prongs, overlapping the U's ends and buried in their tips
    tips.append(tip(base + u * L_P))
tips.append(tip(prong_top + u * 0.8))
trident = unary_union(parts + tips).difference(HAND_FILL)                             # the shaft passes behind the hand: nothing of it inside the loop's outline
comp4 = unary_union([comp, trident]).buffer(0)
mir4 = mirror(comp4)
write_svg('step4.svg', mir4)
pickle.dump(dict(mir=mir, mir4=mir4, trident=mirror(trident), cx=cx, F=F_, notes=[f'trident: standard shaft ({STD}) {L_DOWN:.0f} u below the fist and {L_UP:.0f} u above it to a U of radius {R_U}, three prongs {L_P} u with filled tips {TIP_B} x {TIP_H}; the shaft is knocked out inside the hand loop']), open(os.path.join(HERE, 'step4.pkl'), 'wb'))
print('fist centre', F_.round(2), 'trident bounds', [round(v, 2) for v in trident.bounds])
