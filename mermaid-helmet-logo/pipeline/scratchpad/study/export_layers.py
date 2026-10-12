"""export_layers.py: Step 4 as a layered SVG for hand-off. Groups: figure (Step 3 without the trident hand),
trident-hand (the current open hand, to be redrawn as a grip), trident (shaft + fork, knocked out inside the hand).
Writes step4-layers.svg in the same viewBox as step4.svg. Reads final2.pkl, step4.pkl, compose4.pkl."""
import os, sys, pickle
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely import affinity
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'art'))
from lib import svgpath, fmt
from mat import MAT, compare
F = pickle.load(open(os.path.join(HERE, 'final2.pkl'), 'rb')); S4 = pickle.load(open(os.path.join(HERE, 'step4.pkl'), 'rb'))
C4 = pickle.load(open(os.path.join(HERE, 'compose4.pkl'), 'rb'))
cx = S4['cx']
def mirror(g): return affinity.scale(g, xfact=-1, yfact=1, origin=(cx, 0))
pieces = {p['id']: p for p in C4['pieces']}
hand_fill = mirror(Polygon(MAT.render_clean(None, [pieces['handR_oval']], []).exterior))   # the pre-flip right hand's loop, filled to its outer edge
mir = S4['mir']; trident = S4['trident']
figure = mir.difference(hand_fill); hand = mir.intersection(hand_fill)
check = compare(S4['mir4'], unary_union([figure, hand, trident]))
print(f'layers vs step4: IoU {check["iou"]:.6f}')
b = S4['mir4'].bounds; vb = (b[0] - 2, b[1] - 2, b[2] - b[0] + 4, b[3] - b[1] + 4)
W = 540; H = int(W * vb[3] / vb[2]); INK = '#1c1b19'
def grp(gid, g, note): return f'<g id="{gid}"><desc>{note}</desc><path d="{svgpath(g, prec=3, tol=0.003)}" fill="{INK}" fill-rule="evenodd"/></g>'
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{" ".join(fmt(v, 2) for v in vb)}" width="{W}" height="{H}">'
       + grp('figure', figure, 'Step 3 mirrored figure without the trident hand')
       + grp('trident-hand', hand, 'current open hand around the shaft: to be redrawn as a closed grip')
       + grp('trident', trident, 'current trident: standard shaft, U fork R 2.6, three prongs, triangle tips') + '</svg>')
open(os.path.join(HERE, 'step4-layers.svg'), 'w').write(svg); print('step4-layers.svg', len(svg), 'bytes; hand bounds', [round(v, 2) for v in hand.bounds])
