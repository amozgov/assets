"""dbg.py name x y w h [W]: Step 1 ink (light grey), Step 2 composite outline (black), kept axes (blue) / dropped axes (red)
with ids, construction points from final2.pkl (teal). Reads final2.pkl for comp, keep and points."""
import pickle, sys, os, numpy as np, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'art')); from lib import svgpath, fmt
from shapely.geometry import box, Point
HERE = os.path.dirname(os.path.abspath(__file__))
name = sys.argv[1]; x0, y0, w, h = map(float, sys.argv[2:6]); W = int(sys.argv[6]) if len(sys.argv) > 6 else 700
H = int(round(W * h / w)); win = box(x0, y0, x0 + w, y0 + h)
S1 = pickle.load(open(os.path.join(HERE, 'step1.pkl'), 'rb')); F = pickle.load(open(os.path.join(HERE, 'final2.pkl'), 'rb'))
kept_ids = F['kept_ids']; comp = F['comp']; ink = S1['ink']
parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0} {y0} {w} {h}" width="{W}" height="{H}"><rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#fff"/>',
         f'<path d="{svgpath(ink.intersection(win.buffer(1)), prec=3, tol=0.003)}" fill="#e4e1da" fill-rule="evenodd"/>',
         f'<path d="{svgpath(comp.intersection(win.buffer(1)), prec=3, tol=0.002)}" fill="none" stroke="#111" stroke-width="{w / W * 2.2}" fill-rule="evenodd"/>']
for st in S1['strokes']:
    ax = np.asarray(st['axis'], float)
    if not win.buffer(0.5).intersects(Point(ax[0]).buffer(0.01).union(Point(ax[-1]).buffer(0.01))) and not any(win.contains(Point(p)) for p in ax[::5]): continue
    col = '#2a6fd4' if st['id'] in kept_ids else '#d6473b'
    parts.append(f'<path d="M{" L".join(f"{fmt(a, 3)} {fmt(b, 3)}" for a, b in ax)}" fill="none" stroke="{col}" stroke-width="{w / W * 1.6}"/>')
    k = len(ax) // 2; parts.append(f'<text x="{ax[k][0]:.2f}" y="{ax[k][1]:.2f}" font-size="{w / W * 14}" fill="{col}">{st["id"]}</text>')
    parts.append(f'<circle cx="{ax[0][0]:.3f}" cy="{ax[0][1]:.3f}" r="{w / W * 3}" fill="{col}"/>')
for key in ('c_k', 'c5', 'C4', 'c3', 'P2', 'cj', 'A_', 'B_'):
    if key in F:
        p = np.asarray(F[key]); parts.append(f'<circle cx="{p[0]:.3f}" cy="{p[1]:.3f}" r="{w / W * 4}" fill="none" stroke="#1f8f80" stroke-width="{w / W * 1.5}"/><text x="{p[0] + w / W * 6:.2f}" y="{p[1]:.2f}" font-size="{w / W * 14}" fill="#1f8f80">{key}</text>')
parts.append('</svg>')
tmp = os.path.join(HERE, f'_dbg_{name}.svg'); open(tmp, 'w').write(''.join(parts))
subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=dict(os.environ, NODE_PATH='/opt/node22/lib/node_modules'), check=True, capture_output=True)
os.replace(tmp[:-4] + '.png', os.path.join(HERE, name + '.png')); print(name, W, H)
