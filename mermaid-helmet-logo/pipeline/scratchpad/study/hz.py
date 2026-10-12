"""hz.py name x y w h [W]: aligned raster (nearest-neighbour) | our clean render, for a window in vector units"""
import pickle, numpy as np, re, subprocess, os, sys
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
A = pickle.load(open(os.path.join(HERE, 'align2.pkl'), 'rb'))['mat']; PX = 8.816
im = Image.open(os.path.join(HERE, '..', '..', 'images', '3.png')).convert('RGBA'); bg = Image.new('RGBA', im.size, (255, 255, 255, 255)); bg.alpha_composite(im); im = bg.convert('RGB')
a, b_, c, d, xo, yo = A; Minv = np.linalg.inv(np.array([[a, b_], [c, d]]))
name = sys.argv[1]; x, y, w, h = map(float, sys.argv[2:6]); W = int(sys.argv[6]) if len(sys.argv) > 6 else 600
H = int(round(W * h / w))
corners = np.array([[x, y], [x + w, y], [x, y + h], [x + w, y + h]]) - [xo, yo]; rc = corners @ Minv.T
box = (rc[:, 0].min() * PX, rc[:, 1].min() * PX, rc[:, 0].max() * PX, rc[:, 1].max() * PX)
ras = im.transform((W, H), Image.EXTENT, box, Image.NEAREST)
s = open(os.path.join(HERE, 'step2-clean.svg')).read()
s = re.sub(r'viewBox="[^"]*"', f'viewBox="{x} {y} {w} {h}"', s, count=1); s = re.sub(r'width="[^"]*" height="[^"]*"', f'width="{W}" height="{H}"', s, count=1)
tmp = os.path.join(HERE, '_hz.svg'); open(tmp, 'w').write(s)
subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=dict(os.environ, NODE_PATH='/opt/node22/lib/node_modules'), check=True, capture_output=True)
ours = Image.open(tmp[:-4] + '.png').convert('RGB')
out = Image.new('RGB', (2 * W + 20, H), (190, 190, 190)); out.paste(ras, (0, 0)); out.paste(ours, (W + 20, 0)); out.save(os.path.join(HERE, name + '.png')); print(name, out.size)
