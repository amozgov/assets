"""s12.py name x y w h [W]: Step 1 clean | Step 2 clean, same window (vector units)"""
import re, subprocess, os, sys
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
name = sys.argv[1]; x, y, w, h = map(float, sys.argv[2:6]); W = int(sys.argv[6]) if len(sys.argv) > 6 else 600
H = int(round(W * h / w)); panels = []
for f in ('redraw_v-clean.svg', 'step2-clean.svg'):
    s = open(os.path.join(HERE, f)).read()
    s = re.sub(r'viewBox="[^"]*"', f'viewBox="{x} {y} {w} {h}"', s, count=1); s = re.sub(r'width="[^"]*" height="[^"]*"', f'width="{W}" height="{H}"', s, count=1)
    tmp = os.path.join(HERE, f'_s12_{name}.svg'); open(tmp, 'w').write(s)
    subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=dict(os.environ, NODE_PATH='/opt/node22/lib/node_modules'), check=True, capture_output=True)
    panels.append(Image.open(tmp[:-4] + '.png').convert('RGB'))
out = Image.new('RGB', (2 * W + 20, H), (190, 190, 190))
for i, p in enumerate(panels): out.paste(p, (i * (W + 20), 0))
out.save(os.path.join(HERE, name + '.png')); print(name, out.size)
