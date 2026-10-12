"""ccrop.py name x y w h [W]: crop of step2-construction.svg (grey ink, kept strokes black, new lines orange)"""
import re, subprocess, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
name = sys.argv[1]; x, y, w, h = map(float, sys.argv[2:6]); W = int(sys.argv[6]) if len(sys.argv) > 6 else 600
H = int(round(W * h / w))
s = open(os.path.join(HERE, 'step2-construction.svg')).read()
s = re.sub(r'viewBox="[^"]*"', f'viewBox="{x} {y} {w} {h}"', s, count=1); s = re.sub(r'width="[^"]*" height="[^"]*"', f'width="{W}" height="{H}"', s, count=1)
s = s.replace('stroke-width="0.28"', 'stroke-width="0.04"').replace('stroke-width="0.22"', 'stroke-width="0.04"').replace('r="2.2"', 'r="0.3"').replace('stroke-width="0.18"', 'stroke-width="0.03"')
tmp = os.path.join(HERE, f'_cc_{name}.svg'); open(tmp, 'w').write(s)
subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=dict(os.environ, NODE_PATH='/opt/node22/lib/node_modules'), check=True, capture_output=True)
os.replace(tmp[:-4] + '.png', os.path.join(HERE, name + '.png')); print(name)
