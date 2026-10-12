"""assets for the Study 2 board (vector-only version): Step 1 and Step 2 renders, construction, Step 1 | Step 2 zoom grid"""
import os, re, pickle, subprocess, json
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'board4'); os.makedirs(OUT, exist_ok=True)
ENV = dict(os.environ, NODE_PATH='/opt/node22/lib/node_modules')
def render(svg_path, out_png, width, viewbox=None):
    s = open(svg_path).read()
    if viewbox: s = re.sub(r'viewBox="[^"]*"', f'viewBox="{viewbox}"', s, count=1)
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', s).group(1).split()]
    h = int(round(width * vb[3] / vb[2]))
    s = re.sub(r'width="[^"]*" height="[^"]*"', f'width="{width}" height="{h}"', s, count=1)
    tmp = os.path.join(OUT, '_tmp.svg'); open(tmp, 'w').write(s)
    subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=ENV, check=True, capture_output=True)
    os.replace(tmp[:-4] + '.png', out_png); return out_png
s2 = open(os.path.join(HERE, 'step2-clean.svg')).read(); vb2 = re.search(r'viewBox="([^"]+)"', s2).group(1)
render(os.path.join(HERE, 'step2-clean.svg'), os.path.join(OUT, 'step2-clean.png'), 1080)
render(os.path.join(HERE, 'step2-construction.svg'), os.path.join(OUT, 'step2-construction.png'), 1080)
render(os.path.join(HERE, 'redraw_v-clean.svg'), os.path.join(OUT, 'step1-clean.png'), 1080, viewbox=vb2)
im = Image.open(os.path.join(HERE, '..', '..', 'images', '3.png')).convert('RGBA'); bg = Image.new('RGBA', im.size, (255, 255, 255, 255)); bg.alpha_composite(im); bg.convert('RGB').save(os.path.join(OUT, 'source2.png'))
ZOOMS = [('bun, spear removed', (25.2, 9.2), 4.4), ('left shoulder', (30.3, 16.6), 3.6), ('left armpit', (25.9, 21.0), 3.6), ('breast fold', (23.4, 25.2), 3.6),
         ('right shoulder', (34.2, 16.9), 3.6), ('right armpit, thorns', (39.2, 22.4), 5.6), ('left hand', (10.9, 3.8), 7.6), ('right hand', (52.5, -2.6), 7.6)]
W = 360
try: font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf', 20)
except Exception: font = ImageFont.load_default()
tiles = []
for label, c, S in ZOOMS:
    x, y = c[0] - S / 2, c[1] - S / 2; panels = []
    for f in ('redraw_v-clean.svg', 'step2-clean.svg'):
        s = open(os.path.join(HERE, f)).read()
        s = re.sub(r'viewBox="[^"]*"', f'viewBox="{x} {y} {S} {S}"', s, count=1); s = re.sub(r'width="[^"]*" height="[^"]*"', f'width="{W}" height="{W}"', s, count=1)
        tmp = os.path.join(OUT, '_z.svg'); open(tmp, 'w').write(s)
        subprocess.run(['node', os.path.join(HERE, '..', 'art', 'render.js'), tmp], env=ENV, check=True, capture_output=True)
        panels.append(Image.open(tmp[:-4] + '.png').convert('RGB'))
    tile = Image.new('RGB', (2 * W + 8, W + 34), (232, 229, 222)); tile.paste(panels[0], (0, 34)); tile.paste(panels[1], (W + 8, 34))
    ImageDraw.Draw(tile).text((6, 6), f'{label} ({S:g} u): Step 1 | Step 2', fill=(28, 27, 25), font=font)
    tiles.append(tile)
COLS, GAP = 4, 24; rows = [tiles[i:i + COLS] for i in range(0, len(tiles), COLS)]; tw, th = tiles[0].width, tiles[0].height
strip = Image.new('RGB', (COLS * tw + (COLS - 1) * GAP, len(rows) * th + (len(rows) - 1) * GAP), (232, 229, 222))
for r, row in enumerate(rows):
    for ci, t in enumerate(row): strip.paste(t, (ci * (tw + GAP), r * (th + GAP)))
strip.save(os.path.join(OUT, 'zoom4.png'))
F = pickle.load(open(os.path.join(HERE, 'final2.pkl'), 'rb')); C2 = pickle.load(open(os.path.join(HERE, 'compose2.pkl'), 'rb'))
stats = dict(kept=len(F['kept_ids']), total=len(C2['kept_recs']) + len(C2['dropped']), notes=F['notes'])
json.dump(stats, open(os.path.join(OUT, 'stats.json'), 'w'), indent=1); print(stats['kept'], stats['total'])
print({f: Image.open(os.path.join(OUT, f)).size for f in os.listdir(OUT) if f.endswith('.png')})
