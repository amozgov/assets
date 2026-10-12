"""write the Study 3-4 board html from board34/blobs.json"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, 'board34', 'blobs.json')))
def fig(src, alt, title, cap, h=560):
    return f'''    <figure style="margin: 0; display: flex; flex-direction: column; gap: 8px; background: #ffffff; padding: 18px 18px 14px; border-radius: 6px">
      <div style="height: {h}px; display: flex; align-items: center; justify-content: center"><img src="{src}" alt="{alt}" style="max-height: {h}px; max-width: 100%; display: block"></div>
      <figcaption style="font-size: 12px; line-height: 1.45; color: #4d4942"><strong style="color: #1c1b19">{title}</strong> · {cap}</figcaption>
    </figure>
'''
FIGS = (fig(B['step2'], 'Step 2: the lightning pose built from the Step 1 vector', 'Step 2', 'the finished lightning pose: every join cut square where the kept line is clean and the new piece started at the ink\'s own edge points; no step, spur or notch left at any zoom.')
        + fig(B['step3'], 'Step 3: the Step 2 figure mirrored, tail to the left', 'Step 3 · mirrored', 'the same vector flipped about its own centre line: tail to the left, the sketch hand now on the left. Nothing redrawn.')
        + fig(B['step4'], 'Step 4: the mirrored figure holding a trident in the left hand', 'Step 4 · trident', 'a trident in the left hand, drawn by the rules in place of the lightning bolt: standard straight shaft, a U of radius 2.6 with three prongs and filled tips; the shaft passes behind the hand as the spear did in the original.')
        + fig(B['head'], 'The trident head and the fist at 50 px per unit', 'The head and the grip', 'the fork and the fist. The shaft runs up through the loop\'s outline and stops at it on both sides, so the hand reads in front of it; the side prongs leave the U tangentially; every tip is one filled triangle 2.3 wide, 2.8 high, with the prong buried in its base.'))
html = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Study 3-4: mirrored, trident</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&amp;display=swap" rel="stylesheet">
<style>
body{{margin:0}}
a{{color:#57534b}}a:hover{{color:#1c1b19}}
</style>
</helmet>
<div style="width: 1840px; height: 1400px; box-sizing: border-box; padding: 28px 40px 26px; display: flex; flex-direction: column; gap: 18px; background: #e8e5de; color: #1c1b19; font-family: 'IBM Plex Mono', ui-monospace, Menlo, monospace">
  <div style="display: flex; justify-content: space-between; align-items: baseline; gap: 16px">
    <span style="font-size: 14px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase">Study 3-4 · the finished pose mirrored, and the trident</span>
    <span style="font-size: 12px; letter-spacing: 0.1em; color: #5f5a51">STEPS 3 AND 4 · FRONT LOGO CANDIDATE</span>
  </div>
  <div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 24px">
{FIGS}  </div>
  <div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 24px; font-size: 12px; line-height: 1.5; color: #4d4942">
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">Step 2, checked</span>
      <span>The whole outline was swept for steps, spurs and notches with a turn detector (every pair of opposite turns within 0.2 units, every corner over 150°) against Step 1's own outline as the baseline, and every join was looked at between 150 and 1300 px per unit.</span>
      <span>What remains that Step 1 does not have: the thorn tips, the finger-gap tips of the hands and the junction of the four hand curves. All are drawn shapes, not artefacts. Tail, fin and scales are Step 1's to the last vertex.</span>
    </div>
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">Step 3, the mirror</span>
      <span>A flip about the figure's own centre line (x = 31.2). Line weights, curves and joins are untouched; the figure now faces right with the tail to the left, the sketch hand on the left.</span>
    </div>
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">Step 4, the trident</span>
      <span>Shaft: one standard line (1.38), straight up through the fist, 34 units below it and 6 above it to the fork. Fork: a half-circle U of radius 2.6 at the shaft's weight; the shaft runs on as the middle prong; two side prongs 3 units long rise from the U's ends, tangent to it.</span>
      <span>Tips: three filled isosceles triangles, 2.3 across and 2.8 high, the prongs buried 0.3 in their bases; the middle prong stands 0.8 higher. The shaft stops at the hand loop's outline on both sides: the hand is in front, the knockout the original used for the spear.</span>
      <span>The lightning bolt is dropped. Nothing else of Step 3 changes.</span>
    </div>
  </div>
</div>
</x-dc>
<script type="text/x-dc" data-dc-script data-props='{{"$preview":{{"width":1840,"height":1400}}}}'>
class Component extends DCLogic {{
  renderVals() {{
    return {{}};
  }}
}}
</script>
</body>
</html>
'''
out = os.path.join(HERE, '..', 'artifact-files', '64ae74cb-f239-4963-a2e9-fc125fe724cf', 'project', 'Study-3-4-Trident.dc.html')
open(out, 'w').write(html); print('wrote', len(html), out)
