"""write the Study 2 board html from board4/blobs.json"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, 'board4', 'blobs.json'))); S = json.load(open(os.path.join(HERE, 'board4', 'stats.json')))
def fig(src, alt, title, cap):
    return f'''    <figure style="margin: 0; display: flex; flex-direction: column; gap: 8px; background: #ffffff; padding: 18px 18px 14px; border-radius: 6px">
      <div style="height: 520px; display: flex; align-items: center; justify-content: center"><img src="{src}" alt="{alt}" style="max-height: 520px; max-width: 100%; display: block"></div>
      <figcaption style="font-size: 12px; line-height: 1.45; color: #4d4942"><strong style="color: #1c1b19">{title}</strong> · {cap}</figcaption>
    </figure>
'''
FIGS = (fig(B['source'], 'The lightning mermaid screenshot, used only to read the arm pose', 'Pose reference', 'the screenshot, used only to read where the arms, hands and thorns go. Nothing is traced from it any more.')
        + fig(B['step1'], 'Step 1: the original mermaid re-drawn 1:1', 'Step 1', 'the clean 1:1 re-draw. Head, torso, breast, tail and fin of Step 2 are taken from here unchanged.')
        + fig(B['step2'], 'Step 2: the Step 1 drawing re-posed with raised arms, drawn by the rules', 'Step 2', 'Step 1 with the old arms, hands and spear removed, and the raised arms, hands, thorns and the five joins drawn by the Step 1 rules.')
        + fig(B['construction'], 'Construction: black Step 1 lines kept, orange new lines, purple thorns, teal circles at the rebuilt joins', 'Construction', 'black = Step 1 lines kept, <span style="color: #c96a16">orange</span> = new lines, <span style="color: #6b34a6">purple</span> = thorns, <span style="color: #1f8f80">teal</span> = the joins rebuilt.'))
html = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Study 2: lightning mermaid</title>
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
    <span style="font-size: 14px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase">Study 2 · the lightning pose, built from the Step 1 vector (v5, exact joins)</span>
    <span style="font-size: 12px; letter-spacing: 0.1em; color: #5f5a51">STEP 2 · NOT FOR PRINT</span>
  </div>
  <div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 24px">
{FIGS}  </div>
  <figure style="margin: 0; display: flex; flex-direction: column; gap: 8px; background: #ffffff; padding: 14px 18px 12px; border-radius: 6px">
    <img src="{B['zoom']}" alt="Eight details, Step 1 beside Step 2: bun, left shoulder, left armpit, breast fold, right shoulder, right armpit with thorns, left hand, right hand" style="width: 100%; height: auto; display: block">
    <figcaption style="font-size: 12px; line-height: 1.45; color: #4d4942"><strong style="color: #1c1b19">Details at 100 px per unit</strong> · Step 1 on the left, Step 2 on the right of each pair; the window size is in the label. Where Step 1 had the spear or the old arms, Step 2 shows the rebuilt join.</figcaption>
  </figure>
  <div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 24px; font-size: 12px; line-height: 1.5; color: #4d4942">
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">Kept from Step 1</span>
      <span>{S['kept']} of {S['total']} lines, as they are: the whole head and hair, neck, breast with its fold, nub and inner lines, torso, tail, scales and fin.</span>
      <span>Removed: both old arms and hands, the spear shaft, its blade and barb, the old armpit thorn, the fold's old stub and the knot stumps. Every kept line that the removed ink touched keeps its own two edges, read off the original ink sample by sample and carried straight across the old junction, so nothing of the old ink survives and nothing of the kept line is nibbled.</span>
      <span>The bun: the knot where the shaft entered is gone; the outline is one arc (R 2.7) tangent to the clean outline above and below it, its width running from the one side's to the other's; the strand that ran under the shaft stops under that arc.</span>
    </div>
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">Drawn by the rules</span>
      <span>Arms: two parallel standard lines (1.38), one smooth curve each, no kinks between segments, square ends.</span>
      <span>Hands: one light loop (0.9) around the palm; inside it four single curves leaving one point, so their square caps bury each other; their ends run into the loop.</span>
      <span>Armpit hair: one filled shape seated in the arm line, rising from the crease corner, with three straight-sided thorns 3.4 / 2.8 / 2.3 units long and 1-unit valleys between them. Breast point: Step 1's knot of three stubs is one round dot (R 0.9). The only filled marks beside Step 1's own.</span>
    </div>
    <div style="display: flex; flex-direction: column; gap: 4px">
      <span style="font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: #1c1b19; margin-bottom: 4px">The joins</span>
      <span>Every join is cut square across the kept line where its ink is clean, and the new piece starts exactly at the ink's own edge points on that cut. Left shoulder: the back-of-neck line is carried straight down past the old join as one band with a square end; the arm's top line runs into it and stops flush inside it.</span><span>Forehead: the hair outline turns into the face line in a mitred corner where the old arm used to cover the join.</span>
      <span>Left armpit: Step 1's light chest line runs on into the new arm line; the arm line's end face is flush with the chest line's outer edge.</span>
      <span>Breast fold: the contour turns up a short upright at its own weight and into the light chest line, two mitred corners (47° and 44°); one weight change, at the fold corner; the chest line runs on past the dot as one band.</span>
      <span>Right shoulder: throat line and arm line meet in a mitred L.</span>
      <span>Right armpit: the contour's own two edges run straight on to one mitred corner with the arm's lower line, which leaves the corner straight for 2.8 units before its curve; the hair bush rises from the crease corner.</span>
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
out = os.path.join(HERE, '..', 'artifact-files', '64ae74cb-f239-4963-a2e9-fc125fe724cf', 'project', 'Study-2-Lightning.dc.html')
open(out, 'w').write(html); print('wrote', len(html))
