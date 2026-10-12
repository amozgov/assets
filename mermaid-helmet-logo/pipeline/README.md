# Pipeline

The scripts and data that produced Steps 1 to 4. The folder mirrors the original session layout, so every script runs unchanged:

- `scratchpad/art/`: `lib.py` (drawing kit, `svgpath()` writes shapely geometry as SVG path data) and `render.js` (SVG to PNG).
- `scratchpad/study/`: the scripts, their data (`*.pkl.xz`) and `redraw_v-clean.svg` (Step 1, read by the diagnostics).
- `scratchpad/logo.svg`: the original logo. `images/3.png`: the lightning mermaid screenshot used for alignment.

## Setup

```sh
cd pipeline
xz -dk scratchpad/study/*.pkl.xz
python3 -m venv venv
venv/bin/pip install shapely==2.2.0 numpy==2.5.4 scipy==1.18.1 scikit-image==0.26.0 pillow==12.3.0
```

The `.pkl` files are Python pickles of numpy arrays and shapely geometries, written by these scripts. They need shapely 2.x to load.

`render.js` needs Node with Playwright's Chromium. It writes `file.png` next to `file.svg`, sized by the SVG's `width` and `height` attributes:

```sh
NODE_PATH=/path/to/global/node_modules node scratchpad/art/render.js file.svg
```

## Rebuild Steps 2 to 4

About ten seconds in total. Run from `scratchpad/study`:

```sh
../../venv/bin/python assemble5.py      # Step 2: final2.pkl, step2-clean.svg, step2-construction.svg
../../venv/bin/python step34.py         # Steps 3 and 4: step3.svg, step4.svg, step3.pkl, step4.pkl
../../venv/bin/python export_layers.py  # step4-layers.svg (groups figure / trident-hand / trident)
```

Run from a fresh copy of this folder, these reproduce `../svg/step2-*.svg`, `step3-mirrored.svg`, `step4-trident.svg` and `step4-trident-layers.svg` byte for byte.

## Data flow

```
logo.svg ─ vector.py ─> step1.pkl, redraw_v-clean.svg                      Step 1: exact medial-axis re-draw
3.png + logo.svg ─ align.py ─> align.pkl ─ prep.py ─> prep.pkl ─ align2.py ─> align2.pkl
step1 + prep + align2 ─ compose2.py ─> compose2.pkl ─ hands.py ─> compose3.pkl ─ spikes.py ─> compose4.pkl
compose4 ─ handfit.py ─> handfit.pkl
all of the above ─ assemble5.py ─> final2.pkl ─ step34.py ─> step3/step4 ─ export_layers.py
```

The upstream scripts only need re-running if their inputs change. `compose2.py`, `hands.py` and `spikes.py` also write older Step 2 SVGs as a side effect, so run `assemble5.py` after them.

`step1.pkl` holds, per stroke of the original: the medial axis (`axis`), its width profile (`prof`), primitives (`prims`: line / arc / bezier / poly), `width` and `caps`; the classified stroke polygons (`geoms`); the exact ink polygon (`ink`); and the five fills. `compose4.pkl` holds the new arm and hand pieces read off the lightning raster.

## What assemble5.py does

1. Splits Step 1's strokes into kept and removed. Removed: both old arms and hands, the spear, and the ids in `FORCE_DROP`: 1430 (the old armpit thorn), 320 (the fold's old stub), 540 (the bun knot's stump) and 880 (the neck's old continuation).
2. Protects every kept stroke that the removed ink touched with its exact band: a ribbon through the original ink's own edge points along the stroke's axis. Where one edge is inflated by removed ink, only that edge is carried across, as the smooth continuation of its clean course on either side (`band_info`).
3. Cuts each kept stroke square where its ink is clean, neither swollen by a junction nor tapering (`clean_len`), reads the band's edge points on the cut line (`ink_cut`), and starts the new piece exactly there.
4. Builds the joins: the bun arc, forehead mitre, breast fold (two mitres), left armpit, right armpit corner with the bush, throat L and the neck band with the shoulder cut. It also replaces the breast point's junction stubs with one dot.
5. Draws the arms (smooth single curves) and hands (one junction and four cubics each, from `handfit.pkl`), then composes: Step 1 ink minus the removed strokes, plus the new ink.

Its printed notes give each join's position and angle.

## Diagnostics

The first four take a window in drawing units and write `name.png` next to the script; the last two print to the terminal.

| Command | Shows |
| --- | --- |
| `s12.py name x y w h [W]` | Step 1 and Step 2 side by side |
| `dbg.py name x y w h [W]` | Step 1 ink, Step 2 outline, stroke axes with ids (blue kept, red removed) |
| `ccrop.py name x y w h [W]` | the Step 2 construction view |
| `hz.py name x y w h [W]` | the aligned lightning raster beside Step 2 |
| `notches.py [AMIN DMAX CMIN]` | steps and spikes in Step 2's outline that Step 1's outline does not have |
| `spikes_check.py [AMAX]` | sharp outline features, tiny parts and holes |

Defaults for `notches.py` are 10° turns within 0.12 u and corners over 100°. The last validation used `notches.py 8 0.2 150`: everything it reports is a thorn tip, a finger-gap tip or a hand junction.

`board4.py`, `board4_html.py` and `board34_html.py` built the canvas boards. The HTML scripts write into the session's artifact folder and need the uploaded asset ids in `board4/blobs.json` and `board34/blobs.json`, so treat them as templates.
