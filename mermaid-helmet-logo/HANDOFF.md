# Handoff: mermaid front logo for a Triple Eight helmet

State as of 2026-10-12. Everything referenced here is inside this folder, except the live design canvas (see "Canvas").

## The job

The user is making two vinyl stickers to personalise a Triple Eight skate helmet: a small round one over the front logo, and a wide patch over the "triple eight" wordmark on the back. The user said the back is solved. Work only on the front.

The front is a mermaid in the style of an existing artist's monoline mermaid logo (`reference/original-mermaid.svg`). The user wants that style ("a bit sarcastic, childish, yet stylish and very artistic") but not a 1:1 copy. The remix: the raised-arms pose of the same artist's lightning mermaid, mirrored, holding a trident, with a bush of armpit hair.

The agreed plan has four steps:

1. Re-draw the original 1:1 to learn its rules. **Done and accepted.**
2. Re-pose the Step 1 vector into the lightning mermaid's pose: raised arms, open hands, armpit bush. **Done.**
3. Mirror it, tail to the left. **Done.**
4. Put a trident in the hand that is on the left after the flip, in place of the lightning bolt. **Done, but rejected on two points (next section).**

## What to do next

The user's last feedback, verbatim:

> 1. Hand holding trident needs to be updated. How does mermaid hold it with open hand?
> 2. Trident needs a bit more sophistication it doesn't look "magnificent" now

Start from `svg/step4-trident-layers.svg`. It splits Step 4 into three groups: `figure` (everything except the trident hand), `trident-hand` and `trident`. The three groups reproduce `svg/step4-trident.svg` exactly.

### 1. A hand that grips

Today the trident hand is the lightning mermaid's open hand: a light loop (0.9 u wide) with four inner curves leaving one junction. The shaft is simply knocked out inside the loop, so it reads as passing behind an open palm. Its bounds are x 7.23 to 12.59, y −5.92 to 1.75.

It must become a closed fist around the shaft, built like the original's spear fist. Study that fist in `svg/step1-original-redraw.svg` in the window x 9.5 to 22.5, y 6.5 to 18.5. `reference/original-spear-head-fist-fletching.png` shows it with stroke ids and widths.

- The fist is one rounded loop (strokes 80, 120, 50 and 70, weights 1.06 to 1.23) around two white counters, one above the shaft and one below it.
- The shaft (strokes 40 and 130, about 1.44) runs straight through the loop and stays visible inside it, between the two counters (stroke 90).
- A single vertical cut at x ≈ 16.8 splits the top of the loop into a higher left part and a lower right part, so one side of the fist overlaps the other.
- The forearm's two lines (strokes 60 and 110) run into the bottom of the loop.

The stroke ids are Step 1's medial-axis pieces, not anatomy: study the shape in the picture, not the id list.

Keep the forearm lines as they are; in the `figure` group they end at the old loop's outer edge. Keep the fist roughly where the loop is, centred near (9.81, −1.97). The original's spear runs at 14.8° and passes through the fist sideways, while this shaft is vertical, so the fist has to be turned to match; decide how and show it.

### 2. A trident that looks magnificent

Today: a straight standard shaft (1.38) from 34 u below the fist to 6 u above it, a half-circle U of radius 2.6, three straight prongs 3.0 u long, plain filled triangle tips 2.3 wide and 2.8 high, the middle prong 0.8 higher. It reads as a garden fork (`previews/step4-trident-head.png`).

Directions that stay inside the style. Combine two or three and show the options side by side at the same scale:

- **Barbed tips** modelled on the original's spear head: a filled arrowhead with swept-back barbs and a small notch where the shaft enters (top panel of the reference image; in Step 1 coordinates x −1 to 12, y 8.5 to 17.5).
- **A dominant centre blade**, longer and wider than the side tips.
- **Side prongs that sweep outward and up** from the shaft, as a crescent or a bracket with mitred corners, instead of a plain U.
- **A collar** under the fork: two or three short bands across the shaft with knockout gaps between them.
- **A finished butt end** echoing the original's fletching, a ladder of slanted bars inside a pointed frame (bottom panel of the reference image), or a pointed ferrule.
- **Orientation.** The original spear crosses the figure at 14.8°; a diagonal shaft may read grander than the vertical one. The user has not commented on orientation, so offer it as an option rather than impose it.

Use only the three pen weights. Fill only tips and blades, as the original fills only its arrowhead. Mitre the corners, keep caps flat, and knock the shaft out wherever it passes behind the fist or the body.

### Then

- Validate before showing: run `notches.py` and look at every new join in zoom crops at 150 to 1300 px per unit (tools in `pipeline/README.md`).
- Publish a new board on the canvas to the right of "Study 3-4" (x 8920, y −2520, 1840 × 1400). Keep the old boards.
- After approval: the final round sticker layout and its size in mm. The earlier concept exports assumed a 50 mm tall die-cut; the final size is not decided.

## The rules

These were measured from the original, and the user enforces them strictly.

- **Three pen weights**, in drawing units: light 1.15, standard 1.38, heavy 1.60. In the original they make up 25, 58 and 17 % of the line length.
- **Flat ends**: square, or cut along the edge of the line they meet. **Mitred corners.** No fillets, no rounded joins, no tapering.
- **Filled shapes are rare and deliberate.** The original has five: the spear head, eye lens, nipple, navel and a hair swirl. Step 2 adds the armpit bush (one shape with three thorns, never separate spikes) and the breast-point dot.
- **Knockout gaps** of ½ to ¾ line width where a line passes behind another.
- **Exact primitives**: every line is a straight segment, a circular arc or one smooth curve. The original's 177 lines are 410 straight segments, 65 arcs and 56 curves. No wobble, steps, spurs or notches at any zoom.
- **Scales** are U-cups.
- **Proportions** of the original: head 8.7 × 10.4 u, figure 6.1 heads tall.

## How this user works

- They review at high zoom and circle defects on screenshots; `feedback/` has every one, oldest first. The two from 19:14 critique Step 1's first re-draw; the rest (21:46, 23:38, 00:11, 00:45) critique successive Step 2 rounds. Typical catches: a step where a new line meets an old one, a spike at a corner, a notch, a bump, thorns drawn as separate objects, a join that looks patched.
- "Everything needs to be perfect. This is not about zoom or not zoom."
- "Simple clean less": build each piece once from the rules. Never patch over a bad join; rebuild it.
- Never trace or fit a raster. Use rasters only to read a pose. The Step 1 vector is the ground truth for everything the original already has; the lightning mermaid is the same drawing at 1.008× scale, and only the arms, hands, armpit and bolt differ.
- Short iterations, visible progress, brief reports. They are often away. When told to continue alone, finish and validate a step before starting the next one, and never move past an unvalidated step.
- When a choice is genuinely theirs, ask once and briefly (example: which hand holds the trident; their answer was the hand on the left after the flip).

## Files

| Path | What |
| --- | --- |
| `svg/step1-original-redraw.svg` | Step 1: the original re-drawn 1:1 (IoU 0.999 with the original) |
| `svg/step1-construction.svg` | Step 1 construction: stroke centrelines coloured by weight over the ink |
| `svg/step2-lightning-pose.svg` | Step 2 final (v5, exact joins) |
| `svg/step2-construction.svg` | Step 2 construction: black = Step 1 lines kept, orange = new lines, purple = bush, teal = rebuilt joins |
| `svg/step3-mirrored.svg` | Step 3: Step 2 mirrored |
| `svg/step4-trident.svg` | Step 4: the trident version the user rejected |
| `svg/step4-trident-layers.svg` | Step 4 in groups `figure`, `trident-hand`, `trident`: start here |
| `svg/concepts/` | Earlier concepts, superseded: fronts A1 Trident Siren, A2 Wake Siren, A3 Reverse Mermaid, B1 Golden Snail, B2 Tongue Wave, B3 Eight Hands (each as circle and die-cut), back patches A and B, helmet mockups; `spec.json` holds their sizes in mm |
| `previews/` | PNG renders of everything above, plus Step 1 vs Step 2 detail crops |
| `reference/original-mermaid.svg` | The original logo as the user pasted it |
| `reference/lightning-mermaid-hires.webp` | The lightning mermaid, the pose reference (1724 × 2560, transparent) |
| `reference/original-spear-head-fist-fletching.png` | The original's spear head, fist and fletching, with stroke ids |
| `feedback/` | The user's annotated screenshots, named by UTC time |
| `pipeline/` | The code and data that produce Steps 1 to 4; see `pipeline/README.md` |

The four clean step SVGs are one compound path each, filled `#1c1b19` with `fill-rule="evenodd"`. The layered file holds three such paths, one per group.

## Geometry and coordinates

- Drawing units (u) are the original logo's: viewBox 0 0 54 64, y pointing down. A standard line is 1.38 u.
- Step 2 bounds: x 7.23 to 55.09, y −5.92 to 64.00.
- Step 3 is Step 2 mirrored about x = 31.1595, so x' = 62.319 − x.
- Step 4 adds the trident: bounds x 6.06 to 55.09, y −17.17 to 64.00.
- The trident hand is the lightning mermaid's right-hand-side hand before the flip, centred near (52.51, −1.97); after the flip it is near (9.81, −1.97). The free hand is on the right after the flip.

## Canvas

The design canvas is the user's private artifact https://claude.ai/artifact/DS6Bz4hEMiWgTeZvYU4b26. Its boards:

- **Artboard 1**: the user's moodboard (the helmet's front and back, the original mermaid, an engraved variant, the lightning mermaid, two sets of the artist's other work).
- **Reference · Original mermaid**, **Study · Original re-drawn 1:1**, **Study · Zoom proof**.
- **Study 2 · Lightning mermaid re-posed** (Step 2, v5).
- **Study 3-4 · Mirrored, with trident** (Steps 3 and 4, the rejected trident).

Each board is `project/<Name>.dc.html` (1840 × 1400, IBM Plex Mono, background `#e8e5de`, images uploaded as artifact assets and referenced as `/_blob/<id>`), listed in `project/canvas.json` with position, size, title and order. The user edits the canvas, so re-read `canvas.json` before changing it. The board generators are `pipeline/scratchpad/study/board4_html.py` and `board34_html.py`. The earlier concept boards are no longer on the canvas; their SVGs are in `svg/concepts/`.

## History in brief

- The user first saw six front concepts and two back patches, and liked A2, the helmeted Wake Siren, but found it messy. Their notes on it carry over where relevant: trident in the spare hand, breasts like the original, short distinct armpit hair, narrow elegant arms and hands, a curvier body, more sophisticated tail and scales.
- After several A2 iterations the user stopped it ("we need fresh start") and set Step 1. A first re-draw at IoU 0.967 was "very messy"; the accepted one is an exact medial-axis reconstruction at IoU 0.999.
- Step 2 began as a vectorisation of the lightning raster and was rejected as messy several times. It was rebuilt from the Step 1 vector: the original's lines kept, only the arms, hands and armpit drawn new, and every join cut where the kept line's own ink is clean.

## Known issues

- The Study 2 board's text gives the fold corners as 47° and 44°; the final geometry has 50° and 44°.
- The Study 3-4 board still presents the open-hand trident the user rejected.
