# Showcase asset manifest

This plan precedes the example implementations. All inputs are procedural Python
constants. No downloads, fonts, textures, models or image libraries are required.
Run `python -m examples.showcase.regenerate` from the checkout to rebuild the four
SVGs and their numerical manifest. Add `--verify` to compare a rebuild with the
committed artifacts without overwriting them. Output is confined to the named
files in `examples/showcase/output/`; `--output-dir` selects another directory.

SVGs use a 960 × 560 viewBox, CSS sRGB colors and text with a generic sans-serif
font. Geometry and labels are deterministic; font rasterization depends on the
viewer. SVGs and `measurements.json` are committed. There is no random sampling.
The existing HDR PNGs are reused by link, never copied or overwritten.

| ID / purpose | Procedural input | Exact output / format | Generation and numerical checks | Visual acceptance |
| --- | --- | --- | --- | --- |
| vectors / displacement and orientation | a=(3,1,0), b=(-1,2,0), origin | `examples/showcase/output/vectors.svg`, SVG 960×560 | Vector sum, normalization, dot and cross; sum=(2,3,0), dot=-1, cross=(0,0,7) | Labeled arrows, displaced b, unit a and positive-Z orientation readable |
| transforms / row-vector composition | four XY rectangle points, scale(2,1), +90° Z rotation, translation(3,-1) | `examples/showcase/output/transforms.svg`, SVG 960×560 | Matrix4 products; combined (x,y) → (3-y,2x-1); translation-before-rotation shown separately | Original and transformed silhouettes distinguished; order reversal visibly differs |
| quaternions / interpolated coordinate frames | identity and +120° Z unit quaternion, t=0,.25,.5,.75,1 | `examples/showcase/output/quaternions.svg`, SVG 960×560 | SLERP and quaternion vector rotation; frames at 0°,30°,60°,90°,120°; norm=1 | Five legible frames with increasing counterclockwise orientation |
| bezier / uniform and adaptive sampling | quadratic [(0,0),(2,4),(4,0)], cubic [(0,0),(1,4),(3,-2),(4,1)] | `examples/showcase/output/bezier.svg`, SVG 960×560 | Core evaluation and BezierPath at squared tolerances .16 and .0025; quadratic midpoint=(2,2), cubic midpoint=(2,.875), endpoints and dense-reference distance checks | Control polygons, sample markers and tolerance-related point counts distinguishable |
| lighting / active HDR SH rotation | existing procedural directional environment, canonical RGB L2 SH | `examples/output/sh_original.png` and `examples/output/sh_rotated.png`, PNG RGB8 sRGB 192×192; existing JSON/GLSL beside them | Existing `python -m examples.hdr_sh.regenerate --output-dir /tmp/gem-showcase-hdr`; core projection, analytical +90° Z rotation, one cosine convolution, Lambertian evaluation | Same camera/material/exposure; highlight moves right to up; no golden file changes |

`examples/showcase/output/measurements.json` records inputs, numerical results,
SVG dimensions, byte lengths and SHA-256 hashes. Hashes aid reproducibility;
independent numerical and structural checks provide correctness evidence.

No extra Plane/Ray/Viewport gallery is planned: the four diagrams and existing
lighting reference form a coherent progression without adding local intersection
algorithms or suggesting that the legacy viewport helper is a camera API.
