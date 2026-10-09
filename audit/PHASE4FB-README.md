# README introduction and visual refresh

Base: latest master `ce8e6883d1498c2768cff4485d2498929f087512`.

The landing page now starts with movement, rotation, curves and lighting instead
of technical categories. The existing HDR pair and four SVG diagrams appear near
the top, with short practical captions and links to larger diagrams and source.
Feature tables become use-case bullets; the two longer Python demonstrations
become one six-line movement example. The README is 850 whitespace-delimited words
instead of 1,107, and 139 lines instead of 198. These counts describe size, not a
measured comprehension or reading-time claim.

Installation uses the environment's Python directly, avoiding activation and a
separate import check in the opening instructions. The full installation guide
now labels Linux/macOS and PowerShell sections and supplies Command Prompt commands.
Windows/macOS instructions are documented, not newly verified platform support.

Release status, six, current source versus historical PyPI code, unverified Python
2.7 support, ownership rules, canonical/legacy imports, author and BSD attribution
are retained. The old claim that no documentation site exists is replaced with an
link to the local website build and navigation instructions. The existing Pages
publishing workflow is retained; public deployment availability was not verified. Proposed geometry,
noise, distance fields, grids, volumetric lighting and language projects remain
clearly unreleased. Root ROADMAP remains authoritative and unchanged.

## Executed checks

- Full suite: **2,260 passed**, zero failures, expected failures or skips (11.06s).
- Source and isolated installed-package documentation checks: **42 executable
  Python blocks, 565 source links, 268 declarations and seven constants** checked.
  README contains one Python example; it ran against the installed development
  wheel, with gem's installed path asserted.
- Independent quick-start answer: (1,2,3)+(4,0,-1)=(5,2,2). Original Vector/list
  data remains unchanged and result storage is independent.
- Strict MkDocs build and emitted-site checker: **65 HTML pages, 6,674 local links,
  15 image references, 16 MathML expressions**; no missing references/external assets.
- README preview: all six images and local links resolved. Chromium
  151.0.7922.173 inspected 1200×1000 desktop and 390×844 mobile layouts in light/dark
  appearance. Original image aspect ratios were preserved; no horizontal page
  overflow was detected. The preview used GitHub-compatible Markdown and limited
  inline HTML with local GitHub-style CSS, **not the live GitHub renderer**.
- All 93 protected pre-existing files remain byte-identical. Mathematical source,
  runtime package, dependencies, release metadata, license, tests, benchmark data
  and generated example images are unchanged.

Commands from the checkout:

```sh
/workspace/.venvs/pyGameMath/bin/python -m pytest -q
/tmp/phase4f-docs-env/bin/python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase4fb-source.json
/tmp/phase4f-wheel-env/bin/python -I tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --package-root /tmp/phase4f-wheel-env/lib/python3.12/site-packages --output /tmp/phase4fb-installed.json
/tmp/phase4f-docs-env/bin/python -m mkdocs build --strict
/tmp/phase4f-docs-env/bin/python tools/check_site.py --output /tmp/phase4fb-site.json
```

The local visual preview rendered README with Python Markdown's fenced-code/tables
extensions, preserved the actual image and anchor HTML, and checked each relative
href/src against the repository. A loopback-only Chromium session loaded all six
images, verified natural/rendered aspect ratios and page width, and saved the four
appearance/viewport screenshots under `/tmp/phase4fb-browser`. This was presentation
verification, with no generated artwork, replaced mathematical image or remote asset.

Changed files are README, the installation guide and this report. Nothing was
published, deployed, released or merged by this phase. No new API decision is needed.
