# Phase 4F-C — Documentation and website UX

Base: merged PR #59, `3e714fe949b7a6b7724d5c0da3395ee92483265f`.
Branch: `docs/phase4fc-visual-overhaul`.

## Summary

The MkDocs site has a yellow/white wordmark, consistent light/dark styling,
grouped desktop tabs, clearer module entry points and practical learning paths.
All **268 public API declarations, nine tutorials and 42 original executable
examples** are preserved. A new landing-page example brings executed examples to
43. Production mathematics, public signatures, metadata, dependencies and all
existing library tests/examples/benchmarks are unchanged.

## Design and learning

- Yellow `py` and white `GameMath` share a dark header in both themes. System fonts,
  measured text contrast, spacing, code-block surfaces and wrapping API signatures
  improve readability without remote fonts or decorative effects.
- Desktop tabs separate Getting Started, Tutorials, API Reference, Examples,
  Conventions and Development. Tutorial prerequisites group foundations, graphics
  and lighting; API topics retain their existing page paths.
- The landing page supplies installation, an independently checked Vector example,
  module cards, three learning paths and project/release status. The two HDR sphere
  images are paired with captions describing the unchanged camera/material/display
  settings and active +90-degree Z lighting rotation.
- Existing transformation, quaternion and Bezier diagrams retain their verified
  assets. New Bezier/Legendre equations join the notation guide. The Legendre plot
  evaluates `gem.legendre` and checks three analytical Fraction values; its source
  and full-size zoom link are included. No planned spline or other API is implied.
- All nine tutorials gain concise comments without changing executable syntax
  trees. Shared conventions identify the current repaired master. TeX becomes
  native MathML through the existing Arithmatex/latex2mathml hook, entirely offline.

The notation page contains **24 rendered MathML expressions**. Bezier tutorial
and associated-Legendre reference equations also render through this setup. The
site retains source/reconstruction, coefficient signs, row-vector composition,
quaternion ordering and numerical-limit explanations.

## Validation improvements

The former checker exits 2 on unmodified current master with
`out-of-scope modification: gem/bezier.py`: its byte freeze used the obsolete
Phase 3E checkpoint `dc418923c692b609a9fc611c66433e5950e0a321`.

The guard now pins merged PR #59 explicitly rather than following mutable HEAD.
All previous declaration, constant/reexport, runtime-signature, link, example and
byte-preservation checks remain. Additional checks reject unexpected protected
files, changes to original executable examples and loss of old section anchors.
Only the new documentation-tool regression module is allowed into the protected
paths. Original tests and mathematical examples remain immutable.

Source checks preserve **219 baseline heading anchors**, with compatibility aliases
for renamed Wiki/publication headings. Twelve new deterministic negative/ownership
cases check modified runtime/metadata/tests, unexpected additions, changed/deleted
examples, anchor aliases and guarded Wiki replacement. Pages deployment now runs
the complete documentation checker before upload, with its existing branch and
permission constraints retained.

## Executed verification

| Check | Result |
| --- | --- |
| Focused documentation-tool tests | **15 passed**, zero failed/xfail/skipped |
| Full library suite | **3,959 passed**, zero failed/xfail/skipped/errors |
| Runtime API coverage | **268 declarations**, seven constants plus compatibility identities |
| Source and installed examples | **43 executed** per environment; original 42 ASTs preserved |
| Tutorials and historical anchors | **9 tutorials**, **219 anchors** preserved |
| Strict MkDocs build | **66 HTML pages**, successful exit |
| Emitted links and resources | All internal links/fragments checked; zero external asset dependencies |
| Browser layout sweep | **56 combinations**: seven routes × four widths × two themes |
| Browser interactions | Search, theme toggle, mobile drawer, project-prefix navigation and keyboard skip link pass |
| Requests | Zero missing resources, JavaScript errors or external requests |
| Legendre plot regeneration | Byte-identical SVG with three independent analytical checks |

Widths are 320, 390, 768 and 1440 pixels. Routes include home, Vector/Quaternion/
Legendre API tables, notation, Bezier tutorial and the gallery, under `/pyGameMath/`.
Equations scroll locally, images load and no tested page exceeds viewport width.
Screenshots were inspected; mobile drawer branding was corrected before the final
pass. This is targeted usability/accessibility verification, not a full WCAG
certification or a claim about every browser/platform.

Measured contrast: ordinary links **6.73:1 light / 9.79:1 dark**; header yellow
`py` **11.97:1**, white `GameMath` **15.74:1**. Paragraphs, captions and primary
buttons also meet the checked 4.5:1 floor in both settled themes. Palette
transitions settle before measurement.

The clean wheel/sdist environments from Phase 4G-5R were reused for isolated
`-I` execution of the current documentation checker. All **17 installed runtime
files match current source byte for byte**; neither merged PR #59 nor this phase
changes those artifact bytes. Both environments pass all 43 examples and API
signatures with their installed `gem`, without importing runtime from the checkout.
No packaging rebuild or metadata migration is needed for documentation changes.

Environment: CPython 3.12.14, Linux x86_64/glibc 2.41, pytest 9.1.1, six 1.17.0;
MkDocs 1.6.1, Material 9.7.7, pymdown-extensions 12.1, latex2mathml 3.81.0;
Playwright 1.62.0 with Chromium 151.0.7922.173. Documentation/browser tooling
remains optional and separate from the installed core. Other interpreter/browser
versions are not newly certified.

Machine-readable counts, fingerprints, routes and observations are in
[phase4fc-verification.json](phase4fc-verification.json). Browser screenshots
include the [light homepage](phase4fc-screenshots/home-light-viewport.png),
[dark homepage](phase4fc-screenshots/home-dark-viewport.png) and
[mobile drawer](phase4fc-screenshots/mobile-drawer.png). Image hashes supplement
DOM, resource, contrast, mathematical-example and library tests; they are not
portable byte-only correctness oracles.

## Wiki and publication boundaries

The live Wiki was cloned at `400e42712494836ed016a155f2ab7202ae61ba0e` and already
contains the published historical-preserving portal. The proposed update changes
only **Home, Documentation and the sidebar**, linking the configured documentation
website and retaining GitHub source fallbacks. Guarded dry-run/apply succeeds in a
local clone; all **13 other pages** match their observed hashes. The v2 manifest
pins every current page and permits only reviewed replacements; old v1 guard
behavior remains covered by existing tests. No historical material is deleted.

This Wiki update is **prepared, not pushed**. It can be reviewed with the PR and
published separately after review; the currently published portal remains intact.
The website address comes from repository homepage metadata. Live availability
could not be verified: the environment HTTP proxy rejects
`alexmarinescu.github.io` with **CONNECT 403**. No proxy bypass was attempted and
no website outage is inferred. Local root/project-prefix routes pass independently.
No website deployment, Pages-setting change, PR merge or release publication was
performed.

## Reproduction commands

Run from the checkout. A separate documentation environment uses the unchanged
pinned requirements; Chromium is needed only for optional browser QA.

```sh
python3.12 -m venv .venv-docs
.venv-docs/bin/python -m pip install -r requirements-docs.txt
.venv-docs/bin/python -m pip install -r requirements-docs-browser.txt
.venv-docs/bin/python -m pip install -e .
.venv-docs/bin/python -m mkdocs build --strict --site-dir /tmp/phase4fc-site
.venv-docs/bin/python tools/check_site.py --site-dir /tmp/phase4fc-site --output /tmp/phase4fc-site.json
.venv-docs/bin/python tools/check_site_browser.py --site-dir /tmp/phase4fc-site --browser-executable /usr/bin/chromium --output-dir /tmp/phase4fc-browser
python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase4fc-docs.json
python -m pytest tests/test_documentation_overhaul.py tests/test_documentation_infrastructure.py -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4fc-focused.xml
python -m pytest -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4fc-full.xml
python tools/site/legendre_diagram.py
```

Executed build/browser interpreter: `/tmp/phase4f-docs-env/bin/python`; test/source
checker interpreter: `/workspace/.venvs/pyGameMath/bin/python`. Browser invocation
required sandbox network permission for its loopback server. Exact isolated
installed-check commands:

```sh
/tmp/phase4g5r-final-wheel-env/bin/python -I tools/check_architecture_docs.py --package-root /tmp/phase4g5r-final-wheel-env/lib/python3.12/site-packages --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase4fc-wheel-docs.json
/tmp/phase4g5r-final-sdist-env/bin/python -I tools/check_architecture_docs.py --package-root /tmp/phase4g5r-final-sdist-env/lib/python3.12/site-packages --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase4fc-sdist-docs.json
python tools/prepare_wiki.py --wiki-dir /tmp/phase4fc-wiki
python tools/prepare_wiki.py --wiki-dir /tmp/phase4fc-wiki --apply
git -C /tmp/phase4fc-wiki diff
python tools/site/collect_verification.py --output audit/phase4fc-verification.json
git diff --check
```

The collector uses these executed `/tmp/phase4fc-*` results, checks original file
scope and installed runtime bytes, regenerates the plot and confirms Wiki hashes.
Temporary log locations and screenshots can be recreated; exact image bytes may
vary with browser, font or platform differences.
