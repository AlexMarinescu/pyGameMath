# Documentation stack and maintenance

The website uses **MkDocs 1.6.1** with **Material for MkDocs 9.7.6**. These pinned
versions were installed and built on CPython 3.12.14. MkDocs provides static files,
a local preview server and project-page-compatible relative links; Material adds
search, responsive navigation, light/dark palettes and accessible theme controls.
Both use permissive BSD/MIT licensing (MkDocs BSD, Material MIT). Existing gem
attribution and its BSD 2-Clause license remain unchanged. Theme assets retain
upstream notices. No paid theme extensions are required.

This is a documentation-only selection, not a runtime or interpreter-support
change. `requirements-docs.txt` pins the complete resolved build environment.
Review dependency updates together with strict builds, link/search/browser checks
and the mathematical documentation checks; do not automatically float major
versions. Re-evaluate compatibility before a future MkDocs/theme major upgrade.
Material 9.7.6 emits an upstream advisory about potential MkDocs 2.0 compatibility
and licensing changes. That advisory describes the upstream future direction, not
the installed BSD-licensed MkDocs 1.6.1. This build deliberately stays pinned to
the verified 1.x stack; upgrades need fresh maintenance and license review.
A static Markdown site was preferred to a JavaScript application or custom renderer
because it preserves the reviewed source documents and has a small build surface.

## Offline equations and assets

Pymdown arithmatex identifies inline/display TeX. A local build hook converts it
with latex2mathml 3.81.0 into native MathML. Equations therefore render without a
MathJax/KaTeX CDN, client-side math script or downloaded fonts. Current Chromium,
Firefox and Safari support MathML; older browser engines may display it differently.
Review matrix and equation output in target browsers rather than claiming universal
presentation parity. Conversion errors fail the build. MathML also exposes equation
structure to assistive technologies; actual screen-reader behavior needs platform testing.

Theme font downloads are disabled. Text uses system fonts. Search, theme scripts,
styles and gallery assets are served locally. There are no analytics, web fonts,
remote image embeds or online rendering services. The gallery keeps its original
light diagram backgrounds in either palette for contrast.

## Source of truth

`docs/` and root `ROADMAP.md` remain canonical. `tools/site/hooks.py` presents the
roadmap, benchmark guides and existing example READMEs as virtual build pages,
and copies external gallery assets once into build output. It rewrites source
links to actual repository files where appropriate, without editing their content.
No generated roadmap copy is committed. See [building and publication](website.md)
and [Wiki migration](wiki.md).

## Third-party notices

The built site includes the upstream [Material MIT license](../assets/licenses/material.txt)
and [MkDocs BSD license](../assets/licenses/mkdocs.txt), plus the theme's
[Material icon notice](../assets/licenses/material-icons.txt) and
[Font Awesome notice](../assets/licenses/fontawesome.txt). Theme controls and the
GitHub link use bundled icons; they are not an invented gem logo. No external
font files are loaded. These notices apply to documentation presentation assets,
not a relicensing of gem or its mathematical implementations.
