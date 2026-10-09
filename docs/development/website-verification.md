# Executed website verification

The branch starts at merged PR #43, master
`4243ab0644806b9fa6ba7629e6ad6a3afcafc192`. The website was built locally and
inspected in Chromium; no website, Wiki or package was published.

## Results and states

| Deliverable/check | Executed result |
| --- | --- |
| MkDocs Material website | Strict build passes; 65 HTML pages including 404 |
| Built HTML validation | Local links/fragments, assets, navigation and search pass; no external asset dependencies |
| Search | 42 API index entries; browser search for `rotate_coefficients` returns the actual SH reference |
| Mathematics and code | 16 native MathML expressions; matrices, vectors, quaternion and SH formulas inspected beside highlighted Python |
| Browser | Chromium 151.0.7922.173; 1440×1000 desktop and 390×844 mobile |
| Themes and gallery | Light/dark screenshots inspected; all images load with original aspect ratios |
| Mobile and hosting prefix | Drawer sections and installation link work under `/pyGameMath/`; no horizontal page overflow |
| Browser resources | Zero external requests, missing resources or JavaScript errors |
| Basic accessibility | Language, main landmark, skip link, image alt text, named controls, keyboard search and focus styling checked |
| Mathematical pytest suite | 2,260 passed; 0 failed, 0 xfailed, 0 skipped (12.07s) |
| Phase 4A–4E plus website source checks | 43 Python blocks; 268 declarations and seven constants; all 93 protected files unchanged |
| Wheel/sdist | Clean builds and isolated installed documentation examples pass; all 17 wheel runtime files match source bytes |
| Reproducible examples | Showcase outputs and both HDR PNG hashes match; no canonical assets modified |
| Local preview | MkDocs serve started on loopback; HTTP 200 and homepage content verified |
| GitHub Pages workflow | Prepared manual master-only artifact build; not dispatched; no deployment job or permission changes |
| Wiki | Live inventory read; ten-file migration package tested on a local clone; no remote writes |

The [machine-readable report](website-verification.json) records final link counts,
source checks, browser measurements and publication state. This is not a full
WCAG certification or cross-browser/support-matrix claim. Screen-reader behavior
and other target browsers require platform testing.

## Executed commands

Documentation tools were installed into `/tmp/phase4f-docs-env` with the versions
now pinned in `requirements-docs.txt`. From the repository:

```sh
/tmp/phase4f-docs-env/bin/python -m mkdocs build --strict
/tmp/phase4f-docs-env/bin/python tools/check_site.py --site-dir site --output /tmp/phase4f-site.json
/tmp/phase4f-docs-env/bin/python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase4f-source-docs.json
python tools/check_site_browser.py --site-dir site --browser-executable /usr/bin/chromium --output-dir /tmp/phase4f-browser-evidence
/workspace/.venvs/pyGameMath/bin/python -m pytest -q
/workspace/.venvs/pyGameMath/bin/python -m examples.showcase.regenerate --verify
```

Browser QA serves a loopback-only static site, exercises real theme/search/drawer
controls and checks a project-prefix mount. Screenshots inspected locally:
`home-light.png`, `home-dark.png`, `notation-dark.png`, `gallery-dark.png`,
`home-mobile.png`, `mobile-drawer.png`. They are review evidence, not new mathematical
reference images, and are not required to regenerate the site.

The source was copied into a clean local checkout at `/tmp/phase4f-source`; strict
build and emitted-site validation passed there. Distribution checks used:

```sh
/tmp/phase4f-wheel-env/bin/python setup.py sdist bdist_wheel
/tmp/phase4f-wheel-env/bin/python -m pip install --no-index dist/gem-0.1.12-py3-none-any.whl
/tmp/phase4f-wheel-env/bin/python -I /workspace/pyGameMath/tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --package-root /tmp/phase4f-wheel-env/lib/python3.12/site-packages --output /tmp/phase4f-wheel-docs.json
/tmp/phase4f-sdist-env/bin/python -I /workspace/pyGameMath/tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --package-root /tmp/phase4f-sdist-env/lib/python3.12/site-packages --output /tmp/phase4f-sdist-docs.json
```

The sdist was extracted and installed into a separate clean venv with six and build
tools from local wheels. Isolated wheel execution added only copied example folders
to `sys.path`, asserted gem's installed location, regenerated SVGs to
`/tmp/phase4f-installed-gallery` and HDR outputs to `/tmp/phase4f-installed-hdr`.
No source-tree gem path was added. Existing PNG SHA-256 values remain
`dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a` and
`c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7`.

## Wiki and publication boundaries

`tools/prepare_wiki.py` dry-run and apply were executed against separate local
clones of the inspected Wiki tip. All six class/function pages remain byte-identical;
original Home is preserved verbatim as Historical Home. The live Wiki was untouched.
The connected tools expose no Wiki writer and the environment had no configured
CLI write identity, so the reviewed package is the deliverable. See
[Wiki publishing instructions](wiki.md).

GitHub Pages settings, custom domains, repository permissions and deployment were
not changed. A maintainer must approve publication and the separate deployment
configuration after review. See [build/preview commands](website.md) and the
[stack/license decision](documentation-stack.md), including Material's upstream
MkDocs 2.0 advisory. The verified pinned MkDocs 1.6.1 build remains BSD-licensed;
future tooling upgrades need renewed compatibility and license review.

The runtime, algorithms, signatures, numerical behavior, runtime dependencies,
release metadata and gem license remain unchanged. Root ROADMAP gains only the
long-term 2.0/cross-language direction and is consumed directly by the website;
there is no separately maintained site roadmap copy.
