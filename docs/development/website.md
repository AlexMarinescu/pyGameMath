# Build, preview and publish documentation

The site uses the [pinned documentation stack](documentation-stack.md). These
packages belong in a separate build environment; they are not gem dependencies.
From the repository root with Python 3.12:

```sh
python3.12 -m venv .venv-docs
.venv-docs/bin/python -m pip install -r requirements-docs.txt
.venv-docs/bin/python -m pip install -e .
.venv-docs/bin/python -m mkdocs build --strict
.venv-docs/bin/python tools/check_site.py --site-dir site
.venv-docs/bin/python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website
.venv-docs/bin/python -m examples.showcase.regenerate --verify
.venv-docs/bin/python -m mkdocs serve --dev-addr 127.0.0.1:8000
```

Open the preview at `http://127.0.0.1:8000/`; stop with Ctrl+C. On Windows substitute
`.venv-docs\Scripts\python.exe`. Windows commands are documented, not executed here.
Build emits `site/`, ignored by Git. It needs no NumPy, GPU, rendering package,
remote font or CDN. Installing dependencies is the only online setup step.
The build hook reads repository files directly; no manual copying is needed.

`tools/check_site.py` checks actual emitted HTML links, fragment targets, image/script
references, navigation, search entries and MathML. Browser review is separate;
see [executed verification](website-verification.md). The source checker retains
all earlier immutability, declaration and executable-block checks.

## GitHub Pages preparation

The proposed workflow builds and uploads a static-site artifact on manual dispatch
**only from master**. It has read-only repository permissions and no deploy job,
Pages settings changes, write token, custom domain or feature-branch publication.
The locally built output uses relative paths and was tested under a project prefix.
No production documentation URL is configured or claimed.

After review and merge, a maintainer can run the artifact workflow and inspect its
output. Enabling GitHub Pages then needs an explicit repository decision: choose
GitHub Actions as its source and add a separate reviewed deployment job with
`pages: write`, `id-token: write`, the `github-pages` protected environment and a
master-only condition. Confirm the actual resulting URL before adding `site_url`
or linking it from README/Wiki. Do not run `mkdocs gh-deploy` as part of normal builds.

Until then, documentation is browsable on GitHub and in local preview. The
[Wiki package](wiki.md) points to canonical repository documents and does not
assume a published website exists.

## Optional browser verification

Static builds require no browser dependency. For repeatable visual/search/mobile
QA, install `requirements-docs-browser.txt` in a separate environment with an
available Chromium executable, then run:

```sh
python tools/check_site_browser.py --site-dir site --browser-executable /usr/bin/chromium --output-dir /tmp/gem-site-browser
```

The tool serves only loopback, checks root and `/pyGameMath/` mounts, exercises
search and the mobile drawer, and records requests and screenshots. It fails on
external requests, missing resources or JavaScript errors. The Linux executable
path is a verification example, not a new cross-platform support claim.
