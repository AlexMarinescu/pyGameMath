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

<div id="github-pages-preparation" aria-hidden="true"></div>

## GitHub Pages publication

The existing `publish-docs.yml` workflow builds and deploys documentation on
matching master pushes or manual dispatch. Its build runs only on master;
deployment uses the `github-pages` environment and the existing Pages/OIDC
permissions. Strict builds, emitted-link checks and the complete documentation
checker must pass before artifact upload. Feature-branch builds stay local.

The separate manual artifact workflow produces a preview for inspection without
deploying it. The site uses relative paths and is verified under `/pyGameMath/`.
The repository homepage identifies
[the documentation website](https://alexmarinescu.github.io/pyGameMath/).
Live availability could not be checked from this environment because its HTTP
proxy rejects that domain; this is not evidence of a website failure.
See [current verification](visual-overhaul.md) and the [Wiki portal update](wiki.md).
Normal builds never invoke `mkdocs gh-deploy` or alter repository Pages settings.

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
