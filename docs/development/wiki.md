# GitHub Wiki migration package

The live Wiki was inspected at commit
`715e5039c75e080814a12e957f5148c35cdf8bda`, unchanged since the Phase 1B snapshot.
Its seven pages are Home, Matrix Class, Vector Class, Quaternion Class, Plane
Class, Ray Class and Common Functions. Matrix/Vector contain substantive historical
material; four class/function pages are placeholders. All original source and
attribution remain in the [audited snapshot](../../audit/wiki-snapshot/manifest.json).

## Prepared changes, not published

The [migration package](../wiki-migration/README.md) replaces Home with a concise
navigation entry, adds seven companion pages and a sidebar, and preserves the
original Home verbatim as Historical Home. Existing class/function pages are left
unchanged. No historical page is deleted, truncated or silently corrected. The
new Historical and Legacy Notes page explains that old formulas/examples are
historical evidence rather than current contracts.

Canonical docs and root ROADMAP remain authoritative. Wiki summaries use GitHub
source links; no unverified website URL is inserted. The connected repository tools
do not offer a Wiki write operation, and no CLI write identity was configured in
the observed environment. A reviewed package is provided rather than claiming a
successful Wiki publication.

## Review and publish separately

Clone the Wiki, inspect the [prepared pages](../wiki-migration/pages/Home.md),
and run the guarded preparation tool from the main repository:

```sh
git clone https://github.com/AlexMarinescu/pyGameMath.wiki.git /tmp/pyGameMath-wiki
python tools/prepare_wiki.py --wiki-dir /tmp/pyGameMath-wiki
python tools/prepare_wiki.py --wiki-dir /tmp/pyGameMath-wiki --apply
git -C /tmp/pyGameMath-wiki diff
```

Dry run makes no writes. Apply requires the inspected tip, a clean worktree and
no conflicting target pages; it only copies the reviewed files. If the Wiki has
changed, refresh the inventory and review conflicts instead of overriding the guard.
The tool never commits, pushes, deletes pages or changes remotes.

After maintainer review, commit those specific pages with an ordinary descriptive
message and push normally using authorized Wiki credentials. Do not force push.
Verify live pages and links afterward, then record the actual published commit.
Website deployment is a separate [publication decision](website.md).
