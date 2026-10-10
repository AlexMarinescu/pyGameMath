<div id="github-wiki-migration-package" aria-hidden="true"></div>

# GitHub Wiki navigation portal

The live Wiki was inspected at commit
`400e42712494836ed016a155f2ab7202ae61ba0e`. It already contains the
lightweight navigation portal, six original class/function pages and Historical
Home. Historical Matrix/Vector material and attribution remain intact. The
[audited snapshot](../../audit/wiki-snapshot/manifest.json) preserves the earlier
source at `715e5039c75e080814a12e957f5148c35cdf8bda`.

<div id="prepared-changes-not-published" aria-hidden="true"></div>

## Focused portal update

The [reviewed package](../wiki-migration/README.md) changes only Home,
Documentation and the sidebar to link to the repository's configured
[documentation website](https://alexmarinescu.github.io/pyGameMath/).
Repository source fallbacks remain available. Current API contracts and tutorials
stay in the main documentation; the Wiki does not duplicate them.

The website address comes from repository homepage metadata. This execution
environment's HTTP proxy rejects the GitHub Pages domain with CONNECT 403, so
live website availability was not independently confirmed. The local build and
all matching project-prefix routes are verified separately.

<div id="review-and-publish-separately" aria-hidden="true"></div>

## Review and apply

```sh
git clone https://github.com/AlexMarinescu/pyGameMath.wiki.git /tmp/pyGameMath-wiki
python tools/prepare_wiki.py --wiki-dir /tmp/pyGameMath-wiki
python tools/prepare_wiki.py --wiki-dir /tmp/pyGameMath-wiki --apply
git -C /tmp/pyGameMath-wiki diff
```

Dry run makes no writes. Apply requires the inspected tip, a clean worktree and
matching hashes for every observed Wiki page. All reviewed sources are validated
before any write. Only the three portal files are replaced; all historical and
other companion pages are preserved. The tool never commits, pushes or deletes.
If the Wiki changes, refresh the inventory and review conflicts before applying.

This update is prepared for review alongside the documentation PR; no live Wiki
push is claimed. After review, commit only those portal changes, push normally
with authorized Wiki credentials, and verify the published links. No force push
or second API specification is needed.
