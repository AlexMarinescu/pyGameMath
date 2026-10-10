<div id="non-destructive-wiki-migration" aria-hidden="true"></div>

# Wiki portal update

The current Wiki already has its historical-preserving navigation portal. This
package proposes a focused update to Home, Documentation and the sidebar, using
the documentation website listed in repository homepage metadata and retaining
GitHub source fallbacks.

The manifest pins the observed Wiki commit and every page hash. The guarded
[preparation tool](../../tools/prepare_wiki.py) supports this reviewed replacement
without removing historical pages. See the [Wiki workflow](../development/wiki.md)
for exact dry-run and apply commands, current observations and publication limits.

The remaining files under `pages/` preserve the original migration package's
historical Home and companion content. They are not copied by this update's
three-file manifest. Main documentation remains authoritative.
