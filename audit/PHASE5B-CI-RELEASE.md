# Phase 5B — CI and release safeguards

Base: master `deea819172cf86f6f7aabeaa2693af1e33afb1b7` (PR #61).
Branch: `release/phase5b-ci`. Verification results are recorded below when completed.

## Changes

Reuse Phase 5A checks through a portable harness. Expand the matrix to CPython
3.10–3.14 on Ubuntu, Windows and macOS; retain failure diagnostics. Replace
floating action references with reviewed full-SHA pins. Consolidate strict
MkDocs/API/tutorial/site checking before Pages deployment, with write/OIDC
permissions limited to the deployment job. No Wiki writes occur.

The manual candidate workflow validates canonical master/SHA/version/tag,
clean-source provenance, archive contents and hashes, and isolated installation.
Protected review requires server-side reviewers, prevented self-review, no
administrator bypass and protected-branch restrictions. Missing protection
fails closed. Repository protection settings are not changed by this phase.

GitHub-first release preparation is separate from PyPI authority. The GitHub
workflow defaults to validation-only, consumes the exact previously validated
candidate assets, and never rebuilds them. Publication additionally requires
`GEM_RELEASE_ENABLED=true`, exact owner authorization and protected review.
No publication mode is executed in Phase 5B. PyPI remains disabled; current
publishing authority is unverified and blocks PyPI upload only.

## Verification

Full interpreter, platform, archive and documentation evidence is maintained in
`phase5b-verification.json`. Workflow configuration alone does not establish
platform support. Manual master-only dispatch and real protected approval cannot
be exercised on an unmerged PR; deterministic negative tests validate those
boundaries without creating tags, releases or changing server protections.

## Release checklist

See [CI and release guide](../docs/development/ci-release.md) for exact commands,
review protection requirements, GitHub asset installation, failure handling and
PyPI authority prerequisites. Phase 5C validates the final release commit;
Phase 5D publication requires separate explicit owner approval. Preserve the
verified official archives for eventual PyPI upload, without silent rebuilding.
