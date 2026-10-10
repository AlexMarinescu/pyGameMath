# CI and release-candidate review

Ordinary pull requests and master pushes validate the package; they never upload
it. The matrix covers CPython 3.10–3.14 on Ubuntu, Windows and macOS. A configured
job is not evidence that a platform passed: consult its completed run and the
[verification report](../../audit/PHASE5B-CI-RELEASE.md).

## Reproduce validation

Install CI tools separately from gem's runtime dependency:

```sh
python -m pip install -r requirements-ci.txt
python tools/check_workflows.py
python tools/ci_validate.py --output-dir /tmp/gem-ci --artifacts-dir dist
```

Use an empty artifact directory and an output directory outside the checkout.
The portable harness runs the complete suite, checks 268 API declarations and
43 executable examples, builds both distributions, runs strict Twine checks,
and installs each archive into a fresh isolated environment. It checks installed
runtime bytes and repeats the API/examples checks. Skips and expected failures
fail the candidate gate. Diagnostics remain available when a matrix job fails.
After the full matrix passes, a second installation matrix verifies the same
canonical Ubuntu/CPython 3.12 wheel and sdist on all 15 environments, without
rebuilding them. It reuses the shared installation checks and verifies manifest
provenance and SHA-256 digests before installation.

The archive inspector rejects unexpected wheel contents, links, unsafe paths,
caches, temporary credential files and recognizable credential-shaped content.
This leak heuristic cannot prove the absence of every possible secret.

## Manual candidate workflow

Run **Release candidate validation** from master, supplying its full current
commit SHA, version `1.0.0`, proposed tag `v1.0.0` and mode `validate` or `review`.
The supplied SHA must equal the workflow's dispatched master SHA. A missing
proposed tag is allowed; no tag is created. An existing tag must resolve to the
same commit. Forks and other branches cannot enter candidate validation.

`validate` runs the complete matrix and verifies the resulting candidate bundle.
`review` additionally requires the `gem-release-review` environment to have:

- Required human reviewers with self-review prevented.
- Administrative approval bypass disabled.
- Deployment limited to protected branches, with master protected server-side.

Configure these protections in repository settings before using review mode.
The preflight reads and validates the server configuration; a missing or
inaccessible protection fails closed. It does not create an environment or
change repository settings. Human approval occurs only after the matrix and
artifact integrity checks pass. Approval still enables **no publication**.

Candidate artifacts contain the wheel, sdist, `SHA256SUMS`, installation instructions
and a manifest with source commit,
tree and SHA-256 digests. Only artifacts from the same workflow run are consumed.
Download digest verification is mandatory; the separate gate verifies archive
hashes, clean-source provenance, exact version and contents against the checked-out
candidate. Hashes detect corruption; a hash supplied by an untrusted producer
is not an independent authenticity guarantee.

## GitHub-first release policy and publication checklist

GitHub Release v1.0.0 preparation is independent of PyPI account recovery.
Publishing that release and creating/pushing its tag require final validation
and explicit owner authorization for each action. Validation mode performs neither.
Preserve the approved wheel and sdist unchanged for eventual PyPI publication;
never silently rebuild official 1.0.0 artifacts.

Current authority to publish the historical **gem** PyPI distribution remains
unverified. Repository access does not establish PyPI authority. Before a future,
separately reviewed publishing workflow is enabled, establish:

- The authorized PyPI Owner/Maintainer account and its authority for `gem`.
- Maintainer approval of the exact artifacts, version and release commit.
- Account security requirements, including PyPI's applicable two-factor policy.
- An explicitly authorized trusted-publisher binding for this repository,
  workflow and protected environment, if OIDC is selected.
- Protected branch/review settings and an approved release/tag procedure.

No PyPI token is required or hardcoded. There is no PyPI/TestPyPI upload action
or OIDC package publishing permission. The GitHub publish job below remains
separately gated; ordinary CI and validation mode create no tag or release.
Package identity must not be changed to bypass ownership verification.

If validation fails, retain diagnostics, discard the candidate and rerun from a
reviewed clean commit after correction. Never upload a partially validated
bundle. If publication is enabled in a future phase, its failure/rollback policy
must account for immutable PyPI versions; do not assume deleting and reusing a
version is safe.

## Documentation publication

PR previews and Pages use the same strict documentation build, API/tutorial
checks and emitted-site link/anchor checks. Only the master deployment job has
Pages write/OIDC permissions. Documentation publication is independent of
package release. None of these workflows writes the GitHub Wiki.

Action versions are pinned to reviewed full commit SHAs. Refresh pins through a
reviewed PR with action metadata and syntax verification, not floating tags.

## Phase 5D GitHub publication workflow

**GitHub release preparation** defaults to `validate`. Supply the final master
SHA and successful manual candidate-validation run ID. The workflow downloads
those exact assets rather than rebuilding, verifies run provenance and hashes,
checks the proposed tag, and rejects an existing release. Archive retention is
14 days; if expired, validate a new candidate before official publication.
Preserve official assets durably once released.

Phase 5B tests preparation only. Phase 5C performs final validation. In Phase 5D,
separate explicit owner authorization must precede setting repository variable
`GEM_RELEASE_ENABLED=true` and dispatching mode `publish` with the exact
phrase `publish v1.0.0`. The protected environment requires human approval.
Only this publish job has contents-write permission. It rechecks provenance and
protections before creating the tag/release at the exact reviewed commit and
attaching the original archives, checksums and installation instructions.
Prepared [release notes](release-notes-1.0.0.md) describe compatibility and limits.

Leave the switch unset/false during preparation. A dry run creates no tag or
release. If a future authorized upload fails partway, inspect the existing
release/tag and assets before any recovery; this workflow refuses overwrites.
Never silently replace official artifacts. PyPI remains separately disabled,
regardless of GitHub publication approval.
