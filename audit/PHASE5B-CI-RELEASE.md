# Phase 5B — CI and release safeguards

Base: master `deea819172cf86f6f7aabeaa2693af1e33afb1b7` (PR #61).
Branch: `release/phase5b-ci`. PR: [#62](https://github.com/AlexMarinescu/pyGameMath/pull/62).
Verified implementation: `0e085b4bfd532935bb1267bb3eae95f04de3ed15`.
Final report changes follow that commit; they do not change the tested implementation.

## Changes

Reuse Phase 5A checks through a portable harness. Run CPython 3.10–3.14 on
Ubuntu, Windows and macOS, retaining diagnostics even when the full suite fails.
Pin official Actions to reviewed full SHAs, disable persisted checkout credentials,
and keep default permissions read-only. Strict documentation checks precede
Pages deployment; only its deployment job receives Pages/OIDC permissions.
No workflow writes to the Wiki.

The manual release-validation workflow checks canonical repository, master,
commit, version, tag, clean-source provenance, archive contents and hashes.
Protected review requires reviewers, prevented self-review, no administrator
bypass and protected-branch restrictions. Missing protection fails closed.

The separate GitHub Release workflow defaults to validation-only and consumes
the exact assets from a successful manual validation run without rebuilding.
Publication additionally requires `GEM_RELEASE_ENABLED=true`, the exact owner
confirmation phrase and protected review. It checks remote tag binding before
attaching the wheel, sdist, SHA-256 checksums and installation instructions.
Existing releases cannot be overwritten. PyPI publishing remains disabled.

Portable changes address LF checkout behavior, fixture encoding, path separators,
virtual-environment invocation and benchmark environment metadata. Benchmark
workloads and comparison rules are unchanged. The documentation checker retains
its frozen API baseline with exact-hash allowances for three supporting files.
All 18 runtime Python files, packaging metadata, runtime dependencies and
reference assertions remain unchanged.

## Verification

[Machine-readable evidence](phase5b-verification.json) records interpreter and
platform details, exact commands, artifact digests, installed checks and job URLs.
Hosted [package run 38014899622](https://github.com/AlexMarinescu/pyGameMath/actions/runs/38014899622)
completed all 15 jobs at the verified implementation commit:

| CPython | Ubuntu passed / failed | macOS passed / failed | Windows passed / failed |
| --- | --- | --- | --- |
| 3.10 | 4036 / 0 | 4036 / 0 | 4035 / 1 |
| 3.11 | 4036 / 0 | 4036 / 0 | 4035 / 1 |
| 3.12 | 4036 / 0 | 4036 / 0 | 4035 / 1 |
| 3.13 | 4036 / 0 | 4036 / 0 | 4035 / 1 |
| 3.14 | 4036 / 0 | 4036 / 0 | 4035 / 1 |

Every job has zero expected failures and zero skips. All 15 jobs passed clean
builds, metadata/content checks, `twine check --strict`, isolated wheel and sdist
installation, mathematical smoke checks, 268 API declarations and 43 executable
documentation examples: 30 isolated installations in total. Each job tests its
own built artifacts. A failed suite still prevents candidate approval.

Local Linux verification passed 4036 tests on each of CPython 3.10.21, 3.11.16,
3.12.14, 3.13.15 and 3.14.7. The host was Linux 6.18.44, x86-64, glibc 2.41.
An initial interpreter lacking ensurepip was replaced with an installable CPython
3.13 environment before final verification.

The permanent release-safeguard suite has 70 passing tests, separate from
mathematical regressions. Eight deterministic tampering cases cover an extra
asset, dirty source, wrong commit, wrong tree, wrong version, modified digest,
modified checksum listing and empty installation instructions. Additional cases
exercise forks/events, protection settings, workflow permissions, failed test
results, archive paths, candidate provenance and remote tag binding. Publication
commands are mocked; validation-only execution performs no external mutations.

Actionlint 1.7.12 and the workflow policy checker pass. Shellcheck was unavailable
and was not run. Hosted documentation
[run 38014899541](https://github.com/AlexMarinescu/pyGameMath/actions/runs/38014899541)
passed. Local strict MkDocs/site checks validated 69 pages, 7724 local links and
anchors, 24 MathML equations, 16 image references, seven navigation checks and
42 API search entries, with no external asset dependencies. All nine tutorials
remain present.

### Windows blocker: CI-WIN01

The remaining failure on all five Windows interpreters is
`tests.test_showcase_examples::test_committed_artifacts_match_generation`.
The quaternion scene in `examples/showcase/scenes.py` generates a measurement
of `0.8660254037844389` where the committed JSON contains
`0.8660254037844387`. The strict dictionary comparison rejects this last-bit
difference. Numerical reference checks pass; SVG references match after the LF
checkout correction. This is a cross-platform exact-reference limitation, not
an independently confirmed mathematical defect.

The assertion and committed reference are intact: no rounding, skip, xfail or
relaxed tolerance was introduced. Windows support is not fully verified. The
full matrix is red and release approval remains blocked until a reviewed change
resolves the portable-reference policy. Initial Windows encoding/path failures
and macOS affinity metadata failures were corrected without weakening assertions.

### Reproduction commands

From a clean checkout, install `requirements-ci.txt`, then run:

```sh
python -m pytest -q
python -m pytest -q tests/test_release_workflows.py
python tools/check_workflows.py
python tools/ci_validate.py --output-dir /tmp/gem-ci --artifacts-dir /tmp/gem-assets
python -m mkdocs build --strict --site-dir /tmp/gem-site
python tools/check_site.py --site-dir /tmp/gem-site --output /tmp/gem-site.json
```

The final local candidate used:

```sh
/tmp/phase5a-py312/bin/python tools/ci_validate.py --output-dir /tmp/phase5b-complete --artifacts-dir /tmp/phase5b-complete-assets --offline --wheelhouse /tmp/phase5b-wheelhouse
/tmp/phase5a-py312/bin/python tools/release_gate.py verify-artifacts --directory /tmp/phase5b-complete-assets --commit 0e085b4bfd532935bb1267bb3eae95f04de3ed15 --tag v1.0.0 --output /tmp/phase5b-complete-integrity.json
/tmp/phase5b-tools/actionlint -color=false -shellcheck='' .github/workflows/*.yml
```

These preparation artifacts passed installation and integrity checks; they are
not official release assets:

| Artifact | SHA-256 |
| --- | --- |
| Wheel | `f4d6b096ab99c7c0c8f3b5bd71239a307ec0d467bc819d6d01e045f2c0b2d141` |
| Sdist | `230fca9becbcd0afb2b5e520eeae099fe60fc00213abf292d2bbcfdc1c134de8` |

## Release safeguards and remaining prerequisites

- CI-WIN01 blocks full-matrix release validation.
- Server-side release-environment configuration could not be inspected through
  the available connector. Guards reject missing protection; configuration and
  a real protected approval still require maintainer verification. Master reports
  protected, but that alone does not prove all required status policies.
- A moderate default-branch Dependabot alert remains untriaged in this review:
  [alert #2](https://github.com/AlexMarinescu/pyGameMath/security/dependabot/2).
  Alert details were inaccessible; no vulnerability-free claim is made.
- Publishing authority for historical `gem` remains unverified. This blocks PyPI
  upload only, not GitHub-first preparation or separately authorized publication.
- Real master-only release dispatch and protected approval were not executed on
  this unmerged branch. Repository protection/enablement settings were not changed.

See the [CI and release guide](../docs/development/ci-release.md) and
[prepared release notes](../docs/development/release-notes-1.0.0.md) for installation,
protection settings, artifact retention and failure handling. Phase 5C must
validate the final reviewed release commit. Phase 5D publication requires separate
explicit owner authorization. Preserve official verified archives for eventual
PyPI upload rather than silently rebuilding them.

No tags, GitHub Releases, PyPI/TestPyPI uploads, Wiki writes or merges were performed.
