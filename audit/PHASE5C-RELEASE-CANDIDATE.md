# Phase 5C — Release candidate verification

Base: master `a21273c857e6ed9e948fa92246226100efb80355` (PR #62).
Branch: `release/phase5c-candidate`. Target: gem 1.0.0, GitHub-first.

## Windows reference investigation

An unchanged-baseline trace identifies the first difference in `math.sin` of
`0x1.921fb54442d18p-1` (the binary64 approximation of pi/4). Linux/macOS return
`0x1.6a09e667f3bccp-1`; Windows returns `0x1.6a09e667f3bcdp-1`, one ULP higher.
A 90-digit Decimal series using the exact input float gives
0.707106781186547502751942956217516746261543239537492789524366119137482021518043926217856922.
The Linux/macOS result is nearest; Windows differs by one ULP. Inputs and all
other traced math calls match on the representative CPython 3.12 jobs.

This changes SLERP's scalar weight for t=0.25, then the quaternion's scalar
component from `0x1.ee8dd4748bf15p-1` to `0x1.ee8dd4748bf16p-1`. The Hamilton
sandwich computes X as w*w-z*z: the rotated axis changes from
`0x1.bb67ae8584cabp-1` (0.8660254037844387) to
`0x1.bb67ae8584cadp-1` (0.8660254037844389), two ULP apart.
The mathematical result sqrt(3)/2 is
0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093;
its nearest float is `0x1.bb67ae8584caap-1`. The two results are respectively
one and three ULP above that reference. Absolute error remains below 2^-51.
The same sine call affects the t=0.75 imaginary component and the small
cancellation residual in the 90-degree rotated axes.

This is a benign platform libm difference and an overly strict full-precision
JSON serialization requirement, not a defect in the frozen mathematics.
Verification now permits differences no larger than 2^-52 between only eight
audited quaternion measurement fields. Both manifests must independently satisfy
closed-form values within 2^-51. SVG bytes, all other fields, key/type structure,
serialization, metadata and artifact hashes remain exact. Existing numerical
assertions and committed references are unchanged. Regression tests replay the
specific libm difference and reject larger, unrelated and shared numerical errors.

## Cross-platform verification

[Machine-readable verification](phase5c-release-verification.json) records every
job, interpreter/platform/compiler, command, artifact digest and installation
result. [Package run 38016538984](https://github.com/AlexMarinescu/pyGameMath/actions/runs/38016538984)
passed all 30 jobs: 15 full validations followed by 15 installations of the same
canonical archives. [Documentation run 38016538948](https://github.com/AlexMarinescu/pyGameMath/actions/runs/38016538948)
also passed.

| CPython | Linux | Windows | macOS | Canonical wheel/sdist installation |
| --- | --- | --- | --- | --- |
| 3.10 | 4050 passed | 4050 passed | 4050 passed | All three platforms passed |
| 3.11 | 4050 passed | 4050 passed | 4050 passed | All three platforms passed |
| 3.12 | 4050 passed | 4050 passed | 4050 passed | All three platforms passed |
| 3.13 | 4050 passed | 4050 passed | 4050 passed | All three platforms passed |
| 3.14 | 4050 passed | 4050 passed | 4050 passed | All three platforms passed |

Every full suite has zero failures, errors, expected failures and skips. All
jobs checked 268 declarations and 43 executable documentation examples.
Each primary job built and installed its own wheel and sdist; the second matrix
installed the identical canonical wheel and sdist on all 15 environments.
These are 60 isolated installations in total. Every installed runtime matched
all 18 source Python files and passed independent mathematical smoke checks.

The unchanged baseline trace run
[38015993880](https://github.com/AlexMarinescu/pyGameMath/actions/runs/38015993880)
reproduced the original failure on all five Windows interpreters, with 4035
passes and one failure each; Linux/macOS passed 4036. The same first sine
variation occurs on every Windows version. Intermediate policy validation run
38016288200 passed 4047 tests on every platform before the final three workflow
guard cases were added. No tests were skipped or converted to expected failures.

Representative CPython 3.12 trace environments were Linux x86-64, glibc 2.39,
GCC 13.3.0 / Python 3.12.15; Windows Server 2025 build 26100, AMD64,
MSC v.1943 / Python 3.12.10; and macOS 26.6.2 ARM64, Clang 13.0.0 /
Python 3.12.10. These are tested runners, not claims about every OS version.

Local CPython 3.12.14 passed the final 4050-test suite and both clean archive
installations. Four other local interpreters passed the preceding 4047-test
implementation. The hosted final matrix supplies complete final-version evidence.
The safeguard suite now has 74 tests, including the eight permanent tampering
cases and canonical-install dependency/digest/matrix guards. Existing mathematical
assertions, fixtures, APIs and all 18 runtime files remain unchanged from PR #62.

Strict MkDocs/site checking passed: 69 pages, 7724 links/anchors, 24 MathML
equations, 16 image references and seven navigation checks. All nine tutorials
and 43 examples remain present. Actionlint 1.7.12 and workflow policy checks
passed; Shellcheck was unavailable and was not executed.

## Preserved candidate and provenance

Verified branch implementation: `fb71f4feb2adcf20533afdeab961fb0db0a4a499`.
The hosted PR checkout and preserved candidate source commit is
`775fcafddcc6c930aa24e0ba9d963b9f8885c10d`, tree
`51cc45700514e8c4f11fad7b0a67c0ed882b3b9d`.
This is the tested merge ref of that branch into master
`a21273c857e6ed9e948fa92246226100efb80355`, not a final approved master release.
The following report/notes commit is separate from the tested source snapshot.

The canonical artifact is `gem-candidate-38016538984`, artifact ID 11656447664,
ZIP SHA-256 `3aef3e7232fe82ddb10cbacd419a91260e82393886ef12735b1841925126b4bd`.
It is retained by Actions until 2026-10-24 02:20:09 UTC. An unchanged copy is
preserved outside the source checkout at
`/workspace/release-candidates/phase5c-fb71f4f/` with its manifest,
`SHA256SUMS` and `INSTALL.md`.

| File | SHA-256 |
| --- | --- |
| gem-1.0.0-py3-none-any.whl | `ee06ecc4ca27ea5934fed62149bf2cc4a62c029c33abcd3f059da1df39b32c95` |
| gem-1.0.0.tar.gz | `b63ca118bbe3f409db042a28952374b45d5c92c4e0ffccf6f7bbfef795b42027` |
| SHA256SUMS | `ac228c890608bc8c87b67b57b29f7b5c82d2f59d5dd8402490dcd0aa7e2a309c` |
| INSTALL.md | `ff0a41a420db7bb41a40772703cd5b360555167fb7188985aa6ac6e3c138c91e` |

The manifest records clean-source provenance, distribution `gem`, version
`1.0.0` and exact archive hashes. The wheel has exactly 18 runtime modules plus
five metadata/license entries. Strict metadata, contents and Twine checks
passed. The sole mandatory dependency remains six; packaging metadata is intact.
A detached checkout of the recorded merge commit independently passed the
artifact integrity gate against the downloaded bundle.

Two local builds at the same branch commit produced identical member payloads
(23 wheel files, 393 sdist files), but different archive hashes: five wheel
metadata timestamps, 36 tar member mtimes and gzip timestamps differed.
Build reproduction is therefore content-equivalent, not byte-identical with
these existing build tools. No archive was rewritten to disguise this limitation.
Preserve the exact validated bytes; rebuilding cannot substitute for their
verified hashes. This phase changes no build metadata or archive format.

## Reproduce verification

Use a clean Git checkout with history and the CI-only requirements. Release
checks require that checkout; an extracted sdist is an installation artifact,
not a substitute for Git provenance and the audit baseline.

```sh
python -m pip install -r requirements-ci.txt
python -m pytest -q
python tools/trace_showcase.py --output /tmp/showcase-trace.json
python tools/check_workflows.py
python tools/ci_validate.py --output-dir /tmp/gem-ci --artifacts-dir /tmp/gem-candidate
python tools/ci_validate.py --output-dir /tmp/gem-exact-ci --candidate-directory /tmp/gem-candidate
python -m examples.showcase.regenerate --verify
python -m mkdocs build --strict --site-dir /tmp/gem-site
python tools/check_site.py --site-dir /tmp/gem-site --output /tmp/gem-site.json
```

Exact local candidate pipeline:

```sh
/tmp/phase5a-py312/bin/python tools/ci_validate.py --output-dir /tmp/phase5c-final-ci --artifacts-dir /tmp/phase5c-final-assets --offline --wheelhouse /tmp/phase5b-wheelhouse
/tmp/phase5a-py312/bin/python tools/ci_validate.py --output-dir /tmp/phase5c-exact-ci --candidate-directory /tmp/phase5c-final-assets --offline --wheelhouse /tmp/phase5b-wheelhouse
/tmp/phase5a-py312/bin/python /tmp/phase5c-verified-source/tools/release_gate.py verify-artifacts --directory /workspace/release-candidates/phase5c-fb71f4f --commit 775fcafddcc6c930aa24e0ba9d963b9f8885c10d --tag v1.0.0 --output /tmp/phase5c-canonical-integrity.json
```

The local build is separate from the hosted canonical build; their hashes are
not interchangeable. Machine-readable evidence distinguishes both.

## GitHub release recommendation: NOT READY for publication

Numerical and cross-platform candidate verification passed. CI-WIN01 is closed.
Publication still requires the following human-controlled prerequisites:

| Prerequisite | Status and action |
| --- | --- |
| Reviewed final master | PR #63 is unmerged. Review and merge, then run manual release validation at the exact resulting master SHA; the PR candidate cannot pass the publication provenance guard. |
| Protected release environment | Connector cannot read `gem-release-review`. Verify required reviewers, prevented self-review, disabled admin bypass and protected-branch deployment; exercise protected review successfully. |
| Release enablement | Connector cannot inspect `GEM_RELEASE_ENABLED`. Keep unset/false during preparation and verify its actual setting manually. |
| Explicit owner authorization | No authorization for tag creation or GitHub Release publication has been granted. Obtain it separately in Phase 5D. |

The public master ruleset is active and requires PRs/linear history, but has zero
required approving reviews and no required status-check rule. Administrator
bypass is present. The release-tag ruleset prevents updates/deletion/force pushes
and also has administrator bypass. These settings were inspected, not changed;
`protected: true` alone is not proof of enforced CI or protected human approval.
Review the effective merge policy and permissions before final signoff.

Two moderate default-branch Dependabot alerts were reported by Git push. Details
are inaccessible through the connector; owner triage remains a documented risk,
not a claim of a newly confirmed runtime vulnerability. The release environment,
branch-protection endpoint and enablement variable could not be fully inspected.
Missing release protection fails closed at runtime. No real manual master dispatch
or approval was simulated as a successful protected release.

Prepared release title: **gem 1.0.0**. See
[release notes](../docs/development/release-notes-1.0.0.md) and
[installation/release guide](../docs/development/ci-release.md). PyPI publication
is pending account recovery and verified publishing authority; this blocks PyPI
only. There is no claim that gem 1.0.0 is available on PyPI.

Final approval checklist:

- Review this PR and verify the environment, merge policy and outstanding alerts.
- Validate the exact final master SHA using the existing manual workflow;
  complete both matrices, integrity verification and protected review.
- Preserve that run's exact wheel/sdist/checksums/instructions before retention
  expires, and identify them explicitly in the publication approval.
- Obtain separate authorization for tag creation and GitHub Release publication;
  ensure v1.0.0 targets that same reviewed and tested commit.
- Publish only in the separately authorized Phase 5D. Preserve official bytes
  for later PyPI use; never silently rebuild or replace them.

No tags, GitHub Releases, package uploads, repository-setting changes or merges
were performed. The 1.0 mathematics, API, dependencies and feature scope remain frozen.
