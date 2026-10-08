# Prepared Phase 1 review units

GitHub API access is blocked by the environment's egress proxy (CONNECT 403 for api.github.com; `gh api` Forbidden). No remote PR was created, and no branch was pushed just to bypass that blocker. Native Git reads from the existing checkout succeed. Opening drafts requires allowing `api.github.com` in environment network settings and then rechecking access with the existing authentication; a new token is not established as necessary.

Two local commits separate the executable test baseline from documentation and measured benchmarks. They can be reviewed as stacked draft PRs on the existing separate audit branch, or as one tests-and-audit PR with two reviewable commits. None contains library corrections. Do not merge or publish as part of this audit.

## Draft 1: Add Phase 1 mathematical regression and property tests

The previous CI launcher prints examples without checking answers. Add 271 pytest cases covering every core and experimental source module, independent Fraction matrix oracles, deterministic algebraic properties, known geometry/rotation/projection answers, and edge cases. Keep reproduced defects as strict xfails so `--runxfail` exposes the unchanged source's failures and future fixes cannot silently leave stale markers.

Validation on Python 3.12.14: 179 passed / 92 xfailed normally; 179 passed / 92 failed with `--runxfail`, exit 1 as expected. Audit tools are development-only requirements. Library source, runtime requirements, namespace, and mathematical conventions are unchanged.

Review scope: `tests/`, `requirements-audit.txt`.

## Draft 2: Document conventions, correctness findings, and measured benchmarks

Document 44 reproduced finding groups with source locations and severity, distinguish contract-dependent expectations and source-review-only issues, and record matrix row-vector conventions, quaternion ordering, mixed angle units, projection depth, and coordinate/SH sign inconsistencies. Include a reproducible standard-library benchmark harness, 14 operation baselines with raw samples, and a compatibility assessment.

Validation: 82.24% statement / 62.99% branch coverage; local wheel builds but omits experimental modules. Benchmark harness completed five timed repetitions per operation against the unchanged source. No fixes, API changes, modernization implementation, remote publication, or merges.

Review scope: `audit/`, `benchmarks/`. Stack after Draft 1 because test-results and reproduction instructions reference its tests.

## Future review sequence (not implemented)

After Phase 1 findings are accepted, prioritize separate small correctness PRs: equality; inverse2; quaternion inverse; paired division/inverse4 correction; angle helpers; refraction; plane representation/normalization; projection/transform composition; remaining quaternion/ray defects; then usable experimental algorithms. Contract-dependent changes need explicit decisions. Propose packaging/CI/documentation/typing/performance modernization only after the correctness phase is established.
