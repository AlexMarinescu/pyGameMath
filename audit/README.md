# Reproduce the Phase 1 audit

> This section records the original Phase 1 baseline. On `audit/phase1b-wiki`, use the updated counts and commands in the Phase 1B section below. Original raw Phase 1 results and benchmarks are retained.

Run from the existing `/workspace/pyGameMath` checkout. Keep its separate audit branch; no extra worktree is required. A fresh installation can use any chosen environment directory outside the source tree:

```bash
cd /workspace/pyGameMath
python3 -m venv /workspace/.venvs/pyGameMath-audit
source /workspace/.venvs/pyGameMath-audit/bin/activate
python -m pip install -e .
python -m pip install -r requirements-audit.txt
python -m pip check
```

The baseline used Python 3.12.14 with six 1.17.0. For matching the runtime dependency version use `python -m pip install six==1.17.0`. `requirements-audit.txt` pins direct audit tools, not the core package's runtime metadata or supported Python versions. Tests and benchmarks use no NumPy dependency. Coverage is optional development tooling; omit it when only running pytest.

## Tests

```bash
python -m pytest -q -p no:cacheprovider -o junit_family=legacy
```

Expected baseline: **179 passed, 92 xfailed**, exit 0. Each known failing regression has a finding ID and strict expected-failure marking. Unexpected passes fail the normal run, forcing review of the marker when a fix arrives.

To expose all confirmed baseline failures as real failures:

```bash
python -m pytest --runxfail -q --tb=short -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/pygamemath-audit-failures.xml
```

Expected baseline: **92 failed, 179 passed**, exit 1. No skips. This command is expected to fail because Phase 1 makes no library corrections. Preserve its exit status in automation; do not interpret a piped log viewer's status as pytest success.

Focused reproduction examples:

```bash
python -m pytest --runxfail tests/test_vector_common.py::test_equality_all_components -q
python -m pytest --runxfail tests/test_matrix.py::test_inverse2_known_answer -q
python -m pytest --runxfail tests/test_quaternion.py::test_inverse_identity -q
```

Each case uses small known answers or independent oracles. Random properties use local fixed seeds, not shared global state. Sample-generation tests temporarily fix `random.random` through pytest monkeypatch and restore it automatically. Irradiance fixtures are generated in pytest temporary directories, not fetched image assets.

## Coverage

```bash
COVERAGE_FILE=/tmp/pygamemath-audit.coverage python -m coverage run \
  --branch --source=gem -m pytest -q -p no:cacheprovider -o junit_family=legacy
COVERAGE_FILE=/tmp/pygamemath-audit.coverage python -m coverage report -m
```

Baseline statement coverage: 82.24%; branch coverage: 62.99%. A combined percentage of 78.34% is not branch coverage. Every source module is imported; tests cannot reach many deeper experimental paths due to earlier errors. No claim of exhaustive input coverage or formal proof is made.

## Benchmarks and package check

```bash
python benchmarks/baseline.py --output /tmp/pygamemath-benchmark-results.json
python -m pip wheel --no-deps . -w /tmp/pygamemath-audit-wheel
```

The benchmark warms each operation, automatically calibrates iteration count, then records five samples with GC disabled by timeit. It writes full metadata, iteration counts, samples, and medians. Do not overwrite the committed baseline merely to run a comparison. Run without coverage instrumentation and under similar machine load. See [BENCHMARKS.md](BENCHMARKS.md) for limitations.

The wheel build succeeded, but inspecting its ZIP entries found no `gem/experimental/` paths. This is a packaging defect distinct from editable-source functionality; packaging remains unchanged in Phase 1.

## Review artifacts

- [Prioritized source audit](PHASE1.md)
- [Observed mathematical conventions](CONVENTIONS.md)
- [Compatibility impact assessment](COMPATIBILITY.md)
- [Benchmark summary](BENCHMARKS.md) and [raw samples](benchmark-results.json)
- [Machine-readable test outcomes and module coverage](test-results.json)
- [Prepared draft PR descriptions](PULL_REQUESTS.md)

Reports refer to the unchanged library source at commit `5257291431bb45db0274dc48edf24694ecfe2e2d`. Regeneration of reports, future correctness fixes, broader interpreter testing, and remote PR creation must be recorded separately from these baseline observations.

## Phase 1B wiki validation

See [complete wiki/source report](PHASE1B-WIKI.md), [frozen seven-page wiki](wiki-snapshot/manifest.json), [API decision list and implementation order](PHASE2-DECISIONS.md), and [updated machine-readable results](phase1b-test-results.json). The library and benchmark baseline are unchanged.

Run the same pytest commands above from the Phase 1B branch. Current expected outcomes are 254 passed and 94 strict xfailed (86 confirmed-defect cases plus 8 contract questions). With `--runxfail`, expect 94 failed and 254 passed, exit 1. There are no skips. Use `python -m pytest tests/test_wiki_contracts.py -q -p no:cacheprovider -o junit_family=legacy` to run only the historical-documentation additions: 75 passed, 2 xfailed.

The marker registry and JUnit properties distinguish `defect(id)` from `contract_question(id)`; the latter are proposed policies awaiting review, not approved implementation requirements. Both kinds execute assertions. Questions should not be promoted to confirmed defect requirements merely because a normal run marks them xfailed.

For a fresh Phase 1B coverage run, use a separate path such as `COVERAGE_FILE=/tmp/pygamemath-phase1b.coverage`; do not overwrite retained baseline result files. Measured Phase 1B statement coverage is 83.80%, branch coverage 64.71% on Python 3.12.14.
