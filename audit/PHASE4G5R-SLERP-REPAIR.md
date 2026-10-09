# Phase 4G-5R — SLERP numerical repair

Base: `4e48dc68cc23d47abd7501897c6dfbb0c05156c3` (master, merged PR #58).
Branch: `repair/phase4g5-slerp-endpoints`.

4G5-A01 is repaired. All 12 original endpoint reproductions pass with their
inputs, assertions and tolerances unchanged. The full suite passes **3,947 tests**,
with zero failures, expected failures, errors or skips. The Bézier performance
concern remains outside this repair.

## Arithmetic and repair

Let `u = 2**-1074`, the smallest positive binary64 value. The public axis-angle
constructor can represent `[1,u,0,0]`. Its Hamilton sandwich rotates `[0,1,0]`
to `[0,1,2u]`; the transverse result is representable. The former SLERP angle
calculation obtains difference norm `u` and sum norm `2`, then calculates
`2*atan2(u,2)`. The intermediate half-angle rounds to zero before multiplication
by two. The zero-angle branch therefore copied the starting quaternion even
at `t=1`, losing available endpoint information. This is avoidable arithmetic
loss, distinct from an interior component such as `u/2` being unrepresentable.

Only `quat_slerp` changes in production:

- At `t=0/1`, copy the selected endpoint's components into fresh storage.
  Existing strict-negative-dot sign selection still chooses the shortest path;
  a zero dot product retains the supplied half-turn branch. Preserve signed
  zeros. The existing angle calculation precedes this selection, and nonfinite
  norm results retain the previous arithmetic path.
- If a nonzero separation produces a zero or subnormal spherical angle, use
  the continuous weighted linear limit for `t` in `[0,1]`. The boundary is the
  smallest **normal binary64** value, not a geometric equality epsilon.
- Blend component weights directly, retaining floating-point weight dispatch
  for already accepted Fraction parameters. Do not route through public LERP's
  narrower scalar-operand dispatch.

For unit inputs in this interval, spherical weights differ from the linear
limit by `O(theta**2)`. When theta is subnormal, that correction is far below
the binary64 grid. Dividing separately quantized subnormal sine values instead
can introduce a large relative weight error: an endpoint with component `3u`
at `t=3/8` previously produced `2u`, rather than the representable rounded `u`.
The limit avoids this additional loss without normalizing inputs or results.

`squad4` inherits the repair through its existing nested SLERP calls; its formula
and implementation are unchanged. Legacy three-control SQUAD, no-invert SLERP
and LERP are unchanged.

## Compatibility and limits

Signatures, return types, component order, shortest-path rotation conventions
and caller ownership are preserved. Endpoints return fresh Quaternions and
independent component lists. There is no new normalization, domain validation,
parameter clamping, dependency or quaternion representation.

Ordinary interior results are bit-identical in 960 matched SLERP/squad4 calls;
480 legacy-SQUAD calls and 320 ordinary endpoint calls also match. These include
observations of the existing nonunit behavior, not an expansion of its supported
domain. Twenty-five invalid/nonfinite observations retain their previous results
or exceptions. These checks characterize specific historical inputs, rather
than defining a new general invalid-input policy.

Separately rounded weighted products can still differ from a single-rounded
exact blend. For example, nested tiny-angle SQUAD products `1.25u` and `2.25u`
round separately to `u` and `2u`. The resulting `3u` is within their derived
one-ulp combined budget of the exact `3.5u`; it need not equal its single-rounded
`4u`. Tests distinguish this unavoidable evaluation rounding from angle-ratio
loss. Unrepresentable orientations are not recovered. Arbitrary extrapolation,
ill-conditioned/nonunit inputs and nonfinite-input policies are not redefined.

## Independent verification

The 82 additional cases use exact Fraction Hamilton expansions, dyadic
small-angle references and independently constructed great-circle rotations.
They cover minimum subnormals, adjacent minimum-normal boundaries, both signs,
nonidentity controls, signed-zero endpoint copies, antipodal representations,
half-turn ties, Fraction parameters, input storage, ordinary rotations and
rotation/matrix agreement. SQUAD references include a separately derived angle
polynomial, rather than only comparisons with nested gem calls.

| Executed check | Passed | Failed | Xfail | Skipped |
| --- | ---: | ---: | ---: | ---: |
| Unmodified master, full suite | 3,853 | 0 | 12 | 0 |
| Original reproduction with markers disabled | 0 | 12 | 0 | 0 |
| Original cases after repair | 12 | 0 | 0 | 0 |
| Focused quaternion/integration suites | 1,120 | 0 | 0 | 0 |
| Final full suite | 3,947 | 0 | 0 | 0 |

The reproduction failure is intentional evidence from the unmodified baseline.
All other listed final checks have zero errors. Full-suite execution took
16.34 seconds on this host; counts come from executed JUnit records.

Clean wheel and sdist installations each pass 12 endpoint cases, six additional
subnormal/Fraction checks, two independent rotation/matrix checks, 268 API
declarations plus seven constants/reexports, and **42 executable documentation
examples**. Source examples also pass. Both installations use isolated `-I`
execution outside the checkout and pass `pip check`.

All 17 runtime Python files in both artifacts match the final source byte for
byte. Only `gem/quaternion.py` changes among them. The sole declared dependency
remains `six`; there are no compiled extensions. Build-generated sdist
`setup.cfg` comment removal and `[egg_info]` settings are recorded separately
from unchanged source packaging metadata. The documentation checker remains
unchanged: verification reuses its API/signature helpers and runs the existing
snippets without its historical documentation-only scope gate.

Machine-readable evidence: [verification](phase4g5r-verification.json).
The collector checks original test-body preservation, production scope,
source/artifact fingerprints, installed results and benchmark source hashes.

## Matched performance

Three completed runs use the immutable master quaternion source and the final
repair in the same interpreter/process, with deterministic data, nine trials
and 20,000 calls per implementation per workload. Implementations alternate
order; workload order uses seed `475102+trial`. Each side has 200 warmup calls.
Process CPU time excludes setup/reference checks and includes call-loop,
allocation and wrapper costs; GC is disabled only during each timed batch.

The table shows run 3 medians and median absolute deviations (MAD), in
microseconds per call. Paired ratios are after/before medians for runs 1/2/3;
they need not equal ratios of unpaired median latencies.

| Workload | Before µs ± MAD | After µs ± MAD | Paired ratios, runs 1 / 2 / 3 |
| --- | ---: | ---: | --- |
| Raw ordinary SLERP | 2.948 ± 0.058 | 3.118 ± 0.147 | 1.012 / 1.028 / 1.072 |
| Method ordinary SLERP | 2.903 ± 0.046 | 2.863 ± 0.053 | 1.007 / 1.025 / 0.967 |
| Nearby rotations | 2.900 ± 0.138 | 2.839 ± 0.091 | 0.995 / 1.046 / 0.975 |
| Tiny normal separation | 2.753 ± 0.140 | 2.772 ± 0.073 | 1.053 / 1.036 / 0.998 |
| Identical controls | 2.343 ± 0.054 | 2.315 ± 0.050 | 0.993 / 1.059 / 0.993 |
| Antipodal representations | 2.700 ± 0.112 | 2.590 ± 0.075 | 1.037 / 1.016 / 0.940 |
| SLERP t=0 | 2.870 ± 0.042 | 2.472 ± 0.031 | 0.852 / 0.907 / 0.883 |
| SLERP t=1 | 2.937 ± 0.137 | 2.486 ± 0.066 | 0.869 / 0.873 / 0.850 |
| Minimum-subnormal endpoint | 2.824 ± 0.104 | 3.022 ± 0.171 | 1.027 / 1.017 / 1.032 |
| Subnormal interior blend | 3.552 ± 0.061 | 3.287 ± 0.071 | 0.911 / 0.907 / 0.925 |
| Ordinary squad4 | 8.966 ± 0.333 | 8.729 ± 0.153 | 0.985 / 1.003 / 1.018 |
| squad4 t=0 | 8.910 ± 0.243 | 7.603 ± 0.205 | 0.871 / 0.847 / 0.853 |
| squad4 t=1 | 8.757 ± 0.122 | 7.780 ± 0.281 | 0.886 / 0.842 / 0.876 |
| Unchanged legacy SQUAD control | 3.714 ± 0.048 | 3.610 ± 0.056 | 0.935 / 1.022 / 0.969 |

Endpoint paths consistently avoid sine-weight calculations: ordinary SLERP
endpoint latency is 9–15% lower and squad4 endpoint latency 11–16% lower in
these runs. Raw ordinary SLERP shows 1–7% paired overhead. Method SLERP and
ordinary squad4 changes vary around baseline; nearby/tiny normal cases also
vary. Added comparisons cost work, but shared-host variance prevents a precise
general overhead claim. The unchanged legacy control varies from −6.5% to
+2.2%, illustrating measurement noise. Minimum-subnormal endpoint comparisons
also compare a previously incorrect result with the corrected operation.

These are local measurements, not universal speed guarantees or timing gates.
Full trials, ranges, MADs and matched inputs are retained in
[run 1](phase4g5r-performance.json), [run 2](phase4g5r-performance-repeat.json)
and [run 3](phase4g5r-performance-third.json).

Environment: CPython 3.12.14, Linux x86_64/glibc 2.41, Intel Xeon Platinum 8370C
on a shared host; pytest 9.1.1, six 1.17.0. Build tools: setuptools 84.0.0,
wheel 0.48.0, packaging 26.3. Other Python versions were not verified.

## Commands

From the checkout, `PY` below denotes `/workspace/.venvs/pyGameMath/bin/python`.
The first two commands ran before the repair; the reproduction exits 1.

```sh
PY=/workspace/.venvs/pyGameMath/bin/python
$PY -m pytest -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5r-baseline.xml > /tmp/phase4g5r-baseline.log 2>&1
$PY -m pytest tests/test_final_integration_audit.py --runxfail -k subnormal_interpolation_endpoint -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5r-reproductions.xml > /tmp/phase4g5r-reproductions.log 2>&1
$PY -m pytest tests/test_final_integration_audit.py -k subnormal_interpolation_endpoint -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5r-original-passing.xml > /tmp/phase4g5r-original-passing.log 2>&1
$PY -m pytest tests/test_slerp_numerical_repair.py tests/test_final_integration_audit.py tests/test_quaternion_interpolation.py tests/test_vector_quaternion_optimization.py tests/test_quaternion_numerical_repairs.py tests/test_quaternion.py tests/test_quaternion_contracts.py tests/test_quaternion_operations.py tests/test_quaternion_matrix.py tests/test_quaternion_powers.py -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5r-focused.xml > /tmp/phase4g5r-focused.log 2>&1
$PY -m pytest -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5r-full.xml > /tmp/phase4g5r-full.log 2>&1
$PY benchmarks/slerp_repairs.py --output audit/phase4g5r-performance.json
$PY benchmarks/slerp_repairs.py --output audit/phase4g5r-performance-repeat.json
$PY benchmarks/slerp_repairs.py --output audit/phase4g5r-performance-third.json
$PY benchmarks/slerp_wheel_smoke.py --source --output /tmp/phase4g5r-source-smoke.json
```

Benchmark runs are sequential, with local test/build work finished before timing.
The comparison requires the base commit in the local Git object database.

Distribution staging copied only `gem/` (excluding caches), unchanged setup and
manifest files, README/license, legacy CI configuration and the two
manifest-listed migration guides into `/tmp/phase4g5r-final-dist-build`.
The build command ran there:

```sh
/tmp/phase4g5-build-env/bin/python setup.py sdist bdist_wheel
```

Clean environments were created with the same Python 3.12 interpreter. Cached
standard build wheels and six were installed offline; replace the cache location
with a locally prepared wheelhouse when reproducing elsewhere. Exact installation
and isolated check commands, run from `/tmp`:

```sh
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 -m venv /tmp/phase4g5r-final-wheel-env
/tmp/phase4g5r-final-wheel-env/bin/python -m pip install --no-index --find-links /tmp/phase4d-dependencies six
/tmp/phase4g5r-final-wheel-env/bin/python -m pip install --no-index --no-deps --no-build-isolation /tmp/phase4g5r-final-dist-build/dist/gem-0.1.12-py3-none-any.whl
/tmp/phase4g5r-final-wheel-env/bin/python -m pip check
/tmp/phase4g5r-final-wheel-env/bin/python -I /workspace/pyGameMath/benchmarks/slerp_wheel_smoke.py --output /tmp/phase4g5r-wheel-smoke.json
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 -m venv /tmp/phase4g5r-final-sdist-env
/tmp/phase4g5r-final-sdist-env/bin/python -m pip install --no-index --find-links /tmp/phase4d-dependencies six setuptools==84.0.0 wheel==0.48.0 packaging==26.3
/tmp/phase4g5r-final-sdist-env/bin/python -m pip install --no-index --no-deps --no-build-isolation /tmp/phase4g5r-final-dist-build/dist/gem-0.1.12.tar.gz
/tmp/phase4g5r-final-sdist-env/bin/python -m pip check
/tmp/phase4g5r-final-sdist-env/bin/python -I /workspace/pyGameMath/benchmarks/slerp_wheel_smoke.py --output /tmp/phase4g5r-sdist-smoke.json
```

Logs and JUnit records use `/tmp/phase4g5r-*`. Collect their final evidence:

```sh
$PY audit/phase4g5r_collect.py --output audit/phase4g5r-verification.json
git diff --check
```

Changed files: the SLERP implementation, removal of one marker covering the
12 original cases, the independent repair test module, quaternion API notes,
two benchmark/smoke scripts, this report, its evidence collector, verification
JSON and three measured performance JSON files. Packaging and version metadata,
all other runtime modules and all other existing tests are unchanged.
