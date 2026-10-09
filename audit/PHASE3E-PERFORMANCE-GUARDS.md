# Phase 3E — Performance regression safeguards

Base: `45b708ec44da8be8352e6ce3b70abfb4c2737773` (after merged PR #37).
Timing captures use `553f1664987ab45b4212a4e5c105ebbf92493e81`; the only intervening
master change updates LICENSE, with identical mathematical source.
Mathematical implementations, packaging metadata, dependencies and public APIs
are unchanged. Benchmark tools remain outside the installed runtime package.
Python 2.7 compatibility has not been evaluated.

## Workflow

`benchmarks/run_core.py` reuses Phase 3A's deterministic datasets and timing
routine. Its 103-workload quick suite covers all supported core categories;
five extended cases retain costly adaptive subdivision and 1024-sample SH
projection. Raw kernels and wrappers stay distinct. Input construction, calibration
and optional profiling are excluded from warm trial timings. Returned allocation
and ctypes synchronization remain included. Ray intersection is unavailable in
the supported API and is explicitly identified, not fabricated.

Versioned workload hashes describe inputs, driver code and units rather than the
mathematical implementation. JSON retains every trial, calibrated operation count,
round, environment, source fingerprint and adaptive output count. See
[the workflow and CI policy](../benchmarks/REGRESSION_POLICY.md) for clean-checkout
commands, dataset ranges, baseline versioning and refresh rules. Original Phase
3A–3D scripts and measurements remain available without modification.

`benchmarks/compare.py` reads repeated versioned captures or opposite sides of
one historical paired artifact. It refuses performance classifications for unknown
or mismatched environments, changed workloads or insufficient repetitions.
Reports include absolute timings, medians/MAD, per-block variability, reciprocal
ratios and seeded bootstrap intervals. Repeated copies of one artifact cannot
inflate evidence. Historical files do not gain invented host/workload provenance.

Review heuristics default to a 15% relative margin and .05 µs absolute floor,
increased by observed within/between-block noise. Classification also requires
three blocks, five trials per block, an interval beyond the margin and 80% block
agreement. These are advisory flags, not pass/fail timing assertions or accuracy
contracts. All valid comparisons exit successfully. Inconclusive means insufficient
evidence for a flag, not proof of equivalent performance.

## Measurements

CPython 3.12.14, Linux x86_64/glibc 2.41, Intel Xeon Platinum 8573C, six 1.17.0.
Full host identity hash, affinity and CPU count are recorded in each new artifact.
The two quick captures use the same core source, three rounds of seven trials,
and .005-second calibration. All 103 comparisons are inconclusive; there are no
regression or improvement flags. These unchanged-code runs demonstrate ordinary
environmental variation, not an optimization result.

| Workload | Capture A µs | Capture B µs | B/A |
|---|---:|---:|---:|
| Vector3 dot | .286 | .270 | .944 |
| Matrix3 raw inverse | 11.859 | 11.550 | .974 |
| Matrix4 raw inverse | 20.176 | 21.834 | 1.082 |
| squad4 | 7.174 | 7.381 | 1.029 |
| Cubic Bezier | 1.044 | 1.084 | 1.038 |
| L2 RGB SH rotation | 27.654 | 23.563 | .852 |
| Plane normalization | 2.628 | 2.668 | 1.015 |
| Ray translation | 6.183 | 5.852 | .946 |

Extended medians: 16,395-point subdivision 287.7 ms; depth-16, 65,537-point
subdivision 1152.0 ms; 1024-sample SH projection 0.953/3.666/8.699 ms for
1/3/5 bands. Depth and work counts are preserved. Separate profiles cover raw
Matrix4 inversion, squad4 and RGB SH rotation, with traced single-call peaks of
3,433/616/3,872 bytes respectively. These are instrumented allocation observations,
not RSS or benchmark latency.

The Phase 3D known-data adapter reports nine possible improvements and 24
inconclusive comparisons from the existing paired artifact. It does not compare
that historical host with the new captures. Synthetic reporting fixtures correctly
identify a 2× slowdown, a 2× improvement, a 2% inconclusive difference, a host
mismatch and missing/changed workloads. They are labeled synthetic, not measured
performance. Raw captures, profiles and reports are `audit/phase3e-*.json`.

## Verification

Commands executed with `/workspace/.venvs/pyGameMath/bin/python`:

```sh
python -m pytest -q --junitxml=/tmp/phase3e-tests.xml -o junit_family=legacy
python benchmarks/run_core.py --output audit/phase3e-quick.json --rounds 3 --trials 7 --target-seconds .005
python benchmarks/run_core.py --output audit/phase3e-repeat.json --rounds 3 --trials 7 --target-seconds .005
python benchmarks/run_core.py --suite extended --output audit/phase3e-extended.json --rounds 3 --trials 7 --target-seconds .005
python benchmarks/run_core.py --case matrix4_raw_inverse --case quaternion_squad4 --case sh_rotate_l2_rgb --profiles --output audit/phase3e-profiles.json --rounds 3 --trials 7 --target-seconds .005
python benchmarks/compare.py --baseline audit/phase3e-quick.json --current audit/phase3e-repeat.json --output audit/phase3e-comparison.json
python benchmarks/compare.py --baseline audit/phase3d-benchmarks.json --baseline-side before --current audit/phase3d-benchmarks.json --current-side after --output audit/phase3e-historical-comparison.json
python benchmarks/verify_hdr_reference.py --output audit/phase3e-hdr-reference.json --output-dir /tmp/phase3e-hdr-reference
```

Full regression suite: **2,249 passed, zero expected failures, zero unexpected
failures, zero skips**. The 55 new reporting tests use independent synthetic
timings and executable CLI checks; they contain no wall-clock performance assertions.
Existing numerical stability, ownership, ctypes, Bezier termination, SH rotation,
linear HDR reference and installed-wheel regressions remain intact.

A fresh local checkout of the base plus the staged patch ran all 103 quick and
five extended workloads with one round, three trials and .001-second calibration.
This checks executable reproducibility, not calibrated performance equivalence.
It also ran reporting tests and the clean wheel/sdist isolated import test:
**56 passed**. The wheel is built from a clean source copy, installed offline into
an isolated target, and tested outside the source tree, including compatibility
imports and representative numerical operations. No packaging changes are needed.

HDR regeneration produced seven byte-identical files on this host. Coefficients,
selected linear values and decoded pixels also matched. Independent numerical
tests remain the portable correctness oracle; matching PNG hashes alone do not
establish mathematical correctness across platforms. Verification counts and
unchanged runtime/packaging fingerprints are recorded separately in JSON.

## CI recommendation

Keep ordinary PR checks functional. Optional smoke jobs may execute a few cases
and upload JSON without timing gates. Calibrated comparisons belong on a pinned,
dedicated runner through manual or scheduled jobs, with multiple captures and
reviewed baseline refreshes. Matching machine metadata cannot rule out frequency,
load or virtualization changes. Review flagged differences with profiles and
independent correctness evidence before accepting or rejecting a change.
