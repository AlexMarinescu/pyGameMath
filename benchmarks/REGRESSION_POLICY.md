# Performance monitoring

Benchmark tools run on Python 3, separately from the installed `gem` runtime.
They use the standard library and the existing `six` dependency. Runtime Python
version support is unchanged; these measurements do not verify Python 2.7.

From a clean checkout:

```sh
python3 -m venv /tmp/gem-perf
/tmp/gem-perf/bin/python -m pip install six pytest
/tmp/gem-perf/bin/python -m pytest -q
/tmp/gem-perf/bin/python benchmarks/run_core.py --list
/tmp/gem-perf/bin/python benchmarks/run_core.py --output /tmp/quick-a.json
/tmp/gem-perf/bin/python benchmarks/run_core.py --output /tmp/quick-b.json
/tmp/gem-perf/bin/python benchmarks/compare.py --baseline /tmp/quick-a.json --current /tmp/quick-b.json --output /tmp/comparison.json --markdown /tmp/comparison.md
```

Run sequentially on an otherwise idle host. Default quick runs have three rounds
of seven trials, calibrated to at least .01 seconds per trial. Alternating rounds
reverse workload order. Setup, imports, calibration and optional instrumented
profiles are outside the reported trial samples. All raw samples, operation counts,
block medians and variability are retained. Returned allocation and ctypes export
costs are included; raw kernels and public wrappers remain separate cases. The
original Phase 3A harness remains available for cold-import and historical inverse
measurements. Its dataset ranges and allocation caveats are documented in README.md.

The quick suite covers Vector2/3/4, Matrix2/3/4 including stable inverses,
Quaternion rotation/interpolation/conversion, Bezier evaluation and ordinary paths,
SH basis/projection/convolution/rotation/reconstruction, Plane and Ray operations.
There is no supported Ray intersection method to benchmark. It is explicitly
reported as unavailable rather than substituting an invented implementation.

## Stress and profiling

```sh
python benchmarks/run_core.py --suite extended --list
python benchmarks/run_core.py --suite extended --output /tmp/extended.json
python benchmarks/run_core.py --case matrix4_raw_inverse --case quaternion_squad4 --case sh_rotate_l2_rgb --profiles --output /tmp/profile.json
python benchmarks/verify_hdr_reference.py --output /tmp/hdr-verification.json --output-dir /tmp/gem-hdr-reference
```

Extended cases retain Bezier's full depth-16 workload and 1024-sample SH projection
with 1/3/5 bands. A slow single call can exceed calibration's target duration.
`--case` selects any named workload independently of the suite filter.
`--profiles` adds separate cProfile and tracemalloc observations; those durations
are not timing samples. High cumulative profile times overlap, and traced Python
allocation is not RSS. Existing Phase 3B–3D paired scripts provide targeted before/
after diagnostics and independently loaded historical implementations.

## Comparability and review

Each new result records the interpreter, OS, architecture, CPU, CPU count,
affinity, six version and hashed host identity. Unknown or different identity,
environment or timing methodology prevents a performance classification. Ratios
may still be displayed descriptively, with `environment_not_comparable` status.
Hardware metadata cannot detect changing CPU frequency, shared-host contention,
virtualization scheduling or all power policies. Matching metadata is necessary,
not sufficient evidence of controlled conditions.

Workload fingerprints cover deterministic input/driver source, metadata, units
and `core-v1` definition version. They deliberately exclude `gem` implementations
so an optimization can be compared. Changed fingerprints, units or observed Bezier
point counts are reported as changed work, never a speedup. Missing workloads are
explicit. Calibration may choose different iteration counts; comparisons use
per-call timings and retain those counts. Repeated files can be supplied with
repeated `--baseline` and `--current` options; duplicate artifacts are rejected.

At least three blocks per side and five trials per block are required to classify
a change. The tool reports medians, MAD, extrema, block medians, within-block
relative MAD, absolute changes, reciprocal ratios and a reproducible seeded
bootstrap interval. Unpaired captures resample blocks independently; only matching
blocks inside one historical paired artifact are resampled together. The interval
is descriptive conditional evidence, not a population confidence guarantee.

Default review flags require a ratio interval beyond a 15% relative margin,
an absolute difference exceeding .05 microseconds, and agreement from at least
80% of current blocks. The margin also grows with observed within/between-block
MAD. Both margins are configurable reporting heuristics, not mathematical accuracy
limits. Small, noisy or contradictory changes are `inconclusive`, which does not
establish equivalence. All timing reports exit successfully: there is no automatic
performance gate. Invalid files/options produce a normal command error.

Historical Phase 3B–3D paired JSON can be inspected without inventing missing
host/workload provenance:

```sh
python benchmarks/compare.py --baseline audit/phase3d-benchmarks.json --baseline-side before --current audit/phase3d-benchmarks.json --current-side after --output /tmp/historical-comparison.json
```

Only the identical paired artifact is comparable through this adapter. Different
historical files or new versioned captures cannot establish identical workloads
automatically. Use the original phase-specific summarizers for their documented
multi-run matched comparisons. Phase 3A's single-run unversioned JSON remains
exploratory historical data, not a calibrated regression baseline.

## CI and baseline refresh

Keep ordinary PR checks functional: the complete mathematical suite, reporting
contract tests and installed-wheel checks contain no timing assertions. An optional
performance smoke job can execute a few workloads with one round and three trials
to check the harness, upload JSON and avoid timing gates. It intentionally lacks
enough evidence for a regression classification on a noisy shared runner.

Run long calibrated measurements manually or on a scheduled dedicated runner with
pinned interpreter, dependencies, CPU allocation and documented machine settings.
Collect multiple before/after executions, alternate their order, preserve raw
artifacts and investigate flags with targeted profiles and independent correctness
tests. Compare historical implementations in the same process using existing
paired tools where feasible. Do not compare unrelated hosts as controlled data.
Manual review decides whether a flagged change is material.

Baselines belong to a specific workload version, host and environment. Refresh
after a reviewed change or intentional environment migration, record the commit,
reason and correctness evidence, and retain earlier raw captures. Never refresh
merely to silence a slowdown. Changes to input values, units or driver semantics
require a new workload definition version and a fresh measured baseline; checksum
changes already prevent accidental comparisons. Schema changes require a versioned
reader/migration, not reinterpretation of old records. Performance improvements
must remain subordinate to numerical, ownership and compatibility regressions.

The installed-wheel test builds with the existing setuptools/wheel configuration
using the base interpreter, installs without dependencies into an isolated target,
and exercises core and compatibility imports outside the source tree. Build tools
must be available to that interpreter. Packaging metadata is unchanged.
