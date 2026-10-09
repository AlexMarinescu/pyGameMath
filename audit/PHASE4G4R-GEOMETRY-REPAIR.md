# Phase 4G-4R — Geometry numerical repair

Equal-index refraction now preserves grazing directions, and Newell polygon
normals retain local geometry at large coordinate offsets. All 14 original audit
failures pass with unchanged inputs, assertions and tolerances.

Base: `495dc3d5b4d5f8bd4969b690b3a207a9f3339023` (merged PR #56).
Branch: `repair/phase4g4-geometry`.

## Repairs and compatibility

### 4G4-A01: grazing-angle refraction

For opposing unit inputs at equal refractive indices, the normal correction is
`d+sqrt(d*d)=d+abs(d)=0`, where `d=normal.dot(incident)`. The former implementation
recovered `d*d` as `1-(1-d*d)`. At d=-1e-9, the inner subtraction rounds to 1,
the discriminant becomes zero, and the correction incorrectly removes the normal
component. At d=-1e-100, squaring and subsequent subtraction lose it as well.

The equal-index path uses the exact zero correction when the calculated dot is
in the unit/opposing interval [-1,0]. All other paths retain their existing
arithmetic. Dot calculation and Vector multiplication/subtraction still execute,
preserving operand dispatch, dimension errors, fresh result storage and normal
orientation behavior. This is an exact identity, not an epsilon, near-critical
clamp, implicit normalization or normal flip. Ordinary Snell refraction and the
zero-Vector total-internal-reflection sentinel remain unchanged.

### 4G4-A02: translated polygon normals

Newell's sum combines coordinate sums and edge differences. At a world origin
of 2^52, sums around 2^53 lose unit offsets before cancellation. An exactly
represented, well-conditioned triangle can acquire a zero or incorrect normal.

The wrapped calculation now subtracts the first vertex before forming these
terms. Translating every point by the same anchor preserves Newell's oriented
area mathematically. A rolling tuple of local coordinates avoids allocating
temporary Vectors or a second polygon. The final exact-zero guard and existing
stable Vector normalization are unchanged. Empty, collinear, coincident and
zero-area self-cancelling polygons retain ZeroDivisionError.

The helper remains a Newell area normal for nonplanar polygons, not a least-squares
fit. Reversal changes orientation; cyclic starts and repeated closing vertices
are supported. `bestFitD` and its per-vertex mean weighting are unchanged.

Only `refract` and `Plane.bestFitNormal` change in runtime code. Public names,
signatures, dimensions, return types, ownership, conventions, compatibility
imports, dependencies and package/version metadata are preserved.

## Independent verification

The original reproductions ran before either implementation changed:
**14 failed, 227 deselected** (six A01 and eight A02). Baseline master had
**3,540 passed, 14 xfailed, zero failures and zero skips**.

After correction, all 14 cases passed with strict markers disabled before their
two defect decorators were removed. An AST comparison verifies that no audit
test body, parameter, assertion or tolerance changed. The new module adds
**160 passing cases**, using exact Fraction triangle-fan areas and 110-digit
Decimal normalization rather than another gem geometry operation as an oracle:

- Vector2/3/4 grazing components through 5e-324, both travel directions and
  normal orientations, non-axis-aligned interfaces, exact tangency and normal
  incidence, input and storage preservation.
- Independently calculated ordinary Snell transmission and TIR, retained
  unsupported scalar dispatch and historical wrong-oriented-normal behavior.
- Triangles, quads, concave and nonplanar polygons, cyclic starts, both winding
  directions, open/closed lists and positive/negative/mixed coordinate offsets.
- Uniform scales 2^-200 to 2^200; world offsets/edges from 1e-140/1e-150 to
  1e150/1e140, using the exact stored coordinates in the references.
- Empty, one/two-point, collinear, coincident and self-cancelling polygons.

Final full suite: **3,714 passed; 0 failed, 0 xfailed, 0 skipped, 0 errors**
in 19.13 seconds. The geometry/refraction/plane/numerical modules account for
580 passing cases in that run. The initial focused run passed 572 cases before
the final eight exponent-boundary cases were added.

Executed environment: CPython 3.12.14, Linux x86_64, pytest 9.1.1, six 1.17.0.
No other interpreter compatibility claim is made.

Machine-readable [verification](phase4g4r-verification.json) records every original
failure and its repaired status, artifact hashes and preservation checks. Both
clean installed artifacts pass all 14 original configurations and additional
core/import/degeneracy checks from `/tmp`, under isolated `python -I`. All 17
runtime modules in each wheel/sdist match the source bytes. The other 15 runtime
modules and packaging metadata match master byte-for-byte. AST checks also verify
that no adjacent vector or plane mathematics changed.

## Matched performance comparison

Two independent runs load immutable master source and the repair into the same
interpreter. Each workload has nine alternating before/after trials, 100 warmup
calls, process CPU timing and GC disabled only during timing. Setup, imports,
reference calculations and data generation are excluded; loop overhead and result
allocations are included. Refraction and short polygons have 20,000 calls/trial;
64-vertex polygons have 1,250 and 256-vertex polygons have 312.

Refraction inputs are unit Vector3 directions with normal +Y and ratios 2/3,
1, 1.25 and 1.5. The grazing component is -1e-100. Polygon datasets have 3/4/8/64/256
vertices with X in [-3,3], Y in [-2,2], Z in approximately [-1.12,1.12], generated
deterministically from an elliptical planar loop. Correctness references are
untimed and independent; the historically incorrect grazing result is never an
expected answer. These timings measure completed calculations, not a caught
baseline degeneracy exception versus a repaired normal.

CPU: Intel Xeon Platinum 8370C, shared host. Results, complete trials, median
absolute deviations, ranges, environment and source hashes are in
[run one](phase4g4r-performance.json) and
[run two](phase4g4r-performance-repeat.json). A latency ratio above 1 is slower.
Ratios are medians of paired ratios, not ratios of separate medians.

| Workload | Run 1 before → after (µs) | Paired ratio | Run 2 before → after (µs) | Paired ratio |
| --- | ---: | ---: | ---: | ---: |
| Equal indices, ordinary | 2.518 → 2.431 | 0.969 | 2.517 → 2.484 | 0.979 |
| Equal indices, grazing | 2.505 → 2.466 | 1.002 | 2.451 → 2.699 | 1.041 |
| Air to glass | 2.483 → 2.488 | 0.995 | 2.482 → 2.493 | 0.994 |
| Glass to air | 2.679 → 2.931 | 0.970 | 2.667 → 2.544 | 1.005 |
| TIR | 0.940 → 0.951 | 1.041 | 0.745 → 0.763 | 1.014 |
| Critical boundary | 2.616 → 2.668 | 1.020 | 2.553 → 2.585 | 1.002 |
| Polygon 3 | 4.533 → 4.860 | 1.078 | 4.532 → 5.095 | 1.113 |
| Polygon 4 | 4.938 → 5.266 | 1.078 | 4.972 → 5.863 | 1.134 |
| Polygon 8 | 6.193 → 6.932 | 1.141 | 6.488 → 7.073 | 1.122 |
| Polygon 64 | 22.784 → 27.513 | 1.229 | 21.973 → 28.980 | 1.319 |
| Polygon 256 | 76.802 → 112.777 | 1.449 | 77.649 → 105.476 | 1.340 |

Refraction ratios fall within the observed paired variation; there is no supported
speed claim. Local-coordinate polygon arithmetic adds measurable cost, especially
for larger polygons: approximately 23–45% in the 64/256-vertex workloads across
the two runs. Triangle paired medians increase 8–11%. Shared-host frequency and
scheduling variation are uncontrolled; these are measured trade-offs on this
environment, not universal performance bounds.

Untimed deterministic comparisons use seed 4404. Other-index refraction is
bit-identical in all 1,000 cases. Equal-index ordinary results are bit-identical
in 977/1,000, with maximum absolute change 1.11e-16. Random four-vertex Newell
results are bit-identical in 53/1,000, with maximum absolute change 9.53e-15.
Those arbitrary polygons can have area cancellation; this comparison supplements,
rather than replaces, the exact geometric references.

## Reproduction commands

From the repository root, using the development interpreter:

```sh
python -m pytest tests/test_geometry_audit.py --runxfail -q -k 'equal_media_grazing or translated_polygon'
python -m pytest tests/test_geometry_audit.py tests/test_geometry_repair.py tests/test_angles_refraction.py tests/test_planes.py tests/test_numerical_robustness.py -q
python -m pytest -q -o junit_family=xunit1 --junitxml=/tmp/phase4g4r-full.xml
python benchmarks/geometry_repair_compare.py --output audit/phase4g4r-performance.json
python benchmarks/geometry_repair_compare.py --output audit/phase4g4r-performance-repeat.json
```

The first command fails all 14 cases on the recorded base and passes after repair.
Before-edit baseline and reproduction JUnit reports were saved separately as
`/tmp/phase4g4r-baseline.xml` and `/tmp/phase4g4r-original.xml`. Exact executed
commands used `/workspace/.venvs/pyGameMath/bin/python` in place of `python`;
benchmark defaults were nine trials and 20,000 iterations.

Clean packaging used a staging copy in `/tmp/phase4g4r-dist-build`, containing
`gem/` without caches, the unchanged packaging files and manifest-listed docs:

```sh
cd /tmp/phase4g4r-dist-build
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 setup.py sdist bdist_wheel
```

Two newly created venvs installed six 1.17.0 and their respective artifacts via
`pip install --no-index --no-deps --no-build-isolation`. The sdist build environment
also contained setuptools 84.0.0, wheel 0.48.0 and packaging 26.3. Offline build
tools are verification-only, not added gem dependencies. Executed installed checks:

```sh
cd /tmp
/tmp/phase4g4r-wheel-env/bin/python -I /workspace/pyGameMath/audit/geometry_repair_smoke.py > /tmp/phase4g4r-wheel-smoke.json
/tmp/phase4g4r-sdist-env/bin/python -I /workspace/pyGameMath/audit/geometry_repair_smoke.py > /tmp/phase4g4r-sdist-smoke.json
cd /workspace/pyGameMath
/workspace/.venvs/pyGameMath/bin/python audit/verify_geometry_repair.py --output audit/phase4g4r-verification.json
git diff --check
```

The verification collector checks actual JUnit outcomes, unchanged audit AST,
source preservation and built artifact contents; it does not substitute for
running tests or building/installing the artifacts.

## Remaining numerical limitations

Neither repair recovers geometric differences rounded away before input reaches
gem. Newell area products can still overflow/underflow or cancel for extreme or
nearly zero-area geometry; this is not a universal polygon error bound. The
nonplanar normal retains the historical oriented-area meaning.

Refraction retains the existing discriminant for unequal indices, including its
critical-angle rounding behavior and extreme-IOR limitations. Unit, matching
vectors and an opposing normal remain prerequisites. No broader malformed or
nonfinite-input policy is defined. Dot/cross/barycentric range issues, mean-offset
overflow, ray endpoint semantics and other findings classified separately in the
audit remain outside these two repairs.

## Changed files

- `gem/vector.py`, `gem/plane.py`: the two numerical repairs.
- `tests/test_geometry_audit.py`: remove only the two corrected defect markers.
- `tests/test_geometry_repair.py`: independent boundary/ownership regressions.
- `benchmarks/geometry_repair_compare.py`: immutable-base paired measurements.
- `audit/geometry_repair_smoke.py`, `audit/verify_geometry_repair.py`: isolated
  artifact checks and evidence/preservation collection.
- `audit/phase4g4r-performance.json`, `audit/phase4g4r-performance-repeat.json`,
  `audit/phase4g4r-verification.json`: executed machine-readable results.
- `audit/CONVENTIONS.md`, `audit/COMPATIBILITY.md`, this report: numerical and
  compatibility documentation.
