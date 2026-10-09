# Bezier and spherical-harmonics performance

Base: master `09fb57a31fd1b62fad027822cb4ab1eba52ca158`, merged PR #36.
Branch: `perf/phase3d-bezier-sh`.

## Focused changes

Native, equally sized Vector controls use the same quadratic/cubic Bernstein
weights and left-associated component additions in one fresh result Vector.
The previous evaluator constructed five/quadratic or seven/cubic wrappers and
component lists. Scalar, mixed-representation/dimension, subclass and other
parameter cases retain the original operator expressions. The native shortcut
does not add accepted representations, dimensions or parameter policies.

Adaptive subdivision shares the endpoint chord, its hypot length and its unit
direction across each polygon's interior controls. Offsets, ordered sum, clamped
projection and residual hypot are unchanged. Coincident endpoints still use
distance to that endpoint. This is the same finite-segment flatness criterion,
including overshoot/backtracking; squared tolerance and maximum depth 16 remain.
The iterative stack, midpoint de Casteljau split, endpoint/join ordering, path
validation and source builders are untouched. No points are dropped for speed.

SH basis arrays cache immutable direction-independent normalization layouts,
with a 16-entry LRU limit. Each layout preserves the original root2*K grouping.
One cos(theta) and one deterministic associated Legendre evaluation per
(degree,absolute order) replace repeated calls for paired +/- orders. Results
remain fresh arrays. SPH, K, Factorial, Legendre and GenerateSamples stay
unchanged; no polynomial recurrence or historical basis is redesigned.

Analytical rotation supplies the same three products in the same order to
builtin sum using a fixed tuple instead of a generator. It deliberately retains
sum's interpreter-specific accumulation behavior, including CPython 3.12's
improved floating summation; explicit naive additions would not preserve it.
All rotation matrices, validation and coefficient extraction remain unchanged.
Radiance projection retains math.fsum; angular projection retains its existing
compensated update and O(coefficient count) accumulation storage.

## Correctness and compatibility

94 new independent regressions cover rational de Casteljau evaluation including
empty/generic dimensions and extrapolation; subclass and mixed-dimension
fallbacks; finite-segment/coincident flatness known answers; dense polynomial
approximation, ordered sampling, fresh storage and shared joins; and a depth-16
parabola with exactly 65,537 samples and independent parameter references.

Cartesian SH polynomials verify basis/reconstruction, known weighted constant
and asymmetric environments, cosine factors, cancellation/channel isolation,
active rotations, noncommuting composition and band energies. Cache tests verify
independent repeated arrays, bounded layout count and agreement at higher orders.
Existing quaternion boundary, invalid-input, extreme-scale, integration, shader
and ownership tests retain all original bounds. No expected failure is removed.

The preservation tool checks **1,313 bit-identical result groups** against
independently loaded baseline modules: 600 Vector evaluations, 150 flatness
values, 150 full sampled curves, 180 bases, 200 scalar/RGB rotations and all 33
benchmark workloads. Full depth-limited sample coordinates are included. All
20 unaffected top-level definitions, including BezierPath and SH classes, are
AST-identical. No deterministic rounding difference is intended or observed.
Returning results own independent storage; builders and compatibility imports
keep established semantics. Coefficient signs/order, weights, normalization,
unit-direction prerequisites and the legacy conversion boundary are unchanged.

## Reproduction

```sh
python benchmarks/bezier_sh.py --output audit/phase3d-benchmarks.json
python benchmarks/bezier_sh.py --output audit/phase3d-repeat.json
python benchmarks/summarize_bezier_sh.py audit/phase3d-performance-summary.json audit/phase3d-benchmarks.json audit/phase3d-repeat.json
python benchmarks/verify_bezier_sh.py audit/phase3d-equivalence.json
python benchmarks/measure_sh_cache.py audit/phase3d-cache-memory.json
python benchmarks/verify_hdr_reference.py --output audit/phase3d-hdr-reference.json --output-dir /tmp/phase3d-hdr-reference
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase3d-tests.xml
```

These commands ran with `/workspace/.venvs/pyGameMath/bin/python`. Baseline
comparison needs the merged PR #36 git object. Tools remain outside gem;
functools and all example tools use the standard library. No runtime dependency
or packaging/API change is introduced.

## Methodology and input ranges

Two consecutive runs each use three alternating before/after rounds and seven
trials per case, with a calibrated minimum .02-second timed batch. Warm calls
and count-doubling calibration precede samples. Timeit disables GC only during
measurement. Slow depth-limited single calls exceed the target. No tests or other
benchmark workloads run concurrently. Input generation is excluded; fresh result
allocation and public validation are included. JSON retains individual counts,
samples, ranges, MADs, source hashes and separate profiles.

Phase 3A data supply 33 cases: 3D quadratic/cubic controls [0,0,0], [.2,2,.5],
[.8,-2,-.5], [1,0,0]; scalar cubic [0,2,-2,1]; parameter .37; adaptive distance
tolerances .1,.01,1e-4,1e-8,1e-15 and an eight-segment path at .001. The public
field remains tolerance squared. Raw helpers include splitting, flatness and
subdivision. No conversion or setup time is hidden in returned sample counts.

SH projection uses 64/256/1024 deterministic sphere samples, bands 1/3/5,
RGB [1+4*max(X,0)^8,.5+max(Z,0),.25] and uniform 4*pi/N weights. These APIs
project RGB; scalar coefficient rotation is separately supported, not a new
scalar projection API. SPH controls include degrees 0,2,8,12. Full basis arrays
use 3/5/9 bands at theta=.73, phi=1.27. Rotation uses the existing 73-degree
axis [1,2,3], with scalar/RGB L2 rows and a +Z reconstruction normal. Added
angular projection uses a 32x16 deterministic RGB gradient [2+(col+.5)/32,
1+(row+.5)/16,.25]. Its basis is computed inside the timed call.

Warm basis-layout timings exclude cache initialization. Cache-miss and retained
storage are measured separately below. cProfile/tracemalloc runs are separate;
their instrumented times are not latency and cumulative frames overlap.
Traced Python peaks are not RSS or a count of all native allocations.

Environment: CPython 3.12.14, Linux-6.18.44-x86_64-with-glibc2.41,
CPU INTEL(R) XEON(R) PLATINUM 8573C, six 1.17.0; pytest 9.1.1.
No affinity/frequency or shared-host load isolation is imposed. Descriptive
intervals use 10,000 seeded paired-round resamples of six blocks; independence
and stationarity are not guaranteed. They are conditional observations, not
universal speedups or population confidence guarantees. No cross-host baseline
numbers are treated as controlled comparisons.

## Paired latency results

Microseconds per named operation. Times are medians of round medians; paired
speedup is the median of before/after block ratios, not necessarily the quotient
of the displayed aggregate times. MAD is median within-block relative MAD.

| Operation | Before us | After us | Paired speedup | Descriptive 95% interval | MAD before/after % |
| --- | ---: | ---: | ---: | --- | ---: |
| `bezier_quadratic` | 3.015 | 0.797 | 3.745x | 3.590-4.204 | 2.4/2.5 |
| `bezier_cubic` | 4.332 | 1.001 | 4.245x | 4.066-4.640 | 4.2/1.7 |
| `bezier_cubic_scalar` | 0.252 | 0.280 | 0.902x | 0.891-0.918 | 1.9/1.9 |
| `bezier_subdivide_0.1` | 102.013 | 88.848 | 1.151x | 1.061-1.168 | 2.3/3.4 |
| `bezier_subdivide_0.01` | 320.446 | 285.536 | 1.130x | 1.093-1.168 | 1.5/2.5 |
| `bezier_subdivide_0.0001` | 2935.561 | 2553.518 | 1.150x | 1.123-1.165 | 2.8/3.1 |
| `bezier_subdivide_1e-08` | 309718.921 | 268531.013 | 1.143x | 1.065-1.166 | 1.7/1.7 |
| `bezier_subdivide_1e-15` | 1265676.082 | 1133282.381 | 1.142x | 0.911-1.160 | 3.0/1.7 |
| `bezier_path_8_segments` | 7453.205 | 6558.480 | 1.152x | 1.113-1.265 | 2.8/4.2 |
| `sh_basis_0_0` | 0.707 | 0.661 | 1.008x | 0.982-1.098 | 3.3/1.6 |
| `sh_basis_2_1` | 1.118 | 1.179 | 0.997x | 0.950-1.013 | 3.3/1.9 |
| `sh_basis_8_3` | 2.210 | 2.263 | 0.988x | 0.948-1.013 | 2.9/3.9 |
| `sh_basis_12_6` | 3.031 | 2.965 | 1.011x | 0.968-1.100 | 3.2/2.9 |
| `sh_project_n64_b1` | 59.917 | 60.372 | 0.994x | 0.967-1.037 | 1.9/1.9 |
| `sh_project_n64_b3` | 223.988 | 221.689 | 1.013x | 0.946-1.053 | 3.4/2.7 |
| `sh_project_n64_b5` | 576.374 | 567.470 | 1.019x | 0.851-1.046 | 2.4/1.9 |
| `sh_project_n256_b1` | 226.026 | 231.273 | 0.993x | 0.812-1.016 | 2.2/2.8 |
| `sh_project_n256_b3` | 872.088 | 835.959 | 1.035x | 0.995-1.080 | 2.1/2.0 |
| `sh_project_n256_b5` | 2313.416 | 2154.250 | 1.009x | 0.984-1.403 | 3.2/1.7 |
| `sh_project_n1024_b1` | 898.272 | 895.695 | 1.010x | 0.913-1.024 | 2.5/1.8 |
| `sh_project_n1024_b3` | 3507.099 | 3408.834 | 1.053x | 0.991-1.125 | 4.5/4.4 |
| `sh_project_n1024_b5` | 8770.120 | 8499.514 | 1.035x | 1.002-1.089 | 2.2/3.2 |
| `sh_convolve_l2` | 6.241 | 6.386 | 0.968x | 0.948-1.001 | 2.2/2.9 |
| `sh_rotate_l2_scalar` | 18.070 | 11.771 | 1.535x | 1.518-2.225 | 2.0/3.1 |
| `sh_rotate_l2_rgb` | 42.731 | 24.977 | 1.746x | 1.629-1.814 | 2.0/3.8 |
| `sh_reconstruct_l2` | 17.479 | 13.181 | 1.358x | 1.311-1.394 | 2.5/1.9 |
| `bezier_raw_split` | 5.223 | 5.169 | 1.023x | 0.988-1.036 | 2.4/2.9 |
| `bezier_raw_subdivide` | 307.423 | 265.202 | 1.150x | 0.866-1.240 | 3.1/2.4 |
| `bezier_raw_flatness` | 6.129 | 5.719 | 1.197x | 1.016-1.349 | 2.9/6.7 |
| `sh_raw_basis_b3` | 9.673 | 5.582 | 1.725x | 1.700-1.950 | 2.9/4.4 |
| `sh_raw_basis_b5` | 33.116 | 17.146 | 1.921x | 1.777-2.016 | 1.7/3.6 |
| `sh_raw_basis_b9` | 146.365 | 64.690 | 2.287x | 2.171-2.464 | 2.4/3.6 |
| `sh_angular_probe_b3` | 6163.675 | 4370.977 | 1.415x | 1.350-1.457 | 2.0/2.5 |

Measured gains exceed descriptive intervals for native Vector polynomial
evaluation (3.745x quadratic, 4.245x cubic), ordinary adaptive sampling
(about 1.13-1.15x), full SH bases (1.725-2.287x), scalar/RGB analytical rotation
(1.535x/1.746x), reconstruction (1.358x) and angular projection (1.415x).
The depth-cap and raw subdivision intervals include 1: no reliable improvement
is claimed for those individual noisy cases. Public SPH, weighted projection,
splitting and convolution are unchanged controls; small observed timing shifts
do not establish implementation gains/regressions.

Scalar cubic evaluation is a measured trade-off: .252 to .280 us, roughly
28 ns (11%) additional latency from the native-control dispatch check. Its
paired interval excludes 1. The check is necessary to select the allocation
shortcut while keeping scalar, mixed and subclass behavior. This is an
explained narrow slowdown, not measurement noise or a claim that every path
improves. Mathematical output is bit-identical. No other material slowdown is
demonstrated. The report retains this cost rather than hiding the scalar case.

## Adaptive work and memory

| Case | Baseline count | Optimized count |
| --- | ---: | ---: |
| `bezier_subdivide_0.1` | 7 | 7 |
| `bezier_subdivide_0.01` | 19 | 19 |
| `bezier_subdivide_0.0001` | 161 | 161 |
| `bezier_subdivide_1e-08` | 16395 | 16395 |
| `bezier_subdivide_1e-15` | 65537 | 65537 |
| `bezier_path_8_segments` | [51, 50, 50, 50, 50, 50, 50, 50] | [51, 50, 50, 50, 50, 50, 50, 50] |
| `bezier_raw_subdivide` | 19 | 19 |

Output grows sharply near the depth limit: 161 at 1e-4, 16,395 at 1e-8 and
65,537 at 1e-15 for the fixed curve. The geometric criterion/depth cap is
unchanged; tolerance may remain unsatisfied at depth 16. Timing gains come
from less repeated chord work, not fewer splits or weaker approximation.
The original de Casteljau and whole-path validation overhead remain visible.

| Workload | Traced peak before bytes | Traced peak after bytes |
| --- | ---: | ---: |
| `bezier_cubic` | 848 | 112 |
| `bezier_subdivide_0.0001` | 55416 | 55224 |
| `bezier_path_8_segments` | 107760 | 107568 |
| `sh_reconstruct_l2` | 1312 | 1400 |
| `sh_rotate_l2_rgb` | 4512 | 3344 |
| `sh_angular_probe_b3` | 4584 | 4880 |
| `sh_project_n256_b3` | 28000 | 28504 |

Cubic wrapper construction falls from seven to one per evaluated point;
quadratic falls from five to one. Adaptive output allocation still dominates
peak memory, so the chord reuse mostly reduces arithmetic/call overhead.
A basis array evaluates b*(b+1)/2 Legendre objects instead of b*b; L2 is six
instead of nine, with one theta cosine instead of nine. Warm normalization
layout hits avoid repeated factorial/root calculations. RGB rotation retains
63 three-term sums per call but removes their generator/range/zip iteration.

The basis dictionary/immutable cache trades some memory for arithmetic savings:
reconstruction and angular projection traced peaks increase in the primary run.
The unchanged weighted-projection peak also varies with allocator state; these
single peak observations are not allocation guarantees. Dedicated cache tracing
starts empty, discards results and records initialization/resident cost:
- Cold L2 layout plus basis peak: 1464 bytes; retained layout trace: 904 bytes.
- Layouts 1-16 retained trace: 160648 bytes; fill peak: 174536 bytes.

There are at most 16 layouts, each O(bands^2); the entry count does not cap
band order or guarantee a fixed byte bound independent of order. No direction,
sample, radiance, coefficient or caller-owned object is cached. Warm primary
timings do not imply zero initialization cost or no additional resident memory.

## HDR showcase verification

The verifier executes the actual core pipeline with:

```sh
python -m examples.hdr_sh.regenerate --output-dir /tmp/phase3d-hdr-reference
```

All seven generated files match committed bytes. Original/rotated radiance
coefficients, convolved irradiance rows, selected pre-tone-map irradiance and
Lambertian values, manifest statistics and decoded RGB pixels are identical.
The synthetic environment, camera, albedo, exposure, color conversion, shader
formulas and golden images are untouched. Independent Cartesian/pixel/shader
tests remain passing; file equality is supplemental evidence, not the sole
mathematical reference.

| Image | SHA-256 (committed and regenerated) |
| --- | --- |
| `sh_original.png` | `dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a` |
| `sh_rotated.png` | `c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7` |

Recorded zlib build/runtime: 1.3.2/1.3.2.
Byte identity is verified on this host; Python/libm, zlib and platform changes
can affect floating or compressed bytes. The JSON stores both file and decoded
RGB hashes and numerical pixel properties. No golden output is replaced.

## Full verification and scope

Full pytest: **2,194 passed, 0 expected failures, 0 unexpected failures, 0 skips**.
Baseline: 2,100 passes; 94 new independent cases account for the increase. Clean
offline wheel/sdist builds and isolated installed-wheel checks outside the
source tree pass, including polynomial evaluation, adaptive order, SH
reconstruction/rotation and fresh output storage. A fresh baseline clone plus
the staged patch executes all 33 before/after workloads, reproduces preservation
JSON and regenerates the same seven HDR reference files independently. No tests are deleted, skipped
or weakened. Results are saved in `phase3d-test-results.json`.

Production changes are confined to gem/bezier.py and gem/spherical_harmonics.py.
Benchmark/proof tools, tests, raw measurements, cache/HDR verification and audit
notes support the change. This phase does not include Phase 3E, documentation
modernization, merging or package publication.
