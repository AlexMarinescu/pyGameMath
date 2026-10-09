# Core performance baseline

Run from a clean checkout with Python 3.12 and the existing `six` dependency:

```
python benchmarks/core_baseline.py --output /tmp/gem-baseline.json --trials 7 --target-seconds .02 --profiles
```

Add `--historical` to compare Matrix3/Matrix4 inversion with Phase 2F-1 in the
same interpreter. This option requires the historical commit object
`2cbd899a47546b8c40f8244f57799969c86b0d87`; shallow checkouts/source archives
can run every current-core benchmark without it. No runtime package changes,
optional benchmark dependency, RNG seed or downloaded asset is required.

Inputs are constants in `prepare()`. Ordinary vectors range from -3.75 to 3.75;
matrices are nonsymmetric diagonally dominant, with entries
`3*identity + .13*(i+1)/(j+1)`, and a separate nonsymmetric multiplication operand.
Extreme norms use 1e-300/1e300, mixed norms span 600 decimal orders, and raw
inverse cases use uniform 1e-200/1e200 scales. Historical known answers cover
1e-300 through 1e300 using a well-conditioned triangular matrix. Quaternions
represent unit rotations around nontrivial axes; interpolation uses t=.37,
including a .001-degree near-angle case. No singular/unsupported operation is
timed as successful mathematics.

Bezier controls are 3D, bounded by eight units in X, two in Y and .5 in Z;
a single-segment distance-tolerance sweep covers .1, .01, 1e-4, 1e-8 and 1e-15.
The core stores this as squared tolerance. Eight segments use tolerance .001.
The sweep records emitted point counts and reaches the depth-16 cap.
SH directions use a deterministic midpoint-Z/golden-ratio-azimuth lattice;
projection uses 64/256/1024 directions and 1/3/5 bands, precomputing basis values
outside the timed calls. RGB radiance is finite and asymmetric, approximately
[1,5], [.5,1.5], [.25,.25]. Basis-only cases cover degrees 0, 2, 8 and 12.
L2 rotation/convolution/reconstruction uses fixed signed coefficients; negative
coefficients are valid and do not imply negative radiance input.

Each case executes once before calibration. Calibration doubles its operation
count until a trial reaches the requested wall time; seven trials then report
all samples, median, median absolute deviation and extrema. GC is disabled only
inside timeit and restored by it. Data preparation and import timing are separate.
Returned object/ctypes allocation is included in public-operation timing;
constructor and raw-kernel cases help locate overhead but are not a causal
subtraction experiment. Ray transform cases include duplication to avoid
accumulated mutable state, and report duplicate cost separately.

Fresh interpreter import timing excludes process launch and uses warm filesystem
caches; it is not a disk-cold machine startup measurement. cProfile and single-call
tracemalloc runs occur after timings and their instrumented durations are not
benchmark timings. Profile cumulative times overlap; never sum them. Peak traced
bytes measure Python allocation pressure, not RSS or total native allocation.
Machine load, CPU frequency and timer overhead are uncontrolled. A no-op baseline
is reported without subtracting it. Repeat runs on the same machine, compare
medians and variability, and validate numerical behavior before performance claims.

Ray intersection is unavailable in the supported implementation and deliberately
not invented for benchmarking. Core source and harness SHA-256 fingerprints
identify the measured code; JSON contains environment and operation counts.


If a full repeat shows large outliers, recheck selected cases across interleaved
rounds without a core change:

```
python benchmarks/recheck_variability.py --output /tmp/gem-recheck.json
```

This executes nine cases in three rounds, seven trials each, with .05-second
calibration. It isolates run-to-run variability from mathematical changes; it
does not automatically declare performance gains or impose a noise threshold.
