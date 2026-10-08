# Adaptive Bezier sampling

E04 is repaired in `gem.bezier`. The experimental module and former private
legacy module reexport the same core class, with no duplicate algorithms.

## Historical behavior and compatibility

The historical wiki has no Bezier page. The source establishes
`minimum_sqr_distance = 0.01`, squared-distance builder arguments, nested
lists from `getDrawingPoints`, and `None` returns from builders.
`interpolate(segmentPoints, scale)` explicitly appends controls. That
accumulation behavior is retained, including the possibility of malformed
layouts if complete independent paths are appended without arranging their
shared controls. Sampling rejects malformed nonempty layouts rather than
silently dropping trailing controls.

`samplePoints` historically invokes interpolation and therefore accumulates
stale generated controls. It now builds a fresh control list on each
successful call. This is a compatibility change for callers relying on that
accumulation; explicit accumulation remains available through `interpolate`.
Calls with fewer than two source points remain no-ops. Existing controls are
retained during interpolation, but the containing list is detached before
appending so a caller-owned list supplied by `setControlPoints` is not
extended. Generated Vector controls and sampled points have fresh storage.

## Sampling criterion

Midpoint de Casteljau subdivision uses maximum interior-control distance to
the finite endpoint chord. This is perpendicular distance when the projection
falls on the chord, and distance to the nearest endpoint otherwise. The
endpoint extension is necessary to detect collinear overshoot/backtracking.
Coincident endpoints use distance to that endpoint, avoiding division by zero.
The same criterion applies to scalar, Vector2 and Vector3 controls.

The public `minimum_sqr_distance` remains a positive finite **squared**
coordinate distance; effective distance tolerance is its square root.
Default 0.01 therefore means distance 0.1. Subdivision preserves increasing
parameter order, includes both standalone endpoints and emits each connected
segment boundary once. `getDrawingPoints()` retains its historical nested
per-curve list shape; later lists omit their first shared endpoint. Fully
coincident standalone curves retain two separately allocated endpoint samples.

The control-polygon criterion is conservative in exact arithmetic, but no
absolute floating-point approximation guarantee is made. Maximum depth is
16 per segment, bounding output to 65,536 chords / 65,537 samples. At the cap,
output is best effort and may exceed tolerance. Extreme arithmetic scales and
conditioning are not covered by a new numerical-robustness policy.
`segments_per_curve` and misspelled `divison_threshold` remain available for
compatibility; adaptive sampling no longer uses the historical angle heuristic
or those tuning fields.

`findDrawingPointsAdded(curveIndex,t0,t1,pointList,insertionIndex)` restricts
the cubic to the requested interval, inserts ordered fresh interior samples,
and returns the number inserted. The list's existing interval endpoints and
other elements remain intact. Invalid intervals raise ValueError; invalid
curve/insertion indices raise IndexError. Sampling requires complete 3*k+1
layouts and finite matching scalar/Vector2/Vector3 controls; malformed layouts,
representations and invalid tolerances raise ValueError. Empty paths return [].
Evaluation APIs retain their existing domains and signatures.

## Source-point construction

`samplePoints` retains ordered source endpoints. An interior vertex is kept
when its squared distance from the last retained vertex reaches
`minSqrDistance`, or skipping it would put the next vertex beyond
`maxSqrDistance`. Comparisons use squared distances directly. Thresholds must
be finite with 0 <= min <= max and max > 0. These are thinning heuristics,
not strict gap limits or curve-error bounds; sparse source gaps can exceed
max. These thresholds are independent of adaptive subdivision tolerance.

Interpolation retains historical endpoint/interior tangent formulas and
scale interpretation, using nonmutating local arithmetic. Coincident
neighbours give a zero tangent. Both builders return None and preserve source
objects/storage. Scale must be finite.

## Verification

Latest master baseline: 6916a52caab384b1623611c690ef34ff8cdf08fc (PR #25).
Baseline: **1441 passed, 15 xfailed**. Both E04 cases failed with `--runxfail`.
Diagnostic bypasses exposed scalar-left Vector multiplication and invalid
list.append operations. Source inspection additionally confirms missing
sqrMagnitude, `(t0-t1)/2`, math.abs, inappropriate one-element buffer lookbacks,
and no recursion bound. The old tail adjustment could move source geometry;
source thinning now retains original vertices instead.

Final suite: **1475 passed, 0 failed, 13 xfailed** (Python 3.12.14,
pytest 9.1.1). Exact remaining identities equal baseline minus the two E04
cases: ten unrelated defect cases and three unresolved contract questions.
See `phase2f3b-test-results.json`.

Independent polynomial references cover cubic and degree-elevated quadratic
curves in 2D/3D, ordered parameters, dense-reference chord error and tolerance
reduction. Tests cover scale sensitivity, coincidence, endpoint loops,
collinear overshoot, connected boundaries, subinterval insertion, ownership,
rebuilding, intentional append behavior, equal/zero thresholds and sparse gaps.
An actual depth-limited curve emits 65,537 samples and demonstrates error
above its requested tolerance. The two original E04 probes now pass through
core imports. A built wheel verifies core and compatibility class identity,
adaptive sampling and source construction in an isolated import directory.

No unrelated experimental algorithms are promoted or repaired. The dedicated
experimental retirement cleanup remains separate.
