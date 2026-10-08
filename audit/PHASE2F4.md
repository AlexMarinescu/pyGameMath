# Legendre recurrence and core promotion

E05 is corrected in `gem.legendre`. `gem.experimental.legendre` is a thin
reexport of the identical `Legendre` class. The existing spherical-harmonics
consumer imports the core class; no basis or normalization formula changes.

## Mathematics and domains

The supported associated domain is integer 0 <= m <= l with real x in [-1,1].
The source establishes unnormalized P_l^m with Condon–Shortley phase:
P_m^m(x) = (-1)^m (2m-1)!! (1-x²)^(m/2).
The next seed is P_(m+1)^m = x(2m+1)P_m^m. Each subsequent recurrence uses
the two preceding degrees, without recomputing either seed. Ordinary m=0
polynomials also evaluate outside [-1,1]. Negative associated orders are not
a supported Legendre domain; SPH handles negative harmonic orders through
positive-order Legendre functions and sine terms, as before.

Spherical harmonics retain their separate normalization K(l,m), real cosine
terms for m>0, sine terms for m<0 and polar/azimuth conventions. No implicit
normalization, sign change, clipping, epsilon or extreme-order policy is added.

Invalid inputs remain outside the accuracy contract. Existing arithmetic
continues to expose ValueError from sqrt for m>0 and |x|>1, and TypeError
from range for noninteger orders requiring iteration. An empty recurrence
with l<m retains its historical PML result. Negative degrees/orders and
nonfinite inputs receive no new standardized policy. Extreme degree/order
can overflow, underflow or lose accuracy; m=l=200, x=.2 still overflows to
infinity. NaN propagation is not replaced by a new validation contract.

## State and compatibility

`run()` now uses local recurrence state and leaves l, m, x and scratch fields
P, PM1, PML untouched. Repeated evaluation is deterministic even if scratch
fields contain previous results. Code relying on incidental scratch-field
updates from run must instead use the explicit helper methods.

Public constructor and helper signatures and return types are retained.
`mGreaterThan0` resets P to the diagonal seed; `calculatePM1` populates P
and PM1; `calculatePML(i)` computes the requested degree into PML. Helpers
retain their explicit mutable scratch-state behavior and None returns, but
repeated initialization no longer compounds prior state. Repository consumers
use only run; the historical wiki has no Legendre page or helper-state
contract. Existing fields remain available rather than being removed.

## Verification

Base: df8d96e207d5489ca6bf9cb2bc2eafa262a1aaf4, merged PR #26.
Phase 2F-3B baseline: **1475 passed, 13 xfailed**.
Running experimental defects with --runxfail reproduces all six E05 failures:
three recurrence cases, repeated nonzero-order evaluation and two higher-order
addition-theorem cases. P3(.2) historically returns -.6 rather than -.28.

Final full suite: **1642 passed, 0 failed, 7 xfailed** on Python 3.12.14,
pytest 9.1.1. Exact remaining identities equal the baseline minus six E05
cases: four unrelated defect cases and three unresolved contract questions.
See `phase2f4-test-results.json`.

Independent references differentiate explicit Rodrigues polynomial
coefficients with exact rational arithmetic. Tests cover degrees 0–12, all
valid orders, signed arguments and boundaries -1/0/1, recurrence/parity,
repeatability, scratch-state preservation and deterministic helpers.
Addition-theorem tests reach degree 12 across multiple orientations including
poles; known harmonic values retain the established sign convention.

An actual wheel build and pip installation into an isolated target directory
verify both import paths, class identity, polynomial values and degree-12 SH
interoperability. Existing packaging already includes gem and the transitional
experimental package, so no setup or dependency changes are needed.

## Performance

Five median 30,000-call trials compare constructor+run against the pre-fix
implementation at x=.37. New/old runtime ratios are about 1.23 for (l,m)=(0,0),
1.35 for (2,2), .45 for (12,0) and .44 for (12,5). Low-degree overhead rises
by roughly .09–.24 microseconds per call; higher degrees avoid redundant seed
work. The old high-degree answers are incorrect, so this compares runtime
rather than equivalent accuracy. These local measurements are not an
application performance guarantee. Evaluation is O(l+m) time and O(1)
auxiliary storage. Details are in `phase2f4-benchmarks.json`.

No spherical sampling, irradiance, transport or Bezier algorithm changes are
included. Further experimental promotion and final retirement remain separate.
