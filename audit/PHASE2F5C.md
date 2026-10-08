# Analytical spherical-harmonics rotation

`rotate_coefficients(coefficients, orientation)` actively rotates canonical
scalar or RGB coefficients with 1, 4 or 9 entries. L0 is unchanged; L1/L2
rotate independently. Radiance and already convolved irradiance use the same
transformation. Historical probe coefficients require explicit conversion.

## Mathematical construction

For a column direction d, f_rotated(d)=f_original(R^-1 d). R is the active
right-handed rotation associated with the Hamilton Quaternion [w,x,y,z].
In gem's row-vector notation the inverse direction is d*R_row^T. Positive
90-degree Z rotation moves a +X feature to +Y. Quaternion composition acts
right operand first; rotating by a then b is represented by b*a.

Canonical L1 is the linear form a*(-c3*X-c1*Y+c2*Z), where
a=sqrt(3/(4*pi)). Its vector rotates by R. L2 is d^T T d with symmetric
traceless T: Txx=b*c8-a2*c6, Tyy=-b*c8-a2*c6, Tzz=2*a2*c6,
Txy=b*c4, Txz=-b*c7, Tyz=-b*c5, where a2=sqrt(5/(16*pi)) and
b=sqrt(15/(16*pi)). Active rotation is T'=R*T*R^T, followed by the inverse
basis mapping. Symmetric entries are averaged during extraction; c6 uses
the traceless diagonal combination to suppress roundoff trace contamination.
This is analytical algebra, with no sampled directions, integration or
projection in the implementation. Per-band orthonormal energy is preserved
within ordinary floating-point precision.

## Domain and ownership

Only Quaternion orientations are accepted. Four finite components must have
abs(hypot(w,x,y,z)-1) <= 1e-12. A temporary component copy is normalized;
the supplied Quaternion and its storage are untouched. Larger deviations,
zero/nonfinite values, malformed coefficient arrays and unsupported orientation
types raise ValueError. The tolerance is specific to this additive API and
does not modify other quaternion operations. q and -q represent the same
rotation. RGB/scalar arrays cannot be mixed. Outputs and RGB rows are fresh.
Matrix orientation adapters and higher-order rotations are separate work.
Extreme coefficients whose intermediate results are unrepresentable are
outside the ordinary finite-input accuracy contract.

## Verification

Base: b6d503adb73a322bec8ddfcd6b3ed7a96b7a2322, merged PR #28.
Baseline: 1687 passed, 5 xfailed. Final: **1752 passed, 0 failed, 5 xfailed**,
Python 3.12.14 / pytest 9.1.1. Exact expected-failure identities are unchanged.
See `phase2f5c-test-results.json`.

65 new tests verify independent Cartesian SH polynomials evaluated at inverse
Rodrigues-rotated directions, axis quarter/half turns, arbitrary axes,
noncommuting Hamilton composition, inverse recovery, L0 invariance, band
energy, RGB independence, ownership, sign equivalence, legacy conversion and
commutation with diffuse convolution. Tolerance tests locate adjacent binary64
norms on either side of both bounds; they do not assume decimal 1+1e-12 is
exactly representable. Invalid/degenerate cases remain explicit.

A built wheel is actually installed into an isolated target and passes a
known active Z quarter-turn computation. Representative scalar/RGB L2 timings
are recorded in `phase2f5c-benchmarks.json`: median of five 10,000-call trials,
including validation and fresh storage. Work and temporary storage are bounded
for the supported bands; no HDR reprojection cost is incurred. Measurements
are local CPU timings, not application performance guarantees.

No experimental, projection, sampling, quaternion or matrix algorithm is
changed. No mandatory dependency or publication is introduced.
