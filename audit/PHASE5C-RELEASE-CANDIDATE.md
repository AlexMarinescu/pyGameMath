# Phase 5C — Release candidate verification

Base: master `a21273c857e6ed9e948fa92246226100efb80355` (PR #62).
Branch: `release/phase5c-candidate`. Target: gem 1.0.0, GitHub-first.

## Windows reference investigation

An unchanged-baseline trace identifies the first difference in `math.sin` of
`0x1.921fb54442d18p-1` (the binary64 approximation of pi/4). Linux/macOS return
`0x1.6a09e667f3bccp-1`; Windows returns `0x1.6a09e667f3bcdp-1`, one ULP higher.
A 90-digit Decimal series using the exact input float gives
0.707106781186547502751942956217516746261543239537492789524366119137482021518043926217856922.
The Linux/macOS result is nearest; Windows differs by one ULP. Inputs and all
other traced math calls match on the representative CPython 3.12 jobs.

This changes SLERP's scalar weight for t=0.25, then the quaternion's scalar
component from `0x1.ee8dd4748bf15p-1` to `0x1.ee8dd4748bf16p-1`. The Hamilton
sandwich computes X as w*w-z*z: the rotated axis changes from
`0x1.bb67ae8584cabp-1` (0.8660254037844387) to
`0x1.bb67ae8584cadp-1` (0.8660254037844389), two ULP apart.
The mathematical result sqrt(3)/2 is
0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093;
its nearest float is `0x1.bb67ae8584caap-1`. The two results are respectively
one and three ULP above that reference. Absolute error remains below 2^-51.
The same sine call affects the t=0.75 imaginary component and the small
cancellation residual in the 90-degree rotated axes.

This is a benign platform libm difference and an overly strict full-precision
JSON serialization requirement, not a defect in the frozen mathematics.
Verification now permits differences no larger than 2^-52 between only eight
audited quaternion measurement fields. Both manifests must independently satisfy
closed-form values within 2^-51. SVG bytes, all other fields, key/type structure,
serialization, metadata and artifact hashes remain exact. Existing numerical
assertions and committed references are unchanged. Regression tests replay the
specific libm difference and reject larger, unrelated and shared numerical errors.

## Verification status

Final platform validation, artifact inventory and publication prerequisites are
being collected. No readiness claim applies until the completed report records
all required evidence. No tag, release, package publication or merge is performed.
