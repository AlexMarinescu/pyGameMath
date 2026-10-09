# Maintainer decisions and documented limits

This reference records existing behavior; it does not resolve a new runtime
contract. The [central development decisions](../development/roadmap.md#open-decisions)
and [compatibility policy](../architecture/compatibility.md) remain authoritative.
Settled row-vector, unit, basis, ownership and legacy-helper choices are not reopened.

| Area | Verified behavior / remaining question | Recommended separate review |
|---|---|---|
| Malformed sizes/storage | Constructors retain raw storage; many operations index directly, ignore extra entries or fail natively | Define shape/dimension validation and numerical protocols before uniform exceptions or reflected operators |
| Sequence protocols | Vector/Matrix/Quaternion wrappers lack indexing/iteration; storage lists supply it | Decide any additive protocol separately; do not document a nonexistent interface |
| Clamp explicit size | i_clamp replaces storage but retains receiver size even if explicit size differs | Review malformed-size policy; normal use supplies matching size |
| Matrix2 rotation | Wrapper extracts the origin-only linear block of the homogeneous pivot builder through declared-size multiplication | Any future pivot overload or public builder-shape change needs a separate compatibility decision |
| Translation dimensions | translate3 is last-row replacement, and returning/in-place wrappers have different unsupported-argument paths | Specify unsupported dimensions and validation before changing exceptions or reinterpreting translate3 |
| Nonfinite/extreme inputs | Stable finite norms and Matrix3/4 scaled inverses are narrow guarantees; legacy paths retain native arithmetic | Define any nonfinite policy or broader conditioning strategy separately, without arbitrary epsilon rejection |
| Unit Quaternion powers/log | Unit-only contract is a prerequisite, not enforced normalization/tolerance; general nonunit results unsupported | Decide general-domain support separately, preserving current branch/return contracts |
| SH unit direction | reconstruct checks finite/nonzero XYZ and Z range, not full norm | Do not claim norm validation; any broader direction validation needs compatibility review |
| High-order Legendre/SH | Integer-domain recurrence and ordinary tests are established, but extreme-order stabilization is absent | Define supported numerical range and overflow/error policy before guarantees |
| Plane field synchronization | Constructors/normalization synchronize coefficients and normal; direct field edits can stale the snapshot | Review property/update mechanisms separately; no invented automatic tracking |
| Ray end state | Retained placeholder/hit state has no validity flag, and transforms leave it unchanged | Define hit validity, ownership and transformation before adding intersections |
| Bezier mutable fields | Getter exposes controls; caller edits can stale curveCount; retained fields do not all control adaptive sampling | Consider explicit rebuild/property policies independently, preserving append-only interpolation |
| Alias removal/stability | Experimental shims remain importable, not deprecated by a new deadline | Choose release boundary/window and public-stability policy before removing paths |
| Release/support/docs distribution | Historical version/classifiers are not modern support evidence; only CPython 3.12/Linux executed here | Release engineering must choose versions, matrix and doc inclusion; no packaging changes here |

Matrix2 origin-only rotation is established behavior, not an implicit new pivot
API. Invalid inputs are not turned into universal exception guarantees.

See [coverage and verification](coverage.md), [API index](index.md) and
[canonical roadmap](../../ROADMAP.md) for boundaries.
