# API questions and recommended Phase 2 order

This is a review checklist, not an implementation plan already in progress. The retrieved wiki resolves vector axis identities and several ordinary API requirements; it does not supply the missing policies below. Phase 2 requires a new instruction.

## Questions requiring approval before their behavior changes

| Question | Existing evidence | Decision needed |
| --- | --- | --- |
| QD01: Supported dimensions and errors | Wiki requires equal dimensions for operators and describes Vector sizes 1…N; source accepts empty/malformed shapes | Is Vector0 supported? Should mixed-size equality return false or reject? What exception should arithmetic mismatches use? Preserve equal-size component semantics. |
| QD02: Input ownership | Wiki distinguishes new result from changed receiver, but says nothing about other mutable inputs | Should returning clamp preserve its separate value list? Should axis conversion and Ray construction copy or mutate caller inputs? Must ray duplication preserve distance/end and fully isolate storage? |
| QD03: Projection API | Wiki requires project to return a 3D Vector; source uses lists for project and wrappers for unproject | Should both functions accept Matrix wrappers, nested lists, or both? Confirm window Z=[0,1], viewport origin, homogeneous input rules, and zero-W behavior. Three-component return is already documented. |
| QD04: Ray dimensional policy | Ray page is a placeholder; Matrix multiplication requires matching dimensions | Should ray methods explicitly support Matrix4 with Vector3 via internal position/direction promotion, or require matching input representations? Define distance/end semantics. Do not add general implicit promotion to Matrix*Vector without approval. |
| QD05: Planes | Placeholder page; scalar consumers use n·p+d=0 but constructors/bestFitD differ | Confirm scalar coefficient representation and d sign; raw versus unit normal field; bestFitD's meaning; degenerate points/polygons and open versus repeated-closure polygon inputs. |
| QD06: Quaternion domains and returns | Placeholder page; source implies unit rotations; helper docstring/implementation disagree | Confirm general versus unit-only inverse/log/power/rotation domains, arbitrary-axis helper return semantics, and existing log list return. Which three-control spline behavior should legacy SQUAD represent? No signature replacement is approved. |
| QD07: Accuracy and degeneracy | No wiki tolerances; SLERP deliberately uses a nearby linear approximation | Should SLERP guarantee unit norm and at what tolerance? Define zero normalization, NaN/Inf, singular/near-singular inversion, antipodal interpolation and invalid frusta policy. Preserve LERP's existing component/extrapolation behavior unless explicitly changed. |
| QD08: Mixed axes and units | Vector −Z front is explicitly documented; Quaternion +Z forward and mixed degrees/radians are established source behavior | Retain both forward helpers and mixed units with documentation, or add new explicit alternatives? No sign/unit unification should occur by default. Angle conversion helper names can be corrected separately. |
| QD09: 2D rotation, transform, and shear | Wiki names methods but does not specify pivot, argument representation, or shear displacement source | Define 2D pivot versus axis contract, transform homogeneous behavior, and shear-plane mapping. Preserve row-vector storage/conventions. |
| QD10: Refraction and viewport | Wiki lists refraction; common viewport page is a placeholder | Confirm IOR ratio direction, unit-vector prerequisites and TIR output. Define what getViewPort is meant to compute instead of choosing a formula from its name. |
| QD11: Experimental scope | No experimental wiki; transport incomplete and irradiance map assumes square angular probes | Support rectangular probes or reject them clearly? Confirm raw float endianness and SH sign basis interoperability. Decide usable algorithm scope without inventing shadow-transport behavior. |
| QD12: Public numeric/export protocol | Scalar examples use float; ctypes export and sequence access are undocumented | Preserve c_matrix snapshot/export? Add integer/scalar reflected operators or indexing as additive APIs? Decide before widening protocol requirements. Do not add mandatory compiled dependencies. |

The default Matrix identity is **not** a proposal to change behavior: it is preserved, supported by contemporaneous source despite a wrong wiki output. Likewise, `i_*` returning self is compatible with “without returning a new object”; no change to None is recommended. `indetity` is a documentation typo, not an API requiring a compatibility alias.

### Phase 2C decision update

The user explicitly approved the refraction portion of QD10 before V02 implementation: IOR=n1/n2 (incident/transmitted indices), normalized incident and normal vectors with matching dimensions, normal opposing incidence and pointing into the incident medium, and preservation of the historical zero-Vector result for total internal reflection. No automatic normalization or normal flipping is added. See [current conventions](CONVENTIONS.md#phase-2c-current-angle-and-refraction-conventions). The viewport portion of QD10 and all other unresolved questions remain pending; approval of refraction does not resolve them.

### Plane representation decision

QD05's representation is settled: scalar `a,b,c,d` satisfy `n·p+d=0`, and `.normal` stores `[a,b,c]` at coefficient scale. Three-point construction produces unit coefficients and a negative dot-product offset; normalization scales all four coefficients and synchronizes `.normal`. Polygon normals wrap edges and support repeated-first closure. `bestFitD` retains signed geometric D, with coefficient `d=-D`. See [plane conventions](CONVENTIONS.md#plane-representation). Broader degeneracy, validation, and numerical error policies remain unresolved.

### Vector transformation decision

QD09's transform portion is settled: same-size inputs use row-vector multiplication; an N-component position with an (N+1)×(N+1) matrix receives local w=1 promotion and returns N components. Explicit homogeneous inputs compute every output component normally, including w; no perspective divide is performed. Implicit promotion is intended for affine positions. Raw position and nested matrix lists, returning versus receiver mutation, and the general multiplication operator's matching-dimension requirement are preserved. See [transformation conventions](CONVENTIONS.md#vector-transformations). Pivot rotation, shear mapping, projection/unprojection, and unsupported-shape error policy remain separate questions.

### Projection and unprojection decision

QD03 is settled for ordinary valid inputs: both functions accept Matrix4 wrappers or raw 4×4 nested lists, including mixed forms; project requires explicit Vector4 input. The viewport uses `[x,y,width,height]`, lower-left origin, and upward Y. Window depth is `(ndc.z+1)/2` without clamping. Both functions return fresh Vector3 values and preserve inputs. Row-vector composition is modelview*projection, with homogeneous division after forward or inverse transformation. Project raises ZeroDivisionError for zero clip W; unproject preserves its zero-Vector3 output-W sentinel and singular-inverse ZeroDivisionError. Projection does not require invertibility. See [projection conventions](CONVENTIONS.md#projection-and-unprojection). Unsupported-shape validation, invalid viewports, near-zero thresholds, and broader numerical policies remain separate questions.

### Pivot rotation and shear decision

QD09's remaining mathematical mappings are settled. `rotate2(point,theta)` uses a 2D pivot, degrees, and positive counterclockwise row-vector rotation; its 3×3 homogeneous matrix implements `pivot+(p-pivot)R`. Matrix2 origin-only rotation and Matrix3/Matrix4 axis-angle dispatch are preserved. Shear XY adds x*X+y*Y to Z; YZ adds y*Y+z*Z to X; XZ adds x*X+z*Z to Y. Size-3 forms operate on XYZ, and size-4 forms preserve W. Existing YZ/XZ mappings, signatures, argument order, postmultiplication, and ctypes synchronization are preserved. See [pivot/shear conventions](CONVENTIONS.md#pivot-rotation-and-shear). Broader unsupported-input and numerical policies remain separate questions.

### Quaternion/matrix conversion scope

Q02/Q07 are corrected within the established unit-rotation domain: [w,x,y,z], proper row-vector rotation matrices, sign-equivalent q/-q orientations, and synchronized Matrix4 exports. No implicit normalization, orthogonalization, canonical-sign requirement, or new input representation is introduced. Existing zero/nonunit quaternion-to-matrix arithmetic is preserved. Other quaternion-domain, return, validation, and numerical questions in QD06/QD07 remain unresolved. See [conversion conventions](CONVENTIONS.md#quaternionmatrix-conversions).

### Quaternion axis ownership

QD02 is settled for `quat_from_axis_angle` and `quat_rotate_from_axis_angle`: normalize a temporary axis and preserve caller Vector/list storage and values. Both retain degrees and existing supported representations; the latter retains its legacy normalized-axis Quaternion sandwich result, leaving Q11 unresolved. Other ownership questions remain separate. Q09 exposes the legacy float-only division protocol to Python 3 without extending QD12's scalar/reflected-operator policy. See [axis/arithmetic conventions](CONVENTIONS.md#quaternion-axes-and-arithmetic).

### Quaternion power and logarithm domains

QD06/QD07 are settled for Q03/Q04: unit inputs without normalization or a
new norm tolerance; finite real powers; principal atan2 imaginary/scalar
angle; fresh Quaternion powers and fresh four-element list logarithms with
zero scalar. Zero raises ValueError. Negative identity permits integer
powers by parity and rejects fractional powers/logarithms. Tiny imaginary
directions are retained without a cutoff; quaternion signs are preserved.
General nonunit domains, broader accuracy/invalid-input policies, SQUAD,
SLERP and the arbitrary-axis helper's Q11 return contract remain separate.
See [power/log conventions](CONVENTIONS.md#quaternion-powers-and-logarithms).

### Quaternion interpolation contracts

QD06/QD07 are settled for Q05/Q10: preserve the legacy three-control nested
no-invert blend and its sign-sensitive approximation/degeneracy behavior;
use accurate shortest-path spherical SLERP for unit inputs, with 1e-12 norm
and known-axis component accuracy over [0,1]. No input normalization or
parameter clamping is introduced. Exact half-turn ties retain supplied sign
branches. `squad4(q0,q1,s0,s1,t)` is a separate conventional four-control
blend using accurate SLERP and explicit unit controls. Control generation,
nonunit/nonfinite domains and antipodal no-invert policy remain separate.
See [interpolation conventions](CONVENTIONS.md#quaternion-interpolation).

### Quaternion API compatibility

Q11 is settled by preserving `quat_rotate_from_axis_angle` as the legacy
pure Quaternion axis-rotation result and recommending the existing
`quat_from_axis_angle` constructor for conventional rotations. Vector/list
support, degrees, caller ownership, numerical sandwich and signatures are
retained. No redundant constructor or runtime warning is added. Public API
documentation separates return types, units, forward axes and domains;
Q01-Q11 are resolved within the previously defined scopes. General numerical
robustness, malformed/nonfinite inputs, nonunit rotation/power/log policies
and no-invert antipodal behavior remain separate. See the
[quaternion API guide](../docs/QUATERNIONS.md).

### Ray copying, rotation and translation

QD02/QD04 are settled for ordinary 3D ray copying and rigid transforms:
retain constructor references/in-place direction normalization; duplicate
all stored Vectors independently and copy distance exactly without
construction; rotate about the coordinate origin with unit-quaternion
conjugate sandwiches; retain existing Matrix3 rotation; locally promote
Vector3/Matrix4 positions/directions with w=1/w=0 for pure translation.
Public signatures, None transform returns and general Matrix*Vector rules
are unchanged.

Ray `.end` remains the historical intersection placeholder/state. It is
copied by duplication and left unchanged by transforms. A zero placeholder
cannot be distinguished from a valid hit at the origin without a validity
contract; nonzero values do not establish validity either. Hit-state
validity, ownership and transformation are a separate unresolved API design
question. No geometric endpoint is assigned automatically. Scale/shear,
projective and wider invalid-input policies remain separate. See the
[ray guide](../docs/RAYS.md).

### Finite numerical robustness

QD07 is settled for N01-N03: direct zero Vector/Quaternion normalization
returns zero/identity; narrow exact-zero core-caller guards preserve geometric
ZeroDivisionError; finite norms/normalizations use scaled hypot; overflowing
norms may be infinity. 3x3/4x4 cofactor inverses use power-of-two scaling,
exact binary64 singularity checks and signed infinity for exponent-rescaling
overflow. No epsilon/condition cutoff is introduced. Legacy NaN/Inf input
paths, 2x2 inversion, public determinants and experimental algorithms remain
unchanged. Severe conditioning, unrepresentable results and broader invalid
input policies remain outside the accuracy guarantee. Experimental zero
normalization inherits the direct fallback without new caller semantics.
See [numerical conventions](CONVENTIONS.md#numerical-robustness).

## Recommended implementation order after review

1. **Small ordinary-math fixes with clear contracts:** equal-size equality/inequality (V01 component cases), inverse2 (M01), quaternion inverse (Q01). Preserve exact comparisons, return types, and component order; defer unsupported-dimension policy changes.
2. **Coupled division fix:** M02/M03/M07, updating inverse4's current transpose dependency together. Preserve legacy special methods, add Python 3 division, and check both inverse identities plus ctypes state.
3. **Named angle conversions and refraction:** C01, then V02 once IOR/prerequisite policy is accepted. Retain every other API's established units. Avoid blanket reflected-operator refactoring just to make local formulas work.
4. **Representation/state fixes after decisions:** G01–G04 after QD05; M04 translation consistency; Q07/Q08 and Q06 runtime issues with explicit ownership/type contracts. Keep commits scoped to one defect or tightly coupled cause.
5. **Transform/projection fixes after QD03/QD09:** V03/P01/P02 and M05/M06. Preserve row-vector composition, documented 3D project result, and explicit operator dimension restrictions. Verify noncommuting transforms.
6. **Ray and remaining quaternion work after QD02/QD04/QD06/QD07:** R01/R02/R03; Q02 matrix conversion; Q03/Q04 powers/log; Q05 SQUAD. Contract-dependent Q10/Q11 changes require explicit approval rather than being bundled with crash fixes.
7. **Numerical robustness after domain/error policy:** N01–N03, SLERP accuracy if approved, and conditioning/overflow behavior. Use meaningful bounds and pure-Python implementations; avoid a broad inversion redesign without measured need.
8. **Experimental algorithm repairs last:** E01/E02/E03 Bezier scalar/operand/Python-3 fixes, then independently reproduce masked E04 issues; E05 recurrence, E06 generation, and explicitly scoped E07/E08 work. Do not expand unfinished shadow transport implicitly.
9. **Only after correctness is established:** propose Phase 3 packaging, supported Python versions, CI, documentation/typing, and performance changes. Existing benchmarks stay as baseline; no PyPI publication or merge is included.

Each later review unit must show the pre-fix failure, retain unrelated xfails/questions, remove only the relevant strict defect markers, and include a compatibility note. Wiki examples with mistakes should receive documentation-only corrections separately from mathematical corrections. Existing API/import names and pure-Python design remain constraints throughout.

## Experimental retirement roadmap

Validated Bezier evaluation is supported under `gem.bezier`; experimental
imports retain compatibility reexports and legacy sampling extensions.
E04 remains separate from E01–E03. Later validated Legendre, spherical
harmonics sampling and irradiance work will migrate to coherent core modules.
Incomplete shadow transport is reviewed separately. Final removal of the
experimental directory requires a dedicated cleanup after all migrations and
compatibility decisions; no unrelated modules move in Phase 2F-3A.

## Bezier adaptive sampling and builders

E04 is supported in core with midpoint subdivision, squared-distance tolerance
and maximum depth 16. The finite-chord flatness criterion handles coincidence
and collinear overshoot; capped output is best effort. Nested path output and
ordered endpoints are preserved. `interpolate` remains append-only;
`samplePoints` rebuilds from ordered source vertices using distinct squared
thinning thresholds. Experimental modules reexport the core class. See
[contract and verification](PHASE2F3B.md).

## Legendre core support

E05 is corrected under `gem.legendre` with a compatibility class reexport.
The historical unnormalized Condon–Shortley convention and existing real SH
basis are retained. Supported associated inputs are integer 0 <= m <= l and
x in [-1,1]. `run` is state-preserving; named scratch helpers remain explicit
mutators. Invalid/extreme domains receive no new generalized policy. See
[verification and compatibility](PHASE2F4.md).

## Spherical-harmonics sampling and probe contracts

E06/E08 are supported under gem.spherical_harmonics with transitional imports.
QD11's probe contract is angular-disk pixel-center quadrature stretched to
rectangular images, native-endian raw float32 RGB and explicit legacy-basis
conversion. No latitude-longitude or mirrored-ball interpretation is added.
Canonical radiance projection is separate from first-three-band diffuse
convolution and unit-direction reconstruction. Global RNG stratification is
retained. Recalculation rebuilds coefficients; direct updates accumulate.
E07 transport and coefficient rotation remain separate. See
[API and numerical limits](../docs/SPHERICAL_HARMONICS.md).
