# Existing mathematical conventions

Baseline: `5257291431bb45db0274dc48edf24694ecfe2e2d`. These describe the observed implementation, not a redesigned API. No convention was changed in Phase 1 or Phase 1B. Historical wiki validation uses snapshot `715e5039c75e080814a12e957f5148c35cdf8bda`; see [complete reconciliation](PHASE1B-WIKI.md).

## Matrices and vectors

`Matrix.matrix[i][j]` is a nested list addressed as row `i`, column `j`. `matrix_multiply` implements the ordinary product `C[i][j] = sum(A[i][k] B[k][j])` ([matrix.py:17](../gem/matrix.py#L17)). Flattening and the ctypes buffer preserve row-major list order ([common.py:21](../gem/common.py#L21), [common.py:34](../gem/common.py#L34)). This is storage order, not a claim that all consumers use the same interpretation.

Despite the expression `M * v`, the vector kernel computes `out[j] = sum(v[i] M[i][j])`: mathematically **row-vector multiplication `v M`** ([matrix.py:28](../gem/matrix.py#L28)). A concrete probe is `[[1,2],[3,4]] * [5,6] -> [23,34]`, not `[17,39]`. Accordingly, `(A * B) * v` means apply `A`, then `B`; it equals `B * (A * v)` in the library's syntax. Matrix products themselves are ordinary products.

3D homogeneous translation occupies the final row: `translate4([tx,ty,tz])` returns rows `[..., [tx,ty,tz,1]]` ([matrix.py:119](../gem/matrix.py#L119)). Use four-component positions `[x,y,z,1]` and directions `[x,y,z,0]`. There is no automatic promotion of a Vector3 to Vector4. `translate2` similarly uses a 3×3 homogeneous matrix. `translate3` is a separate, ambiguous 3×3 operation that overwrites the last row, including its diagonal; it cannot encode general 3D affine translation.

For OpenGL column-vector consumers, the contiguous row-major bytes describe the transpose when interpreted as column-major data. This can represent the equivalent transform, but the upload call's transpose flag and matrix composition must agree. No external integration was assumed or tested. `c_matrix` is a float32 snapshot, whereas `.matrix` arithmetic uses Python numbers; it is not a live view and direct list mutation can leave it stale.

## Coordinates, orientation, and projections

Cross products use the right-handed formula: X cross Y = Z ([vector.py:43](../gem/vector.py#L43)). Positive rotation around +Z sends +X toward +Y for the row-vector representation ([matrix.py:144](../gem/matrix.py#L144)). Vector helpers identify +X as right, +Y as up, and **−Z as front** ([vector.py:413](../gem/vector.py#L413)). These support a right-handed, negative-Z viewing convention. The library does not enforce one universal world-coordinate convention.

`lookAt` builds a right-handed view matrix with a negative-Z viewing direction ([matrix.py:657](../gem/matrix.py#L657)). `perspective`, `perspectiveX`, and `orthographic` use OpenGL-style NDC depth **[−1,+1]**; points at camera Z = −near and −far map to −1 and +1 respectively. `perspective` takes vertical FOV, `perspectiveX` horizontal FOV; aspect = width / height. FOV is in degrees.

`unproject` expects window depth [0,1], mapped by `2*z−1`, and window X/Y mapped relative to `[x,y,width,height]`. Its current multiplication order is incorrect for noncommuting model/projection matrices. `project` is unusable; its unreachable return also leaves NDC depth unremapped and declares a size-3 vector with four values. The crash and size-3/four-value mismatch are defects. Phase 1B confirms a three-dimensional return from the wiki; window depth [0,1] is a proposed compatibility choice consistent with unproject, not an explicit wiki specification. Project input types and homogeneous behavior require approval.

Quaternion `getForward()` uses **+Z**, whereas Vector `front()` uses −Z ([quaternion.py:418](../gem/quaternion.py#L418)). Preserve both until a compatibility policy is approved; the mismatch must be explicit in documentation.

## Quaternions and angles

Quaternions are **[w,x,y,z]**, identity `[1,0,0,0]`, with Hamilton multiplication ([quaternion.py:6](../gem/quaternion.py#L6), [quaternion.py:18](../gem/quaternion.py#L18)). Vector rotation is `q * (0,v) * conjugate(q)` for unit `q`, implemented by `quat_rotate_vector`. `q * Vector` alone returns a Quaternion product, not a rotated Vector. Rotation preserves length only for unit quaternions; nonunit inputs scale the rotated vector by the squared norm. Rotation-matrix conversion also assumes unit quaternions.

`toMatrix()` uses a row-vector matrix consistent with `Matrix.rotate`; conversion round trips must allow `q` and `−q` to represent the same rotation. Quaternion composition applies the right operand first in `q1*q2`; the equivalent row matrices appear in reversed order.

Angle units vary:

| API | Observed unit |
| --- | --- |
| `rotate2`, `rotate3`, `rotate4`, `Matrix.rotate` | degrees |
| `rotate_origin2` | radians |
| `quat_from_axis_angle`, `quat_rotate`, `quat_rotate_from_axis_angle` | degrees |
| `quat_rotate_x/y/z_from_angle` | radians |
| `toAngle` | returns radians |
| `SPH(theta, phi)` / sample theta, phi | radians |

`common.radiansToDegrees` and `degreesToRadians` implement the opposite conversions using 3.14, contrary to their names. Keep the established angle units in the other APIs; do not globally switch to radians.

## Planes, rays, interpolation, and spherical harmonics

`Plane.dot` and `point_location` express `a*x+b*y+c*z+d=0`, with positive values on the normal side ([plane.py:82](../gem/plane.py#L82), [plane.py:115](../gem/plane.py#L115)). `fromPoints` instead stores points in the coefficient fields and sets `d=normal.dot(point)`, incompatible with that equation. `bestFitD` returns the positive average dot product, corresponding to `normal.dot(point)=D`; its sign differs from coefficient `d`. This is inconsistent representation, not evidence to silently choose a new equation.

`Ray` stores a mutable origin, normalizes the caller's direction in place, and records its original magnitude as `distance`. The stored `end` remains a zero vector: intersections are unfinished. Preserve the misspelled public method `roateUsingMatrix`; any corrected spelling should be an additive alias. Translation should distinguish homogeneous positions and directions, but currently changes only the direction and fails on common 4×4 inputs.

Vector and scalar LERP use `a+t*(b−a)` without clamping `t`, so extrapolation is supported. Quaternion LERP is linear component interpolation, **not** normalized LERP. SLERP uses shortest-path sign correction, with an unnormalized linear approximation for nearby inputs. `slerp_no_invert` deliberately omits sign correction; its antipodal midpoint is the zero quaternion, an ambiguous rotation rather than a unique mathematical answer. SQUAD's signature has three quaternions, unlike the usual four-control implementation: choosing its intended spline contract requires approval.

`Legendre` and `SPH` use associated Legendre polynomials with the **Condon–Shortley phase** and real SH: `m>0` cosine terms, `m<0` sine terms, `m=0` zonal terms. `theta` is polar angle from +Z, `phi` azimuth from +X toward +Y. SH sample index is `l*(l+1)+m`, number of coefficients = bands². Low orders 0–2 pass orthonormality checks. Higher-order recurrence is defective.

The irradiance map's nine hard-coded SH polynomials instead use positive X/Y first-order terms, which differ in sign from `SPH(1,±1,...)`. Its intended input is raw native-endian 32-bit RGB floats, not a general HDR decoder; the integration assumes a square angular light probe. Rectangular images crash. Basis/sign interoperability and file-endian semantics need documentation before mathematical corrections.

## Historical wiki evidence added in Phase 1B

All seven current wiki pages are archived in [wiki-snapshot](wiki-snapshot/manifest.json). Vector and Matrix pages contain API lists/examples; Quaternion, Plane, Ray, and Common Functions remain placeholders, including their earlier revisions. No wiki page documents experimental functions.

The wiki explicitly confirms Vector front=−Z, back=+Z, right=+X, up=+Y; same-size operands for vector operators; matching matrix/vector dimensions for multiplication; vertical versus horizontal FOV APIs; project returning a 3D Vector; and the ordinary/new versus i-prefixed/in-place receiver distinction. “Without returning a new object” permits returning self. It does not specify ownership of a separate mutable value/axis argument.

The wiki does **not** establish storage order, multiplication/composition orientation, quaternion components, angle units, window depth, plane offset sign, ray homogeneous promotion, unit versus general quaternion domains, or numerical tolerances. All corresponding source observations above remain observations rather than newly approved contracts. The equivalent transposed column-vector interpretation of contiguous data does not authorize changing observable multiplication.

Matrix default identity remains established source behavior. The wiki prints zeros, but the contemporaneous 2015 constructor already used identity; its output is a documentation error. Likewise, the wiki's in-place scale output incorrectly shows an unchanged identity while reporting determinant 24. Neither mistake justifies changing constructors or in-place mutation.

Decision-sensitive Phase 1 expectations are separated into contract-question markers: mixed/empty equality, clamp value-list copying, project argument/depth choices, ray Vector3/Matrix4 promotion, SLERP unit tolerance, and arbitrary-axis quaternion return semantics. See [decision list](PHASE2-DECISIONS.md). Core equal-size equality, inverse identities, and dimensional consistency keep their confirmed defect regressions.

## Phase 2C current angle and refraction conventions

The angle-helper defect described in the baseline above is corrected: `radiansToDegrees` multiplies by `180/math.pi`, and `degreesToRadians` multiplies by `math.pi/180`. Their original keyword parameter names are retained for compatibility. Every other API's angle units in the table above remain unchanged.

The user explicitly approved the following refraction contract for Phase 2C; it is a new recorded decision, not a claim that the historical wiki specified it. `refract(IOR, incidentVec, normal)` uses **IOR = n1/n2**, where n1 is the incident medium's refractive index and n2 is the transmitted medium's. Both vectors must be normalized and have matching dimensions. The incident vector points toward the interface; the normal points into the incident medium and opposes the incident direction (`normal.dot(incidentVec) <= 0`). The caller supplies this orientation and normalization; the function does not normalize, flip normals, or infer which medium the ray is entering.

For air-to-glass with illustrative indices 1.0 and 1.5, pass `IOR=1.0/1.5` and a normal pointing into air. For glass-to-air, pass `IOR=1.5/1.0` and a normal pointing into glass. When reversing travel across the same interface, reverse the normal and exchange the indices. Always pass the ratio rather than the transmitted medium's index alone.

With `d=normal.dot(incidentVec)` and `k=1-IOR*IOR*(1-d*d)`, the transmitted direction is `IOR*incidentVec-(IOR*d+sqrt(k))*normal`. This obeys Snell's law `n1*sin(theta1)=n2*sin(theta2)` for unit inputs. If `k<0`, preserve the historical **fresh zero Vector** result with the input dimension for total internal reflection; this is a sentinel, not a reflected direction. If `k=0`, the transmitted direction is tangent to the interface. Neither input is mutated. The existing floating-point `k<0` comparison is retained, without a new near-critical clamping tolerance or unsupported-input policy.

## Plane representation

Planes use scalar coefficients in `a*x+b*y+c*z+d=0`. `.normal` contains `[a,b,c]` at the same scale as the coefficients. `fromCoeffs` preserves supplied values; `fromPoints` uses the unit cross product `(b-a) cross (c-a)` and sets `d=-normal.dot(a)`. Reversing point order reverses orientation. Normalization divides all four coefficients by the original normal magnitude and refreshes `.normal`. For a unit-normal plane, evaluating the equation gives signed perpendicular distance; positive values lie on the normal side.

`bestFitNormal` computes a unit Newell normal from ordered polygon vertices, wrapping the final edge to the first vertex. Open lists and lists with a repeated first vertex are supported. `bestFitD` retains `D=average(normal.dot(point))`, representing `normal.dot(point)=D`; D can be negative. Construct a coefficient plane using `d=-D`:

```python
p = plane.Plane()
n = p.bestFitNormal(vertices)
D = p.bestFitD(vertices, n)
p.fromCoeffs(n.vector[0], n.vector[1], n.vector[2], -D)
```

For planar polygons this plane contains every vertex. For nonplanar input, the Newell normal and mean offset describe an approximation, not a least-squares fit. Each supplied vertex contributes to the mean, including a repeated endpoint. Constructors retain their None returns; helpers do not mutate the receiver or input vertices. Public fields remain mutable snapshots: changing coefficients or `.normal` directly does not automatically synchronize the other representation.

The historical plane wiki contains only “Coming soon.” These conventions resolve the representation inconsistencies identified in G01–G04. Zero-normal normalization, collinear three-point construction, and empty best-fit helpers retain their existing ZeroDivisionError behavior. Validation and exception policy for malformed, degenerate, nonfinite, or extreme-scale inputs remain unresolved under QD05/QD07.

## Vector transformations

`vector.transform(size, position, matrix)` accepts a raw position list and a square nested matrix list. `Vector.transform(position, matrix)` uses the receiver's size and returns a fresh Vector; `Vector.i_transform(position, matrix)` replaces only the receiver's vector list and returns self. The supplied position list and matrix rows are preserved, including when the position list is the receiver's old storage. The free helper returns a fresh list. Wrapper arguments are not added to these APIs.

For a size N position and an N×N matrix, transformation is the ordinary row-vector product: `out[j]=sum(position[i]*matrix[i][j])`. An explicitly supplied homogeneous component participates normally. For example, `[x,y,z,0]` is a direction under affine Matrix4 transforms, and `[x,y,z,w]` receives translation weighted by w. Output w is computed from the matrix's final column; it is not forcibly preserved.

For an N-component position and an (N+1)×(N+1) matrix, the transform helper locally supplies w=1, computes `out[j]=sum(position[i]*matrix[i][j])+matrix[N][j]`, and returns N components. This promotion is intended for **affine position transforms**. No automatic perspective divide occurs. For a projective matrix, the returned N components are homogeneous numerators rather than perspective-correct Cartesian coordinates. Use the [projection/unprojection path](#projection-and-unprojection) when those coordinates are needed. Use explicit homogeneous input when the computed output w is needed.

Row-major storage, translation in the final row, and ordinary matrix products are unchanged. In row-vector order, `T*R` applies translation then rotation; `R*T` applies rotation then translation. In library operator syntax, `(T*R)*v` equals `R*(T*v)` for matching dimensions. Local promotion in the transform helper does not extend general `Matrix*Vector`: its established matching-dimension precondition remains, and mismatched Vector3/Matrix4 multiplication still raises IndexError.

`Matrix(4).i_translate(Vector3)` now builds a 4×4 translation matrix, matching the returning method and refreshing the float32 ctypes snapshot through existing in-place multiplication. Vector4 offsets retain their existing ignored fourth component. Matrix3/Vector2 translation remains homogeneous 2D translation; Matrix3/Vector3 retains the legacy last-row replacement helper, not general affine 3D translation. Matrix2 translation remains unsupported. Shape validation and unsupported-input exception policies remain unresolved; this phase defines ordinary valid transformation shapes without introducing a broader validation API.

## Projection and unprojection

Both functions accept Matrix4 wrappers and raw 4×4 nested lists, including mixed representations. `project(obj, model, proj, viewport)` requires an explicit Vector4 `[x,y,z,w]`; a Cartesian position normally uses w=1. Other supplied w values participate in multiplication. Vector3 promotion and new object-input forms are not supported. Both functions return a fresh Vector3 with exactly three stored components and preserve caller vectors, matrix rows, ctypes snapshots, and viewport data.

The viewport is `[x,y,width,height]`, with a lower-left origin and upward-increasing Y. For row vectors, clip coordinates are `obj * (model * proj)` mathematically. Divide clip X/Y/Z by the computed clip W to obtain NDC, then map:

```text
winx = viewport.x + (ndc.x + 1) * viewport.width / 2
winy = viewport.y + (ndc.y + 1) * viewport.height / 2
winz = (ndc.z + 1) / 2
```

OpenGL NDC depth [-1,1] becomes window depth [0,1]; near/far map to 0/1 for the existing perspective and orthographic cameras. No coordinate or depth clamping occurs, so points outside the frustum can produce values outside the viewport or depth interval. No top-left-origin flip or graphics-driver state is inferred.

`unproject(winx, winy, winz, modelview, projection, viewport)` reverses the viewport mapping, uses NDC `[x,y,2*winz-1,1]`, and applies `(modelview * projection).inverse()` in row-vector mathematics. Divide the resulting X/Y/Z by the computed homogeneous W to obtain object Cartesian coordinates. This composition is essential for noncommuting modelview and projection matrices.

Zero clip W in `project` raises ZeroDivisionError. Zero homogeneous output W in `unproject` retains the historical fresh zero-Vector3 sentinel, which cannot distinguish an invalid finite inverse image from a genuine object origin. Singular combined matrices retain the inverse routine's ZeroDivisionError. Projection requires no inverse and can map through a singular matrix when clip W is nonzero. No near-zero-W tolerance, singularity threshold, invalid-viewport validation, or nonfinite/extreme-scale policy is introduced; broader numerical and error policies remain under QD07.

## Pivot rotation and shear

`rotate2(point, theta)` takes raw pivot coordinates and an angle in degrees, positive counterclockwise. It returns a fresh 3×3 nested matrix for homogeneous row-vector positions `[X,Y,1]`. With `c=cos(theta)` and `s=sin(theta)`, its final row is `[px*(1-c)+py*s, py*(1-c)-px*s, 1]`. Thus `p' = pivot + (p-pivot)R`, or equivalently `T(-pivot)*R*T(pivot)` in row-vector composition order. The pivot stays stationary; directions `[X,Y,0]` receive rotation without pivot translation. The input pivot is preserved.

Use `Matrix(3, data=rotate2([px,py], theta))` for this 2D affine transform. Matrix2 `rotate`/`i_rotate` retain origin-only rotation and their historical Vector argument, whose pivot components do not affect the result. A 2×2 matrix cannot represent translation about an arbitrary pivot. Matrix3/Matrix4 rotation methods retain their axis-angle dispatch; no pivot overload or new method is introduced.

The shear argument names identify the coordinates that supply displacement to the remaining coordinate:

| API | Row-vector mapping |
| --- | --- |
| `shearXY(x,y)` | `Z' = Z + x*X + y*Y`; X/Y unchanged |
| `shearYZ(y,z)` | `X' = X + y*Y + z*Z`; Y/Z unchanged |
| `shearXZ(x,z)` | `Y' = Y + x*X + z*Z`; X/Z unchanged |

The size-3 helpers operate on XYZ vectors, not homogeneous 2D positions. The size-4 helpers apply the same linear mapping to XYZ and preserve the supplied W, including zero and nonunit values. Each simple shear has determinant 1 and inverse given by negating both factors. Named shears generally do not commute with other shear planes or translations.

Returning Matrix methods postmultiply the receiver by the shear matrix and return a fresh Matrix; in-place methods postmultiply, replace receiver storage, return self, and refresh `c_matrix` through existing multiplication. This leaves separately supplied rows and vector inputs unchanged. Matrix2 shear remains unsupported. Malformed inputs, nonfinite factors, extreme-scale accuracy, and general exception policies remain outside this ordinary-input contract.

## Quaternion/matrix conversions

`quat_to_matrix(quat)` and `Quaternion.toMatrix()` preserve `[w,x,y,z]` ordering and return a fresh row-vector Matrix4 with a synchronized float32 `c_matrix` snapshot. For a unit quaternion, the upper-left 3×3 block is an orthogonal rotation with determinant +1, consistent with `q*(0,v)*conjugate(q)`. The homogeneous row/column remain identity. Positive +Z rotation sends +X toward +Y. Neither the quaternion nor its component list is mutated, and each returned matrix owns separate rows and export storage.

`quat_from_matrix(matrix)` reads the upper-left 3×3 rotation block of the Matrix wrapper, returning a fresh Quaternion. Existing Matrix3 and Matrix4 extraction behavior is retained; values outside that block do not participate. Raw nested-list inputs are not added. The conversion assumes a proper rotation block: orthogonal, determinant +1. It selects the largest of the four squared-component candidates and reconstructs the remaining components using the established row-vector signs. Half-turns therefore use a nonzero axis component rather than dividing by the vanishing scalar component.

The quaternions q and -q represent the same rotation and produce the same matrix. Conversion round trips must compare sign-equivalent orientations; neither a fixed component sign nor continuity of the returned sign is promised. Tests use independent literal matrices and Rodrigues vector rotation, supplemented by round trips.

No implicit quaternion normalization or matrix orthogonalization is performed. Quaternion-to-matrix formulas for nonunit inputs, including the zero quaternion's identity matrix result, retain their existing numerical behavior; these do not establish a general rotation-domain policy. Scale, shear, reflection, malformed/nonfinite matrices, and extreme numerical conditioning remain outside the proper-rotation conversion domain, with validation/error policies unresolved under QD06/QD07. Forward-axis helpers and other quaternion operations are unchanged. Direct external mutation of matrix rows still does not refresh ctypes snapshots.

## Quaternion axes and arithmetic

`quat_from_axis_angle(axis,theta)` and `quat_rotate_from_axis_angle(axis,theta)` accept their existing Vector or Python-list axis representations. Ordinary inputs use a finite, nonzero three-component axis. Both helpers normalize a temporary axis, preserving the caller's Vector component-list identity, values, and list inputs. Angles remain degrees. The X/Y/Z angle-only helpers retain radians, and `quat_rotate` keeps its separate existing axis handling.

`quat_from_axis_angle` returns the rotation Quaternion `[cos(theta/2), n*sin(theta/2)]` for unit axis n. `quat_rotate_from_axis_angle` preserves its legacy result: the Quaternion sandwich rotating normalized n about itself, approximately `[0,n]`. This is not a resolution of Q11's proposed rotation-quaternion return contract. The repaired list branch uses a temporary Vector for that sandwich instead of adding Quaternion/list multiplication. Unsupported axis representations still return NotImplemented; zero-axis normalization still raises ZeroDivisionError. No broader shape, nonfinite, or extreme-scale policy is introduced.

Quaternion/Vector multiplication is the Hamilton product `q*(0,v)`, returning a Quaternion rather than a rotated Vector. In-place multiplication uses the same component kernel, replaces receiver storage, and returns the receiver. For `q=(w,u)`, the result is `(-u dot v, w*v + u cross v)`. Neither form normalizes operands or mutates the Vector. Rotation still requires the conjugate sandwich in the existing rotation helper.

Python 3 `/` and `/=` delegate to the retained `__div__`/`__idiv__` implementations. Quaternion wrappers accept floats, including float subclasses, and return NotImplemented for other operand types. Integer, boolean, complex, vector, and quaternion divisors are not added, and no reflected operators are introduced. Returning division creates a fresh Quaternion; in-place division returns the receiver and replaces its component list. A zero float divisor raises ZeroDivisionError without replacing receiver storage. The raw `quat_div_float` helper retains its existing Python arithmetic independently of wrapper operand acceptance.

## Quaternion powers and logarithms

`pow(e)` / `quat_pow(q,e)` and `log()` / `quat_log(q)` support unit
quaternions in [w,x,y,z] order. They do not normalize inputs; general
nonunit quaternions are outside the supported domain. No norm tolerance is
introduced. Zero inputs raise ValueError.

The principal quaternion angle is `atan2(|imaginary|,w)` in [0,pi]. Powers
use `[cos(e*angle), axis*sin(e*angle)]` for finite real exponents, returning
a fresh Quaternion. Zero and one powers return identity and a fresh exact
copy, respectively. Negative powers retain Hamilton inverse semantics.
Signs are not canonicalized: q and -q can have different fractional powers.

Logarithms return a fresh four-element list `[0,axis*angle]`, with identity
mapping to four zeros. Negative identity has no unique imaginary axis:
integer powers use parity, while fractional powers and logarithms raise
ValueError. Exactly zero imaginary parts are handled separately. Stable
hypotenuse calculations preserve tiny nonzero imaginary directions without
an arbitrary cutoff. Very large exponents remain subject to floating-point
angle-reduction accuracy; exponent reduction avoids intermediate overflow.

## Quaternion interpolation

`Quaternion.slerp(other,t)` / `quat_slerp(q0,q1,t)` perform accurate
shortest-path spherical interpolation of unit quaternions. Inputs are not
normalized and t is not clamped. For t in [0,1], unit norm and known-axis
components are checked within 1e-12, subject to ordinary floating-point
limitations. Difference/sum norms compute the spherical angle stably, even
when dot rounds to one. Equal orientations return a fresh copy.

Negative endpoint dot selects the sign-equivalent shortest path. At exactly
zero dot (a spatial half-turn), both paths are equally short; supplied signs
retain the existing tie branch. There is no global sign canonicalization.
Accuracy requirements do not extend to arbitrary extrapolation or nonunit
inputs. Quaternion LERP and slerp_no_invert retain their historical policies.

Legacy `Quaternion.squad(q1,q2,t)` / `quat_squad(q0,q1,q2,t)` use q0 as
start, q2 as end and q1 as an additional blend control:
`N(N(q0,q2,t),N(q0,q1,t),2*t*(1-t))`, where N is slerp_no_invert.
That helper uses spherical interpolation for -0.95 < dot < 0.95 and
unnormalized LERP otherwise. The legacy blend is sign-sensitive, may be
nonunit, and may degenerate to zero for antipodal inputs. Independently
negating controls can change the path; negating all controls preserves
represented rotations. No new antipodal-axis policy is introduced.

`gem.quaternion.squad4(q0,q1,s0,s1,t)` is separate conventional SQUAD.
q0/q1 are endpoint unit quaternions; s0/s1 are intermediate unit SQUAD
controls, not simply neighbouring animation keyframes. It computes
`S(S(q0,q1,t),S(s0,s1,t),2*t*(1-t))` using accurate shortest-path SLERP S.
It returns fresh Quaternion storage and preserves every input. Endpoints
match q0/q1 as rotations, with unit norm within 1e-12 over [0,1]. Sign
equivalence follows SLERP, including half-turn ties. Control generation and
an exponential helper are not part of this API.

### Examples

```python
from gem.quaternion import Quaternion, quat_from_axis_angle, squad4

start = Quaternion()
end = quat_from_axis_angle([0, 0, 1], 90)
legacy_control = quat_from_axis_angle([0, 0, 1], 180)
legacy = start.squad(legacy_control, end, 0.5)  # 67.5 degrees about Z

# Explicit SQUAD controls for this example, not neighbouring keyframes.
s0 = quat_from_axis_angle([0, 0, 1], 30)
s1 = quat_from_axis_angle([0, 0, 1], 120)
curve = squad4(start, end, s0, s1, 0.5)        # 60 degrees about Z
```
