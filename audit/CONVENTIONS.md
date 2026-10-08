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
