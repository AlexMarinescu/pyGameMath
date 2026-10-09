# Real spherical harmonics and radiometry

Import from `gem.spherical_harmonics`. `Legendre` is a retained imported class,
not a second implementation; its canonical reference is [gem.legendre](legendre.md).
[Source](../../gem/spherical_harmonics.py), [projection/basis tests](../../tests/test_spherical_harmonics.py),
[analytical rotation](../../tests/test_sh_rotation.py), [HDR reference](../../tests/test_hdr_sh_example.py).

## Basis and coefficient layout

Degree l≥0 and order −l≤m≤l are integers; theta is polar angle from +Z and
phi azimuth from +X toward +Y, both radians. P_l^m includes Condon–Shortley phase.
K(l,m)=sqrt((2l+1)(l−m)!/(4pi(l+m)!)). Real orthonormal basis:
Y_l0=K(l,0)P_l^0(cos theta);
Y_lm=sqrt(2)K(l,m)cos(m phi)P_l^m(cos theta) for m>0;
Y_lm=sqrt(2)K(l,−m)sin(−m phi)P_l^(−m)(cos theta) for m<0.

Index=l(l+1)+m; complete numBands bands contain numBands² coefficients.
Canonical L2 coordinate polynomials and scales are:

| Index | (l,m) | Basis at unit (X,Y,Z) |
|---|---|---|
| 0 | (0,0) | 1/sqrt(4pi) |
| 1 | (1,−1) | −sqrt(3/(4pi))*Y |
| 2 | (1,0) | sqrt(3/(4pi))*Z |
| 3 | (1,1) | −sqrt(3/(4pi))*X |
| 4 | (2,−2) | sqrt(15/(4pi))*XY |
| 5 | (2,−1) | −sqrt(15/(4pi))*YZ |
| 6 | (2,0) | sqrt(5/(16pi))*(3Z²−1) |
| 7 | (2,1) | −sqrt(15/(4pi))*XZ |
| 8 | (2,2) | sqrt(15/(16pi))*(X²−Y²) |

RGB APIs use coefficient-by-channel lists [[R,G,B],...]. Only analytical rotation
also accepts scalar coefficient arrays. Reconstruction is RGB-only. No public
module constants or cache-management APIs are exposed; imported math/random/etc
and private basis-layout cache are implementation details.

## Functions and containers

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `Factorial(n)` | `n`: nonnegative integer prerequisite; scalar n!, iterative. Historical n≤1 returns 1 even outside domain; no broad validation. |
| `K(l, m)` | Integers `l,m`, 0≤m≤l; positive scalar normalization K above. No high-order overflow stabilization or uniform invalid-domain errors. |
| `SPH(l, m, theta, phi)` | Integer `l,m`, −l≤m≤l, finite numeric `theta,phi` radians; scalar real Y_lm. Uses unnormalized Legendre and K, preserving phase. No generalized angular/domain validator. |
| `SPHSample` | Mutable sample container; precomputed values drive projection, not re-evaluation of dir. |
| `SPHSample.__init__(self, theta, phi, dirc, sampleNumber)` | `theta,phi`: stored radians; `dirc`: Vector reference retained, otherwise fresh zero Vector3; `sampleNumber`: basis-array length (normally bands²), initialize values to zeros; return None. |
| `GenerateSamples(sqrtNumSamples, numBands)` | Positive integer `sqrtNumSamples` and `numBands` in supported use; fresh list of N=sqrtNumSamples² SPHSamples with complete basis values. Global-RNG jittered equal-solid-angle strata, default weight 4pi/N; no progress printing. Invalid arguments retain legacy arithmetic (zero causes ZeroDivisionError). |
| `project_radiance(samples, radiances, weights=None)` | Nonempty matching `samples`: SPHSamples and `radiances`: finite RGB triples. `weights=None`: 4pi/N, suitable for uniform sphere samples; otherwise matching finite nonnegative solid angles. Return fresh canonical RGB coefficients via weighted basis sums. Complete matching finite values required (ValueError); inputs preserved. |
| `project_angular_probe(hdr, numBands=3)` | `hdr`: nonempty rectangular row/column/finite-RGB lists; `numBands=3`: positive integer, bool rejected. Fresh canonical RGB radiance coefficients from angular-disk pixel-center quadrature; malformed layouts/nonfinite RGB raise ValueError on supported sequence forms. No image decoding. |
| `reconstruct(coefficients, direction)` | `coefficients`: finite RGB complete nonempty bands; `direction`: unit Vector3 or XYZ triple prerequisite. Fresh RGB sum c_i*Y_i; no normalization, clamping or convolution. Checks finite/nonzero triple and Z∈[−1,1], not full unit norm. Invalid checked data raises ValueError; malformed foreign containers may raise native errors. |
| `convolve_diffuse(radiance_coefficients)` | `radiance_coefficients`: finite RGB complete bands, 1/4/9 entries only; fresh RGB irradiance coefficients, degree factors pi,2pi/3,pi/4. Invalid checked layouts/nonfinite data or higher bands raise ValueError; no mutation. |
| `legacy_to_canonical(coefficients)` | Exactly nine finite legacy RGB radiance rows; fresh canonical rows, flipping odd-order signs and correcting rounded scale constants. Invalid checked layout/nonfinite data ValueError. Import migration alone does not convert a basis. |
| `SPH_IrradianceMapCoeff` | Historical file-backed angular-probe class; .coeffs is legacy radiance despite class name, not convolved irradiance. |
| `SPH_IrradianceMapCoeff.__init__(self, fileU, width, height)` | `fileU`: binary file path; `width,height`: positive ints excluding bool. Store file/dimensions, read native-endian float32 RGB, build hdr and nine legacy radiance rows; returns None. Invalid dimensions/short data ValueError; file errors propagate. |
| `SPH_IrradianceMapCoeff.load(self)` | Read exactly width*height*3 floats from file (trailing bytes ignored); replace row-major hdr and recompute legacy coeffs; return None. Missing file/OSError and short-data ValueError propagate. |
| `SPH_IrradianceMapCoeff.calculateCoefficients(self)` | Clear coeffs then accumulate hdr pixel-center quadrature, return None; repeated calls rebuild rather than double accumulation. Uses angular validation; .hdr preserved. |
| `SPH_IrradianceMapCoeff.updateCoefficients(self, hdr, domega, x, y, z)` | `hdr`: one RGB triple, `domega`: weight, `x,y,z`: supplied direction coordinates. Accumulate nine rounded legacy polynomial radiance terms into existing coeffs, return None. No normalization/convolution or new validation; input triple preserved. |
| `SPH_IrradianceMapCoeff.output(self)` | Pretty-print current coeffs to stdout; return None, no coefficient mutation. |
| `rotate_coefficients(coefficients, orientation)` | Canonical finite scalar arrays or RGB rows of length 1/4/9; `orientation`: finite Quaternion with norm drift ≤1e−12. Normalize temporary data, return independent active analytic rotation. Reject invalid orientation/layout with ValueError; no Matrix adapter, resampling or implicit legacy conversion. |

## Sampling, projection and reconstruction

`SPHSample.theta`, `.phi`, `.dir` and `.values` are public mutable state. Its
constructor retains a supplied Vector; GenerateSamples creates separate Vectors.
Set the global random seed for repeatability (this changes application RNG state):
u=(i+random())/sqrtNumSamples, v=(j+random())/sqrtNumSamples,
theta=2acos(sqrt(1−u)), phi=2pi*v. Projection uses stored values and caller RGB
samples, so keeping dir/angles/values consistent is the caller's responsibility.
Explicit weights are steradians and are not forced to sum to 4pi.

```python
import math
import random
from gem.spherical_harmonics import GenerateSamples, project_radiance, convolve_diffuse, reconstruct

random.seed(17)
samples = GenerateSamples(4, 1)
colors = [[1.0, 2.0, 3.0] for _ in samples]
radiance = project_radiance(samples, colors)
assert len(radiance) == 1
assert all(abs(a - math.sqrt(4 * math.pi) * b) < 1e-14
           for a, b in zip(radiance[0], [1, 2, 3]))
irradiance = convolve_diffuse(radiance)
value = reconstruct(irradiance, [0, 0, 1])
assert all(abs(a - math.pi * b) < 1e-13 for a, b in zip(value, [1, 2, 3]))
assert radiance[0] is not irradiance[0]
assert colors[0] == [1.0, 2.0, 3.0]
```

## Angular probes and the legacy coefficient class

Pixels use `hdr[row][column][channel]`, native float32 RGB for the historical file
class. For width w/height h, centers map to
u=2(col+.5)/w−1, v=1−2(row+.5)/h, r=hypot(u,v);
ignore r>1; theta=pi*r, phi=atan2(v,u). Right/up/center map to +X/+Y/+Z.
Rectangular images stretch the disk. Weight=4pi²/(w*h)*sin(theta)/theta, with
center limit 1; weights are not renormalized. Finite quadrature approximates
integrals and converges with resolution; no universal finite-resolution error
bound exists. This is not mirrored-ball photographic mapping or lat-long.

Legacy `.file`, `.width`, `.height`, `.hdr`, `.coeffs` are mutable. The nine
legacy radiance basis constants are 0.282095 (L0), +0.488603*(Y,Z,X),
1.092548*(XY,YZ,XZ), 0.315392*(3Z²−1), 0.546274*(X²−Y²), in indexed order
[1,Y,Z,X,XY,YZ,3Z²−1,XZ,X²−Y²]. `legacy_to_canonical` scales by the ratio of
exact to rounded constants and flips indices 1,3,5,7. Never pass these raw legacy
coeffs directly to canonical rotation/reconstruction. File errors are not hidden.

The following small probe deliberately samples just the disk center: its weight
is 4pi², not a renormalized 4pi sphere integral. This independently checks the
pixel-center mapping and illustrates coarse quadrature error, not an accurate
constant-environment approximation.

```python
import math
from gem.spherical_harmonics import SPH, project_angular_probe

assert abs(SPH(1, 1, math.pi / 2, 0) + math.sqrt(3 / (4 * math.pi))) < 1e-14
hdr = [[[1.0, 0.0, 0.0]]]
c = project_angular_probe(hdr, 1)
expected = 4 * math.pi * math.pi / math.sqrt(4 * math.pi)
assert abs(c[0][0] - expected) < 1e-14
assert c[0][1:] == [0.0, 0.0] and hdr == [[[1.0, 0.0, 0.0]]]
```

## Analytical active rotation

f_rotated(d)=f_original(R^-1 d): a +90° rotation about +Z moves a +X feature to +Y.
L0 stays fixed, L1 is a Cartesian linear-form rotation, and L2 a symmetric traceless
tensor transform T′=R*T*R^T. No directional sampling or reintegration is used.
Scalar/RGB channels transform independently. Inputs/storage are preserved;
orientation normalization is temporary and specific to this API, not a new
Quaternion-wide policy.

```python
from gem.quaternion import quat_from_axis_angle
from gem.spherical_harmonics import rotate_coefficients

# Canonical negative-X basis means coefficient -1 at index 3 points toward +X.
original = [0.0, 0.0, 0.0, -1.0]
q = quat_from_axis_angle([0, 0, 1], 90)
saved = q.data[:]
rotated = rotate_coefficients(original, q)
assert all(abs(a - b) < 1e-14 for a, b in zip(rotated, [0, -1, 0, 0]))
assert original == [0.0, 0.0, 0.0, -1.0] and q.data == saved
assert rotated is not original
```

## Radiance, irradiance and limits

Convolve radiance once with the clamped-cosine kernel to get irradiance. Reconstruct
irradiance directly; Lambertian reflected radiance is albedo*irradiance/pi.
All coefficients and calculations stay linear until final display conversion.
Truncated SH can ring and give negative values; APIs do not silently clamp them.
Basis order, quadrature error and ill-conditioned/extreme-value arithmetic are
separate limits. Historical scalar helpers do not validate all invalid domains;
newer checked finite/layout requirements are not a blanket guarantee for every
foreign object/type or overflowing sum.

[HDR/GLSL workflow](../../examples/hdr_sh/README.md) uses example-only lat-long/RGBE
adapters and a CPU sphere renderer. Its GLSL is a reference formula, not evidence
of GPU execution. See [SH guide](../SPHERICAL_HARMONICS.md), [legacy](legacy.md),
[decisions](decisions.md) and [index](index.md).

See the [graphics gallery example](../examples/gallery/lighting.md) for an executable visualization.
