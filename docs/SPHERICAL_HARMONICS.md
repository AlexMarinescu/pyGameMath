# Spherical harmonics and diffuse environment lighting

Use `gem.spherical_harmonics` for validated basis evaluation, directional
sampling, RGB radiance projection and low-order diffuse convolution.
The experimental `sph`, `sph_sample` and `sph_irradiance_map` modules reexport
the same supported functions/classes. The unfinished `sph_object` transport module is retired with no core replacement.
These compatibility reexports remain pending a release-boundary removal decision;
see [import migration](EXPERIMENTAL_MIGRATION.md).

## Basis and coefficient layout

`SPH(l,m,theta,phi)` uses integer l >= 0 and -l <= m <= l. The associated
Legendre function includes Condon–Shortley phase. Theta is polar angle from
+Z and phi is azimuth from +X toward +Y, in radians. Normalization is
sqrt((2l+1)(l-|m|)!/(4*pi*(l+|m|)!)); nonzero orders also multiply sqrt(2),
with cosine for positive m and sine for negative m. `K(l,m)` takes
nonnegative m. `Factorial(n)` retains its historical nonnegative-integer API.

Coefficients use index `l*(l+1)+m`, with bands² entries. Each entry is an
RGB list. The first nine canonical basis functions are proportional to
`[1,-Y,Z,-X,XY,-YZ,3Z²-1,-XZ,X²-Y²]`, with orthonormal scale factors.
For example Y_1,1(+X)=-sqrt(3/(4*pi)). Analytical coefficient rotation through L2 is documented below.

## Directional samples and radiance

`SPHSample(theta,phi,dirc,sampleNumber)` retains a supplied Vector reference
and allocates sampleNumber basis slots. Its ownership behavior is unchanged.
`GenerateSamples(sqrtNumSamples,numBands)` returns sqrtNumSamples² samples,
with numBands² precomputed values per sample. Jittered strata use
u=(i+random())/sqrtNumSamples, phi=2*pi*v and cos(theta)=1-2*u. Thus each sample
represents equal solid angle 4*pi/N. The global RNG is unchanged: use a fixed
seed or reset its state for repeatability. Generation no longer prints status.
Positive integer grid/band counts are the supported generation domain;
broader legacy invalid-input behavior is not standardized here.

```python
import random
from gem.spherical_harmonics import GenerateSamples, project_radiance

random.seed(7)
samples = GenerateSamples(32, 3)
radiances = [[1.0, 2.0, 3.0] for sample in samples]
radiance_coeffs = project_radiance(samples, radiances)
```

`project_radiance(samples,radiances,weights=None)` integrates RGB radiance
against each sample's precomputed basis. Default weights assume uniform
sphere sampling. Explicit finite nonnegative weights are solid angles in
steradians; they are not renormalized. Basis arrays must contain matching
complete bands. Inputs and storage are preserved; output is a fresh RGB list.
The integration does not include cosine convolution.

## Angular probes

`project_angular_probe(hdr,numBands=3)` takes nonempty rectangular
`hdr[row][column][channel]` RGB data. It supports **angular disks only**,
not latitude-longitude maps or mirrored-ball photographs. A rectangular
probe stretches the disk independently along the two image axes.

At each pixel center:

- u=2*(column+0.5)/width-1, v=1-2*(row+0.5)/height.
- r=hypot(u,v); pixels with r>1 are ignored.
- theta=pi*r, phi=atan2(v,u).
- Direction=(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)).

Image right corresponds to +X, up to +Y; the center points toward +Z and
the disk perimeter toward -Z. Data outside the disk do not contribute.
The Jacobian is pi²*sinc(theta), because dOmega=sin(theta)dtheta dphi and
du dv=r dr dphi. Pixel solid angle is therefore
`4*pi²/(width*height)*sinc(theta)`, with sinc(0)=1.

This is midpoint quadrature, without weight renormalization. Finite-resolution
integrals approximate the continuum; tiny images can have large errors.
Increasing resolution improves smooth-probe accuracy, with no universal
finite-image error guarantee. RGB values are finite, independently integrated
and not clipped. No mandatory image-decoding dependency is introduced.

## Historical raw-probe class

`SPH_IrradianceMapCoeff(fileU,width,height)` preserves its name/signature,
reading width*height native-endian float32 RGB triples in row-major order.
It is not a Radiance HDR/image decoder. Positive integer dimensions and
sufficient binary data are required; trailing bytes remain ignored. Short
files and malformed dimensions raise ValueError; ordinary file I/O errors
propagate. Foreign-endian/image formats need explicit external conversion.

Despite the class name, `.coeffs` contains **radiance** coefficients in the
historical nine-polynomial basis, not irradiance coefficients. The legacy
basis has positive X/Y first-order terms and rounded normalization constants.
`updateCoefficients(hdr,domega,x,y,z)` retains those formulas and accumulates.
`calculateCoefficients()` and `load()` rebuild rather than double accumulated
coefficients; `.hdr` remains unchanged by calculation. `output()` still prints
coefficients.

`legacy_to_canonical(coefficients)` explicitly converts nine legacy RGB
entries into fresh canonical coefficients. Indices 1,3,5,7 change sign;
all entries also account for the historical rounded scale factors. Existing
`.coeffs` values are never silently reinterpreted.

```python
from gem.spherical_harmonics import (
    SPH_IrradianceMapCoeff, legacy_to_canonical,
    convolve_diffuse, reconstruct,
)

probe = SPH_IrradianceMapCoeff("probe.float", 512, 256)
radiance_coeffs = legacy_to_canonical(probe.coeffs)
irradiance_coeffs = convolve_diffuse(radiance_coeffs)
irradiance_rgb = reconstruct(irradiance_coeffs, [0.0, 0.0, 1.0])
```

## Reconstruction and diffuse convolution

`reconstruct(coefficients,direction)` returns RGB at a caller-supplied unit
Vector3 or XYZ triple. It never normalizes the direction or convolves the
coefficients. Unit length is a prerequisite, not a tolerance-based validation
policy. Empty/invalid coefficient layouts, nonfinite RGB, zero directions,
wrong component counts or Z outside [-1,1] raise ValueError.

`convolve_diffuse(radiance_coefficients)` supports one, two or three complete
bands. It returns fresh irradiance coefficients with degree factors pi,
2*pi/3 and pi/4. Reconstruction must use these directly, without another
convolution. Constant radiance L produces irradiance pi*L. Lambertian outgoing
radiance is albedo*irradiance/pi; rendering and albedo are outside this API.

The SH representation is a low-frequency approximation; all three RGB channels
remain independent. Extreme orders, ill-conditioned numerical integrations,
and unrepresentable binary64 results are outside this phase's accuracy contract.
No new NaN/Infinity policy is applied to historical scalar basis helpers.

## Analytical coefficient rotation

`rotate_coefficients(coefficients, orientation)` accepts complete canonical
L0, L0–L1 or L0–L2 arrays (1, 4 or 9 entries), as scalar values or RGB rows.
It returns independent output storage and never mixes bands or changes L0.
Use the same API for radiance and already cosine-convolved irradiance;
rotation commutes with the per-degree diffuse factors and never convolves.

Rotation is active and right-handed: `f_rotated(d)=f_original(R^-1 d)`.
Positive 90-degree rotation about Z moves a +X lighting feature toward +Y.
Orientation is a gem Quaternion in [w,x,y,z] order; q and -q are equivalent.
Applying a then b corresponds to Hamilton product b*a. The equivalent inverse
direction uses the transpose of gem's row-vector rotation matrix.

```python
import math
from gem.quaternion import Quaternion
from gem.spherical_harmonics import rotate_coefficients

orientation = Quaternion([math.sqrt(0.5), 0.0, 0.0, math.sqrt(0.5)])
# Canonical scalar L1: negative c3 corresponds to a positive-X feature.
coefficients = [0.0, 0.0, 0.0, -1.0]
rotated = rotate_coefficients(coefficients, orientation)
# Approximately [0, -1, 0, 0]: the feature now points toward +Y.
```

Only Quaternion orientations are supported. The finite quaternion norm may
deviate from 1 by at most 1e-12; a temporary copy is normalized to remove that
small drift. Zero, nonfinite or larger deviations raise ValueError, without
mutating the supplied object. Existing quaternion normalization rules are
unchanged. Matrix adapters and higher bands remain outside this API.

Historical probe arrays must first pass through `legacy_to_canonical`.
The rotation implementation uses L1 linear forms and L2 symmetric traceless
tensors; it does not rotate sample directions or reintegrate an environment.
See `audit/PHASE2F5C.md` for the derivation, verification and benchmark method.
