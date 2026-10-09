# Sampled environment lighting with real spherical harmonics

**Advanced.** Prerequisites: [vectors](vectors.md), [rotation](quaternions.md),
[numerical accuracy](numerical.md) and the idea of integrating a function.
Learn to represent low-frequency RGB directional radiance, rotate its coefficients
analytically and evaluate diffuse irradiance without confusing the three radiometric
quantities. This is a CPU lighting reference, not a renderer or visibility model.

## Directional functions and real basis coordinates

An environment assigns linear radiance L(d) to a unit direction on the sphere.
Real orthonormal SH basis functions Y_lm give coefficients
c_lm=integral_sphere L(d)*Y_lm(d) dOmega. Reconstruction sums c_lm*Y_lm(d).
The coefficient index is l(l+1)+m, -l≤m≤l. One, two and three bands contain
1, 4 and 9 entries; L2 means maximum degree 2, thus three bands.
Higher bands describe finer angular variation but increase work and may not
capture a sharp light adequately at low order.

Canonical real SH uses the Condon–Shortley phase. L0 is 1/sqrt(4pi), and L1
indices 1/2/3 are sqrt(3/(4pi))*[-Y,Z,-X]. Theta is polar from +Z and phi azimuth
+X toward +Y, in radians. The negative X/Y signs matter for projection and rotation.
RGB arrays store [[R,G,B],...] per coefficient; scalar rotation arrays are a
separate supported representation. `reconstruct` expects RGB, not scalar arrays.

## A field with independent coefficients

Use L(d)=[1+X/2,2,3]. It is positive everywhere and asymmetric along +X.
Integrals of X and XY vanish, while integral X² dOmega=4pi/3. Thus
c0=sqrt(4pi)*[1,2,3], c3=[-0.5*sqrt(4pi/3),0,0], and other L1 entries are zero.
Six axis samples with equal weights 4pi/6 integrate these degree≤2 products
exactly in real arithmetic. This special quadrature is not an arbitrary high-band
environment integration rule; the hand-derived coefficients validate this example.

```python
import math
from gem.quaternion import quat_from_axis_angle
from gem.spherical_harmonics import SPHSample, SPH, project_radiance, rotate_coefficients, convolve_diffuse, reconstruct
from gem.vector import Vector

directions = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
samples = []
for xyz in directions:
    theta, phi = math.acos(xyz[2]), math.atan2(xyz[1], xyz[0])
    sample = SPHSample(theta, phi, Vector(3, xyz[:]), 4)
    sample.values = [SPH(l, m, theta, phi) for l in range(2) for m in range(-l, l + 1)]
    samples.append(sample)
colors = [[1 + xyz[0] * 0.5, 2, 3] for xyz in directions]
coefficients = project_radiance(samples, colors, [4 * math.pi / 6] * 6)
expected = [[math.sqrt(4 * math.pi) * x for x in [1, 2, 3]], [0, 0, 0], [0, 0, 0], [-0.5 * math.sqrt(4 * math.pi / 3), 0, 0]]
assert all(abs(a - b) < 1e-13 for row, ref in zip(coefficients, expected) for a, b in zip(row, ref))
rotated = rotate_coefficients(coefficients, quat_from_axis_angle([0, 0, 1], 90))
irradiance = convolve_diffuse(rotated)  # exactly once
on_y = reconstruct(irradiance, [0, 1, 0])
on_x = reconstruct(irradiance, [1, 0, 0])
assert all(abs(a - b) < 1e-13 for a, b in zip(on_y, [4 * math.pi / 3, 2 * math.pi, 3 * math.pi]))
assert all(abs(a - b) < 1e-13 for a, b in zip(on_x, [math.pi, 2 * math.pi, 3 * math.pi]))
albedo = [0.5, 0.5, 0.5]
reflected = [a * e / math.pi for a, e in zip(albedo, on_y)]
assert all(abs(a - b) < 1e-14 for a, b in zip(reflected, [2.0 / 3, 1, 1.5]))
assert coefficients is not rotated and rotated[0] is not irradiance[0]
assert colors[0] == [1.5, 2, 3] and samples[0].dir.vector == [1, 0, 0]
print([round(x, 6) for x in on_y], [round(x, 6) for x in reflected])
```

Output: `[4.18879, 6.283185, 9.424778] [0.666667, 1.0, 1.5]`.
The +90° Z rotation moves the red +X feature to +Y, while constant green/blue
remain unaffected. Active rotation satisfies L_rotated(d)=L_original(R^-1*d).
The implementation transforms L1 linear forms/L2 tensors analytically; it does
not rotate sample directions and reintegrate. Rotation accepts finite Quaternion
norm drift up to 1e-12 and normalizes only temporary storage.

## Radiance, irradiance and reflection

Clamped-cosine convolution has degree factors pi, 2pi/3, pi/4 through L2.
It creates separate irradiance coefficients; reconstruction does not convolve again.
Our red irradiance after rotation is pi+(pi/3)*Y. Reflected Lambertian radiance
is albedo*irradiance/pi, explaining the printed 2/3 rather than 4pi/3.
Applying cosine convolution twice or adding another 1/pi to irradiance itself
would mix these quantities. There is no scene visibility or occlusion term here.

Constant radiance is an especially useful reference: irradiance must be pi times
radiance, independent of normal. Sampling cannot justify that result more strongly
than an exact L0 coefficient check:

```python
import math
from gem.spherical_harmonics import convolve_diffuse, reconstruct, rotate_coefficients
from gem.quaternion import quat_from_axis_angle

radiance = [[math.sqrt(4 * math.pi), 0.0, 0.0]]
e = convolve_diffuse(rotate_coefficients(radiance, quat_from_axis_angle([1, 0, 0], 90)))
for normal in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
    rgb = reconstruct(e, normal)
    assert abs(rgb[0] - math.pi) < 1e-14 and rgb[1:] == [0.0, 0.0]
assert radiance[0] == [math.sqrt(4 * math.pi), 0.0, 0.0]
```

## From a real environment to the reference pipeline

The existing [HDR-to-SH workflow](../../examples/hdr_sh/README.md) generates a
linear asymmetric environment, projects canonical L2 RGB coefficients, rotates,
convolves once, reconstructs on a sphere and exports GLSL coefficients/reference
expressions. Run from the checkout with gem installed:

```sh
python -m examples.hdr_sh.regenerate --output-dir /tmp/gem-tutorial-hdr
```

This generates CPU PNGs, a manifest and shader coefficient files; it does not
exercise a GPU shader. Keep computations linear until final tone mapping/display
encoding. The generated images supplement independent linear pixel checks.
Latitude-longitude/RGBE adapters are example-only and have distinct solid-angle
mapping from angular-disk probes.

Core angular probes use row-major RGB pixel centers, right/up/center → +X/+Y/+Z,
theta=pi*r and Jacobian weighting, ignoring centers outside the disk. They are
not mirrored-ball photographs. Historical SPH_IrradianceMapCoeff stores rounded
legacy **radiance** coefficients: call legacy_to_canonical before using this
pipeline. Updating an import path alone does not convert their basis.

Finite sample integration, SH truncation and floating arithmetic are separate
error sources. Sharp lights may ring or reconstruct negative values; no silent
clamp or universal approximation bound is promised. GenerateSamples uses global
RNG jittered strata: seed externally for repeatability. Uniform weights assume
uniform-sphere sampling; arbitrary samples require justified solid-angle weights.
Reconstruction requires supplied unit directions and does not normalize them.

See [SH API](../api/spherical-harmonics.md), [SH guide](../SPHERICAL_HARMONICS.md),
[analytical rotation tests](../../tests/test_sh_rotation.py), [accuracy](numerical.md),
[quaternions](quaternions.md) and [tutorial index](index.md).
