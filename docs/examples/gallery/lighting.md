# From HDR radiance to diffuse lighting

![Original sphere](../../../examples/output/sh_original.png)
![Rotated sphere](../../../examples/output/sh_rotated.png)

These existing 192×192 PNGs come from the complete CPU-only
[HDR reference workflow](../../../examples/hdr_sh/README.md). An asymmetric
procedural linear HDR environment has a prominent feature along (0.8,0,0.6).
Pixel-center latitude-longitude quadrature projects it into canonical RGB L2 SH
coefficients. Analytical active +90° Z rotation moves the feature toward +Y:
f_rotated(d)=f_original(R⁻¹d). No directions are resampled or reintegrated for rotation.

Diffuse convolution is applied exactly once with factors π,2π/3,π/4 by band.
Core reconstruction evaluates irradiance at visible sphere normals. Reflected
Lambertian radiance is albedo×irradiance/π. Both images share camera, material,
exposure and Reinhard tone mapping followed by sRGB encoding. The right-side
brightness moves upward under rotation. Negative L2 truncation lobes are clipped
only for display; linear coefficient/reference data remains available.

Canonical ordering is l(l+1)+m, with Condon–Shortley signs. Historical positive-X/Y
probe coefficients require explicit `legacy_to_canonical` conversion. Angular
disk probes and latitude-longitude maps use different mappings; mirrored-ball
photographs are not an interchangeable supported format.

This sphere is reconstructed on the CPU. The accompanying GLSL formula and
coefficient exports illustrate shader integration but do not constitute a GPU
renderer. The gallery links the golden outputs without duplicating or changing them.

[Manifest](../../../examples/output/visualization.json) ·
[GLSL evaluation](../../../examples/hdr_sh/diffuse.glsl) ·
[Lighting tutorial](../../tutorials/lighting.md) ·
[SH API](../../api/spherical-harmonics.md) · [Gallery](../index.md)
