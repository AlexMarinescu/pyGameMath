# Spherical-harmonics sampling and irradiance

E06/E08 are corrected and validated functionality is promoted into the single
core `gem.spherical_harmonics` module. Experimental basis/sample/probe modules
are thin reexports. Legendre remains the existing canonical dependency.

## Corrections and compatibility

GenerateSamples fixes math.PI and the downstream `.vec` write, constructing
real Vector directions. Stratification, global RNG and constructor direction
ownership are retained; the progress print is removed. Seeded calls are
repeatable, rather than silently replacing random sampling with a fixed grid.

The raw probe loader uses row-major height/width indexing and supports angular
disks stretched to rectangular images. Pixel centers replace historical edge
sampling. The solid-angle Jacobian is preserved mathematically and uses both
image dimensions. Recalculation clears accumulated coefficients; direct
updateCoefficients continues accumulating. Native-endian float32 RGB is still
the file format. Dimensions and truncated data now produce clear ValueError
rather than unsafe indexing/unpacking failures.

The historical `.coeffs` basis and rounded update constants are preserved.
Explicit legacy_to_canonical conversion changes the odd-|m| signs and accounts
for rounded scales. New radiance projection uses canonical Condon–Shortley real
SH. Diffuse convolution produces a separate coefficient list and is never
implicitly applied by projection or reconstruction. First-three-band cosine
factors are pi, 2*pi/3, pi/4. No coefficient rotation, renderer, visibility or
shadow transport is included.

See [API, mapping, weights and examples](../docs/SPHERICAL_HARMONICS.md).
Existing packaging automatically includes the new core module and transition
shims; no mandatory dependency or unrelated experimental move is introduced.

## Reproduction and verification

Base: f03b0f525400a9e238a08f20853e53dec12eef51, merged PR #27.
Baseline full suite: **1642 passed, 7 xfailed**. Both E06/E08 probes fail with
--runxfail. Temporarily bypassing math.PI demonstrates zero direction lengths
and misplaced `.vec` storage. A square probe's DC coefficient doubles on
recalculation; +X legacy/canonical basis values have opposite signs.

Final full suite: **1687 passed, 0 failed, 5 xfailed**, Python 3.12.14,
pytest 9.1.1. Exact remaining identities match baseline minus E06/E08:
E07 transport, C02 viewport and three unresolved Vector contract cases.
The two corrected defect probes exercise core imports. See
`phase2f5-test-results.json` for identities and numerical measurements.

Independent Cartesian basis formulas verify all nine signs, indexing and
normalization. Equal-solid-angle midpoint references verify constant/asymmetric
RGB integrals, low-frequency reconstruction, orthogonality, explicit weights
and diffuse factors. Tests cover seeded repeatability, storage preservation,
2D rectangular probe orientation, individual pixel Jacobians, outside-disk
pixels, legacy rounded-scale conversion and malformed inputs. Mixed-scale RGB
values 1e200/1e-200 remain representable without cross-channel contamination.

An end-to-end constant angular probe (RGB=[1,2,3]) approaches irradiance
pi*RGB. Maximum error over +X/+Y/+Z is approximately .00635, .000928, .000300
and .0000437 for resolutions 16x8, 32x16, 64x32 and 128x64. These demonstrate
quadrature convergence, not a universal error bound. A one-pixel test explicitly
confirms weights are not renormalized to force a constant answer.

A wheel is built and actually pip-installed into an isolated target. Smoke
tests verify canonical/compatibility identities, sample generation, RGB
projection, convolution, reconstruction, rectangular loading and conversion.
Angular projection uses compensated coefficient accumulation with O(bands²)
auxiliary storage, rather than storing an image-sized sample object array.
The ordinary sample projection uses fsum; high-order/extreme numerical policy
is not redesigned.

Experimental retirement remains staged. E07's unfinished object/transport code
stays experimental with its existing expected failure. SH rotation is a
separate Phase 2F-5B review.
