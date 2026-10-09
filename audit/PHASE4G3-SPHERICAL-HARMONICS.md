# Independent spherical-harmonics review

Base: master `01ca2fc9f84e45e2c6978f3f5f9a173831533d97`, including merged PR #53.
Branch: `audit/phase4g3-spherical-harmonics`.

One numerical defect is confirmed: low-order SH loses nonzero transverse
components near the poles. Ten independent cases reproduce its basis and
reconstruction manifestations. Ordinary signs, normalization, projection,
convolution, analytical rotation, channel independence and ownership pass the
new checks. Production mathematics is unchanged. This is a bounded audit,
not a universal accuracy guarantee.

## Implementation and contract inventory

Review covered the actual source, [API inventory](../docs/architecture/api-inventory.md),
[API reference](../docs/api/spherical-harmonics.md), [SH guide](../docs/SPHERICAL_HARMONICS.md),
CONVENTIONS, COMPATIBILITY, PHASE2-DECISIONS, Phase 1/1B, 2F-4/5/5C/5D,
3D and 4G-2/2R evidence, existing SH/Legendre/optimization/HDR tests, and the
headless example adapters and shader formula. The historical wiki snapshot
`715e5039c75e080814a12e957f5148c35cdf8bda` contains no SH page. Frozen source
`5257291431bb45db0274dc48edf24694ecfe2e2d` already converts theta to cos(theta)
before associated-Legendre evaluation; this numerical weakness is inherited.

| Existing implementation | Contract and scope |
|---|---|
| `Factorial`, `K`, `SPH` | Nonnegative integer factorial; integer degree l≥0, order −l≤m≤l; unnormalized Condon–Shortley Legendre multiplied by orthonormal K. Polar theta from +Z, azimuth phi from +X toward +Y, radians. Historical scalar helpers have no uniform invalid-domain validator. |
| `gem.legendre.Legendre` | Canonical imported dependency, not a duplicate SH polynomial implementation. Associated 0≤m≤l, x∈[−1,1]. Its factored seed and recurrence operate on the supplied x; `run()` preserves scratch fields. |
| Private `_basis_layout`, `_basis` | Direction-independent normalization cache, at most 16 layouts; fresh basis arrays, shared +/-order polynomial evaluations. Private implementation details, not public cache APIs. |
| `SPHSample`, `GenerateSamples` | Mutable precomputed values, angles and direction. Container retains a supplied Vector reference; generation creates independent storage and jittered equal-solid-angle strata using global RNG. Default sample weight 4pi/N. |
| `project_radiance` | Finite RGB samples, complete matching finite basis arrays, optional nonnegative finite steradian weights. Uses stored values rather than regenerating from directions; fresh RGB coefficient rows, `fsum` accumulation. |
| `project_angular_probe` | Rectangular row-major linear RGB angular disk, pixel centers, area Jacobian, no weight renormalization. Compensated accumulation with O(bands²) auxiliary storage. Not lat-long or mirrored-ball projection. |
| `reconstruct` | Canonical complete-band RGB coefficients and caller-supplied unit Vector3/XYZ triple. Fresh RGB; no normalization, clamping or convolution. Checks finite/nonzero components and Z∈[−1,1], not full unit norm. |
| `convolve_diffuse` | Separate RGB irradiance coefficients for one, two or three bands. Degree factors pi, 2pi/3, pi/4. Lambertian reflected radiance additionally multiplies by albedo/pi. No general spherical-kernel convolution API. |
| `legacy_to_canonical` | Exactly nine legacy radiance rows; explicit odd-order sign conversion and exact/rounded scale ratios. Merely changing an import does not convert coefficient conventions. |
| `SPH_IrradianceMapCoeff` | Native-endian raw float32 angular probe, historical rounded positive-X/Y basis; `.coeffs` is radiance despite the class name. Load/recalculation rebuild; direct updates accumulate; output pretty-prints. |
| `rotate_coefficients` | Analytical active scalar/RGB rotation with 1/4/9 coefficients, finite Quaternion norm drift ≤1e−12 and temporary normalization. L1 linear forms and L2 traceless tensors; f_rotated(d)=f_original(R^-1 d). No directional reintegration, Matrix adapter or higher-band rotation. |
| `examples/hdr_sh` | Example-only lat-long adapter, limited RGBE decoder, linear radiometric pipeline, GLSL reference and CPU renderer. SH mathematics delegates to core. Inputs/formats are explicit; no GPU execution is claimed. |
| `gem.experimental.sph`, `sph_sample`, `sph_irradiance_map` | Thin compatibility imports of identical core objects. No independent algorithms. Shadow transport is retired and not part of this audit. |

Index is l(l+1)+m. Degree l occupies indices l² through (l+1)²−1;
complete B bands contain B² rows. Canonical L2 polynomials are
`[1,-Y,Z,-X,XY,-YZ,3Z²−1,-XZ,X²−Y²]` with their distinct orthonormal scales.
RGB storage is coefficient-major, and only rotation accepts scalar arrays.
Normalization means basis normalization, not normalization of a coefficient
array or an implicit direction/quaternion policy across the library.

**Order range:** rotation and diffuse convolution explicitly stop at degree 2
(three bands). Basis evaluation, sample generation, projection and reconstruction
accept general complete bands without a fixed public upper cap, but high-order
range/accuracy is explicitly not guaranteed. This review checks selected degrees
through 85, every order at nine selected degrees via the addition theorem,
orthogonality through degree 6, and projection/reconstruction through seven bands.
It does not establish a new supported-order ceiling. Degree 86 diagonal basis
already encounters factorial conversion overflow; some lower-order terms remain
representable at higher degrees. There is no single numerical ceiling shared by
all degree/order/direction combinations.

## Independent mathematical validation

[New regressions](../tests/test_spherical_harmonics_audit.py): **191 cases**, of
which 181 pass and ten are strict expected failures for 4G3-A01. Standard-library
references and the existing pytest dependency are used; no dependency is added.

| Area | Independent reference and evidence |
|---|---|
| Basis values and K | Explicit differentiated Rodrigues coefficients in Fraction; 180-digit Decimal associated factors and normalization before final float conversion. 280 signed-order comparisons across degrees 0,1,2,3,8,12,20,32,64,85 and five signed/equatorial arguments. |
| Low bands and coordinates | Closed-form Cartesian polynomials, including poles/equator, coordinate signs and all nine indices. These also supply the near-pole references without using cos(theta) to recover sine. |
| Orthonormality | Twelve Gauss–Legendre nodes solved from explicit Rodrigues polynomials by Decimal Newton iteration, times 32 Fourier azimuth nodes: 384 directions and 1,225 integrals over 49 terms through degree 6. No gem polynomial generates quadrature nodes. |
| Addition theorem/parity | Sum_m Y_lm²=(2l+1)/(4pi), antipodal parity (−1)^l and pole values at selected degrees through 85. |
| Projection | Exact-degree Gauss/Fourier quadrature of independent RGB signals `[2+.3X+.2XY, 3−.7Y+.8Z², 4+.4Z+.5(X²−Y²)]`; analytic coefficient moments, including absent higher bands. Stored sample values are independently calculated. |
| Reconstruction/diffuse | Closed-form linear/quadratic signals, plus direct independent hemisphere integrals of L(d)(n·d). First-three-band factors and albedo/pi are separately checked. |
| Directional lighting | A weighted directional sample and addition-theorem cosine-kernel polynomial `1/4+u/2+5(3u²−1)/32`; no nested implementation call is the oracle. |
| Rotation | Inverse Rodrigues directions and independent Cartesian SH for identity, X/Y/Z quarter/half turns, arbitrary axis, q/−q equivalence, inverse and noncommuting compositions. Per-band norms, scalar/RGB isolation and fresh storage. |
| Numerical/ownership boundaries | Scales 1e−280…1e280 in rotation, 1e−200/1e200 in reconstruction/projection, compensated cancellation, zero energy/weights, repeated calls, aliased input rows, malformed layouts/nonfinite data and actual quaternion norm-tolerance boundaries. |
| Sampling/probes | Seed 473002 restored after use; independent values and strata. Rectangular angular pixels, disk exclusion, solid-angle weights, raw loading, rounded legacy conversion and rebuild ownership. |
| HDR workflow | Independent continuous zonal moments for the procedural sharp feature, both original and actively rotated lighting, selected pre-tone-map pixels and Lambertian normalization. Existing shader/display tests remain intact. |

Other deterministic inputs use seed 473001+coefficient count, coefficients in
[−2,2], fixed axes/angles and explicitly separated scales. No randomized numerical
integration or uncontrolled RNG state is introduced. No gem matrix conversion,
quaternion rotation helper or SH basis is the sole oracle for another algorithm.

Measured reference errors, excluding the confirmed near-pole cases:

| Check | Maximum observed error |
|---|---:|
| K normalization, relative | 2.16e−16 |
| Rodrigues SH comparison, absolute | 1.27e−14 |
| Rodrigues SH comparison, relative where reference magnitude >1e−12 | 1.51e−13 |
| Orthonormality integral, absolute | 1.12e−15 |
| Addition theorem, relative | 1.72e−14 |

These are observations over the stated datasets. Tests allow ordinary rounding,
trigonometric and polynomial-conditioning error; they do not impose impossible
exactness at cancellation zeros or establish a universal high-order bound.

## Confirmed finding, ranked for repair

### 4G3-A01 — Near-pole SH loses nonzero transverse components

**Severity: low (P3), numerical correctness.** A valid, well-conditioned low-order
basis value can become exactly zero near either pole. Absolute errors are small
for ordinary coefficient magnitudes, so this is not evidence of widespread HDR
lighting failure. Relative error in the affected transverse basis contribution
is 100%, and finite coefficient scaling can make the lost contribution visible.

Locations in unchanged source:

- [spherical_harmonics.py:43](../gem/spherical_harmonics.py#L43): public SPH converts theta to cos(theta), including its nonzero-order branches.
- [spherical_harmonics.py:132](../gem/spherical_harmonics.py#L132): private basis repeats that conversion before associated-Legendre evaluation.
- [spherical_harmonics.py:243](../gem/spherical_harmonics.py#L243): reconstruction uses acos(Z), discarding the available transverse magnitude in X/Y.
- [legendre.py:27](../gem/legendre.py#L27): the associated seed correctly reflects its already-rounded x, but cannot recover information lost by SH coordinate conversion.

Minimal reproductions:

```python
from gem.spherical_harmonics import SPH, reconstruct
print(SPH(1, 1, 1e-9, 0.0))
c = [[0.0, 0.0, 0.0] for _ in range(4)]
c[3] = [1.0, 0.0, 0.0]
print(reconstruct(c, [1e-9, 0.0, 1.0]))
```

| Case | Actual | Independent expected |
|---|---:|---:|
| Y_1,1(theta=1e−9,phi=0) | −0.0 | −4.8860251190292e−10 |
| Y_1,1(theta=1e−100,phi=0) | −0.0 | −4.8860251190292e−101 |
| Y_1,1(theta=pi−1e−9,phi=0) | −0.0 | −4.886026121666232e−10 |
| Y_2,2(theta=1e−9,phi=.3) | 0.0 | About 4.51e−19 |
| Reconstruct c3=[1,0,0] at [1e−9,0,1] | [0,0,0] | [−4.8860251190292e−10,0,0] |

Derivation: Y_1,1=−sqrt(3/(4pi))·sin(theta)·cos(phi), whereas the implementation
recovers the imaginary factor from sqrt((1−x)(1+x)), x=cos(theta). For theta=1e−9,
x rounds to 1 in binary64, although sin(theta) and the answer are ordinary
representable nonzero numbers. Relative condition with respect to theta tends
to one. This is avoidable loss from the coordinate conversion, not unavoidable
underflow or severe polynomial conditioning. At theta=1e−7, the same mechanism
already causes about 0.04% relative error before the answer becomes zero.

Reconstruction has independent X/Y data available but uses only Z to recover
theta. `[1e−9,0,1]` is the binary64-rounded representation of a unit direction;
its hypot is exactly one and the exact-length discrepancy is below binary64
precision. Treating such rounded directions as invalid would reject normal unit
vector representations rather than preserve their geometry. The tests do not
extend reconstruction to arbitrary nonunit inputs.

**Contract:** the public API specifies real SH in polar angles and canonical
unit-direction reconstruction. The Cartesian polynomials in that reference
establish the expected values. This low-degree angular loss is separate from the
explicitly documented high-order and extreme-coefficient arithmetic limits.

**Coverage gap:** previous independent Cartesian tests used poles, equator and
moderate transverse components; higher-order tests compared against the same
rounded cos(theta) argument. Neither retained small nonzero transverse geometry.
Six new basis cases cover positive/negative orders, both poles, L1/L2 and a very
small yet normal answer. Four reconstruction cases cover X/Y, north/south and
independent RGB storage. With xfail disabled, all ten fail the numerical assertion.

**Recommended repair:** retain angle-aware sine magnitude when seeding associated
functions inside SH; preserve the existing public x-based Legendre contract.
For reconstruction, retain hypot(X,Y) and use a stable polar conversion, together
with a basis path that actually preserves that transverse factor. Merely replacing
acos(Z) with atan2(hypot(X,Y),Z) and then taking cos(theta) again is insufficient.
Check both poles, signs, exact endpoints, higher orders, ordinary rounding and
near-subnormal products. No epsilon, arbitrary transverse axis, whole-direction
normalization or new SH API is needed.

**Compatibility:** repaired tiny numerical contributions would replace erroneous
zeros while signatures, basis/order, input ownership and return types remain.
Historical compatibility imports would share the corrected implementation.
Keep the pole/angle domain and existing invalid-input policies explicit rather
than broadening them during repair.

**Performance:** an angle-aware seed/basis path may add work to SPH, cached basis
arrays and reconstruction, or may avoid costly inverse trig for L2. Compare
ordinary and near-pole paths with the existing benchmark tooling before selecting
an implementation. No performance gain or repair implementation is claimed here.

## Graphics reference and approximation checks

For the procedural feature L(d)=ambient+peak·max(axis·d,0)^32, axis=[.8,0,.6],
Funk–Hecke yields c_lm=2pi·I_l·Y_lm(axis), with
I0=1/33, I1=1/34, I2=(3/35−1/33)/2. These exact moments independently validate
the original and +90° Z-rotated axes, canonical signs and RGB channels.
At the committed 128×64 resolution, maximum radiance-coefficient discrepancy
from this continuous reference is about 4.13e−4. Selected pre-tone-map irradiance
pixels differ by at most 8.00e−4. The tests budget integration error separately
from roundoff. The GLSL remains a reference formula, not a tested GPU program.

Constant RGB=[1,2,4] maps approach irradiance pi·RGB. Angular-disk quadrature
and lat-long center radiance with exact pixel-area weights have different errors:

| Image size | Angular max irradiance error | Lat-long max irradiance error |
|---|---:|---:|
| 16×8 | .00846836 | .10496828 |
| 32×16 | .00123723 | .02547838 |
| 64×32 | .00040018 | .00632348 |
| 128×64 | .00005823 | .00157801 |

Errors are over +X/+Y/+Z and all three channels. Exact lat-long area weights
make the constant DC coefficient exact here, but midpoint sampling of the other
basis functions still produces L2 error. Angular weights are intentionally not
renormalized. New tests compare each finite quadrature with independently summed
Cartesian moments; an overly tight continuous-sphere tolerance would incorrectly
label coarse quadrature a defect. The table is convergence evidence for these
images, not a universal error guarantee.

All seven HDR files regenerate without overwriting references. Coefficients,
selected linear pixel values, decoded RGB and file bytes match. The two PNG hashes
remain:

| File | SHA-256 |
|---|---|
| sh_original.png | dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a |
| sh_rotated.png | c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7 |

Byte identity supplements the independent mathematical references; it is not
an accuracy oracle across platforms. Python/libm and zlib differences may alter
floating or compressed bytes on other systems.

## Other observations and classification

These are not additional strict expected failures or expanded contracts.

| Classification | Reproduction/evidence | Disposition |
|---|---|---|
| Documented high-order limitation | `SPH(86,86,pi/2,0)` raises OverflowError although normalized Y is finite; `SPH(170,0,.73)` returns infinity although its normalized answer is finite. Factorials and unnormalized functions are separate intermediates. | Already excluded by API high-order stabilization limits. Future scaled/normalized recurrence review; no band expansion or new guaranteed order range here. |
| Documented extreme rotation limit | Identity L2 rotation with only c6=1e308 returns infinity; finite coefficient 1e308 is the independent answer. Tensor extraction forms an overflowing intermediate numerator. | Phase 2F-5C explicitly excludes unrepresentable intermediate extreme arithmetic. Record range cost; do not manufacture an ordinary-domain defect or change the policy. |
| Documented extreme product limits | One angular center pixel at RGB=[1e307,0,0] overflows color×weight before multiplication by Y00 despite finite final c0≈1.11e308. A radiance product 1e−300×Y_1,-1(about −4.89e−101) underflows before weight 1e300 can recover it. | Supported finite input checks do not promise range-safe products in every extreme combination. Record exact reproductions in JSON; potential separately scoped scaling review. |
| Documented truncation | L2 directional cosine reconstruction is 17/16 at the light direction and 1/16 opposite it, rather than exact max(u,0). Negative reconstructed lighting can occur for other truncated signals. | Ordinary low-band approximation/ringing, not a sign, energy or convolution defect. No clamp added. |
| Documented quadrature | Coarse constant/angular/lat-long images have different residual L2 moments despite correct mapping/weights. | Preserve pixel-center quadrature and convergence qualification; not a universal finite-resolution error bound. |
| Unsupported behavior | Matrix orientation inputs, higher-band rotation/convolution, scalar reconstruction/projection, mirrored-ball interpretation and new visibility/probe systems. | No planned API assumed or added. |
| Compatibility questions | Historical Factorial/K/SPH/GenerateSamples do not uniformly reject invalid degrees/orders; reconstruction's unit norm is a prerequisite rather than a checked tolerance. Malformed foreign containers may raise native errors. | No generalized validation, epsilon or invalid-input exception policy inferred. Existing checked finite/layout errors pass. |
| Intentional ownership | SPHSample retains its Vector reference; direct legacy updates accumulate, rebuilding calculation replaces rows. | Preserve state contracts; new output-row alias and repeated-call tests pass. |
| False positives ruled out | Odd-order legacy sign changes require explicit conversion; isqrt(index) correctly identifies each complete degree band; q and −q rotations agree; noncommuting composition follows active inverse-direction order. | Independent low-band, directional, rotation and hemisphere references pass. No convention mismatch found in tested domains. |
| Unverified concerns | Arbitrarily high degree/order, large ill-conditioned cancellations, generalized finite-product scaling, other interpreters/platforms and actual GPU behavior. | No exhaustive guarantee or untested compatibility claim. |

## Verification and scope

Environment: CPython 3.12.14, Linux x86_64, pytest 9.1.1, six 1.17.0,
zlib 1.3.2. Only this interpreter/platform was executed.

| Run | Passed | Xfailed | Unexpected failures/errors | Skips |
|---|---:|---:|---:|---:|
| Unchanged merged-master baseline | 2,984 | 0 | 0 | 0 |
| New audit tests | 181 | 10 | 0 | 0 |
| Complete suite | 3,165 | 10 | 0 | 0 |

The reproduction-only `--runxfail -k near_pole` command fails exactly ten
assertions, deselects 181 cases and intentionally exits 1. All strict failures
are identified as 4G3-A01. There are no collection errors, unexpected passes,
retired tests, removed tests or weakened existing tests. Previous repaired
quaternion/curve regressions and existing isolated-wheel packaging/import tests
continue passing in the full suite.

[Machine-readable verification](phase4g3-test-results.json) records exact case
identities, actual/expected values, numerical metrics, environment, HDR equality
and 38 protected-file hashes. All tracked runtime/compatibility modules, HDR
sources/fixtures/outputs, dependency files and packaging metadata are byte-identical
to base. Public SH signatures and every existing test are unchanged.

Reproduce from the checkout using the existing audit environment:

```sh
python -m pytest tests/test_spherical_harmonics_audit.py -q --tb=short -o junit_family=legacy --junitxml=/tmp/phase4g3-focused.xml
python -m pytest tests/test_spherical_harmonics_audit.py --runxfail -k near_pole -q --tb=short -o junit_family=legacy --junitxml=/tmp/phase4g3-reproductions.xml
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase4g3-full.xml
python benchmarks/verify_hdr_reference.py --output /tmp/phase4g3-hdr-reference.json --output-dir /tmp/phase4g3-hdr-reference
git diff --exit-code 01ca2fc9f84e45e2c6978f3f5f9a173831533d97 -- gem examples setup.py setup.cfg MANIFEST.in requirements-audit.txt requirements-docs.txt requirements-docs-browser.txt
git diff --check
```

Commands used `/workspace/.venvs/pyGameMath/bin/python`; the system interpreter
does not have pytest installed. Only three files are added: this report, the
independent regression module and verification JSON. The gem 1.0 feature scope
remains frozen. No implementation repair, release or subsequent audit phase is
included.
