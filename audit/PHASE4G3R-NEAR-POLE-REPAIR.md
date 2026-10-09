# Phase 4G-3R — SH near-pole numerical repair

The 4G3-A01 repair retains transverse directional information that was lost when
`cos(theta)` rounded to ±1. All ten audit failures now pass, together with 138 new
boundary cases. The runtime change is confined to `gem/spherical_harmonics.py`.

Base: `72fb9dac9beadc6c7ca1eae8a615b4ad86441254` (merged PR #54).
Branch: `fix/phase4g3r-near-pole-sh`.

## Mathematics and implementation

The former associated seed recovered sine as `sqrt((1-cos(theta))*(1+cos(theta)))`.
At theta = 1e-9, cosine rounds to 1 and that calculation produces zero even though
sine is accurately represented. The angle-aware private seed now uses `sin(theta)`
inside the documented polar interval. It feeds the existing Legendre recurrence;
the public Legendre API, implementation and mutable scratch-state semantics remain
unchanged. Exact angular endpoints 0 and `math.pi` retain zero transverse terms,
including the historical treatment of the small library sine residual at pi.
Outside the documented polar interval, the historical seed arithmetic is retained.

Unit-direction reconstruction uses `hypot(X,Y)` and Z directly, avoiding the
`acos(Z)` round trip. Reconstructing theta through `atan2(hypot(X,Y),Z)` and then
recomputing sine would still lose tiny south-pole components when theta rounded to
pi. The direction's first- and second-order azimuth factors likewise come from
X/hypot(X,Y), Y/hypot(X,Y) and double-angle identities. This prevents a second loss:
`atan2(Y,X)` can round onto a cardinal azimuth, whose trigonometric residual would
overwrite a much smaller supplied transverse component. These ratios represent
azimuth, not implicit normalization of the supplied direction.

The shared associated seed is necessarily used for all existing degrees. First-
and second-order azimuth identities are shared across degrees; higher-order phase
calculations, recurrence, normalization, factorial limits and supported band counts
are otherwise unchanged. There is no new threshold, clipping, direction-domain
validation or public API. Condon–Shortley signs, `l*(l+1)+m` indexing, RGB behavior,
independent result storage and transitional imports remain established contracts.

Ordinary-angle results can differ in their last bits. Finite unit directions remain
a prerequisite; nonunit directions are not normalized. Subnormal output rounding,
unrepresentable products and documented high-order/extreme-coefficient limitations
remain binary64 limitations rather than new accuracy guarantees.

## Independent verification

The original reproductions were executed before changes with markers disabled:
**10 failures, 181 deselected**, as expected. Unmodified master had **3,165 passes,
10 expected failures, zero unexpected failures and zero skips**.

The six basis and four reconstruction test bodies and tolerances are unchanged;
only their two strict defect markers and audit-module description changed. Initial
repaired SH/rotation/HDR/optimization/audit checks passed all 414 cases with markers
disabled before marker removal. The new module adds 138 cases using 180-digit
Decimal Cartesian L0/L1/L2 expressions and differentiated Fraction Rodrigues
polynomials. It covers:

- Adjacent binary64 angles across cosine's rounding-to-one boundary, both poles,
  signed orders, and angles/components through subnormal magnitudes.
- Signed and mixed-scale X/Y reconstruction, independent RGB channels, exact
  endpoints, fresh output, cached basis ownership and input storage preservation.
- The unchanged recurrence at degrees 3, 8, 12 and 32, existing numeric angle
  representations, and independent near-pole rotation/convolution/projection values.

Error budgets are specified in ulps for these new low-order/subnormal tests; the
original audit's relative tolerances are retained. Existing independent rotation,
convolution, quadrature, shader, projection and coefficient-indexing tests also pass.
Machine-readable counts, reproduction identities, wheel results and preservation
checks are in [phase4g3r-verification.json](phase4g3r-verification.json).

### Ordinary comparison with master

Deterministic matched comparisons use the immutable base source in the same
interpreter. They supplement independent references, rather than serve as the sole
correctness oracle.

| Workload | Components | Bit-identical | Maximum absolute difference |
| --- | ---: | ---: | ---: |
| Ordinary `SPH`, degrees 0–12 | 4,732 | 2,083 | 4.22e-15 |
| Ordinary cached basis, 1/3/5/9 bands | 3,248 | 1,566 | 3.33e-15 |
| Ordinary RGB reconstruction, 1/4/9/25 coefficients | 1,536 | 687 | 3.22e-15 |
| Exact angular poles, degrees 0–31 | 2,048 | 2,048 | 0 |
| Analytical rotation, convolution, legacy conversion | 27 each | 27 each | 0 |

Ten matched historical invalid/out-of-domain observations retain the same numerical
results or exception classes. These observations do not define a broader supported
input domain. Public signatures are unchanged. All other runtime modules, existing
examples/fixtures, dependencies and packaging metadata are byte-identical to base.

## HDR regeneration

Executed from the repository root using CPython 3.12.14:

```sh
/workspace/.venvs/pyGameMath/bin/python -m examples.hdr_sh.regenerate --output-dir /tmp/phase4g3r-hdr-output
```

All seven reference files were generated through the actual core. Both 192×192
PNG files retain identical bytes and decoded pixels: zero changed channel bytes.

| Image | SHA-256, committed and regenerated |
| --- | --- |
| `sh_original.png` | `dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a` |
| `sh_rotated.png` | `c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7` |

JSON and GLSL files are **not** byte-identical. Improved sine/azimuth arithmetic
changes last bits and tiny quadrature cancellation residuals. Shapes and nonnumeric
fields are unchanged; the largest recorded linear irradiance difference is 1.07e-14.

| File | Changed numeric fields/constants | Maximum absolute difference |
| --- | ---: | ---: |
| `reference.json` | 57 | 8.88e-16 |
| `visualization.json` | 66 | 1.07e-14 |
| `irradiance_coefficients.glsl` | 16 | 3.67e-17 |
| `sh_original_coefficients.glsl` | 9 | 4.03e-18 |
| `sh_rotated_coefficients.glsl` | 9 | 4.03e-18 |

The committed goldens were not replaced. The comparison data stores both sets of
hashes and every changed JSON value. The older performance-only
`benchmarks/verify_hdr_reference.py` requires all seven files to be byte-identical;
that criterion cannot describe this numerical repair. Its implementation is
unchanged. The repair harness instead reports byte differences, checks decoded PNG
pixels and compares JSON/GLSL numerically; existing independent pre-tone-mapping
and shader tests remain the correctness checks. Same-host PNG equality does not
promise byte equality across libm/zlib/platform versions.

## Performance

Two independent nine-trial runs compare the base and repaired modules in the same
CPython 3.12.14 interpreter, Linux x86_64, Intel Xeon Platinum 8370C. Measurements
use process CPU time, alternating paired implementation order, deterministic
shuffled workload order, warmup and GC disabled only during timing. Setup is
excluded; loop overhead and returned allocations are included. RNG is reset outside
timing for paired sample-generation trials and restored afterward. Each normal
workload executes 20,000 calls/trial, B9 basis 5,000, angular projection 39, and
sample generation 625. All individual trials, medians, MADs, environment/source
hashes and separate profiles are stored in
[run one](phase4g3r-performance.json) and [run two](phase4g3r-performance-repeat.json).

The table below is generated from those measurements. Latency ratios are medians
of paired after/before ratios, not ratios of separately computed medians. A ratio
above 1 means more time. Shared-host scheduling/frequency effects are uncontrolled.

| Workload | Run 1 before → after (µs) | Paired ratio | Run 2 before → after (µs) | Paired ratio |
| --- | ---: | ---: | ---: | ---: |
| `SPH_L0` | 0.794 → 0.794 | 1.016 | 0.741 → 0.777 | 1.028 |
| `SPH_L1_ordinary` | 1.097 → 1.212 | 1.094 | 1.091 → 1.300 | 1.191 |
| `SPH_L1_near_north` | 1.106 → 1.236 | 1.120 | 1.094 → 1.235 | 1.121 |
| `SPH_L2_ordinary` | 1.270 → 1.372 | 1.130 | 1.356 → 1.359 | 1.085 |
| `SPH_L2_near_south` | 1.305 → 1.397 | 1.126 | 1.346 → 1.432 | 1.091 |
| `SPH_L12_m6` | 3.052 → 3.143 | 1.033 | 3.098 → 3.204 | 1.052 |
| `basis_B3_ordinary` | 6.185 → 6.262 | 1.022 | 6.273 → 6.374 | 1.014 |
| `basis_B3_near` | 6.082 → 6.108 | 1.029 | 6.256 → 6.278 | 0.994 |
| `basis_B9_ordinary` | 69.899 → 66.731 | 0.961 | 68.136 → 65.276 | 0.944 |
| `reconstruct_L2_ordinary` | 14.093 → 13.900 | 0.998 | 13.668 → 14.012 | 1.025 |
| `reconstruct_L2_near_north` | 14.191 → 14.237 | 0.994 | 14.358 → 14.126 | 1.010 |
| `reconstruct_L2_near_south` | 14.911 → 14.548 | 0.948 | 14.703 → 14.427 | 0.963 |
| `reconstruct_B5` | 37.773 → 36.380 | 0.970 | 37.496 → 38.064 | 1.017 |
| `rotate_L2_RGB_control` | 27.537 → 27.730 | 1.038 | 28.043 → 26.623 | 0.983 |
| `convolve_L2_control` | 6.538 → 6.741 | 1.026 | 6.594 → 6.819 | 1.014 |
| `project_samples_16_control` | 70.648 → 69.196 | 0.973 | 73.360 → 74.384 | 0.996 |
| `project_angular_32x16` | 4619.759 → 4693.094 | 1.026 | 5262.251 → 4942.604 | 1.005 |
| `GenerateSamples_4x4_B3` | 198.445 → 211.792 | 1.042 | 206.081 → 231.554 | 1.094 |

Scalar nonzonal L1/L2 calls consistently take more time: paired increases span
8.5–19.1% across the two runs, roughly 0.1–0.2 µs per call. Extra private dispatch,
angle-domain handling and retained sine calculations are the cost of the repair.
L0 is unchanged. Cached L2 basis and ordinary/near-north reconstruction differences
are small or reverse between runs; they do not establish a speed improvement.
B9 basis and near-south reconstruction show lower paired times in both runs, but
variation and unchanged controls preclude a general speedup claim. Sample generation
shows a 4–9% paired increase, with substantial variability in the second run.

Individual ordinary/near L1/L2 before/after MADs range from 1.2–8.6% of their
respective medians. Unchanged rotation, convolution and precomputed projection
controls shift by roughly -3–4%, illustrating noise; the second angular-projection
run has a 14% baseline MAD. All samples are retained, including outliers.

Separate profiling of 1,000 cached L2 basis calls confirms that six Legendre objects
per call and recurrence counts are retained. Associated seeds now avoid three
`run`/square-root seed paths per call, sharing one transverse magnitude. Direct
reconstruction retains two azimuth ratios and first/second-order identities instead
of recovering them through cosine. The cache size and O(bands²) result storage are
unchanged. Profile instrumentation is not used as a latency oracle.

Reproduce the two measurement runs after HDR regeneration:

```sh
python benchmarks/sh_near_pole.py --output audit/phase4g3r-performance.json --hdr-directory /tmp/phase4g3r-hdr-output --trials 9 --iterations 20000
python benchmarks/sh_near_pole.py --output audit/phase4g3r-performance-repeat.json --hdr-directory /tmp/phase4g3r-hdr-output --trials 9 --iterations 20000
```

The harness requires the base commit in local git history. Fixed inputs include
angles .73/1.27, near-pole 1e-9/1e-100 components, unit XYZ [.3,.4,sqrt(.75)],
deterministic signed RGB coefficients, 16 precomputed sample bases, and a 32×16
angular image with radiance ranges R=[1,1.96875], G=[2,2.9375], B=.5.
The broader ordinary comparisons use seed 473101 and finite RGB values in [-2,2].

## Full suite and installed artifacts

Exact final full-suite command:

```sh
/workspace/.venvs/pyGameMath/bin/python -m pytest -q --junitxml=/tmp/phase4g3r-final.xml
```

Result: **3,313 passed, zero failed/errors/xfails/xpasses/skips**, 15.68 s
(pytest elapsed time). The count is 3,165 existing passes + 10 repaired regressions
+ 138 new cases. JUnit xunit2 output reports 3,313 warnings from the existing
`record_property` autouse fixture's format incompatibility; these are report-format
warnings, not runtime numerical warnings. Fixture and pytest configuration are
unchanged. No tests were deleted, weakened or converted to skips.

The existing installed-wheel tests pass in the full suite. A separate clean copy
of the package and unchanged build files was built outside the checkout:

```sh
cd /tmp/phase4g3r-package/source
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 setup.py sdist bdist_wheel
/workspace/.venvs/pyGameMath/bin/python -m venv /tmp/phase4g3r-wheel-env
/tmp/phase4g3r-wheel-env/bin/python -m pip install --no-index --no-deps /tmp/phase4d-dependencies/six-1.17.0-py2.py3-none-any.whl /tmp/phase4g3r-package/source/dist/gem-0.1.12-py3-none-any.whl
cd /tmp
/tmp/phase4g3r-wheel-env/bin/python -I -
```

The isolated stdin smoke program imports core/compatibility SH, asserts installation
under the venv rather than the source tree, checks four analytical angular L1 and
four signed/mixed-scale reconstruction cases, two exact poles and independent output
storage. It passes. Wheel and sdist Python module bytes match the current source.
The unchanged historical distribution version is not a release claim. Verification
records dependency versions; only CPython 3.12.14 was exercised here.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/spherical_harmonics.py` | Private angle-aware seed and stable direction basis |
| `tests/test_spherical_harmonics_audit.py` | Ten findings become ordinary regressions |
| `tests/test_sh_near_pole_repair.py` | 138 independent numerical/ownership boundary cases |
| `benchmarks/sh_near_pole.py` | Matched timing, rounding, profiles and HDR comparisons |
| `audit/phase4g3r-performance.json` | First complete paired run and HDR evidence |
| `audit/phase4g3r-performance-repeat.json` | Independent repeated measurement |
| `audit/phase4g3r-verification.json` | Test counts, original reproductions and preservation checks |
| `docs/api/spherical-harmonics.md` | Numerical behavior of the retained transverse basis |
| `audit/COMPATIBILITY.md` | Numerical differences and unchanged contracts |
| `audit/PHASE4G3R-NEAR-POLE-REPAIR.md` | Repair, compatibility and measured tradeoffs |
