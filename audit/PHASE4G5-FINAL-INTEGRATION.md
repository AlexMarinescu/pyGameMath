# Phase 4G-5 — Final cross-module integration audit

**Release signoff is withheld.** One new supported-contract defect is confirmed:
SLERP loses a representable endpoint rotation at the smallest subnormal
separation, and `squad4` inherits it. Its absolute effect is extremely small,
but it violates endpoint interpolation and remains a correctness gate.
No implementation is repaired in this audit.

Base: `fe453893b2c3a68b6619206cc3214036bf8c24c1`, merged PR #57.
Branch: `audit/phase4g5-final-integration`. The unchanged baseline is exactly
**3,714 passing tests**. The final suite has **3,853 passes, 12 strict expected
failures, zero unexpected failures/errors and zero skips**. All 12 expected
failures describe the single new defect; every original test is unchanged.

[Machine-readable evidence](phase4g5-verification.json) contains case identities,
counts, independent actual/expected endpoint values, distribution hashes,
installed checks, documentation inventory, two performance runs, HDR numeric
differences and preservation fingerprints. This is bounded deterministic
validation, not a proof of every possible input or platform.

## Inventory and established contracts

The review covers the focused [core audit](PHASE4G1-CORE-REVIEW.md) and
[repair](PHASE4G1R-QUATERNION-REPAIRS.md), [curves/Legendre audit](PHASE4G2-CURVES-LEGENDRE.md)
and [repair](PHASE4G2R-NUMERICAL-REPAIRS.md), [SH audit](PHASE4G3-SPHERICAL-HARMONICS.md)
and [repair](PHASE4G3R-NEAR-POLE-REPAIR.md), [geometry audit](PHASE4G4-GEOMETRY.md)
and [repair](PHASE4G4R-GEOMETRY-REPAIR.md), the Phase 3 performance reports,
current source/tests, README, API pages and
[conventions](CONVENTIONS.md), [compatibility](COMPATIBILITY.md) and
[decisions](PHASE2-DECISIONS.md).

AST inventory and installed runtime signatures agree with all **268 documented
declarations**: 112 functions, nine classes, 143 methods including dunders, and
four division aliases. Seven exposed constants/buffers and compatibility
reexports are checked separately. The full declaration list is in the JSON.

| Canonical module | Declarations | Implemented scope |
| --- | ---: | --- |
| `gem.vector` | 76 | Dimensioned Vector, arithmetic, stable norms, transform, reflection/refraction, barycentric and direction helpers |
| `gem.matrix` | 57 | Dimensioned Matrix, multiplication, 2/3/4 determinants/inverses, transformations, camera/projection helpers |
| `gem.quaternion` | 65 | Quaternion arithmetic, rotation, matrix conversion, unit-domain powers/log, interpolation |
| `gem.bezier` | 12 | Quadratic/cubic evaluation and BezierPath construction/adaptive sampling |
| `gem.legendre` | 6 | Unnormalized ordinary/associated Legendre functions and historical scratch helpers |
| `gem.spherical_harmonics` | 18 | Real SH basis, samples, RGB projection/reconstruction/convolution, scalar/RGB analytical L0–L2 rotation, legacy probes |
| `gem.plane` | 15 | Plane coefficients, points/polygon normals, normalization, evaluation, orientation helpers |
| `gem.ray` | 7 | Ray state, deep duplication, rigid rotation and pure translation |
| `gem.common` | 12 | Numeric, angle and ctypes/list helpers |

There are 17 runtime Python files: nine nonempty core modules, the empty root
initializer and seven experimental package/shim files. `Vector2` and `Matrix4`
are descriptions of dimensions, not separate public classes. There is no core
ray-intersection query, plane projection/distance method, public Bezier derivative
API, new spline family, collision primitive or spatial index. Tutorial-local
analytical functions are not promoted into the public API. Retired shadow
transport remains absent; core runtime imports do not depend on experimental shims.

The integrations retain these specific boundaries:

- Matrices store rows; a wrapper's `Matrix*Vector` calculates a row-vector
  product. Composition applies the left matrix first. Hamilton `q2*q1` applies
  q1 then q2, so its row matrix is `M(q1)*M(q2)`.
- Positions get local w=1 only in documented transform/ray helpers. Directions
  use w=0; general matrix multiplication requires matching dimensions.
  Plane normals under an affine linear map use inverse transpose, not a position
  transform. Project/unproject perform homogeneous division separately.
- Quaternions use [w,x,y,z]. Axis-angle constructors use degrees; historical
  axis-specific helpers retain their own documented units. Quaternion forward
  +Z and Vector front −Z are intentional. Quaternion/Vector multiplication is
  not a substitute for the conjugate-sandwich rotation helper.
- Returning operations allocate independent results where promised, while
  constructors may retain caller lists/Vectors. Ray construction retains and
  normalizes direction storage; transforms preserve distance and `.end` as
  intersection state. This audit does not redefine it as an endpoint.
- Canonical real SH uses Condon–Shortley signs and index l*(l+1)+m. Historical
  probe coefficients require explicit conversion. Rotation is active:
  f_rotated(d)=f_original(R^-1 d); diffuse convolution is applied once.
  Reconstruction uses unit directions without implicit normalization.
- Ordinary Legendre extrapolation, associated-domain restrictions and explicit
  scratch mutations remain distinct. Legacy SQUAD is sign-sensitive and need
  not be unit length; `squad4` and shortest-path SLERP retain their unit-input
  contracts. No arbitrary branch canonicalization is imposed.

## Independent integration coverage

[New tests](../tests/test_final_integration_audit.py) add **139 passing cases and
12 strict failures**. Many cases exercise multiple operations and deterministic
inputs. Independent oracles include Fraction Bernstein polynomials, rational
row products and Gauss–Jordan inverse, explicit Rodrigues rotation,
Cartesian L0/L1/L2 SH formulas, analytical Snell refraction and camera mapping.
Cross-module agreement supplements these references rather than replacing them.

| Integration | Independent evidence and boundaries |
| --- | --- |
| Quaternion → Vector/Matrix → quaternion and composition | Seed 475000+case; axes scaled 1e−300/1/1e300, signed angles and noncommuting rotations; rational inverse and float32 export |
| Rotation → reflection/refraction/Ray | Analytical rotated geometry, ordinary Snell cases, origin rotation, local translation, distance and deep-copy ownership |
| Affine transform → Plane | Nonuniform scale/shear and signed determinant; exact dual-normal/offset, incidence, polygon winding |
| Camera pose → perspective/orthographic projection | Analytically calculated window values; noncommuting pose, offsets/non-square viewport, four raw/wrapper combinations, independent unprojection inputs |
| Vector curves → affine maps/adaptive paths | Quadratic/cubic Fraction references, Vector2/3/4, extrapolation/endpoints, transformed samples, ordered parameters, repeated output and control ownership |
| Extreme inverse → Vector application | Matrix3/4 scales 2^-996 through 2^996, opposite-scale vectors, nonsymmetric inputs, inverse application, rational inverse, in-place/returning buffers |
| Interpolation → Matrix/SH | Known-axis angle polynomials including t=0, 1e−12 and nextafter(1,0), sign-equivalent controls, active rotation/reconstruction and convolution |
| Near poles → Quaternion/Matrix/SH | Both poles, transverse components 1e−9 through minimum subnormal, exact half-turn, ulp-aware representability |
| Legendre → SH | Explicit P0–P3 and P2^1 polynomials, phase/normalization, preserved scratch state; prior high-degree independent suites retained |
| Projection → rotation → diffuse lighting | Independent six-axis cubature/basis, RGB moment references, analytic irradiance, explicit native-float32 legacy conversion/pixel-center quadrature |
| Shared Vector/Quaternion storage | Constructor alias retention, returning independence, receiver mutation, immutable ctypes snapshot |
| Exact degeneracy across modules | Direct zero-norm fallbacks versus preserved geometric ZeroDivisionError and checked SH ValueError behavior |

All prior extreme-range, near-singular, repeated-transform, shader, sampling,
exception and ownership regressions remain in the complete suite. No test,
tolerance, existing marker, implementation, signature or dependency is changed.

## Confirmed defect: 4G5-A01 — subnormal interpolation endpoint

**Severity P3; release blocker under the final correctness gate.** Source:
[quat_slerp lines 289–293](../gem/quaternion.py#L289), propagated by
[squad4 lines 331–333](../gem/quaternion.py#L331).
Public free SLERP, `Quaternion.slerp` and conventional SQUAD are affected.

Minimal reproduction, using a supported axis-angle constructor:

```python
import math
from fractions import Fraction
from gem.quaternion import Quaternion, quat_from_axis_angle, quat_rotate_vector
from gem.vector import Vector

u = math.ulp(0.0)  # 2^-1074
end = quat_from_axis_angle([1.0, 0.0, 0.0], math.degrees(2*u))
assert end.data == [1.0, u, 0.0, 0.0]
assert math.hypot(*end.data) == 1.0
out = Quaternion().slerp(end, 1.0)
expected_z = float(2*Fraction(end.data[0])*Fraction(end.data[1]))
print(out.data)  # actual: [1.0, 0.0, 0.0, 0.0]
print(expected_z)  # 9.881312916824931e-324, nonzero and representable
print(quat_rotate_vector(out, Vector(3, [0.0, 1.0, 0.0])).vector)
# actual [0.0, 1.0, 0.0]; expected [0.0, 1.0, expected_z]
```

The endpoint is an ordinary binary64 representation of a unit rotation, produced
by the public constructor. Its cosine rounds to one. No input domain expansion
or whole-quaternion normalization is needed. Direct endpoint rotation preserves
Z=2*w*x; Fraction establishes that value independently. Negating the quaternion
does not change the expected rotation, and either rotation sign fails.

SLERP computes difference norm u and sum norm 2, then theta=2*atan2(u,2).
The half-angle u/2 is below the binary64 subnormal grid and rounds to zero;
multiplication afterward cannot recover it. The theta==0 branch returns the
start orientation even at t=1. This is premature intermediate underflow, not
an unrepresentable final result, a half-turn tie or ill-conditioned geometry.
Absolute error is about 1e−323, with negligible practical visual effect; the
available endpoint component is nevertheless lost completely.

Coverage gap: prior nearby-angle tests stop far above this angle; powers/log
and near-pole SH repairs test subnormals but do not feed them into interpolation.
Twelve new strict cases cover three entry points × two rotation directions ×
two quaternion signs. They also verify fresh result storage and preserved inputs.
With markers disabled, all 12 fail exactly the final independent rotated-component
assertion. The same reproduction fails in the separately installed wheel and sdist.

Recommended separate repair: preserve exact shortest-path endpoint rotations
with independent copies at t=0/1, and investigate a scale-aware continuous small-
angle blend so an intermediate half-angle cannot erase a representable result.
Verify interior t, nonidentity controls, signs and SQUAD propagation as well as
endpoint cases. Do not infer that every subnormal midpoint must remain nonzero;
genuine output underflow is unavoidable. Preserve the branch convention and avoid
epsilon-based equality or implicit normalization. Benchmark the affected paths
after selecting a repair; no repair-cost estimate is claimed here.

## Performance review

Prior matched repair reports already quantify the costs of corrected arithmetic.
These are their own same-environment comparisons, not a new combined benchmark:

| Repair report | Recorded tradeoff |
| --- | --- |
| 4G-1R | Exact cyclic power at exponent 8 costs 3.92–4.09×; subnormal half-power 1.45–1.50× and log 1.61–1.63× |
| 4G-2R | Ordinary scalar quadratic/cubic roughly 16–28% slower; tiny-parameter repaired paths roughly 4.2–5.5×; P2 overflow retry about 2.5× including construction |
| 4G-3R | Nonzonal scalar L1/L2 about 8.5–19.1% slower; cached ordinary L2 reconstruction not distinguishable consistently from noise |
| 4G-4R | Triangle normal 8–11% slower; 64/256 vertices 23–45%; refraction differences inconclusive within paired variation |

These changes pay for validated results, rather than comparable historical
accuracy. There is no evidence here to justify reverting the numerical repairs.
Cross-host timing records are not treated as controlled comparisons.

### 4G5-P01 — whole-path validation scaling

**P3, post-1.0 improvement.** The earlier curves audit identified repeated
validation as a concern. This audit measures the specific hypothesis: even when
every cubic is straight and needs only endpoints, `getDrawingPoints` validates
all 3N+1 controls once initially and again for each of N segments.
[Call sites](../gem/bezier.py#L169) and
[validation](../gem/bezier.py#L219) cause (N+1)(3N+1) control visits.

The [dedicated harness](../benchmarks/final_integration_audit.py) excludes
construction/reference checks, warms once, disables GC only in timing, uses
process CPU time, shuffles workload order with seed 475005 and performs seven
trials. Counts per trial are 16/8/2/1 calls for N=8/32/128/512. Separate cProfile
calls verify N+1 full-coordinate validations; outputs are independently checked
against [3*i,0]. No curved/depth-cap growth confounds this measurement.

| Straight segments | Output points | Control visits | Run 1 median ± MAD (ms) | Run 2 median ± MAD (ms) |
| --- | ---: | ---: | ---: | ---: |
| 8 | 9 | 225 | 0.130 ± 0.006 | 0.141 ± 0.013 |
| 32 | 33 | 3,201 | 1.268 ± 0.010 | 1.301 ± 0.066 |
| 128 | 129 | 49,665 | 17.273 ± 0.147 | 18.618 ± 1.837 |
| 512 | 513 | 788,481 | 275.380 ± 3.720 | 295.773 ± 29.666 |

Increasing 128→512 segments multiplies output by about four and measured time
by about sixteen in both runs. Full-path coordinate validation accounts for
about 98% of the separate 512-segment profiled time in run one. Profiling adds
overhead and its absolute times are not benchmark latencies. Raw trials, MAD,
range and profile counts are retained in JSON. This is a measured scalability
concern for rebuilding long paths, not a new mathematical defect or timing gate.
Future work can validate once per invocation while respecting caller-owned
mutable controls and historical invalid-input behavior. No caching or optimization
is introduced now. Shared-host scheduling/frequency and Python differences prevent
universal latency claims.

## Documentation and distribution checks

Both isolated installations pass seven mathematical/ownership smoke groups,
268 runtime declarations, seven constants and compatibility identity checks.
The source and each installed environment also execute all **42 Python snippets**
and check **569 local links/anchors** across the selected public documentation.
Installed execution uses `python -I` from `/tmp`; imported `gem` paths are under
the respective venv site-packages, not the checkout. The audit tool reads source
for byte/signature comparison but does not import source-tree gem in that mode.

The original `tools/check_architecture_docs.py` command exits 2 with
`out-of-scope modification: gem/bezier.py`. It freezes checkpoint
`dc418923c692b609a9fc611c66433e5950e0a321`, before the legitimate merged repairs.
The [audit adapter](phase4g5_verify.py) supplies this audit's explicit base to the
unchanged checker in memory. API/link/snippet and preservation checks all pass
at that checkpoint: 119 protected runtime/test/example/benchmark/metadata files.
The collector additionally verifies **all 358 original tracked files** are
byte-identical to base. The original checker is not edited and its original
command is not reported as passing. A reusable explicit-base option is a
non-blocking release-tooling follow-up.

Fresh staging builds produce `gem-0.1.12-py3-none-any.whl` and
`gem-0.1.12.tar.gz`. Both contain all 17 runtime files with source-identical bytes,
canonical promoted modules and transitional shims. There are no native extensions,
orphaned retired imports or additional runtime requirements: `six` remains the
sole dependency. Both installed environments pass `pip check`.
Setuptools normalizes the sdist's setup.cfg by removing comments and adding empty
tag_build/tag_date=0 under egg_info; the declared metadata and repository file
are unchanged. Artifact hashes/member lists and the normalization are recorded.
Build-tool/timestamp differences may change archive hashes on another run.

Existing packaging is not release-ready metadata: version 0.1.12 and historical
URLs/classifiers remain deliberately unchanged. README.rst still advertises
unverified PyPy/Python2 support and obsolete feature descriptions, unlike the
current Markdown documentation. The sdist includes only the two selected
migration guides, not the complete reference or docs/QUATERNIONS.md linked by
README.rst; the wheel contains runtime modules and metadata, not tutorials.
These are previously documented release-engineering risks, not a runtime import
failure. Python 2.7 has known source blockers, including raise-from syntax and
unavailable APIs; it was not executed. Other interpreters/platforms remain
unverified. No compatibility claim is expanded by this audit.

### HDR reference regeneration

The existing CPU example was executed into a separate directory. All seven
outputs are generated from core SH, and committed references remain untouched.
Both 192×192 PNGs match bytes and decoded RGB exactly:

| Image | Committed and regenerated SHA-256 |
| --- | --- |
| sh_original.png | `dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a` |
| sh_rotated.png | `c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7` |

JSON/GLSL differ in the already recorded last bits from 4G-3R:
57 numeric fields in reference.json (max 8.88e−16), 66 in visualization.json
(max 1.07e−14), and 16/9/9 GLSL constants (max 3.67e−17 / 4.03e−18 / 4.03e−18).
Shapes and nonnumeric content match. PNG channels range 0–211, mean 102.349465.
Linear lighting centroids move from (122.8861,95.5) to (95.5,68.1139), consistent
with active +90° lighting rotation. Existing independent linear-pixel/GLSL tests
pass; hashes are supplemental host-specific reproducibility evidence. The older
byte-only HDR verifier is not a portable mathematical oracle.

## Release-readiness matrix and remaining limitations

| Area / concern | Classification | Evidence and disposition |
| --- | --- | --- |
| Mathematical correctness: interpolation endpoint | **Release blocker** | 4G5-A01, P3; 12 strict failures and installed reproductions. Resolve and rerun before signoff |
| Cross-module order/normal/basis mismatch | Not reproducible / false positive | Independent integration references pass; row/Hamilton composition and inverse-transpose normals agree |
| Numerical boundaries promised by prior repairs | Not reproducible / false positive | Original extreme-scale inverse/norm, powers/log, Bezier/Legendre, near-pole and grazing/Newell tests remain passing |
| Broader extreme arithmetic | Documented limitation | Unscaled Matrix2 inverse/determinants, quaternion inverse/dot/cross, severely conditioned cofactors, extreme barycentric/plane/reflection products retain explicit limits |
| Very small nonzero projection W | Release risk, non-blocking | Prior 4G-1 policy question: reciprocal may overflow although component quotients are finite; zero-W exceptions/sentinel are established, broader policy is not resolved here |
| Public API/ownership inconsistency | Not reproducible / false positive | 268 declarations and independent ownership tests agree; constructor aliases, unit prerequisites and legacy return types are intentional |
| Current reference/docs execution | Not reproducible / false positive | Links, declarations, aliases and 42 snippets pass in source and both installed environments with explicit current checkpoint |
| Historical docs checker checkpoint | Release risk, non-blocking | Original CLI fails against obsolete base; audited adapter passes without changing production or historical tooling |
| Metadata/distribution documentation | Release risk, non-blocking | Historical version, URLs/classifiers/README.rst and incomplete packaged docs require release-engineering reconciliation; no unverified interpreter claim |
| Runtime dependencies/imports/distributions | Not reproducible / false positive | Pure wheel, six only, 17 source-identical runtime files, shims equivalent, E07 absent, no orphan runtime imports |
| Long Bezier paths | Post-1.0 improvement | 4G5-P01: two-run quadratic validation evidence; numerical output and finite subdivision remain correct |
| SH truncation, quadrature and extreme orders | Documented limitation | L2 ringing, finite-resolution integration, factorial/unnormalized intermediate range and extreme coefficient products are previously bounded |
| Other versions/platforms and actual GPU execution | Documented limitation | CPython3.12/Linux only; CPU/GLSL reference checks do not establish graphics-driver or other-platform support |

Specific retained limitations include Matrix3/4 cofactor cancellation on severely
ill-conditioned represented matrices; inverse2 and quaternion inverse at squared
extreme ranges; dot/cross/barycentric underflow; unrepresentable Newell products;
overflow of polygon offset means; unequal-index near-critical refraction policy;
huge arbitrary quaternion phase error; high-order SH factorial/recurrence limits
and extreme color×basis×weight products. The public operation-specific boundaries
and prior exact reproductions remain authoritative. None is silently converted
into an epsilon policy or a new unsupported-domain regression.

Adaptive sampling's depth-16 cap may emit a best approximation beyond tolerance,
and source thinning thresholds are heuristics. Direct list edits need not update
ctypes snapshots. Ray `.end` transformation/actual-hit interpretation remains an
API design question. Nonunit powers/logarithms, nonunit rotation inputs, matrix
SH orientations, higher-band analytical rotation and new geometry APIs remain
outside supported contracts. These are not manufactured new findings.

The confirmed endpoint issue is the only new mathematical blocker found in this
audit. Packaging policy and broader numerical questions remain visible for later
review; passing bounded tests is not unconditional release readiness.

## Executed verification and reproduction commands

Environment: CPython **3.12.14**, Clang22.1.3, Linux6.18.44 x86_64/glibc2.41,
Intel Xeon Platinum8370C shared host, pytest9.1.1, six1.17.0.
Build tools: setuptools84.0.0, wheel0.48.0, packaging26.3.
No other interpreter or OS was tested.

| Run | Passed | Xfailed | Failed/errors | Skipped |
| --- | ---: | ---: | ---: | ---: |
| Unchanged baseline | 3,714 | 0 | 0 | 0 |
| New integration file | 139 | 12 | 0 | 0 |
| Confirmed reproduction, markers disabled | 0 | 0 | 12 intentional failures | 0 |
| Final complete suite | 3,853 | 12 | 0 | 0 |

Baseline console time was 16.67s; final 16.15s. These single-suite wall times are
not a performance comparison. Reproduction exits 1 and deselects 139 passing
cases. There are no unexpected passes, retired/deleted tests or weakened checks.

Exact test, documentation and measured-performance commands, from the checkout:

```sh
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -o junit_family=xunit1 --junitxml=/tmp/phase4g5-baseline.xml > /tmp/phase4g5-baseline.log 2>&1
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_final_integration_audit.py -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5-focused.xml > /tmp/phase4g5-focused.log 2>&1
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_final_integration_audit.py --runxfail -k subnormal_interpolation_endpoint -q --tb=short -o junit_family=xunit1 --junitxml=/tmp/phase4g5-reproductions.xml > /tmp/phase4g5-reproductions.log 2>&1
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -o junit_family=xunit1 --junitxml=/tmp/phase4g5-full.xml > /tmp/phase4g5-full.log 2>&1
/workspace/.venvs/pyGameMath/bin/python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --output /tmp/phase4g5-source-docs.json > /tmp/phase4g5-source-docs.log 2>&1
/workspace/.venvs/pyGameMath/bin/python audit/phase4g5_verify.py docs --output /tmp/phase4g5-docs.json > /tmp/phase4g5-docs.log 2>&1
/workspace/.venvs/pyGameMath/bin/python benchmarks/final_integration_audit.py --output /tmp/phase4g5-path-cost-1.json
/workspace/.venvs/pyGameMath/bin/python benchmarks/final_integration_audit.py --output /tmp/phase4g5-path-cost-2.json
/workspace/.venvs/pyGameMath/bin/python -m examples.hdr_sh.regenerate --output-dir /tmp/phase4g5-hdr-output > /tmp/phase4g5-hdr.log 2>&1
```

Run the baseline before adding the audit test, or in an unmodified checkout of
the explicit base. Paths identify the executed environment; substitute the
checkout/interpreter paths on another host. The first docs CLI is expected to
fail for the checkpoint reason above; the following audit adapter is the verified
current-base command. No baseline or implementation monkeypatch is used for tests.

Clean artifact staging copied `gem/` excluding pycache/egg-info, unchanged
setup.py/setup.cfg/MANIFEST.in/README.rst/LICENSE/.travis.yml/.landscape.yaml and
the two manifest-selected migration guides into `/tmp/phase4g5-dist-build`.
The staging command was run from the repository root before building:

```sh
python - <<'PY'
from pathlib import Path
import shutil
root = Path('/workspace/pyGameMath')
dest = Path('/tmp/phase4g5-dist-build')
dest.mkdir()
shutil.copytree(root/'gem', dest/'gem', ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.egg-info'))
for name in ['setup.py', 'setup.cfg', 'MANIFEST.in', 'README.rst', 'LICENSE',
             '.travis.yml', '.landscape.yaml', 'docs/EXPERIMENTAL_MIGRATION.md',
             'docs/VECTOR_VIEWPORT_CONTRACTS.md']:
    target = dest/name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root/name, target)
PY
```

Offline build wheels were already cached in `/tmp/phase4d-dependencies`:

```sh
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 -m venv /tmp/phase4g5-build-env
/tmp/phase4g5-build-env/bin/python -m pip install --no-index --find-links /tmp/phase4d-dependencies setuptools wheel packaging six
cd /tmp/phase4g5-dist-build
/tmp/phase4g5-build-env/bin/python setup.py sdist bdist_wheel > /tmp/phase4g5-build.log 2>&1

cd /tmp
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 -m venv /tmp/phase4g5-wheel-env
/tmp/phase4g5-wheel-env/bin/python -m pip install --no-index --find-links /tmp/phase4d-dependencies /tmp/phase4g5-dist-build/dist/gem-0.1.12-py3-none-any.whl
/tmp/phase4g5-wheel-env/bin/python -I /workspace/pyGameMath/audit/phase4g5_verify.py smoke --output /tmp/phase4g5-wheel-smoke.json > /tmp/phase4g5-wheel-smoke.log 2>&1
/tmp/phase4g5-wheel-env/bin/python -m pip check

/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 -m venv /tmp/phase4g5-sdist-env
/tmp/phase4g5-sdist-env/bin/python -m pip install --no-index --find-links /tmp/phase4d-dependencies setuptools wheel packaging six
/tmp/phase4g5-sdist-env/bin/python -m pip install --no-index --no-deps --no-build-isolation /tmp/phase4g5-dist-build/dist/gem-0.1.12.tar.gz
/tmp/phase4g5-sdist-env/bin/python -I /workspace/pyGameMath/audit/phase4g5_verify.py smoke --output /tmp/phase4g5-sdist-smoke.json > /tmp/phase4g5-sdist-smoke.log 2>&1
/tmp/phase4g5-sdist-env/bin/python -m pip check
```

On a clean host, obtain the same build wheels first; six remains the only runtime
dependency. Install/build tooling is not added to gem metadata. Both isolated
smokes intentionally record the blocker as reproduced, not as a passing repair.
After these commands, the standard-library collector inspects archives, checks
all original tracked files, compares regenerated HDR values and gathers the
JUnit/installed/performance evidence:

```sh
cd /workspace/pyGameMath
/workspace/.venvs/pyGameMath/bin/python audit/phase4g5_collect.py --output audit/phase4g5-verification.json
git diff --check
```

## Changed files

| New file | Purpose |
| --- | --- |
| `tests/test_final_integration_audit.py` | Independent integration regressions and 12 strict defect cases |
| `audit/PHASE4G5-FINAL-INTEGRATION.md` | Findings, readiness matrix and executed methodology |
| `audit/phase4g5-verification.json` | Machine-readable measurements/verification |
| `audit/phase4g5_verify.py` | Isolated mathematical/API/documentation smoke adapter |
| `audit/phase4g5_collect.py` | Asserted collection of test/artifact/scope evidence |
| `benchmarks/final_integration_audit.py` | Single-hypothesis Bezier scaling benchmark/profile |

All changes are additive audit/test/development tooling. No production math,
existing tests, examples/assets, packaging metadata, release settings or dependency
files change. The draft audit PR stops at findings for review; no repair or next
phase is started.
