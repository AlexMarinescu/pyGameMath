# Experimental package and shadow-transport audit

Base: master `d3bbb68330946e782f60ef9eeca8184511bae721` (merged PR #30).
This phase changes documentation and diagnostic tests only. No implementation,
import path, packaging declaration or expected-failure marker is changed.

## Inventory and destinations

All eight remaining modules are listed below. Reexports have no independent
algorithm or state. Their implementation, validation and ownership contracts
are those of the canonical modules; no further promotion is necessary.

| Experimental module | Complete public surface and purpose | Dependencies / completeness | Recommendation |
| --- | --- | --- | --- |
| `__init__` | Package marker; no public definitions | Empty | Retire with package after compatibility decision |
| `_bezier_legacy` | `BezierPath` compatibility alias (private module but importable) | `gem.bezier`; no implementation | Promote already complete; retire alias after sunset |
| `bezier` | `cubicBezierPoint`, `quadraticBezierPoint`, `BezierPath`: evaluation and adaptive paths | `gem.bezier`; validated core implementation | Promote already complete; retire shim after sunset |
| `legendre` | `Legendre`: ordinary/associated polynomial evaluation | `gem.legendre`; validated core implementation | Promote already complete; retire shim after sunset |
| `sph` | `Factorial`, `K`, `SPH`, incidental `Legendre`: real canonical SH basis | `gem.spherical_harmonics` (and its Legendre import); complete | Promote already complete; retire shim after sunset |
| `sph_sample` | `SPHSample`, `GenerateSamples`: stratified directional quadrature | `gem.spherical_harmonics`; complete, global RNG contract | Promote already complete; retire shim after sunset |
| `sph_irradiance_map` | `SPH_IrradianceMapCoeff`: legacy raw angular-probe radiance coefficients | `gem.spherical_harmonics`; complete under documented probe/basis contract | Promote already complete; retire shim after sunset |
| `sph_object` | `SPHVertex(position, normal)`, `SPHObject(indices, vertices)`, misspelled `GenereateCoeffs(numSamples, numBands, samples, objects)` | `math`, `six.moves`; only remaining algorithm implementation; unfinished and failing | Retire scaffold; Requires decision on symbol removal. Extract any future scene transport into a separate application/library |

Reexported class methods are also part of the compatibility surface:

- `BezierPath`: `__init__`, `setControlPoints`, `getControlPoints`,
  `calculateBezerPoint`, `interpolate`, `samplePoints`, `getDrawingPoints`,
  `findDrawingPoints`, `findDrawingPointsAdded`.
- `Legendre`: `__init__`, `mGreaterThan0`, `calculatePM1`, `calculatePML`, `run`.
- `SPHSample`: `__init__`; stores theta, phi, direction and basis values.
- `SPH_IrradianceMapCoeff`: `__init__`, `load`, `calculateCoefficients`,
  `updateCoefficients`, `output`. Its legacy RGB coefficients require explicit
  basis conversion; the class name does not imply cosine convolution.
- `SPHVertex` and `SPHObject`: only `__init__`. These alias caller inputs;
  coefficient fields begin as `None`. No copying, validation, collision or
  scene management methods exist.

There are no duplicate algorithm bodies left in the shims. The commented
allocation, ray and collision blocks in `sph_object` are unfinished scaffolding,
not callable functionality. Its object/vertex containers have no core equivalent
and provide no reusable geometry operations worth promoting.

## Evidence and consumers

Reviewed historical findings in PHASE1/PHASE1B-WIKI, archived wiki pages,
CONVENTIONS, COMPATIBILITY and PHASE2-DECISIONS (especially QD11), current source
and tests. No historical wiki page documents experimental functions. The
available path history includes `4253839` (2015-10-29, package-path update);
current transport retains the originally audited algorithm and placeholders.
This evidence does not establish a working visibility contract.

`rg -n 'gem.experimental|from gem import experimental' gem tests examples README.rst`
finds test consumers but no core or example dependency. `sph_object` is exercised
by the original E07 test and the new diagnostics only. External callers cannot
be determined from repository searches. Existing imports are therefore a real
compatibility consideration even for code currently unable to produce valid
transport.

Bezier, Legendre and SH tests cover canonical imports and identity-equivalent
compatibility aliases. The E07 test remains a strict xfail. Legacy import
removal could affect qualified imports, incidental exports and downstream
serialization assumptions. Do not infer that unknown consumers do not exist.

## E07: actual mathematical intention

The active inner loop sums `max(normal·direction, 0) * sample.values[l]`.
The intended uniform-sphere scale is `4*pi/numSamples`. This is scalar
**clamped-cosine transfer**, with an unfinished visibility mask suggested by
comments:

`T_lm(p,n) = integral V(p,omega) max(n·omega,0) Y_lm(omega) domega`.

It is not pure visibility projection: cosine is already included. There is no
radiance input, material, bounce, geometry intersection, distance falloff or
actual visibility query. Positions and indices never enter the calculation.
The code is PRT-like scaffolding, not a complete precomputed-radiance-transport
system. No shadowed irradiance result is currently implemented.

If transfer were implemented, irradiance would be the dot product of canonical
RGB **radiance** coefficients with scalar transfer coefficients. Applying cosine
convolution again would double-count it. For unoccluded lighting, the existing
`convolve_diffuse` plus `reconstruct` already evaluates normal-dependent
irradiance; another scene-container API is unnecessary.

### Reproduction and masked faults

Original command:

```
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_experimental.py -k generate_object --runxfail -q -p no:cacheprovider
```

Result: **1 failed, 16 deselected**, `TypeError: type 'object' is not
subscriptable`. Nine new passing diagnostic cases reproduce current faults;
they assert observations, not a desired future transport contract. A test-only
module-global alias bypasses the typo; test fixtures allocate arrays to reach
later stages. Production code is not patched.

| Stage | Observed failure or incorrect behavior |
| --- | --- |
| Entry | Builtin `object[i]` instead of `objects[i]` |
| After typo bypass | Coefficient arrays remain `None`; assignment raises TypeError |
| After fixture allocation | Only the last vertex is rescaled; equal normals produce unequal coefficients |
| Shadowed output | No contribution accumulation; last vertex gets additive `4*pi/N` bias, others remain zero |
| Empty vertex object | UnboundLocalError in final rescaling |
| Zero samples | ZeroDivisionError in rescaling |
| Sample count beyond buffer | IndexError; no validated count contract |
| Changed positions/indices | Results unchanged: no visibility/geometric transport |
| Nonunit normals | Doubled normal doubles transfer; comment's '[0,1]' clamp only tests positivity |

For two identical +Z normals and one +Z sample, independently known
`Y00 = 1/sqrt(4*pi)`: first vertex remains `Y00`, last becomes `sqrt(4*pi)`;
shadowed outputs are respectively zero and `4*pi`. Neither is a valid shadow
calculation. Empty scenes return `None`. Constructors retain references to input
vectors/lists. These observations expose missing policies, not invitations to
invent validation or ownership changes in this audit.

## Core/application boundary and unresolved contracts

Supported core math already contains SH basis evaluation, weighted radiance
projection, reconstruction, diffuse convolution and analytical rotation. Keep
these canonical APIs and their independent known-answer tests.

Scene traversal, acceleration structures, ray origins/self-hit handling,
visibility, materials and multiple-bounce GI belong outside this math package.
A future separate transport implementation would need explicit decisions on:

- Geometry/topology, index interpretation, coordinate frames and unit normals.
- Visibility callback or scene interface, occluders, sidedness, origin offsets,
  self-intersection and hit validity. Ray `.end` remains ambiguous intersection
  state and cannot serve as a newly invented endpoint/hit contract.
- Canonical basis/sign/index ordering; scalar versus RGB transfer; cosine versus
  pure visibility and direct-only versus multi-bounce scope.
- Weighted quadrature versus uniform-sphere sampling; sample counts/bands,
  convergence limitations, empty/invalid inputs and finite domains.
- Allocation, caller ownership, mutation, return values, repeated-call behavior,
  batching, performance and cache lifetime.

Recommendation: retire this unfinished scaffold, preserve its diagnostic evidence
and git history, and design any future visibility transport separately. A small
weighted scalar transfer integrator could be evaluated independently if needed;
no new core API is proposed here. These architecture and removal choices require
a separate decision before implementation.

## Packaging and compatibility

A wheel built with `python setup.py bdist_wheel --dist-dir /tmp/phase2f6-wheel`
includes all eight experimental modules and canonical Bezier/Legendre/SH modules.
The original PHASE1 omission is repaired: setup currently lists both `gem` and
`gem.experimental`. Installing with `--no-deps --target /tmp/phase2f6-installed`
and importing from `/tmp` verified all eight modules load from the installed
wheel and representative aliases are the same core objects. `six` remains the
existing dependency; no dependency is added. No current module omission found.

Build warnings identify existing `description-file` metadata, legacy setup.py
invocation/license classifiers and version normalization. Python classifiers
still advertise old interpreter versions despite modern math APIs. These need
separate packaging modernization; this audit changes none of them.

## Verification

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Complete command:

```
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider --junitxml=/tmp/phase2f6.xml
```

Baseline: **1771 passed, 5 xfailed**. Final: **1780 passed, 5 xfailed,
0 unexpected failures, 0 skips**. Nine added diagnostics account for the change.
JUnit identities match the prior phase exactly:

- `tests.test_experimental::test_generate_object_coefficients` — E07.
- `tests.test_vector_common::test_clamp_preserves_input` — V04 ownership.
- `tests.test_vector_common::test_empty_equality` — unsupported empty dimension.
- `tests.test_vector_common::test_equality_dimensions` — dimension contract.
- `tests.test_vector_common::test_viewport_vector` — C02.

## Phase 2F-7 retirement checklist (decision-gated)

1. Choose compatibility duration and release/version policy for removing imports;
   a transitional release may be necessary before physical package deletion.
2. Confirm retirement of `SPHVertex`, `SPHObject`, `GenereateCoeffs`; archive the
   failure evidence and document the absent visibility implementation. Do not
   present a replacement as supported without its own math/application contract.
3. Migrate ordinary regression consumers to canonical modules, keeping dedicated
   shim tests until the compatibility window closes. Include incidental `Legendre`
   and private `_bezier_legacy` imports in the migration note.
4. Preserve existing core signatures and basis/ownership contracts. Review external
   import and serialization risks; repository search alone cannot rule them out.
5. Resolve the E07 test explicitly as retired functionality: archive reproduction
   evidence and add intentional package/API-absence checks at removal. Do not
   silently delete/skip the xfail. Preserve the other four expected failures.
6. Remove the directory and setup package entry together only after that decision.
   Do not simulate old imports through namespace hacks while claiming retirement.
7. Clean build artifacts, build and install a fresh wheel outside the source tree;
   assert canonical modules work and no experimental files ship or import.
8. Update docs/import examples, scan for stale consumers, run full regression and
   installed-wheel tests, and record compatibility/release notes.

No Phase 2F-7 implementation or package removal is included in this phase.
