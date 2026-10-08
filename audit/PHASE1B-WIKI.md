# Phase 1B: historical wiki and API validation

## Scope, evidence, and conclusions

All seven available wiki pages were retrieved from the wiki Git repository at **`715e5039c75e080814a12e957f5148c35cdf8bda`**. The wiki tip dates to 2015-11-01, earlier than the audited library tip (2017-06-21, `5257291431bb45db0274dc48edf24694ecfe2e2d`). Byte-for-byte snapshots, SHA256 hashes, and provenance are saved under [wiki-snapshot](wiki-snapshot/manifest.json). Branch inspection found only the master wiki branch; history inspection found no deleted page paths. The four placeholder pages' histories contain only “Coming soon.” or the earlier spelling “Comming soon.” There is no missing substantive quaternion/plane/ray/common documentation in those histories.

This review starts from `974a01c7db4d00065f90896ae41be41e59ffcb6a` on **`audit/phase1b-wiki`**. The original two audit commits, their raw results, and the benchmark baseline are preserved. No production source, runtime dependency, public namespace, mathematical convention, or PR #9 was changed. No Phase 2 implementation was started.

The wiki strengthens several existing API conclusions but does **not** define storage order, row/column multiplication, quaternion ordering, angle units, plane equations, refraction ratio direction, or numerical tolerances. Matrix/vector conventions remain as established by source and executable probes. “Not documented” is distinct from “wrong.”

### Complete page inventory

| Page | Content examined | Consequence |
| --- | --- | --- |
| [Home](wiki-snapshot/Home.md) | Library purpose, OpenGL origin, graphics/physics scope, README pointer | Context only; no new mathematical contract. “Coming soon” statement is stale relative to existing class pages. |
| [Vector Class](wiki-snapshot/Vector-Class.md) | Operators, free functions, every listed method, swizzles, axis identities, examples, non-i/i rule | Main historical evidence for vector API. |
| [Matrix Class](wiki-snapshot/Matrix-Class.md) | Operators, projection functions, every listed method, identity helpers, examples, dimension restriction, non-i/i rule | Main historical evidence for matrix API. |
| [Quaternion Class](wiki-snapshot/Quaternion-Class.md) | Exactly `Coming soon.` | No historical contract; source/math evidence remains necessary. |
| [Plane Class](wiki-snapshot/Plane-Class.md) | Exactly `Coming soon.` | No plane equation or normal-field policy. |
| [Ray Class](wiki-snapshot/Ray-Class.md) | Exactly `Coming soon.` | No ray promotion, copying, or distance policy. |
| [Common Functions](wiki-snapshot/Common-Functions.md) | Exactly `Coming soon.` | No angle conversion or viewport specification. |

Original website: <https://github.com/AlexMarinescu/pyGameMath/wiki>. Frozen local evidence is used for reproducibility rather than assuming later website content is unchanged.

## Wiki-to-source API compatibility matrix

“Agrees” means the stated requirement is supported on tested ordinary inputs, not that every numerical edge case passes. “Underspecified” means the wiki cannot determine the disputed policy. Source links point to unchanged implementation locations; tests are in [test_wiki_contracts.py](../tests/test_wiki_contracts.py) or the existing Phase 1 suite.

### Vector API: every documented entry

| Documented entry / wiki lines | Intended requirement | Source and audit result |
| --- | --- | --- |
| `Vector(size, data)`; 65–85 | 1D through ND, default zero, supplied list data; dimension-limited methods allowed | [vector.py:160](../gem/vector.py#L160): agrees for 1,2,3,4,8 dimensions. Validation, ownership of supplied data, and size zero are not specified. |
| `+`, `−`; 6–7, 88, 95–97 | Component addition/subtraction of equal-size vectors | [vector.py:171](../gem/vector.py#L171), :191: examples agree. Scalar addition/subtraction are source extensions, not documentation violations. Mismatch precondition is explicit; exception type is not. |
| `*`, `/`; 8–9, 99–101 | Scalar multiplication/division | [vector.py:211](../gem/vector.py#L211), :232: examples agree on Python 3. Operator handedness and scalar type protocol beyond examples are unspecified; do not infer mandatory `scalar*Vector`. |
| Unary `−`, `negate()`; 10, 24, 103 | Negate every component | [vector.py:273](../gem/vector.py#L273), :288: agrees. Returning operations preserve receiver. |
| `==`; 11 | Boolean equality | [vector.py:253](../gem/vector.py#L253): only first component compared (V01). Equal-size equality is a confirmed defect; mixed-size and empty-vector policies are not settled by wiki. `!=` is later API, not documented here. |
| Free `lerp(a,b,time)`; 14–15 | Linear interpolation; return same dimension | [vector.py:37](../gem/vector.py#L37): agrees for 1–8 dimensions. Clamping/extrapolation policy not stated; existing extrapolation is preserved. |
| Free `cross(a,b)`; 16, 107–119 | 3D cross product; example `[-17,−17,17]` | [vector.py:43](../gem/vector.py#L43): agrees, confirms sign orientation for that example. Rejection behavior for non-3D input is unspecified. |
| Free `reflect(incident,norm)`; 17 | 3D reflection | [vector.py:51](../gem/vector.py#L51): ordinary unit-normal example passes. Normalization responsibility unspecified. |
| Free `refract(IOR,incident,normal)`; 18 | 3D refraction with that argument order | [vector.py:55](../gem/vector.py#L55): V02 confirmed formula/runtime failures. Wiki does not specify eta=n1/n2 versus reciprocal or TIR policy. |
| `clone()`; 21, 130–141 | Exact copy; independent normalization example | [vector.py:277](../gem/vector.py#L277): agrees; cloned data independent. |
| `one()`, `zero()`; 22–23, 73–85 | Fill every component with 1 or 0 | [vector.py:280](../gem/vector.py#L280), :284: agrees, returns receiver. Neither has an i-prefixed counterpart; “most functions” prose does not require one. |
| `maxV(other)`, `minV(other)`; 25,27 | Component extrema between two vectors | [vector.py:292](../gem/vector.py#L292), :298: agrees. |
| `maxS()`, `minS()`; 26,28 | Unary component extrema inferred from signatures | [vector.py:295](../gem/vector.py#L295), :301: agrees. “Between two vectors” descriptions conflict with unary signatures (D06), not evidence of a missing second argument. |
| `magnitude()`; 29 | Scalar magnitude | [vector.py:304](../gem/vector.py#L304): ordinary cases pass; N02 extreme-scale instability remains. |
| `clamp(size,value,minS,maxS)`; 30 | Clamp components; generic returning/in-place distinction | [vector.py:307](../gem/vector.py#L307), :311: receiver rule passes. Ownership of the separate `value` argument is unspecified; V04 is a contract question. |
| `normalize()`, `i_normalize()`; 31, 122–141 | New normalized vector versus normalization in place | [vector.py:316](../gem/vector.py#L316), :321: corrected numerical example agrees. Zero normalization policy absent (N01 uses source fallback evidence). “Without returning a new vector” permits returning `self`; it does not require `None`. |
| `dot()`; 32,105,118 | Two-vector scalar dot, example 57 | [vector.py:326](../gem/vector.py#L326): agrees; abbreviated `dot()` listing omits required second-vector argument (D08). |
| `isInSameDirection(other)`, `isInOppositeDirection(other)`; 33–34 | Return booleans | [vector.py:333](../gem/vector.py#L333), :340: returns dot-sign hemisphere checks. Wiki gives no parallel-versus-hemisphere requirement. |
| `barycentric(a,b,c)`; 35 | Coordinates relative to triangle | [vector.py:347](../gem/vector.py#L347): known-answer test passes. Degenerate-domain policy absent. |
| `transform(position,matrix)`; 36 | Matrix-based transformation | [vector.py:371](../gem/vector.py#L371), helper :150: identity fails V03. Input wrapper/list forms, homogeneous behavior, and convention remain unspecified by wiki. |
| `xy`, `yz`, `xz`; 41–43 | Return 2D component groups | [vector.py:383](../gem/vector.py#L383): agrees. |
| `xw`, `yw`, `zw`; 44–46 | Return 2D groups requiring W component | [vector.py:392](../gem/vector.py#L392): agrees for Vector4; page's “depends on case” does not make W available on Vector3. |
| `xyw`, `yzw`, `xzw`; 47–50 | Return 3D groups | [vector.py:401](../gem/vector.py#L401): agrees; `xyw` repeated and `xyz` omitted (D07). No extra duplicate API is required. |
| `right`, `left`, `front`, `back`, `up`, `down`; 55–60 | Explicit +X, −X, **−Z**, +Z, +Y, −Y identities | [vector.py:414](../gem/vector.py#L414): all agree. Historical support for Vector front=−Z is explicit; Quaternion forward=+Z remains separate. |
| General non-i/i convention; 122 | Ordinary method returns new result; i-method changes receiver | Ordinary normalize/clone/negate and i-normalize agree. Do not turn i-method `self` returns into `None`, nor infer blanket copying of every argument. |

### Matrix API: every documented entry

| Documented entry / wiki lines | Intended requirement | Source and audit result |
| --- | --- | --- |
| `Matrix(size,data)`; 41,47–73 | NxN with dimension-specific limits; explicit nested matrix data | [matrix.py:345](../gem/matrix.py#L345): explicit constructors agree. Default-zero printed output conflicts with contemporaneous identity source (D01); preserve default identity. |
| `Matrix*Matrix`; 6,76,87 | Same-size matrices, Matrix result | [matrix.py:356](../gem/matrix.py#L356): matched-size cases pass; mismatch raises ValueError. Storage/composition not defined by identity-only examples. |
| `Matrix*Vector`; 6,76,90 | **Vector dimension must match matrix size**, Vector result | [matrix.py:28](../gem/matrix.py#L28): matched-size examples pass. Wiki does not imply column-vector algebra merely by spelling Matrix first. Implicit promotion is not established. |
| `Matrix/scalar`; 7,93 | Element scalar division; float example 2.0 | [matrix.py:384](../gem/matrix.py#L384): M03 Python 3 operator failure; M02 transpose in helper. Integer support not demonstrated. |
| Free `orthographic`; 13 | 4×4 orthographic Matrix | [matrix.py:604](../gem/matrix.py#L604): type/depth known-answer tests pass. No near/far validation or coordinate contract provided. |
| Free `perspective`, `perspectiveX`; 14–15 | 4×4 perspective; **FOVY versus FOVX** | [matrix.py:616](../gem/matrix.py#L616), :634: agrees; asymmetric aspect=2 verifies difference. Wiki does not state degrees. |
| Free `lookAt`; 16 | 4×4 view Matrix | [matrix.py:657](../gem/matrix.py#L657): ordinary tests pass; degenerate input policy absent. |
| Free `project`; 17 | **3D Vector result**, viewport projection | [matrix.py:670](../gem/matrix.py#L670): P01 crash. Three-component storage requirement clarified; input forms and depth remap still undecided. |
| Free `unproject`; 18 | Screen-to-object coordinates | [matrix.py:693](../gem/matrix.py#L693): identity passes; noncommuting transform fails P02 under established row-vector convention. Depth range is not specified by wiki. |
| `scale(value)`; 21,104–128 | Scale by Vector, returning/new versus i/in-place; determinant 24 example | [matrix.py:398](../gem/matrix.py#L398), :409: corrected 4×4 example agrees. Listing allows dimension limits, so arbitrary-N scaling is not newly required. |
| `det()`; 22 | Scalar determinant | [matrix.py:419](../gem/matrix.py#L419): 2–4 known answers pass; unsupported-N policy not specified. |
| `inverse()`; 23 | Matrix inverse | [matrix.py:445](../gem/matrix.py#L445): M01 2×2 defect strengthened; ordinary 3×3/4×4 pass. General-N support not promised. |
| `rotate(axis,theta)`; 24 | Axis-angle rotation | [matrix.py:473](../gem/matrix.py#L473): ordinary 3D/4D tests pass; M06 pivot helper error remains source/math-derived. Units and 2D pivot interpretation absent. |
| `translate(vecA)`; 25 | Vector translation; returning/new and i/in-place versions | [matrix.py:488](../gem/matrix.py#L488), :510: M04 same-input inconsistency confirmed. Wiki's dimension rule concerns multiplication, not the translation offset size; Vector3 offset in a 4×4 transform is not disproved. |
| `transpose()`; 26 | Matrix transpose | [matrix.py:535](../gem/matrix.py#L535): agrees. |
| `shearXY`, `shearYZ`, `shearXZ`; 27–29 | Named shear operations with listed argument order | [matrix.py:539](../gem/matrix.py#L539): M05 3×3 crash; other helpers callable. Wiki's short wording does not settle which coordinate supplies each shear displacement. |
| Free `identity(size)`, `zero_matrix(size)`; 33–36 | Nested matrices passed into Matrix constructor | [matrix.py:6](../gem/matrix.py#L6), :10: agree; not class methods. |
| General non-i/i convention; 96 | New Matrix vs changed receiver | All eight documented mutator pairs tested on supported 4×4 inputs; returning `self` from i-methods is compatible. Output printed after i-scale is wrong (D04). |

Undocumented extensions such as scalar vector addition, `!=`, `xyz`, ctypes export, and quaternion helpers are not removed or renamed merely because the wiki omits them. Experimental algorithms have no wiki page and no newly inferred historical requirements.

## Documentation mistakes, kept separate from implementation findings

| ID | Frozen evidence | Classification and handling |
| --- | --- | --- |
| D01 | Matrix 60–63 prints zeros for `Matrix(4)` | Documentation/source conflict. Source from 2015 commit `63f3630c1a2973add885be4d04a734a65243ae69` already defaults to identity, as does audited source. Treat printed output as an error with strong historical evidence; no default change authorized. |
| D02 | Matrix 82–83,103 calls `matrix.indetity` | Typo; actual API is `identity`. Correct only the local test transcription, not the API or archived page. |
| D03 | Matrix 103 constructs size 3 from identity(4), then prints 4×4 output | Inconsistent example dimensions. Corrected test uses Matrix4 to match data/output. Enforcing constructor shape is a separate API decision. |
| D04 | Matrix 109–119 applies i-scale then prints original identity, yet determinant is 24 | Output contradicts in-place prose and determinant. Corrected expected receiver is diag(2,3,4,1). Existing in-place tests should not be changed to preserve identity. |
| D05 | Vector 128 uses `data[2.0,4.0,8.0]` | Incorrect expression/missing keyword `=`; corrected test uses `data=[...]`. As written it may parse but references/indexes `data`, not a valid documented constructor invocation. |
| D06 | Vector 26,28 says maxS/minS work “between two vectors” | Description conflicts with zero-argument signatures; unary extrema matches source. Do not add a mandatory vector argument. |
| D07 | Vector 47,50 repeats `xyw`; omits `xyz` | Duplicate/omission. Existing `xyz` stays public; no mathematical defect. |
| D08 | Vector 32 lists `dot()` without argument, but example uses `.dot(vectorB)` | Abbreviated/incomplete signature, resolved by the correct example. |
| D09 | Home says individual class pages are “Coming soon” | Stale status text; two pages do exist. Four actual placeholder pages remain incomplete, not hidden API contracts. |

Rounded decimal outputs for vector division/normalization are display approximations, not demands for truncating mathematical results. “Without returning a new vector/matrix” is not a promise to return `None`. Corrected examples are labeled in tests; frozen snapshots are not edited.

## Revisit of all 44 Phase 1 finding groups

Evidence strength is distinct from severity. “Math/source” findings remain valid without wiki detail; “contract question” means a failing assertion alone does not establish an approved behavior change. This table supersedes any earlier implication that all 44 groups are equally confirmed implementation defects.

| ID | Historical evidence and revised classification | Test/implementation consequence |
| --- | --- | --- |
| V01 | Vector boolean equality documented; **confirmed** on equal-size ordinary vectors. Dimension/empty policies not specified; same-size precondition explicit. | Keep two equality/inequality component regressions as defects; move dimension and empty hypotheses to contract-question marking. |
| V02 | Refraction explicitly listed; math/runtime defect **confirmed**. Ratio direction and TIR return policy unspecified. | Keep reproductions; approve ratio semantics before documenting a correction. |
| V03 | Transform documented as matrix transformation; identity violation **confirmed**. | Keep identity reproduction; do not use wiki to switch conventions or choose wrapper inputs. |
| V04 | Clamp plus receiver non-i/i rule documented, but ownership of separate `value` list unspecified. **Reclassify as ownership contract question**, not an approved defect fix. | Keep observation under contract-question marker; add passing receiver-ownership test. |
| M01 | Inverse explicitly documented; 2×2 result **confirmed defect**. | Keep independent inverse regressions. |
| M02 | Scalar division documented; transpose violates scalar division: **confirmed defect**. | Correct with inverse4 dependency in the same later PR. |
| M03 | Literal `matrixB / 2.0` example strengthens **confirmed Python 3 operator defect**. | Add exact documented division reproduction after fixing only `indetity` typo. |
| M04 | Returning vs in-place semantics support **confirmed same-input inconsistency**. Multiplication dimension rule does not forbid a 3D offset for a 4×4 transform. | Keep i-translate regression; do not misclassify its offset as mismatched matrix-vector multiplication. |
| M05 | Shear APIs listed; 3×3 indexing crash **confirmed**. Exact shear-plane mapping not resolved. | Keep crash test; approve mapping before rewriting XY4. |
| M06 | Rotation named, but no 2D pivot semantics. Pivot invariance defect remains **math/source-derived**, not wiki-proven. | Keep helper probe; require explicit 2D rotation contract before public API correction. |
| M07 | No ctypes contract in wiki; stale secondary state remains **source consistency defect**. | Keep buffer-sync regression; preserve export API. |
| N01 | Zero normalization not covered; guard/fallback show **source-level zero-division defect**. | Keep reproductions; existing zero/identity fallback remains source evidence, not wiki evidence. |
| N02 | Magnitude named, no scale guarantees; overflow/underflow **confirmed numerical issue**. | Retain robustness tests; tolerance/domain policy remains separate. |
| N03 | Inverse named, no scaling guarantees; scaled identity **confirmed numerical issue**. | Retain scale probes; no new inverse algorithm approved. |
| P01 | Wiki explicitly says project returns **3D Vector**. Crash and size-3/four-value return are **confirmed**; argument types and [0,1] depth are not specified. | Add return-type/three-storage-components regression. Move old center/depth/type hypotheses to contract questions. |
| P02 | Unproject objective documented; inverse-composition defect **confirmed under preserved row vectors**. | Keep noncommuting test; don't use identity example to infer column vectors. |
| Q01 | Quaternion page placeholder; inverse violation **confirmed by Hamilton identity**, no wiki conclusion. | Keep inverse tests and ordering. |
| Q02 | Placeholder; valid rotation matrix failures **confirmed by rotation round trip**. | Keep half-turn tests; q/−q equivalence preserved. |
| Q03 | Placeholder; q^1 and identity-power failures **confirmed mathematical defects**. | Keep tests; general versus unit-only domain undecided. |
| Q04 | Placeholder; unit log scalar component **confirmed mathematical defect**. | Keep unit test; return type/general log domain not settled. |
| Q05 | Placeholder; callable-float exception **confirmed implementation defect**. Three-control spline semantics remain ambiguous. | Keep crash regression; approve intended SQUAD contract before mathematics. |
| Q06 | Placeholder; explicit source list branch crashes **confirmed implementation defect**. | Keep list probes; axis mutation remains ownership question. |
| Q07 | Placeholder; stale matrix buffer **confirmed consistency defect**, not wiki-proven. | Keep snapshot test; public ctypes state preserved. |
| Q08 | Placeholder; `.data` on Vector **confirmed implementation defect**. | Keep in-place product test; product is not a Vector rotation. |
| Q09 | Placeholder; **confirmed Python 3 compatibility gap** relative to legacy division API. | Keep division regression; do not claim wiki documents Quaternion `/`. |
| Q10 | Placeholder; nearby branch deliberately approximates interpolation. **Reclassify strict unit accuracy as numerical contract question**. | Retain measured drift under question marker, not a normative 1e−12 accuracy promise. |
| Q11 | Placeholder; docstring suggests rotation quaternion but code rotates its axis. **Reclassify as unresolved return/semantic contract**. | Question marker; wiki cannot authorize return change. |
| G01 | Plane page placeholder; scalar-coefficient crash **confirmed implementation defect**. Normal length policy absent. | Narrow test to coefficient preservation and normal direction, avoiding an unapproved raw-vs-unit norm requirement. |
| G02 | Placeholder; coefficient-plane normalization and stale normal **confirmed math/state defects**. | Keep offset/normalization tests; decide field policy before implementation. |
| G03 | Placeholder; point fields conflict with coefficient consumers: **confirmed internal representation inconsistency**. | Keep incidence test; decide plane equation and D sign explicitly. |
| G04 | Placeholder; final polygon vertex indexing **confirmed runtime defect**. | Keep polygon probe; closed-list conventions to document. |
| R01 | Ray page placeholder; shared clone and lost distance **source-level copy inconsistencies**, not wiki-proven. | Keep probes; approve broader ownership/distance contract before implementation. |
| R02 | Placeholder; Quaternion replaces expected Vector state: **confirmed type/math defect**. | Keep ordinary quaternion-rotation probe. |
| R03 | Placeholder; origin is not translated, but Phase 1 test also mixes Matrix4 with Vector3. **Split source defect from unapproved implicit promotion requirement**. | Existing failing test moved to contract-question marking. Fix ray dimensional policy before treating automatic promotion as required. |
| C01 | Common page placeholder; conversions disagree with names: **confirmed named-math defect**, no wiki resolution. | Keep tests; do not change other APIs' units. |
| C02 | Placeholder; normalize then subscript incompatible types: **confirmed runtime defect**. | Keep probe; intended viewport formula needs approval. |
| E01 | No experimental page; cubic endpoint violation **confirmed mathematical defect**. | Keep endpoint/interior tests; no historical additions. |
| E02 | No experimental page; curve/Vector operand incompatibility **confirmed integration defect**. Wiki documents Vector*scalar but does not require scalar*Vector. | Prefer a later local operand-order correction over assuming reflected multiplication is mandatory. |
| E03 | No experimental page; float range count **confirmed Python 3 defect**. | Keep reproduction; no inferred curve input-validation policy. |
| E04 | No experimental page; sampling errors **confirmed/source-inspected**, downstream faults still masked. | Keep end-to-end probes; reproduce masked faults before fixing them. |
| E05 | No experimental page; recurrence/idempotence **confirmed mathematical defects**. | Keep low/high-order tests and sign basis. |
| E06 | No experimental page; math.PI and wrong direction field **confirmed source defects**, only first reproduced end-to-end. | Keep generated-sample probe; no additional unit policy from wiki. |
| E07 | No experimental page; object subscripting, unallocated arrays, unfinished shadow transport **confirmed/incomplete implementation**. | Keep failure; do not invent shadowing requirements. |
| E08 | No experimental page; rectangular indexing **confirmed runtime limitation**, rectangular map projection unsupported/undefined. | Keep limitation probe; decide square-only versus rectangular support before correction. |

After this review there are **40 finding IDs with confirmed/source-math defect probes** and **four groups dominated by unresolved contracts** (V04, Q10, Q11, R03). V01 and P01 also contain mixed definite-defect and contract-question subcases. Every original group is retained; no unsupported behavior correction has been silently approved.

## Existing expectations versus historical documentation

No test should adopt the wiki's erroneous default-zero output, identity-after-i-scale output, misspelled function, or incompatible constructor size. Existing default-identity, in-place mutation, and exact arithmetic tests are retained. Matching dimensions are an explicit supported-domain requirement, so `test_equality_dimensions` and the ray Matrix4/Vector3 test cannot be advertised as wiki-required behavior.

Specific changes to existing tests:

| Test | Change and rationale |
| --- | --- |
| `test_equality_dimensions`, `test_empty_equality` | `defect(V01)` -> explicit dimension/empty contract questions. Keep observations; no input support policy chosen. |
| `test_clamp_preserves_input` | `defect(V04)` -> value-list ownership question. Wiki's new-result/receiver rule does not promise deep copying of other arguments. |
| `test_slerp_nearby_unit_length` | `defect(Q10)` -> unit-accuracy question; tolerance is not historical specification. |
| `test_arbitrary_axis_helper_returns_rotation_quaternion` | `defect(Q11)` -> return-semantics question; placeholder wiki adds no authority. |
| `test_ray_translation_moves_origin` | `defect(R03)` -> homogeneous-promotion question; raw multiplication dimension precondition conflicts with its promotion assumption. |
| `test_project_center` (two input forms) | `defect(P01)` -> input/depth questions. Keep newly wiki-supported three-component return test as definite P01 regression. |
| `test_plane_from_coefficients` | Keep G01 crash regression but assert aligned normal, not specifically raw `[0,2,0]`; wiki does not choose unit-versus-coefficient normal. |

Both marker classes execute assertions. Normal runs use strict xfails; `--runxfail` exposes both kinds as failures. Contract questions are explicitly distinguishable in JUnit properties and [Phase 1B results](phase1b-test-results.json); a green-with-xfails run does not mean every proposed expectation is approved.

## Added tests and verification

Added **77 cases** in `test_wiki_contracts.py`: **75 pass**, two reproduce existing M03/P01 failures. They cover all unambiguous numerical examples, constructor dimensions, vector axis identities, listed swizzles, LERP return dimensions, clone/normalize behavior, unary extrema, receiver ownership, Matrix operator dimensions/results, corrected scaling and determinant, non-i/i pairs, projection return types, and FOVY/FOVX distinction. Nontrivial existing math tests supplement descriptions too vague to yield new known answers.

Combined suite: **348 cases, 254 passed, 86 confirmed-defect xfails, 8 contract-question xfails, no skips**. Normal exit 0. `--runxfail`: **94 failed / 254 passed**, expected exit 1. Python 3.12.14, pytest 9.1.1. Statement coverage 83.80%, branch coverage 64.71%. Coverage does not certify unsupported domains or unexecuted downstream experimental code.

No benchmarks were rerun: production implementation did not change, and the existing Phase 1 baseline remains applicable. No production corrections, package modernization, merge, or publication occurred.

## Decisions and subsequent work

See [API questions and recommended Phase 2 order](PHASE2-DECISIONS.md). These are recommendations and decision points, not authorization to start implementation. Primary constraints remain namespace/API preservation, existing row-vector and storage conventions, mixed established angle units, and a pure-Python core.
