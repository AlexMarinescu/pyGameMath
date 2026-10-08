# API questions and recommended Phase 2 order

This is a review checklist, not an implementation plan already in progress. The retrieved wiki resolves vector axis identities and several ordinary API requirements; it does not supply the missing policies below. Phase 2 requires a new instruction.

## Questions requiring approval before their behavior changes

| Question | Existing evidence | Decision needed |
| --- | --- | --- |
| QD01: Supported dimensions and errors | Wiki requires equal dimensions for operators and describes Vector sizes 1…N; source accepts empty/malformed shapes | Is Vector0 supported? Should mixed-size equality return false or reject? What exception should arithmetic mismatches use? Preserve equal-size component semantics. |
| QD02: Input ownership | Wiki distinguishes new result from changed receiver, but says nothing about other mutable inputs | Should returning clamp preserve its separate value list? Should axis conversion and Ray construction copy or mutate caller inputs? Must ray duplication preserve distance/end and fully isolate storage? |
| QD03: Projection API | Wiki requires project to return a 3D Vector; source uses lists for project and wrappers for unproject | Should both functions accept Matrix wrappers, nested lists, or both? Confirm window Z=[0,1], viewport origin, homogeneous input rules, and zero-W behavior. Three-component return is already documented. |
| QD04: Ray dimensional policy | Ray page is a placeholder; Matrix multiplication requires matching dimensions | Should ray methods explicitly support Matrix4 with Vector3 via internal position/direction promotion, or require matching input representations? Define distance/end semantics. Do not add general implicit promotion to Matrix*Vector without approval. |
| QD05: Planes | Placeholder page; scalar consumers use n·p+d=0 but constructors/bestFitD differ | Confirm scalar coefficient representation and d sign; raw versus unit normal field; bestFitD's meaning; degenerate points/polygons and open versus repeated-closure polygon inputs. |
| QD06: Quaternion domains and returns | Placeholder page; source implies unit rotations; helper docstring/implementation disagree | Confirm general versus unit-only inverse/log/power/rotation domains, arbitrary-axis helper return semantics, and existing log list return. Which three-control spline behavior should legacy SQUAD represent? No signature replacement is approved. |
| QD07: Accuracy and degeneracy | No wiki tolerances; SLERP deliberately uses a nearby linear approximation | Should SLERP guarantee unit norm and at what tolerance? Define zero normalization, NaN/Inf, singular/near-singular inversion, antipodal interpolation and invalid frusta policy. Preserve LERP's existing component/extrapolation behavior unless explicitly changed. |
| QD08: Mixed axes and units | Vector −Z front is explicitly documented; Quaternion +Z forward and mixed degrees/radians are established source behavior | Retain both forward helpers and mixed units with documentation, or add new explicit alternatives? No sign/unit unification should occur by default. Angle conversion helper names can be corrected separately. |
| QD09: 2D rotation, transform, and shear | Wiki names methods but does not specify pivot, argument representation, or shear displacement source | Define 2D pivot versus axis contract, transform homogeneous behavior, and shear-plane mapping. Preserve row-vector storage/conventions. |
| QD10: Refraction and viewport | Wiki lists refraction; common viewport page is a placeholder | Confirm IOR ratio direction, unit-vector prerequisites and TIR output. Define what getViewPort is meant to compute instead of choosing a formula from its name. |
| QD11: Experimental scope | No experimental wiki; transport incomplete and irradiance map assumes square angular probes | Support rectangular probes or reject them clearly? Confirm raw float endianness and SH sign basis interoperability. Decide usable algorithm scope without inventing shadow-transport behavior. |
| QD12: Public numeric/export protocol | Scalar examples use float; ctypes export and sequence access are undocumented | Preserve c_matrix snapshot/export? Add integer/scalar reflected operators or indexing as additive APIs? Decide before widening protocol requirements. Do not add mandatory compiled dependencies. |

The default Matrix identity is **not** a proposal to change behavior: it is preserved, supported by contemporaneous source despite a wrong wiki output. Likewise, `i_*` returning self is compatible with “without returning a new object”; no change to None is recommended. `indetity` is a documentation typo, not an API requiring a compatibility alias.

### Phase 2C decision update

The user explicitly approved the refraction portion of QD10 before V02 implementation: IOR=n1/n2 (incident/transmitted indices), normalized incident and normal vectors with matching dimensions, normal opposing incidence and pointing into the incident medium, and preservation of the historical zero-Vector result for total internal reflection. No automatic normalization or normal flipping is added. See [current conventions](CONVENTIONS.md#phase-2c-current-angle-and-refraction-conventions). The viewport portion of QD10 and all other unresolved questions remain pending; approval of refraction does not resolve them.

## Recommended implementation order after review

1. **Small ordinary-math fixes with clear contracts:** equal-size equality/inequality (V01 component cases), inverse2 (M01), quaternion inverse (Q01). Preserve exact comparisons, return types, and component order; defer unsupported-dimension policy changes.
2. **Coupled division fix:** M02/M03/M07, updating inverse4's current transpose dependency together. Preserve legacy special methods, add Python 3 division, and check both inverse identities plus ctypes state.
3. **Named angle conversions and refraction:** C01, then V02 once IOR/prerequisite policy is accepted. Retain every other API's established units. Avoid blanket reflected-operator refactoring just to make local formulas work.
4. **Representation/state fixes after decisions:** G01–G04 after QD05; M04 translation consistency; Q07/Q08 and Q06 runtime issues with explicit ownership/type contracts. Keep commits scoped to one defect or tightly coupled cause.
5. **Transform/projection fixes after QD03/QD09:** V03/P01/P02 and M05/M06. Preserve row-vector composition, documented 3D project result, and explicit operator dimension restrictions. Verify noncommuting transforms.
6. **Ray and remaining quaternion work after QD02/QD04/QD06/QD07:** R01/R02/R03; Q02 matrix conversion; Q03/Q04 powers/log; Q05 SQUAD. Contract-dependent Q10/Q11 changes require explicit approval rather than being bundled with crash fixes.
7. **Numerical robustness after domain/error policy:** N01–N03, SLERP accuracy if approved, and conditioning/overflow behavior. Use meaningful bounds and pure-Python implementations; avoid a broad inversion redesign without measured need.
8. **Experimental algorithm repairs last:** E01/E02/E03 Bezier scalar/operand/Python-3 fixes, then independently reproduce masked E04 issues; E05 recurrence, E06 generation, and explicitly scoped E07/E08 work. Do not expand unfinished shadow transport implicitly.
9. **Only after correctness is established:** propose Phase 3 packaging, supported Python versions, CI, documentation/typing, and performance changes. Existing benchmarks stay as baseline; no PyPI publication or merge is included.

Each later review unit must show the pre-fix failure, retain unrelated xfails/questions, remove only the relevant strict defect markers, and include a compatibility note. Wiki examples with mistakes should receive documentation-only corrections separately from mathematical corrections. Existing API/import names and pure-Python design remain constraints throughout.
