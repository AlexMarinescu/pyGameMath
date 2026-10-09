# Development roadmap entry

The [root ROADMAP.md](../../ROADMAP.md) is the single authoritative roadmap.
This page links to it rather than reproducing milestone scope or status.
Use the [architecture charter](../architecture/philosophy.md),
[API inventory](../architecture/api-inventory.md),
[conventions](../architecture/conventions.md) and
[compatibility policy](../architecture/compatibility.md) for design constraints.

## Open decisions

These recommendations need maintainer review before implementation or release;
they do not change existing APIs or block this documentation-only review.

| Decision | Recommended next step | Compatibility consequence |
|---|---|---|
| Release interpreter/platform matrix and Python 2.7 legacy target | Assess the concrete source/tool blockers, choose modern minimum/test matrix, and determine how the 2.7 target can be engineered/tested | Do not advertise unverified interpreters or remove six in advance |
| Experimental alias removal boundary/window | Retain current shims through architecture/release review; explicitly announce a future approved boundary before deletion | Complete package elimination would break documented imports and some serialized references |
| 1.0 public stability/versioning promise | Review raw kernels, historical helpers and mutable fields separately from well-documented contracts; confirm proposed semantic-version policy | Passing tests alone cannot freeze every incidental public attribute |
| Wider shape/scalar/error policy | Keep current narrow contracts; separately design mismatch, nonfinite, degeneracy and numeric-protocol rules | Uniform validation/reflected operators would change current exceptions/accepted operands |
| Future Ray intersection results | Define hit validity, distance versus parameter and ownership before adding queries; keep legacy end state untouched | Zero end cannot identify an unset hit; silently reinterpreting it would break state semantics |
| Future module organization and documentation/release tooling | Review cohesive conceptual areas and navigation, then select framework/Wiki migration/build work in later phases | Do not rename current imports or treat provisional API sketches as shipped |

Previously settled row-vector, quaternion, refraction, plane, sampling and legacy
helper conventions are not reopened here. Their historical decision trail is in
[PHASE2-DECISIONS.md](../../audit/PHASE2-DECISIONS.md); current summaries describe
the implemented outcomes. See [verification](verification.md) for the documentation
review evidence and unresolved scope boundaries.
