# Tutorial boundaries and open decisions

These guides use established contracts. The local movement rule, scalar comparison,
Bezier derivative and ray-plane query are educational calculations, not new gem
functions. Their assumptions are stated beside the code; none changes public API
or chooses an application-wide policy.

| Subject | Tutorial treatment | Remaining review |
|---|---|---|
| Ray hit state | Separate query result; no automatic assignment to end; distance is stored constructor state | Hit validity, transform ownership and future intersection-result contracts remain unresolved |
| Near-parallel queries | Exact-zero denominator classification for finite demo data | Application tolerance/range policy requires explicit units and accuracy choices |
| Camera/window integration | Lower-left viewport and OpenGL NDC depth; explicit homogeneous coordinates | External top-left/depth conventions, clipping and zero-W sentinel interpretation require application handling |
| Bezier motion | Analytic derivative locally; parameter-time traversal; sampled chord-length discussion | No constant-speed traversal or automatic orientation frame is claimed |
| Interpolation controls | Unit endpoints/controls supplied explicitly to squad4 | No control generation or new general-domain validation is introduced |
| Nonfinite/conditioning | Narrow finite-domain examples and existing fallbacks/errors | Wider shape/scalar/nonfinite and numerical-range policies remain separate |
| Foreign matrix use | Verify memory layout/lifetime; no OpenGL call | Binding, shader, transpose, endian and platform behavior require actual external integration testing |
| Interpreter/release support | Executed CPython 3.12/Linux only; current repository package | No Python 2.7 claim, new support matrix, packaging change or release designation |

Use [API decision notes](../api/decisions.md) and the
[central development decisions](../development/roadmap.md#open-decisions) for
maintainer review rather than treating these examples as a new compatibility
policy. Future noise, SDF, frustum, voxel, area-light and general collision APIs
remain planned, not demonstrated or importable. Existing experimental shims stay
intact. See [verification](verification.md), [tutorial index](index.md) and
[canonical roadmap](../../ROADMAP.md).
