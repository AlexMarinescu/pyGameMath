# Architecture charter

gem is a lightweight, portable, pure-Python mathematics library for computer
graphics, computational geometry, game development, simulation, procedural
generation and numerical research. Its foundation is small mathematical objects
and independently testable algorithms, with reference examples showing how to
use them. The distribution and import name remain `gem`; the project is
pyGameMath. This charter sets direction, not a declaration that proposed features
or a 1.0 release already exist.

## Principles

- Implement mathematical algorithms in Python. No NumPy runtime dependency,
  Cython, compiled mathematical extensions, native SIMD backend or mandatory
  numerical framework.
- Use the standard library, including `math` and `ctypes`. Retain `six` for now;
  minimize third-party runtime dependencies over time.
- Put correctness, numerical stability, predictable ownership and lightweight
  imports ahead of unverified speed claims. Use independent mathematical
  references and repeatable performance measurements.
- Preserve established imports, signatures and conventions. Breaking changes
  need an explicit compatibility decision, migration guidance and release boundary.
- Keep engine-specific image loading, shaders and rendering demonstrations outside
  core. Core lighting mathematics must not require an OpenGL context or GPU.
- Verify interpreter/platform support explicitly. Python 2.7 is unsupported;
  the tested release matrix is CPython 3.10–3.14 on Linux x86_64.

gem does not replace NumPy's large-array processing. It is not a game engine,
GPU renderer, shader compiler, scene graph, ECS, asset manager, general-purpose
physics engine or dependency-heavy scientific computing distribution. Future
sampling, BRDF and volume algorithms are mathematical building blocks rather
than a complete rendering or simulation framework.

## Present architecture

| Layer | Existing modules | Dependencies within gem |
|---|---|---|
| Component arithmetic and interoperability | `gem.vector`, `gem.common` | Neither imports another core module |
| Matrices and transforms | `gem.matrix` | vector, common |
| Quaternion algebra and rotations | `gem.quaternion` | vector, matrix, common |
| Supporting geometry | `gem.plane`, `gem.ray` | plane: vector; ray: vector, quaternion |
| Curves | `gem.bezier` | vector |
| Polynomial basis | `gem.legendre` | None |
| Lighting basis and integration | `gem.spherical_harmonics` | vector, legendre; quaternion loaded inside coefficient rotation |
| Import compatibility | `gem.experimental` shims | Reexports the canonical core implementations |

`gem` itself exposes version metadata, not a facade reexporting mathematical classes.
Use explicit module imports. `Vector(size, data)` and `Matrix(size, data)` are
the existing constructors; there are no separate `Vector3` or `Matrix4` classes.
The experimental directory still exists solely for compatibility. Unfinished
shadow-transport code has been retired without a core replacement.

The [API inventory](api-inventory.md) distinguishes documented contracts from
historical conveniences and unreviewed edge domains. Passing tests is evidence
for the cases exercised, not a blanket stability certificate. The current
package version and metadata are not changed by this charter.

## Proposed mathematical foundation

The following are conceptual areas, **not existing modules or public APIs**.
Existing module names remain canonical. Prefer adding cohesive modules after
review rather than reorganizing working imports or splitting every algorithm
into a separate package.

| Conceptual area | Mathematical scope | Prerequisites |
|---|---|---|
| Algebra and numerical utilities | Small vectors/matrices, rotations, numerical primitives | Existing core contracts; verified numeric domains |
| Geometry and robust predicates | Triangles, barycentrics, intersections, orientation tests | Algebra; exact/filtered predicate policy |
| Curves, surfaces and transformations | Additional splines, frames, derivatives, surface parameterization | Bezier; algebra; continuity and parameter contracts |
| Spatial mathematics and acceleration | Bounds, frusta, BVH, octrees, spatial hashing | Geometry, bounds and ray-query contracts |
| Procedural noise and sampling | Noise families, fractals, low-discrepancy and Poisson sampling | RNG ownership, dimensions, distribution/reference tests |
| Implicit geometry | SDF primitives, composition, gradients and sphere tracing | Geometry, transforms and distance/Lipschitz conventions |
| Volumetric mathematics | Voxel coordinates, DDA, interpolation, volume integration | Frames, bounds, sampling and spatial queries |
| Lighting and radiometry | BRDFs, Fresnel, solid angles, SH and integration | Unit directions, probability measures, geometry and sampling |
| Advanced geometry and topology | Triangulation, Voronoi, half-edge connectivity, differential geometry | Robust predicates, spatial queries and explicit topology invariants |

The dependency direction should flow from algebra and numeric contracts toward
geometry/sampling, then spatial, implicit and lighting systems. Visibility-dependent
transport needs an explicitly supplied scene/query abstraction; it must not be
smuggled into the SH basis module. Mesh connectivity must not become a dependency
of simple Vector imports. Future module names and interfaces require review before
implementation; this document does not authorize a new namespace hierarchy.

## Engineering boundaries

Use independent known answers, invariants and input-preservation tests for each
algorithm. State units, dimensions, degeneracy policies and approximation limits
before implementing ambiguous behavior. Keep exact predicates distinct from
epsilon heuristics. Distinguish integration error, truncation error, ill-conditioning
and binary64 representability.

Measure public wrapper costs as well as raw kernels. Preserve necessary object
allocation and float32 ctypes snapshots; removing them would change compatibility.
The [performance policy](../../benchmarks/REGRESSION_POLICY.md) separates functional
CI from calibrated, manually reviewed timing evidence. It does not promise one
machine's speedup on every interpreter.

The [canonical roadmap](../../ROADMAP.md) sets priorities. The
[compatibility policy](compatibility.md) and [conventions](conventions.md) constrain
future work. Outstanding choices are listed on the [development page](../development/roadmap.md).
