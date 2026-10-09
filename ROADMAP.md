# pyGameMath roadmap

This is the **single authoritative development roadmap** for gem. Architecture
pages describe principles and current behavior; historical audit plans record
completed work and do not supersede this roadmap. Future features and API sketches
below are proposals, not shipped interfaces or guaranteed release commitments.
Assignments beyond 1.0 are provisional, with no release dates.

gem is a pure-Python graphics and computational mathematics library. Keep existing
imports, signatures, row-vector transforms, quaternion ordering and ownership
contracts. No NumPy, Cython, native mathematical backend or mandatory numerical
framework is planned. Standard-library integrations and current six support remain.
See the [architecture charter](docs/architecture/philosophy.md),
[current API inventory](docs/architecture/api-inventory.md),
[conventions](docs/architecture/conventions.md) and
[compatibility evidence](docs/architecture/compatibility.md).

## Present foundation

Implemented core: Vector/Matrix/Quaternion algebra, transforms/projection,
Plane/Ray geometry, quadratic/cubic Bezier evaluation and adaptive sampling,
Legendre functions and real SH projection/convolution/analytical L2 rotation.
`Vector.barycentric` already exists; future triangle work will broaden it into
coherent geometry/query contracts. The headless HDR/SH example is a reference
workflow, not a renderer. Ray intersections, robust predicates, noise, spatial
trees and the other future systems below have not shipped.

Confirmed audit correctness fixes and Phase 3 performance work are complete
within their stated domains. Phase 3E records 2,249 passing tests, zero expected
or unexpected failures, with performance safeguards separate from mathematical
assertions. Passing tests do not settle every malformed/extreme input or confer
a blanket stable API designation. Transitional experimental aliases remain;
unfinished E07 shadow transport is removed with no replacement. The current
package declaration is v0.1.12, not 1.0.

## Version 1.0 — Stable foundation

**Purpose/status:** prepare a maintainable, accurately documented foundation.
Correctness and performance phases are completed; documentation modernization
and release engineering are in progress/planned. 1.0 is not released.

**Scope:** consolidate current correctness contracts, performance monitoring,
complete API reference, installation/quick start, tutorials and runnable examples;
verify builds/distributions, Python compatibility and public release preparation.
Review existing public interfaces without unnecessary breaking changes. No new
mathematical feature is required merely to call the foundation 1.0.

**Prerequisites:** architecture/inventory review, approved stability/compatibility
policy, chosen interpreter/platform matrix, shim removal/retention decision and
resolved release-blocking documentation gaps. **API direction:** preserve current
modules; proposed reference organization is not a renamed runtime hierarchy.

**Numerics/performance:** retain independent known answers, extreme-scale and
ownership regressions, binary32 export contracts and approximation limits.
Use the versioned performance suite with manual review; never sacrifice numerical
correctness to recover historical timings. Do not infer cross-machine speedups.

**Python:** test chosen modern interpreters on real environments. Assess the
deliberate Python 2.7 legacy target, source blockers and compatible tools before
making support claims. Retain six until that decision is engineered and verified.

**Applications/docs:** small transform/geometry examples, curve sampling and the
existing headless HDR lighting pipeline. Provide reference pages for every intended
public interface, units/ownership/error domains, practical tutorials, migration
guidance and reproducible benchmark/build commands.

**Acceptance:** no unexpected failures; reviewed public API/domain inventory;
working documented examples; clean wheel/sdist and outside-source imports on the
declared matrix; verified metadata/URLs/licenses; maintained functional CI; reviewed
compatibility notes and release artifacts. Packaging/documentation framework
choices and any publishing require subsequent phases; this roadmap changes neither.

## Future 1.x — Geometry foundations

**Purpose/status:** coherent primitive geometry and query building blocks;
proposed, unimplemented expansion beyond existing Plane/Ray/barycentric helpers.

**Scope:** triangles and barycentrics, AABB/OBB/sphere/capsule bounds, ray/primitive
intersections, frustum extraction/culling, robust geometric predicates, coordinate
frames and tangent-space mathematics.

**Prerequisites:** Vector/Matrix conventions, hit-record and ray distance/parameter
contracts, boundary/inclusion policies, predicate accuracy model. **Provisional API
sketches:** `triangle_barycentric`, `intersect_ray_triangle`, `AABB`, `orientation2d`,
`tangent_frame`; these names do not exist as new public APIs and require review.

**Numerics:** degenerate triangles, parallel rays, signed/zero extents, closed/open
boundaries, nearly coplanar predicates and valid near-singular transforms. Separate
exact/filtered predicates from heuristic tolerances; specify handedness and winding.
**Performance:** benchmark primitive queries and allocation; use simple operations
as building blocks for spatial acceleration rather than requiring it for every query.
**Python:** pure-Python small-object routines; verify the chosen support matrix and
numeric type behavior before using newer interpreter-only helpers.

**Applications/docs:** picking, conservative culling, collision-query reference math
and procedural mesh frames, not a physics engine. Document diagrams, query result
ownership, parameter units, degeneracy policies and independently worked examples.
**Acceptance:** independent incidence/containment/barycentric tests, boundary and
degenerate cases, orientation invariants, transform covariance, caller preservation
and measured representative throughput without weakening correctness.

## Future 1.x — Procedural mathematics

**Purpose/status:** deterministic procedural fields and useful sample distributions;
proposed, not implemented.

**Scope:** gradient/simplex-style and value noise, Worley/cellular noise, fBm,
ridged fractals/domain warping, Poisson disk and low-discrepancy sampling, additional
spline families. Existing Bezier evaluation remains canonical and is not replaced.

**Prerequisites:** seed/RNG ownership, coordinate/dimension contracts, noise period
and output-range definitions, curve continuity/parameter rules. **Provisional API
sketches:** `gradient_noise`, `cellular_noise`, `fbm`, `poisson_disk`, `halton`,
`evaluate_bspline`; concepts only, with no final signatures or module names.

**Numerics:** interpolation derivatives, negative coordinates, boundary periodicity,
octave amplitude/frequency, distribution bias and minimum-distance checks. Avoid
claiming exact uniformity from empirical samples. **Performance:** quantify lattice
lookup, dimension/octave growth and Poisson neighbour-query cost; bound attempts
and memory. **Python:** explicit seeds/state rather than incidental global RNG where
a new contract permits it; distinguish reproducibility within an implementation
from identical random output across interpreter versions.

**Applications/docs:** terrain textures, procedural displacement, sampling and
animation curves. Explain range/seed/period semantics, composition examples,
continuity and distributions. **Acceptance:** independent lattice/known-spline
answers, reproducibility/ownership tests, continuity or derivative checks,
statistical diagnostics without brittle random assertions, termination and scaling
measurements for stated dimensions.

## Future 1.x — Spatial and implicit mathematics

**Purpose/status:** reusable spatial query and implicit-geometry math; proposed,
not implemented.

**Scope:** SDF primitives/operations/gradients/transforms, sphere tracing, voxel
coordinate transforms, 3D DDA, trilinear sampling, Morton encoding, octrees,
spatial hashing and BVH. Algorithms should expose clear data/query contracts
without an engine scene graph or mandatory array framework.

**Prerequisites:** bounds/intersections, frames, sampling, distance-sign and ray
parameter conventions, finite traversal/error policies. **Provisional API sketches:**
`sdf_sphere`, `sdf_union`, `sphere_trace`, `voxel_dda`, `trilinear`, `morton_encode`,
`BVH`; unimplemented names requiring architecture review.

**Numerics:** signed-distance versus approximate fields, nonuniform scale and
Lipschitz bounds, gradient degeneracy, voxel-face ties, zero ray components,
integer bit range and empty bounds. **Performance:** finite traversal limits;
construction/query tradeoffs, coherent memory ownership, pathological tree depth,
cell load and sample count. Pure Python targets small/reference workloads;
large voxel/mesh processing is not promised to match native array libraries.
**Python:** document integer semantics and serialization, avoid platform-sized
assumptions and test the approved interpreter matrix.

**Applications/docs:** sparse procedural worlds, implicit modelling, picking and
spatial search references. Explain grid/index coordinates, distance validity,
traversal tie rules and scaling limits. **Acceptance:** analytic SDF known answers,
finite-difference gradient references, exhaustive small-grid traversal/Morton cases,
brute-force spatial-query comparisons, deterministic termination and memory/time
scaling on adversarial inputs.

## Future 1.x — Lighting and volumetric mathematics

**Purpose/status:** radiometric and integration building blocks; expansion proposed.
Existing SH/environment reference functionality is already implemented and remains
distinct from the future features below.

**Scope:** radiometric/photometric helpers, BRDF evaluation/sampling, Fresnel and
microfacet math, rectangular/disk/spherical area-light geometry, solid angles,
Monte Carlo/numerical integration, cone geometry and voxel-cone-tracing reference
math, volume integration/transmittance and SH-based volumetric lighting utilities.

**Prerequisites:** unit/frame conventions, geometric sampling, probability measures,
integration error policy and explicit visibility/grid query interfaces. **Provisional
API sketches:** `fresnel_dielectric`, `ggx_brdf`, `sample_brdf`, `solid_angle_disk`,
`integrate_transmittance`; proposed functions, not shipped APIs or a renderer.

**Numerics:** radiance versus irradiance, steradians, projected-area PDFs, energy
conservation/reciprocity, grazing angles, delta limits, positive extinction and
underflow. SH convolution must occur once and sign/basis conversion stay explicit.
Cone tracing is a reference approximation, not a guaranteed visibility solution.
**Performance:** integration/sample counts, variance, quadrature convergence and
volume step growth; avoid hiding runtime-heavy transport inside simple basis APIs.
**Python:** independent RGB channels, linear working space, core free of image/GPU
dependencies; optional rendering adapters stay outside gem.

**Applications/docs:** material/reference shading, environment probes and homogeneous
volume attenuation. Document physical units, PDFs, linear/display color, integration
limits and analytic comparisons. **Acceptance:** independent solid-angle/constant
integrals, BRDF normalization/energy tests, PDF correctness, exponential homogeneous
transmittance, SH/channel/frame invariants, bounded integration and calibrated
runtime/error measurements. No complete GI probe-volume or visibility tracer is implied.

## Advanced computational geometry

**Purpose/status:** deeper topology and surface algorithms; research-backed future
milestone, with version assignment undecided.

**Scope:** Delaunay triangulation, bounded Voronoi diagrams, Lloyd relaxation,
half-edge connectivity, surface differential geometry, advanced topology and
subdivision. Package organization is conceptual until reviewed.

**Prerequisites:** robust orientation/in-circle predicates, clipping/bounds, spatial
queries, explicit winding/manifold/topological invariants and ownership model.
**Provisional API sketches:** `delaunay2d`, `bounded_voronoi`, `lloyd_relax`,
`HalfEdgeMesh`; none is implemented or approved as a final interface.

**Numerics:** duplicate/coincident points, cocircular ties, exact predicate escalation,
bounded/unbounded cells, nonmanifold input, curvature discretization and mesh
degeneracy. **Performance:** asymptotic versus pathological cases, Python object
memory, incremental updates and bounded iteration; establish realistic input-size
limits before claiming production-scale mesh processing. **Python:** explicit integer/
index representation and portable serialization; no compiled geometry dependency.

**Applications/docs:** planar meshing, partitioning, procedural topology and surface
analysis. Provide algorithm derivations, topology invariants, failure domains,
visual reference examples and complexity/memory guidance. **Acceptance:** independent
small known triangulations, empty-circle/tessellation properties, Euler/connectivity
invariants under supported topology, controlled relaxation convergence, degeneracy
tests and measured size scaling. Reject unsupported input explicitly once a policy
is chosen, rather than concealing it as a valid mesh.

## Research & experimental

Speculative work has no release commitment: higher-band analytical SH rotation,
quaternion exponential/control generation, improved high-order basis numerics,
adaptive integration, advanced implicit/volume approximations and visibility-aware
lighting interfaces. First define mathematical domains, independent references,
ownership, performance feasibility and engine boundaries. Retired E07 scaffolding
is historical evidence, not a ready transport implementation.

Research belongs in clearly labeled external prototypes/examples or design notes
until validated and approved for core. It must not recreate unfinished mathematics
inside transitional `gem.experimental` shims. Retaining compatibility imports does
not reserve that namespace for new implementations.

## Review and progress policy

Milestones depend on foundations, not dates: algebra/contracts support geometry
and sampling; those support spatial/implicit/lighting systems; advanced topology
depends on robust predicates and spatial queries. Future releases can select a
small coherent subset rather than deliver every 1.x proposal together.

Changes to proposed scope, API names, support policy or compatibility boundaries
require reviewed roadmap updates. Record completion only with implementation,
independent tests, examples, numeric limits and build evidence. The
[development page](docs/development/roadmap.md) lists open decisions without
duplicating this roadmap; the [documentation index](docs/README.md) provides
navigation. Architecture review does not begin feature implementation or authorize
release publication.
