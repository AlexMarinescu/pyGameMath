# Documentation navigation

Start with the [project charter](architecture/philosophy.md),
[current API inventory](architecture/api-inventory.md) and
[canonical roadmap](../ROADMAP.md). The [website homepage](index.md) and [build instructions](development/website.md)
provide the static-site entry. The source pages remain canonical.

| Section | Current entry point | Future documentation work |
|---|---|---|
| Installation and quick start | [README](../README.md), [getting started](getting-started/README.md) | Release-specific installation and broader interpreter verification |
| API reference | [Topic reference](api/index.md), [audited inventory](architecture/api-inventory.md) | Continued examples and reviewed domain clarifications |
| Practical tutorials | [Learning paths](tutorials/index.md), [Quaternion guide](QUATERNIONS.md), [Ray guide](RAYS.md), [SH guide](SPHERICAL_HARMONICS.md) | Further worked applications and reviewed integration examples |
| Examples/use cases | [Headless HDR/SH workflow](../examples/hdr_sh/README.md), [launcher](../launcher.py) | More executable, independently verified small examples |
| Mathematical conventions | [Current conventions](architecture/conventions.md) | Keep reference examples synchronized with tests |
| Benchmarks/performance | [Benchmark guide](../benchmarks/README.md), [monitoring policy](../benchmarks/REGRESSION_POLICY.md) | Maintained calibrated baselines and review guidance |
| Compatibility/migration | [Policy/evidence](architecture/compatibility.md), [import migration](EXPERIMENTAL_MIGRATION.md), [Vector/viewport](VECTOR_VIEWPORT_CONTRACTS.md) | Release support matrix and compatibility announcements |
| Contributing | [Verification workflow](development/verification.md) | Contribution/review process, code style and release-specific tooling guide |
| Roadmap/development status | [Development entry](development/roadmap.md) | Reviewed progress updates in the root roadmap |

Future sections in the table are proposals, not empty published APIs or nonexistent
linked pages. The [selected documentation stack](development/documentation-stack.md) and
[non-destructive Wiki migration](development/wiki.md) are prepared for review. Existing audit reports remain chronological evidence rather
than a second public roadmap. The complete architecture documentation is not yet
included by current distribution manifests; packaging changes are separate work.

Explore the [headless graphics gallery](examples/index.md) for reproducible visual examples.
