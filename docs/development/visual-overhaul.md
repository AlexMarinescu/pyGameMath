# Documentation design and verification

The website uses a local wordmark, yellow **py** on a dark header, high-contrast
light/dark palettes and system fonts. Desktop tabs group Getting Started,
Tutorials, API Reference, Examples, Conventions and Development. Mobile uses the
same content in the theme's keyboard-accessible drawer. Existing page paths and
API declaration tables are retained.

## Learning and mathematical presentation

Begin with the [quick start](../getting-started/quick-start.md), then follow one
of the [nine tutorial paths](../tutorials/index.md). Exact signatures live in the
[API reference](../api/index.md); shared contracts live in
[conventions](../architecture/conventions.md).

The [notation guide](../architecture/notation.md) connects row-vector composition,
quaternion interpolation, Bezier subdivision, Legendre functions and spherical
harmonics. Arithmatex and the build hook convert TeX to native MathML, using the
pinned documentation stack with no external scripts, fonts or CDN requests.
Equations can scroll independently on narrow screens.

The [visual gallery](../examples/index.md) retains its verified generated assets.
The two HDR spheres share camera, albedo, exposure and display settings; only
analytical lighting rotation changes. The new Legendre plot can be regenerated:

```sh
python tools/site/legendre_diagram.py
```

It evaluates the supported core module and verifies analytical reference values;
it does not add an algorithm or mathematical API.

## Validation workflow

Run the commands in the [website build guide](website.md). The documentation
checker now freezes runtime, existing tests, examples, benchmarks, metadata and
dependencies at merged PR #59. It continues checking all 268 declarations and
runtime signatures, and additionally preserves the syntax trees of all 42
baseline executable examples and all nine tutorials. Comment-only clarification
is allowed; removing or changing their executable calculations is rejected.
The landing page adds one independently executable introductory example.

Browser verification checks light/dark contrast, equation dimensions, rendered
images, search, mobile navigation, keyboard interaction and hosting under
`/pyGameMath/`. It records screenshots and fails on external resources,
JavaScript errors or missing requests. These checks supplement the emitted HTML
link/anchor checks; they are not a complete accessibility certification.

The [Phase 4F-C report](../../audit/PHASE4FC-DOCUMENTATION.md) records executed
results, commands and remaining publication limitations. Production mathematics,
public signatures, package metadata and core dependencies remain unchanged.
