# Compatibility imports and public-interface status

The `gem` initializer is empty: import classes/functions from their modules.
There are no root Vector/Matrix aliases, separate Vector2/3/4 classes or a
root version API. `gem.experimental` is an empty retained package marker.

| Status | Meaning here |
|---|---|
| Canonical core | Implemented algorithms under gem.vector/common/matrix/quaternion/plane/ray/bezier/legendre/spherical_harmonics |
| Legacy but retained | Preserved helpers and field spellings described on the topic pages; not blanket 1.0 stability |
| Transitional import | Thin reexport of the canonical object; no duplicate implementation |
| Retired | Removed unfinished E07; no core transport replacement |
| Experimental algorithm | None currently retained in the experimental directory; only compatibility shims remain |
| Private | Leading-underscore helpers and implementation imports/cache; not supported mathematical entry points |

No runtime deprecation, removal date or new support promise is declared here.
See [migration policy](../EXPERIMENTAL_MIGRATION.md) and
[compatibility assessment](../architecture/compatibility.md).

## Retained paths and exact exports

| Historical module | Exported objects | Canonical module |
|---|---|---|
| gem.experimental.bezier | cubicBezierPoint, quadraticBezierPoint, BezierPath | [gem.bezier](bezier.md) |
| gem.experimental._bezier_legacy | BezierPath (private historical path retained) | [gem.bezier](bezier.md) |
| gem.experimental.legendre | Legendre | [gem.legendre](legendre.md) |
| gem.experimental.sph | Factorial, K, SPH, Legendre | [gem.spherical_harmonics](spherical-harmonics.md); prefer gem.legendre for Legendre |
| gem.experimental.sph_sample | SPHSample, GenerateSamples | [gem.spherical_harmonics](spherical-harmonics.md) |
| gem.experimental.sph_irradiance_map | SPH_IrradianceMapCoeff | [gem.spherical_harmonics](spherical-harmonics.md) |

These are identical objects with identical signatures, basis conventions and
ownership, not adapters. `__all__` in experimental.bezier names its three exports;
experimental.legendre names Legendre. Other shims use direct imports without a
new `__all__`; incidental attributes are not an additional mathematical contract.
Canonical `gem.spherical_harmonics.Legendre` is the same imported core class.
Changing a probe-class import does not change legacy coefficients into canonical
SH; use `legacy_to_canonical` explicitly.

```python
from gem.bezier import BezierPath, cubicBezierPoint
from gem.legendre import Legendre
from gem.spherical_harmonics import SPH, SPHSample
from gem.experimental.bezier import BezierPath as OldPath, cubicBezierPoint as old_cubic
from gem.experimental.legendre import Legendre as OldLegendre
from gem.experimental.sph import SPH as old_sph
from gem.experimental.sph_sample import SPHSample as OldSample

assert OldPath is BezierPath and old_cubic is cubicBezierPoint
assert OldLegendre is Legendre and old_sph is SPH and OldSample is SPHSample
```

## Retired and non-core interfaces

`gem.experimental.sph_object`, including SPHVertex, SPHObject and misspelled
GenereateCoeffs, was removed. Import fails; no compatibility transport container
or serialization shim exists. [Archived audit evidence](../../audit/PHASE2F6.md)
and [packaging/retirement tests](../../tests/test_core_packaging.py) distinguish
retired diagnostics from fixed mathematics. External consumers are unknown.

`launcher.py` is a source-tree print demonstration; its ignored `test` argument
is not an assertion runner. `examples/hdr_sh`, benchmark scripts and documentation
tools are development/example code, not installed core mathematical APIs. No
future roadmap module is implied to be importable. The historical Wiki is audit
evidence and may describe corrected defects; it is not the current API reference.

The removal window for retained aliases, broader public stability and invalid-input
policy remain [maintainer decisions](decisions.md). See [index](index.md).
