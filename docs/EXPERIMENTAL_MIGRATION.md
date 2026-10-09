# Experimental import migration

Validated mathematics has one implementation in the supported core. New code
should use core imports. The experimental directory still exists: it contains
only minimal compatibility reexports and a package marker. Complete removal
would break previously documented imports, so the remaining paths await an
explicit release-boundary and compatibility-window decision. No runtime
warnings or duplicated algorithms are introduced.

| Historical import | Canonical replacement | Compatibility status | Removal status / migration action |
| --- | --- | --- | --- |
| `gem.experimental.bezier` | `gem.bezier` | Identical functions and `BezierPath` class | Retained shim; update imports |
| `gem.experimental._bezier_legacy.BezierPath` | `gem.bezier.BezierPath` | Identical class; private path remains importable | Retained shim; update imports |
| `gem.experimental.legendre.Legendre` | `gem.legendre.Legendre` | Identical class | Retained shim; update imports |
| `gem.experimental.sph` | `gem.spherical_harmonics` | Identical `Factorial`, `K`, `SPH`; incidental `Legendre` retained | Retained shim; prefer `gem.legendre` for `Legendre` |
| `gem.experimental.sph_sample` | `gem.spherical_harmonics` | Identical `SPHSample`, `GenerateSamples` | Retained shim; update imports |
| `gem.experimental.sph_irradiance_map` | `gem.spherical_harmonics` | Identical `SPH_IrradianceMapCoeff` | Retained shim; update imports; legacy basis still requires explicit conversion |
| `gem.experimental.sph_object` | None | Unfinished transport retired | Removed; `SPHVertex`, `SPHObject`, `GenereateCoeffs` imports fail |
| `gem.experimental` | Core modules above | Package marker retained for shims | Not fully removed; no new experimental implementations |

E07 never provided implemented visibility transport. Its removal is an import
break for consumers of its containers or unfinished routine; it is not a
mathematical correction or a successful replacement of shadow transport.
Repository consumers consisted of its expected-failure and diagnostic tests.
External usage is unknown. Historical source remains in git; the original test
and nine diagnostic cases are archived under `audit/retired-tests`, with fault
results in `audit/PHASE2F6.md`. They are deliberately retired from collection,
not skipped or relabeled as fixed. New tests assert this module's absence.

The retained paths continue to ship in both wheels and source distributions.
Changing an import to its core equivalent preserves signatures, storage,
ownership, basis signs, angle units and numerical behavior. Updating the import
does not convert historical probe coefficients into canonical SH: use
`legacy_to_canonical` explicitly. Existing serialized references to retired
transport classes may require application-specific migration; no serialization
shim is provided.

Future package deletion must choose a release boundary, announce all removed
paths, migrate compatibility tests and change packaging together. Scene
visibility, acceleration and transport remain outside core math. No new
transport API or Phase 3 optimization is included.
