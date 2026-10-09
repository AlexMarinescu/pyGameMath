# Tutorial verification and coverage

The nine subject guides contain self-contained executable Python blocks, each
with independently calculated values or a meaningful invariant/ownership check.
No companion script, new binary asset or external download is necessary. Longer
HDR rendering remains delegated to the existing source-tree example. This keeps
teaching code close to its explanation without duplicating algorithms in core.

## Run all checks

From the checkout with the existing audit requirements installed:

```sh
python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials
python -m pytest -q
python -m examples.hdr_sh.regenerate --output-dir /tmp/gem-tutorial-hdr
```

For installed examples, build the unchanged wheel/sdist in a clean copy, install
in fresh environments, and run the checker from outside the source tree. Supply
that environment's actual site-packages path with `--package-root`; the
[installation verification](../getting-started/verification.md) explains the
isolated `python -I` workflow. Do not rely on an editable/source-tree import as
proof of distribution execution.

The checker includes all tutorial pages, records block counts by page and executes
each block with an independent namespace. It also retains Phase 4A–4C link,
source-signature, API-coverage, alias-identity and protected-file checks. Assertions
validate unrounded numerical values; printed rounded values are explanatory.
Code examples alone do not establish wider API guarantees.

## Coverage and evidence

| Guide | Executed Python blocks | Independent checks |
|---|---|---|
| Vectors | 1 | 3–4–5 distance, facing angle, movement and ownership |
| Transformations | 1 | Explicit transformed coordinates, composition order, inverse recovery |
| Camera | 1 | Camera/clip/NDC/window values, near/far depth and zero W |
| Quaternions | 2 | Known rotations, noncommuting composition, sign equivalence and interpolation |
| Bezier | 1 | Polynomial points/derivative, ordered samples and ownership |
| Geometry | 2 | Analytic plane hit/classification and near/far picking positions |
| Lighting | 2 | Hand-derived L0/L1 RGB projection, active rotation, convolution and reflection |
| Interoperability | 1 | Offsets, snapshot replacement and retained-buffer lifetime |
| Numerical accuracy | 2 | Extreme norms, zero fallbacks, diagonal inverses and conditioning |

**13 new tutorial blocks** pass in source, isolated wheel and installed-sdist runs.
Combined Phase 4A–4D checks execute **42 blocks**, validate **418 internal links**,
check 268 source declarations/seven type aliases or reference buffers, and confirm
93 protected files unchanged. The complete suite reports **2,249 passed**, zero
failures, expected failures or skips. Clean wheel/sdist builds and prior Phase
4A–4C checks pass. Both HDR images regenerated with matching committed hashes in
the documented command and a separate isolated installed-core run.


[Machine-readable results](verification-results.json) record the per-guide block
coverage, command contexts, environment, full-suite counts and installed-wheel/
sdist checks. The HDR command is run headlessly and regenerated hashes compared
against existing reference images, alongside numerical reference tests. Image bytes
are corroborating same-environment evidence, not the sole cross-platform oracle.

Verification targets current repository code on CPython 3.12/Linux with six.
Python 2.7 and other platform/interpreter support remain unverified. Mathematical
implementations, existing tests/examples, benchmarks, runtime dependencies,
packaging/release metadata and license are unchanged. New tutorials remain
source-tree Markdown under the current manifest; they are not promised to ship
inside the installed package. [Open decisions](decisions.md) records assumptions
requiring separate integration or API review. Return to [index](index.md).
