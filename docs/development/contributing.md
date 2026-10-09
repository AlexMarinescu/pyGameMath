# Contributing

Start with the [architecture](../architecture/philosophy.md),
[conventions](../architecture/conventions.md), [compatibility policy](../architecture/compatibility.md)
and [root roadmap](../../ROADMAP.md). Keep changes focused and explain the concrete
problem, numerical behavior and verification in a pull request to master.

## Mathematical changes

Use independent known-answer and property checks, including ownership and extreme
values where relevant. Preserve the pure-Python core and established row-vector,
quaternion and public API contracts. Resolve ambiguous behavior before changing it;
do not broaden invalid-input policy as an incidental cleanup. Benchmarks inform
review but do not replace mathematical correctness tests.

```sh
python -m pip install -e .
python -m pip install -r requirements-audit.txt
python -m pytest -q
```

These commands run from a source checkout. Use clean installed-wheel checks for
packaging/import changes. See the [verification workflow](verification.md) and
[performance policy](../../benchmarks/REGRESSION_POLICY.md).

## Documentation changes

Edit canonical Markdown, not generated site files. Keep executable examples,
source links, formulas and generated assets consistent. Run the
[documentation build and checks](website.md); use existing gallery regeneration
rather than manually editing diagrams. Concise engineering notes should describe
what changed, why, compatibility consequences and the tests performed.

For defects or proposed features, open a repository issue with reproduction inputs,
expected mathematical result, actual result and relevant environment information.
Keep API proposals distinct from shipped interfaces. Publication and release
settings require a separate maintainer review; merging a documentation change
must not silently publish a website or change the release version.
