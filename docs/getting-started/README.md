# Getting started

Follow this sequence with current repository code:

1. [Install in a fresh environment](installation.md) and confirm which package is imported.
2. [Run the quick start](quick-start.md): vectors, matrix transforms/inversion, quaternion rotation, Bezier evaluation, SH and ctypes.
3. Follow [practical tutorial learning paths](../tutorials/index.md), or explore the [headless HDR/SH reference](../../examples/hdr_sh/README.md), [quaternion guide](../QUATERNIONS.md) or [ray guide](../RAYS.md).
4. Check [conventions](../architecture/conventions.md), [ownership/compatibility](../architecture/compatibility.md) and the [API reference](../api/index.md) before broader integration.

These examples need no optional renderer, GPU or numerical framework. six remains
the runtime dependency. Verified reference environment: CPython 3.12.14, Linux
x86_64; other interpreter/platform combinations are unverified here. Cross-version
syntax in an example does not establish Python 2.7 support.

The [root roadmap](../../ROADMAP.md) describes planned expansion separately from
implemented APIs. Return to the [documentation index](../README.md) for navigation.
See [example/build verification](verification.md) for the exact evidence and limits.
