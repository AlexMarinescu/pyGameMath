# Vector comparison, clamp ownership and the legacy viewport helper

Vector equality is exact component equality for equal declared dimensions.
Different dimensions compare unequal in both operand orders; two empty Vectors
compare equal and empty/nonempty compare unequal. Inequality is the logical
complement of equality for supported Vectors. Both special methods return
`NotImplemented` for other operand types. No approximate comparison is added.
Malformed component storage remains outside this cleanup's domain.

`clamp(size, value, minS, maxS)` and `Vector.clamp(size, value, minS, maxS)`
return fresh Vector/list storage. Value and bound lists remain unchanged,
including when multiple Vectors share them. The method receiver is preserved.
`Vector.i_clamp(size, value, minS, maxS)` still returns and changes the receiver,
but replaces its component list without changing a separate caller list or
another Vector sharing the old list. Clamping retains max-bound comparison
followed by min-bound comparison; no new validation of reversed bounds is added.
Aliased bounds are read from their original data rather than incidentally changed
by clamping the value list.

```python
from gem.vector import Vector, clamp

values = [-2, 2, 10]
a = Vector(3, values)
b = Vector(3, values)
result = clamp(3, values, [0]*3, [5]*3)  # [0, 2, 5], fresh storage
assert values == [-2, 2, 10]
a.i_clamp(3, values, [0]*3, [5]*3)
assert a.vector == [0, 2, 5]
assert b.vector == values == [-2, 2, 10]
```

## getViewPort

`getViewPort(coords, width, height)` retains the historical source formula.
For a finite nonzero Vector with at least two components, it normalizes the
**whole Vector**, including Z/W where supplied, using the existing stable
normalization routine. If L is its mathematical norm:

```
x_result = (coords.x / L + 1) * width / 2 + coords.x
y_result = (coords.y / L + 1) * height / 2 + coords.y
result = [x_result, y_result, width, height]
```

Vector2/3/4 are covered by known-answer regressions. Original coordinates are
both normalized direction components and additive XY offsets; this unusual
formula is preserved, not replaced with conventional NDC-to-window mapping.
It returns a new four-element list without changing the Vector or its storage.
Repeated calls are independent. Zero Vectors retain `ZeroDivisionError`.
Width/height are used arithmetically without a new positivity rule. Unsupported
shapes/types and nonfinite inputs receive no new validation or error policy.

```python
from gem.common import getViewPort
from gem.vector import Vector

assert getViewPort(Vector(2, [3, 4]), 100, 200) == [83, 184, 100, 200]
# Including Z changes the norm, even though only XY are returned.
assert getViewPort(Vector(3, [3, 4, 12]), 130, 260) == [83, 174, 130, 260]
```

This helper neither sets an OpenGL viewport nor projects world/object points.
For modelview/projection matrices, perspective division and window depth,
use the existing `project`/`unproject` API.

## Compatibility

Cross-dimension comparisons previously returned inconsistent values or raised
IndexError; empty comparisons returned None. They now return the documented
booleans. Returning clamp no longer changes caller data; code relying on that
side effect must explicitly assign the result or use `i_clamp` on its receiver.
External references to old receiver storage remain unchanged after `i_clamp`.
Viewport calls with supported Vector inputs now execute instead of raising the
Vector-subscripting TypeError. The normalized-coordinate formula, zero error,
function names/signatures and namespace remain unchanged.
