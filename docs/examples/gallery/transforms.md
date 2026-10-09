# Transformation order

![Transformation panels](../../../examples/showcase/output/transforms.svg)

The dashed rectangle is the original; the solid rectangle is the result.
Positions are explicit Vector4 values with w=1. Scale doubles X, rotation is
+90° about Z, and translation is (3,-1,0).

For gem's row-vector convention, S*R*T applies scale, then rotation, then
translation. Independently, (x,y) → (2x,y) → (-y,2x) → (3-y,2x-1).
S*T*R instead yields (1-y,2x+3), visibly moving the rectangle to another region.
`Matrix * Vector` is wrapper syntax for the mathematical row product vM.

The six panels use current Matrix4 APIs; no implicit promotion, perspective
division or changes to multiplication order are involved.

[Source](../../../examples/showcase/scenes.py) ·
[Transform tutorial](../../tutorials/transformations.md) ·
[Matrix API](../../api/matrix.md) · [Gallery](../index.md)
