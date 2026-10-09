# Vector directions and displacement

![Vector geometry](../../../examples/showcase/output/vectors.svg)

The diagram uses a=(3,1,0), b=(-1,2,0). Moving b's tail to a's head reaches
(2,3,0), the vector sum. Normalization produces (3/√10,1/√10,0) without changing
a. The dot product is 3×(-1)+1×2=-1, so the angle is obtuse. The cross product's
Z component is 3×2-1×(-1)=7: +Z points out of the XY diagram.

Arrow coordinates and labels come from `Vector`, `.normalize()`, `.dot()` and
`cross()`. The drawing utility only maps XY coordinates into screen positions.

[Source](../../../examples/showcase/scenes.py) ·
[Measurements](../../../examples/showcase/output/measurements.json) ·
[Vector tutorial](../../tutorials/vectors.md) · [Vector API](../../api/vector.md) ·
[Gallery](../index.md)
