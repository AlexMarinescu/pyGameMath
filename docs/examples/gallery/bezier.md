# Bezier curves and sample placement

![Bezier curves](../../../examples/showcase/output/bezier.svg)

Dashed polygons connect control points; blue polylines trace evaluated curves.
The quadratic midpoint is (2,2), the cubic midpoint (2,0.875). Orange markers on
the quadratic use uniform parameter steps of 1/8, which do not give equal arc
length steps or constant travel speed.

The cubic uses existing BezierPath subdivision at `minimum_sqr_distance` values
0.16 and 0.0025, corresponding to distances 0.4 and 0.05 in curve coordinates.
Large orange markers show the coarser path; small green markers show the finer
path. The current sampling criterion measures control-point distance to the
endpoint chord segment. The maximum depth is 16; depth-limited sampling can
exceed tolerance. This demonstration does not claim a universal error guarantee.

Independent checks evaluate the expanded cubic x=3t+3t²-2t³,
y=12t-30t²+19t³ at 1,001 parameter values and measure distance to each generated
polyline. Both paths include endpoints, increase in parameter order (X is
strictly increasing here), and the finer path reduces the measured error.

[Source](../../../examples/showcase/scenes.py) ·
[Measurements](../../../examples/showcase/output/measurements.json) ·
[Bezier tutorial](../../tutorials/bezier.md) ·
[Bezier API](../../api/bezier.md) · [Gallery](../index.md)
