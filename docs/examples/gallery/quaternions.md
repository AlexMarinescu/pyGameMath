# Interpolating orientation

![SLERP frames](../../../examples/showcase/output/quaternions.svg)

Identity and +120° about Z define the endpoints. Shortest-path SLERP at
0, 1/4, 1/2, 3/4 and 1 yields rotations of 0°,30°,60°,90° and120°.
Blue is the rotated local X axis; orange is local Y. Active positive rotation
moves +X toward +Y, preserving perpendicularity and unit length.

Independent references use q=(cos(θ/2),0,0,sin(θ/2)), X'=(cosθ,sinθ,0) and
Y'=(-sinθ,cosθ,0). Quaternion ordering remains [w,x,y,z]. The example uses
`quat_from_axis_angle`, `quat_slerp` and `quat_rotate_vector` without changing
caller-owned inputs. This is SLERP; neither legacy three-control SQUAD nor
four-control squad4 is needed for this two-endpoint demonstration.

[Source](../../../examples/showcase/scenes.py) ·
[Quaternion tutorial](../../tutorials/quaternions.md) ·
[Quaternion API](../../api/quaternion.md) · [Gallery](../index.md)
