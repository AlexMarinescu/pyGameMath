# Camera, projection and unprojection

Import all functions below from `gem.matrix`.
[Source](../../gem/matrix.py), [known-answer and error tests](../../tests/test_projection.py).

Cameras use a right-handed frame looking along negative Z; near/far distances
are positive in ordinary use. OpenGL-style NDC depth is [-1,1]. Functions return
fresh Matrix4/Vector3 objects and preserve inputs. Invalid frusta/viewports have
legacy arithmetic/index errors, not a uniform validator; exact zero denominators
raise ZeroDivisionError. No clipping or automatic depth clamping occurs.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `orthographic(left, right, bottom, top, zNear, zFar)` | Numeric `left,right,bottom,top,zNear,zFar`; fresh Matrix4 with diagonal [2/(right−left),2/(top−bottom),−2/(zFar−zNear),1] and final-row offsets. Camera depths −zNear/−zFar map to NDC −1/+1. |
| `perspective(fov, aspect, znear, zfar)` | `fov`: vertical degrees, `aspect`: width/height, numeric `znear,zfar`; fresh Matrix4. XY diagonals [1/(aspect*tan(fov/2)),1/tan(fov/2)], M22=−(far+near)/(far−near), M23=−1, M32=−2*far*near/(far−near); clip W=−camera Z. |
| `perspectiveX(fov, aspect, znear, zfar)` | `fov`: horizontal degrees, aspect width/height, `znear,zfar`: numeric; fresh Matrix4, XY diagonals [cot(fov/2),aspect*cot(fov/2)] and same depth/W convention as perspective. No new frustum validation. |
| `lookAt(eye, center, up)` | `eye,center,up`: Vector3; fresh view Matrix4 using f=normalize(center−eye), s=normalize(f×normalize(up)), u=s×f; final row [−s·eye,−u·eye,f·eye,1]. Zero direction/up or parallel frame raises ZeroDivisionError. Inputs preserved. |
| `project(obj, model, proj, viewport)` | Explicit `obj`: Vector4 [x,y,z,w], no Vector3 promotion. `model,proj`: Matrix4 or raw 4×4 lists, mixed allowed. `viewport`: [x,y,width,height]. Return fresh window Vector3 after model*proj and clip-W division; zero clip W raises ZeroDivisionError. No inverse required. |
| `unproject(winx, winy, winz, modelview, projection, viewport)` | Numeric `winx,winy,winz`; `modelview,projection`: Matrix4/raw 4×4, mixed allowed; viewport [x,y,width,height]. Fresh object Vector3 from reversed viewport and inverse(modelview*projection), with homogeneous division. Singular combined matrix raises ZeroDivisionError; exact zero output W returns historical zero Vector3 sentinel. |

## Window mapping

Viewport origin is lower-left, +Y upward. From clip coordinates c, NDC=c.xyz/c.w:
window XY=viewport.xy+(NDC.xy+1)*viewport.wh/2;
window Z=(NDC.z+1)/2. Visible near/far map to 0/1; outside values remain unclamped.
Unprojection applies NDC XY=2*(window−origin)/extent−1 and NDC Z=2*windowZ−1,
then the combined inverse. These functions are distinct from
[common.getViewPort](common.md#historical-viewport) and OpenGL state setting.

```python
from gem.matrix import Matrix, perspective, orthographic, project, unproject
from gem.vector import Vector

viewport = [10, 20, 200, 100]
p = perspective(90, 2, 1, 9)
near = project(Vector(4, [0, 0, -1, 1]), Matrix(4), p.matrix, viewport)
far = project(Vector(4, [0, 0, -9, 1]), Matrix(4).matrix, p, viewport)
assert all(abs(a - b) < 1e-14 for a, b in zip(near.vector, [110, 70, 0]))
assert all(abs(a - b) < 1e-14 for a, b in zip(far.vector, [110, 70, 1]))
world = unproject(110, 70, 0, Matrix(4), p, viewport)
assert all(abs(a - b) < 1e-14 for a, b in zip(world.vector, [0, 0, -1]))
o = orthographic(-2, 2, -1, 1, 1, 9)
assert project(Vector(4, [2, 1, -5, 1]), Matrix(4), o, viewport).vector == [210.0, 120.0, 0.5]
```

Near-zero W handling, malformed dimensions and invalid frusta remain
[open policies](decisions.md), not new checks. See [transforms](transformations.md) and [index](index.md).
