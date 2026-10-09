# Mathematical notation

These equations restate established contracts. TeX in the source becomes native
MathML at website build time; there is no external equation-rendering service.

## Vectors and row matrices

For $\mathbf{v}=(v_x,v_y,v_z)$, the Euclidean norm is

$$
\lVert\mathbf{v}\rVert=\sqrt{v_x^2+v_y^2+v_z^2}.
$$

This defines the mathematics, not the floating-point algorithm: gem uses stable
hypot-style lengths and scaled normalization. See [numerical accuracy](../tutorials/numerical.md).

A homogeneous translation with row-vector inputs is

$$
T=\begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&1&0\\t_x&t_y&t_z&1\end{pmatrix},\qquad
\mathbf{p}'=\mathbf{p}T.
$$

Thus $[x,y,z,1]T=[x+t_x,y+t_y,z+t_z,1]$. Directions use w=0.
The Matrix wrapper expression `T * p` evaluates this row product. A*B applies
A before B. [Transformation tutorial](../tutorials/transformations.md).

## Quaternions

A unit axis $\mathbf{u}$ and angle $\theta$ define

$$
q=\left[\cos(\theta/2),\,\mathbf{u}\sin(\theta/2)\right],\qquad
[0,\mathbf{v}']=q[0,\mathbf{v}]q^*.
$$

Storage is [w,x,y,z], and products use Hamilton multiplication. API units remain
function-specific; `quat_from_axis_angle` takes degrees.
[Quaternion reference](../api/quaternion.md).

## Spherical harmonics

Canonical coefficient index is $i=l(l+1)+m$. For canonical real orthonormal basis
functions, reconstruction and active rotation obey

$$
f(\mathbf{d})=\sum_{l=0}^{2}\sum_{m=-l}^{l}c_{lm}Y_{lm}(\mathbf{d}),\qquad
f_R(\mathbf{d})=f(R^{-1}\mathbf{d}).
$$

Condon–Shortley signs are retained. The first-order XYZ functions are
$Y_{1,-1}=-\sqrt{3/(4\pi)}y$, $Y_{1,0}=\sqrt{3/(4\pi)}z$,
and $Y_{1,1}=-\sqrt{3/(4\pi)}x$. Coefficients may be RGB.
Cosine convolution applies band factors $\pi$, $2\pi/3$, $\pi/4$ exactly once for
L0,L1,L2; reflected Lambertian radiance is $L_o=\rho E/\pi$.
[Lighting tutorial](../tutorials/lighting.md).

```python
import math
from gem.quaternion import quat_from_axis_angle, quat_rotate_vector
from gem.vector import Vector
q = quat_from_axis_angle([0, 0, 1], 90)
v = quat_rotate_vector(q, Vector(3, [1, 0, 0]))
assert math.isclose(v.vector[1], 1.0, abs_tol=1e-12)
```

The adjacent Python block follows the same active +Z convention. Equations,
code and [verified diagrams](../examples/index.md) describe the same mathematics.
