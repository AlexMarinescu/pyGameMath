|ScreenShot|

pyGameMath |Build Status| |Code Health| |Codacy Badge|
======================================================

| This is a math library written in python for 2D/3D game development
  which is also compatible with pypy. I made it while I was learning
  more about the math used in graphics development and for personal use
  in OpenGL related projects.
| It’s still a work in progress.

Dependencies:
-------------
It uses six to allow support between python2.x and python3.x.

Install:
--------
To install the library just do

.. code:: Python

    pip install gem

It will install the dependicies automatically.

Documentation and Examples:
---------------------------
The examples on how to use the library and more info are maintained on the github wiki:

`Wiki Link <https://github.com/explosiveduck/pyGameMath/wiki>`_

Supported features:
~~~~~~~~~~~~~~~~~~~

NxN Matrices:
'''''''''''''

-  Transpose
-  Scale
-  NxN Matrix Multiplication
-  NxN Matrix \* N Dimensions Vector Multiplication
-  4x4 Perspective Projection Matrix
-  lookAt 4x4 Matrix
-  Translation (3x3, 4x4)
-  Rotation (2x2, 3x3, 4x4)
-  Shear (2x2, 3x3, 4x4)
-  Project
-  Unproject
-  Orthographic Projection
-  Perspective Projection
-  lookAt 4x4 matrix
-  Determinant 2x2, 3x3, 4x4
-  Inverse 2x2, 3x3, 4x4

N Dimensions Vectors:
'''''''''''''''''''''

-  Dot Product
-  Cross Product (3D, No 7D as of now)
-  2D get angle of vector
-  2D -90 degree rotation
-  2D +90 degree rotation
-  Refraction
-  Reflection
-  Negate
-  Normalize
-  Linear Interpolation
-  Max Vector/Scalar
-  Min Vector/Scalar
-  Clamp
-  Transform 
-  Barycentric 
-  isInSameDirection test
-  isInOppositeDirection test
-  3D Vector swizzling, similar to GLSL
-  3D Vector idenitities

Quaternions:
''''''''''''

-  Normalize
-  Dot Product
-  Rotation
-  Conjugate
-  Inverse
-  Negate
-  Rotate X, Y, Z
-  Arbitary Axis Rotation
-  From angle Rotation
-  To Rotation Matrix (4x4)
-  From Rotation Matrix (4x4)
-  Vector3D, Scalar Multiplication
-  Logarithm
-  Power
-  Liner Interpolation (LERP)
-  Spherical Interpolation (SLERP)
-  Spherical Interpoliaton No Invert
-  Quaternion Splines (SQUAD)

See the `quaternion API guide <docs/QUATERNIONS.md>`_ for helper return
types, angle units, input ownership and rotation examples.

Plane:
''''''

-  Define using

   -  3 Vectors
   -  Point and Normal
   -  Manual input

-  Dot Product
-  Normalize
-  Best fit normal and D value
-  Distance from plane to a point
-  Point location
-  Output
-  Flip

Ray:
''''

-  Rotate using Matrix
-  Rotate using Quaternions
-  Translate
-  Output

Legendre Polynomial (Experimental, not complete):
'''''''''''''''''''''''''''''''''''''''''''''''''

-  For spherical harmonics
-  (l - m)PML(x) = x(2l - 1)PML-1(x
-  Irradiance maps

.. |ScreenShot| image:: https://raw.github.com/AlexMarinescu/pyGameMath/master/data/pyGameMathLogo.png
.. |Build Status| image:: https://travis-ci.org/explosiveduck/pyGameMath.svg?branch=master
   :target: https://travis-ci.org/explosiveduck/pyGameMath
.. |Code Health| image:: https://landscape.io/github/explosiveduck/pyGameMath/master/landscape.svg?style=flat
   :target: https://landscape.io/github/explosiveduck/pyGameMath/master
.. |Codacy Badge| image:: https://api.codacy.com/project/badge/907e4230379f40a8bedcfc0a9a0ed43c
   :target: https://www.codacy.com
Bezier evaluation
-----------------

Use ``gem.bezier`` for supported scalar and Vector Bezier evaluation::

    from gem.bezier import cubicBezierPoint, BezierPath
    from gem.vector import Vector

    midpoint = cubicBezierPoint(0.5, Vector(2, [0, 0]),
                                Vector(2, [1, 2]), Vector(2, [2, 2]),
                                Vector(2, [3, 0]))  # [1.5, 1.5]
    path = BezierPath()
    path.setControlPoints([0, 1, 2, 3])
    midpoint = path.calculateBezerPoint(0, 0.5)  # 1.5

Quadratic evaluation is available as ``quadraticBezierPoint(t, p0, p1, p2)``.
Controls are preserved and parameters are not clamped. Cubic paths use
``3*k+1`` controls for ``k`` segments. Experimental Bezier imports remain compatible; sampling is supported in core. See
``audit/PHASE2F3A.md`` for migration and compatibility details.


Adaptive Bezier sampling
-----------------------

Sample a cubic path without changing its control points::

    from gem.bezier import BezierPath
    from gem.vector import Vector

    path = BezierPath()
    path.setControlPoints([Vector(2, [0, 0]), Vector(2, [1, 2]),
                           Vector(2, [2, 2]), Vector(2, [3, 0])])
    path.minimum_sqr_distance = 0.0001  # distance tolerance 0.01
    points = path.findDrawingPoints(0)  # ordered, includes both endpoints
    curves = path.getDrawingPoints()   # nested lists; shared joins appear once

The geometric flatness test applies to Vector2/Vector3 and scalar controls.
Depth is capped at 16; capped output may exceed tolerance. Source-point
``samplePoints(sourcePoints, minSqrDistance, maxSqrDistance, scale)`` uses
separate squared-distance thinning heuristics and rebuilds its generated path.
``interpolate(segmentPoints, scale)`` retains append-only behavior. Both
builders return None and preserve source storage. See ``audit/PHASE2F3B.md``
for limits, validation and compatibility details.

Legendre functions
------------------

Use the core module for unnormalized Legendre and associated Legendre values::

    from gem.legendre import Legendre

    value = Legendre(3, 0, 0.2).run()  # -0.28
    associated = Legendre(2, 2, 0.5).run()  # 2.25

Associated inputs use integer ``0 <= m <= l`` and ``-1 <= x <= 1``, with the
Condon–Shortley phase. ``run()`` preserves internal scratch state and is
repeatable. The experimental import remains a compatibility alias.
Spherical-harmonics normalization is separate. See ``audit/PHASE2F4.md`` for
numerical-domain and compatibility details.
