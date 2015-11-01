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
-  Cross Product
-  Vector3D, Scalar Multiplication
-  Logarithm
-  Exponential
-  Power
-  Liner Interpolation (LERP)
-  Spherical Interpolation (SLERP)
-  Spherical Interpoliaton No Invert
-  Quaternion Splines (SQUAD)

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