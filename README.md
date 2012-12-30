pyGameMath
==========
This is a math library written in python for 2D/3D game development which is also compatible with pypy. I made it while I was learning more about the math used in graphics development and for personal use in OpenGL related projects.
It's still a work in progress. However, it will be finished soon.

Supported features:
===================
* NxN Matrices
  * Transpose
  * Scale
  * NxN Matrix Multiplication
  * NxN Matrix * N Dimensions Vector Multiplication
  * NxN Matrix Inverse Calculation (Gauss Jordan Elimination)
  * 4x4 Perspective Projection Matrix
  * lookAt 4x4 Matrix
  * Translation (3x3, 4x4)
  * Rotation (2x2, 3x3, 4x4)
  * Shear (2x2, 3x3, 4x4)
  * Output
  * Stack
  
* N Dimensions Vectors
  * Dot Product
  * Cross Product (No 7D)
  * Refraction
  * Reflection
  * Invert
  * Normalize
  * Output
  * Stack
  
* Quaternions
  * Normalize
  * Dot Product
  * Rotation
  * Conjugate
  * Inverse
  * Negate
  * Rotate X, Y, Z
  * Arbitary Axis Rotation
  * To Rotation Matrix (4x4)
  * Cross Product
  * Vector3D, Scalar Multiplication
  * Logarithm
  * Exponential
  * Power
  * Liner Interpolation (LERP)
  * Spherical Interpolation (SLERP)
  * Quaternion Splines (SQUAD)
  * Output
  
* Plane
  * Define using
     * 3 Vectors
     * Point and Normal
     * Manual input
  * Dot Product
  * Normalize
  * Best fit normal and D value
  * Distance from plane to a point
  * Point location
  * Output
  
* Point
  >Basically this is a wrap around the vector class. It adds a location parameter so you can classify the point. 
  >
  >However, the location param could have been in the vector class itself but meh.
  >
  >This will probably change in the future.
  
* Ray
  * Rotate using Matrix
  * Rotate using Quaternions
  * Translate
  * Output

* Legendre Polynomial
  * For spherical harmonics
  * (l - m)PML(x) = x(2l - 1)PML-1(x) - (l + m -1)PML-2(x)
  * PMM(x) = (-1)^m * (2m - 1)!!(1 - x^2)^m/2
  * PMM+1(x) = x(2m + 1)PMM(x)
  
* Spherical Harmonics (Not complete yet)
  * Normalization Constant (K)
  * Sample a Spherical Harmonic function Y(l, m)
  * Sample
  * Generate samples
  
TO DO:
===============
* AABB/OBB/Sphere/Cylinder (The intersections between a ray and each, will be inside the ray class)
* Matrix Normalizing and Finding the Determinant.
* Tutorial on how to use everything.

Note:
=====
* Vectors Cross Product is only for 3D and 7D vectors, but currently only the 3D version is implemented.
* [Multivectors](http://en.wikipedia.org/wiki/Multivector) and [Bivector](http://en.wikipedia.org/wiki/Bivector) operations are not yet implemented.
* The only matrix transformation that will be implemented for now are going to be the ones up to 4x4 matrices.
* While the Vector and Matrix Stacks are not necessary and can be created on the go by anyone, they simplify and organize things a bit better.
* Optimizations are going to be done more indepth once the whole thing is finished and perhaps when I notice that a certain operations/transformation takes too long to execute.
* Also do report any bugs or speed related issues because some functions did not get a chance to be tested.