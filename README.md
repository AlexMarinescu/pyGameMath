[![Build Status](https://travis-ci.org/AlexMarinescu/pyGameMath.svg?branch=master)](https://travis-ci.org/AlexMarinescu/pyGameMath)

![ScreenShot](https://raw.github.com/AlexMarinescu/pyGameMath/master/data/pyGameMathLogo.png)

pyGameMath
==========
This is a math library written in python for 2D/3D game development which is also compatible with pypy. I made it while I was learning more about the math used in graphics development and for personal use in OpenGL related projects.
It's still a work in progress. However, it will be finished soon.

####Update 1:####
* The old version which is fully OOP is inside the folder "oop math lib", some bugs are still present.
* The new version has most of the features functional based for performance, it is inside the folder "src" and all the bugs from the previous have been nuked.
* There is a new folder called "experimental" inside the "src" that stores the features that are unfinished.
* Spherical harmonics were moved under "experimental".

####Update 2:####
* Updated the matrix and vector, fixed a lot of bugs, still dynamic class that allows Nth dimensions
* A lot of the specific 2D, 3D, 4D functions are outside the dynamic class
* Now it allows in-place transformations
* Some math operations do not support Nth dimensions and will return exceptions
* Removed LU decomposition
* Temporary Removed NxN Matrix Normalization
* Temporary Removed NxN Matrix Determinant
* Temporary Removed NxN Matrix Inverse Calculation
* Quaternion and the rest need to be updated to match the new style.

##Supported features:##

###NxN Matrices###
* Transpose
* Scale
* NxN Matrix Multiplication
* NxN Matrix * N Dimensions Vector Multiplication
* 4x4 Perspective Projection Matrix
* lookAt 4x4 Matrix
* Translation (3x3, 4x4)
* Rotation (2x2, 3x3, 4x4)
* Shear (2x2, 3x3, 4x4)
* Project
* Unproject
  
###N Dimensions Vectors###
* Dot Product
* Cross Product (No 7D)
* Refraction
* Reflection
* Invert
* Normalize
  
###Quaternions###
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
  
###Plane###
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
  
###Ray###
* Rotate using Matrix
* Rotate using Quaternions
* Translate
* Output

###Legendre Polynomial (Experimental)###
* For spherical harmonics
* (l - m)PML(x) = x(2l - 1)PML-1(x) - (l + m -1)PML-2(x)
* PMM(x) = (-1)^m * (2m - 1)!!(1 - x^2)^m/2
* PMM+1(x) = x(2m + 1)PMM(x)

###Spherical Harmonics (Experimental)###
* Normalization Constant (K)
* Sample a Spherical Harmonic function Y(l, m)
* Sample
* Generate samples
* Coefficients from Irradiance map
  
###TO DO:###
* AABB/Sphere/Cylinder  (Everything is done, just needs to be cleaned up before pushed.)
* GJK (Done, needs to be tested)
* Primitives/Supports for GJK (Basic ones are done, needs to be updated)
* Spherical Harmonics.
* Tutorial on how to use everything.

###Note:###
* Vectors Cross Product is only for 3D and 7D vectors, but currently only the 3D version is implemented.
* While the Vector and Matrix Stacks are not necessary and can be created on the go by anyone, they simplify and organize things a bit better.
* Optimizations are going to be done more indepth once the whole thing is finished and perhaps when I notice that a certain operations/transformation takes too long to execute.
* Also do report any bugs or speed related issues because some functions did not get a chance to be tested.
