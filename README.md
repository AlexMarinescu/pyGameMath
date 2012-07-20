pyGameMath
==========
This is a math library written in python for 2D/3D game development which is also compatible with pypy. I made it while I was learning more about the math used in graphics development and for personal use in OpenGL related projects.
It's still a work in progress. However, it will be finished soon (possibly a week or so).

Supported features:
===================
* NxN Matrices
* 4x4 Perspective Projection Matrix
* lookAt 4x4 Matrix
* Translation and Rotation
* NxN Matrix Inverse Calculation (Gauss Jordan Elimination)
* N Dimensions Vectors
* Vectors Dot Product
* Vectors Cross Product
* Matrix Stack
* Vector Stack
* Quaternions

TO DO:
===============
* Some more 3D Specific Vector Operations.
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