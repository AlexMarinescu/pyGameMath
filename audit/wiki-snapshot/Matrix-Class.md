## Matrix Class features breakdown:

***

**Operators:**
* multiplication (vector or matrix)
* division (scalar)

**Special cases:**

These functions are not members of the Matrix class, they return Matrix and Vector type.

* orthographic(left, right, bottom, top, zNear, zFar) (Orthographic Projection 4x4 matrix)
* perspective(fov, aspect, znear, zfar) (Perspective projection matrix 4x4. FOVY)
* perspectiveX(fov, aspect, znear, zfar) (Perspective projection matrix 4x4. FOVX)
* lookAt(eye, center, up) (returns a view matrix 4x4)
* project(obj, model, proj, viewport) (returns a 3D vector, projects to viewport)
* unproject(winx, winy, winz, modelview, projection, viewport) (Unproject a point from the screen and return the object coordinates.)

**Functions:**
* scale(value) (scale matrix by a vector)
* det() (returns the determinant of the matrix)
* inverse() (returns the inverse of the matrix)
* rotate(axis, theta) (rotate the matrix by axis and angle)
* translate(vecA) (translate the matrix using a vector)
* transpose() (Transpose the matrix)
* shearXY(x, y) (shear the matrix along the x, y)
* shearYZ(y, z) (shear the matrix along the z, y)
* shearXZ(x, z) (shear the matrix aling the x, z)

**Identities:**

These are not part of the Matrix class, so their output needs to be passed as input to the Matrix constructor.

* identity(size) (returns a matrix identity)
* zero_matrix(size) (returns a matrix filled with zeros)


## Examples:

It is a NxN Matrix class. However, some functions only apply for certain dimensions only so for example a 8x8 matrix will be limited in functionality.

**Declaration of a matrix:**

    from gem import matrix

    matrixA = matrix.Matrix(4)

    mdata = [[1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0]]

    matrixB = matrix.Matrix(4, data=mdata)

    matrixC = matrix.Matrix(4, data=matrix.identity(4))


    # Output:
    # MatrixA:  0.0, 0.0, 0.0, 0.0
    #           0.0, 0.0, 0.0, 0.0
    #           0.0, 0.0, 0.0, 0.0
    #           0.0, 0.0, 0.0, 0.0

    # MatrixB:  1.0, 0.0, 0.0, 0.0
    #           0.0, 1.0, 0.0, 0.0
    #           0.0, 0.0, 1.0, 0.0
    #           0.0, 0.0, 0.0, 1.0

    # MatrixC:  1.0, 0.0, 0.0, 0.0
    #           0.0, 1.0, 0.0, 0.0
    #           0.0, 0.0, 1.0, 0.0
    #           0.0, 0.0, 0.0, 1.0

**Operators usage:**
They have to be the same size matrices in order for the operators to work. The vector size needs to match that of the matrix, for example a 3D vector can only be multiplied by a 3x3 matrix, a 4D vector by a 4x4 matrix and so on.

    from gem import matrix
    from gem import vector

    # Initialize some matrices and vectors
    matrixA = matrix.Matrix(4, data=matrix.indetity(4))
    matrixB = matrix.Matrix(4, data=matrix.indetity(4))
    vectorA = vector.Vector(4, data=[2.0, 1.0, 3.0, 5.0])

    # Matrix Multiplication (NxN * NxN)
    matrixC = matrixA * matrixB

    # Vector Multiplication (NxN * ND)
    matrixD = matrixA * vectorA

    # Matrix Division (NxN / scalar)
    matrixE = matrixB / 2.0

**Function usage:**
Inside the matrix class most functions have two declartion, regular function name and with a "i_" in front of it. For example, the scale function is declared as "scale" and "i_scale", the former will return a new scaled matrix, the latter will scale the matrix in its place. So any function with "i_" in front of the function name will do the changes in place without returning a new matrix.

Example:
    
    from gem import matrix
    from gem import vector

    matrixA = matrix.Matrix(3, data=matrix.indetity(4))
    vectorS = vector.Vector(3, data=[2.0, 3.0, 4.0])

    # This will scale the matrix by creating a new one
    matrixS = matrixA.scale(vectorS)

    # This will scale the matrix in-place
    matrixA.i_scale(vectorS)

    # Return the determinant of the matrix
    matrixADet = matrixA.det()

    #Output:
    # matrixA:  1.0, 0.0, 0.0, 0.0
    #           0.0, 1.0, 0.0, 0.0
    #           0.0, 0.0, 1.0, 0.0
    #           0.0, 0.0, 0.0, 1.0

    # vectorS:  2.0, 3.0, 4.0

    # matrixS:  2.0, 0.0, 0.0, 0.0
    #           0.0, 3.0, 0.0, 0.0
    #           0.0, 0.0, 4.0, 0.0
    #           0.0, 0.0, 0.0, 1.0

    # marixADet: 24

