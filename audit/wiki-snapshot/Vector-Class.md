## Vector Class features breakdown:

***

**Operators:**
* adding (Addition of two vectors)
* subtracting (Subtraction of two vectors)
* multiplication (scalar only)
* division (scalar only)
* negation (Negates every value inside the vector)
* equal boolean test (==)

**Special cases:**
These functions are not members of the Vector class.
* lerp(vecA, vecB, time) (linear interpolation between two vectors, returns vectors same size as the input)
* cross(vecA, vecB) (Cross product for 3D vectors only)
* reflect(incidentVec, norm) (reflect a 3D vector)
* refract(IOR, incidentVec, normal) (refract a 3D vector)

**Functions:**
* clone() (Returns exact copy of a vector)
* one() (Returns a vector with all its components set to 1)
* zero() (Returns a vector with all its components set to 0)
* negate() (Negates a vector)
* maxV(vectorB) (Returns the max vector between two vectors)
* maxS() (Returns the max component between two vectors)
* minV(vectorB) (Returns the min vector between two vectors)
* minS() (Returns the min component between two vectors)
* magnitude() (Returns magnitude of a vector)
* clamp(size, value, minS, maxS) (Clamps a vector)
* normalize() (Normalized the vector)
* dot() (dot product of two vectors)
* isInSameDirection(otherVec) (Returns a boolean)
* isInOppositeDirection(otherVec) (Returns a boolean)
* barycentric(a, b, c) (Get a coordinates with respect to a triangle)
* transform(position, matrix) (matrix multiplication for transformation)

**Swizzling:**

This is similar to GLSL and will work with 3D or 4D vectors. Returns a 2D or 3D vector depends on the case.
* xy()
* yz()
* xz()
* xw()
* yw()
* zw()
* xyw()
* yzw()
* xzw()
* xyw()

**Identities:**

These are only for 3D vectors.
* right() (Returns [1.0, 0.0, 0.0])
* left() (Returns [-1.0, 0.0, 0.0])
* front() (Returns [0.0, 0.0, -1.0])
* back() (Returns [0.0, 0.0, 1.0])
* up() (Returns [0.0, 1.0, 0.0])
* down() (Returns [0.0, -1.0, 0.0])


## Examples:

It is a 1D, 2D, ... , ND Vector class. However, some functions only apply for certain dimensions only so for example a 8D vector will be limited in functionality.

**Declaration of a vector:**

    from gem import vector

    vectorA = vector.Vector(3, data=[1.0, 2.0, 3.0])
    vectorB = vector.Vector(3)
    vectorC = vector.Vector(3).one()
    vectorD = vector.Vector(3).zero()
    
    vectorData = [5.0, 6.0, 9.0]
    vectorE = vector.Vector(3, vectorData)


    # Output:
    # VectorA: 1.0, 2.0, 3.0
    # VectorB: 0.0, 0.0, 0.0
    # VectorC: 1.0, 1.0, 1.0
    # VectorD: 0.0, 0.0, 0.0
    # VectorE: 5.0, 6.0, 9.0

**Operators usage:**
They have to be the same size vectors in order for the operators to work.

	from gem import vector

	vectorA = vector.Vector(3, data=[5.0, 3.0, 8.0])
	vectorB = vector.Vector(3, data=[1.0, 4.0, 5.0])

	vectorAaddB = vectorA + vectorB

	vectorAsubB = vectorA - vectorB

	vectorAMulScalar = vectorA * 2.0

	vectorBDivScalar = vectorB / 3.0

	vectorANeg = -vectorA

	vectorAdotB = vectorA.dot(vectorB)

	# Since we are using 3D vectors here we can use the cross product
        # Note this function is not a member of the Vector class
	vectorAcrossB = vector.cross(vectorA, vectorB)

	# Output:
	# VectorAaddB = 6.0, 7.0, 13.0
	# VectorAsubB = 4.0, -1.0, 3.0
	# VectorAMulScalar = 10.0, 6.0, 16.0
	# VectorBDivScalar = 0.33333, 1.33333, 1.66666
	# VectorANeg = -5.0, -3.0, -8.0

	# vectorAdotB = 57 (This returns a scalar)
	# vectorAcrossB = -17, -17, 17

**Function usage:**
Inside the vector class most functions have two declartion, regular function name and with a "i_" in front of it. For example, the normalize function is declared as "normalize" and "i_normalize", the former will return a new normalized vector, the latter will normalize the vector in its place. So any function with "i_" in front of the function name will do the changes in place without returning a new vector.

Example:
	
	from gem import vector

	vectorA = vector.Vector(3, data[2.0, 4.0, 8.0])

	# This will clone vectorA to vectorB
	vectorB = vectorA.clone()

	# Normalize vectorA
	vectorN = vectorA.normalize()

	# Normalize vectorB in-place, same result as vectorN
	vectorB.i_normalize()

	#Output:
	# vectorN = 0.218218, 0.436436, 0.872872
	# vectorB = 0.218218, 0.436436, 0.872872
