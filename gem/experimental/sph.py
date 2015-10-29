import math

from gem.experimental.legendre import Legendre


# n! where n >= 0
def Factorial(n):
    if n <= 1:
        return 1

    result = n

    while n > 1:
        n -= 1
        result *= n

    return result

# Normalization constant for a Spherical Harmonic function
def K(l, m):
    K = ((2.0 * l + 1.0) * Factorial(l - m)) / ((4.0 * math.pi) * Factorial(l + m))
    return math.sqrt(K)


# Sample a Spherical Harmonic function Y(l, m) at a point on the unit sphere
def SPH(l, m, theta, phi):
    root2 = math.sqrt(2.0)
    if m == 0:
        return K(l, 0) * Legendre(l, m, math.cos(theta)).run()
    elif m > 0:
        return root2 * K(l, m) * math.cos(m * phi) * Legendre(l, m, math.cos(theta)).run()
    elif m < 0:
        return root2 * K(l, -m) * math.sin(-m * phi) * Legendre(l, -m, math.cos(theta)).run()
    else:
        print ("WTF... The m is ...")
        return 0
