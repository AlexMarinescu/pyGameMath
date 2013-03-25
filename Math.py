import math
import sys
import random

# Constants
PI = 3.14159265358979323846
E  = 2.71828182845904523536

# Sinc function
def sinc(x):
    if math.fabs(x) < 1.0e-4:
        return 1.0
    else:
        return math.sin(x) / x

# Scalar lerp
def scalarLerp(a,b,time):
    return a + time*(b-a)
