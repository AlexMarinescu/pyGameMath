# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import math

# Sinc function
def sinc(x):
    if math.fabs(x) < 1.0e-4:
        return 1.0
    else:
        return math.sin(x) / x

# Scalar lerp
def scalarLerp(a, b, time):
    return a + time * (b - a)
