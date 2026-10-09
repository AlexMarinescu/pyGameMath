"""Unnormalized Legendre functions with the Condon–Shortley phase.

The supported associated domain is integer 0 <= m <= l, -1 <= x <= 1.
No clipping, normalization, extreme-order stabilization or new invalid-input
policy is imposed. Ordinary polynomials (m=0) also evaluate outside [-1,1].
"""
import math


class Legendre(object):
    """Evaluate P_l^m(x); run() does not change polynomial scratch fields.

    The historical helpers explicitly update P, PM1 and PML. Repeated helper
    initialization is deterministic rather than multiplying previous state.
    """
    def __init__(self, l, m, x):
        self.l = l
        self.m = m
        self.x = x
        self.P = 1.0
        self.PM1 = 0.0
        self.PML = 0.0

    def _seed(self):
        value = 1.0
        if self.m > 0:
            root = math.sqrt((1.0-self.x)*(1.0+self.x))
            for order in range(1, self.m+1):
                value *= -(2.0*order-1.0)*root
        return value

    def _evaluate(self, degree, previous, current):
        if degree == self.m:
            return previous
        if (self.m == 0 and (type(self.x) is float or type(self.x) is int)
                and (self.x > 1 or self.x < -1)):
            # For integer l <= 64 and |x| <= 2, each recurrence step grows
            # its absolute state by less than 8. Even the undivided products
            # remain below 2**(3*l+8), far inside binary64's finite range.
            if not (type(degree) is int and degree <= 64 and -2 <= self.x <= 2):
                return self._evaluate_extrapolation(degree, previous, current)
        for index in range(self.m+2, degree+1):
            following = (self.x*(2.0*index-1.0)*current
                         - (index+self.m-1.0)*previous)/(index-self.m)
            previous, current = current, following
        return current

    def _evaluate_extrapolation(self, degree, previous, current):
        for index in range(self.m+2, degree+1):
            following = (self.x*(2.0*index-1.0)*current
                         - (index+self.m-1.0)*previous)/(index-self.m)
            if (math.isinf(following) and math.isfinite(self.x)
                    and math.isfinite(previous) and math.isfinite(current)):
                # Align product exponents before subtraction and division.
                # Scaling by powers of two avoids an overflowing numerator
                # without changing the ordinary finite recurrence's rounding.
                x, xe = math.frexp(self.x)
                c, ce = math.frexp(current)
                p, pe = math.frexp(previous)
                exponent = max(xe+ce, pe)
                numerator = (math.ldexp(x*c*(2.0*index-1.0), xe+ce-exponent)
                             - math.ldexp(p*(index-1.0), pe-exponent))
                try:
                    following = math.ldexp(numerator/index, exponent)
                except OverflowError:
                    # A genuinely unrepresentable result retains infinity.
                    following = math.copysign(float('inf'), numerator)
            previous, current = current, following
        return current

    def mGreaterThan0(self):
        self.P = self._seed()

    def calculatePM1(self):
        self.mGreaterThan0()
        self.PM1 = self.x*(2.0*self.m+1.0)*self.P

    def calculatePML(self, i):
        self.calculatePM1()
        self.PML = self._evaluate(i, self.P, self.PM1)

    def run(self):
        # Retain the historical empty-recurrence result outside m <= l.
        if self.l < self.m:
            return self.PML
        previous = self._seed()
        current = self.x*(2.0*self.m+1.0)*previous
        return self._evaluate(self.l, previous, current)
