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
        for index in range(self.m+2, degree+1):
            following = (self.x*(2.0*index-1.0)*current
                         - (index+self.m-1.0)*previous)/(index-self.m)
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
