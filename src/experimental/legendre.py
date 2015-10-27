import math
import six.moves as sm

# Legendre Polynomial Class
class Legendre (object):

    def __init__(self, l, m, x):
        self.l = l
        self.m = m
        self.x = x
        self.P = 1.0
        self.PM1 = 0.0
        self.PML = 0.0

    def mGreaterThan0(self):
        if self.m > 0:
            sqrtOneMinusXSquared = math.sqrt((1.0 - self.x) * (1.0 + self.x))
            f = 1.0
            for _ in sm.range(1, self.m + 1):
                self.P *= (-f) * sqrtOneMinusXSquared
                f += 2.0
        else:
            self.P = 1.0

    def calculatePM1(self):
        self.mGreaterThan0()
        self.PM1 = (self.x) * (2.0 * self.m + 1.0) * self.P

    def calculatePML(self, i):
        self.calculatePM1()
        self.PML = (((self.x) * (2.0 * i - 1.0) * self.PM1 - (i + self.m - 1.0) * self.P)) / (i - self.m)

    def run(self):
        if self.l == self.m:
            self.mGreaterThan0()
            return self.P

        elif self.l == (self.m + 1):
            self.calculatePM1()
            return self.PM1

        else:
            for i in sm.range(self.m + 2, self.l + 1):
                self.calculatePML(i)
                self.P = self.PM1
                self.PM1 = self.PML
            return self.PML
