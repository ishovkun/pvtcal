#!/usr/bin/env python
import utils

class PVDGtable:
    p = []             # pressure vector
    B = []             # volume factor
    mu = []            # viscosity
    p_atm = 101325.0/100000. # atmospheric


    def __init__(self, data):
        """
        data is a list of floats
        """
        self.readTable_(data)
        if (self.p[0] > self.p_atm): self.extrapolateToStandardConditions()
        self.checkConsistency_()

    def checkConsistency_(self):
        for i in range(1, len(self.p)):
            assert self.B[i] <= self.B[i-1], "Bg must decrease with pressure"
            assert self.mu[i] >= self.mu[i-1], "mug must increase with pressure"

    def getVolumeFactor(self, p):
        return self.getProperty_(self.B, p)

    def getViscosity(self, p):
        return self.getProperty_(self.mu, p)

    def getProperty_(self, prop, p):
        pos, i, j = utils.findSurroundingElements(p, self.p)
        assert pos == utils.Position.BETWEEN, "extrapolation not allowed"
        y = utils.interp_lin( self.p[i], self.p[j], prop[i], prop[j], p)
        return y

    def readTable_(self, data):
        assert len(data) % 3 == 0, "Wrong entry in PVDG"
        for i in range(0, len(data), 3):
            self.p.append( data[i] )
            self.B.append( data[i+1] )
            self.mu.append( data[i+2] )

    def extrapolateToStandardConditions(self):
        # this interpolation gives a relatively good
        # extrapolation on SPE1
        mu_atm = utils.interp_lin( self.p[0], self.p[1],
                                   self.mu[0], self.mu[1], self.p_atm)
        B_atm = 1.0

        assert mu_atm > 0, "linear extrapolation failed"
        self.p.insert(0, self.p_atm)
        self.B.insert(0, B_atm)
        self.mu.insert(0, mu_atm)

    def extrapolateProperties(self, p):
        """
        currently does linear extrapolation
        """
        assert p > self.p[-1] or p < self.p[0], "this is not extrapolation"
        B = (self.B[-1] - self.B[-2]) / (self.p[-1] - self.p[-2]) * (p - self.p[-2])
        mu = (self.mu[-1] - self.mu[-2]) / (self.p[-1] - self.p[-2]) * (p - self.p[-2])
        if (p > self.p[-1]):
            self.p.append(p)
            self.B.append(B)
            self.mu.append(mu)
        else: # if p < self.p[0]
            self.p.insert(0, p)
            self.B.insert(0, B)
            self.mu.insert(0, mu)

