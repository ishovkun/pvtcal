#!/usr/bin/env python
import numpy as np
import utils

Position = utils.Position


class PVTOtable:
    # saturated data
    p_bub = np.array(0)
    Rs_sat = np.array(0)
    B_sat = np.array(0)         # formation volume factor
    mu_sat = np.array(0)        # oil viscosity

    # unsaturated data
    # NOTE: p_usat = pressure - p_bub ; this table always starts with 0
    p_usat = np.array(0)        # 2d numpy array
    B_usat = np.array(0)        # volume factor
    mu_usat = np.array(0)       # volume factor
    Rs_usat = []

    # constants
    p_atm = 101325.0/100000.;           # atmospheric
    min_B_limit = 1e-6

    # simple parameters
    p_max = 0
    p_min = 0
    p_usat_rel_max = 0              # max relative pressure

    def __init__(self, data):
        """
        data is a list of lists of float values
        parsing is happening right here
        """
        self.readTable_(data)
        # self.convertUnits_()
        self.checkConsistency_()
        self.makeStandardConditionsBranch_()
        self.extrapolateUSatToMaximumRelativePressure_()
        self.checkConsistency_()


    def computeSaturatedPressure_(self, sat_lower, sat_upper, Rs):
        if (sat_lower < sat_upper):
            b1 = sat_lower
            b2 = sat_upper
            p_sat = (self.p_bub[b2] - self.p_bub[b1]) / \
                (self.Rs_sat[b2] - self.Rs_sat[b1]) * \
                (Rs - self.Rs_sat[b1]) + self.p_bub[b1]
        else: #  lower == upper never >
            p_sat = self.p_bub[sat_lower]
        return p_sat

    def getProperty(self, Rs, p, sat_values, usat_values):
        # find lower and higher rs value indices
        # pos, sat_lower, sat_upper = utils.findSurroundingElements(Rs, self.Rs_sat)
        # print(sat_lower, sat_upper)
        # print(self.Rs_sat)
        sat_lower, sat_upper = utils.binarySearchExtrapolationIndices(Rs, self.Rs_sat)
        # exit(0)
        # print(pos, self.Rs_sat[sat_lower], self.Rs_sat[ sat_upper ], Rs)

        # if (pos != pos.BETWEEN):
        #     raise LookupError("out of bounds and extrapolation not supported here")

        p_rel = p - self.computeSaturatedPressure_(sat_lower, sat_upper, Rs)

        # find undersaturated branches to interpolate between
        usat_lower, usat_upper = utils.binarySearchExtrapolationIndices(Rs, self.Rs_usat)
        # print(pos, self.Rs_usat[usat_lower], self.Rs_usat[ usat_upper ], Rs)
        # if (pos == Position.BETWEEN):
        # find usaturated relative pressures to get B values on
        # each of the under-saturated branches
        i11, i12 = utils.binarySearchExtrapolationIndices(p_rel, self.p_usat[usat_lower],
                                                 shift=self.p_usat[usat_lower][0])
        i21, i22 = utils.binarySearchExtrapolationIndices(p_rel, self.p_usat[usat_upper],
                                                 shift=self.p_usat[usat_upper][0])

        # interpolate on each usat branch to get two Bo values
        usat_value1 = (usat_values[usat_lower][i12] - usat_values[usat_lower][i11]) / \
                (self.p_usat[usat_lower][i12] - self.p_usat[usat_lower][i11]) * \
                (p_rel - (self.p_usat[usat_lower][i11] - (self.p_usat[usat_lower][0]))) +\
                usat_values[usat_lower][i11]

        usat_value2 = (usat_values[usat_upper][i22] - usat_values[usat_upper][i21]) / \
            (self.p_usat[usat_upper][i22] - self.p_usat[usat_upper][i21]) * \
            (p_rel - (self.p_usat[usat_upper][i21] - (self.p_usat[usat_upper][0]))) +\
            usat_values[usat_upper][i21]

        # NOTE: we should interpolate B wrt Rs between branches since
        # this causes bad interpolatiopn (trust me i checked)
        # B = (B2 - B1) / (self.Rs_usat[usat_upper] - self.Rs_usat[usat_lower]) * \
        #     (Rs - self.Rs_usat[usat_lower]) + B1

        # that's why we interpolate along the saturated B curve
        # first interpolate B_sat
        sat_value = (sat_values[sat_upper] - sat_values[sat_lower]) / \
                (self.Rs_sat[sat_upper] - self.Rs_sat[sat_lower]) * \
                (Rs - self.Rs_sat[sat_lower]) + sat_values[sat_lower]

        # use B_sat as x in interpolation
        usat_value = (usat_value2 - usat_value1) / \
            (usat_values[usat_upper][0] - usat_values[usat_lower][0]) * \
            (sat_value - usat_values[usat_lower][0]) + usat_value1

        return usat_value

    def getViscosity(self, Rs, p):
        return self.getProperty(Rs, p, self.mu_sat, self.mu_usat)

    def getVolumeFactor(self, Rs, p):
        return self.getProperty(Rs, p, self.B_sat, self.B_usat)

    def convertUnits_(self):
        # convert pressure from bars to Pa
        self.p_bub = [x*100000. for x in self.p_bub]
        for i in range(len(self.p_usat)):
            for j in range(len(self.p_usat[i])):
                self.p_usat[i][j] *= 100000.
        self.p_atm *= 100000.

    def makeStandardConditionsBranch_(self):
        """
        check if there is an unsaturated branch and
        make it if not
        """
        assert( max(self.p_bub) > self.p_atm )
        # at standard conditions Rs = 0
        if (self.p_bub[0] >= self.p_atm):

            # interpolate viscosity -- different from Alex' calcs
            mu_atm = interp_log( self.p_bub[0], self.p_bub[1],
                                 self.mu_sat[0], self.mu_sat[1], self.p_atm)

            # this is what Alex has
            # log_mu = (log(self.mu_sat[1]) - log(self.mu_sat[0])) / \
            #          ( log(self.p_bub[1]) - log(self.p_bub[0]) ) * \
            #          (log(self.p_atm) - log(self.p_bub[1])) + log(self.mu_sat[0])
            # mu_atm = np.exp(log_mu)

            self.p_bub.insert (0, self.p_atm)  # atmospheric
            self.Rs_sat.insert(0, 0.0)        # no gas dissolved on surface
            self.B_sat.insert (0, 1.0)        # B on surface = 1 (no gas dissolved)
            self.mu_sat.insert(0, mu_atm)     # B on surface = 1 (no gas dissolved)

            # create usat branch;
            self.Rs_usat.insert( 0, 0.0)
            self.p_usat.insert( 0, [ self.p_atm ])
            self.B_usat.insert( 0, [ 1.0 ])
            self.mu_usat.insert( 0, [ mu_atm ])
            # NOTE: the second point will be extrapolated later
        else:    # assign the first branch
            # create usat branch;
            self.Rs_usat.insert( 0, self.Rs_sat[0])
            self.p_usat.insert( 0, [ self.p_bub[0] ])
            self.B_usat.insert( 0, [ self.B_sat[0] ])
            self.mu_usat.insert( 0, [ self.mu_sat[0] ])

    def readTable_(self, data):
        Rs_sat = []; Rs_usat = []; p_bub = [];
        B_sat = []; mu_sat = []
        p_usat = []; B_usat = []; mu_usat = []
        for row in data:
            Rs_sat.append(row[0])
            p_bub.append(row[1])
            B_sat.append(row[2])
            mu_sat.append(row[3])

            if (len(row) > 4): # unsaturated
                rs = 0; p = []; B = []; mu = []
                rs = row[0]
                # p.append( row[1] - p_bub[-1])
                p.append( row[1] )
                B.append( row[2] )
                mu.append( row[3] )

                for i in range(4, len(row), 3): # next columns go with a skip
                    # p.append( row[i] - p_bub[-1] )
                    p.append( row[i] )
                    B.append( row[i+1] )
                    mu.append( row[i+2] )

                Rs_usat.append(rs)
                p_usat.append(p)
                B_usat.append(B)
                mu_usat.append(mu)

        self.Rs_sat = Rs_sat
        self.p_bub = p_bub
        
        self.B_sat = B_sat
        self.mu_sat = mu_sat

        self.Rs_usat = Rs_usat
        self.p_usat = p_usat
        self.B_usat = B_usat
        self.mu_usat = mu_usat

        # find maximum relative pressure
        self.p_usat_rel_max = max( self.p_bub )
        for i in range(len(self.Rs_usat)):
            # self.p_usat_rel_max = max( self.p_usat_rel_max, max(self.p_usat[i]) ) #
            self.p_usat_rel_max = max( self.p_usat_rel_max, max(self.p_usat[i]) - self.p_usat[i][0] )


    def extrapolateUSatToMaximumRelativePressure_(self):
        """
        extrapolate to the maximum unsaturated pressure:
        I want tables for all rs span the same p_usat range
        """
        create_sc_branch = False
        for i in range(len(self.p_usat)):
            # only a standard conditions branch previously
            # at this point can have a single value
            # remember that and treat it later
            if (len(self.p_usat[i]) < 2):
                create_sc_branch = True
            else:
                if (self.p_usat[i][-1] < self.p_usat[i][0] + self.p_usat_rel_max):
                    # extrapolate linearly the unsaturated properties to the
                    # maximum relative pressure
                    p = self.p_usat[i]
                    B = self.B_usat[i]
                    mu = self.mu_usat[i]
                    B_ext = utils.interp_lin( p[-2], p[-1], B[-2], B[-1], p[0] + self.p_usat_rel_max)
                    mu_ext = utils.interp_lin( p[-2], p[-1], mu[-2], mu[-1], p[0] + self.p_usat_rel_max)
                    p.append( p[0] + self.p_usat_rel_max )
                    B.append( B_ext )
                    mu.append( mu_ext )

        if (create_sc_branch):
            # p1 = self.p_usat[1][-1]
            dp = self.p_usat[1][-1] - self.p_usat[1][-2]
            dB = self.B_usat[1][-1] - self.B_usat[1][-2]
            dmu = self.mu_usat[1][-1] - self.mu_usat[1][-2]
            p_ext = self.p_usat[0][0] + self.p_usat_rel_max
            B_ext = self.B_usat[0][0] + dB / dp * self.p_usat_rel_max
            B_ext = max(B_ext, self.min_B_limit)
            mu_ext = self.mu_usat[0][0] + dmu / dp * self.p_usat_rel_max
            self.p_usat[0].append(p_ext)
            self.B_usat[0].append(B_ext)
            self.mu_usat[0].append(mu_ext)


    def checkConsistency_(self):
        assert(len(self.p_usat[-1]) > 1), "Unsaturated extra data must be specified for the highest Rs in PVTO"

        # check the saturated table
        for i in range(1, len(self.p_bub)):
            assert( self.p_bub[i] >= self.p_bub[i-1] ), "Bubble-point pressure saturated values must" +\
                                                       "be in ascending order"
            assert( self.Rs_sat[i] >= self.Rs_sat[i-1] ), "Rs must increase with bubble point pressure"
            assert( self.mu_sat[i] <= self.mu_sat[i-1] ), "Saturated viscosity must increase with bubble pressure"

        # check the under-saturated table
        for i in range(len(self.p_usat)):
            for j in range ( 1, len(self.p_usat[i]) ):
                assert(self.p_usat[i][j] >= self.p_usat[i][j-1] ), \
                    "Pressure must increase in unsaturated region"
                assert(self.B_usat[i][j] <= self.B_usat[i][j-1] ), \
                    "Volume facto must descrease with P in usat region"
                assert(self.mu_usat[i][j] >= self.mu_usat[i][j-1] ), \
                "Viscosity must increase with pressure in undersaturated conditions"

            # check for crossing of branches
            for j in range ( i+1, len(self.p_usat) ):
                assert ( not utils.curves_intersect(self.p_usat[i], self.p_usat[j],
                                              self.B_usat[i], self.B_usat[j]) ), \
                    "Volume factor under-saturated branches must not intersect"
                assert ( not utils.curves_intersect(self.p_usat[i], self.p_usat[j],
                                              self.mu_usat[i], self.mu_usat[j]) ), \
                    "Viscosity under-saturated branches must not intersect"
