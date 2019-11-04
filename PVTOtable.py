#!/usr/bin/env python
import numpy as np
log = np.log

def interp_lin(x1, x2, y1, y2, x):
    """
    simple linear interpolation
    y1 and y2 are the function values in the points
    x1 and x2, respectively
    x is the point where we want to evaluate the
    function
    """
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1

def interp_log(x1, x2, y1, y2, x):
    """
    simple logarithmic interpolation
    no checks
    """
    return np.exp( interp_lin(x1, x2, log(y1), log(y2), x) )

def curves_intersect(x1, x2, y1, y2):
    """
    return true if two curves y1 = y1(x1) and y2 = y2(x2)
    intersect
    x1, x2, y1, and y2 are all either lists or np arrays
    """
    assert(len(x1) == len(y1))
    assert(len(x2) == len(y2))
    # split curves into segments : y = a1 x + b1 and y = a2 x + b2
    # and check all pairs for intersection
    for i in range(1, len(x1)):
        # build segment line y = a1 x + b1
        assert(x1[i] != x1[i-1])
        a1 = (y1[i] - y1[i-1]) / (x1[i] - x1[i-1])
        b1 = y1[i] - a1 * x1[i]
        for j in range(1, len(x2)):
            # second line : y = a2 x + b2
            assert(x2[j] != x2[j-1])
            a2 = (y2[j] - y2[j-1]) / (x2[j] - x2[j-1])
            b2 = y2[j] - a2 * x2[j]

            # compute intersection
            if (abs(a1 - a2) < 1e-6):
                if (abs (b1 - b2) < 1e-6):
                    return True # conicide

            # intersection
            x = - (b2 - b1) / (a2 - a1)
            if ( x >= x1[i-1] and x <= x1[i] and x >= x2[j-1] and x <= x2[j]):
                return True
    return False


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
        self.makeStandardConditionsBranch_()
        self.extapolateUSatToMaximumRelativePressure_()
        self.checkConsistency_()

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

            self.Rs_sat.insert( 0, 0 )
            self.p_bub.insert (0, self.p_atm)  # atmoscpheric
            self.Rs_sat.insert(0, 0.0)        # no gas dissolved on surface
            self.B_sat.insert (0, 1.0)        # B on surface = 1 (no gas dissolved)
            self.mu_sat.insert(0, mu_atm)     # B on surface = 1 (no gas dissolved)

            # create usat branch;
            self.Rs_usat.insert( 0, 0.0)
            # self.p_usat.insert( 0, [ 0 ])
            self.p_usat.insert( 0, [ self.p_atm ])
            self.B_usat.insert( 0, [ 1.0 ])
            self.mu_usat.insert( 0, [ mu_atm ])
            # NOTE: the second point will be extrapolated later

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


    def extapolateUSatToMaximumRelativePressure_(self):
        """
        extrapolate to the maximum unsaturated pressure:
        I want tables for all rs span the same p_usat range
        """
        extrapolate_into_sc = False
        for i in range(len(self.p_usat)):
            # only a standard conditions branch previously
            # at this point can have a single value
            # remember that and treat it later
            if (len(self.p_usat[i]) < 2):
                extrapolate_into_sc = True
            else:
                if (self.p_usat[i][-1] < self.p_usat[i][0] + self.p_usat_rel_max):
                    # extrapolate linearly the unsaturated properties to the
                    # maximum relative pressure
                    p = self.p_usat[i]
                    B = self.B_usat[i]
                    mu = self.mu_usat[i]
                    B_ext = interp_lin( p[-2], p[-1], B[-2], B[-1], p[0] + self.p_usat_rel_max)
                    mu_ext = interp_lin( p[-2], p[-1], mu[-2], mu[-1], p[0] + self.p_usat_rel_max)
                    p.append( p[0] + self.p_usat_rel_max )
                    B.append( B_ext )
                    mu.append( mu_ext )

        # if (extrapolate_into_sc):
        #     # two closest branches to the SC branch
        #     p1 = self.p_usat[1][-1]
        #     p2 = self.p_usat[2][-1]
        #     rs1 = self.Rs_usat[1]
        #     rs2 = self.Rs_usat[2]
        #     B1 = self.B_usat[1][-1]
        #     B2 = self.B_usat[2][-1]
        #     mu1 = self.mu_usat[1][-1]
        #     mu2 = self.mu_usat[2][-1]
        #     p_ext = self.p_usat[0][0] + self.p_usat_rel_max
        #     # B_ext = interp_lin( rs1, rs2, B1, B2, self.Rs_usat[0] )
        #     # B_ext = interp_lin( p1, p2, B1, B2, p_ext )
        #     # mu_ext = interp_lin( rs1, rs2, mu1, mu2, self.Rs_usat[0] )
        #     B_ext = 1.5
            
        #     mu_ext = interp_log( p1, p2, mu1, mu2, p_ext )
        #     self.p_usat[0].append( p_ext )
        #     self.B_usat[0].append( B_ext )
        #     self.mu_usat[0].append( mu_ext )

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
                    "Pressure must decrease in unsaturated region"
                assert(self.B_usat[i][j] <= self.B_usat[i][j-1] ), \
                    "Volume facto must descrease with P in usat region"
                assert(self.mu_usat[i][j] >= self.mu_usat[i][j-1] ), \
                "Viscosity must increase with pressure in undersaturated conditions"

            # check for crossing of branches
            for j in range ( i+1, len(self.p_usat) ):
                assert ( not curves_intersect(self.p_usat[i], self.p_usat[j],
                                              self.B_usat[i], self.B_usat[j]) ), \
                    "Volume factor under-saturated branches must not intersect"
                assert ( not curves_intersect(self.p_usat[i], self.p_usat[j],
                                              self.mu_usat[i], self.mu_usat[j]) ), \
                    "Viscosity under-saturated branches must not intersect"
