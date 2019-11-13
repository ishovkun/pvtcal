#!/usr/bin/env python

fname = "examples/gprs.data"

import KeywordIterator
from PVTOtable import PVTOtable
from PVDGtable import PVDGtable

import numpy as np
import random
import matplotlib.pyplot as plt
plt.style.use('dark_background')

kwd_it  = KeywordIterator.KeywordIterator(fname)

for kwd, value in kwd_it:
    if (kwd == "PVDG"):
        value = np.array([float(x) for x in value])
        # value = value.reshape( int(len(value) / 3), 3 )
        # exit(0)
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=[13, 5])
        pvdg = PVDGtable(value)
        ax1.plot(pvdg.p, pvdg.B)
        ax2.plot(pvdg.p, pvdg.mu)

        # sprinkle some random points
        for i in range(5):
            x = random.random()
            p = (max(pvdg.p) - min(pvdg.p)) * x + min(pvdg.p)
            ax1.plot( p, pvdg.getVolumeFactor(p), "ro" )
            ax2.plot( p, pvdg.getViscosity(p), "ro" )

        ax1.set_xlabel("Pressure [bar]")
        ax1.set_ylabel("Gas volume factor B$_g$")
        ax2.set_xlabel("Pressure [bar]")
        ax2.set_ylabel("Gas viscosity [cP]")
        # exit(0)

    elif (kwd == "PVTO"):
        data = []
        for i in range(len(value)):
            data.append([float(x) for x in value[i]])
        pvto = PVTOtable(data)

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=[13, 5])
        
        # volume factor
        ax1.plot(pvto.p_bub, pvto.B_sat, "bo-")
        for i in range(len(pvto.Rs_usat)):
            ax1.plot( pvto.p_usat[i] , pvto.B_usat[i] , 'o-')

        # plot viscosity
        ax2.plot(pvto.p_bub, pvto.mu_sat, "o-")
        for i in range(len(pvto.Rs_usat)):
            ax2.plot( pvto.p_usat[i] , pvto.mu_usat[i], "ro-" )

        ax1.set_xlabel("Pressure [bar]")
        ax1.set_ylabel("Oil volume factor B$_g$")
        ax2.set_xlabel("Pressure [bar]")
        ax2.set_ylabel("Oil viscosity [cP]")
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50,500,640, 545)

        # plot in relative pressure coordinates
        # fig = plt.figure()
        # plt.plot(np.zeros(len(pvto.p_bub)), pvto.B_sat, "bo-")
        # for i in range(len(pvto.Rs_usat)):
        #     plt.plot( np.array(pvto.p_usat[i]) - pvto.p_usat[i][0] , pvto.B_usat[i] , 'o-')

        # let's fucking build more usat branches
        for i in range(len(pvto.Rs_sat)):
            if(pvto.Rs_sat[i] not in pvto.Rs_usat):
                # get a couple of pressures just to
                # know what range we interpolate
                # print(pvto.Rs_sat[i], pvto.Rs_usat)
                p1 = pvto.p_bub[i] + 2
                p2 = pvto.p_bub[i] + 400
                B1 = pvto.getVolumeFactor(pvto.Rs_sat[i] + 1e-4, p1)
                B2 = pvto.getVolumeFactor(pvto.Rs_sat[i], p2)
                ax1.plot([pvto.p_bub[i], p1, p2], [pvto.B_sat[i], B1, B2], "g*-")
                mu1 = pvto.getViscosity(pvto.Rs_sat[i], p1)
                mu2 = pvto.getViscosity(pvto.Rs_sat[i], p2)
                ax2.plot([pvto.p_bub[i], p1, p2], [pvto.mu_sat[i], mu1, mu2], "g*-")


# if (__name__ == "__main__"):
mngr = plt.get_current_fig_manager()
mngr.window.setGeometry(50,100,640, 545)
plt.show()
