#!/usr/bin/env python

fname = "/home/ishovkun/sim/blackoil/1d/gprs.data"

import KeywordIterator
from PVTOtable import PVTOtable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kwd_it  = KeywordIterator.KeywordIterator(fname)

for kwd, value in kwd_it:
    if (kwd == "PVDG"):
        value = np.array([float(x) for x in value])
        value = value.reshape( int(len(value) / 3), 3 )
        # df = pd.DataFrame(value, columns=["p", "Bg", "mug"])
        # plt.plot(df["p"], df["Bg"])
        # plt.plot(df["p"], 1./df["mug"])
        # plt.show()

    elif (kwd == "PVTO"):
        data = []
        for i in range(len(value)):
            data.append([float(x) for x in value[i]])
        pvto = PVTOtable(data)

        # df = pd.DataFrame(tsat, columns=["Rs", "p", "Bo", "muo"])
        # dfu = pd.DataFrame(tusat, columns=["Rs", "p", "Bo", "muo"])
        # plt.plot(df["p"], df["Bo"], 'bo-')
        # plt.plot(dfu["p"], dfu["Bo"], 'ko')
        # plt.show()
