#!/usr/bin/env python

fname = "/home/ishovkun/sim/blackoil/1d/gprs.data"

import WordIterator
import KeywordIterator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kwd_it  = KeywordIterator.KeywordIterator(fname)

for kwd, value in kwd_it:
    if (kwd == "PVDG"):
        value = np.array([float(x) for x in value])
        value = value.reshape( int(len(value) / 3), 3 )
        df = pd.DataFrame(value, columns=["p", "Bg", "mug"])
        # plt.plot(df["p"], df["Bg"])
        # plt.plot(df["p"], 1./df["mug"])
        # plt.show()
    elif (kwd == "PVTO"):
        print(value)
        # value = np.array([float(x) for x in value])
        # value = value.reshape( int(len(value) / 3), 3 )
        # df = pd.DataFrame(value, columns=["p", "Bg", "mug"])
