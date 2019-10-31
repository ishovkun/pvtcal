#!/usr/bin/env python

fname = "/home/ishovkun/sim/Blackoil/1d/gprs.data"

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
        # print(value.shape)
        # df = pd.DataFrame(value, columns=["p", "Bg", "mug"])
        # plt.plot(df["p"], df["Bg"])
        # plt.plot(df["p"], 1./df["mug"])
        # plt.show()
    elif (kwd == "PVTO"):
        # print(value)
        tsat = []
        tusat = []
        for row in value:
            if (len(row) == 4): # saturated
                tsat.append([float(x) for x in row])
                # print(row)
            else:               # unsaturated
                row = [float(x) for x in row]
                tsat.append(row[:4])
                tusat.append(row[:4])
                for i in range(4, len(row), 4):
                    line = [tsat[-1][0]] + row[i:i+4]
                    tusat.append(line)
                pass

        tsat = np.array(tsat)
        tusat = np.array(tusat)
        df = pd.DataFrame(tsat, columns=["Rs", "p", "Bo", "muo"])
        dfu = pd.DataFrame(tusat, columns=["Rs", "p", "Bo", "muo"])
        plt.plot(df["p"], df["Bo"])
        plt.plot(dfu["p"], dfu["Bo"], 'ko')
        plt.show()
