#!/usr/bin/env python
import numpy as np

class PVTOtable:
    tsat = np.array(0)
    tusat = np.array(0)

    def __init__(self, data):
        """
        data is a list of lists of float values
        parsing is happening right here
        """
        tsat = []
        tusat = []
        for row in data:
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

        self.tsat = np.array(tsat)
        self.tusat = np.array(tusat)
