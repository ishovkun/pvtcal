#!/usr/bin/env python


class PVTPhase:
  # Constructors and etc
  def __init__(self, Table_len : int = 0):
      pass

  def __init__(self, other : PVTPhase):
      pass

  def computeViscosity(self, vars, ip: int) -> float:
      pass


  def computeFugacity(xcp, pres, temp, ip, fug):
      pass

  def computeEnthalpy(vars, ip):
      pass

  def readTable(self, table, keyword):
      pass

  def getPressure(self, i):
      pass

  # double        getSurfaceDensity()     const { return dSurfaceDensity; }
  # double        getMolarWeight()        const { return dMolarWeight; }
  # const double* getXcpTbl(int number)   const { return &(vAllTables[XCP_SATURATED + number][0]);}
  # double        getFvf (int i)          const { return vAllTables[FVF_ORIGINAL][i]; }
  # double        getCompressibility ()   const { return dCompressibility; }
  # double        getViscosibility ()     const { return dViscosibility; }
  # double        getRatio (int i)        const { return vAllTables[RCP_ORIGINAL][i]; }
  # double        getViscosity (int i)    const { return vAllTables[VSC_ORIGINAL][i]; }
  # double        getViscosityRef ()    const { return dViscosityReference; }
