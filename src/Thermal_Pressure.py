"""
  Phasego -- Automatic calculation and plot of phase diagram

  Copyright (C) 2012-2016 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""

import numpy as np
from numpy import linspace


def thermal_pressure(tdata, free_v, pv_t, outpath, T_accurate):
    filethermalpressure_v = open(outpath + "/ThermalP_V-%s.dat" %(str(T_accurate)), "w")
    filethermalpressure_v.write("# T (K)     Thermal pressure (GPa)\n")
    nvolumes = 100
    thermalpressure_v = np.zeros(shape=(len(tdata), nvolumes, 2))

    # Extrplolate volumes and calculate thermal pressure by: P(T) - P(0)
    # at fixed V
    for nvol, vol in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
        filethermalpressure_v.write("# V= %17.13f\n" % vol)
        for nT in range(len(free_v)):
            thermalpressure_v[nT][nvol][0] = tdata[nT]
            thermalpressure_v[nT][nvol][1] = pv_t[nT][nvol][0] - pv_t[0][nvol][0]
            filethermalpressure_v.write("%7.1f%22.13f\n" % (tdata[nT], pv_t[nT][nvol][0] - pv_t[0][nvol][0]))

    # at fixed T.
    filethermalpressure_t = open(outpath + "/ThermalP_T-%s.dat" %(str(T_accurate)), "w")
    filethermalpressure_t.write("# V (a.u.^3)     Thermal pressure (GPa)\n")
    for nT, T in enumerate(tdata):
        filethermalpressure_t.write("# T= %7.1f\n" % T)
        for nvol, vol in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
            filethermalpressure_t.write("%17.13f%22.13f\n" % (vol, thermalpressure_v[nT][nvol][1]))

    filethermalpressure_v.close()
    filethermalpressure_t.close()
    # print "Thermal pressure done..."

    return thermalpressure_v
