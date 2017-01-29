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
from scipy import interpolate


# Transiform PV_T to VT_P
def pv_t2vt_p(tdata, pdata, free_v, pv_t, outpath, T_accurate):
    vt_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    for nT in range(len(free_v)):
        # interpolation
        tck = interpolate.splrep(pv_t[nT][:, 0], pv_t[nT][:, 1])
        for nP, P in enumerate(pdata):
            vt_p[nP][nT][0] = tdata[nT]
            vt_p[nP][nT][1] = interpolate.splev(float(P), tck)

    filevt_p = open(outpath + "/VT_P-%s.dat" %(str(T_accurate)), "w")
    filevt_p.write("#T (K)" + 6 * " " + "V (a.u.^3)\n")
    for nP, P in enumerate(pdata):
        filevt_p.write("# P= %f GPa\n" % P)
        for nT in range(len(free_v)):
            # Write to files
            filevt_p.write("%7.1f%22.13f\n" % (tdata[nT], vt_p[nP][nT][1]))

    filevt_p.close()
    # print "VT_P done..."
    return vt_p