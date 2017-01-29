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


def thermal_expan(tdata, pdata, vt_p, outpath, T_accurate):
    filealpha = open(outpath + "/Alpha-%s.dat" %(str(T_accurate)), "w")
    filealpha.write("#T (K)    alpha (10^(-5) K^(-1))\n")
    alpha = np.zeros(shape=(len(pdata), len(tdata), 2))

    for nP, P in enumerate(pdata):
        filealpha.write("#P= %f GPa\n" % P)
        tck = interpolate.splrep(vt_p[nP][:, 0], vt_p[nP][:, 1])

        for nT, T in enumerate(tdata):
            # The first derivative
            alpha[nP][nT][0] = T
            alpha[nP][nT][1] = 100000.0 * interpolate.splev(T, tck, der=1) / vt_p[nP][nT][1]
            # Write to files
            if T == 0:
                filealpha.write("%7.1f%22.13f\n" % (0, 0))
            else:
                filealpha.write("%7.1f%22.13f\n" % (alpha[nP][nT][0], alpha[nP][nT][1]))

        filealpha.write("\n")

            # print "Alpha done..."
    filealpha.close()
    return alpha