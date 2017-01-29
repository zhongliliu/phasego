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
from ReadInput import indict
from scipy import interpolate


def getb_t_t(tdata, pv_t, outpath, T_accurate):
    if indict['If_Incl_Phonon'][0] == 'no':
        tdata = [0]

    fileb_t_t = open(outpath+"/B_T_T-%s.dat" %(str(T_accurate)), "w")

    # Number of V for interpolation
    nvolumes = 100

    # At fixed T
    b_t_t = np.zeros(shape=(len(tdata), nvolumes, 2))
    for nT, T in enumerate(tdata):
        fileb_t_t.write("#T= %f K" % T + '\n')

        tck = interpolate.splrep(pv_t[nT][:, 1][::-1], pv_t[nT][:, 0][::-1])
        for nvol, vol in enumerate(linspace(max(pv_t[nT][:, 1]), min(pv_t[nT][:, 1]), nvolumes)):
            b_t_t[nT][nvol][0] = interpolate.splev(vol, tck)
            b_t_t[nT][nvol][1] = - vol * interpolate.splev(vol, tck, der=1)
            fileb_t_t.write("%17.13f     %17.13f" % (b_t_t[nT][nvol][0], b_t_t[nT][nvol][1]) + "\n")

    # print "B_T_T finished"
    return b_t_t
