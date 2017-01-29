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

import Constants
from scipy import integrate
from scipy import interpolate
from Constants import R


# Calculate Debye temperature
def debye_temp(num_formula_units, tdata, cvph_t, cvph_v, inppath, outpath, T_accurate):
    filedebye = open(outpath + "/DebyeT-%s.dat" %(str(T_accurate)), 'w')
    filedebye.write("# T/K      Debye temperature\n")

    for nV, P in enumerate(cvph_v[0][:, 0]):
        filedebye.write("# V=%18.13f\n" % cvph_v[0][:, 0][nV])
        for nT, t in enumerate(tdata):
            xlist = []
            ydifflist = []
            theta_d = 1.0
            while theta_d < 1000:
                xlist.append(theta_d)

                if t == 0:
                    t = 1

                x0 = theta_d / t
                yint, err = integrate.quad(Constants.y, 0.0, x0)
                cv_debye = 9 * R * x0 ** (-3) * yint

                # Compare C_V of the lattice vib. with C_V from the Debye model.
                ydifflist.append(cvph_t[nV][nT][1] / num_formula_units - cv_debye)
                theta_d += 10

            tck = interpolate.splrep(xlist, ydifflist)
            root = interpolate.sproot(tck, mest=100)
            if list(root):
                filedebye.write("%8.3f%20.1f\n" % (t, root[0]))

        filedebye.write("\n")
    filedebye.close

    return None
