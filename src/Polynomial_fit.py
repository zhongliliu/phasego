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
from ReadInput import indict
from Constants import RyperBohr3toGPa


if indict['Eos_Name'][0][0] == 'p' or indict['Eos_Name'][0][0] == 'P':
    order = int(indict['Eos_Name'][0][1:])
elif indict['Names_of_Strs'][0] == 'liquid':
    order = 1


# ========================== fitting using the selected eos ==========================
class Polyfit:
    def __init__(self, v=[], e=[]):
        self.order = order
        self.poly_pars = []
        self.v = np.array(v)
        self.e = np.array(e)

    def getparas(self):
        poly_pars = np.polyfit(self.v, self.e, self.order)
        self.poly_pars = poly_pars

        return poly_pars

    def getenergy(self, vol):
        p = np.poly1d(self.poly_pars)
        return p(vol)

    def getpressure(self, vol):
        dv = 0.0001
        p = np.poly1d(self.poly_pars)
        pressure = -RyperBohr3toGPa*(p(vol + dv) - p(vol - dv))/dv/2

        return pressure