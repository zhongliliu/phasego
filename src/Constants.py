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

# Some constants
ry = 13.6058
R = 8.314472                    # J / (mol*K)
h = 6.62606876e-34 * 2.99792458e10   # not hbar
kb = 1.3806504e-23              # J/K
e0 = 1.602176462e-19
NA = 6.02214199e23              # 1/mol

# Some transformation factors
RyperBohr3toGPa = 14710.531917404165
RyperAtomperK2JperMolperK = 131.22837394961041 * 10000

# Some functions
y = lambda z: z**4.0 * np.exp(z) / (np.exp(z)-1.0)**2.0
