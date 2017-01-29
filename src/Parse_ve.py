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


# Parse ve file
# ve file can have any number of blank lines and comment lines
def parse_ve(vefile):
    ve_file = open(vefile, "r")
    vvee = {}

    lines = ve_file.readlines()

    for l in lines:
        line = l.strip().split()
        # The lines contained no numerical data in the first two columns are ignored
        if len(line) >= 2:
            try:
                line0 = float(line[0])
                line1 = float(line[1])
            except:
                pass

            else:
                vvee[line0] = line1

    vsort = vvee.keys()
    vsort.sort()
    vsort.reverse()

    if indict['Units_VE'][0] == 'Bohr3':
        v_conv = 1.0
    elif indict['Units_VE'][0] == 'A3':
        v_conv = 1 / 0.529177**3.0

    if indict['Units_VE'][1] == 'Ry':
        e_conv = 1.0
    elif indict['Units_VE'][1] == 'Hartree':
        e_conv = 2.0
    elif indict['Units_VE'][1] == 'eV':
        e_conv = 2 / 27.21138

    ve = np.array([[i * v_conv, vvee[i] * e_conv] for i in vsort])

    ve_file.close()

    return ve

if __name__ == '__main__':
    ve = parse_ve("ve-0")
    print ve
