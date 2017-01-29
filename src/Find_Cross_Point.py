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

from scipy import interpolate


def find_cross_point(y1, y2):
    ydiff = y2 - y1

    xlist = [i[0] for i in y1]
    ydifflist = [i[1] for i in ydiff]

    tck = interpolate.splrep(xlist, ydifflist)
    root = interpolate.sproot(tck, mest=100)

    return root
