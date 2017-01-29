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

import math
from ReadInput import indict
from scipy import interpolate
import matplotlib.pyplot as plt


def plot_pt(ptdatap, ptdatat, lowest):
    ifplot = True
    ifhasplotted = False
    plt.figure(figsize=(8, 6))

    for i in ptdatap.keys():
        if len(ptdatap[i]) >= 4:
            ifhasplotted = True
            k = i.split("---")
            times = (float(indict['Tdata_Write'][1]) - float(indict['Tdata_Write'][0])) / \
                    (float(indict['Pdata'][1]) - float(indict['Pdata'][0]))

            plocation = ptdatap[i][int(len(ptdatap[i]) / 2)]
            tlocation = ptdatat[i][int(len(ptdatat[i]) / 2)]

            tck = interpolate.splrep(ptdatat[i], ptdatap[i])
            rad = 6 / interpolate.splev(tlocation, tck, der=1) / times / 8

            angle = math.degrees(math.atan(rad))

            pshiftdist = 0.03 * (float(indict['Pdata'][1]) - float(indict['Pdata'][0]))

            deltap = pshiftdist * abs(math.sin(math.radians(abs(angle))))
            deltat = pshiftdist * abs(math.cos(math.radians(abs(angle)))) * times
            if angle >= 0:
                plt.text(plocation - deltap, tlocation + deltat, k[0], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')
                plt.text(plocation + deltap, tlocation - deltat, k[1], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')

            if angle < 0:
                plt.text(plocation - deltap, tlocation - deltat, k[0], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')
                plt.text(plocation + deltap, tlocation + deltat, k[1], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')

            plt.plot(ptdatap[i], ptdatat[i], 'o-', linewidth=1.0)

        else:
            print "Warning: too few data to plot by scanning temperature!"
            ifplot = False
            pass

    if ifplot:
        plt.xlabel("Pressure (GPa)")
        plt.ylabel("Temperature (K)")
        plt.xlim(float(indict['Pdata'][0]), float(indict['Pdata'][1]))
        plt.ylim(float(indict['Tdata_Write'][0]), float(indict['Tdata_Write'][1]))

        ppos = (float(indict['Pdata'][0]) + float(indict['Pdata'][1])) / 2
        tpos = (float(indict['Tdata_Write'][0]) + float(indict['Tdata_Write'][1])) / 2

        if not ifhasplotted:
            plt.text(ppos, tpos, lowest)

        plt.savefig('Phase-PT/P-T.eps', dpi=600)
        plt.savefig('Phase-PT/P-T.pdf', dpi=600)
        plt.show()

    return None


def plot_tp(ptdatap, ptdatat, lowest):
    ifplot = True
    ifhasplotted = False
    plt.figure(figsize=(8, 6))
    for i in ptdatap.keys():
        if len(ptdatap[i]) >= 4:
            ifhasplotted = True
            k = i.split("---")
            times = (float(indict['Tdata_Write'][1]) - float(indict['Tdata_Write'][0])) / \
                    (float(indict['Pdata'][1]) - float(indict['Pdata'][0]))

            plocation = ptdatap[i][int(len(ptdatap[i]) / 2)]
            tlocation = ptdatat[i][int(len(ptdatat[i]) / 2)]

            tck = interpolate.splrep(ptdatap[i], ptdatat[i])
            rad = 6 * interpolate.splev(plocation, tck, der=1) / times / 8

            angle = math.degrees(math.atan(rad))

            pshiftdist = 0.03 * (float(indict['Pdata'][1]) - float(indict['Pdata'][0]))

            deltap = pshiftdist * abs(math.sin(math.radians(abs(angle))))
            deltat = pshiftdist * abs(math.cos(math.radians(abs(angle)))) * times

            if angle >= 0:
                plt.text(plocation + deltap, tlocation - deltat, k[0], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')
                plt.text(plocation - deltap, tlocation + deltat, k[1], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')

            if angle < 0:
                plt.text(plocation - deltap, tlocation - deltat, k[0], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')
                plt.text(plocation + deltap, tlocation + deltat, k[1], rotation=angle,
                         horizontalalignment='center', verticalalignment='center')

            plt.plot(ptdatap[i], ptdatat[i], 'o-', linewidth=1.0)

        else:
            print "Warning: too few data to plot by scanning pressure!"
            print
            ifplot = False
            pass

    if ifplot:
        plt.xlabel("Pressure (GPa)")
        plt.ylabel("Temperature (K)")
        plt.xlim(float(indict['Pdata'][0]), float(indict['Pdata'][1]))
        plt.ylim(float(indict['Tdata_Write'][0]), float(indict['Tdata_Write'][1]))

        ppos = (float(indict['Pdata'][0]) + float(indict['Pdata'][1])) / 2
        tpos = (float(indict['Tdata_Write'][0]) + float(indict['Tdata_Write'][1])) / 2

        if not ifhasplotted:
            plt.text(ppos, tpos, lowest)

        plt.savefig('Phase-PT/T-P.eps', dpi=600)
        plt.savefig('Phase-PT/T-P.pdf', dpi=600)
        plt.show()

    return None