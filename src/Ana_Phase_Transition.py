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

import os
import numpy as np
from ReadInput import indict
from Plot_Matplotlib import plot_pt, plot_tp
from Find_Cross_Point import find_cross_point


def ana_phasetran(bigdict_p, bigdict_t):
    tdata = np.arange(int(indict['Tdata_Write'][0]), int(indict['Tdata_Write'][1]) + int(indict['Tdata_Write'][2]),
                      int(indict['Tdata_Write'][2]))
    pdata = np.arange(int(indict['Pdata'][0]), int(indict['Pdata'][1]) + int(indict['Pdata'][2]),
                      int(indict['Pdata'][2]))

    # if Phase-PT dir does not exist, then mkdir Phase-PT
    if not os.path.isdir("%s" % "Phase-PT"):
        os.mkdir("%s" % "Phase-PT")

    # Scan pressure at fixed temperature
    ldict = {}
    filept = open("Phase-PT/P-T.dat", "w")
    filept.write("# at fixed T/K    |    phase1    |   transits to   |     phase2   |   when P increases to P/GPa\n")

    p_cross = []
    ptdatap = {}
    ptdatat = {}
    for nt, t in enumerate(tdata):
        for phasename in bigdict_t.keys():
            ldict[phasename] = [i[1] for i in bigdict_t[phasename][t]]

        phasecross = []
        for nP in range(len(pdata)):
            ptmplist = []
            for phasename in bigdict_t.keys():
                ptmplist.append(ldict[phasename][nP])

            for phasename in bigdict_t.keys():
                if min(ptmplist) == ldict[phasename][nP] and phasename not in phasecross:
                    # print Pdata[nP],phasename
                    phasecross.append(phasename)

        if len(phasecross) >= 2:
            for i in range(len(phasecross) - 1):

                p_cross = find_cross_point(bigdict_t[phasecross[i]][t], bigdict_t[phasecross[i + 1]][t])
                if phasecross[i] + "---" + phasecross[i + 1] in ptdatap:
                    ptdatap[phasecross[i] + "---" + phasecross[i + 1]].append(p_cross[0])
                    ptdatat[phasecross[i] + "---" + phasecross[i + 1]].append(t)
                else:
                    ptdatap[phasecross[i] + "---" + phasecross[i + 1]] = [p_cross[0]]
                    ptdatat[phasecross[i] + "---" + phasecross[i + 1]] = [t]

        if indict['If_Incl_Phonon'][0] == 'no':
            break

    for i in ptdatap.keys():
        k = i.split("---")
        for j in range(len(ptdatap[i])):
            # print "T=",t,phase1,"===>",phase2,"at:",p
            filept.write("%14.2f %14s %18s %14s %18.2f\n" % (ptdatat[i][j], k[0], '--------->', k[1], ptdatap[i][j]))

        filept.write("\n")

    filept.close()

    lowest = ''
    if not list(p_cross):
        lowe = []
        for phasename in ldict.keys():
            lowe.append(ldict[phasename][0])
        for phasename in ldict.keys():
            if min(lowe) == ldict[phasename][0]:
                lowest = phasename

        print lowest, "is the most stable structure."
        print
        print "At fixed T: No phase transition is found\n" + " " * 12 + "in the predefined PT regime."
        print

    if indict['If_Plot'][0] == 'yes' and indict['If_Incl_Phonon'][0] == 'yes':
        plot_pt(ptdatap, ptdatat, lowest)

    # Scan temperature at fixed pressure
    lowest = ''
    if indict['If_Incl_Phonon'][0] == 'yes':
        ldict = {}
        filetp = open("Phase-PT/T-P.dat", "w")
        filetp.write(
            "# at fixed P/GPa   |     phase1   |    transits to    |    phase2   |   when T increases to T/K\n")

        t_cross = []
        ptdatap = {}
        ptdatat = {}
        for p in pdata:
            for phasename in bigdict_p.keys():
                ldict[phasename] = [i[1] for i in bigdict_p[phasename][p]]

            phasecross = []
            for nT in range(len(tdata)):
                ttmplist = []
                for phasename in bigdict_p.keys():
                    ttmplist.append(ldict[phasename][nT])

                for phasename in bigdict_p.keys():
                    if min(ttmplist) == ldict[phasename][nT] and phasename not in phasecross:
                        # print Tdata[nT],phasename
                        phasecross.append(phasename)

            if len(phasecross) >= 2:
                for i in range(len(phasecross) - 1):
                    t_cross = find_cross_point(bigdict_p[phasecross[i]][p], bigdict_p[phasecross[i + 1]][p])
                    if phasecross[i] + "---" + phasecross[i + 1] in ptdatap:
                        ptdatap[phasecross[i] + "---" + phasecross[i + 1]].append(p)
                        ptdatat[phasecross[i] + "---" + phasecross[i + 1]].append(t_cross[0])
                    else:
                        ptdatap[phasecross[i] + "---" + phasecross[i + 1]] = [p]
                        ptdatat[phasecross[i] + "---" + phasecross[i + 1]] = [t_cross[0]]

        if not list(t_cross):
            lowe = []
            for phasename in ldict.keys():
                lowe.append(ldict[phasename][0])
            for phasename in ldict.keys():
                if min(lowe) == ldict[phasename][0]:
                    lowest = phasename

            print lowest, "is the most stable structure."
            print
            print "At fixed P: No phase transition is found\n" + " " * 12 + "in the predefined PT regime."
            print

    for i in ptdatap.keys():
        k = i.split("---")
        for j in range(len(ptdatap[i])):
            # print "P=",p,phase1,"===>",phase2,"at:",t
            filetp.write("%16.2f %14s %18s %14s %18.2f\n" % (ptdatap[i][j], k[0], '--------->', k[1], ptdatat[i][j]))

        filetp.write("\n")

    filetp.close()

    if indict['If_Plot'][0] == 'yes' and indict['If_Incl_Phonon'][0] == 'yes':
        plot_tp(ptdatap, ptdatat, lowest)

    return None