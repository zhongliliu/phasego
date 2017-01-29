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
import Constants, warnings
from EOSfit import EOSfit
from Polynomial_fit import Polyfit
from numpy import linspace
from ReadInput import indict
from scipy import interpolate

warnings.filterwarnings("ignore")


# Calculate Gibbs free energy: F + PV
def gibbsfree_en(tdata, pdata, free_v, outpath, T_accurate):
    # Prepare output files
    filehelmfree_t = open(outpath + "/HelmFreeE_T-%s.dat" %(str(T_accurate)), "w")
    filepvt = open(outpath + "/PV_T-%s.dat" %(str(T_accurate)), "w")
    filehelmfree_v = open(outpath + "/HelmFreeE_V-%s.dat" %(str(T_accurate)), "w")
    filehelmfree_v.write("# T (K)" + 15 * " " + "Helmholtz Free Energy (Ry.)\n")
    filegibbsfree = open(outpath + "/G_T-%s.dat" %(str(T_accurate)), "w")
    filehfree = open(outpath + "/FittedHelmFreeE_T-%s.dat" %(str(T_accurate)), "w")
    filehelmfree_t.write("# V (a.u.^3)  Helmholtz Free Energy (Ry.)\n")
    filepvt.write("# P (GPa)" + 13 * " " + "V (a.u.^3)\n")
    filegibbsfree.write("# P (GPa)" + 13 * " " + "Gibbs free energy (Ry.)\n")
    filehfree.write("# V (a.u.^3)" + 10 * " " + "Helmholtz Free Energy (Ry.)\n")
    fileenthalpy = open(outpath + "/Enthalpy-%s.dat" %(str(T_accurate)), "w")
    fileenthalpy.write("# P (GPa)" + 6 * " " + "H (Ry.)\n")

    # Number of volumes for interpolation
    nvolumes = 100

    # Prepare numpy array to contain P,V,G,T data
    pv_t = np.zeros(shape=(len(tdata), nvolumes, 2))
    gibbst = np.zeros(shape=(len(tdata), nvolumes, 2))
    gibbs_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    helmholtz_v = np.zeros(shape=(nvolumes, len(tdata), 2))

    # Calculate G at each T
    for nT in range(len(free_v)):
        filehelmfree_t.write("# T= %f K\n" % tdata[nT])
        filepvt.write("# T= %f K\n" % tdata[nT])
        filehfree.write("# T= %f K\n" % tdata[nT])

        for nP in range(len(free_v[nT])):
            filehelmfree_t.write("%9.4f %16.9f\n" % (free_v[nT][nP][0], free_v[nT][nP][1]))

        # Fitting Polynomial or EOS to F-V data before Gibbs free energy calculation
        if indict['Eos_Name'][0][0] == 'p' or indict['Eos_Name'][0][0] == 'P':
            func_fit = Polyfit(v=free_v[nT][:, 0], e=free_v[nT][:, 1])
            func_fit.getparas()
        else:
            func_fit = EOSfit(v=free_v[nT][:, 0], e=free_v[nT][:, 1])
            func_fit.getparas()

        # Extrapolation of volumes for each single structure
        for nvol, vol in enumerate(linspace(max(free_v[nT][:, 0]), min(free_v[nT][:, 0]), nvolumes)):
            filepvt.write("%17.13f%22.13f\n" % (func_fit.getpressure(vol), vol))
            filehfree.write("%17.13f%22.13f\n" % (vol, func_fit.getenergy(vol)))
            helmholtz_v[nvol][nT][0] = tdata[nT]
            helmholtz_v[nvol][nT][1] = func_fit.getenergy(vol)
            pv_t[nT][nvol][0] = func_fit.getpressure(vol)
            pv_t[nT][nvol][1] = vol

            # Obtain G by: F + PV
            gibbst[nT][nvol][0] = func_fit.getpressure(vol)
            gibbst[nT][nvol][1] = func_fit.getenergy(vol) + func_fit.getpressure(vol) * vol / Constants.RyperBohr3toGPa

    # Helmholtz free energy at each fixed volume
    for nvol, vol in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
        filehelmfree_v.write("# V= %17.13f\n" % vol)
        for nT, T in enumerate(tdata):
            filehelmfree_v.write("%8.2f%22.13f\n" % (helmholtz_v[nvol][nT][0], helmholtz_v[nvol][nT][1]))

    gibbs_t = np.zeros(shape=(len(tdata), len(pdata), 2))
    for nT in range(len(free_v)):
        # print gibbst[nT][:, 0], gibbst[nT][:, 1]
        tck = interpolate.splrep(gibbst[nT][:, 0], gibbst[nT][:, 1])

        for nP in range(len(pdata)):
            # Gibbs free energy at each fixed P
            gibbs_p[nP][nT][0] = tdata[nT]
            gibbs_p[nP][nT][1] = interpolate.splev(float(pdata[nP]), tck)

            # Gibbs free energy at each fixed T
            gibbs_t[nT][nP][0] = float(pdata[nP])
            gibbs_t[nT][nP][1] = gibbs_p[nP][nT][1]

    # Build big Gibbs free energy dictionary
    gibbsdict_p = {}
    for nP in range(len(pdata)):
        gibbsdict_p[pdata[nP]] = gibbs_p[nP]

    gibbsdict_t = {}
    for nT in range(len(tdata)):
        gibbsdict_t[tdata[nT]] = gibbs_t[nT]

    # Write Gibbs data into files
    filegibbsfree_p = open(outpath + "/G_P-%s.dat" %(str(T_accurate)), "w")
    filegibbsfree_p2 = open(outpath + "/G_P2-%s.dat" %(str(T_accurate)), "w")
    filegibbsfree_p.write("# T (K)     Gibbs free energy (Ry.)\n")
    filegibbsfree_p2.write("# T (K)     Gibbs free energy (Ry.)\n")

    # Write enthalpy data into files
    for nP in range(len(pdata)):
        fileenthalpy.write("%9.2f%22.13f\n" % (pdata[nP], gibbs_p[nP][0][1]))

    # Write Gibbs data into files
    for nP in range(len(pdata)):
        filegibbsfree_p.write("# P= %f GPa\n" % pdata[nP])
        for nT in range(len(free_v)):
            filegibbsfree_p.write("%7.1f%22.13f\n" % (tdata[nT], gibbs_p[nP][nT][1]))

    for nP in range(len(pdata)):
        for nT in range(len(free_v)):
            filegibbsfree_p2.write("%7.1f%22.13f\n" % (pdata[nP], gibbs_p[nP][nT][1]))

    # Write Gibbs data into files
    for nT in range(len(tdata)):
        filegibbsfree.write("# T= %f K\n" % tdata[nT])
        for nP in range(len(pdata)):
            filegibbsfree.write("%7.1f%22.13f\n" % (pdata[nP], gibbs_p[nP][nT][1]))

        # print "HelmholtzFreeEnergy_T done..."
        #    print "PV_T done..."
        #    print "GibbsFreeEnergy_T done..."
        #    print "GibbsFreeEnergy_V done..."
        #    print "FittedHelmholtzFreeEnergy done..."
        #    print "Enthalpy done..."
        #    print "GibbsFreeEnergy_P done..."

    filehelmfree_v.close()
    filehelmfree_t.close()
    filepvt.close()
    filegibbsfree.close()
    filehfree.close()
    fileenthalpy.close()
    filegibbsfree_p.close()
    return pv_t, gibbs_p, helmholtz_v, gibbsdict_p, gibbsdict_t
