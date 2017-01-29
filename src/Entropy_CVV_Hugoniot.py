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
import numpy as np
from numpy import linspace
from ReadInput import indict
from scipy import interpolate


def entropy_cvv_hugoniot(tdata, free_v, pv_t, helmholtz_v, outpath, T_accurate):
    fileentropy_v = open(outpath + "/Entropy_V-%s.dat" %(str(T_accurate)), 'w')
    fileentropy_v.write("# T (K)      Entropy (J K^-1 mol^-1)" + '\n')
    filefittede_v = open(outpath + "/fittedE_V-%s.dat" %(str(T_accurate)), 'w')
    filefittede_v.write("# T (K)      fittedE (Ry.)      P (GPa)" + '\n')
    filec_v_v = open(outpath + "/C_V_V-%s.dat" %(str(T_accurate)), 'w')
    filec_v_v.write("# T (K)      C_V_V (J K^-1 mol^-1)" + '\n')
    nvolumes = 100
    entropy_v = np.zeros(shape=(nvolumes, len(tdata), 2))
    fittede_v = np.zeros(shape=(nvolumes, len(tdata), 2))
    c_v_v = np.zeros(shape=(nvolumes, len(tdata), 2))


    # print V0

    # interpolation
    for nV, V in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
        fileentropy_v.write("# V = %f (a.u.^3)\n" % V)
        filefittede_v.write("# V = %f (a.u.^3)\n" % V)
        filec_v_v.write("# V = %f (a.u.^3)\n" % V)
        filefittedc_v_v = open(outpath + "/fittedC_V_V-%s.dat" %(str(T_accurate)), 'w')
        filefittedc_v_v.write("# T (K)      fittedC_V_V (J K^-1 mol^-1)\n")

        t_data = np.array([helmholtz_v[nV][nT][0] for nT in range(len(tdata))])
        helmholtz_data = np.array([helmholtz_v[nV][nT][1] for nT in range(len(tdata))])

        tck = interpolate.splrep(t_data, helmholtz_data)

        for nT, T in enumerate(tdata):
            # Entropy at fixed V
            entropy_v[nV][nT][0] = T
            entropy_v[nV][nT][1] = -Constants.RyperAtomperK2JperMolperK * interpolate.splev(T, tck, der=1)
            fileentropy_v.write("%7.1f%22.13f\n" % (entropy_v[nV][nT][0], entropy_v[nV][nT][1]))

            # interal energy at fixed V
            fittede_v[nV][nT][0] = T
            fittede_v[nV][nT][1] = helmholtz_data[nT] + T * entropy_v[nV][nT][1] / Constants.RyperAtomperK2JperMolperK
            filefittede_v.write("%7.1f%22.13f%22.13f\n" % (fittede_v[nV][nT][0], fittede_v[nV][nT][1], pv_t[nT][nV][0]))

            # C_V at fixed V
            c_v_v[nV][nT][0] = T
            c_v_v[nV][nT][1] = -T * Constants.RyperAtomperK2JperMolperK * interpolate.splev(T, tck, der=2)
            filec_v_v.write("%7.1f%22.13f\n" % (c_v_v[nV][nT][0], c_v_v[nV][nT][1]))



    if indict['If_Hugoniot'][0] == 'yes':
        pdata = []
        vdata = []

        for nV, V in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
            for nT, T in enumerate(tdata):
                if T == 300:
                    pdata.append(pv_t[nT][nV][0])
                    vdata.append(pv_t[nT][nV][1])

        tckpv = interpolate.splrep(pdata, vdata)
        v0 = interpolate.splev(0, tckpv, der=0)

        # For Hugoniot relation: p(V0-V) - 2(E-E0) = 0
        # Fix V: increasing T to find the solution of this equation
        vvdata = []
        eedata = []
        for nV, V in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
            for nT, T in enumerate(tdata):
                if T == 300:
                    vvdata.append(V)
                    eedata.append(fittede_v[nV][nT][1])

        tckev = interpolate.splrep(vvdata[::-1], eedata[::-1])
        e0 = interpolate.splev(v0, tckev, der=0)
        
        filehogniotptv = open(outpath + "/HugoniotPTV-%s.dat" %(str(T_accurate)), 'w')
        filehogniotptv.write("# Hugoniot P(GPa)-T(K)-V(Bohr^3)\n")

        for nV, V in enumerate(linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
            pvedata = []
            t_data = []
            pdata = []
            for nT, T in enumerate(tdata):
                pvedata.append(pv_t[nT][nV][0] * (v0 - V) - 2 * (fittede_v[nV][nT][1] - e0) * Constants.RyperBohr3toGPa)
                t_data.append(T)
                pdata.append(pv_t[nT][nV][0])

            # print T_data,pvedata

            tckpve = interpolate.splrep(t_data, pvedata)
            root = interpolate.sproot(tckpve, mest=100)

            if list(root):
                ht = root[0]
                hv = V
                tckpt = interpolate.splrep(t_data, pdata)
                hp = interpolate.splev(ht, tckpt, der=0)
                if hp >= 0:
                    filehogniotptv.write("%8.2f%12.2f%12.4f\n" % (hp, ht, hv))
                    # print HP,HT,HV

        filehogniotptv.close()

    # For C_V_V
    fittedc_v_v = np.zeros(shape=(len(tdata), nvolumes, 2))
    for nT, T in enumerate(tdata):
        v_data = []
        cv_data = []
        n_p = len(free_v[0][:, 0]) - 1
        while n_p >= 0:
            v_data.append(free_v[0][:, 0][n_p])
            cv_data.append(c_v_v[n_p][nT][1])
            n_p -= 1

        tck = interpolate.splrep(v_data, cv_data)
        for nvol, vol in enumerate(np.linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
            fittedc_v_v[nT][nvol][0] = vol
            fittedc_v_v[nT][nvol][1] = interpolate.splev(vol, tck)

    for nvol, vol in enumerate(np.linspace(max(free_v[0][:, 0]), min(free_v[0][:, 0]), nvolumes)):
        filefittedc_v_v.write("# V= %f (a.u.^3)\n" % vol)
        for nT, T in enumerate(tdata):
            filefittedc_v_v.write("%7.1f%22.13f\n" % (T, fittedc_v_v[nT][nvol][1]))

        # print "Entropy_V done..."
        #    print "C_V_V done..."

    fileentropy_v.close()
    filec_v_v.close()
    filefittedc_v_v.close()

    return c_v_v
