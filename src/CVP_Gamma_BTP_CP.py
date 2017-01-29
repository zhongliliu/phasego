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
from scipy import interpolate


def cvp_gamma_btp_cp(tdata, pdata, free_v, gibbs_p, pv_t, vt_p, c_v_v, b_t_t, outpath, T_accurate):
    # Prepare files to write in
    filec_v_p = open(outpath + "/C_V_P-%s.dat" %(str(T_accurate)), 'w')
    filec_v_p.write("# T (K)      C_V_P (J K^-1 mol^-1)" + '\n')

    filegamma_t = open(outpath + "/gamma_P-%s.dat" %(str(T_accurate)), "w")
    filegamma_t.write("# T (K)       gamma" + '\n')
    fileb_s = open(outpath + "/B_S-%s.dat" %(str(T_accurate)), 'w')
    fileb_s.write("# T (K)      B_S (GPa)" + '\n')
    fileb_t_p = open(outpath + "/B_T_P-%s.dat" %(str(T_accurate)), 'w')
    fileb_t_p.write("# T (K)      B_T_P (GPa)" + '\n')
    filec_p = open(outpath + "/C_P-%s.dat" %(str(T_accurate)), 'w')
    filec_p.write("# T (K)      C_P (J K^-1 mol^-1)" + '\n')

    gamma = np.zeros(shape=(len(pdata), len(tdata), 2))
    b_s_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    b_t_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    c_p_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    c_v_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    nvolumes = 100
    for nP, P in enumerate(pdata):
        filegamma_t.write("# P = %f (GPa)" % P + '\n')
        fileb_s.write("# P = %f (GPa)" % P + '\n')
        fileb_t_p.write("# P = %f (GPa)" % P + '\n')
        filec_p.write("# P = %f (GPa)" % P + '\n')
        filec_v_p.write("# P = %f (GPa)" % P + '\n')

        # Scan T
        for nT, T in enumerate(tdata):

            if T == 0:
                T = 1

            pressure_data = np.array([pv_t[nT][nvol][0] for nvol in range(nvolumes)])
            volume_data = np.array([pv_t[nT][nvol][1] for nvol in range(nvolumes)])
            tckk = interpolate.splrep(pressure_data, volume_data)

            p_data = np.array([b_t_t[nT][nv][0] for nv in range(nvolumes)])
            bt_data = np.array([b_t_t[nT][nv][1] for nv in range(nvolumes)])
            tckl = interpolate.splrep(p_data, bt_data)
            v = interpolate.splev(float(P), tckk)
            bt = interpolate.splev(float(P), tckl)

            v_data = []
            cv_data = []
            n_p = len(free_v[0][:, 0]) - 1

            while n_p >= 0:
                v_data.append(free_v[0][:, 0][n_p])
                cv_data.append(c_v_v[n_p][nT][1])
                n_p -= 1

            # interpolation
            tckc = interpolate.splrep(v_data, cv_data)
            cv = interpolate.splev(v, tckc)

            tt = np.array([vt_p[nP][n_T][0] for n_T in range(len(tdata))])
            vv = np.array([vt_p[nP][n_T][1] for n_T in range(len(tdata))])
            tckt = interpolate.splrep(tt, vv)

            alp = interpolate.splev(T, tckt, der=1) / vt_p[nP][nT][1]

            # Gamma constant
            gamma[nP][nT][0] = T

            if cv == 0:
                cv = 1

            gamma[nP][nT][1] = 89.207089645992284 * alp * bt * v / cv
            filegamma_t.write("%7.1f%22.13f\n" % (gamma[nP][nT][0], gamma[nP][nT][1]))

            # B_S and C_P
            b_s_p[nP][nT][0] = T
            b_s_p[nP][nT][1] = bt * (1 + alp * gamma[nP][nT][1] * T)

            b_t_p[nP][nT][0] = T
            b_t_p[nP][nT][1] = bt

            c_p_p[nP][nT][0] = T
            c_p_p[nP][nT][1] = cv + 89.207089645992284 * alp ** 2. * bt * v * T
            c_v_p[nP][nT][0] = T
            c_v_p[nP][nT][1] = cv

            # Write to files
            fileb_t_p.write("%7.1f%22.13f\n" % (T, bt))
            filec_v_p.write("%7.1f%22.13f\n" % (T, cv))
            fileb_s.write("%7.1f%22.13f\n" % (T, b_s_p[nP][nT][1]))
            filec_p.write("%7.1f%22.13f\n" % (T, c_p_p[nP][nT][1]))

        gamma[nP][0][1] = 0

    fileentropy_p = open(outpath + "/Entropy_P-%s.dat" %(str(T_accurate)), 'w')
    fileentropy_p.write("# T (K)" + 4 * " " + "Entropy (J K^-1 mol^-1)\n")

    entropy_p = np.zeros(shape=(len(pdata), len(tdata), 2))
    for nP, P in enumerate(pdata):
        fileentropy_p.write("# P = %f (GPa)\n" % P)
        t_data = np.array([gibbs_p[nP][nT][0] for nT in range(len(tdata))])
        gibbs_data = np.array([gibbs_p[nP][nT][1] for nT in range(len(tdata))])
        tck = interpolate.splrep(t_data, gibbs_data)

        # Entropy at fixed P
        for nT, T in enumerate(tdata):
            entropy_p[nP][nT][0] = T
            entropy_p[nP][nT][1] = - Constants.RyperAtomperK2JperMolperK * interpolate.splev(T, tck, der=1)
            fileentropy_p.write("%7.1f%22.13f\n" % (entropy_p[nP][nT][0], entropy_p[nP][nT][1]))

    fileentropy_p.close()
    filec_p.close()
    filec_v_p.close()
    filegamma_t.close()
    fileb_t_p.close()
    fileb_s.close()

    return b_t_p, b_s_p, c_p_p, c_v_p, gamma

