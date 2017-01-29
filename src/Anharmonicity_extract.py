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
from numpy import linspace

# Extracting anharmonic data from the dos
def anharm_extract1(alpha_anh, gdict_p, gdict_t, alpha, free_v, pv_t, b_t_t, gibbsdict_t, pdata, tdata_read,
                    tdata_write, gibbsdict_p, T_accurate, fileb_t_t, filepvt_anh, filegibbsfree):

    g_p = np.zeros(shape=(len(pdata), len(tdata_read), 2))

    nvolumes = 100

    for nP in range(len(pdata)):
        for nT in range(len(tdata_write)):
            if tdata_write[nT] == T_accurate:
                alpha_anh[nP][list(tdata_read).index(T_accurate)][0] = alpha[nP][:, 0][nT]
                alpha_anh[nP][list(tdata_read).index(T_accurate)][1] = alpha[nP][:, 1][nT]
                g_p[nP][list(tdata_read).index(T_accurate)][0] = gibbsdict_p[pdata[nP]][:, 0][nT]
                g_p[nP][list(tdata_read).index(T_accurate)][1] = gibbsdict_p[pdata[nP]][:, 1][nT]

    for nT in range(len(tdata_write)):
        if tdata_write[nT] == T_accurate:
            fileb_t_t.write("#T= %f K" % T_accurate + '\n')
            filepvt_anh.write("# T= %f K\n" % T_accurate)
            for nvol, vol in enumerate(linspace(max(free_v[nT][:, 0]), min(free_v[nT][:, 0]), nvolumes)):
                filepvt_anh.write("%17.13f%22.13f\n" % (pv_t[nT][nvol][0], pv_t[nT][nvol][1]))
                fileb_t_t.write("%17.13f     %17.13f" % (b_t_t[nT][nvol][0], b_t_t[nT][nvol][1]) + "\n")

        gdict_p[pdata[nP]] = g_p

    gdict_t[T_accurate] = gibbsdict_t[T_accurate]

    for nT in range(len(tdata_write)):
        if tdata_write[nT] == T_accurate:
            filegibbsfree.write("# T= %f K\n" % tdata_write[nT])
            for nP in range(len(pdata)):
                filegibbsfree.write("%7.1f%22.13f\n" % (gibbsdict_t[T_accurate][nP][0], gibbsdict_t[T_accurate][nP][1]))

    return alpha_anh, gdict_p, gdict_t

def anharm_extract2(alpha_anh, pdata, tdata_read, tdata_write, dict_cvph_t, dict_sph_t, dict_eph_t, dict_c_v_v,
                    dict_free_v, dict_helmholtz_v, dict_thermalpressure_v, dict_vt_p, gibbsdict_p, b_t_p, b_s_p,
                    c_p_p, c_v_p, gamma, ve, Num_Formula_Units_Comp, free_v, outpath, filealpha):

    filehelmfree_t = open(outpath + "/HelmFreeE_T_anh.dat", "w")
    filehelmfree_t.write("# V (a.u.^3)  Helmholtz Free Energy (Ry.)\n")
    filehelmfree_v = open(outpath + "/HelmFreeE_V_anh.dat", "w")
    filehelmfree_v.write("# T (K)" + 15 * " " + "Helmholtz Free Energy (Ry.)\n")

    filegibbsfree_p = open(outpath + "/G_P_anh.dat", "w")
    filegibbsfree_p.write("# T (K)     Gibbs free energy (Ry.)\n")
    fileb_t_p = open(outpath + "/B_T_P_anh.dat", 'w')
    fileb_t_p.write("# T (K)      B_T_P (GPa)" + '\n')
    fileb_s = open(outpath + "/B_S_anh.dat", 'w')
    fileb_s.write("# T (K)      B_S (GPa)" + '\n')
    filec_p = open(outpath + "/C_P_anh.dat", 'w')
    filec_p.write("# T (K)      C_P (J K^-1 mol^-1)" + '\n')

    filec_v = open(outpath + "/C_V_P_anh.dat", 'w')
    filec_v.write("# T (K)      C_V (J K^-1 mol^-1)" + '\n')
    filegamma_t = open(outpath + "/gamma_P_anh.dat", "w")
    filegamma_t.write("# T (K)       gamma" + '\n')

    cvdirectfile = open(outpath + "/C_V-ph-direct_anh.dat", "w")
    cvdirectfile.write("# T (K)      C-ph_V_V (J K^-1 mol^-1)" + '\n')
    sdirectfile = open(outpath + "/Entropy-ph_V_anh.dat", "w")
    sdirectfile.write("# T (K)      Entropy-ph (J K^-1 mol^-1)" + '\n')
    eintfile = open(outpath + "/E-ph_V_anh.dat", "w")
    eintfile.write("# T (K)       E-ph (Ry.)" + '\n')

    filec_v_v = open(outpath + "/C_V_V_anh.dat", 'w')
    filec_v_v.write("# T (K)      C_V_V (J K^-1 mol^-1)" + '\n')
    filethermalpressure_v = open(outpath + "/ThermalP_V_anh.dat", "w")
    filethermalpressure_v.write("# T (K)     Thermal pressure (GPa)\n")
    filevt_p = open(outpath + "/VT_P_anh.dat", "w")
    filevt_p.write("#T (K)" + 6 * " " + "V (a.u.^3)\n")

    nvolumes = 100

    for nP, P in enumerate(pdata):

        fileb_t_p.write("# P = %f (GPa)" % pdata[nP] + '\n')
        fileb_s.write("# P = %f (GPa)" % P + '\n')
        filec_p.write("# P = %f (GPa)" % P + '\n')
        filec_v.write("# P = %f (GPa)" % P + '\n')
        filegamma_t.write("# P = %f (GPa)" % P + '\n')
        filevt_p.write("# P= %f GPa\n" % P)
        filealpha.write("#P= %f GPa\n" % P)
        filegibbsfree_p.write("# P= %f GPa\n" % P)

        for nT, T in enumerate(tdata_read):
            # Write to files
            if T == 0:
                filealpha.write("%7.1f%22.13f\n" % (0, 0))
            else:
                filealpha.write("%7.1f%22.13f\n" % (alpha_anh[nP][nT][0], alpha_anh[nP][nT][1]))
        filealpha.write("\n")

        for T1 in tdata_read:
            for nT2, T2 in enumerate(tdata_write):
                if T2 == T1:
                    filegibbsfree_p.write("%7.1f%22.13f\n" % (gibbsdict_p[pdata[nP]][:, 0][nT2], gibbsdict_p[pdata[nP]][:, 1][nT2]))
                    fileb_t_p.write("%7.1f%22.13f\n" % (b_t_p[nP][nT2][0], b_t_p[nP][nT2][1]))
                    fileb_s.write("%7.1f%22.13f\n" % (b_s_p[nP][nT2][0], b_s_p[nP][nT2][1]))
                    filec_p.write("%7.1f%22.13f\n" % (c_p_p[nP][nT2][0], c_p_p[nP][nT2][1]))
                    filec_v.write("%7.1f%22.13f\n" % (c_v_p[nP][nT2][0], c_v_p[nP][nT2][1]))
                    filegamma_t.write("%7.1f%22.13f\n" % (gamma[nP][nT2][0], gamma[nP][nT2][1]))
                    filevt_p.write("%7.1f%22.13f\n" % (dict_vt_p[T1][nP][nT2][0], dict_vt_p[T1][nP][nT2][1]))

    for n in range(len(ve[0])):
        cvdirectfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / Num_Formula_Units_Comp))
        sdirectfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / Num_Formula_Units_Comp))
        eintfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / Num_Formula_Units_Comp))

        for T1 in tdata_read:
            for nT2, T2 in enumerate(tdata_write):
                if T2 == T1:
                    cvdirectfile.write("%7.1f%22.13f\n" % (dict_cvph_t[T1][n][nT2][0], dict_cvph_t[T1][n][nT2][1]))
                    sdirectfile.write("%7.1f%22.13f\n" % (dict_sph_t[T1][n][nT2][0], dict_sph_t[T1][n][nT2][1]))
                    eintfile.write("%7.1f%22.13f\n" % (dict_eph_t[T1][n][nT2][0], dict_eph_t[T1][n][nT2][1]))

    for nvol, vol in enumerate(linspace(max(free_v[nT][:, 0]), min(free_v[nT][:, 0]), nvolumes)):
        filec_v_v.write("# V = %f (a.u.^3)\n" % vol)
        filehelmfree_v.write("# V= %17.13f\n" % vol)
        filethermalpressure_v.write("# V= %17.13f\n" % vol)
        for T1 in tdata_read:
            for nT2, T2 in enumerate(tdata_write):
                if T2 == T1:
                    filec_v_v.write("%7.1f%22.13f\n" % (dict_c_v_v[T1][nvol][nT2][0], dict_c_v_v[T1][nvol][nT2][1]))
                    filehelmfree_v.write("%17.13f%22.13f\n" % (dict_helmholtz_v[T1][nvol][nT2][0], dict_helmholtz_v[T1][nvol][nT2][1]))
                    filethermalpressure_v.write("%7.1f%22.13f\n" % (dict_thermalpressure_v[T1][nT2][nvol][0], dict_thermalpressure_v[T1][nT2][nvol][1]))

    for T1 in tdata_read:
        for nT2, T2 in enumerate(tdata_write):
            if T2 == T1:
                filehelmfree_t.write("# T= %f K\n" % T1)
                for n in range(len(ve[0])):
                    filehelmfree_t.write("%9.4f %16.9f\n" % (dict_free_v[T1][nT2][n][0], dict_free_v[T1][nT2][n][1]))

    fileb_s.close()
    filec_p.close()
    filec_v.close()
    filevt_p.close()
    filealpha.close()
    fileb_t_p.close()
    cvdirectfile.close()
    filehelmfree_t.close()
    filehelmfree_v.close()
    filegibbsfree_p.close()
    filethermalpressure_v.close()
