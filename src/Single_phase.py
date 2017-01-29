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
from GetB_T_T import getb_t_t
from PV_T2VT_P import pv_t2vt_p
from Helm_free import helm_free
from Debye_temp import debye_temp
from GibbsFree_en import gibbsfree_en
from Thermal_Expan import thermal_expan
from Thermal_Pressure import thermal_pressure
from CVP_Gamma_BTP_CP import cvp_gamma_btp_cp
from Anharmonicity_extract import anharm_extract1, anharm_extract2
from Entropy_CVV_Hugoniot import entropy_cvv_hugoniot


# Do single phase calculations
def single_phase(num_formula_units, inppath, outpath):
    tdata_write = np.arange(int(indict['Tdata_Write'][0]), int(indict['Tdata_Write'][1]) + int(indict['Tdata_Write'][2]),
                      int(indict['Tdata_Write'][2]))

    pdata = np.arange(int(indict['Pdata'][0]), int(indict['Pdata'][1]) + int(indict['Pdata'][2]),
                      int(indict['Pdata'][2]))

    if indict['If_Incl_Anharm_Phonon'][0] == 'yes':
        tdata_read = np.arange(int(indict['Tdata_Read'][0]),
                               int(indict['Tdata_Read'][1]) + int(indict['Tdata_Read'][2]),
                               int(indict['Tdata_Read'][2]))

        T_accurate = tdata_read[0]

        gdict_p = {}
        gdict_t = {}
        dict_cvph_t = {}
        dict_sph_t = {}
        dict_eph_t = {}
        dict_c_v_v = {}
        dict_free_v = {}
        dict_helmholtz_v = {}
        dict_thermalpressure_v = {}
        dict_vt_p = {}

        alpha_anh = np.zeros(shape=(len(pdata), len(tdata_read), 2))

        filealpha = open(outpath + "/Alpha_anh.dat", "w")
        filealpha.write("#T (K)    alpha (10^(-5) K^(-1))\n")
        fileb_t_t = open(outpath+"/B_T_T_anh.dat", "w")
        filepvt_anh = open(outpath + "/PV_T_anh.dat", "w")
        filepvt_anh.write("# P (GPa)" + 13 * " " + "V (a.u.^3)\n")
        filegibbsfree = open(outpath + "/G_T_anh.dat", "w")
        filegibbsfree.write("# P (GPa)" + 13 * " " + "Gibbs free energy (Ry.)\n")

        print
        print "Anharmonic effects analysis..."
        print
        print "Extracting anharmonic data"
        print "from the phonon DOS at:"
        print

        for T_accurate in tdata_read:

            print "T = %s K ..." %(T_accurate)

            outpath_T = "out-" + indict['Names_of_Strs'][0] + "-" + indict['Formula_Name'][0] + "/" + str(T_accurate)
            if not os.path.isdir("%s" % outpath_T):
                os.mkdir("%s" % outpath_T)

            # Calculate Helmholtz free energy
            ve, free_v, free_t, cvph_v, cvph_t, sph_t, eph_t = helm_free(num_formula_units, inppath, outpath_T, tdata_write, T_accurate)
            dict_cvph_t[T_accurate] = cvph_t
            dict_sph_t[T_accurate] = sph_t
            dict_eph_t[T_accurate] = eph_t
            dict_free_v[T_accurate] = free_v

            # Calculate Gibbs free energy
            pv_t, gibbs_p, helmholtz_v, gibbsdict_p, gibbsdict_t = gibbsfree_en(tdata_write, pdata, free_v, outpath_T, T_accurate)
            dict_helmholtz_v[T_accurate] = helmholtz_v

            # Calculate thermal pressure
            thermalpressure_v = thermal_pressure(tdata_write, free_v, pv_t, outpath_T, T_accurate)
            dict_thermalpressure_v[T_accurate] = thermalpressure_v

            # Calculate bulk moduli
            b_t_t = getb_t_t(tdata_write, pv_t, outpath_T, T_accurate)

            # Change PV_T to VT_P
            vt_p = pv_t2vt_p(tdata_write, pdata, free_v, pv_t, outpath_T, T_accurate)
            dict_vt_p[T_accurate] = vt_p
            if indict['If_Incl_Phonon'][0] == 'yes' and indict['Tdata_Write'][0] != indict['Tdata_Write'][1]:
                alpha = thermal_expan(tdata_write, pdata, vt_p, outpath_T, T_accurate)
                c_v_v = entropy_cvv_hugoniot(tdata_write, free_v, pv_t, helmholtz_v, outpath_T, T_accurate)
                dict_c_v_v[T_accurate] = c_v_v
                b_t_p, b_s_p, c_p_p, c_v_p, gamma = cvp_gamma_btp_cp(tdata_write, pdata, free_v, gibbs_p, pv_t,
                                                                     vt_p, c_v_v, b_t_t, outpath_T, T_accurate)
                if indict['If_Calc_Debye_Temp'][0] == 'yes':
                    debye_temp(num_formula_units, tdata_write, cvph_t, cvph_v, inppath, outpath_T, T_accurate)

            if indict['Tdata_Write'][0] != indict['Tdata_Write'][1]:

                alpha_anh, gdict_p, gdict_t = anharm_extract1(alpha_anh, gdict_p, gdict_t, alpha, free_v, pv_t, b_t_t,
                                                              gibbsdict_t, pdata, tdata_read, tdata_write, gibbsdict_p,
                                                              T_accurate, fileb_t_t, filepvt_anh, filegibbsfree)
        if indict['Tdata_Write'][0] != indict['Tdata_Write'][1]:
            anharm_extract2(alpha_anh, pdata, tdata_read, tdata_write, dict_cvph_t, dict_sph_t, dict_eph_t, dict_c_v_v,
                            dict_free_v, dict_helmholtz_v, dict_thermalpressure_v, dict_vt_p, gibbsdict_p, b_t_p, b_s_p,
                            c_p_p, c_v_p, gamma, ve, num_formula_units, free_v, outpath, filealpha)

        fileb_t_t.close()
        filepvt_anh.close()
        filegibbsfree.close()

        print
        print "=========\/========="
        print

    else:
        T_accurate = 0 # tdata_write[0]
        tdata = tdata_write
        # Calculate Helmholtz free energy
        ve, free_v, free_t, cvph_v, cvph_t, sph_t, eph_t = helm_free(num_formula_units, inppath, outpath, tdata, T_accurate)
        print "-----\\" + 8 * "-" + "/-----"

        # Calculate Gibbs free energy
        pv_t, gibbs_p, helmholtz_v, gdict_p, gdict_t = gibbsfree_en(tdata, pdata, free_v, outpath, T_accurate)
        print "------\\" + 6 * "-" + "/------"

        # Calculate thermal pressure
        thermal_pressure(tdata, free_v, pv_t, outpath, T_accurate)
        print "-------\\" + 4 * "-" + "/-------"

        # Calculate bulk moduli
        b_t_t = getb_t_t(tdata, pv_t, outpath, T_accurate)
        print "--------\\" + 2 * "-" + "/--------"

        # Change PV_T to VT_P
        vt_p = pv_t2vt_p(tdata, pdata, free_v, pv_t, outpath, T_accurate)
        print "---------\/---------"

        if indict['If_Incl_Phonon'][0] == 'yes':
            thermal_expan(tdata, pdata, vt_p, outpath, T_accurate)
            print "------\\" + 6 * "-" + "/------"
            c_v_v = entropy_cvv_hugoniot(tdata, free_v, pv_t, helmholtz_v, outpath, T_accurate)
            print "-------\\" + 4 * "-" + "/-------"
            cvp_gamma_btp_cp(tdata, pdata, free_v, gibbs_p, pv_t, vt_p, c_v_v, b_t_t, outpath, T_accurate)
            print "--------\\" + 2 * "-" + "/--------"
            if indict['If_Calc_Debye_Temp'][0] == 'yes':
                debye_temp(num_formula_units, tdata, cvph_t, cvph_v, inppath, outpath, T_accurate)

        print "=========\/========="
        print

    return gdict_p, gdict_t
