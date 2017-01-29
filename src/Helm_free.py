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

import copy
import numpy as np
from scipy import integrate
from ReadInput import indict
from Parse_ve import parse_ve
from Constants import h, kb, e0, R, ry
from Products_Parse import name_parse


# Read in ve and phonon dos info and calculate Helmholtz free energy
def helm_free(num_formula_units, inppath, outpath,  tdata, T_accurate):
    ve = {}
    free_t = {}
    cvph_v = {}
    cvph_t = {}
    sph_t = {}
    eph_t = {}
    if indict['If_Incl_Electronic_Excitation'][0] == 'yes':
        for nT, T in enumerate(tdata):
            # Parse static volume-energy data in unit: Bohr^3: Ryd.
            ve[T] = parse_ve(inppath + "/" + indict['VE_data_File_Name'][0] + str(T))
    else:
        ve[0] = parse_ve(inppath + "/" + indict['VE_data_File_Name'][0] + "0")

    # if include thermal effect, then prepare output files
    if indict['If_Incl_Phonon'][0] == 'yes':
        cvdirectfile = open(outpath + "/C_V-ph-direct-%s.dat" % (str(T_accurate)), "w")
        cvdirectfile.write("# T (K)      C-ph_V_V (J K^-1 mol^-1)" + '\n')
        sdirectfile = open(outpath + "/Entropy-ph_V-%s.dat" % (str(T_accurate)), "w")
        sdirectfile.write("# T (K)      Entropy-ph (J K^-1 mol^-1)" + '\n')
        eintfile = open(outpath + "/E-ph_V-%s.dat" % (str(T_accurate)), "w")
        eintfile.write("# T (K)       E-ph (Ry.)" + '\n')

        if indict['Read_what'][0] == 'dos_data':
            filehelmfree_v = open(outpath + "/tf.dat", "w")
            filehelmfree_v.write("# T (K)" + 15 * " " + "Helmholtz Free Energy (Ry.)\n")

        if indict['Read_what'][0] == 'TF_data':
            fhelmfree_v = open(inppath + "/" + indict['TF_data_File_Name'][0], "r")
            fhelmfree_v.readline()

    # Read in phonon dos and calculate ve plus phonon contribution
    if indict['If_Incl_Phonon'][0] == 'yes':
        free_t = np.zeros(shape=(len(ve[0]), len(tdata), 2))
        free_v = np.zeros(shape=(len(tdata), len(ve[0]), 2))
        eph_t = np.zeros(shape=(len(ve[0]), len(tdata), 2))
        eph_v = np.zeros(shape=(len(tdata), len(ve[0]), 2))
        cvph_t = np.zeros(shape=(len(ve[0]), len(tdata), 2))
        cvph_v = np.zeros(shape=(len(tdata), len(ve[0]), 2))
        sph_t = np.zeros(shape=(len(ve[0]), len(tdata), 2))

        for n in range(len(ve[0])):
            if indict['Read_what'][0] == 'TF_data':
                fhelmfree_v.readline()

            if indict['If_Incl_Phonon'][0] == 'yes' and indict['Read_what'][0] == 'dos_data':
                cvdirectfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / num_formula_units))
                sdirectfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / num_formula_units))
                eintfile.write("# V = %f (a.u.^3)\n" % (float(ve[0][n][0]) / num_formula_units))

            if indict['Units_VE'][0] == 'Bohr3':
                v_conv = 1.0
            elif indict['Units_VE'][0] == 'A3':
                v_conv = 1 / 0.529177 ** 3.0

            if indict['Read_what'][0] == 'dos_data':
                filehelmfree_v.write("# V= %17.13f\n" % ve[0][n][0])
                phdos = np.loadtxt(inppath + "/%s%.4f-%s" % (indict['Ph_Dos_File_Base_Name'][0], ve[0][n][0]/v_conv,
                                                             str(T_accurate)))

                if indict['Unit_of_Freq'][0] == 'cm-1':
                    freq = copy.deepcopy(phdos[:, 0])
                    dos = copy.deepcopy(phdos[:, 1])

                if indict['Unit_of_Freq'][0] == 'THz':
                    freq = copy.deepcopy(phdos[:, 0]) * 33.3564
                    dos = copy.deepcopy(phdos[:, 1]) / 33.3564

                if indict['Unit_of_Freq'][0] == 'meV':
                    freq = copy.deepcopy(phdos[:, 0]) * 8.065541
                    dos = copy.deepcopy(phdos[:, 1]) / 8.065541

                if indict['Unit_of_Freq'][0] == 'eV':
                    freq = copy.deepcopy(phdos[:, 0]) * 8065.541
                    dos = copy.deepcopy(phdos[:, 1]) / 8065.541

#                ds = list(dos)
#                fr = list(freq)
#                del ds[:list(dos).index(max(dos))]
#                del fr[:list(dos).index(max(dos))]
                integral = integrate.simps(list(dos), list(freq))

                if int(round(integral)) == 1:
                    # change to 3*num_formula_units*Num_Atoms for PHON code in which dos is normalized to 1.0.
                    Formula_Name_dict = name_parse(indict['Formula_Name'][0])
                    num_atoms = sum(Formula_Name_dict.values())
                    # print num_atoms
                    dos *= 3 * num_formula_units * num_atoms

                for nT, T in enumerate(tdata):

                    if T == 0:
                        T = 10

                    for i, f in enumerate(freq):
                        if f < 1.5 * 2.2e-16:
                            freq[i] = 0.001
                            dos[i] = 0

                    y0 = 0.5 * h * freq * dos + kb * T * np.log(1.0 - np.exp(-h * freq / (kb * T))) * dos
                    # y0 = 0.5 * h * freq * dos + kb * T * np.log(2.0 * np.sinh(h * freq / (2.0 * kb * T))) * dos
                    # free energy:
                    if indict['Names_of_Strs'][0] == 'liquid':
                        y1 = 1.055 * y0
                        y2 = kb * T * np.log(h*freq/(kb*T))*dos
                    else:
                        y1 = y0

                    f_t = integrate.trapz(y1, freq) / e0 / ry
                    # print f_t
                    # f_t += -0.5 * alpha * T * T*6.333624*10**(-06)

                    if indict['Names_of_Strs'][0] == 'liquid':
                        f_t2 = integrate.trapz(y2, freq) / e0 / ry
                        print f_t2,f_t
                    #else:
                        # f_t3 = integrate.trapz(y3, freq) / e0 / ry
                        # f_t3 = y3
                        # print T,f_t, f_t3
                    #    pass
                        
                    # print T,F_T
                    if T == 10:
                        T = 0

                    if indict['If_Incl_Electronic_Excitation'][0] == 'yes':
                        free_t[n][nT] = [T, ve[T][n][1] / num_formula_units + f_t / num_formula_units]
                    else:
                        free_t[n][nT] = [T, ve[0][n][1] / num_formula_units + f_t / num_formula_units]

                    filehelmfree_v.write("%8.2f%22.13f\n" % (free_t[n][nT][0], free_t[n][nT][1]))
                    # direct calc. of Cv, S and E.
                    # from C. Lee, X. Gonze, Phys. Rev. B 51, 8610 (1995)
                    if indict['If_Incl_Phonon'][0] == 'yes':
                        if T == 0:
                            T = 5

                        x = h * freq / (2.0 * kb * T)
                        # C_V:
                        #  csch(x) = 1 / sinh(x)
                        y2 = dos * x ** 2.0 / np.sinh(x) ** 2.0
                        cv = integrate.trapz(y2, freq) * R
                        cvph_t[n][nT] = [T, cv / num_formula_units]

                        cvdirectfile.write("%7.1f%22.13f\n" % (T, cv / num_formula_units))

                        # Entropy:
                        # coth(x) = 1 / tanh(x)
                        y3 = dos * (x / np.tanh(x) - np.log(2.0 * np.sinh(x)))
                        s = integrate.trapz(y3, freq) * R
                        sph_t[n][nT] = [T, s  / num_formula_units]

                        sdirectfile.write("%7.1f%22.13f\n" % (T, s / num_formula_units))
                        # print S

                        # Internal energy:
                        # coth(x) = 1 / tanh(x)
                        y4 = 0.5 * h * freq * dos / np.tanh(x)
                        e = integrate.trapz(y4, freq) / e0 / ry  # Ry.
                        eph_t[n][nT] = [T, e / num_formula_units]

                        eintfile.write("%7.1f%22.13f\n" % (T, e / num_formula_units))
                        # print E

            if indict['Read_what'][0] == 'TF_data':
                for nT, T in enumerate(tdata):
                    line = fhelmfree_v.readline().split()
                    free_t[n][nT][0] = float(line[0])
                    free_t[n][nT][1] = float(line[1])

                    # print free_t[n][nT][0],free_t[n][nT][1]

        for nT, T in enumerate(tdata):
            for n in range(len(ve[0])):
                if indict['If_Incl_Electronic_Excitation'][0] == 'yes':
                    free_v[nT][n] = np.array([ve[T][n][0] / num_formula_units, free_t[n][nT][1]])
                    eph_v[nT][n] = np.array([ve[T][n][0] / num_formula_units, eph_t[n][nT][1]])
                    cvph_v[nT][n] = np.array([ve[T][n][0] / num_formula_units, cvph_t[n][nT][1]])
                else:
                    free_v[nT][n] = np.array([ve[0][n][0] / num_formula_units, free_t[n][nT][1]])
                    eph_v[nT][n] = np.array([ve[0][n][0] / num_formula_units, eph_t[n][nT][1]])
                    cvph_v[nT][n] = np.array([ve[0][n][0] / num_formula_units, cvph_t[n][nT][1]])

    else:
        free_v = np.zeros(shape=(1, len(ve[0]), 2))
        for n in range(len(ve[0])):
            free_v[0][n] = np.array([ve[0][n][0], ve[0][n][1]])

        free_v /= num_formula_units

    if indict['If_Incl_Phonon'][0] == 'yes':
        cvdirectfile.close()
        sdirectfile.close()
        eintfile.close()

    return ve, free_v, free_t, cvph_v, cvph_t, sph_t, eph_t


if __name__ == '__main__':
    Tdata = np.arange(int(indict['Tdata'][0]), int(indict['Tdata'][1]) + int(indict['Tdata'][2]),
                      int(indict['Tdata'][2]))
    num_formula_units = int(indict['Num_Formula_Units'][0])
    inppath = 'inp-fcc'
    outdir = './'
    free_V, free_T = helm_free(num_formula_units, inppath, outdir, Tdata)
    print free_V, free_T
