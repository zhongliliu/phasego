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

import re
from Find_Comp_Ratio import *
from itertools import combinations


def name_parse(Name):
    chemical_symbols = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
                       'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
                       'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                       'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
                       'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                       'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
                       'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                       'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                       'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
                       'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
                       'Fm', 'Md', 'No', 'Lr']

    tuple_list = re.findall(r'([A-Z][a-z]*|[a-z][a-z]*)(\d*)', Name)
    # print tuple_list
    name_dict = {}
    for i in tuple_list:
        if i[0].capitalize() in chemical_symbols:
            if i[1]:
                name_dict[i[0].capitalize()] = int(i[1])
            else:
                name_dict[i[0].capitalize()] = 1
        else:
            print "Something in Formula_Name cannot be parsed. Please check!!!"
            exit(0)

    return name_dict


def products_parse(Formula_Name, Names_Components):
    Formula_Name_dict = name_parse(Formula_Name)
    # print Formula_Name_dict
    # print
    name_comp_list = []
    for Name in Names_Components:
        name_dict = name_parse(Name)
        name_comp_list.append(name_dict)
        # print name_dict

    # print name_comp_list

    Comp_Ratio_Dict = {}
    for n in range(len(name_comp_list) + 1):
        if n > 1:
            comb = list(combinations(range(len(name_comp_list)), n))
        # print n
        if n == 2 and len(name_comp_list) >= 2:
            # print comb
            for cc in comb:
                n1, n2, n3 = find2(name_comp_list[cc[0]], name_comp_list[cc[1]], Formula_Name_dict)
                # print n1, n2, n3
                if n3 != 1000:
                    # print cc
                    Comp_Ratio_Dict["%s:%sx%s" %(n1, n2, n3)] = cc

        if n == 3 and len(name_comp_list) >= 3:
            for cc in comb:
                n1, n2, n3, n4 = find3(name_comp_list[cc[0]], name_comp_list[cc[1]],
                                       name_comp_list[cc[2]], Formula_Name_dict)
                # print n1, n2, n3, n4
                if n4 != 1000:
                    # print cc
                    Comp_Ratio_Dict["%s:%s:%sx%s" % (n1, n2, n3, n4)] = cc

        if n == 4 and len(name_comp_list) >= 4:
            for cc in comb:
                n1, n2, n3, n4, n5 = find4(name_comp_list[cc[0]], name_comp_list[cc[1]], name_comp_list[cc[2]],
                                           name_comp_list[cc[3]], Formula_Name_dict)
                # print n1, n2, n3, n4, n5
                if n5 != 1000:
                    # print cc
                    Comp_Ratio_Dict["%s:%s:%s:%sx%s" % (n1, n2, n3, n4, n5)] = cc

        if n == 5 and len(name_comp_list) >= 5:
            for cc in comb:
                n1, n2, n3, n4, n5, n6 = find5(name_comp_list[cc[0]], name_comp_list[cc[1]], name_comp_list[cc[2]],
                                               name_comp_list[cc[3]], name_comp_list[cc[4]], Formula_Name_dict)
                # print n1, n2, n3, n4, n5, n6
                if n6 != 1000:
                    # print cc
                    Comp_Ratio_Dict["%s:%s:%s:%s:%sx%s" % (n1, n2, n3, n4, n5, n6)] = cc

    # print Comp_Ratio_Dict
    if not Comp_Ratio_Dict:
        print "No proper product ratio, please check input parameters!!!"
        # exit(0)

    return Comp_Ratio_Dict

if __name__ == '__main__':

    Formula_Name = "Al2S3O12"

    # Names_Components = ["Mg3O", "Mg2O3", "MgO", "MgO2", "MgO4"]
    # Names_Components = ["Mg", "MgO"]
    # Comp_Ratio_Dict = products_parse(Formula_Name, Names_Components)
    # print Comp_Ratio_Dict

    Formula_Name_dict = name_parse(Formula_Name)
    print Formula_Name_dict
