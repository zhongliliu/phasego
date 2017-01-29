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


def plus_dict(dict1, dict2):
    for it in dict1:
        if it in dict2:
            dict2[it] += dict1[it]
        if it not in dict2:
            dict2[it] = dict1[it]

    return dict2


def multi_values(dictin, n):
    for i in dictin:
        dictin[i] *= n
    return dictin


def divid_values(dictin, n):
    for i in dictin:
        dictin[i] /= float(n)
    return dictin


def find2(sp1, sp2, Name_sp):
    ni = 1
    n3 = 1000
    if_find = False
    while ni <= 10:
        # print "ni=", ni
        nj = 1
        n1 = ni
        while nj <= 10:
            # print "nj=", nj
            n2 = nj
            # print multi_values(copy.deepcopy(sp1), ni)
            new_dict = plus_dict(multi_values(copy.deepcopy(sp1), ni), multi_values(copy.deepcopy(sp2), nj))
            # print "new_dict=", new_dict

            if new_dict == Name_sp:
                n3 = 1
                if_find = True
                # print ni, nj
                break
            else:
                nk = 1
                while nk <= 10:
                    # print new_dict
                    if divid_values(copy.deepcopy(new_dict), nk) == Name_sp:
                        if_find = True
                        n3 = nk
                        # print nk
                        break
                    nk += 1

            if if_find:
                break
            nj += 1

        if if_find:
            break
        ni += 1

    return n1, n2, n3


def find3(sp1, sp2, sp3, Name_sp):
    ni = 1
    n4 = 1000
    if_find = False
    while ni <= 10:
        # print "ni=", ni
        nj = 1
        n1 = ni
        while nj <= 10:
            # print "nj=", nj
            n2 = nj
            # print multi_values(copy.deepcopy(sp1), ni)
            nm = 1
            while nm <= 10:
                n3 = nm
                dict1 = plus_dict(multi_values(copy.deepcopy(sp1), ni), multi_values(copy.deepcopy(sp2), nj))
                new_dict = plus_dict(copy.deepcopy(dict1), multi_values(copy.deepcopy(sp3), nm))
                # print "new_dict=", new_dict

                if new_dict == Name_sp:
                    n4 = 1
                    if_find = True
                    # print ni, nj
                    break
                else:
                    nk = 1
                    while nk <= 10:
                        # print new_dict
                        if divid_values(copy.deepcopy(new_dict), nk) == Name_sp:
                            if_find = True
                            n4 = nk
                            # print nk
                            break
                        nk += 1

                if if_find:
                    break
                nm += 1

            if if_find:
                break
            nj += 1

        if if_find:
            break
        ni += 1

    return n1, n2, n3, n4


def find4(sp1, sp2, sp3, sp4, Name_sp):
    ni = 1
    n5 = 1000
    if_find = False
    while ni <= 10:
        # print "ni=", ni
        nj = 1
        n1 = ni
        while nj <= 10:
            # print "nj=", nj
            n2 = nj
            # print multi_values(copy.deepcopy(sp1), ni)
            nm = 1
            while nm <= 10:
                n3 = nm
                nn = 1
                while nn <= 10:
                    n4 = nn

                    dict1 = plus_dict(multi_values(copy.deepcopy(sp1), ni), multi_values(copy.deepcopy(sp2), nj))
                    dict2 = plus_dict(copy.deepcopy(dict1), multi_values(copy.deepcopy(sp3), nm))
                    new_dict = plus_dict(copy.deepcopy(dict2), multi_values(copy.deepcopy(sp4), nn))
                    # print ni, nj, nm, nn
                    # print "new_dict=", new_dict

                    if new_dict == Name_sp:
                        n5 = 1
                        if_find = True
                        # print ni, nj
                        break
                    else:
                        nk = 1
                        while nk <= 10:
                            # print new_dict
                            if divid_values(copy.deepcopy(new_dict), nk) == Name_sp:
                                if_find = True
                                n5 = nk
                                # print nk
                                break
                            nk += 1

                    if if_find:
                        break
                    nn += 1

                if if_find:
                    break
                nm += 1

            if if_find:
                break
            nj += 1

        if if_find:
            break
        ni += 1

    return n1, n2, n3, n4, n5


def find5(sp1, sp2, sp3, sp4, sp5, Name_sp):
    ni = 1
    n6 = 1000
    n_max = 10
    if_find = False
    while ni <= n_max:
        # print "ni=", ni
        nj = 1
        n1 = ni
        while nj <= n_max:
            # print "nj=", nj
            n2 = nj
            # print multi_values(copy.deepcopy(sp1), ni)
            nm = 1
            while nm <= n_max:
                n3 = nm
                nn = 1
                while nn <= n_max:
                    n4 = nn
                    nl = 1
                    while nl <= n_max:
                        n5 = nl
                        dict1 = plus_dict(multi_values(copy.deepcopy(sp1), ni), multi_values(copy.deepcopy(sp2), nj))
                        dict2 = plus_dict(copy.deepcopy(dict1), multi_values(copy.deepcopy(sp3), nm))
                        dict3 = plus_dict(copy.deepcopy(dict2), multi_values(copy.deepcopy(sp4), nn))
                        new_dict = plus_dict(copy.deepcopy(dict3), multi_values(copy.deepcopy(sp5), nl))
                        # print ni, nj, nm, nn, nl
                        # print "new_dict=", new_dict

                        if new_dict == Name_sp:
                            n6 = 1
                            if_find = True
                            break
                        else:
                            nk = 1
                            while nk <= n_max:
                                # print new_dict
                                if divid_values(copy.deepcopy(new_dict), nk) == Name_sp:
                                    if_find = True
                                    n6 = nk
                                    # print nk
                                    break
                                nk += 1

                        if if_find:
                            break
                        nl += 1

                    if if_find:
                        break
                    nn += 1

                if if_find:
                    break
                nm += 1

            if if_find:
                break
            nj += 1

        if if_find:
            break
        ni += 1

    return n1, n2, n3, n4, n5, n6


if __name__ == '__main__':
    Name_sp = {"Mg": 3, "O": 2}

    sp1 = {"Mg": 3, "O": 1}
    sp2 = {"Mg": 2, "O": 3}
    sp3 = {"Mg": 1, "O": 1}
    sp4 = {"Mg": 1, "O": 2}
    sp5 = {"Mg": 1, "O": 4}

    n1, n2, n3, n4, n5 = find4(sp1, sp2, sp3, sp4, Name_sp)

    print Name_sp
    print
    print sp1
    print sp2
    print sp3
    print sp4
    # print sp5
    print n1, n2, n3, n4, n5