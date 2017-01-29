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

from ReadInput import indict


def products_gibbs(Comp_Ratio_Dict, Comp_Ratios, Names_Components, BigDict_Comp_P, BigDict_Comp_T, BigDict_P, BigDict_T):
    for ii in Comp_Ratios:
        # print ii
        n = ii.count(":") + 1
        if n == 2:
            for it in BigDict_Comp_P.keys():
                for jt in BigDict_Comp_P.keys():
                    if it.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][0]] and \
                                    jt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][1]]:
                        BigDict_P[it + " + " + jt] = {}
                        BigDict_T[it + " + " + jt] = {}

                        ni, nj = tuple(ii.split("x")[0].split(":"))
                        multi = int(ii.split("x")[1])
                        print str(multi)+indict['Formula_Name'][0], "<-->", ni+it, "+", nj+jt
                        for P in BigDict_Comp_P[it].keys():
                            new_p = (int(ni) * BigDict_Comp_P[it][P] +
                                     int(nj) * BigDict_Comp_P[jt][P]) / multi
                            new_p[:, 0] *= float(multi) / float((int(ni) + int(nj)))
                            BigDict_P[it + " + " + jt][P] = new_p
                            # print new_p

                        for T in BigDict_Comp_T[it].keys():
                            new_t = (int(ni) * BigDict_Comp_T[it][T] +
                                     int(nj) * BigDict_Comp_T[jt][T]) / multi
                            new_t[:, 0] *= float(multi) / float((int(ni) + int(nj)))
                            BigDict_T[it + " + " + jt][T] = new_t
                            # print new_t

        elif n == 3:
            for it in BigDict_Comp_P.keys():
                for jt in BigDict_Comp_P.keys():
                    for kt in BigDict_Comp_P.keys():
                        if it.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][0]] and \
                                        jt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][1]] and \
                                        kt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][2]]:
                            BigDict_P[it + " + " + jt + " + " + kt] = {}
                            BigDict_T[it + " + " + jt + " + " + kt] = {}

                            ni, nj, nk = tuple(ii.split("x")[0].split(":"))
                            multi = int(ii.split("x")[1])
                            print str(multi) + indict['Formula_Name'][0], "<-->", ni+it, "+", nj+jt, "+", nk+kt
                            # print nj,nk
                            for P in BigDict_Comp_P[it].keys():
                                new_p = (int(ni) * BigDict_Comp_P[it][P] +
                                        int(nj) * BigDict_Comp_P[jt][P] +
                                        int(nk) * BigDict_Comp_P[kt][P]) / multi
                                new_p[:, 0] *= float(multi) / float((int(ni) + int(nj) + int(nk)))
                                BigDict_P[it + " + " + jt + " + " + kt][P] = new_p
                                # print new_p

                            for T in BigDict_Comp_T[it].keys():
                                new_t = (int(ni) * BigDict_Comp_T[it][T] +
                                        int(nj) * BigDict_Comp_T[jt][T] +
                                        int(nk) * BigDict_Comp_T[kt][T]) / multi
                                new_t[:, 0] *= float(multi) / float((int(ni) + int(nj) + int(nk)))
                                BigDict_T[it + " + " + jt + " + " + kt][T] = new_t
                                # print new_t

        elif n == 4:
            for it in BigDict_Comp_P.keys():
                for jt in BigDict_Comp_P.keys():
                    for kt in BigDict_Comp_P.keys():
                        for lt in BigDict_Comp_P.keys():
                            if it.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][0]] and \
                                            jt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][1]] and \
                                            kt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][2]] and \
                                            lt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][3]]:
                                BigDict_P[it + " + " + jt + " + " + kt + " + " + lt] = {}
                                BigDict_T[it + " + " + jt + " + " + kt + " + " + lt] = {}

                                ni, nj, nk, nl = tuple(ii.split("x")[0].split(":"))
                                multi = int(ii.split("x")[1])
                                print str(multi)+indict['Formula_Name'][0], "<-->", \
                                    ni+it, "+", nj+jt, "+", nk+kt, "+", nl+lt
                                # print nj,nk
                                for P in BigDict_Comp_P[it].keys():
                                    new_p = (int(ni) * BigDict_Comp_P[it][P] +
                                             int(nj) * BigDict_Comp_P[jt][P] +
                                             int(nk) * BigDict_Comp_P[kt][P] +
                                             int(nl) * BigDict_Comp_P[lt][P]) / multi
                                    new_p[:, 0] *= float(multi) / float((int(ni) + int(nj) + int(nk) + int(nl)))
                                    BigDict_P[it + " + " + jt + " + " + kt + " + " + lt][P] = new_p
                                    # print new_p

                                for T in BigDict_Comp_T[it].keys():
                                    new_t = (int(ni) * BigDict_Comp_T[it][T] +
                                             int(nj) * BigDict_Comp_T[jt][T] +
                                             int(nk) * BigDict_Comp_T[kt][T] +
                                             int(nl) * BigDict_Comp_T[lt][T]) / multi
                                    new_t[:, 0] *= float(multi) / float((int(ni) + int(nj) + int(nk) + int(nl)))
                                    BigDict_T[it + " + " + jt + " + " + kt + " + " + lt][T] = new_t
                                    # print new_t

        elif n == 5:
            for it in BigDict_Comp_P.keys():
                for jt in BigDict_Comp_P.keys():
                    for kt in BigDict_Comp_P.keys():
                        for lt in BigDict_Comp_P.keys():
                            for mt in BigDict_Comp_P.keys():
                                if it.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][0]] and \
                                                jt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][1]] and \
                                                kt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][2]] and \
                                                lt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][3]] and \
                                                mt.split("-")[1] == Names_Components[Comp_Ratio_Dict[ii][4]]:
                                    BigDict_P[it + " + " + jt + " + " + kt + " + " + lt + " + " + mt] = {}
                                    BigDict_T[it + " + " + jt + " + " + kt + " + " + lt + " + " + mt] = {}

                                    ni, nj, nk, nl, nm = tuple(ii.split("x")[0].split(":"))
                                    multi = int(ii.split("x")[1])
                                    print str(multi)+indict['Formula_Name'][0], "<-->", \
                                        ni+it, "+", nj+jt, "+", nk+kt, "+", nl+lt, "+", nm+mt
                                    for P in BigDict_Comp_P[it].keys():
                                        new_p = (int(ni) * BigDict_Comp_P[it][P] +
                                                 int(nj) * BigDict_Comp_P[jt][P] +
                                                 int(nk) * BigDict_Comp_P[kt][P] +
                                                 int(nl) * BigDict_Comp_P[lt][P] +
                                                 int(nm) * BigDict_Comp_P[mt][P]) / multi
                                        new_p[:, 0] *= float(multi) / float((int(ni) + \
                                                        int(nj) + int(nk) + int(nl) + int(nm)))
                                        BigDict_P[it + " + " + jt + " + " + kt + " + " + lt + " + " + mt][P] = new_p
                                        # print new_p

                                    for T in BigDict_Comp_T[it].keys():
                                        new_t = (int(ni) * BigDict_Comp_T[it][T] +
                                                 int(nj) * BigDict_Comp_T[jt][T] +
                                                 int(nk) * BigDict_Comp_T[kt][T] +
                                                 int(nl) * BigDict_Comp_T[lt][T] +
                                                 int(nm) * BigDict_Comp_T[mt][T]) / multi
                                        new_t[:, 0] *= float(multi) / float((int(ni) + \
                                                        int(nj) + int(nk) + int(nl) + int(nm)))
                                        BigDict_T[it + " + " + jt + " + " + kt + " + " + lt + " + " + mt][T] = new_t
                                        # print new_t

    return BigDict_P
