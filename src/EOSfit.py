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
from pylab import polyfit
from ReadInput import indict
from scipy.optimize import leastsq
from Constants import RyperBohr3toGPa


def residuals_murnaghan(pars, y, x):
    return y - murnaghan(pars, x)


def murnaghan(parameters, vol):
    """
    From PRB 28,5480 (1983)
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    e = e0 + b0*vol/bp*(((v0/vol)**bp)/(bp-1)+1) - v0*b0/(bp-1)

    return e


def residuals_birch(pars, y, x):
    return y - birch(pars, x)


def birch(parameters, v):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Princples
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    paper downloaded from Web

    case where n=0
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    e = (e0 + 9.0/8.0*b0*v0*((v0/v)**(2.0/3.0) - 1.0)**2
            + 9.0/16.0*b0*v0*(bp-4.)*((v0/v)**(2.0/3.0) - 1.0)**3)

    return e


def residuals_birch_murnaghan3(pars, y, x):
    return y - birch_murnaghan3(pars, x)


def birch_murnaghan3(parameters, vol):
    """
    BirchMurnaghan equation from PRB 70, 224107
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    eta = (v0/vol)**(1./3.)
    e = e0 + 9.*b0*v0/16.*(eta**2-1)**2*(6 + bp*(eta**2-1.) - 4.*eta**2)

    return e


def residuals_birch_murnaghan4(pars, y, x):
    return y - birch_murnaghan4(pars, x)


def birch_murnaghan4(parameters, vol):
    """
    BirchMurnaghan equation from PRB 70, 224107
    """
    e0 = parameters[0]
    b0 = parameters[1]
    b0p = parameters[2]
    v0 = parameters[3]
    b0pp = parameters[4]

    t1 = (v0/vol)**(1./3.)
    t2 = t1**2
    t3 = t2-1.
    t4 = t3**2/4.
    t5 = b0p**2

    e = e0 + 3./8.*b0*v0*t4*(9.*t4*b0*b0pp+9.*t4*t5-63.*t4*b0p+143.*t4+6.*b0p*t3-24.*t2+36.)

    return e


def residuals_vinet(pars, y, x):
    return y - vinet(pars, x)


def vinet(parameters, vol):
    """
    Vinet equation from PRB 70, 224107
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    eta = (vol/v0)**(1./3.)

    e = (e0 + 2.*b0*v0/(bp-1.)**2 * (2. - (5. + 3. * bp*(eta-1.)-3.*eta)*np.exp(-3.*(bp-1.)*(eta-1.)/2.)))

    return e


def residuals_universal(pars, y, x):
    return y - universal(pars, x)


def universal(parameters, vol):
    """
    Universal equation of state(Vinet P et al., J. Phys.: Condens. Matter 1, p1941 (1989))
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    t1 = b0*v0
    t2 = bp-1.
    t3 = (vol/v0)**(1./3.)
    t4 = np.exp(-3./2.*t2*(-1.+t3))
    t5 = t2**2
    t6 = 1./t5
    e = e0 - 2. * t1 * t4 * (3. * t3 * bp-3. * t3+5.-3. * bp) * t6+4. * t1 * t6

    return e


def residuals_natural_strain3(pars, y, x):
    return y - natural_strain3(pars, x)


def natural_strain3(parameters, vol):
    e0 = parameters[0]
    b0 = parameters[1]
    b0p = parameters[2]
    v0 = parameters[3]

    t1 = b0*v0
    t2 = np.log(v0/vol)
    t3 = t2**2
    t4 = t3*t2
    e = e0 + t1*t3/2.+t1*t4*b0p/6.-t1*t4/3.

    return e


def residuals_natural_strain4(pars, y, x):
    return y - natural_strain4(pars, x)


def natural_strain4(parameters, vol):
    e0 = parameters[0]
    b0 = parameters[1]
    b0p = parameters[2]
    v0 = parameters[3]
    b0pp = parameters[4]

    t1 = b0*v0
    t2 = np.log(v0/vol)
    t3 = t2**2
    t4 = t3**2
    t5 = b0**2
    t6 = b0p**2
    t7 = t3*t2
    e = e0 + t1*t4/8.+t5*v0*t4*b0pp/24.-t1*t4*b0p/8.+t1*t4*t6/24.+t1*t7*b0p/6.-t1*t7/3.+t1*t3/2.

    return e


if indict['Eos_Name'][0][0] == 'p' or indict['Eos_Name'][0][0] == 'P':
    Eos_Name = 1
else:
    Eos_Name = int(indict['Eos_Name'][0])

if Eos_Name == 1:
    eos = murnaghan
    fresiduals = residuals_murnaghan
elif Eos_Name == 2:
    eos = birch
    fresiduals = residuals_birch
elif Eos_Name == 3:
    eos = birch_murnaghan3
    fresiduals = residuals_birch_murnaghan3
elif Eos_Name == 4:
    eos = birch_murnaghan4
    fresiduals = residuals_birch_murnaghan4
elif Eos_Name == 5:
    eos = vinet
    fresiduals = residuals_vinet
elif Eos_Name == 6:
    eos = universal
    fresiduals = residuals_universal
elif Eos_Name == 7:
    eos = natural_strain3
    fresiduals = residuals_natural_strain3
elif Eos_Name == 8:
    eos = natural_strain4
    fresiduals = residuals_natural_strain4


# ======================= End of Equation of state definitions =======================

# ========================== fitting using the selected eos ==========================


class EOSfit:
    def __init__(self, v=[], e=[], residuals=fresiduals):
        self.eos = eos
        self.eos_pars = []
        self.v = np.array(v)
        self.e = np.array(e)
        self.residuals = residuals

    def getparas(self):
        a, b, c = polyfit(self.v, self.e, 2)  # this is from pylab
        # initial guesses.
        v0 = -b/(2*a)
        e0 = a*v0**2 + b*v0 + c
        b0 = 2*a*v0
        bp = 4
        bpp = 1.0
        x0 = [e0, b0, bp, v0, bpp]
        eos_pars, ier = leastsq(self.residuals, x0, args=(self.e, self.v))
        self.eos_pars = [eos_pars[0], eos_pars[1], eos_pars[2], eos_pars[3], eos_pars[4]]

        return eos_pars

    def getenergy(self, vol):
        return self.eos(self.eos_pars, vol)

    def getpressure(self, vol):
        dv = 0.0001
        pressure = -RyperBohr3toGPa*(self.eos(self.eos_pars, vol + dv) - self.eos(self.eos_pars, vol - dv))/dv/2

        return pressure
