"""
  Phasego -- Automatic calculation and plot of phase diagram

  Copyright (C) 2012-2015 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

py_m=["Ana_Phase_Transition","Anharmonicity_extract", "Constants", "CVP_Gamma_BTP_CP",\
      "Debye_temp", "Entropy_CVV_Hugoniot", "EOSfit", "Find_Comp_Ratio", "Find_Cross_Point", "GetB_T_T",\
      "GibbsFree_en", "Helm_free", "Parse_ve", "Plot_Matplotlib", "Polynomial_fit",\
      "Products_Parse", "Products_Gibbs", "PV_T2VT_P", "ReadInput", "setup", "Single_phase", "Thermal_Expan", \
      "Thermal_Pressure"]

setup(
      name="phasego",
      version="3.0.0",
      description="A toolkit for automatic calculation and plot of phase diagram",
      author="Zhong-Li Liu",
      author_email="zl.liu@163.com",
      url="https://sourceforge.net/projects/phasego/",
      license="GNU GPL version 3",
      py_modules=py_m,
      package_dir = {'':'src'},
      scripts=["src/phasego"]
      )
