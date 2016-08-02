# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import ImportCheck
from PyAstronomy.pyaC import pyaErrors as PE

ic = ImportCheck(["numpy", "scipy"], required=["numpy"]) 

from .periodBase import *
from .lombScargle import *
from .gls import *

if ic.check["scipy"]:
  from .fourier import *
else:
  PE.warn(PE.PyARequiredImport("No Fourier transform will be available, because scipy cannot be imported.", \
                               where="Import of PyPeriod", solution="Install scipy."))