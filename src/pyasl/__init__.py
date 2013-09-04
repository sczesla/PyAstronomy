# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import ImportCheck

_ic = ImportCheck(["numpy", "quantities", "scipy", "matplotlib", "pyfits"], required=["numpy"])

from PyAstronomy.pyasl.asl import *
from PyAstronomy.pyasl.resBased import *