# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import ImportCheck

_ic = ImportCheck(["numpy", "quantities", "scipy", "matplotlib", "pyfits", "astropy.io.fits", "pandas"], required=["numpy", "scipy"])

from PyAstronomy.pyasl.asl import *
from PyAstronomy.pyasl.resBased import *