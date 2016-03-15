from __future__ import print_function, division
import importlib
from PyAstronomy import pyaC

# https://stackoverflow.com/questions/21221358/python-how-to-import-all-methods-and-attributes-from-a-module-dynamically
# def pyaslimp(n):
#   try:
#     n = importlib.import_module("." + n, "PyAstronomy.pyasl.asl")
#     try:
#       to_import = n.__all__
#     except AttributeError:
#       to_import = [name for name in n.__dict__ if not name.startswith('_')]
#       for name in to_import:
#         if name in globals():
#           print("Reassigning: ", name, " from ", n)
#       
#     globals().update({name: n.__dict__[name] for name in to_import})
#   except:
#     print(n)


_modules = ["astroTimeLegacy", "airtovac", "aitoffLegacy", "folding", \
           "keplerOrbit", "inTransit", "statTest", "binning", "smooth", \
           "planet_Teq", "planck", "planetPhase", "fluxConversion", "lineOfSight", \
           "lineOfSight", "dopplerShift", "fluxConversion", "baryvel"]

# from astroTimeLegacy import *
# from airtovac import *
# from aitoffLegacy import *
# from folding import *
# from keplerOrbit import *
# from inTransit import *
# from statTest import *
# from binning import *
# from smooth import *
# from planet_Teq import *
# from planck import *
# from planetPhase import *
# from fluxConversion import *
# from lineOfSight import *
# from dopplerShift import *
# from fluxConversion import *
# from baryvel import *

_modules.extend(["idlMod", "idlCirrange", "rotBroad", "svd", "unred", "estimateSNR", \
                "outlier", "intep", "magnitudes", "eq2hor", "sunpos", "observatory"])

# from idlMod import *
# from idlCirrange import *
# from rotBroad import *
# from svd import *
# from unred import *
# from estimateSNR import *
# from outlier import *
# from intep import *
# from magnitudes import *
# from eq2hor import *
# from sunpos import *
# from observatory import *

_modules.extend(["moonpos", "moonphase", "localtime", "lineWidth", "cardinalPoint", \
                "twilight", "angularDistance", "airmass", "decimalYear", "broad", \
                "crosscorr", "quadextreme", "readFitsSpec", "posAngle", "sindex", \
                "atomicNo", "coordinates", "writeFitsSpec", "aslExt_1"])

# from moonpos import *
# from moonphase import *
# from localtime import *
# from lineWidth import *
# from cardinalPoint import *
# from twilight import *
# from angularDistance import *
# from airmass import *
# from decimalYear import *
# from broad import *
# from crosscorr import *
# from quadextreme import *
# from readFitsSpec import *
# from posAngle import *
# from sindex import *
# from atomicNo import *
# from coordinates import *
# from writeFitsSpec import *

from .aslExt_1 import *

for m in _modules:
  pyaC.pyaimport(m, "PyAstronomy.pyasl.asl", globals())

