from __future__ import print_function, division
from PyAstronomy import pyaC
from PyAstronomy.pyasl import _moduleImportStatus

_modules = ["astroTimeLegacy", "airtovac", "aitoffLegacy", "folding", \
           "keplerOrbit", "inTransit", "statTest", "binning", "smooth", \
           "planet_Teq", "planck", "planetPhase", "fluxConversion", "lineOfSight", \
           "lineOfSight", "dopplerShift", "fluxConversion", "baryvel"]


_modules.extend(["idlMod", "idlCirrange", "rotBroad", "svd", "unred", "estimateSNR", \
                "outlier", "intep", "magnitudes", "eq2hor", "sunpos", "observatory"])


_modules.extend(["moonpos", "moonphase", "localtime", "lineWidth", "cardinalPoint", \
                "twilight", "angularDistance", "airmass", "decimalYear", "broad", \
                "crosscorr", "quadextreme", "readFitsSpec", "posAngle", "sindex", \
                "atomicNo", "coordinates", "writeFitsSpec", "aslExt_1"])


from .aslExt_1 import *

for m in _modules:
  _moduleImportStatus[(__name__, m)] = pyaC.pyaimportallfrom(m, "PyAstronomy.pyasl.asl", globals())
 
