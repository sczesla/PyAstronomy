from PyAstronomy import pyaC
from PyAstronomy.pyasl import _ic, _moduleImportStatus

_modules = ["pizzolato2003", "expCorrRN", "ramirez2005", "stellarAge", "ballesterosBV_T", \
           "abundances", "spectralClass", "mmscetsTable", "roche", "sysrem", "parallactic_angle", \
           "scaleHeight", "solarDirectFlux"]

for m in _modules:
    _moduleImportStatus[(__name__, m)] = pyaC.pyaimportallfrom(m, "PyAstronomy.pyasl.asl.aslExt_1", globals())

if _ic.check["six"]:
    _moduleImportStatus[(__name__, "transCurves")] = pyaC.pyaimportallfrom("transCurves", "PyAstronomy.pyasl.asl.aslExt_1", globals())