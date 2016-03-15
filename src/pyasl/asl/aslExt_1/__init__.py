from PyAstronomy import pyaC

from PyAstronomy.pyasl import _ic

_modules = ["pizzolato2003", "expCorrRN", "ramirez2005", "stellarAge", "ballesterosBV_T", \
           "abundances"]

# from pizzolato2003 import *
# from expCorrRN import *
# from ramirez2005 import *
# from stellarAge import *
# from ballesterosBV_T import *
# from abundances import *

for m in _modules:
  pyaC.pyaimport(m, "PyAstronomy.pyasl.asl.aslExt_1", globals())

if _ic.check["six"]:
#   from transCurves import TransmissionCurves
  pyaC.pyaimport("transCurves", "PyAstronomy.pyasl.asl.aslExt_1", globals())