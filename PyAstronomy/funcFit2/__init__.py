# -*- coding: utf-8 -*-
from __future__ import print_function
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyaC import ImportCheck
import sys

# Check Python version
_majV, _minV = sys.version_info[0:2]
if (_minV < 7) and (_majV == 2):
    PE.warn("funcFit needs 2.7.x or greater. See documentation (Prerequisites) for explanation.")

ic = ImportCheck(["numpy", "scipy", "matplotlib",
                  "matplotlib.pylab", "pyfits", "emcee", "progressbar"])

# Get out if numpy not present
if not ic.check["numpy"]:
    raise(PE.PyARequiredImport("Numpy cannot be imported.", solution="Install numpy (see http://numpy.scipy.org/, you probably should also install SciPy).",
                               addInfo="The numpy package provides array support for Python and is indispensable in many scientific applications."))

# Check whether fitting modules can be imported
_scoImport = ic.check["scipy"]
_emceeImport = ic.check["emcee"]
_mplImport = ic.check["matplotlib"]

from PyAstronomy.funcFit.utils import *
from .fbo2 import *
from .anaMCMCTraces2 import TraceAnalysis2

def __question(tf):
    if tf:
        return "yes"
    return "no"


def status():
    print("Status of funcFit2:")
    print("--------------------------")
    print("Is scipy.optimize available? ", __question(_scoImport))
    if _scoImport:
        print("  Version: ", ic.versions["scipy"])
    if not _scoImport:
        print("  The minimization algorithms from SciPy cannot be used.")
        print("  Install SciPy to make them available (use 'pip install scipy' or see http://www.scipy.org/).")
    print("Is emcee available? ", __question(_emceeImport))
    if _emceeImport:
        print("  Version: ", ic.versions["emcee"])
    if not _emceeImport:
        print("  No Markov-Chain Monte-Carlo sampling via emcee will be possible.")
        print("  Install emcee to enable it (use 'pip install emcee' or see http://dfm.io/emcee/current/).")

