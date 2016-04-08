# -*- coding: utf-8 -*-
from __future__ import print_function
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyaC import ImportCheck
import sys

# Check Python version
_majV, _minV = sys.version_info[0:2]
if (_minV < 7) and (_majV == 2):
  PE.warn("funcFit needs 2.7.x or greater. See documentation (Prerequisites) for explanation.")
  
ic = ImportCheck(["numpy", "scipy", "pymc", "matplotlib", "matplotlib.pylab", "pyfits", "emcee", "progressbar"])

# Get out if numpy not present
if not ic.check["numpy"]:
  raise(PE.PyARequiredImport("Numpy cannot be imported.", solution="Install numpy (see http://numpy.scipy.org/, you probably should also install SciPy).", \
                             addInfo="The numpy package provides array support for Python and is indispensable in many scientific applications."))

# Check whether fitting modules can be imported
_scoImport = ic.check["scipy"]
_pymcImport = ic.check["pymc"]
_mplImport = ic.check["matplotlib"]

from PyAstronomy.funcFit.utils import *
from .modelRebin import turnIntoRebin, _ModelRebinDocu
from .onedfit import *
from .gauss1d import *
from .gauss2d import *
from .params import *
from .sinexp1d import *
from .polyModel import *
from .cauchyLorentz1d import *
# from .anneal import *
from .syncFit import SyncFitContainer
from .coordinateGrid import *
from .circle2d import *
from .extFitter import NelderMead
from .fufDS import FufDS

if ic.check["scipy"]:
  from .voigt1d import *
else:
  PE.warn(PE.PyARequiredImport("No Voigt profile available, because scipy cannot be imported", \
                               solution="Install scipy."))

def __question(tf):
  if tf: return "yes"
  return "no"

def status():
  print("Status of funcFit:")
  print("--------------------------")
  print("Is scipy.optimize available? ", __question(_scoImport))
  if not _scoImport:
    print("  The minimization algorithms from SciPy cannot be used.")
    print("  Install SciPy to make them available (see http://www.scipy.org/).")
  print("Is pymc available? ", __question(_pymcImport))
  if not _pymcImport:
    print("  No Markov-Chain Monte-Carlo sampling via pymc will be possible.")
    print("  Install pymc to enable it (see http://code.google.com/p/pymc/).")
  else:
    import pymc
    print("  pymc is available in version: ", pymc.__version__)
    try:
      majorversion = int(pymc.__version__.split('.')[0])
      if majorversion > 2:
        print("  NOTE: Sampling is only supported with pymc 2.x (not 3.x)!")
    except:
      print("  Could not parse pymc version string...")
  print("Is emcee available? ", __question(ic.check["emcee"]))
  if not ic.check["emcee"]:
    print("  No Markov-Chain Monte-Carlo sampling via emcee will be possible.")
    print("  Install emcee to enable it (see http://dan.iel.fm/emcee/current/).")

__all__ = ["cauchyLorentz1d", "gauss1d", "modelRebin", "onedfit", "params", "voigt1d", "sinexp1d", "polyModel",
           "BoltzmannSAS", "PyAAnneal"]

