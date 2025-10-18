# -*- coding: utf-8 -*-
import numpy
import PyAstronomy.funcFit as fuf
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.modelSuite.XTran import _ZList

try:
  from . import occultnl
  _importOccultnl = True
except ImportError:
  _importOccultnl = False
  
    
class MandelAgolNLLC(_ZList, fuf.OneDFit):
  """
    Calculate and fit analytical transit light-curves using the formulae
    provided by Mandel & Agol 2002. This model uses the non-linear
    limb-darkening prescription by Claret 2011:
    
    .. math :: \\frac{I(\mu)}{I(1)}= 1 - \\sum_{n=1}^{4}{a_n(1-\mu^{n/2})}
    
    .. note :: The computation of transit light curves
               is done using the external *occultnl* FORTRAN library.
               This library must be compiled manually using SciPy's f2py
               wrapper (http://www.scipy.org/F2py). Simply go to the
               *forTrans* directory of the source distribution of PyAstronomy,
               then invoke
               
               f2py -c occultnl.pyf occultnl.f
               
               
    
    *Model parameters*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `a1` - Non-Linear limb-darkening coefficient.
      - `a2` - Non-Linear limb-darkening coefficient.
      - `a3` - Non-Linear limb-darkening coefficient.
      - `a4` - Non-Linear limb-darkening coefficient.
      - `T0` - Time offset of transit center.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).
    
  """
  
  def __init__(self):
    _ZList.__init__(self, "circular")
    if not _importOccultnl:
        raise(PE.PyARequiredImport("Could not import required shared object library 'occultquad.so'",\
                                   solution=["Invoke PyA's install script (setup.py) with the --with-ext option.",
                                             "Go to 'forTrans' directory of PyAstronomy and invoke\n    f2py -c occultnl.pyf occultnl.f"]
                                   ))

    fuf.OneDFit.__init__(self, ["p", "a", "i", "a1", "a2", "a3", "a4", "T0", "per", "b"])
    self.freeze(["p", "a", "i", "a1", "a2", "a3", "a4", "T0", "per", "b"])
    self.setRootName("MandelAgolNL")
    
    self.__zlist=None
    
  def evaluate(self, time):
    """ 
     Calculate a light curve according to the analytical models
     given by Mandel & Agol 2002.
        
     Parameters
     ----------
     time : array
         An array of time points at which the light curve
         shall be calculated.
        
     .. note:: time = 0 -> Planet is exactly in the line of sight (phase = 0).

     Returns
     -------
     Model : array
         The analytical light curve is stored in the property `lightcurve`.
    """

    # Translate the given parameters into an orbit and, finally,
    # into a projected, normalized distance (z-parameter)
    self._calcZList(time - self["T0"])

    # Use occultquad Fortran library to compute flux decrease
    result = occultnl.occultnl(self["p"],self["a1"],self["a2"],self["a3"], \
                               self["a4"],self._zlist[self._intrans])
    
    df = numpy.zeros(len(time))
    df[self._intrans] = (1.0 - result[0])
 
    self.lightcurve = (1.-df)*1./(1.+self["b"]) + self["b"]/(1.0+self["b"])
    
    return self.lightcurve


MandelAgolNLLC_Rebin = fuf.turnIntoRebin(MandelAgolNLLC)