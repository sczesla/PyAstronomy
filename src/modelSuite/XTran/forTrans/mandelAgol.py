# -*- coding: utf-8 -*-
import numpy
import PyAstronomy.funcFit as fuf
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.modelSuite.XTran import _ZList

try:
  import occultquad
  _importOccultquad = True
except ImportError:
  _importOccultquad = False
  
    
class MandelAgolLC(_ZList, fuf.OneDFit):
  """
    Calculate and fit analytical transit light-curves using the formulae
    provided by Mandel & Agol 2002.
    
    .. note :: The computation of transit light curves
               is done using the external *occultquad* FORTRAN library.
               This library must be compiled manually using SciPy's f2py
               wrapper (http://www.scipy.org/F2py). Simply go to the
               *forTrans* directory of the source distribution of PyAstronomy,
               then invoke
               
               f2py -c occultquad.pyf occultquad.f
               
               That's it. 
    
    *Model parameters*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `linLib` - Linear limb-darkening coefficient.
      - `quadLimb` - Quadratic limb-darkening coefficient.
      - `T0` - Time offset of transit center.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).

    This class inherits the functionality of funcFit's OneDFit object.
    You can, therefore, set and get the parameter using the brackets:
      e.g., pallc["p"] = 0.12345

    .. warning::
        Time units have to be consistent.
  """
  
  def __init__(self):
    _ZList.__init__(self, "circular")
    if not _importOccultquad:
        raise(PE.PyARequiredImport("Could not import required shared object library 'occultquad.so'",\
                                   solution=["Invoke PyA's install script (setup.py) with the --with-ext option.",
                                             "Go to 'forTrans' directory of PyAstronomy and invoke\n    f2py -c occultquad.pyf occultquad.f"]
                                   ))

    fuf.OneDFit.__init__(self, ["p", "a", "i", "linLimb", "quadLimb", "T0", "per", "b"])
    self.freeze(["p", "a", "i", "linLimb", "quadLimb", "T0", "per", "b"])
    self.setRootName("Occultquad")
    
    self.__zlist=None
    
  def evaluate(self, time):
    """ 
     Calculate a light curve according to the analytical models
     given by Pal 2008.
        
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
    result = occultquad.occultquad(self._zlist[self._intrans],self["linLimb"],self["quadLimb"], \
                                   self["p"],len(self._intrans))
    df = numpy.zeros(len(time))
    df[self._intrans] = (1.0 - result[0])
    
    self.lightcurve = (1.-df)*1./(1.+self["b"]) + self["b"]/(1.0+self["b"])
    
    return self.lightcurve

MandelAgolLC_Rebin = fuf.turnIntoRebin(MandelAgolLC)