# -*- coding: utf-8 -*-
import numpy
import PyAstronomy.funcFit as fuf
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.modelSuite.XTran import _ZList

try:
  from . import occultquad
  _importOccultquad = True
except ImportError:
  _importOccultquad = False

try:
  from . import occultnl
  _importOccultnl = True
except ImportError:
  _importOccultnl = False

    
class MandelAgolLC(_ZList, fuf.OneDFit):
  """
    Analytical transit light-curves using the formulae provided by Mandel & Agol 2002.
    
    .. note :: The computation of transit light curves
               is done using the external *occultquad* FORTRAN library.
               This library must be compiled manually using SciPy's f2py
               wrapper (http://www.scipy.org/F2py). Simply go to the
               *forTrans* directory of the source distribution of PyAstronomy,
               then invoke
               
               f2py -c occultquad.pyf occultquad.f
               
               f2py -c occultnl.pyf occultnl.f
               
               That's it. 
    
    :Model parameters:
    
    The set of parameters specifying this model depends on: the type
    of orbit chosen (circular or keplerian) and the type of limb darkening
    chosen (quadratic or non-linear).
    
    More information on the Keplerian orbit can be found here: :ref:`keplerorbitpyasl`
    
    *Orbital model parameters (circular orbit)*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `T0` - Time offset of transit center.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).

    *Orbital model parameters (Keplerian orbit)*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).
      - `tau` - Time of periapsis passage.
      - `Omega` - Longitude of the ascending node [deg].
      - `w` - Argument of periapsis [deg]. Note that the longitude if periapsis is given by Omega+w.
      - `e` - Orbital eccentricity (0-1).

    *Limb darkening parameters (quadratic)*:
      - `linLib` - Linear limb-darkening coefficient.
      - `quadLimb` - Quadratic limb-darkening coefficient.
    
    *Limb darkening parameters (non-linear)*:
      - `a1` - Non-Linear limb-darkening coefficient.
      - `a2` - Non-Linear limb-darkening coefficient.
      - `a3` - Non-Linear limb-darkening coefficient.
      - `a4` - Non-Linear limb-darkening coefficient.

    :Limb-darkening laws:
    
    The quadratic limb-darkening law is given by:
    
    .. math :: \\frac{I(\mu)}{I(1)}= 1 - linLimb \\times (1-\mu) - quadLimb \\times (1-\mu)^2
    
    The non-linear limb-darkening law is given by:
    
    .. math :: \\frac{I(\mu)}{I(1)}= 1 - \\sum_{n=1}^{4}{a_n(1-\mu^{n/2})}

    .. warning::
        Time units have to be consistent.
        
    Parameters
    ----------
    orbit : string, {"circular", "keplerian"}, optional
        Determines whether a circular or full keplerian
        planetary orbit is used in the calculations. The
        default is a circular orbit.
    ld : string, {"quad", "nl"}
        Determines whether quadratic or non-linear
        limb darkening is used. The default is quadratic
        limb darkening.
    collCheck : boolean, optional
        If True (default), the model will check whether there
        is a physical collision between the star and the planet
        on evaluating the model and raises an exception when there
        is one.
  """
  
  def __init__(self, orbit="circular", ld="quad", collCheck=True):
    if not orbit in ["circular", "keplerian"]:
      raise(PE.PyAValError("Invalid option for orbit: " + str(orbit), \
                           soltuion="Use either 'circular' or 'keplerian'."))
    if not ld in ["quad", "nl"]:
      raise(PE.PyAValError("Invalid option for orbit: " + str(ld), \
                           soltuion="Use either 'quad' or 'nl'."))      
    _ZList.__init__(self, orbit, collCheck)
    
    if (not _importOccultquad) and (ld == "quad"):
        raise(PE.PyARequiredImport("Could not import required shared object library 'occultquad.so'",\
                                   solution=["Invoke PyA's install script (setup.py) with the --with-ext option.",
                                             "Go to 'forTrans' directory of PyAstronomy and invoke\n    f2py -c occultquad.pyf occultquad.f"]
                                   ))
    
    if (not _importOccultnl) and (ld == "nl"):
        raise(PE.PyARequiredImport("Could not import required shared object library 'occultquad.so'",\
                                   solution=["Invoke PyA's install script (setup.py) with the --with-ext option.",
                                             "Go to 'forTrans' directory of PyAstronomy and invoke\n    f2py -c occultnl.pyf occultnl.f"]
                                   ))
    
    if orbit == "circular":
      plist = ["p", "a", "i", "T0", "per", "b"]
    else:
      # It is an elliptical orbit
      plist = ["p", "a", "i", "per", "b", \
               "e", "Omega", "tau", "w"]
    
    if ld == "quad":
      plist.extend(["linLimb", "quadLimb"])
    else:
      # Non-linear LD
      plist.extend(["a1", "a2", "a3", "a4"])
    
    fuf.OneDFit.__init__(self, plist)
    self.freeze(plist)
    self.setRootName("Occultquad")
    
    if orbit == "keplerian":
      self["w"] = -90.0
      
    self._orbit = orbit
    self._ld = ld
    
  def evaluate(self, time):
    """ 
     Calculate a light curve according to the analytical models
     given by Mandel and Agol.
        
     Parameters
     ----------
     time : array
         An array of time points at which the light curve
         shall be calculated.

     Returns
     -------
     Model : array
         The analytical light curve is stored in the property `lightcurve`.
    """

    # Translate the given parameters into an orbit and, finally,
    # into a projected, normalized distance (z-parameter)
    if self._orbit == "circular":
      self._calcZList(time - self["T0"])
    else:
      # The orbit is keplerian
      self._calcZList(time)

    # Use occultquad Fortran library to compute flux decrease
    df = numpy.zeros(len(time))
    if len(self._intrans) > 0:
      if self._ld == "quad":
        # Use occultquad Fortran library to compute flux decrease
        result = occultquad.occultquad(self._zlist[self._intrans],self["linLimb"],self["quadLimb"], \
                                       self["p"],len(self._intrans))
      else:
        result = occultnl.occultnl(self["p"],self["a1"],self["a2"],self["a3"], \
                                   self["a4"],self._zlist[self._intrans])
      df[self._intrans] = (1.0 - result[0])
      
    self.lightcurve = (1.-df)*1./(1.+self["b"]) + self["b"]/(1.0+self["b"])
    
    return self.lightcurve

MandelAgolLC_Rebin = fuf.turnIntoRebin(MandelAgolLC)