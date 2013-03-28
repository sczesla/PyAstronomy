# -*- coding: utf-8 -*-
import numpy
from onedfit import OneDFit
import scipy.special as sps
from PyAstronomy.pyaC import pyaErrors as PE

class Voigt1d(OneDFit):
  """
    Implements a Voigt profile (convolution of Cauchy-Lorentz
    and Gaussian distribution).
    
    .. note:: The profile is implemented so that `al` is half
              the FWHM.
    
    *Fit parameters*:
      - `A` - Amplitude
      - `al` - Scale parameter of the Cauchy-Lorentz distribution
      - `ad` - The width of the Gaussian (usually called sigma)
      - `mu` - Center
      - `off` - Constant offset
      - `lin` - Linear contribution
    
    Notes
    -----
    
    The Voigt profile V is defined as the convolution
      
    .. math::
      V = A\int G(x')L(x-x')dx'
      
    of a Gaussian distribution
      
    .. math::
      G=1/(2 \pi \ ad)\exp(-(x-mu)^2/(2 \ ad^2))
      
    and a Cauchy-Lorentz distribution
      
    .. math::
      L=al/(\pi ((x-mu)^2+al^2)) + lin \ x + off
      
    The Voigt profile is calculated via the real part of
    the Faddeeva function.
    For details, see http://en.wikipedia.org/wiki/Voigt_profile 
    and http://en.wikipedia.org/wiki/Error_function.
  """

  def __init__(self):
    OneDFit.__init__(self, ["A", "al", "ad", "mu", "lin", "off"])
    self.setRestriction({"al":[0.0,None], "ad":[0.0,None]})
    self.setRootName("Voigt")

  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    try:
        
        
      z = ((x-self["mu"]) + (1.j)*self["al"])/(numpy.sqrt(2.0) * self["ad"])
      y = self["A"] * numpy.real(sps.wofz(z)) 
      y /= (self["ad"] * numpy.sqrt(2.0*numpy.pi))
      y += x * self["lin"] + self["off"]
  
      y[numpy.where(numpy.isnan(y))] = 0.0
    except FloatingPointError, fpe:
      raise(PE.PyAFloatingPointError("The following floating point error occurred:\n  " + str(fpe) + "\n" +
                                     "Current Parameter values:\n"+str(self.parameters()),
                                     solution=["Try to rescale/shift your abscissa. For instance, put" + 
                                               "the spectral line you try to fit at position `0`."]))
    return y
  
  def FWHM(self):
      """
        Calculates an approximation of the FWHM.
        
        The approximation is accurate to 
        about 0.03% (see http://en.wikipedia.org/wiki/Voigt_profile).
        
        Returns
        -------
        FWHM : float
            The estimate of the FWHM
      """
      # The width of the Lorentz profile
      fl =  2.0 * self["al"]
      # Width of the Gaussian [2.35 = 2*sigma*sqrt(2*ln(2))]
      fd = 2.35482 * self['ad']
      return 0.5346*fl + numpy.sqrt(0.2166*(fl**2.) + fd**2.)  

