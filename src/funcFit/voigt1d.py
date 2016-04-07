# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy
from .onedfit import OneDFit
import scipy.special as sps
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

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
        
        
      z = ((x-self["mu"]) + (1.j)*abs(self["al"]))/(numpy.sqrt(2.0) * abs(self["ad"]))
      y = self["A"] * numpy.real(sps.wofz(z)) 
      y /= (abs(self["ad"]) * numpy.sqrt(2.0*numpy.pi))
      y += x * self["lin"] + self["off"]
  
      y[numpy.where(numpy.isnan(y))] = 0.0
    except FloatingPointError as fpe:
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



class MultiVoigt1d(OneDFit):
  """
    Multicomponent Voigt with a single linear continuum component.
    
    The parameters are the same as for *Voigt1d*,
    except that all are extended by a number specifying the
    component to which they belong. Therefore, they read, e.g.,
    `A1`, `mu2`, and so on; only `off` and `lin`
    remain unnumbered.

    Parameters
    ----------
    n : int
        The number of Voigt components.
  """
  
  def __init__(self, n):
    # Number of components
    self.n = n
    # Building parameter list
    params = ["off", "lin"]
    for i in range(n):
      p = str(i+1)
      params.extend(["A"+p, "mu"+p, "al"+p, "ad"+p])
    OneDFit.__init__(self, params)
    self.setRootName("MultiVoigt")
    # Use Voigt1d for evaluation
    self._v1d = Voigt1d()
    self._v1d["off"] = 0.0
    self._v1d["lin"] = 0.0
  
  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    # Assign 'continuum'
    y = self["off"] + self["lin"]*x
    # Add Voigt lines
    for i in smo.range(self.n):
      p = str(i+1)
      self._v1d.assignValues({"A":self["A"+p], "al":self["al"+p], "ad":self["ad"+p], "mu":self["mu"+p]})
      y += self._v1d.evaluate(x)
    return y
  
  def evalComponent(self, x, p):
    """
      Evaluate the model considering only a single component.
      
      Parameters
      ----------
      x : array
          The abscissa.
      p : int
          Component number (starts with one).
      
      Returns
      -------
      Single component model : array
          The model considering a single component. Note that the
          linear continuum is included.
    """
    if p > 0 and p <= self.n:
      p = str(p)
      y = self["off"] + self["lin"] * x
      self._v1d.assignValues({"A":self["A"+p], "al":self["al"+p], "ad":self["ad"+p], "mu":self["mu"+p]})
      y += self._v1d.evaluate(x)
      return y
    else:
      raise(PE.PyAValError("No such component (no. "+str(p)+")", where="MultiVoigt1d::evalComponent", \
            solution="Use value between 1 and "+str(self.n)))

