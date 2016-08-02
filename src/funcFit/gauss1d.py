# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy
from .onedfit import OneDFit
from PyAstronomy.pyaC import pyaErrors as PE

class GaussFit1d(OneDFit):
  """
    Implements a one dimensional Gaussian.
    
    The functional form is:
    
    .. math:: \\frac{A}{\\sqrt{2\\pi\\sigma^2}}e^{-(x-\\mu)^2/(2\\sigma^2)} + x \\times lin + off
    
    Here, `lin` and `off` denote the linear and the offset term.
    
    *Fit parameters*:
     - `A` - Amplitude (the area of the Gaussian)
     - `mu` - Center of the Gaussian
     - `sig` - Standard deviation
     - `off` - Offset
     - `lin` - Linear term
  """
  
  def __init__(self):
    OneDFit.__init__(self, ["A", "mu", "sig", "off", "lin"])
    self.setRootName("Gaussian")

  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    if self["sig"] == 0.0:
      raise(PE.PyAValError("Width of Gaussian must be larger than zero.", \
                           solution="Change width ('sig')."))
    y = self["A"]/numpy.sqrt(2.0*numpy.pi*self["sig"]**2) * \
        numpy.exp(-(self["mu"] - x)**2/(2.0 * self["sig"]**2)) \
        + self["off"] + (self["lin"] * x)
    return y


class MultiGauss1d(OneDFit):
  """
    A multicomponent Gaussian with a single linear continuum component.
    
    The parameters are the same as for the *GaussFit1d*,
    except that all receive a number specifying the Gaussian
    component to which they belong. Therefore, they are, e.g.,
    named `A1`, `mu2`, and so on, only `off` and `lin`
    remain unnumbered.

    Parameters
    ----------
    n : int
        The number if Gaussian components.
  """
  
  def __init__(self, n):
    self.n = n
    params = ["off", "lin"]
    for i in range(n):
      p = str(i+1)
      params.extend(["A"+p, "mu"+p, "sig"+p])
    OneDFit.__init__(self, params)
    self.setRootName("MultiGauss")
  
  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    y = self["off"] + self["lin"]*x
    for i in range(self.n):
      p = str(i+1)
      y += self["A"+p]/numpy.sqrt(2.0*numpy.pi*self["sig"+p]**2) * \
          numpy.exp(-(self["mu"+p] - x)**2/(2.0 * self["sig"+p]**2))
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
      single component model : array
          The model considering a single component.
    """
    if p > 0 and p <= self.n:
      p = str(p)
      y = self["off"] + self["lin"] * x
      y += self["A"+p] / numpy.sqrt(2.0*numpy.pi*self["sig"+p]**2) * \
           numpy.exp(-(self["mu"+p] - x)**2/(2.0 * self["sig"+p]**2))
      return y
    else:
      raise(PE.PyAValError("No such component (no. "+str(p)+")", where="MultiGauss1d::evalComponent", \
            solution="Use value between 1 and "+str(self.n)))

