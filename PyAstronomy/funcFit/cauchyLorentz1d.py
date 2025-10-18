# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy
from .onedfit import OneDFit

class CauchyLorentz1d(OneDFit):
  """
    A Cauchy-Lorentz profile
    
    *Fit parameters*:
     - `A` - Amplitude
     - `g` - Scale parameter (usually gamma)
     - `mu` - Center
     - `off` - A constant offset
     - `lin` - A linear contribution
  """

  def __init__(self):
    OneDFit.__init__(self, ["A", "g", "mu", "off", "lin"])
    self.setRootName("CauLor")

  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    y = self["A"] * (self["g"]/numpy.pi / ((x-self["mu"])**2 + self["g"]**2)) + \
        (self["lin"] * x) + self["off"]
    return y

