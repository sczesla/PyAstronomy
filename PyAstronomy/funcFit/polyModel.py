# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy
from .onedfit import OneDFit
import six.moves as smo

class PolyFit1d(OneDFit):
  """
    Implements a polynomial fit.
    
    *Fit parameters*:
      - cn - Here n is a number indicating degree
             (e.g., c0 + c1*x + c2*x**2 ...)
    
    If `xoff` is specified, the polynomial will be
    evaluated at the points `x-xoff`. This can be
    useful to suppress correlation.
  """
  
  def __init__(self, degree, xoff=0.0):
    coeffs = []
    self.degree = degree
    self.xoff = xoff
    for i in smo.range(degree+1):
      coeffs.append("c"+str(i))
    OneDFit.__init__(self, coeffs)
    self.setRootName("Poly")

  def asNPPoly(self):
    """
      Construct a numpy.poly1d object from the coefficients.
      
      Returns
      -------
      Polynomial : numpy.poly1d
          A numpy polynomial with the coefficients stored
          by the fit object. 
    """
    poly = []
    for i in smo.range(self.degree, -1, -1):
      poly.append(self["c"+str(i)])
    return numpy.poly1d(poly)

  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    return self.asNPPoly()(x-self.xoff)


class ConstantFit1d(OneDFit):
  """
    Implements a constant model.
    
    *Fit parameters*:
      - `c` - The constant.
  """
  def __init__(self):
    OneDFit.__init__(self, ["c"])
    self.setRootName("Constant")
  
  def evaluate(self, x):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    return numpy.ones(x.size) * self["c"]