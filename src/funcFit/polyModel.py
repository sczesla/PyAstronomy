# -*- coding: utf-8 -*-
import numpy
from onedfit import OneDFit

class PolyFit1d(OneDFit):
  """
    Implements a polynomial fit.
    
    *Fit parameters*:
      - cn - Here n is a number indicating degree
             (e.g., c0 + c1*x + c2*x**2 ...)
  """
  
  def __init__(self, degree):
    coeffs = []
    self.degree = degree
    for i in xrange(degree+1):
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
    for i in xrange(self.degree, -1, -1):
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
    return self.asNPPoly()(x)


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