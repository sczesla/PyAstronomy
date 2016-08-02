from __future__ import print_function, division
import numpy as np
from .onedfit import OneDFit

class Circle2d(OneDFit):
  """
    Implements a simple, two-dimensional circular orbit.
    
    The functional form is:
    
    .. math:: (x,y) = (r \cos(\omega (t-t_0)), r \sin(\omega (t-t_0)))
    
    where :math:`\omega` is :math:`2\pi/Period`.
    
    *Fit parameters*:
     - `r` - Radius of the circle
     - `per` - Period of orbit
     - `t0` - Starting time of orbit
  """
  
  def __init__(self):
    OneDFit.__init__(self, ["r", "per", "t0"])
    self.setRootName("Circle2d")

  def evaluate(self, t):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      x : array
          Specifies the points at which to evaluate the model.
    """
    result = np.zeros((len(t), 2))
    w = 2.0*np.pi / self["per"]
    result[::,0] = abs(self["r"]) * np.cos(w*(t-self["t0"]))
    result[::,1] = abs(self["r"]) * np.sin(w*(t-self["t0"]))
    return result