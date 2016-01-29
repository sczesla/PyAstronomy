from __future__ import division
import re
import os.path
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class _FDs:
  
  def __init__(self, fn):
    self.fds = {}
    with open(fn) as f:
      for line in f:
        r = re.match("^\s*(\d+)\s+(\d+)\s+(.*)$", line)
        if r is None:
          continue
        derivative = int(r.group(1))
        errOrd = int(r.group(2))
        # Extract finite differences
        fdsStrings = r.group(3).split()
        coeffs = np.zeros(len(fdsStrings))
        offsets = np.arange(len(fdsStrings), dtype=np.int) - len(fdsStrings)//2
        for i, s in enumerate(fdsStrings):
          r = re.match("(-)?(\d+)(/(\d+))?", s)
          coeffs[i] = float(r.group(2))
          if r.group(4) is not None:
            coeffs[i] /= float(r.group(4))
          if r.group(1) is not None:
            coeffs[i] *= -1.0
        self.fds[(derivative, errOrd)] = (offsets, coeffs)
  
  def getFD(self, derivative, errOrd):
    """
    """
    if not (derivative, errOrd) in self.fds:
      raise(PE.PyAValError("No finite differences for derivative " + str(derivative) + \
                           ", error order " + str(errOrd), \
                           where="getFD", \
                           solution="Check derivative and error order."))
    return self.fds[(derivative, errOrd)]
          
          
         

def diffCFD(x, y, derivative, errOrd=2):
  """
    Calculate the numerical derivative using central finite distances.
    
    Parameters
    ----------
    x : array
        The abscissa values.
    y : array
        The ordinate values.
    derivative : int
        The derivative to calculate. Use, e.g., 1 to obtain
        the first derivative.
    errOrd : int
        The order of the error term. The default is 2. Other
        possibilities are 4, 6, and 8. Note the finite
        differences may not be available for all combinations.
    
    Returns
    -------
    indices : array of int
        The indices with respect to the input abscissa to which
        the derivative applies.
    der : array
        The derivative.
  """
  cfd = _FDs(os.path.join(os.path.dirname(__file__), "centralFiniteDiff.dat"))
  offs, fds = cfd.getFD(derivative, errOrd)
  if len(x) < len(offs):
    raise(PE.PyAValError("The length of the input is insufficient.", \
                         where="diffCFD"))
  h = (x[1]-x[0])
  if np.any(np.abs((x[1:]-x[0:-1])/(x[1]-x[0]) - 1.0) > 1e-6):
    raise(PE.PyAValError("The abscissa needs to be equally spaced.", \
                         where="diffCFD"))
  
  imin = abs(min(offs))
  imax = len(x) - imin
  der = np.zeros(imax-imin)
  for i in smo.range(imin, imax):
    der[i-imin] = np.sum(y[i+offs] * fds) / h**derivative
  
  return list(range(imin, imax)), der
  
  
  
  