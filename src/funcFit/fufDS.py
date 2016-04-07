from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import copy
import numpy as np

class FufDS:
  """
    FuncFit data set.
    
    Parameters
    ----------
    x : array
        Specify the abscissa.
    y : array
        Specify the data.
    yerr : array, optional
        The error of the data.
  
  """
  
  def __init__(self, x, y, yerr=None):
    if not isinstance(y, np.ndarray):
      raise(PE.PyAValError("`y` must be a numpy.ndarray. Current type is '" + str(type(y) + "'.", \
            solution="Make `y` a numpy array.")))
    self.y = y.copy()
    if yerr is not None:
      if not isinstance(yerr, np.ndarray):
        raise(PE.PyAValError("`yerr` must be a numpy.ndarray. Current type is '" + str(type(y) + "'.", \
              solution="Make `yerr` a numpy array.")))
        if y.shape != yerr.shape:
          raise(PE.PyAValError("`y` and `yerr` must have the same shape.", \
                solution="Check and adapt shapes of arrays."))
      self.yerr = yerr.copy()
    else:
      self.yerr = None
    if isinstance(x, np.ndarray):
      self.x = x.copy()
    else:
      self.x = copy.copy(x)
  
  def yerrDefined(self):
    """
      Check whether yerr is defined or None.
      
      Returns
      -------
      Defined : boolean
          True, if yerr is not None.
    """
    return not (self.yerr is None)
  
  def xyyerrDefined(self):
    """
      Check whether x, y, and yerr are not None:
      
      Returns
      -------
      Defined : boolean
          True, if x, y, and yerr are not None. 
    """
    xnn = not (self.x is None)
    ynn = not (self.y is None)
    yenn = self.yerrDefined()
    return (xnn & ynn & yenn)