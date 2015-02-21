import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

class BallesterosBV_T:
  """
    Black-body based conversion between effective temperature and B-V color.
    
    Ballesteros 2012 (EPL 97, 34008) present a conversion between
    effective temperature and B-V color index based on a black body
    spectrum and the filter functions.
  """
  
  def __init__(self):
    self._a = 0.92
    self._b = 1.7
    self._c = 0.62
    self._T0 = 4600.0
  
  def bv2T(self, bv):
    """
      Convert B-V color into temperature.
      
      Parameters
      ----------
      bv : float
          B-V color index [mag]
      
      Returns
      -------
      T : float
          Temperature [K]
    """
    T = self._T0*(1.0/(self._a*bv + self._b) + 1.0/(self._a*bv + self._c))
    return T
  
  def t2bv(self, T):
    """
      Convert temperature into B-V color.
      
      Parameters
      ----------
      T : float
          Temperature in K.
      
      Returns
      -------
      bv : float
          B-V color index [mag].
    """
    z = T/self._T0
    ap = z*self._a**2
    bp = self._a*self._c*z + self._b*self._a*z - 2.0*self._a
    cp = self._b*self._c*z -self._c - self._b
    
    sqrtarg = bp**2 - 4.0*ap*cp
    # By comparison with a BB it can be verified that
    # the physical solution is this one
    bv1 = (-bp + np.sqrt(sqrtarg))/(2.0*ap)
    return bv1
    
    
def bv2T_Ballesteros(bv):
  b = BallesterosBV_T()
  return b.bv2T(bv)
bv2T_Ballesteros.__doc__ = BallesterosBV_T.bv2T.__doc__


def t2bv_Ballesteros(T):
  b = BallesterosBV_T()
  return b.t2bv(T)
t2bv_Ballesteros.__doc__ = BallesterosBV_T.t2bv.__doc__