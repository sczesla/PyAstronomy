from __future__ import division
import numpy

class SaniBase:
  
  def mrd(self, x, y):
    """
      Returns maximum relative deviation = max(abs(x-y)/x)
      
      Parameters:
       - `x`, `y` - float or array 
    """
    return numpy.max(numpy.abs(x-y)/x)