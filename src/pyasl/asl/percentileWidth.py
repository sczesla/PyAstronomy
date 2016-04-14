import numpy as np
import copy

def percentileWidth(x, y, rng=[], percentile=0.1):
  """
    Measure the width of an emission line.
    
    Parameters
    ----------
    x, y : arrays
        The abscissa and ordinate of the data.
        x should be of increasing value
    rng : array
        left and right bound of the region to consider
    percentile : float
        Fraction of the peak
        
    Returns
    -------
      w0, w1 - floats
        Can be np.nans if the percentile is not reached within the speicified range
  """
  
  # Check if len(x) == len(y)
  if len(x) != len(y):
    raise(PE.PyAValError("x and y need to have the same length.", \
                         solution="Check the lengths of the input arrays.", \
                         where="percentileWidth"))
  # Check if x is increasing  
  if not all(x[i] <= x[i+1] for i in xrange(len(x)-1)):
    raise(PE.PyAValError("x-axis should be sorted",
                         where="percentileWidth"))
  
  gi = np.where(np.logical_and(x>rng[0], x<=rng[1]))[0]
  mx = np.nanmax(y[gi])
      
  
  i0 = np.where(y == mx)[0]
  if len(i0)!=1:
    return np.nan, np.nan
  i1 = copy.copy(i0)
  x0 = x[i0]      
  x1 = x[i0]      

  while x0>rng[0] and i0>0 and y[i0]>percentile * mx:
      i0-=1
      x0 = x[i0]      
  while x1<rng[1] and i1<len(x)-1 and y[i1] > percentile * mx: 
      i1+=1
      x1 = x[i1]
      
  if y[i0] > percentile * mx:
      x[i0] = np.nan
  if y[i1] > percentile * mx:
      x[i1] = np.nan
  return x[i0], x[i1]
      
          
  
  