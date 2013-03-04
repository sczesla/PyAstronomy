import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

def getCardinalPoint(azimuth):
  """
    Get the cardinal point (North, East, South, or West) for a given azimuth.
    
    Here, the azimuth is measured from North to East.
    N = (315, 45] deg
    E = (45, 135] deg
    S = (135, 225] deg
    W = (225, 315] deg

    Parameters
    ----------
    azimuth : float
        Azimuth of an object in deg (0-360).

    Returns
    -------
    Cardinal point : string, {N,E,S,W}
        Returns the cardinal point (N, E, S, W) corresponding
        to the given azimuth of the object.
  """  
  
  if (azimuth < 0.0) or (azimuth > 360.):
    raise(PE.PyAValError("The azimuth needs to be a number in the range [0,360]. Given azimuth: " + str(azimuth), \
          where="getCardinalPoint", \
          solution="Specify a number between 0 and 360"))
  
  if np.logical_and(azimuth > 315., azimuth <= 360.): return "N"
  if np.logical_and(azimuth >= 0., azimuth <= 45.):   return "N"
  if np.logical_and(azimuth > 45., azimuth <= 135.):  return "E"
  if np.logical_and(azimuth > 135., azimuth <= 225.): return "S"
  if np.logical_and(azimuth > 225., azimuth <= 315.): return "W"