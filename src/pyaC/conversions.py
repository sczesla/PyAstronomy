from __future__ import division
import numpy as np

def degtorad(d):
  """
    Convert degrees into radians.
    
    Parameters
    ----------
    d : float or array
        Angle in degrees.
    
    Returns
    -------
    Angle : float or array
        The angle converted into radians.
  """
  return (d/180.0)*np.pi

def radtodeg(r):
  """
    Convert radians into degrees.
    
    Parameters
    ----------
    d : float or array
        Angle in radians.
    
    Returns
    -------
    Angle : float or array
        The angle converted into degrees.
  """
  return (r/np.pi)*180.0