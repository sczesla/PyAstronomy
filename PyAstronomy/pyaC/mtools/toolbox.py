from __future__ import division
import numpy as np
import scipy.special as ss

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
  return (r/np.pi)*180.


def farat(x, y):
    """
    Compute ratio of factorials.
    
    Computes x!/y! via ln(x!) - ln(y!) to avoid
    numerical overflow.
    
    Parameters
    ----------
    x : int, float
        The factorial of x is the numerator
    y : int, float
        The factorial of y is the denumerator
    
    Returns
    -------
    Ratio : float
        The ratio x!/y! (not the logarithm).
    """
    lnr = ss.gammaln(x+1) - ss.gammaln(y+1)
    return np.exp(lnr)

