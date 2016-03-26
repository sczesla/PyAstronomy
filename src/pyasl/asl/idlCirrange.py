from __future__ import print_function, division
from .idlMod import idlMod
import numpy as np

def cirrange(x, radians=False):
  """
    Force angle into range.
    
    Emulates IDL's `cirrange`.
    
    Parameters
    ----------
    x : float or array
        The angle(s).
    radians : boolean, optional
        If True, angle will be forced into
        0-2pi range. If False (default),
        the 0-360 deg range will be applied.
    
    Returns
    -------
    Modified angle : float or array
        The angle forced into the 0-360 deg
        or 0-2pi radians range.
  """
  if radians:
    m = 2.0 * np.pi
  else:
    m = 360.0
  return x % m
    