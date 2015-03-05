from PyAstronomy.pyaC import _ic
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np

def ibtrapz(x, y, x0, x1, iaout=False):
  """
    Use the trapezoid rule to integrate and interpolate boundary values.
    
    Can be used for integration on tabled data, where the integration boundaries
    are not located on the grid. The values to use at the boundaries are
    determined using linear interpolation.
    
    Parameters
    ----------
    x,y : arrays
        The data.
    x0, x1 : float
        The integration boundaries.
    iaout : boolean, optional
        If True, the arrays used internally for
        integration are also returned. The default
        is False.
    
    Returns
    -------
    Integral : float
        The value of the resulting integral.
    xi, yi : arrays, optional
        Internally used arrays for integration including
        the values derived at the boundaries. Only returned
        if `iaout` is set to True.
  """
  if not _ic.check["scipy"]:
    raise(PE.PyARequiredImport("scipy could not be imported.", \
                               where="ibtrapz", \
                               solution="Please install scipy to use ibtrapz."))
  import scipy.interpolate as sci
  import scipy.integrate as scinteg
  if (x1 < x0):
    raise(PE.PyAValError("x1 must be larger than x0.", \
                         where="ibtrapz", \
                         solution="Exchange boundaries."))
  if np.any(x[1:]-x[0:-1] < 0.0):
    raise(PE.PyAValError("x axis must be in ascending order.", \
                         where="ibtrapz", \
                         solution="Sort values."))
  if (x0 < np.min(x)) or (x1 > np.max(x)):
    raise(PE.PyAValError("One or both integration boundaries are beyond the range of data given.", \
                         where="ibtrapz", \
                         solution="Give more data or adjust integration boundaries."))
  # Interpolate at integration boundaries
  fi = sci.interp1d(x, y)
  y0, y1 = fi(x0), fi(x1)
  # Search data within boundaries
  indi = np.where((x > x0) & (x < x1))[0]
  # Construct arrays to be used in integration
  xi = np.zeros(len(indi)+2)
  yi = np.zeros(len(indi)+2)
  xi[1:-1] = x[indi]
  yi[1:-1] = y[indi]
  xi[0] = x0
  xi[-1] = x1
  yi[0] = y0
  yi[-1] = y1
  
  if iaout:
    return scinteg.trapz(yi, xi), xi, yi
  return scinteg.trapz(yi, xi)
  
    
  
  