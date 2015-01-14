from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyaC import invertIndexSelection
import numpy as np

def intep(x, y, xinter, boundsError=True, fillValue=None):
  """
    The INTEP interpolation algorithm
    
    The INTEP interpolation algorithm is described
    by Hill 1982, PDAO 16, 67 ("Intep - an Effective
    Interpolation Subroutine"). The implementation
    at hand is based on the FORTRAN code stated
    therein.
    
    The aim of the algorithm is to imitate the curve
    "an experienced scientist" would draw through a
    given set of points. 
    
    Parameters
    ----------
    x : array
        Independent values.
    y : array
        Dependent values.
    xinter : array
        Values at which to interpolate the tabulated
        data given by `x` and `y`.
    boundsError : boolean, optional
        If True, an exception will be raised if values
        need to be extrapolated beyond the limits of
        the given tabulated data. Values beyond the
        limits are simply replaced with the closest
        valid value available, which might not be a
        good approximation. Set this flag to False
        suppress the exception.
    fillValue : float, optional
        If given (i.e., not None), this value will be
        used to represent values outside of the given
        bounds. Note that `boundsError` must be set
        False for this to have an effect. For instance,
        use np.NaN.
    
    Returns
    -------
    Interpolated values : array
        Interpolated values at the locations
        specified by `xinter`.
  """
  # Check whether x-array is sorted
  if not np.all((x[1:] - x[0:-1]) > 0.0):
    raise(PE.PyAValError("The array of independent values (`x`) is not sorted in (strictly) " + \
                         "ascending order, but it needs to be.", \
          solution=["Sort the arrays `x` (and `y` accordingly), e.g., using numpy.argsort.", \
                    "Check whether there are duplicate values in `x`."]))
  
  # Create result array
  result = np.zeros(len(xinter))
  # Treat extrapolation points
  ilow = np.where(xinter < min(x))[0]
  if fillValue is None:
    # Use first point beyond limit
    result[ilow] = y[0]
  else:
    result[ilow] = fillValue
  iup = np.where(xinter > max(x))[0]
  if fillValue is None:
    # Use last point beyond limit
    result[iup] = y[-1]
  else:
    result[iup] = fillValue
  
  # No extrapolation (noepo)
  noepo = invertIndexSelection(xinter, np.concatenate((ilow, iup)))
  
  if (len(noepo) < len(xinter)) and boundsError:
    raise(PE.PyAValError("There are " + str(len(xinter) - len(noepo)) + " points, which need to be extrapolated. " + \
          "Attention: Extrapolation is simply done by using the closest valid value.", \
          solution=["Use only interpolation points (`xinter`) within the valid range.",
                    "Set the `boundsError` flag to False (will suppress this exception)."]))
  
  # Loop over all indices of xinter, which were
  # not already subject to extrapolation.
  for i in noepo:
    xp = xinter[i]
    # Index of the first entry in x larger (or equal) than
    # `val`
    infl = np.where(x >= xp)[0][0]
    infl -= 1
    lp1 = 1.0/(x[infl] - x[infl+1])
    lp2 = -lp1
    if infl <= 0:
      # Special treatment for first point
      fp1 = (y[1] -y[0]) / (x[1] - x[0])
    else:
      fp1 = (y[infl+1] - y[infl-1]) / (x[infl+1] - x[infl-1])
    if infl >= (len(x) - 2):
      # Special treatment for last points
      fp2 = (y[-1] - y[-2]) / (x[-1] - x[-2])
    else:
      fp2 = (y[infl+2] - y[infl]) / (x[infl+2] - x[infl])
    
    xpi1 = xp - x[infl+1]
    xpi  = xp - x[infl]
    l1 = xpi1 * lp1
    l2 = xpi * lp2
    result[i] = y[infl]*(1. - 2.*lp1*xpi)*l1**2 + \
                y[infl+1]*(1. - 2.*lp2*xpi1)*l2**2 + \
                fp2*xpi1*l2**2 + fp1*xpi*l1**2
  return result