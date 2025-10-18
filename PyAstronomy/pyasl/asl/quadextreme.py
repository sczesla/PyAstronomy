import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

def quadExtreme(x, y, mode="max", dp=(1,1), exInd=None, fullOutput=False, fullPoint=False):
  """
    Find the extreme (minimum or maximum) by means of a parabolic fit.
    
    This function searches for the maximum (or minimum) in the given
    ordinate values, fits a second-order polynomial to the
    surroundings of that point, and, finally, determines the
    thus approximated abscissa value of the extreme point. 
    
    Parameters
    ----------
    x : array
        Abscissa values.
    y : array
        Ordinate values.
    mode : string, {min, max}, optional
        Determines whether a minimum or a maximum is searched
        for. The default is a maximum.
    dp : tuple with two integers, optional
        Determines the width around the extreme point, which will
        be used in the fitting. The default is one point to the
        right and left, i.e., dp = (1,1).
    extInd : integer, optional
        If given, the function will assume that this is the
        index of the extreme point and not search for it.
    fullOutput : boolean, optional
        If True, the output will also cover the points used in
        the fitting and the resulting polynomial.
    fullPoint : boolean, optional
        If True, the return value `epos` will be a tuple holding
        the abscissa and ordinate values of the extreme point.
        The default is False.
    
    Returns
    -------
    epos : float or tuple
        Position of the extreme point. If `fullPoint` is True,
        a tuple with the abscissa and ordinate values of the
        extreme point.
    mi : int
        The index of the extreme point (maximum or minimum).
    xb : array, optional
        Only returned if `fullOutput` is True. The abscissa
        values used in the polynomial fit. Note the the
        x-value of the extreme point has been subtracted.
    yb : array, optional
        Only returned if `fullOutput` is True. The ordinate
        values used in the polynomial fit.
    p : numpy polynomial, optional
        Only returned if `fullOutput` is True. The best-fit
        polynomial. Note that the fit refers to the `xb` axis,
        where the x-value of the extreme point has been
        subtracted.
  """
  if exInd is None:
    if mode == "max":
      mi = np.argmax(y)
    elif mode == "min":
      mi = np.argmin(y)
    else:
      raise(PE.PyAValError("Unknown mode '" + str(mode) + "'.", \
                           solution="Choose 'min' or 'max'"))
  else:
    mi = exInd

  if ((mi-dp[0]) < 0) or ((mi+dp[0]+1 > len(x))):
    raise(PE.PyAValError("The requested fitting range around the extreme is not covered by the data.", \
                         solution=["Adapt 'dp'.", "Do you need to change the mode?"]))
  
  # Cut out the relevant range
  xb = x[mi-dp[0]:mi+dp[1]+1].copy()
  yb = y[mi-dp[0]:mi+dp[1]+1].copy()
  
  # Subtract the x-value of the extreme
  xm = x[mi]
  xb -= xm
  
  # Fit second-order polynomial
  p = np.polyfit(xb, yb, 2)
  
  # Check whether curvature agrees with mode
  signs = {"max":-1, "min":+1}
  if np.sign(p[0]) != signs[mode]:
    PE.warn(PE.PyAValError("The curvature of the fit does not match the requested mode of '" + str(mode) + "'.", \
                           solution="Did you choose the right mode (min/max)?"))
  
  # Location of the extreme
  epos = -p[1]/(2.0 * p[0])
  # Shift back the location
  epos += xm
  
  if fullPoint:
    epos = (epos, -p[1]**2/(4.0*p[0]) + p[2])
  
  if not fullOutput:
    return epos, mi
  else:
    return epos, mi, xb, yb, p
    
  