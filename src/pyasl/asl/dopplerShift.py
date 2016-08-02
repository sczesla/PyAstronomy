from __future__ import print_function, division
import scipy.interpolate as sci
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

def dopplerShift(wvl, flux, v, edgeHandling=None, fillValue=None):
  """
    Doppler shift a given spectrum.
    
    A simple algorithm to apply a Doppler shift
    to a spectrum. This function, first, calculates
    the shifted wavelength axis and, second, obtains
    the new, shifted flux array at the old, unshifted
    wavelength points by linearly interpolating.
    
    Due to the shift, some bins at the edge of the
    spectrum cannot be interpolated, because they
    are outside the given input range. The default
    behavior of this function is to return numpy.NAN
    values at those points. One can, however, specify
    the `edgeHandling` parameter to choose a different
    handling of these points.
    
    If "firstlast" is specified for `edgeHandling`,
    the out-of-range points at the red or blue edge
    of the spectrum will be filled using the first
    (at the blue edge) and last (at the red edge) valid
    point in the shifted, i.e., the interpolated, spectrum.
    
    If "fillValue" is chosen for edge handling,
    the points under consideration will be filled with
    the value given through the `fillValue` keyword.
    
    .. warning:: Shifting a spectrum using linear
                interpolation has an effect on the
                noise of the spectrum. No treatment
                of such effects is implemented in this
                function.
    
    Parameters
    ----------
    wvl : array
        Input wavelengths in A.
    flux : array
        Input flux.
    v : float
        Doppler shift in km/s
    edgeHandling : string, {"fillValue", "firstlast"}, optional
        The method used to handle the edges of the
        output spectrum.
    fillValue : float, optional
        If the "fillValue" is specified as edge handling method,
        the value used to fill the edges of the output spectrum.
    
    Returns
    -------
    nflux : array
        The shifted flux array at the *old* input locations.
    wlprime : array
        The shifted wavelength axis.
  """
  # Shifted wavelength axis
  wlprime = wvl * (1.0 + v/299792.458)
  
  fv = np.nan
  if edgeHandling == "fillValue":
    if fillValue is None:
      raise(PE.PyAValError("Fill value not specified", where="pyasl.dopplerShift", \
                           solution="If you request 'fillValue' as edge handling method, you need to specify the 'fillValue' keyword."))
    fv = fillValue
      
  f = sci.interp1d(wlprime, flux, bounds_error=False, fill_value=fv)
  nflux = f(wvl)
  
  if edgeHandling == "firstlast":
    firsts = []
    # Search for first non-NaN value save indices of
    # leading NaN values
    for i in smo.range(len(nflux)):
      if np.isnan(nflux[i]):
        firsts.append(i)
      else:
        firstval = nflux[i]
        break
    # Do the same for trailing NaNs
    lasts = []
    for i in smo.range(len(nflux)-1,0,-1):
      if np.isnan(nflux[i]):
        lasts.append(i)
      else:
        lastval = nflux[i]
        break
    # Use first and last non-NaN value to
    # fill the nflux array
    nflux[firsts] = firstval
    nflux[lasts] = lastval
  return nflux, wlprime