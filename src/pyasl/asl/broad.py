import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf


def broadGaussFast(x, y, sigma, edgeHandling=None):
  """
    Apply Gaussian broadening. 
    
    This function broadens the given data using a Gaussian
    kernel.
    
    Parameters
    ----------
    x, y : arrays
        The abscissa and ordinate of the data.
    sigma : float
        The width (i.e., standard deviation) of the Gaussian
        profile used in the convolution.
    edgeHandling : string, {None, "firstlast"}, optional
        Determines the way edges will be handled. If None,
        nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and
        last value at the start or end. Note that this is
        not necessarily appropriate. The default is None.
    
    Returns
    -------
    Broadened data : array
        The input data convolved with the Gaussian
        kernel.
  """
  # Check whether x-axis is linear
  dxs = x[1:] - x[0:-1]
  if abs(max(dxs) - min(dxs)) > np.mean(dxs)*1e-6:
    raise(PE.PyAValError("The x-axis is not equidistant, which is required.", \
                         where="broadGaussFast"))
  
  # To preserve the position of spectral lines, the broadening function
  # must be centered at N//2 - (1-N%2) = N//2 + N%2 - 1
  nx = (np.arange(len(x), dtype=np.int) - sum(divmod(len(x),2)) + 1) * dxs[0]
  gf = fuf.GaussFit1d()
  gf["A"] = 1.0
  gf["sig"] = sigma
  e = gf.evaluate(nx)
  # This step ensured that the 
  e /= np.sum(e)
  
  if edgeHandling == "firstlast":
    nf = len(y)
    y = np.concatenate((np.ones(nf)*y[0], y, np.ones(nf)*y[-1]))
    result = np.convolve(y, e, mode="same")[nf:-nf]
  elif edgeHandling is None:
    result = np.convolve(y, e, mode="same")
  else:
    raise(PE.PyAValError("Invalid value for `edgeHandling`: " + str(edgeHandling), \
                         where="broadGaussFast", \
                         solution="Choose either 'firstlast' or None"))
  return result
  

def instrBroadGaussFast(wvl, flux, resolution, edgeHandling=None, fullout=False):
  """
    Apply Gaussian instrumental broadening. 
    
    This function broadens a spectrum assuming a Gaussian
    kernel. The width of the kernel is determined by the
    resolution. In particular, the function will determine
    the mean wavelength and set the Full Width at Half
    Maximum (FWHM) of the Gaussian to
    (mean wavelength)/resolution. 
    
    Parameters
    ----------
    wvl : array
        The wavelength
    flux : array
        The spectrum
    resolution : int
        The spectral resolution.
    edgeHandling : string, {None, "firstlast"}, optional
        Determines the way edges will be handled. If None,
        nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and
        last value at the start or end. Note that this is
        not necessarily appropriate. The default is None.
    fullout : boolean, optional
        If True, also the FWHM of the Gaussian will be returned.
    
    Returns
    -------
    Broadened spectrum : array
        The input spectrum convolved with a Gaussian
        kernel.
    FWHM : float, optional
        The Full Width at Half Maximum (FWHM) of the
        used Gaussian kernel.
  """
  # Check whether wvl axis is linear
  dwls = wvl[1:] - wvl[0:-1]
  if abs(max(dwls) - min(dwls)) > np.mean(dwls)*1e-6:
    raise(PE.PyAValError("The wavelength axis is not equidistant, which is required.", \
                         where="instrBroadGaussFast"))
  meanWvl = np.mean(wvl)
  fwhm = 1.0/float(resolution) * meanWvl
  sigma = fwhm/(2.0*np.sqrt(2.*np.log(2.)))
  
  result = broadGaussFast(wvl, flux, sigma, edgeHandling=edgeHandling)
  
  if not fullout:
    return result
  else:
    return (result, fwhm)
  

  