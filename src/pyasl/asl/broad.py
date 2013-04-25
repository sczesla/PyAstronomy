import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf

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
  
  # To preserve the position of spectral lines, the broadening function
  # must be centered at N//2 - (1-N%2) = N//2 + N%2 - 1
  nwl = (np.arange(len(wvl), dtype=np.int) - sum(divmod(len(wvl),2)) + 1) * dwls[0]
  gf = fuf.GaussFit1d()
  gf["A"] = 1.0
  gf["sig"] = sigma
  e = gf.evaluate(nwl)
  # This step ensured that the 
  e /= np.sum(e)
  
  if edgeHandling == "firstlast":
    nf = len(flux)
    flux = np.concatenate((np.ones(nf)*flux[0], flux, np.ones(nf)*flux[-1]))
    result = np.convolve(flux, e, mode="same")[nf:-nf]
  elif edgeHandling is None:
    result = np.convolve(flux, e, mode="same")
  else:
    raise(PE.PyAValError("Invalid value for `edgeHandling`: " + str(edgeHandling), \
                         where="instrBroadGaussFast", \
                         solution="Choose either 'firstlast' or None"))
  if not fullout:
    return result
  else:
    return (result, fwhm)

  