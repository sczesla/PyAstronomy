import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic
import os

def read1dFitsSpec(fn, hdu=0, fullout=False, CRPIX1=None):
  """
    Read a simple 1d spectrum from fits file.
    
    Reads a 1d-spectrum from a file and constructs the
    associated wavelength axis. To this end, the expression:
    wvl = ((np.arange(N) + 1.0) - CRPIX1) * CDELT1 + CRVAL1
    will be evaluated, where N is the number of bins and
    CRPIX1, CDELT1, and CRVAL1 are header keywords.
    
    Parameters
    ----------
    fn : string
        Filename
    hdu : int, optional
        The number of the HDU to be used. The default
        is 0, i.e., the primary HDU.
    fullout : boolean, optional
        If True, the header keywords used to construct
        the wavelength axis will be returned. The default
        is False.
    CRPIX1 : int, optional
        Can be used to circumvent missing CRPIX1 entry.
    
    Returns
    -------
    wvl : array
        The wavelength array.
    flx : array
        The flux array.
    
    Examples
    --------
    
    ::
      
      from PyAstronomy import pyasl
      wvl, flx = pyasl.read1dFitsSpec("mySpectrum.fits")
  """
  if (not _ic.check["pyfits"]) and (not _ic.check["astropy.io.fits"]):
    raise(PE.PyARequiredImport("Could neither import module 'pyfits' and 'astropy.io.fits'.", \
                               where="read1dFitsSpec", \
                               solution="Install pyfits: http://www.stsci.edu/institute/software_hardware/pyfits"))
  
  if not os.path.isfile(fn):
    raise(PE.PyAFileError("File " + str(fn) + "does not exist.", \
                          where="read1dFitsSpec", \
                          solution="Check file name."))
  
  if _ic.check["pyfits"]:
    import pyfits
  else:
    import astropy.io.fits as pyfits
  
  
  hl = pyfits.open(fn)
  
  naxis = hl[hdu].header["NAXIS"]
  if hl[hdu].header["NAXIS"] != 1:
    raise(PE.PyAValError("There is more than one axis (NAXIS = " + str(naxis) + ").", \
                         where="read1dFitsSpec", \
                         solution="Check file and HDU."))

  hkeys = {"CRVAL1":None, "CRPIX1":CRPIX1, "CDELT1":None}
  for k in hkeys.keys():
    if hkeys[k] is not None:
      continue
    if not k in hl[hdu].header:
      raise(PE.PyAValError("Header does not contain required keyword '" + str(k) + "'", \
                         where="read1dFitsSpec", \
                         solution="Check file and HDU."))
    hkeys[k] = hl[hdu].header[k]
  
  N = int(hl[hdu].header["NAXIS1"])
  
  # Construct wvl axis
  wvl = ((np.arange(N) + 1.0) - hkeys["CRPIX1"]) * hkeys["CDELT1"] + hkeys["CRVAL1"]
  # Get flux
  flx = np.array(hl[hdu].data)
  
  if fullout:
    # Full output
    return wvl, flx, hkeys
  
  return  wvl, flx
  
  