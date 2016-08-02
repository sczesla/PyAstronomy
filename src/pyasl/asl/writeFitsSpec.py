from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic
import os
import six.moves as smo
import six

def write1dFitsSpec(fn, flux, wvl=None, waveParams=None, fluxErr=None, header=None, clobber=False, refFileName=None, refFileExt=0):
  """
    Write a 1d spectrum with equidistant binning to a fits file.
    
    Write a 1d-spectrum to a file. Wavelength axis and header keywords
    are related through the following expression:
    wvl = ((np.arange(N) + 1.0) - CRPIX1) * CDELT1 + CRVAL1,
    where CRPIX1, CDELT1, and CRVAL1 are the
    relevant header keywords.
    
    The function allows to specify an existing fits extension, from
    which the header will be cloned. Alternatively, an arbitrary,
    user-defined header may be given.
    
    Parameters
    ----------
    fn : string
        Filename
    flux : array
        Flux array
    wvl : array, optional
        Wavelength array. Either the wavelength array or the header
        keywords have to be provided (see `waveParams`). 
    waveParams : dict, optional
        Wavelength information required in the fits-header.
        Required keys are CDELT, CRVAL, CRPIX or (CDELT1, CRVAL1, CRPIX1).
    fluxErr : array, optional
        Flux errors. If given, the error will be stored in an additional extension.
    header : dict, optional
        Dictionary with header information to be transfered to the new file.
        Note that the wavelength information will be overwritten by the information
        by the information given to this routine.
        If both, a reference file to clone the header from and the header parameters
        are given, the cloned header from the reference file will be overwritten and
        extended by the keywords specified here.
    refFileName : string, optional
        Clone header keywords from a reference file.
        Note that the wavelength information will be overwritten by the information
        by the information given to this routine.
    refFileExt : int, optional
        Reference-file extension to be used for cloning the header keywords.
        The default is 0.
    
  """
  
  reservedKeyWords = ["BITPIX", "SIMPLE","NAXIS","NAXIS1","CDELT1","CRPIX1","CRVAL1","EXTEND"]
  
  if (not _ic.check["pyfits"]) and (not _ic.check["astropy.io.fits"]):
    raise(PE.PyARequiredImport("Could neither import module 'pyfits' and 'astropy.io.fits'.", \
                               where="write1dFitsSpec", \
                               solution="Install pyfits: http://www.stsci.edu/institute/software_hardware/pyfits"))

  if wvl is None and waveParams is None:
    raise(PE.PyAValError("The wavelength axis is not defined, i.e., neither \'wvl\' nor \'waveParams\' is specified.", \
                         where="write1dFitsSpec", \
                         solution="Provide wavelength array as \'wvl\' or wavelength information as \'waveParams\'."))                               

  if wvl is not None and waveParams is not None:
    raise(PE.PyAValError("You provided the wavelength axis as \'wvl\' AND wavelength information as \'waveParams\'. Don't know which to use.", \
                         where="write1dFitsSpec", \
                         solution="Provide only wavelength array via \'wvl\' or wavelength information via \'waveParams\'."))                               
                         
  if wvl is not None:
    # Check whether the wavelength axis is equidistant
    dwl = wvl[1:] - wvl[0:-1]
    if (np.max(dwl) - np.min(dwl)) / np.max(dwl) > 1e-6:
      raise(PE.PyAValError("Wavelength axis seems not to be equidistant.", \
                           where="write1dFitsSpec", \
                           solution=["Check wavelength array.", \
                                     "Consider passing wavelength information via `waveParams`."]))
                         
  if os.path.isfile(fn) and not clobber:
    raise(PE.PyAFileError("File " + str(fn) + "already exists.", \
                          where="write1dFitsSpec", \
                          solution="File exists; set clobber=True to overwrite."))
    
  if _ic.check["pyfits"]:
    import pyfits
  else:
    import astropy.io.fits as pyfits
    
  
  # Put the flux into the output file.
  # Generate primary HDU
  hdu = pyfits.PrimaryHDU(flux)
  
  if refFileName:
    # If the is a reference file, its header will be used to populate
    # the header of the file to be written.
    if not os.path.isfile(refFileName):
      raise(PE.PyAFileError("File " + str(refFileName) + "does not exist.", \
            where="write1dFitsSpec", \
            solution="Check file name."))
    ff = pyfits.open(refFileName)[refFileExt]
    
    for k in ff.header.keys():
      if k not in reservedKeyWords and k != "COMMENT":  
        try:
          if len(k) > 8: 
            hdu.header["HIERARCH " + k] = ff.header[k]
          else:  
            hdu.header[k] = ff.header[k]
        except:
          PE.warn(PE.PyAValError("Cannot write keyword <" + str(k) + "> with content " + str(ff.header[k])))
  
  if header is not None:
    # A use-defined header was specified
    for k in header.keys():
      hdu.header[k] = header[k]
  
  # Header keywords relevant for wavelength axis
  hk = {}  
  # The only allowed type here
  hk["CTYPE1"] = "Linear"
  # If wavelength array is provided, create header keywords:
  if wvl is not None:           
    hk["CRPIX1"] = 1
    hk["CRVAL1"] = wvl[0]
    hk["CDELT1"] = float(wvl[-1] - wvl[0]) / (len(wvl)-1)
  elif waveParams is not None:
    requiredKeys = ["CRPIX", "CRVAL", "CDELT"]
    # Counts whether all required keywords are provided
    count = 0
    for k in requiredKeys:
      subCount = 0
      for p in six.iterkeys(waveParams):
        if p.upper() == k or p.upper() == k+"1":
          hk[k+"1"] = waveParams[p]
          subCount += 1
      if subCount == 1:
        count += 1
      else:
        # This occurs if, e.g., both "CRVAL" and "CRVAL1" are present 
        raise(PE.PyAValError("The wavelength parameters provided via \'waveParams\' contain ambiguous parameters. " + \
                             "You provided multiple values for " + k + ".", \
                              where="write1dFitsSpec", \
                              solution="Check content of \'waveParams\' so that only one value for "+ k + \
                                       " is provided (including lower/upper case and presence of a trailing \'1\'."))           
    if count < 3:
      raise(PE.PyAFileError("You provided an incomplete set of waveParams.", \
            where="write1dFitsSpec", \
            solution="Required keywords are CRPIX, CRVAL, and CDELT."))
  
  # (Over-)Write wavelength-related header information
  for k in hk.keys():
          hdu.header[k] = hk[k]
  
  if not fluxErr is None:
    # An error on the flux was given
    hdue = pyfits.ImageHDU(fluxErr)
    for k in hk.keys():
      # Add wavelength information to error header
      hdue.header[k] = hk[k]
    hdulist = pyfits.HDUList([hdu, hdue])
  else:  
    hdulist = pyfits.HDUList([hdu])
    
  hdulist.writeto(fn, clobber=clobber)