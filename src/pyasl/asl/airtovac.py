# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic

def airtovac(wave):
  """
    Convert air wavelengths to vacuum wavelengths 
    
    Parameters
    ----------
    wave : float, array
        The wavelength in air [Angstrom]
   
    Returns
    -------
    Wavelength : array
        Wavelength in vacuum [Angstrom]
   
    Notes
    -----

    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:
    
    NAME:
          AIRTOVAC
    PURPOSE:
          Convert air wavelengths to vacuum wavelengths 
    EXPLANATION:
          Wavelengths are corrected for the index of refraction of air under 
          standard conditions.  Wavelength values below 2000 A will not be 
          altered.  Uses the IAU standard for conversion given in Morton 
          (1991 Ap.J. Suppl. 77, 119)
    
    CALLING SEQUENCE:
          AIRTOVAC, WAVE
    
    INPUT/OUTPUT:
          WAVE - Wavelength in Angstroms, scalar or vector
                 WAVE should be input as air wavelength(s), it will be
                 returned as vacuum wavelength(s).  WAVE is always converted to
                 double precision upon return.
    
    EXAMPLE:
          If the air wavelength is  W = 6056.125 (a Krypton line), then 
          AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
    
    METHOD:
          See Morton (Ap. J. Suppl. 77, 119) for the formula used
    
    REVISION HISTORY
          Written W. Landsman                   November 1991
          Converted to IDL V5.0   W. Landsman   September 1997
  """

  sigma2 = (1.e4/wave)**2.                      #Convert to wavenumber squared

  # Compute conversion factor

  fact = 1. + 6.4328e-5 + 2.94981e-2/(146.-sigma2) + 2.5540e-4/(41.-sigma2)

  fact = fact*(wave >= 2000.) + 1.*(wave < 2000.)

  wave = wave*fact              #Convert Wavelength

  return wave



def vactoair(wave):
  """
    Convert vacuum wavelengths to air wavelengths

    Parameters
    ----------
    wave : float, array
        The wavelength in vacuum [Angstrom]
   
    Returns
    -------
    Wavelength : array,
        Wavelength in air [Angstrom]
   
    Notes
    -----

    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:

    NAME:
          VACTOAIR
    PURPOSE:
          Convert vacuum wavelengths to air wavelengths
    EXPLANATION:
          Corrects for the index of refraction of air under standard conditions.  
          Wavelength values below 2000 A will not be altered.  Accurate to 
          about 0.005 A 

    CALLING SEQUENCE:
          VACTOAIR, WAVE

    INPUT/OUTPUT:
          WAVE - Wavelength in Angstroms, scalar or vector
                  WAVE should be input as vacuum wavelength(s), it will be
                  returned as air wavelength(s).  WAVE is always converted to
                  double precision

    EXAMPLE:
          If the vacuum wavelength is  W = 2000, then 

          IDL> VACTOAIR, W 

          yields an air wavelength of W = 1999.353 Angstroms

    METHOD:
          An approximation to the 4th power of inverse wavenumber is used
          See IUE Image Processing Manual   Page 6-15.

    REVISION HISTORY
          Written, D. Lindler 1982 
          Documentation W. Landsman  Feb. 1989
          Converted to IDL V5.0   W. Landsman   September 1997
  """

  wave2 = wave**2.
  fact = 1. + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2**2.)
  fact = fact * ( wave >= 2000. ) + 1.*( wave < 2000. )

  # Convert wavelengths
  
  wave = wave/fact

  return wave



def specAirVacConvert(wvl, flux, direction):
  """
    Transform spectrum between air and vacuum wavelengths.
    
    In contrast to :py:func:`airtovac` and :py:func:`vactoair`,
    which only transform the wavelength axis, this function
    converts the wavelength axis and interpolates the shifted
    spectrum to retain the input wavelength axis.
    
    Due to the shift, a fraction of the spectrum is lost on one
    side, while there is no valid information on the other.
    Therefore, the function returns both a new flux array obtained
    by linear interpolation and an array of valid indices, i.e.,
    those indices in the new spectrum, which are not NaN.
    
    .. warning:: Linear interpolation of the spectrum can affect its
                 noise properties. This function
                 does *nothing* to prevent or correct for that effect.
    
    Parameters
    ----------
    wvl : array
        Input wavelength axis.
    flux : array
        Input flux array.
    direction : string, {"airtovac","vactoair"}
        If the direction is "vactoair", the input spectrum is
        assumed to be obtained in vacuum conditions. The opposite
        is true otherwise.
    
    Returns
    -------
    New flux : array
        The shifted flux array, defined on the input wavelength axis.
    Valid indices : array
        The array of valid, i.e., not NaN, indices.
  """
  if not _ic.check["scipy"]:
    raise(PE.PyARequiredImport("This function requires scipy.", \
                               solution="Install scipy."))
  import scipy.interpolate as sci
  import numpy as np
  
  if direction == "airtovac":
    newwvl = airtovac(wvl)
  elif direction == "vactoair":
    newwvl = vactoair(wvl)
  else:
    raise(PE.PyAValError("The direction '" + str(direction) + "' is not supported.", \
                         solution="Choose either 'airtovac' or 'vactoair'."))
  
  # Use spectrum with converted wavelength axis for interpolation 
  shiftSpec = sci.interp1d(newwvl, flux, bounds_error=False, fill_value=np.nan)
  # Obtain flux on 'old' wvl axis by interpolation
  sflux = shiftSpec(wvl)
  # Determine 'valid' indices
  vind = np.where(np.bitwise_not(np.isnan(sflux)))[0]
  
  return sflux, vind
  
  