# -*- coding: utf-8 -*-
from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic
import numpy as np

def airtovac(wave, depWarn=True):
  """
    Convert air wavelengths to vacuum wavelengths 
    
    .. warning::
       The conversion implemented here is based on the older formulae
       given by Edlen 1953. Furthermore, it seems that wave numbers in air
       are used, where vacuum wave numbers should be used, which, however,
       produces only a second-order deviation. Consider using
       :py:func:`airtovac2` instead.
    
    Parameters
    ----------
    wave : float, array
        The wavelength in air [Angstrom]
    depWarn : boolean, optional
        If True (default), a deprecation warning will be
        given.
   
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

  if depWarn:
    PE.warn(PE.PyADeprecationError("Note: airtovac is outdated; see documentation for reasons.", \
                                   solution="Consider using 'airtovac2'."))

  sigma2 = (1.e4/wave)**2.                      #Convert to wavenumber squared

  # Compute conversion factor

  fact = 1. + 6.4328e-5 + 2.94981e-2/(146.-sigma2) + 2.5540e-4/(41.-sigma2)

  fact = fact*(wave >= 2000.) + 1.*(wave < 2000.)

  wave = wave*fact              #Convert Wavelength

  return wave



def vactoair(wave, depWarn=True):
  """
    Convert vacuum wavelengths to air wavelengths

    .. warning::
       The conversion implemented here is based on the older formulae
       given by Edlen 1953. Furthermore, it seems that wave numbers in air
       are used, where vacuum wave numbers should be used, which, however,
       produces only a second-order deviation. Consider using
       :py:func:`vactoair2` instead.

    Parameters
    ----------
    wave : float, array
        The wavelength in vacuum [Angstrom]
   
    Returns
    -------
    Wavelength : array,
        Wavelength in air [Angstrom]
    depWarn : boolean, optional
        If True (default), a deprecation warning will be
        given.
   
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

  if depWarn:
    PE.warn(PE.PyADeprecationError("Note: vactoair is outdated; see documentation for reasons.", \
                                   solution="Consider using 'vactoair2'."))

  wave2 = wave**2.
  fact = 1. + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2**2.)
  fact = fact * ( wave >= 2000. ) + 1.*( wave < 2000. )

  # Convert wavelengths
  
  wave = wave/fact

  return wave



def specAirVacConvert(wvl, flux, direction, dontUse2=False):
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
    
    .. warning::
       Linear interpolation of the spectrum can affect its
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
    dontUse2 : boolean, optional
        If True, the now deprecated function airtovac and vactoair
        will be used. Otherwise (default), airtovac2 and vactoair2
        will be used with default settings.
    
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
  
  if direction == "airtovac":
    if dontUse2:
      newwvl = airtovac(wvl)
    else:
      newwvl = airtovac2(wvl)
  elif direction == "vactoair":
    if dontUse2:
      newwvl = vactoair(wvl)
    else:
      newwvl = vactoair2(wvl)
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
  
  


def _iocon(inp):
  """
    Convert input into np float array and create function to convert back.
    
    Parameters
    ----------
    inp : some numeric type
        Input
    
    Returns
    -------
    output : float array
        The input converted into a numpy float array
    recovery function : callable
        Function that converts the output back into a float
        number unless the input was an numpy array. In the latter
        case, it remains an array.
  """
  
  try:
    # Using ravel to avoid generating object with zero length
    out = np.ravel(np.array(inp, dtype=np.float))
  except Exception as e:
    raise(PE.PyAValError("Could not convert input of type: " + str(type(inp)) + " into float array. Error: " + str(e), \
          where="_IOCON", \
          solution="Use an input which can be converted into a numpy float array."))
  
  if isinstance(inp, np.ndarray) or len(out) > 1:
    return out, lambda x:x
  else:
    def recover(x):
      try:
        result = float(x)
      except Exception as e:
        raise(PE.PyAValError("Could not convert result into a float. Error: " + str(e), \
                             where="_iocon", \
                             solution="Use a float or numpy array as input. Are you maybe using a list?"))
      return result
    return out, recover
    


  
class RefractiveIndexAV:
  """
    Refractive index for conversion between air- and vacuum wavelengths.
    
    The conversion is based on the formulae provided either by:
      - Edlen 1953, J. Opt. Soc. Am 43 no. 5
      - Peck and Reeder 1972, J. Opt. Soc. 62 no. 8
      - Ciddor 1996, Applied Optics 35 no. 9
    
    This class provides the routine `refractiveIndex`, `vactoair`,
    and `airtovac` to convert between wavelengths in air and vacuum.
    
    The method `refractiveIndex` points to either of `_edlen53`,
    `_ciddor`, or `_peckReeder` depending on the mode chosen.
    
    To conversion from vacuum to air wavelengths is given by the
    respective formulae provided by the authors of the above papers, who
    all specify the refractive index of "standard air" as a function
    of vacuum wave number. The reverse conversion from air into vacuum
    is done by applying an iterative scheme. Standard air is at 15 degrees
    Celsius, dry, and at a pressure of 750 mmHg = 101325 Pa. Ciddor assumes
    a CO2 concentration of 450 ppm, while Edlen, and Peck and Reeder have
    assumed 300 ppm.
    
    Parameters
    ----------
    mode : string, {edlen53, peckReeder, ciddor}, optional
        The source for the wavelength conversion. The default is
        "ciddor".
      
  """
  
  def _edlen53(self, wvl):
    """
      Refractive index according to Edlen 1953.
      
      Parameters
      ----------
      wvl : float or array
          Vacuum wavelength [A]
      
      Returns
      -------
      n : float or array
          The refractive index.
    """
    if np.any(wvl < 2000.):
      raise(PE.PyAValError("The Edlen 1953 conversion covers only the >2000 A range.",
                           where="RefractiveIndexAV::_edlen53", \
                           solution="Adapt the wavelength range"))
      
    s2 = (1e4/wvl)**2
    n = 1.0 + 6.4328e-5 + (2.94981e-2/(146.0 - s2)) + (2.554e-4/(41. - s2)) 
    return n
  
  def _ciddor(self, wvl, xc=None):
    """
      Refractive index given by Ciddor 1996.
      
      Parameters
      ----------
      wvl : float or array
          Vacuum wavelength [A]
      xc : float, optional
          Concentration of CO2 [ppm]
      
      Returns
      -------
      n : float or array
          Refractive index
    """
    if np.any(wvl < 3000.) or np.any(wvl > 16900.):
      raise(PE.PyAValError("The Ciddor conversion covers only the 3000 - 16900 A range.",
                           where="RefractiveIndexAV::_ciddor", \
                           solution="Adapt the wavelength range"))
    
    # Constants in [(10-6 m)**-2] from Ciddor 1996, Applied Optics, Vol. 35, Issue 9, pp. 1566-1573
    # Appendix A
    k0 = 238.0185
    k1 = 5792105.
    k2 = 57.362
    k3 = 167917.0
    
    s2 = (1e4/wvl)**2
    
    # Eqs. 1 and 2
    nas = (k1/(k0 - s2) + k3/(k2 - s2)) / 1e8 + 1.0
    if xc is not None:
      naxs = (nas - 1.0)*(1.0 + 0.534e-6*(xc - 450.)) + 1.0
      return naxs
    return nas
    
  def _peckReeder(self, wvl):
    """
      Refractive index from Peck and Reeder 1972.
      
      Parameters
      ----------
      wvl : float or array
          Vacuum wavelength [A]
      
      Returns
      -------
      n : float or array
          The refractive index.
    """
    if np.any(wvl < 2300.) or np.any(wvl > 16900.):
      raise(PE.PyAValError("The Peck and Reeder conversion covers only the 2300 - 16900 A range.",
                           where="RefractiveIndexAV::_peckReeder", \
                           solution="Adapt the wavelength range"))
    
    s2 = (1e4/wvl)**2
    
    # Eq. 2
    n = (5791817. / (238.0185 - s2) + 167909./(57.362 - s2)) / 1e8 + 1.0
    return n
  
  def vactoair(self, wvl, **kwargs):
    """
      Convert vacuum wavelength into wavelength in air.
      
      Parameters
      ----------
      wvl : float or array
          Wavelength in vacuum [A].
      kwargs : dict
          Additional keywords eventually needed in the call
          to the refractiveIndex function (depends on the
          specified 'mode').
      
      Returns
      -------
      wvl : float or array
          Wavelength in air.
    """
    wvl, rec = _iocon(wvl)
    result = wvl/self.refractiveIndex(wvl, **kwargs)
    return rec(result)
  
  def airtovac(self, wvl, precision=1e-12, maxiter=30, **kwargs):
    """
      Convert air wavelength into wavelength in vacuum.
      
      An iterative scheme is applied to invert the formula.
      
      Parameters
      ----------
      wvl : float or array
          Wavelength in vacuum [A].
      precision : float, optional
          The target precision beyond which iteration should
          be stopped; default is 1e-12.
      maxiter : int, optional
          Maximum number of iterations used to invert the
          conversion; default is 30.
      kwargs : dict
          Additional keywords eventually needed in the call
          to the refractiveIndex function (depends on the
          specified 'mode').
      
      Returns
      -------
      wvl : float or array
          Wavelength in air.
    """
    # First guess for vacuum wavelength is wavelength in air
    wvl, rec = _iocon(wvl)
    lvac = wvl.copy()
    counter = 0
    while True:
      counter += 1
      dl = lvac/self.refractiveIndex(lvac) - wvl
      lvac -= dl
      # Obtain difference between air wavelength resulting from
      # current estimate of vacuum wavelength 
      dlair = self.vactoair(lvac) - wvl
      if np.abs(np.max(dlair)) < precision:
        return rec(lvac) 
      if counter > maxiter:
        p = np.max(dlair)
        raise(PE.PyAValError("Maximum of iterations reached. Current precision is " + str(p) + " A.", \
                             where="RefractiveIndexAV::airtovac", \
                             solution=["Try to increase 'maxiter'", \
                                       "Try to use subrange of wavelengths", \
                                       "Try to adapt 'precision' to necessary value"]))
      
    return rec(lvac)
  
  def __init__(self, mode="ciddor"):
    
    if mode == "edlen53":
      self.refractiveIndex = self._edlen53
    elif mode == "peckReeder":
      self.refractiveIndex = self._peckReeder
    elif mode == "ciddor":
      self.refractiveIndex = self._ciddor
    else:
      raise(PE.PyAValError("Unknown mode: '" + str(mode) + "'.", \
                           where="RefractiveIndexAV", \
                           solution="Use either of 'edlen53', 'peckReeder', or 'ciddor'."))


def vactoair2(wvl, mode="ciddor", **kwargs):
  """
    Converts wavelength in vacuum into wavelength in air. 
    
    This function is a wrapper around the `RefractiveIndexAV` class.

    Parameters
    ----------
    wvl :float or array
        Wavelength in vacuum [A].
    mode : string, {edlen53, peckReeder, ciddor}, optional
        The source for the wavelength conversion. The default is
        "ciddor".
    kwargs : dict
        Additional parameters handed to the vactoair method
        of RefractiveIndexAV.
    
    Returns
    -------
    wvl : float or array
        Wavelength in air [A].
  """
  c = RefractiveIndexAV(mode=mode)
  return c.vactoair(wvl, **kwargs)

def airtovac2(wvl, mode="ciddor", **kwargs):
  """
    Converts wavelength in air into wavelength in vacuum. 

    This function is a wrapper around the `RefractiveIndexAV` class.

    Parameters
    ----------
    wvl :float or array
        Wavelength in vacuum [A].
    mode : string, {edlen53, peckReeder, ciddor}, optional
        The source for the wavelength conversion. The default is
        "ciddor".
    kwargs : dict
        Additional parameters handed to the airtovac method
        of RefractiveIndexAV.
    
    Returns
    -------
    wvl : float or array
        Wavelength in vacuum [A].
  """
  c = RefractiveIndexAV(mode=mode)
  return c.airtovac(wvl, **kwargs)