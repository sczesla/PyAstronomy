# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np
from PyAstronomy.pyasl import _ic

def planck(T, lam=None, nu=None):
  """
    Evaluate Planck's radiation law.
    
    Depending on whether wavelength or frequency is specified as input, the
    function evaluates:
    
    .. math::
    
        B_{\\nu} = \\frac{2\\pi h \\nu^3}{c^2} \\frac{1}{e^{\\frac{h\\nu}{kT}} - 1}
    
    or
    
    .. math::
    
        B_{\\lambda} = \\frac{2\\pi h c^2}{\\lambda^5} \\frac{1}{e^{\\frac{h c}{\\lambda kT}} - 1} \\; .
    
    If lambda is given (in meters), the output units are W/(m^2 m). To convert into erg/(cm^2 A),
    the output has to be multiplied by a factor of 1e-7.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
    lam : float or array, optional
        Wavelength in meters.
    nu : float or array, optional
        Frequency in Hz.
    
    Returns
    -------
    Spectral radiance : float or array
        Depending on whether `lam` or `nu` were specified, returns the
        spectral radiance per area and wavelength or frequency. The unit (SI)
        will be W/(m^2 m) if `lam` was given and W/(m^2 Hz) if `nu` was
        specified.
  """
  if _ic.check["quantities"]:
    from PyAstronomy import constants
    c = constants.PyAConstants(unitSystem="SI")
  else:
    raise(PE.PyARequiredImport("You need to install 'quantities' to use this function.\n 'quantities' can be obtained from 'http://pypi.python.org/pypi/quantities'."))
    return None
  if (lam is not None) and (nu is not None):
    raise(PE.PyAParameterConflict("Specify either 'lam' OR 'nu', but not both."))
  if (lam is None) and (nu is None):
    raise(PE.PyAParameterConflict("Specify either 'lam' OR 'nu'."))
  # Parameters have been specified properly
  if lam is not None:
    result = 2.*np.pi*c.h*(c.c**2)/(lam**5)
    result /= (np.exp(c.h*c.c/(lam*c.k*T)) - 1.0)
    return result
  elif nu is not None:
    result = 2.*np.pi*c.h*(nu**3)/(c.c**2)
    result /= (np.exp(c.h*nu/(c.k*T)) - 1.0)
    return result