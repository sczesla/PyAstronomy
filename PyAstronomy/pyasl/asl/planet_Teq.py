from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np
from PyAstronomy.pyasl import _ic


def equiTempPlanet(Ts, R, a, albedo=0.0, e=0.0, solarUnits=False):
  """
    Calculate the planetary equilibrium temperature.
    
    The equilibrium temperature is defined, e.g., in
    Laughlin et al. 2011 (ApJ, 729), Eq. 1. Albedo
    was added "by hand".
    
    Note that units have to be consistent for stellar
    radius and large semi-major axis. The `solarUnits`
    flag can be used to use solar units for these
    quantities.
    
    Parameters
    ----------
    Ts : float
        Stellar effective temperature [K].
    R : float
        Radius of the star
    a : float
        Large semi-major axis
    albedo : float
        Albedo of the planetary atmosphere.
        Default is 0.0.
    e : float, optional
        Eccentricity of orbit (default is 0.0).
    solarUnits : boolean, optional
        If this flag is set True, stellar radius and
        semi-major axis are taken to be in solar radii
        and AU.
    
    Returns
    -------
    Equilibrium temperature : float
        The equilibrium temperature of the planet.
  """
  if _ic.check["quantities"]:
    from PyAstronomy import constants as c
  else:
    raise(PE.PyARequiredImport("You need to install 'quantities' to use this function."))
    return None
  if solarUnits:
    a = a * c.AU
    R = R * c.RSun
  if a < R:
    raise(PE.PyAValError("The semi-major axis cannot be smaller than the radius of the star."))
  teq = np.power(1. - albedo, 1./4.) * np.sqrt(R / (2.0 * a)) * Ts \
        / np.power(1.0 - e**2, 1./8.)
  return teq
  