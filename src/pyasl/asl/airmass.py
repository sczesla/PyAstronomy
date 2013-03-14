import numpy as np
from PyAstronomy import pyaC as PC
from PyAstronomy.pyaC import pyaErrors as PE

def airmassPP(zangle):
  """
    Calculate airmass for plane parallel atmosphere.
    
    Parameters
    ----------
    zangle : float or array
        The zenith angle in degrees.
    
    Returns
    -------
    Airmass : float or array
        The airmass assuming a plane parallel
        atmosphere.
  """
  return 1.0/np.cos(PC.degtorad(zangle))
  


def airmassSpherical(zangle, obsAltitude, rearth=6371.0, yatm=10.0):
  """
    Calculate the airmass for a given zenith angle and observer altitude.
    
    This routine uses a geometric formula for a homogeneous, spherical
    atmosphere with an elevated observer.
    
    .. note:: In this model, the airmass is *not* necessarily one
              toward the zenith.
    
    Parameters
    ----------
    zangle : float
        Zenith angle of an object in deg.
    obsAltitude :  float
        Elevation of the observer in meter.
    rearth : float, optional
        Earth's radius in km.
    yatm : float, optional
        Height of the atmosphere in km.

    Returns
    -------
    Airmass : float
        The airmass.
  """

  # Convert observer's altitude to km to have
  # consistent units
  obsAltitude = obsAltitude/1000.0
  r = (rearth/yatm)
  y = (obsAltitude/yatm)

  # Find maximal zenith angle 
  zmax = 180.0 - PC.radtodeg(np.arcsin(rearth/(rearth + obsAltitude)))
  if zangle > zmax:
    raise(PE.PyAValError("Zenith angle is too large. The maximum allowed angle is " + str(zmax) + " deg.", \
          solution="Use an angle within the limit. Check observer's altitude."))
  
  # Convert into deg
  zangle = PC.degtorad(zangle)

  return np.sqrt( ( r + y )**2 * np.cos(zangle)**2 + 2.*r*(1.-y) - y**2 + 1.0 ) - \
         (r+y)*np.cos(zangle)
         

