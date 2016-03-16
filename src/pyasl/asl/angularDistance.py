import numpy as np

def getAngDist(ra1, dec1, ra2, dec2):
  """
    Calculate the angular distance between two coordinates.
    
    Parameters
    ----------
    ra1 : float, array
        Right ascension of the first object in degrees.
    dec1 : float, array
        Declination of the first object in degrees.
    ra2 : float, array
        Right ascension of the second object in degrees.
    dec2 : float, array
        Declination of the second object in degrees.

    Returns
    -------
    Angle : float, array
        The angular distance in DEGREES between the first
        and second coordinate in the sky.
         
  """  
  
  delt_lon = (ra1 - ra2)*np.pi/180.
  delt_lat = (dec1 - dec2)*np.pi/180.
  # Haversine formula
  dist = 2.0*np.arcsin( np.sqrt( np.sin(delt_lat/2.0)**2 + \
         np.cos(dec1*np.pi/180.)*np.cos(dec2*np.pi/180.)*np.sin(delt_lon/2.0)**2 ) )  

  return dist/np.pi*180.
