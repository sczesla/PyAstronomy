
def absMagToPower(am, absMagSun=4.75, absLumSun=3.846e33):
  """
    Convert absolute magnitude to power scale
    
    The default values for the absolute magnitude and luminosity
    of the Sun are adopted from Harmanec and Prsa 2011
    (2011PASP..123..976H).
    
    Parameters
    ----------
    am : float
        Absolute magnitude.
    absMagSun : float, optional
        Absolute magnitude of the Sun.
    absLumSun : float, optional
        Absolute luminosity of the Sun.
        The default is given in units of erg/s.
        
    Returns
    -------
    Power : float
        Total emitted power. Same units as `absLumSun`;
        the default corresponds to erg/s.
  """  
  power = 10.0**((am-absMagSun)/(-2.5)) * absLumSun
  return power


def absModuleToDist(magApp, magAbs):
  """
    Convert apparent and absolute magnitude into distance.
    
    Parameters
    ----------
    magApp : float
        Apparent magnitude of object.
    magAbs : float
        Absolute magnitude of object.
    
    Returns
    -------
    Distance : float
        The distance resulting from the difference in
        apparent and absolute magnitude [pc].
  """
  d = 10.0**(-(magAbs - magApp)/5.0 + 1.0)
  return d