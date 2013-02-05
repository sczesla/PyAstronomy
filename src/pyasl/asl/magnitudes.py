
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