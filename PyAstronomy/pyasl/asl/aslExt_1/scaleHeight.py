from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import constants as PC

def atmosphericScaleHeight(T, mu, g):
    """
    Atmospheric scale height
    
    Parameters
    ----------
    T : float
        Temperature in K
    mu : float
        Mean olecular weight
    g : float
        Gravitational acceleration [m/s**2]
    
    Returns
    -------
    Scale height : float
        Atmospheric scale height in km
    """
    if T <= 0:
        raise(PE.PyAValError("Temperature must be positive.", \
                             where="atmosphericScaleHeight"))
    if mu <= 0:
        raise(PE.PyAValError("Mean molecular weight must be positive.", \
                             where="atmosphericScaleHeight"))        
    if g <= 0:
        raise(PE.PyAValError("Gravitational acceleration must be positive.", \
                             where="atmosphericScaleHeight"))
    pc = PC.PyAConstants()
    pc.setSystem("SI")
    H = pc.k * T / (mu*pc.u*g)
    H /= 1e3
    return H


def atmosphericScaleHeight_MR(T, mu, Mp, Rp, ref):
    """
    Atmospheric scale height
    
    Parameters
    ----------
    T : float
        Temperature in K
    mu : float
        Mean olecular weight
    Mp : float
        Mass of planet [wrt Jupiter or Earth (see `ref`)]
    Rp : float
        Radius of planet [wrt Jupiter or Earth (see `ref`)]
    ref : string, {J, E}
        Determines whether Mp and Rp are considered in units of the Earth (E)
        or Jupiter (J)
    
    Returns
    -------
    Scale height : float
        Atmospheric scale height in km
    """
    if Mp <= 0:
        raise(PE.PyAValError("Planetary mass must be positive.", \
                             where="atmosphericScaleHeight"))
    if Rp <= 0:
        raise(PE.PyAValError("Planetary radius must be positive.", \
                             where="atmosphericScaleHeight"))
    if not ref in ['E', 'J']:
        raise(PE.PyAValError("ref must either be 'E' or 'J'.", \
                             where="atmosphericScaleHeight"))
    
    pc = PC.PyAConstants()
    pc.setSystem("SI")
    if ref == "E":
        r0 = pc.REarth
        m0 = pc.MEarth
    else:
        r0 = pc.RJ
        m0 = pc.MJ
    
    g = pc.G * (Mp*m0) / (Rp*r0)**2
    
    return atmosphericScaleHeight(T, mu, g)
