import numpy as np
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE

def solarDirectFluxMeinel(z, height=0.0, I0=1353, c=0.357, s=0.678, a=0.14):
    """
    Direct (unscattered) solar flux in the Earth's atmosphere according to Meinel 1976.
    
    The height-dependent flux of direct solar radiation was given by
    Meinel & Meinel `` Applied Solar Energy'' (Addison Wesley Publishing Co., 1976)
    on page 46 as
    
    .. math::
    
        I_0 \\left( (1-ah) \exp(-c \\times A(z)^s) + ah  \\right)
    
    where :math:`A(z)` is the airmass (evaluated via :func:`airmassPP` using
    the plane parallel approximation)
    
    Parameters
    ----------
    z : float or array
        The zenith angle [deg]
    height : float
        The height above ground [km]
    I0 : float, optional
        Solar constant (W/m**2)
    c : float, optional
        Coefficient
    s : float, optional
        Coefficient
    a : float, optional
        Coefficient
    
    Returns
    -------
    Direct flux : float or array
        Direct solar flux in W/m**2 facing the Sun
     
    """
    f1 = np.exp(-c)
    am = pyasl.airmassPP(z)
    if height == 0.0:
        return I0 * f1**(am**s)
    else:    
        ah = a*height
        if ah >= 1:
            raise(PE.PyAValError("Height is too large for the calculation: (1-ah) is negative."))
        return I0*((1-ah)*f1**(am**s) + ah)