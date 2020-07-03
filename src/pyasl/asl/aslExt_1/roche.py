import numpy as np
from PyAstronomy import constants as PC
from PyAstronomy.pyaC import pyaErrors as PE

def rochepot_dl(x, y, z, q):
    """
    Dimensionless Roche potential (:math:`\\Phi_n`, synchronous rotation)
    
    More massive component (:math:`m_1`) is centered at (x,y,z) = (0,0,0). Less massive
    component (:math:`m_2`) is at (1,0,0). The unit of length is the distance between
    the objects. Both objects are in the x,y plane (x-axis along the connecting line and
    z perpendicular to the orbital plane). 
    
    Parameters
    ----------
    x, y, z : float or array
        Location(s) at which to calculate the potential. Unit
        of length is the distance between the masses m1 and m2.
    q : float
        Mass ratio (0 <= m2/m1 <= 1)
    
    Returns
    -------
    Potential : float or array
        The potential at the specified location(s)
    """
    if (q < 0) or (q > 1):
        raise(PE.PyAValError("Invalid mass ratio q=m2/m1 of '" + str(q) + "'.", \
                             where="rochepot_dl", \
                             solution="Check masses. Perhaps exchange m1 and m2 so that m1 >= m2."))
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x-1)**2 + y**2 + z**2)
    p = 2/((1+q)*r1) + 2*q/((1+q)*r2) + (x - q/(1+q))**2 + y**2
    return p


def rochepot(x, y, z, m1, m2, a):
    """
    Roche potential (:math:`\\Phi`, synchronous rotation)
    
    More massive component (:math:`m_1`) is at (x,y,z) = (0,0,0). Less massive
    component (:math:`m_2`) is at (a,0,0).
    
    Parameters
    ----------
    x, y, z : float or array
        Position [m]
    a : float
        Separation of masses [m]
    m1, m2 : float
        Masses of objects [kg] (m1 >= m2)
    
    Returns
    -------
    Roche potential : float or array
        At specified locations [m**2/s**2]
    """
    # Gravitational constant (m**3/(kg*s**2))
    G = 6.67430e-11
    M = m1 + m2
    return -(G*M)/(2*a) * rochepot_dl(x/a, y/a, z/a, m2/m1)
    
    
