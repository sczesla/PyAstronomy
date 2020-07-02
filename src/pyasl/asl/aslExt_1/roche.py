import numpy as np
from PyAstronomy import constants as PC

def rochepot_dl(x, y, z, q):
    """
    Dimensionless Roche potential
    
    More massive component (m1) is at (x,y,z) = (0,0,0). Less massive
    component (m2) is at (1,0,0). The unit of length is the distance between
    the objects. Both objects are in the x,y plane. The dimensionless Roche
    potential is given by
    
    .. math::
    
        \\Phi_n &= \\frac{2}{(1+q) r_1} + \\frac{2 q}{(1+q) r_2} + \\left(x - \\frac{q}{1+q}\\right)^2 + y^2 \\\\
        r_1 &= \\sqrt{x^2 + y^2 + z^2} \\\\
        r_2 &= \\sqrt{(x-1)^2 + y^2 + z^2}
    
    Parameters
    ----------
    x, y, z : float or array
        Location(s) at which to calculate the potential. Unit
        of length is the distance between the masses m1 and m2.
    q : float
        Mass ratio (m2/m1 <= 1)
    
    Returns
    -------
    Potential : float or array
        The potential at the specified location(s)
    """
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x-1)**2 + y**2 + z**2)
    p = 2/((1+q)*r1) + 2*q/((1+q)*r2) + (x - q/(1+q))**2 + y**2
    return p


def rochepot(x, y, z, m1, m2, a):
    """
    Roche potential
    
    More massive component (m1) is at (x,y,z) = (0,0,0). Less massive
    component (m2) is at (a,0,0). Both objects are in the x,y plane.
    The Roche potential is given by
    
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
    G = 6.67430
    M = m1 + m2
    return -(G*M)/(2*a) * rochepot_dl(x/a, y/a, z/a, m2/m2)
    
    
