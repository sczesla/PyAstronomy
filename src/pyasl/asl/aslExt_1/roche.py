from __future__ import print_function, division
import numpy as np
import inspect
from PyAstronomy import constants as PC
from PyAstronomy.pyaC import pyaErrors as PE
from scipy import optimize as sco

def _r1r2_dl(x, y, z):
    """
    Calculate r1 and r2 for dimensionless Roche potential
    """
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x-1)**2 + y**2 + z**2)
    return r1, r2

def _checkq(q):
    """ Check validity of q """
    if (q < 0) or (q > 1):
        raise(PE.PyAValError("Invalid mass ratio q=m2/m1 of '" + str(q) + "'.", \
                             where=inspect.stack()[1][3], \
                             solution="Check masses. Perhaps exchange m1 and m2 so that m1 >= m2."))

def _checkm(m):
    if not m in [1,2]:
        raise(PE.PyAValError("Invalid m(ass). Current value is '" + str(m) + "'", \
                             where=inspect.stack()[1][3], \
                             solution="Use 1 or 2."))

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
    _checkq(q)
    r1, r2 = _r1r2_dl(x, y, z)
    p = 2/((1+q)*r1) + 2*q/((1+q)*r2) + (x - q/(1+q))**2 + y**2
    return p

def ddx_rochepot_dl(x, q, y=0, z=0):
    """
    Derive of dimensionless Roche potential along x-axis
    
    Parameters
    ----------
    x, y, z : float or array
        Position (default for y and z is zero)
    q : float
        Mass ratio m2/m1
    
    Returns
    -------
    Derivative: float or array
        d/dx of dimensionless Roche potential
    """
    _checkq(q)
    r1, r2 = _r1r2_dl(x, y, z)
    ddx = 2/(1+q) * (-x/r1**3) + 2*q/(1+q) * (-(x-1)/r2**3) + 2*(x - q/(1+q))
    return ddx

def ddz_rochepot(z, q, x=1, y=0):
    """
    Derive of dimensionless Roche potential along z-axis
    
    Parameters
    ----------
    x, y, z : float or array
        Position (default for x is one and zero for z)
    q : float
        Mass ratio m2/m1
    
    Returns
    -------
    Derivative: float or array
        d/dz of dimensionless Roche potential
    """
    _checkq(q)
    r1, r2 = _r1r2_dl(x, y, z)
    ddz = 2/(1+q) * (-z/r1**3) + 2*q/(1+q) * (-z/r2**3) 
    return ddz

def roche_lobe_radius_eggleton(q, m):
    """
    Approximate dimensionless Roche lobe radius according to Eggelton 1983 (ApJ 268, 368).
    
    .. math::
    
        r_L \\approx \frac{0.49\,q^{2/3}}{0.6\,q^{2/3} + \\ln(1+q^{1/3})}
    
    Parameters
    ----------
    q : float or array
        Mass ratio (m2/m1)
    m : int, {1,2}
        Calculate radius for mass 1 or 2
    """
    _checkq(q)
    _checkm(m)
    if m == 1:
        # Hilditch p161
        q = 1/q
    q23 = q**(2./3)
    return 0.49*q23 / (0.6*q23 + np.log(1+q**(1./3)))

def _get_lagrange_123(q, a, b, getdlrp=True):
    """
    Find root of derivative along x-axis
    
    Parameters
    ----------
    q : float
        Mass ratio (m2/m1)
    a, b : float
        Limits for root search
    getdlrp : boolean, optional
        If True (default), also return potential at that point
    """
    r = sco.brentq(ddx_rochepot_dl, a, b, args=(q))
    if getdlrp:
        dlrp = rochepot_dl(r, 0, 0, q)
        return r, dlrp
    return r

def get_lagrange_1(q, getdlrp=True, eps=1e-4):
    """
    Get location of first Lagrange point
    
    Parameters
    ----------
    q : float
        Mass ratio
    getdlrp : boolean, optional
        If True (default), also the dimensionless Roche potential
        at the first Lagrange point is returned.
    eps : float, optional
        The potential diverges at x=0 and x=1. The search for a root of
        the derivative is performed in the interval [eps, 1-eps].
    
    Returns
    -------
    xL1 : float
        The location of the first Lagrange point (xL1, 0, 0).
    Potential : float, optional
        The dimensionless potential at the first Lagrange point (default is True).
    """
    _checkq(q)
    return _get_lagrange_123(q, eps, 1-eps, getdlrp)

def get_lagrange_2(q, getdlrp=True, eps=1e-4):
    """
    Get location of second Lagrange point
    
    Parameters
    ----------
    q : float
        Mass ratio
    getdlrp : boolean, optional
        If True (default), also the dimensionless Roche potential
        at the first Lagrange point is returned.
    eps : float, optional
        The potential diverges at x=0 and x=1. The search for a root of
        the derivative is performed in the interval [1+eps, 1+10*rL], where rL
        is the estimated (effective) Roche lobe radius according to Eggleton 1983.
    
    Returns
    -------
    xL2 : float
        The location of the second Lagrange point (xL2, 0, 0).
    Potential : float, optional
        The dimensionless potential at the second Lagrange point (default is True).
    """
    _checkq(q)
    rl = roche_lobe_radius_eggleton(q,2)
    return _get_lagrange_123(q, 1+eps, 1+10*rl, getdlrp)

def get_lagrange_3(q, getdlrp=True, eps=1e-4):
    """
    Get location of third Lagrange point
    
    Parameters
    ----------
    q : float
        Mass ratio
    getdlrp : boolean, optional
        If True (default), also the dimensionless Roche potential
        at the third Lagrange point is returned.
    eps : float, optional
        The potential diverges at x=0. The search for a root of
        the derivative is performed in the interval [-10, -eps].
    
    Returns
    -------
    xL3 : float
        The location of the first Lagrange point (xL3, 0, 0).
    Potential : float, optional
        The dimensionless potential at the third Lagrange point (default is True).
    """
    _checkq(q)
    return _get_lagrange_123(q, -10, -eps, getdlrp)

def get_lagrange_4():
    """
    Get location of forth Lagrange point
    
    Orbital angular momentum is supposed to point into +z direction.
    L4 is in direction of motion of m2.
    
    Returns
    -------
    x,y,z : float
        The location of the fourth Lagrange point (xL3, 0, 0).
    """
    return (0.5, np.sin(np.radians(60)), 0)

def get_lagrange_5():
    """
    Get location of fifth Lagrange point
    
    Orbital angular momentum is supposed to point into +z direction.
    L5 is behind motion of m2.
    
    Returns
    -------
    x,y,z : float
        The location of the fifth Lagrange point (xL3, 0, 0).
    """
    return (0.5, -np.sin(np.radians(60)), 0)

def rl_vol_tabint(q, m=2, n=50, eps=1e-4, fullout=True):
    """
    Calculate (dimensionless) Roche lobe volume
    
    Parameters
    ----------
    q : float
        Mass ratio (m2/m1)
    m : int, {1,2}
        Whether to calculate volume for m1 or m2
    n : int, optional
        Number of point to be used in numerical integration. Default is 50.
        Increase to improve accuracy.
    eps : float, optional
        Parameter used to avoid singularities of Roche potential. Default is 1e-4.
    fullout : boolean, optional
        If True, provides more output than volume.
    
    Returns
    -------
    Volume : float
        Dimensionless volume of Roche lobe
    Effective radius : float, optional
        Radius of a sphere with the same radius (dimensionless)
    xx, rx : arrays, optional
        Points on the x-axis as well as extent of body in the z (and equivalently y) direction.
        Volume is calculated as V = pi*int_xx rx**2.
    """
    _checkq(q)
    _checkm(m)
    l1, l1pot = get_lagrange_1(q, eps=eps)
    l2, _ = get_lagrange_2(q, eps=eps)
    l3, _ = get_lagrange_3(q, eps=eps)
    f = lambda x: rochepot_dl(x,0,0,q) - l1pot
    if m == 2:
        # Outer limit (beyond m2) of Roche lobe
        p1 = l1
        p2 = sco.brentq(f, 1+eps, l2)
    elif m == 1:
        # Outer limit (behind m1) of Roche lobe
        p1 = sco.brentq(f, l3, -eps)
        p2 = l1
    # Prepare numerical integration
    xx = np.linspace(p1, p2, n)
    fx = np.zeros_like(xx)
    # Leaves outer point at zero expansion
    for i in range(1,len(xx)-1):
        f = lambda z: rochepot_dl(xx[i],0,z,q) - l1pot
        fx[i] = sco.brentq(f, 0, 1)
    # Get volume of rotational body
    vol = np.pi*np.trapz(fx**2, xx)
    if fullout:
        reff = (vol*3./4/np.pi)**(1./3)
        return vol, reff, xx, fx

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
    
    
