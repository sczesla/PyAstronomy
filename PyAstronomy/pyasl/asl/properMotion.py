from __future__ import division
import numpy as np

def applyProperMotion(ra, dec, pmra, pmdec, dt, fixes=1):
    """
    Apply proper motion to coordinates
    
    Parameters
    ----------
    ra, dec : float or array
        Coordinates at reference epoch
    pmra, pmdec : float or array
        Proper motion in mas/yr.
    dt : float
        Elapsed time [yr]
    fixes : int, optional
        Number of steps in which the proper motion is applied.
        Default is one.
    """
    
    # mas/yr -> deg/yr
    pmra = pmra / 3.6e6
    pmdec = pmdec / 3.6e6
    
    dtt = dt / fixes
    for i in range(fixes):
        decnew = dec + pmdec*dtt
        ra = ra + pmra*dtt / np.cos(np.deg2rad((dec+decnew)/2))
        dec = decnew
        
    return ra, dec
    