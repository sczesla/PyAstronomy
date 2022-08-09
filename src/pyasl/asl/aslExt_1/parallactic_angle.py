import numpy as np

def parallacticAngle(t, gl, dec):
    """
    Calculate the parallactic angle
    
    The parallactic angle is the position angle of the zenith
    at the target.
    
    Notes
    -----
    Useful links
    
    `Parallactic Angle Observations (Keck) <https://www2.keck.hawaii.edu/inst/common/parallactic.html>`_
    
    `What Wiki has to say <https://en.wikipedia.org/wiki/Parallactic_angle>`_
    
    
    Parameters
    ----------
    t : float or array
        Hour angle [deg]
    gl : float
        Latitude of the observer [deg]
    dec : float
        Declination of the target [deg]
    
    Returns
    -------
    q : float or array
        The parallactic angle [deg]
        
    """
    t = np.deg2rad(t)
    gl = np.deg2rad(gl)
    dec = np.deg2rad(dec)
    
    up = np.sin(t)
    down = np.tan(gl)*np.cos(dec) - np.sin(dec)*np.cos(t)
    r = np.arctan2(up,down)
    return np.rad2deg(r)
