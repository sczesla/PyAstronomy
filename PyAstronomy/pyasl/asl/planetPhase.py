import numpy as np

def lambertPhaseFunction(alpha):
  """
    Calculate phase function for a Lambert sphere.
    
    The phase function of the Lambert sphere is given by:
    
    .. math:: \\Phi(\\alpha) = \\frac{\\sin(\\alpha) + (\\pi - \\alpha)\\cos(\\alpha)}{\\pi} \\; .
    
    Here, :math:`\\alpha` is the phase angle, which is defined as the angle
    between the star and the Earth as seen from the planet. Hence, at phase
    angle zero, the planet is found in opposition. Formally, the phase angle
    can be between 0 and 180 degrees. This function accounts for cases in
    which the given phase angle violates these limits by projecting it back
    into the valid range.
    
    Parameters
    ----------
    alpha : float or array
        The phase angle(s) [deg] at which the phase function
        is to be calculated.
    
    Returns
    -------
    Phase function : float or array
        The values of the phase function at the input
        phase angles.
  """
  a = np.abs(alpha/180.0 * np.pi)
  if not isinstance(a, float):
    a = a - np.floor(a/(2.*np.pi))*(2.*np.pi)
    indi = np.where(a > np.pi)[0]
    a[indi] = 2.*np.pi - a[indi]
    return (np.sin(a) + (np.pi - a)*np.cos(a)) / np.pi
  return float(np.sin(a) + (np.pi - a)*np.cos(a)) / np.pi