import numpy as np
from PyAstronomy import funcFit as fuf

class RotBroadProfile(fuf.OneDFit):
  """
    Implements rotational broadening with linear limb-darkening.
    
    *Fit Parameters*:
      - `xmax`  - Maximal extent of the profile
      - `eps`   - Linear limb-darkening coefficient
      - `A`     - Area under the profile (negative for absorption)
      - `off`   - An offset applied to the profile
      - `lin`   - Gradient of a linear term to adjust the 'continuum'

    The profile is given by:
    
    .. math::
        
        G(x) = A \\left(c_1\\sqrt(1-(x/x_{max})^2) + c_2(1-(x/x_{max})^2)\\right) + off + lin \\cdot x
    
    with the constants given by:
    
    .. math::
    
        c_1 = \\frac{2(1-\\epsilon)}{\\pi x_{max} (1-\\epsilon/3)} \\;\\;\\;\\;
        c_2 = \\frac{\\epsilon}{2 x_{max} (1-\\epsilon/3)}
    
    Here, `x` can either denote a velocity or a wavelength shift. Thus, `xmax` can
    be given in km/s or Angstrom depending on the input. 
    
    Note that the profile is normalized, i.e., the area amounts to one. This may, however,
    depend on the sampling. If the profile is undersampled, errors can become significant.
  """

  def __init__(self):
    fuf.OneDFit.__init__(self, ["off", "lin", "xmax", "eps", "A"], rootName="RTB")
    self["A"] = 1.0
    self["xmax"] = 1.0
  
  def evaluate(self, v):
    """
      Calculates the rotational broadening profile according to current parameter values.

      Parameters:
      x : array
          Wavelength or velocity.
    """
    denom = np.pi*self["xmax"]*(1.0-self["eps"]/3.0)
    c1 = 2.0*(1.0-self["eps"])/denom
    c2 = (np.pi*self["eps"]/2.0) / denom
    
    indi = np.where(self["xmax"] - np.abs(v) > 0.0)[0]
    y = np.zeros(v.size)
    
    vv = (v[indi]/self["xmax"])**2
    y[indi] += (c1 * np.sqrt(1.0 - vv) + c2 * (1.0 - vv))
    y[indi] *= self["A"]
    
    y += self["off"] + (self["lin"] * v)
    
    return y