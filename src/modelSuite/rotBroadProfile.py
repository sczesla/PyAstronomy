import numpy as np
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl

class RotBroadProfile(fuf.OneDFit):
  """
    Implements rotational broadening with linear limb-darkening.
    
    *Fit Parameters*:
      - `xmax`  - Maximal extent of the profile
      - `eps`   - Linear limb-darkening coefficient
      - `A`     - Area under the profile (negative for absorption)
      - `off`   - An offset applied to the profile
      - `lin`   - Gradient of a linear term to adjust the 'continuum'
      - `mu`    - Center of the profile (same units as `xmax`)
      - `gsig`  - The standard deviation of a Gaussian with which the
                  rotational profile is convoluted, e.g., to model
                  instrumental resolution.

    The profile is given by:
    
    .. math::
        
        G(x) = A \\left(c_1\\sqrt(1-(x/x_{max})^2) + c_2(1-(x/x_{max})^2)\\right) + off + lin \\cdot x
    
    with the constants given by:
    
    .. math::
    
        c_1 = \\frac{2(1-\\epsilon)}{\\pi x_{max} (1-\\epsilon/3)} \\;\\;\\;\\;
        c_2 = \\frac{\\epsilon}{2 x_{max} (1-\\epsilon/3)}
    
    Here, `x` can either denote a velocity or a wavelength shift. Thus, `xmax` can
    be given in km/s or Angstrom depending on the input; see, e.g., "Stellar Photospheres"
    by D.F. Gray for a derivation.
    
    Note that the profile is normalized, i.e., the area amounts to one. This may, however,
    depend on the sampling. If the profile is undersampled, errors can become significant.
  """

  def __init__(self):
    fuf.OneDFit.__init__(self, ["off", "lin", "xmax", "eps", "A", "mu", "gsig"], rootName="RTB")
    self["A"] = 1.0
    self["xmax"] = 1.0
  
  def _rotProf(self, v):
    """
      Calculate the rotationally broadened line-profile.
      
      Parameters
      ----------
      v : array
          Abscissa at which to evaluate profile.
      
      Returns
      -------
      Profile : array
          The rotationally broadened profile.
    """
    denom = np.pi*self["xmax"]*(1.0-self["eps"]/3.0)
    c1 = 2.0*(1.0-self["eps"])/denom
    c2 = (np.pi*self["eps"]/2.0) / denom
    
    indi = np.where(self["xmax"] - np.abs(v-self["mu"]) > 0.0)[0]
    y = np.zeros(v.size)
    
    vv = ((v[indi]-self["mu"])/self["xmax"])**2
    y[indi] += (c1 * np.sqrt(1.0 - vv) + c2 * (1.0 - vv))
    y[indi] *= self["A"]
    return y
  
  def evaluate(self, v):
    """
      Calculates the rotational broadening profile according to current parameter values.

      Parameters:
      x : array
          Wavelength or velocity.
    """
    y = self._rotProf(v)
    if self["gsig"] != 0.0:
      # Convolve with Gaussian to implement instrumental broadening
      y = pyasl.broadGaussFast(v, y, self["gsig"], edgeHandling="firstlast")
    y += self["off"] + (self["lin"] * v)
    return y